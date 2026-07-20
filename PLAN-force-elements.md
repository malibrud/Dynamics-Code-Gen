# Implementation plan: force elements (TSDA/RSDA) and points of interest

Audience: an implementing agent working in this repo. Read the README's
"Force elements and points of interest (planned)" section first — it states
the settled design decisions; this file turns them into staged, verifiable
work. Do not revisit the decisions.

## Ground rules

- **Work in stages, in order.** Each stage ends with `build.bat` passing
  (all twelve models generate, validate to 1e-9, and benchmark). Commit each
  completed stage directly to `master` (this is a solo repo; no feature
  branches) with a message describing the stage.
- **Conventions are load-bearing.** `q` is the physical coordinate and
  `qd = dq/dt`; the Plücker `E` block is `R^T`, transposed explicitly at the
  one place the transform is formed (see README "Why the transpose is
  explicit"). Never introduce a hidden transpose. The Chrono cross-check
  exists precisely because self-consistent sign errors are invisible to the
  dense validator.
- **Two independent formulations per feature.** The generator applies
  element forces via the symbolic path-walk row; the dense reference must
  apply them via a *different* route (equal-and-opposite spatial forces) so
  agreement actually validates the math, not a shared bug.
- **The dense reference and generator share only `src/model/model.h`.**
  Keep it that way.
- Chrono SDK at `C:\dev\chrono`, Eigen at `C:\dev\eigen-5.0.0`
  (overridable via `CHRONO_DIR`/`EIGEN_DIR`). `/arch:AVX2` is mandatory for
  anything linking Chrono — without it the heap corrupts. `build_chrono.bat`
  handles all of this; don't hand-roll compile lines.
- `C:\dev\chrono\data\vehicle\hmmwv` has real HMMWV suspension data
  (masses, hardpoints, spring/damper rates) — the source for stage 6's
  quartercar numbers, and `models/quartercar.mech` documents which files
  its existing numbers came from.

## Target API (final state)

Generated per-model C file exports:

```c
// Existing, with one new parameter. act may be NULL when the model declares
// no input elements (same contract as fext).
void fwd_dyn_crba( double *q, double *qd, double *tau, double *fext,
                   double *act, double *out_qdd );
const int fwd_dyn_crba_nbod;
const int fwd_dyn_crba_nq;
const int fwd_dyn_crba_nv;
const int fwd_dyn_crba_nact;    // number of act[] slots (input elements)

// New: world kinematics of POIs and element states. out_obs layout:
// all POIs in declaration order (6 doubles each: px py pz vx vy vz),
// then all TSDAs (2 each: l ldot), then all RSDAs (2 each: angle rate).
void fwd_dyn_observe( double *q, double *qd, double *out_obs );
const int fwd_dyn_observe_len;

// New: runtime name lookup. kind is "poi" | "tsda" | "rsda" (returns the
// element's starting offset into out_obs) or "act" (returns the index into
// act[]). Returns -1 if not found. Plain strcmp over generated const
// tables; callers resolve once at init.
int fwd_dyn_find( const char *kind, const char *name );
```

## `.mech` grammar additions

After the bodies (interleaving with loop blocks is fine; a `poi`, `tsda`,
or `rsda` keyword ends any open body/loop block, exactly as `body` and
`loop` end each other today):

```
poi <name> <body> <x> <y> <z>    # one line; named point at (x y z) in body
                                 # <body>'s frame; body -1 = fixed base

tsda <name>                      # translational spring-damper-actuator
  from <poi-name>                # anchor 1 (order defines rho = p_to - p_from)
  to <poi-name>                  # anchor 2
  k <k>                          # spring rate, N/m       (default 0)
  c <c>                          # damping, N/(m/s)       (default 0)
  l0 <l0>                        # free length, m         (default 0)
  input                          # optional: reserve an act[] slot

rsda <name>                      # rotational spring-damper-actuator
  joint <body>                   # acts on the tree joint of body <body>
  k <k>                          # Nm/rad (or N/m on a prismatic joint)
  c <c>
  q0 <q0>                        # rest coordinate        (default 0)
  input
```

Force law (signs are part of the contract — document them in `model.h`):

- TSDA: scalar tension `f = k*(l - l0) + c*ldot + act[slot]`, positive `f`
  pulls the two anchors together. `l = |rho|`, `rho = p_to - p_from`,
  `ldot = (rho . rhodot)/l`.
- RSDA: generalized force on the joint coordinate
  `tau_e = -k*(q - q0) - c*qd + act[slot]`, positive `act` drives `+q`.

Parser validation (all with file/line context, matching existing style):

- POI names unique; element names unique within their kind. Length-limit
  names (`MODEL_NAME_LEN 32`) and count-limit (`MODEL_MAX_POI 16`,
  `MODEL_MAX_ELEM 16`).
- `poi` body in `-1 .. nbod-1`; declared only after all its body's `body`
  block (POIs may simply be required to come after all bodies — simplest).
- `tsda`: both POI names must resolve; `from != to`; reject both anchors on
  the same body at the identical local point (guaranteed l = 0 forever;
  `f/l` divides by `l`, so a coincident-anchor state is a model error —
  also note this hazard in the `model.h` header comment).
- `rsda`: `joint` must be a single-DOF tree joint (revolute or prismatic);
  reject floating/spherical/cylindrical.
- Assign `act` slots in declaration order across both element kinds;
  store `nact` on the Model.

`Model` struct additions (`src/model/model.h`):

```c
typedef struct { char name[32]; int body; double p[3]; } ModelPoi;
typedef enum { ELEM_TSDA, ELEM_RSDA } ModelElemType;
typedef struct {
    char name[32];
    ModelElemType type;
    int poi1, poi2;    // TSDA: indices into poi[]
    int joint;         // RSDA: body index of the tree joint
    double k, c, ref;  // ref = l0 (TSDA) or q0 (RSDA)
    int hasInput;
    int actIdx;        // slot in act[], -1 if !hasInput
} ModelElem;
// Model gains: int npoi; ModelPoi poi[MODEL_MAX_POI];
//              int nelem; ModelElem elem[MODEL_MAX_ELEM]; int nact;
```

## The math (generator side)

A TSDA is the `distance` loop joint's geometry with a known force instead
of an unknown constraint force. `gen_loops.h`'s `J_DISTANCE` branch already
emits, symbolically: anchor world positions (`loopAnchorBase`), `rho`,
anchor velocities `pdot = v0 + w x p`, `rhod`, and the spatial force
`T = (p x rho; rho)`, then walks the two root paths projecting `T` onto
each tree joint's base-coordinate motion subspace with sign `-` on the
predecessor(`from`) path and `+` on the successor(`to`) path. Call that
per-coordinate projection `row[r]`.

`row` is the gradient of `delta = (rho.rho)/2`, i.e. `row[r] = l * dl/dq_r`.
With potential `V = k*(l-l0)^2/2` and Rayleigh damping, the generalized
force of tension `f` is

```
Q[r] = -(f/l) * row[r]        f = k*(l - l0) + c*ldot + act[slot]
```

which is exact, not small-angle. The only new symbolic quantities beyond
what the distance branch computes are `l = sqrt(rho.rho)` (the first sqrt
in generated element code — fine; the floating-body quaternion normalize
already uses sqrt), `linv = 1.0/l`, and `ldot = (rho.rhod)*linv`.

An RSDA needs no geometry at all: `Q[vidx[joint]] += tau_e` with
`q = q[qidx[joint]]`, `qd = qd[vidx[joint]]`.

Application point: the solve (`gen_solve_qdd.h`) currently forms
`tau[r] - C_r`. Element forces join as `tau[r] - C_r + Q_r`, where `Q_r`
temps exist only for coordinates on some element's path (structural zeros
stay out, same discipline as `Kentry`).

## Stages

### Stage 1 — parser, model struct, grammar docs

Files: `src/model/model.h`.

- Add the structs, parse the grammar above, run the validations, assign
  `qidx`-style layout: `obsIdx` per POI/element (offset into `out_obs`)
  and `actIdx` per input element, plus `npoi/nelem/nact` and
  `obsLen = 6*npoi + 2*nelem`.
- Extend the header comment (the format documentation) with the new
  grammar, force-law signs, and the coincident-anchor hazard.
- No consumer reads the new fields yet. Gate: `build.bat` passes untouched;
  hand-test a few malformed files against the parser errors (a throwaway
  `.mech` in the scratchpad, driven through `codegen.exe`, is enough).

### Stage 2 — generator refactor: shared kinematics emission (no behavior change)

Files: `src/codegen/gen_fwd_dyn_crba.h`, new `src/codegen/gen_kinematics.h`
(or similar), `src/codegen/main.cpp`.

- `fwd_dyn_observe` needs the same outward kinematics pass (per-body
  `xFromParent`/`xFromBase`, `S`, `vJ`) inside a *second* generated
  function. Extract the per-body kinematics emission loop of
  `gen_fwd_dyn_crba()` into a helper that emits into whatever function
  body is currently open, so both generated functions can call it. The
  observation function will also need body spatial velocities `sys->v[]`
  (currently populated by `gen_rnea_C`) — extract the velocity recursion
  (the RNEA outward pass's velocity part) into the shared helper too, or
  give the observe generator its own small velocity emission; choose
  whichever keeps `gen_rnea_C` byte-identical in output.
- Gate: regenerate all twelve models and diff against pre-refactor output —
  `dynamics_<model>.c` must be **byte-identical**. Keep the diff step in
  your verification, not in the build.

### Stage 3 — element forces in the generator + `act` in the signature

Files: new `src/codegen/gen_elems.h`, `src/codegen/system.h` (Q-entry
storage, element data), `src/codegen/main.cpp` (populate System from
Model, call the new generator), `src/codegen/gen_fwd_dyn_crba.h`
(signature + `fwd_dyn_crba_nact`), `src/codegen/gen_solve_qdd.h`
(`+ Q_r` terms), and — same stage, so the build stays green — every
caller: `src/validate/main.c`, `src/bench/main.c`, `src/chrono/main.cpp`
(pass `NULL` for now; real support in stages 5–6).

- `gen_elems` runs after `gen_rnea_C` (needs `sys->v[]` and
  `sys->xFromBase[]`), before `gen_solve_qdd`. For each TSDA follow the
  `J_DISTANCE` pattern in `gen_loops.h` closely: `loopAnchorBase` for
  anchors, `emitVec3` discipline for temps, `invertXform` +
  `fwdXformMotion` to get base-coordinate body velocities (note
  `gen_loops.h:197` — the fixed base has `v = 0`, and for elements there
  is no `avp` subtlety because no accelerations are needed). Emit `l`,
  `linv`, `ldot`, `f` (folding `act[ slot ]` in only when the element has
  one; when `k` and `c` are both zero the baked part vanishes
  symbolically), then the path walk producing `Q_<r>` temps accumulated
  into per-coordinate expressions in `sys`.
- RSDA: emit `Q_<r>` for the one slot directly from `q`/`qd` temps already
  unpacked in the kinematics pass.
- `gen_solve_qdd`: rhs entries become `tau[ r ] - C<r> + Q_<r>` where a
  `Q_<r>` exists (both the loop and no-loop branches).
- The generated header comment for `act` mirrors the `fext` one: slots,
  NULL-ability, sign conventions.
- Gate: `build.bat` passes (existing models have no elements — their
  generated output should differ only by the added `act` parameter and
  `fwd_dyn_crba_nact = 0`).

### Stage 4 — dense-reference element forces + test models + validation

Files: `src/validate/dense.h`, `src/validate/main.c`, `src/bench/main.c`,
new models, `build.bat` model list.

- In `fd_dense_full`, apply each element as **equal-and-opposite spatial
  forces in base coordinates** — deliberately NOT the row route: for a
  TSDA with tension `f`, the force on the `to` body is `-f*rho/l` at
  `p_to` (spatial: `[ p_to x F; F ]`), and the negation at `p_from` on the
  `from` body; add them into the same internal per-body applied-force
  accumulation the `fext` path uses (`dense.h:453-469`), bypassing the
  per-body `hasFext` flag (element forces are always on). For an RSDA add
  `tau_e` into the effective tau of its coordinate. Compute anchor
  world positions/velocities from the `Xa` (base->body) transforms the
  dense pass already builds.
- Add two models:
  - `models/springpend.mech` (or similar): a small tree mechanism — e.g.
    two of the `pendulums` bodies — with a `tsda` from a body POI to a
    base POI (nonzero `k`, `c`, `l0`, plus `input`), and an `rsda` with
    `input` on one joint. Cheap, loop-free, exercises every element code
    path including `act`.
  - Extend `models/rrrs.mech` or `models/slidercrank.mech` with a passive
    TSDA (elements + loops together, floating-free), OR add elements to
    `models/quartercar.mech` with placeholder rates (real rates come in
    stage 6). At least one model must combine elements with loops.
- `src/validate/main.c` and `src/bench/main.c`: generate random `act`
  values per trial/state (`5.0 * frand()`, like tau) whenever
  `fwd_dyn_crba_nact > 0`, and pass them to both implementations.
- Add the new model(s) to `build.bat`'s `MODELS` list.
- Gate: `build.bat` green — the two application routes agreeing to 1e-9 on
  1000 random states is the core correctness evidence for the whole
  feature. If they disagree, suspect the sign of the path walk first
  (`-` on `from`, `+` on `to`) and the `f/l` factor second.

### Stage 5 — observation function + name tables + lookup

Files: new `src/codegen/gen_observe.h`, `src/codegen/main.cpp`,
`src/validate/main.c` (spot-check).

- `fwd_dyn_observe( q, qd, out_obs )`: open a second function after
  `fwd_dyn_crba` is closed, run the shared kinematics + velocity emission
  from stage 2 (temp names may repeat — separate function scope), then for
  each POI emit world position (`loopAnchorBase` pattern) and velocity
  (`v0 + w x p`), for each TSDA emit `l`, `ldot`, for each RSDA copy
  `q`/`qd` of its joint; store into `out_obs[]` at the layout offsets from
  stage 1. Emit `fwd_dyn_observe_len`.
- `fwd_dyn_find`: emit const name/offset tables and a strcmp loop; kinds
  `"poi"`, `"tsda"`, `"rsda"` return `out_obs` offsets, `"act"` returns
  the act slot. Unknown kind or name returns -1.
- Validation: in `src/validate/main.c`, for models with elements, at each
  trial call `fwd_dyn_observe` and check (a) TSDA `l`/`ldot` against a
  finite-difference of `l` across a small `q` integration step or against
  a dense recomputation from the `Xa` transforms (dense recomputation is
  cleaner — add a tiny helper in `dense.h`), to 1e-9; (b) `fwd_dyn_find`
  round-trips every declared name and rejects a bogus one.
- Gate: `build.bat` green including the new checks.

### Stage 6 — Chrono cross-check + real quartercar spring/damper

Files: `src/chrono/main.cpp`, `models/quartercar.mech`.

- Map elements onto Chrono: `ChLinkTSDA` (set spring coefficient, damping
  coefficient, rest length; attach at the two anchor world points) and
  `ChLinkRSDA` (about the tree joint's axis between child and parent).
  For `input` elements, register a force/torque functor that returns the
  full law `±(k, c, ref)` plus a per-trial `act` value stored in the
  `ChronoModel` struct — check the exact functor API in
  `C:\dev\chrono\include\chrono\physics\ChLinkTSDA.h` / `ChLinkRSDA.h`
  rather than trusting recalled signatures, and mind Chrono's sign and
  rest-length conventions: get one element working on a trivial model
  first if signs look wrong. The existing acceleration check solves a
  custom KKT system (`solveAccelerations`) — confirm applied link forces
  from TSDA/RSDA links are included in the load assembly it uses; if not,
  fall back to applying the element forces as accumulated body forces on
  the Chrono side (compute them from the same state, exactly as the dense
  reference does).
- `crossCheck` and `runSim`/`stepGenerated`: thread `act` through
  (random per trial in the cross-check when `nact > 0`; zero in the sim
  loop so both sides agree by construction).
- Upgrade `models/quartercar.mech`: replace/augment with the real HMMWV
  spring and damper as a `tsda` between the lower wishbone (or wherever
  the Chrono HMMWV data attaches it — read the data files) and the
  chassis, rates from `C:\dev\chrono\data\vehicle\hmmwv`. Keep `stab`,
  masses, and hardpoints consistent with the comments already in the
  model file. The shaker input stays the driven prismatic post.
- Gate: `build_chrono.bat` green on all models — accelerations to ~1e-9,
  sim-loop position gap halving at `dt/2` (ratio ~2.00). The quartercar
  now validates spring, damper, loops, and a floating body together.

### Stage 7 — documentation close-out

Files: `README.md`, `PLAN-force-elements.md`.

- Move the README's "Force elements and points of interest (planned)"
  section to implemented status: drop "(planned)", fold in anything
  learned (actual model list, measured benchmark deltas, Chrono mapping
  notes), document the new executables' behavior changes in the "Running
  the executables" section (the `act` parameter, `fwd_dyn_observe`,
  `fwd_dyn_find`), and update the Layout section (new `gen_elems.h` /
  `gen_observe.h`, new models).
- Delete this plan file in the same commit (git history keeps it).

## Verification summary (what "done" means)

1. `build.bat`: all models (twelve existing + new element models) generate,
   validate at 1e-9 over 1000 trials with random `act`, and benchmark.
2. The dense reference validates the generator through a *different*
   force-application route (spatial forces vs path-walk rows).
3. `fwd_dyn_observe` agrees with dense-side recomputation at 1e-9;
   `fwd_dyn_find` round-trips all names.
4. `build_chrono.bat`: all models pass all three checks — Chrono shares no
   conventions with this project, so this is the physics-level check on
   the element force laws and their signs.
5. Generated output for pre-existing element-free models differs from
   before the feature only in the `act` parameter and `nact` export.
