# dynamics

An implementation of Roy Featherstone's spatial-vector rigid-body dynamics
algorithms (RBDA), aimed at being fast enough for real-time use. Featherstone's
formulation uses 6x6 spatial transforms, but for single-DOF joints most of
that matrix is structurally zero or one — this project generates C code that
exploits that sparsity directly instead of doing dense 6x6 math at runtime.

The eventual goal is to JIT-compile the generated C with LLVM and call it
directly, rather than shelling out to a C compiler.

## Layout

- `models/` — mechanism definition files (`.mech`), the single source of
  truth for what gets generated and validated. A simple line-based text
  format (documented in `src/model/model.h`) describing gravity and, per
  body: parent, joint type and axis, tree translation, mass, CG, inertia,
  and whether to support an external force. The coordinate and sign
  conventions are spelled out under [Conventions](#conventions) below. The
  grammar covers revolute, prismatic, cylindrical, spherical, and floating
  joints; the pipeline implements all single-DOF joints (revolute and
  prismatic about x, y, or z) plus genuine 6-DOF `floating` bodies (RBDA
  sec 9.4), which have no gimbal lock at any attitude and take 7 position
  coordinates but 6 velocity coordinates (nq ≠ nv; the packed layouts are
  given by the model's qidx/vidx maps). Spherical and cylindrical joints
  are rejected with a clear error until they're implemented. The example
  set spans the standard textbook topologies, every one generated,
  validated, and benchmarked by `build.bat`:
  - `threelink` / `fourlink` — the original planar chains
  - `elbow` — 6R spatial arm with spherical wrist (Spong/Craig); maximally
    serial, mixed z/y/x axes, dense H
  - `stanford` — Stanford arm, the classic serial chain with a prismatic
    boom (R-R-P plus 3R wrist)
  - `gantry` — Cartesian PPP robot; no trigonometry survives generation
    and H is constant (~19 ns/call, ~36x over dense)
  - `pendulums` — four independent spatial pendulums off the base;
    maximally parallel, H is diagonal by branch-induced sparsity
  - `bintree` — Featherstone's balanced binary tree (7 bodies); block
    sparsity between the two subtrees
  - `satellite` — free-flying bus (`joint floating`) carrying two solar
    panel arms and an antenna boom (RBDA ch. 9 floating base); nq = 10
    (quaternion + position + 3 angles), nv = 9. In free fall the generated
    code reproduces |linear qdd| = 9.81 with zero angular and joint
    acceleration at every attitude — including the ones that were
    ill-conditioned for the earlier Euler-chain formulation
  - `fourbar` — planar four-bar linkage, the canonical closed kinematic
    loop (RBDA ch. 8): crank–coupler–rocker spanning tree, loop closed
    back to the ground pivot with a `point-on-line` joint (see
    [Closed loops](#closed-kinematic-loops) below)
  - `rrrs` — spatial RRRS loop: four revolute tree joints with mixed axes
    and 3D offsets, closed by a `spherical` loop joint; properly
    constrained with nc = 3 = rank(K), mobility 1
  - `slidercrank` — crank and slider as separate branches off the base,
    tied by a rigid rod modeled as a `distance` loop joint (nc = 1); the
    only model whose two loop paths meet at the base with both nonempty
  - `quartercar` — HMMWV front quarter-car test rig, built from the HMMWV
    suspension data shipped with Project Chrono: sprung mass on a vertical
    slider, double-wishbone suspension, and a driven shaker post under the
    wheel as the vertical position input. The spindle carries ball joints
    everywhere, so it is a `floating` body pinned by loop constraints —
    two `spherical` closures to the wishbones, a `distance` radius rod,
    and a `distance` "tire" rod down to the post (nq = 11, nv = 10,
    nc = 8, mobility 2) — the first model to combine a floating body with
    kinematic loops
- `src/model` — the shared model-file parser (plain C, header-only), used
  by the generator, the validator, and the benchmark, so all three always
  agree on the mechanism.
- `src/codegen` — the active approach. A small C++ program that loads a
  mechanism file and, instead of computing with numbers, computes with
  strings: each spatial-vector/transform component is a symbolic expression.
  `compose()`, `multFactors()`, `addTerms()` etc. simplify away structural
  zeros and ones as they go, then `emit()` prints each simplified expression
  as a named C temporary and rewires later composition to reference that name
  (so expressions never re-expand across pipeline stages). Running the
  generator produces `build/dynamics.c` — the specialized, sparsity-exploiting
  C implementation for the mechanism it was given.
- `src/validate` — a dense-matrix validation reference. Implements the same
  forward dynamics (RNEA bias force + CRBA mass matrix + solve) with full
  6x6 spatial matrices built straight from the RBDA equations — explicit
  Plucker transforms composed by plain matrix multiplication, `mcI` spatial
  inertias, `crm`/`crf` cross-product operator matrices, and a generic
  Gaussian solver for force coordinate changes. It shares no code with the
  generator, then runs both implementations on 1000 random states (q, qd,
  tau, external forces) and compares qdd to 1e-9. `build.bat` runs it on
  every build and fails the build on disagreement. On its first run it
  caught a real sign error in the generator's spatial-inertia lower-left
  block that the earlier hand-picked tests couldn't see.
- `src/chrono` — an external cross-check against Project Chrono. Reads the
  same `.mech` file and rebuilds the mechanism as a Chrono multibody system
  (maximal-coordinate bodies + joint constraints). Unlike `src/validate`,
  Chrono shares no conventions with this project, so this checks the physics
  of the `.mech` conventions themselves, not just the code generation — and
  it earned its keep immediately by catching a joint-convention bug that had
  been invisible since the start of the project (see
  [Conventions](#conventions)). Each binary links one model's generated
  `dynamics_<model>.c`, so Chrono validates the artifact this project
  actually ships rather than reaching it transitively through the dense
  reference. Three checks run per model:
  - **`qd == dq/dt`** — a pure kinematic assertion, no dynamics and no
    Chrono solver, that the integration rule and the coordinates agree.
    A sign inversion shows up as a scaled error of exactly 2.0.
  - **accelerations** — Chrono vs both the dense reference and the
    generated code, on 200 staged-then-random states (q, qd, tau, fext),
    agreeing to ~1e-9 on every model.
  - **simulation loop** — integrates the generated dynamics with the
    documented rule while time-stepping Chrono alongside it from the same
    initial state, comparing where the bodies physically end up. This is
    the only test that exercises `q`, `qd` and the integration rule
    *together over time*; the acceleration check evaluates `f(q, qd)` at a
    fixed state and never forms `dq/dt`, which is exactly how the old
    convention bug survived. The two sides use different first-order
    integrators, so the run is repeated at `dt/2` to show the gap halves —
    a real discrepancy would not shrink with the step. Typical result:
    ~1e-7 to 1e-6 over 0.25 s at dt = 2e-5, with a ratio of 2.00.

  One Chrono-specific note is encoded in its comments: `DoAssembly()` gets
  accelerations by finite-differencing a small step, which is only ~1e-5
  accurate and drops the centrifugal `Cq_dot*v` terms, so the checker
  assembles and solves the acceleration-level KKT system exactly instead.
  `--sim <T> <dt>` runs just the simulation comparison at a chosen horizon
  and step. Built and run by `build_chrono.bat` (needs the Chrono SDK and
  Eigen; see below), separate from `build.bat` so the main build has no
  Chrono dependency.
- `src/bench` — benchmark comparing the generated sparse code against the
  dense reference (both built with `/O2`), cycling through a pool of random
  states with a checksum to defeat dead-code elimination. Run by `build.bat`
  after validation; on the original dev machine the 3-link chain runs at
  ~50 ns/call generated vs ~800 ns/call dense (~16x).
- `src/numeric` — an earlier, currently non-building attempt at a runtime
  (non-generated) numeric implementation of the same algorithms. Kept for
  reference; `system.h` was removed mid-refactor and never replaced, so this
  does not currently compile. `build.bat` still attempts it and reports
  failure without stopping the rest of the build. It follows spatial_v2
  directly rather than building `E` as a matrix — the transforms are
  hand-written sparse applications with the sin/cos entries baked into the
  arithmetic — and it already carries the `E = Rᵀ`, `S = +a` convention
  documented above, so it needs no sign changes if it is ever revived. Note
  that it is a stub in other ways: the CRBA ancestor loop is empty and
  `x_fapply_TxRz` is a copy of `x_mapply_TxRz` rather than the force
  transform, so it would need real work beyond restoring `system.h`.
- `build/` — all build output: object files, executables, and generated
  sources such as `dynamics.c`. Not checked into git.
- `references/` — course notes, slides, and the `spatial_v2` MATLAB
  reference implementation this project is checked against. Not checked into
  git (see `.gitignore`).

## Closed kinematic loops

Mechanisms may contain kinematic loops (RBDA ch. 8). The bodies always form
a spanning tree; a `loop` block in the `.mech` file adds a loop-closing
joint between two existing bodies (or a body and the base). Loop joints
impose **constraints, never coordinates** — nq/nv are untouched — and the
dynamics becomes method 2 of RBDA sec 8.5: with the loop Jacobian K and the
acceleration-level RHS k + kstab (Table 8.4, computed in base coordinates),
solve

```
A lambda = k + kstab - K qdd0        A = K H^-1 K^T,  qdd0 = H^-1 (tau - C)
H qdd    = tau - C + K^T lambda
```

One Gaussian elimination of H serves qdd0 and all columns of H^-1 K^T
(`solveLinearSystemM` in the generated code), and everything before the
solve — K's epsilon structure (only tree joints on the loop's two root
paths get entries, RBDA eq 8.18–8.20), the Table 8.1 position-error
formulae, the axis alignment — is symbolic, so loop constraints get the
same structural-sparsity treatment as the tree dynamics.

Loop joint types and the constraints they impose: `revolute`/`prismatic`
(5), `cylindrical` (4), `spherical` (3), `point-on-line` (2, the book's
sphere-in-cylinder), `distance` (1, a rigid rod of given `length` between
an anchor point on each body). Two properties worth knowing:

- **Loops must be properly constrained: rank(K) = nc.** Any planar loop
  closed by a revolute is overconstrained (nc = 5, rank 2 — RBDA sec 8.10),
  which makes A singular; the validator computes rank(K) at every random
  trial and fails the build with a pointer to the fix: close planar loops
  with `point-on-line`, which constrains exactly the two translations
  perpendicular to its axis. Isolated variable-mobility singularities
  (e.g. a four-bar with all links collinear) are fine — random states hit
  them with probability zero, which is also why loop-model trials never
  sit exactly at a designed pose like q = 0.
- **Baumgarte stabilization** (RBDA sec 8.3) with alpha = beta = 1/Tstab is
  baked into the generated k + kstab; `stab <Tstab>` in the model file
  (default 0.1 s) sets the error-decay time constant. Torn-open states are
  pulled back onto the closure manifold exponentially; at consistent
  states the terms vanish.
- **Loop models are validated at states near the closure manifold**,
  generated by a random walk shared by both harnesses (`manifoldSeed` /
  `manifoldStep` in `src/validate/dense.h`): keep an on-manifold state
  across trials, take a moderate random tangent step each trial, and
  Gauss-Newton-project back. Fully random q — what tree models use — tears
  the loops open by meters, and at such states A = K H⁻¹ Kᵀ can be almost
  arbitrarily ill-conditioned (the quartercar model hits cond(A) ~ 1e10
  within 1000 draws, turning roundoff into 1e-7 gen-vs-dense noise that
  says nothing about correctness); near the manifold cond(A) stays that of
  the physical mechanism and the 1e-9 gate is meaningful. Each validation
  state then gets a small deliberate tear on top so the Baumgarte terms
  stay exercised (at exactly-closed states they vanish). The projector's
  Newton update lives in velocity space, mapped into q through the
  quaternion tangent map for floating bodies — the same rates `integrateQ`
  uses — so loop models may contain floating bodies (the quartercar's
  spindle is one: a body whose every physical joint is a ball must be
  floating, since spherical tree joints take 3 DOF the single-DOF pipeline
  doesn't have).

`point-on-line`'s constraint-force basis T is fixed to the *predecessor*
frame (its constraint directions follow the predecessor's axis), unlike the
Table 8.1 types whose T is successor-fixed; RBDA eq 8.47 for k is unchanged
either way by cross-product antisymmetry — see `loopJcalc` in
`src/validate/dense.h`, whose header comment carries the derivation. Its
position error is exact (recovered from rho~ = -E^T B), not small-angle.

`distance` is fixed to *neither* body — its constraint direction
rho = p_s − p_p turns with the mechanism — yet it still fits the single-T
machinery: the constraint force is a pure force along the line through both
anchor points, so the one spatial force T = (p_p × rho; rho) =
(p_s × rho; rho) serves predecessor and successor alike and the K path walk
is unchanged. Its position error is carried in the squared form
delta = (|rho|² − d²)/2, which is polynomial in the state (no square roots
in the generated code) and whose gradient is exactly the K row, so closure
projection works unchanged; the acceleration RHS replaces Table 8.4's
vs × vp term with the rod's own direction rate,
rho·(w_s × ṗ_s − w_p × ṗ_p) + |ρ̇|². Both forms of the derivation (direct
differentiation of ½ rho·rho and Ṫ·(v_s − v_p)) agree — see the branch in
`gen_loops.h` / `fd_dense_full`.

The Chrono cross-check builds loop joints as Chrono constraint links
(point-on-line's free axis is marker x there, not z; `distance` maps to
`ChLinkDistance`) and compares at states from the same on-manifold random
walk the validator uses — comparing constrained dynamics at torn-open
states would just compare stabilization schemes, so unlike the validator
it adds no tear. All four loop models agree with Chrono to ~1e-9 on
accelerations and track its trajectories with the expected first-order
truncation ratio.

## Force elements and points of interest (planned)

Designed but not yet implemented — the full implementation plan lives in
[PLAN-force-elements.md](PLAN-force-elements.md). The design was settled
around two principles the project already follows: **declare what you pay
for** (nothing appears in the generated code that the `.mech` file didn't
ask for — the `fext 1` contract, extended), and **the dynamics API stays
the stable contract** while new machinery grows alongside it.

The decisions:

- **Points of interest (POIs)** — named points fixed to a body, declared in
  the `.mech` file (`poi <name> <body> <x y z>`); the generator emits their
  world kinematics. A TSDA anchor *is* a POI, so anchors and observation
  points share one declaration and one symbolic treatment.
- **Force elements are first-class** — `tsda <name>` blocks (translational
  spring-damper-actuator between two named POIs) and `rsda <name>` blocks
  (rotational, acting on a named tree joint). A TSDA is generated
  symbolically by the same machinery as the `distance` loop joint: the
  applied generalized force `Q = -(f/l) * row` reuses the path-walk row, so
  it touches only the joints on the path between the two anchors — no
  6-vector spatial-force detour. An RSDA on a tree joint is a torque into
  that joint's `tau` slot with the angle being `q[i]` itself (RSDA across a
  loop-closing revolute is deferred; it needs an atan2 relative angle).
- **`fext` unchanged** — stays opt-in per body for genuinely external loads
  (tire contact, validation probes); elements never route through it.
- **Element inputs are opt-in** — every element carries a baked linear law
  (`k`, `c`, and a reference `l0`/`q0`) from the `.mech` file; elements
  additionally declared `input` get a slot in a compact `act[]` array, and
  the total force is `baked law + act[slot]`. Passive elements compile to
  pure expressions with no runtime input surface. The dynamics signature
  becomes `fwd_dyn_crba( q, qd, tau, fext, act, out_qdd )`, with `act`
  NULL-able exactly like `fext`. Caller-side nonlinear laws (real damper
  curves) use the pattern `k 0 c 0 input`: compute the force outside from
  the element's observed state, feed it through `act`.
- **Name lookup via runtime tables** — the generated file exports
  `fwd_dyn_find( kind, name )` (`kind` one of `"poi"`, `"tsda"`, `"rsda"`,
  `"act"`) backed by const string tables, resolved once at caller init.
  Runtime lookup rather than a generated header because the eventual JIT
  consumer is a fixed host program that cannot recompile per model.
- **Observation function** — generated `fwd_dyn_observe( q, qd, out_obs )`
  exposes POI world positions and velocities plus each element's state
  (`l`, `ldot` for TSDAs; angle, rate for RSDAs), indexed through the same
  lookup. The temporaries already exist inside the element generation, so
  the marginal cost is writes to an output array — and it keeps geometry
  out of every caller's force-law code.

## Conventions

These are the contract between a `.mech` file, the generated
`fwd_dyn_crba()`, and whatever integrates it. They are stated normatively in
`src/model/model.h`; this is the narrative version.

### Joint coordinates

Joint coordinates are the ordinary physical ones, and `qd` is the honest time
derivative of `q`:

| joint | `q` | `qd` | `tau` |
| --- | --- | --- | --- |
| revolute | right-handed rotation of `+q` about the `+axis` | `dq/dt` | couple about the `+axis` |
| prismatic | displacement along the `+axis` | `dq/dt` | force along the `+axis` |
| floating | 7 slots: unit quaternion (body-to-base attitude) then base-frame position | 6 slots: body's spatial velocity **in its own frame**, `[ω; v]` | `[moment about the body origin; force]`, body frame |

Gravity is the physical vector (e.g. `0 0 -9.81`), and `fext[i]` is a spatial
force on body `i` in **base** coordinates: `[moment about the base origin;
force]`. Bodies declared without `fext 1` never read their slots.

### Integration

Single-DOF joints integrate with plain `q += qd*dt`. A floating body is the
one exception — its `q` is a base-frame quaternion plus position while its
`qd` is a body-frame spatial velocity, so the rates are

```
Qdot = 0.5 * Q x (0, ω)      quaternion product, ω = qd[0..2], body frame
rdot = R(Q) * v              v = qd[3..5], body frame
```

That is the integrator's job, outside the generated function.
`integrateQ()` in `src/chrono/main.cpp` is a worked reference for both cases.

### Why the transpose is explicit

The Plücker `E` block is a **coordinate** transform (child-from-parent), not
the body's rotation: a child physically rotated by `R` has `E = Rᵀ`. That
transpose is forced by the math, so the only real choice is where to put it.
spatial_v2 folds it into `rotx`/`roty`/`rotz`, which is why its "rotation"
matrices read as the transpose of the textbook `R`. This project instead
writes the ordinary right-handed `R` and transposes it **explicitly**, at the
one place the coordinate transform is formed:

- `Xform::setJointRot()` — `src/codegen/xform.h`
- `jointXS()` — `src/validate/dense.h`
- the floating-body quaternion in `src/codegen/gen_fwd_dyn_crba.h`, where the
  transpose is a pure index swap on the emitted temps (free at runtime)

`src/numeric` differs in style but not in convention: it never builds `E` as a
matrix, hand-coding each sparse transform application with the sin/cos baked
into the arithmetic, so its `E = Rᵀ` is implicit in those expressions — called
out in a comment at the top of `src/numeric/xform.h`.

The motion subspace stays on the positive axis (`S = +a`). `E = Rᵀ` together
with `S = +a` is exactly what makes `q` the physical angle and `qd = dq/dt`.
In the generated code the transpose costs nothing: it folds into code
generation, so `dynamics.c` emits the same operation count either way.

This is worth stating because getting it wrong is nearly undetectable. Until
the Chrono cross-check was written, the generator used `R` directly as `E`,
which silently made `q` and `qd` opposite-handed: forward dynamics stayed
perfectly self-consistent — `src/validate` agreed to 1e-14, because the dense
reference and the generator shared the convention — but `q` ran backwards
relative to its own `qd`, so any integrator stepping `q += qd*dt` drove the
mechanism the wrong way. `checkQdConsistency()` in `src/chrono/main.cpp` now
asserts `qd == dq/dt` on every model, independently of Chrono, so the
invariant cannot silently regress. (A sign inversion shows up there as a
scaled error of exactly 2.0.)

## Building

Requires MSVC (Visual Studio with the C++ workload) on Windows. Both build
scripts locate Visual Studio via `vswhere` automatically, so they work from a
plain terminal — no Developer Command Prompt needed.

```
build.bat
```

For every model in its list (all twelve above), this builds
`build\codegen.exe`, generates `build\dynamics_<model>.c`, compiles
`build\validate_<model>.exe` and `build\bench_<model>.exe` against it, and
runs both — the build **fails** if any model's generated code disagrees with
the dense reference.

```
build_chrono.bat [model ...]
```

The Chrono cross-check, separate from `build.bat` so the main build has no
Chrono dependency. Self-contained: builds the code generator, emits
`dynamics_<model>.c` for each model, links it into a per-model
`build\chrono_check_<model>.exe`, and runs the three checks described in
`src/chrono` above. With no arguments it runs every model; pass model names
(e.g. `build_chrono.bat fourbar quartercar`) to check just those.

Requires a Project Chrono SDK and Eigen; defaults to `C:\dev\chrono` and
`C:\dev\eigen-5.0.0`, overridable via the `CHRONO_DIR` and `EIGEN_DIR`
environment variables. Note the SDK's required consumer flags (encoded in the
script): `/arch:AVX2` is mandatory because Chrono was built with AVX2 and
Eigen objects cross the DLL boundary — without it the heap corrupts.
`Chrono_core.dll` (in the SDK's `bin`) must be on `PATH` at run time; the
script sets this itself, but see the note under `chrono_check_<model>.exe`
below for running the binaries by hand.

## Running the executables

Everything lands in `build\`. All commands below are shown from the repo
root; the executables accept the model path as an argument, so the working
directory doesn't matter. (If you omit the model path, `validate` and `bench`
default to `../models/threelink.mech`, which only resolves when run from
inside `build\`.)

A note on the per-model naming: `validate_<model>.exe`,
`bench_<model>.exe`, and `chrono_check_<model>.exe` each have that model's
generated `dynamics_<model>.c` **compiled in**. The `.mech` path on the
command line tells the harness what to compare against, so it must be the
same model the binary was built for — a mismatch in nq/nv is caught with an
error (`... was generated for nq=3 nv=3 but 'models\elbow.mech' has nq=6
nv=6 -- regenerate`), but a different model with coincidentally equal sizes
would just fail validation.

### `codegen.exe` — generate specialized dynamics for a mechanism

```
build\codegen.exe <mechanism-file> [output.c]
```

Reads a `.mech` file and writes the specialized, sparsity-exploiting C
implementation (default output name `dynamics.c`, in the current directory).
Silent on success; parse errors and unsupported joint types (spherical or
cylindrical **tree** joints) are reported on stderr with a nonzero exit.

```
build\codegen.exe models\elbow.mech build\dynamics_elbow.c
```

The generated file exports one function and two size constants:

```c
void fwd_dyn_crba( double *q, double *qd, double *tau, double *fext, double *out_qdd );
extern const int fwd_dyn_crba_nq;   // length of q
extern const int fwd_dyn_crba_nv;   // length of qd, tau, out_qdd
```

`fext` is `6 * nbodies` doubles (base-coordinate spatial force per body; see
[Conventions](#conventions)). To use a mechanism of your own, drop a
`.mech` file in `models\` and add its name to the `MODELS` list at the top
of `build.bat` — it will be generated, validated, and benchmarked like the
others.

### `validate_<model>.exe` — dense reference vs generated code

```
build\validate_threelink.exe models\threelink.mech
```

Runs 1000 random states (near-manifold states for loop models) through both
the compiled-in generated code and the dense 6x6 reference and compares qdd
to 1e-9, also checking rank(K) = nc for loop models at every trial. Exit
code 0 = PASS, 1 = FAIL:

```
validate: models\threelink.mech: 1000 trials, max scaled |qdd_generated - qdd_dense| = 1.699e-14 (worst at trial 969)
validate: PASS
```

### `bench_<model>.exe` — timing, generated vs dense

```
build\bench_threelink.exe models\threelink.mech
```

Times both implementations (compiled `/O2`) over a pool of 256 random
states, batching until the per-call figure is stable, with a running
checksum to defeat dead-code elimination:

```
bench: models\threelink.mech (3 bodies, nq 3, nv 3)
bench: generated sparse      51.2 ns/call  (checksum/call +0.062690)
bench: dense reference      874.2 ns/call  (checksum/call +0.062686)
bench: speedup 17.1x
```

### `chrono_check_<model>.exe` — cross-check against Project Chrono

```
build\chrono_check_<model>.exe <model.mech>                 all three checks
build\chrono_check_<model>.exe <model.mech> --sim <T> <dt>  simulation comparison only
```

Requires `Chrono_core.dll` on `PATH` when run outside `build_chrono.bat`:

```
set PATH=C:\dev\chrono\bin;%PATH%
build\chrono_check_threelink.exe models\threelink.mech
```

The no-flag form runs the qd-vs-dq/dt consistency check, the 200-trial
acceleration cross-check, and the default simulation comparison
(T = 0.25 s, dt = 2e-5, repeated at dt/2 to show the gap is integrator
truncation):

```
chrono_check: models\threelink.mech: qd-vs-dq/dt consistency, max scaled error = 5.799e-10 (worst at trial 11)
chrono_check: models\threelink.mech: 200 trials, max scaled |qdd_chrono - qdd_dense| = 9.261e-10 (worst at trial 192)
chrono_check: models\threelink.mech: 200 trials, max scaled |qdd_chrono - qdd_generated| = 9.261e-10
chrono_check: PASS
chrono_check: models\threelink.mech: simulation loop, generated dynamics vs Chrono, T=0.25 dt=2e-05
    t= 0.031  gen b2( 2.8723  0.7004  0.0000)  chrono b2( 2.8723  0.7004 -0.0000)  max|dp|=1.07e-08
    ...
    t= 0.250  gen b2( 2.9017  0.6093  0.0000)  chrono b2( 2.9017  0.6093 -0.0000)  max|dp|=8.56e-08
chrono_check: max |position gap| = 8.561e-08 at dt, 4.286e-08 at dt/2  (ratio 2.00 -- shrinks with the step, so the gap is integrator truncation)
chrono_check: SIM PASS (tolerance 1.0e-02)
```

`--sim` skips the first two checks and runs just the simulation comparison
at a chosen horizon and step, printing the trajectory table of both sides —
useful for watching a single model evolve:

```
build\chrono_check_elbow.exe models\elbow.mech --sim 0.25 2e-5
```