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

Requires MSVC (Visual Studio with the C++ workload) on Windows.

```
build.bat
```

This builds `build/codegen.exe` from `src/codegen`, runs it on
`models/threelink.mech` to emit `build/dynamics.c`, builds and runs
`build/validate.exe` (dense reference vs generated code — the build fails if
they disagree), repeats generation + validation for `models/fourlink.mech`,
builds and runs `build/bench.exe` (`/O2` timing of both implementations),
and attempts to build `build/numeric.exe` from `src/numeric` (currently
expected to fail — see above).

To generate dynamics for a different mechanism:

```
build\codegen.exe models\mymodel.mech mydynamics.c
build\validate.exe models\mymodel.mech      (after compiling it against mydynamics.c)
```

### Chrono cross-check

```
build_chrono.bat
```

Self-contained: builds the code generator, emits `dynamics_<model>.c` for each
model, links it into a per-model `build\chrono_check_<model>.exe`, and runs the
three checks above. Requires a Project Chrono SDK and Eigen; defaults to
`C:\dev\chrono` and `C:\dev\eigen-5.0.0`, overridable via the `CHRONO_DIR` and
`EIGEN_DIR` environment variables. Note the SDK's required consumer flags
(encoded in the script): `/arch:AVX2` is mandatory because Chrono was built
with AVX2 and Eigen objects cross the DLL boundary — without it the heap
corrupts. `Chrono_core.dll` (in the SDK's `bin`) must be on `PATH` at run time.

To watch a single model run, with a trajectory table of both sides:

```
build\chrono_check_elbow.exe models\elbow.mech --sim 0.25 2e-5
```
