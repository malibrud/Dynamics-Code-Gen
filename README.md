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
  and whether to support an external force. The grammar covers revolute,
  prismatic, cylindrical, spherical, and floating joints; the pipeline
  implements all single-DOF joints (revolute and prismatic about x, y, or
  z) plus genuine 6-DOF `floating` bodies (RBDA sec 9.4): orientation is a
  unit quaternion — normalized on input, no gimbal lock at any attitude —
  so a floating body takes 7 position coordinates but 6 velocity
  coordinates (nq ≠ nv; the packed layouts are given by the model's
  qidx/vidx maps, and the quaternion-rate integration q̇ = G(q)ω is the
  integrator's job, outside the generated function). Spherical and
  cylindrical joints are rejected with a clear error until they're
  implemented. The example set spans the standard textbook topologies,
  every one generated, validated, and benchmarked by `build.bat`:
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
- `src/bench` — benchmark comparing the generated sparse code against the
  dense reference (both built with `/O2`), cycling through a pool of random
  states with a checksum to defeat dead-code elimination. Run by `build.bat`
  after validation; on the original dev machine the 3-link chain runs at
  ~50 ns/call generated vs ~800 ns/call dense (~16x).
- `src/numeric` — an earlier, currently non-building attempt at a runtime
  (non-generated) numeric implementation of the same algorithms. Kept for
  reference; `system.h` was removed mid-refactor and never replaced, so this
  does not currently compile. `build.bat` still attempts it and reports
  failure without stopping the rest of the build.
- `build/` — all build output: object files, executables, and generated
  sources such as `dynamics.c`. Not checked into git.
- `references/` — course notes, slides, and the `spatial_v2` MATLAB
  reference implementation this project is checked against. Not checked into
  git (see `.gitignore`).

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
