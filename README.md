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
    and H is constant (~10 ns/call, ~69x over dense)
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
    only model whose two loop paths meet at the base with both nonempty.
    Also carries a passive return-spring `tsda` from a fixed base point to
    the crank, combining a force element with a kinematic loop
  - `quartercar` — HMMWV front quarter-car test rig, built from the HMMWV
    suspension data shipped with Project Chrono: sprung mass on a vertical
    slider, double-wishbone suspension, and a driven shaker post under the
    wheel as the vertical position input. The spindle carries ball joints
    everywhere, so it is a `floating` body pinned by loop constraints —
    two `spherical` closures to the wishbones, a `distance` radius rod,
    and a `distance` "tire" rod down to the post (nq = 11, nv = 10,
    nc = 8, mobility 2) — the first model to combine a floating body with
    kinematic loops. Carries the real HMMWV front coil spring and shock as
    two `tsda` elements (see [Force elements and points of
    interest](#force-elements-and-points-of-interest) below)
  - `springpend` — two independent pendulums off the base, loop-free;
    carries a `tsda` and an `rsda`, both declared `input`, exercising every
    force-element code path (POI kinematics, path-walk with `k`/`c`/`act`)
    without a kinematic loop's added complexity
  - `halfcar` — HMMWV rear half-car: a sprung half-mass free to heave and
    roll (chained `prismatic z` then `revolute x`, since this format has no
    combined 2-DOF joint — a massless-ish carrier body realizes the pair,
    same trick as `elbow`'s wrist chain), carrying left and right
    double-wishbone corners built the same way as `quartercar`'s (each
    spindle a `floating` body pinned by two `spherical` closures and a
    `distance` "tierod", plus a `distance` "tire" rod to its own driven
    shaker post) — the rear is used rather than the front because its real
    tierod grounds on the chassis at both ends (no steering rack), so it is
    an ordinary fixed-length link. Real HMMWV rear coil springs and shocks
    (two `tsda` per side) connect chassis_roll to each LCA; unlike
    `quartercar`, both TSDA anchors attach to the body that physically
    carries them rather than the fixed base, since with a roll DOF present
    a base-anchored mount would pick up a spurious length change from pure
    heave/roll (nq = 22, nv = 20, nc = 16, mobility 4 — heave, roll, and the
    two posts)
  - `hmmwv` — the whole HMMWV: `quartercar` and `halfcar` assembled into one
    vehicle. A `floating` chassis (the full 2086.52 kg sprung mass, free-flying
    with all six rigid-body DOF — heave, roll, pitch, yaw, and lateral/
    longitudinal drift — quaternion attitude so no gimbal lock) carries four
    independent double-wishbone corners, each built exactly as `quartercar`'s
    (a `floating` spindle pinned by two `spherical` closures, a `distance`
    tierod, and a `distance` "tire" rod to its own driven shaker post). The two
    front tierods reroute from the chassis onto a laterally-sliding
    `prismatic y` steering rack (HMMWV rack-and-pinion), so rack travel steers
    the front wheels about their kingpin axes — the steering geometry; the rear
    tierods stay grounded on the chassis (the rear doesn't steer). All four
    real HMMWV coil springs and shocks are carried as two `tsda` per corner
    (nq = 48, nv = 43, nc = 32, mobility 11 — the free chassis's 6, the
    steering rack, and the four posts; 16 loops, the first model past the old
    8-loop cap, which this raised to 20). One of the two largest models here,
    and among those that most exercise the loop solver — 32 constraints, full
    rank at every trial
  - `hmmwv_pitman` — the same full HMMWV, but with the vehicle's REAL steering:
    a Pitman-arm / idler linkage instead of `hmmwv`'s rack-and-pinion (Chrono's
    default HMMWV uses `HMMWV_PitmanArm.json`). A `revolute z` Pitman arm on
    the chassis (the steering input) swings a relay rod — a `revolute z` body
    hung off the arm's end, idealizing the universal joint — whose far end is
    held by a `distance` "idler" loop (length = the idler-arm length), so the
    Pitman and idler arms form a parallelogram and the relay translates to
    steer both front tie rods. The steering chain is mobility 1 (Pitman 1 +
    relay 1 − idler distance 1); everything else matches `hmmwv` (nq = 49,
    nv = 44, nc = 33, mobility 11; 17 loops — the most of any model). Shows the
    format expressing a genuine spatial steering linkage from tree revolutes
    plus loop closures
  - `track_min` — the smallest model that exercises **road-wheel face
    contact**: one band, one road wheel, 16 shoes, everything at `y = 0` and
    every dimension chosen so the design pose *is* the equilibrium pose to
    within the static deflection. `m113` carries 882 pin pairs and 504 ground
    pads; this carries 48 and 64, with one face pair loaded at a time, so a
    failing pair can be printed by hand. This is the model to reach for when a
    contact gate trips
  - `track_demo` — a tracked chassis (floating hull, sprocket, idler, two road
    wheels) carrying a 40-shoe band: exercises the `track` block end to end
    (see [Track bands](#track-bands) below) with more than one road wheel.
    Placeholder geometry; the point is the expansion. Its band is centred on
    the hull, because a single band offset to one side is a monocycle stood on
    its edge and never reaches a static attitude
  - `m113` — the M113 armored personnel carrier, the first tracked vehicle
    here and by far the largest model: a `floating` hull (7819.24 kg) carrying
    two track assemblies, each a drive sprocket, a translational idler on a
    fore/aft-sliding tensioned carrier, five torsion-bar road-wheel stations
    (trailing arm + road wheel), and a closed band of 63 shoes riding all
    seven rollers. 27 hand-written bodies expand to 153 once the two bands are
    generated (nq = 171, nv = 168, nc = 4 — just the two band closures —
    mobility 164). Real M113 data from Chrono's `data/vehicle/M113/`
    (single-pin track variant), which the mass total corroborates: 11 325 kg
    against the real vehicle's ~11 340 kg combat weight. 18 force elements —
    ten torsion bars as `rsda` on the arm joints, six shocks and two idler
    tensioners as `tsda` — with Chrono's preloads folded into the rest
    coordinate (`q0 = preload/k` for the bars, `l0 = free length + preload/k`
    for the tensioner), which is exact for a linear law. Also the model that
    made the validator's gate size-aware (see `src/validate` below)
- `src/model` — the shared model-file parser (plain C, header-only), used
  by the generator, the validator, and the benchmark, so all three always
  agree on the mechanism. `kinematics.h` alongside it is the shared spatial
  algebra and forward-kinematics layer — `plux`, `jointXS`, `mcI`,
  `xformOrigin`, `denseKinematics` — extracted from `src/validate/dense.h`
  once it acquired a second consumer that has nothing to do with validation
  (the contact engine needs per-body poses and velocities, and a force
  producer must not depend on the validation reference).
- `src/codegen` — the active approach. A small C++ program that loads a
  mechanism file and, instead of computing with numbers, computes with
  strings: each spatial-vector/transform component is a symbolic expression.
  `compose()`, `multFactors()`, `addTerms()` etc. simplify away structural
  zeros and ones as they go, then `emit()` prints each simplified expression
  as a named C temporary and rewires later composition to reference that name
  (so expressions never re-expand across pipeline stages). Running the
  generator produces `build/dynamics.c` — the specialized, sparsity-exploiting
  C implementation for the mechanism it was given. `gen_kinematics.h` factors
  the per-body kinematics/velocity emission out of the main dynamics pass so
  it can be reused verbatim by a second generated function; `gen_elems.h`
  emits each force element's generalized force (folded into the solve as
  `+ Q_r`); `gen_observe.h` emits that second function, `fwd_dyn_observe`,
  plus the `fwd_dyn_find` name-lookup tables (see [Force elements and points
  of interest](#force-elements-and-points-of-interest)). `gen_solve_qdd.h`
  closes the function with the `H qdd = tau - C` solve, whose sparse
  factorization `gen_prelude.h` emits at file scope (see [The H
  solve](#the-h-solve)).
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

  It also carries a second, independent forward dynamics: `fd_dense_aba`,
  Featherstone's ch. 7 articulated-body algorithm, split into a config-only
  factorization reused across right-hand sides and a cheap per-RHS solve —
  the property that makes the method-2 loop solve O(n) for a long band, and
  16x faster than CRBA at 192 shoes. Every trial cross-checks it against the
  CRBA path on a tree-only, element-free copy of the model.

  **The gate is size-aware above nv = 64**: `1e-9 * (nv/64)^4`. Both
  implementations end in a factorization of H — Gaussian elimination on the
  dense side, `L^T D L` on the generated one (see [The H solve](#the-h-solve))
  — so the noise floor is `eps * cond(H)` either way, and for a long serial
  chain cond(H) grows like the fourth
  power of its length — arithmetic, not correctness. On `m113` the dense CRBA
  and the dense ABA, which share no code and never touch the generator,
  disagree by 2.5e-9, and that gap tracks the chain length to the fourth
  across 15 / 31 / 63 shoes (3.0e-12 / 1.5e-10 / 2.5e-9). Holding a fixed
  1e-9 there would gate on conditioning, the same mistake near-manifold
  sampling exists to avoid at the loop level. Every model but `m113` has
  nv <= 44, so all of them are still gated at exactly 1e-9; and a real bug
  shows up as a large *relative* error (the metric divides by |qdd|), not as
  a few ulps of solve noise.
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
    ~1e-7 to 1e-6 over 0.25 s at dt = 2e-5, with a ratio of 2.00. The gap is
    the worse of each body's **position** (m) and **attitude** (rad) error;
    comparing the two on one scale is crude, but both are O(1) here and the
    metric only has to catch a divergence. Position alone was a blind spot: a
    body that rotates about its own frame origin — `springpend`'s pendulums,
    for instance — never moves, so no amount of rotational disagreement
    registered.

  One Chrono-specific note is encoded in its comments: `DoAssembly()` gets
  accelerations by finite-differencing a small step, which is only ~1e-5
  accurate and drops the centrifugal `Cq_dot*v` terms, so the checker
  assembles and solves the acceleration-level KKT system exactly instead.
  `--sim <T> <dt>` runs just the simulation comparison at a chosen horizon
  and step. Built and run by `build_chrono.bat` (needs the Chrono SDK and
  Eigen; see below), separate from `build.bat` so the main build has no
  Chrono dependency.

  Force elements map onto Chrono unevenly, both notes earned empirically
  since the SDK ships as headers and prebuilt binaries with no `.cpp` source
  to read the default force law from (see [Force elements and points of
  interest](#force-elements-and-points-of-interest)): TSDA maps onto
  `ChLinkTSDA` (its `k`/`c`/rest length match this project's directly; its
  actuator force needs negating). RSDA does **not** map onto `ChLinkRSDA` —
  that class's angle only tracks continuity across a stepped simulation, and
  every cross-check trial here places bodies directly at a random `q`, which
  can land Chrono's angle a multiple of 2π away from the true one with no way
  to disambiguate the branch. RSDA is instead computed directly from this
  project's own `q`/`qd` and applied as an ordinary joint force, exactly like
  `tau[]`.

  That last point has a sharp edge, found by `m113`. A TSDA is a `ChLinkTSDA`,
  so Chrono re-evaluates its force law every step on its own; an RSDA is a
  plain accumulator entry, so it stays **frozen** at whatever was last written
  — and the simulation loop wrote it once, before stepping. Chrono was
  therefore integrating a constant torque while the generated side applied the
  real `-k(q - q0) - c*qd`: two different systems, showing up as a trajectory
  gap that does not shrink with the step (`m113` diverged 0.688 m at a ratio of
  1.00, worst at a road wheel hanging off a preloaded torsion bar). The loop
  now recomputes each RSDA torque from **Chrono's own** joint state before
  every step. Reading a joint coordinate back out of Chrono is exact rather
  than approximate because a `.mech` tree joint carries no fixed rotation, so
  the child differs from the parent by a pure rotation (or translation) about
  the axis; the revolute angle comes back wrapped to (-π, π] and is unwrapped
  against the previous step, which is legitimate here — the very branch
  ambiguity that rules `ChLinkRSDA` out of the staged-state checks is harmless
  in a stepped simulation that moves the joint by a hair per step.

  The bug had been latent since RSDAs landed. `springpend`, the only other
  RSDA model, hid it twice over: its `q0` is 0, so the frozen torque happened
  to be the correct one, and its pendulum rotates about its own frame origin,
  which the then position-only gap metric could not see either way. Both holes
  are now closed.
- `src/bench` — benchmark comparing the generated sparse code against the
  dense reference (both built with `/O2`), cycling through a pool of random
  states with a checksum to defeat dead-code elimination. Run by `build.bat`
  after validation; on the original dev machine the 3-link chain runs at
  ~42 ns/call generated vs ~880 ns/call dense (~21x). The full table across
  every model is under [The H solve](#the-h-solve).
- `src/contact` — two contact force producers and their shared test harness:
  `track_contact.h` (shoe pins against the rollers they ride on) and
  `ground_contact.h` (shoe pads against the terrain), the latter reading its
  surface through `terrain.h`. Both are caller-side and stateless: they read a
  state, write `fext`, and the generated dynamics never learns contact exists.
  They are deliberately **independent of each other** — the shared surface is a
  damping formula and a regularized-Coulomb expression, which is not worth
  giving up the property that a mistake in one cannot be present in the other.
  See [Track contact](#track-contact) and [Ground contact](#ground-contact)
  below. The harness links no generated code — it tests the force laws, not the
  dynamics — so `build.bat` runs it on every model, and models without a band or
  without a ground assert they produce no forces.
- `src/sim` — the trajectory harness: a plain-C stepper (`integrate.h`, with
  the quaternion update ported out of the Chrono cross-check so it needs no
  SDK) plus the checks that only a trajectory can make. Before this, nothing in
  `build.bat` stepped time at all — `validate` and `bench` evaluate at fixed
  states, and the only integrator in the repo lived inside the Chrono harness.
  Two schemes live here, semi-implicit Euler and Störmer–Verlet, and tree models
  are gated on **both** — Euler's energy drift must scale as `dt`, Verlet's as
  `dt²`. All eight report 2.000 and 4.000 to three digits, `satellite` included,
  which is the only thing that exercises the quaternion path. Track models get
  the contact-in-the-loop settling run instead.
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
- `tools/` — the one thing here that is not C: `plot_trace.py`, which plots
  the time histories `sim_<model>.exe --trace` writes. See
  [Track contact](#track-contact).
- `references/` — course notes, slides, and the `spatial_v2` MATLAB
  reference implementation this project is checked against. Not checked into
  git (see `.gitignore`).

## The H solve

Forward dynamics ends in `H qdd = tau - C`, and for a while that solve was the
one part of the generated code that ignored everything the generator knew. The
CRBA emitted `H`'s branch-induced structural zeros symbolically and then handed
the matrix to a dense Gaussian elimination with partial pivoting — `O(n^3)`, and
blind to the sparsity. On `m113` that was a 168x168 elimination, ~1.6M
multiply-adds, and it had grown to dominate: the model ran at 433 us/call, only
4.4x over the dense reference, while `threelink` was managing 17x. Essentially
all of the 74k generated lines were optimizing the cheap part.

The solve is now Featherstone's sparse **`L^T D L` factorization** (RBDA sec 6.5,
Tables 6.3–6.5), emitted in `gen_prelude.h` as `ltdlFactor` / `ltdlSolveM` and
called from `gen_solve_qdd.h`. Three things make it fit this project unusually
well:

- **`H` is symmetric positive definite**, so the partial pivoting was pure
  waste — no pivot search, no row swaps, half the flops.
- **The factorization has zero fill-in.** Given the *expanded parent array*
  `lam` (RBDA Table 6.4 — the parent relation over velocity coordinates rather
  than bodies: within a floating body's 6-wide block the coordinates chain, and
  the block's first coordinate points at the parent body's last), `L` carries
  exactly `H`'s sparsity pattern. So the branch-induced zeros are not merely
  skipped, they are provably never touched, and the cost falls from `n^3/3` to
  `sum_i depth(i)^2 / 2`. `lam` is emitted as a `static const int Hlam[]` table
  beside `H`; the generator asserts it reproduces `crbaHasEntry()`'s pattern on
  every coordinate pair, so the two descriptions of the same sparsity cannot
  drift apart.
- **One factorization still serves every right-hand side**, which is what the
  method-2 loop solve needs for its `1 + nc` columns (see [Closed kinematic
  loops](#closed-kinematic-loops)). Only the `nc x nc` constraint system
  `A = K H^-1 K^T` still goes through a dense solver: it has no tree structure
  to exploit, and `nc <= 33` across every model here.

Table 6.3 gives both `L^T L` and `L^T D L`, and the choice between them was
settled by measurement, not preference. `L^T L` needs a square root per
coordinate, and that was enough to make it measurably **slower** than the dense
elimination it replaced on every model with no branch sparsity to win back —
`elbow` 163 -> 172 ns, `rrrs` 194 -> 216 ns. `L^T D L`'s `L` is unit
lower-triangular, so no square roots appear at all, and the regressions became
speedups.

Generated `ns/call`, before and after, across the whole model set (same machine,
same `/O2`, `x` = speedup over the dense reference):

| model | nv | before | after | | model | nv | before | after |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| threelink | 3 | 53.1 | **41.8** (21x) | | slidercrank | 2 | 34.5 | **23.1** (49x) |
| fourlink | 4 | 80.7 | **66.0** (18x) | | springpend | 2 | 25.6 | **15.5** (33x) |
| elbow | 6 | 163.5 | **139.3** (12x) | | fourbar | 3 | 80.6 | **64.1** (29x) |
| stanford | 6 | 158.4 | **137.6** (12x) | | rrrs | 4 | 194.4 | **180.6** (13x) |
| gantry | 3 | 19.2 | **10.1** (69x) | | quartercar | 10 | 759 | **539** (17x) |
| pendulums | 4 | 48.6 | **20.1** (36x) | | halfcar | 20 | 3747 | **1719** (19x) |
| bintree | 7 | 197.7 | **129.1** (15x) | | hmmwv | 43 | 26795 | **9865** (27x) |
| satellite | 9 | 216.8 | **166.0** (7.7x) | | hmmwv_pitman | 44 | 56442 | **11097** (27x) |
| track_demo | 27 | 2933 | **1640** (9.3x) | | m113 | 168 | 432934 | **148075** (13x) |

(The `m113` figure predates [Track contact](#track-contact), which added
`fext 1` to four roller bodies and with it four more inverse-transform blocks
inside an already 74k-line function. Measured A/B on one machine, interleaved,
min of four: 137.6 us without them against 150.7 us with, so **+9.5%** — and the
variance widens too, 137-141 against 150-173, which is what extra straight-line
code does to the instruction cache. That is the standing price of routing
contact reactions through the bearing into the hull rather than lumping them
into `tau`; the rest of the table is untouched and was measured in a single
earlier session, so it is left as it was rather than half-rewritten.)

(The `track_demo` row also predates [Ground contact](#ground-contact), which
raised that model from 12 shoes to 40 — 17 bodies to 45 — because its sprocket
and idler could not otherwise carry load. Its generated call went 1.6 -> 22.9 us
with it. That is the model getting bigger, not the solver getting worse.)

(`track_min` postdates the table entirely: 20 bodies, nv 30, **2.69 us**
generated against 18.9 us dense, 7.0x. Its whole contact-in-the-loop step is
5.8 us — 172k steps/s — which is what makes it usable for the kind of repeated
settle runs that debugging a contact law needs.)

**Where a contact step actually goes.** `bench_<model>.exe` reports a step
budget for any model with a band or a ground, since the producers are
caller-side and a real step costs `denseKinematics` + both producers +
`fwd_dyn_crba`:

| | m113 | track_demo |
| --- | ---: | ---: |
| `denseKinematics` | 14.0 us (7.2%) | 4.0 us (13.2%) |
| `trackContactFext` | 7.1 us (3.6%) | 1.4 us (4.6%) |
| `groundContactFext` | 5.7 us (2.9%) | 2.1 us (6.8%) |
| `fwd_dyn_crba` | 167.5 us (86.2%) | 22.9 us (75.5%) |
| **full step** | **194 us**, contact 13.8% | **30 us**, contact 24.5% |

Two things worth reading off that. **Contact is not the bottleneck** — the
generated solve outweighs both engines and the kinematics they share by better
than 5:1 on `m113`. And **ground contact cost nothing in the generated code**,
unlike track contact's 9.5%: it needed `fext 1` only on the shoes, which already
had it, so no new inverse-transform blocks were emitted.

Both engines run about 8–13 ns per pair or per point. They are all-pairs with no
culling — `m113` tests all 63x7x2 = 882 shoe-pin-against-roller pairs to find
the ~61 actually pressed into a cylinder, and all 504 pad corners to find the
~147 actually on the ground (both counts are pairs and corners, not rollers or
shoes: the model declares only 14 rollers in total) — so the cost is set by the declared counts rather than by
what is touching, which is also why a random state is a fair place to measure
them. A broadphase would save ~4% of the step and cost a persistent contact set,
so it has not been written.

The real lever is **`dt`, not the per-step cost.** Regularized Coulomb has
tangential slope `mu*f_n/v_eps`, an order of magnitude stiffer than the normal
damper, and it holds `m113` to 6.2e-5 s where the band alone ran happily at
1e-4 — so `src/sim` derives its step from that rather than hard-coding one.
`m113` settles at about 0.35x realtime; recovering the full step would take it
to ~0.57x, a bigger win than deleting contact entirely.

Every model got faster; the gains track how branchy each one is. `m113` is 2.9x
(433 -> 148 us, ~6.8 kHz), `hmmwv_pitman` 5.1x. The generated file also shrank
~75 KB, since only `H`'s lower triangle is now written — the factorization and
both substitutions stay below the diagonal, so the upper half is never read.

**What this does not fix:** `L^T D L` exploits *branch* sparsity only. Down a
single serial chain `lam[i] = i - 1`, `H` is full, and the algorithm degenerates
to a dense factorization — which is why `rrrs` and `elbow` gain the least.
`m113` wins because it is three independent branches off the base (each track
band's anchor shoe is its own root, so the two 68-coordinate band blocks are
mutually zero) and not because its 63-long chains got any cheaper. Making a long
chain itself `O(n)` needs the articulated-body algorithm — see `fd_dense_aba` in
`src/validate`, which is already split into a config-only factorization and a
cheap per-RHS solve for exactly that purpose.

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

One factorization of H serves qdd0 and all columns of H^-1 K^T (`ltlFactor`
then a single `ltlSolveM` over all 1 + nc right-hand sides in the generated
code — see [The H solve](#the-h-solve) below), and everything before the
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

## Track bands

A tracked vehicle's band is tens of near-identical shoes per side — far too
many to hand-write, and enough bodies that the algorithm choice starts to
matter. A `track` block declares the band once and the loader expands it:

```
track band_L
root 0                              hull body: the band's geometry frame
plane xz                            band plane (pin axis y); the only option
shoes 63                            N
pitch 0.154                         declared pin-to-pin (advisory, see below)
mass 18.02                          per shoe
cg 0.0752 0 0
inertia 0.22 0.04 0.25 0 0 0
fext 1                              per shoe, for the contact forces
roller 1 0.245 end drive teeth 10   sprocket: a path-defining end pulley
roller 3 0.255 end                  idler: the other one
roller 5 0.305                      road wheel, rides the band
ground z -0.55                      contact plane
kn 2e5   kt 2e5   mu 0.8            contact material
```

Each block expands into N shoe **bodies** — a `floating` anchor shoe plus a
chain of N-1 `revolute y` pins — and **one** `point-on-line` closure loop
tying the last shoe's trailing pin back to the anchor's leading pin. That is
the planar-loop rule from above: a planar loop closed by a revolute would be
overconstrained, so the closure constrains only the two in-plane
translations and leaves the pin rotation free. The chain is parented to the
**base**, not the hull: like the quartercar's spindle, it is held in place by
contact rather than by a mount.

The band cannot sit at q = 0 — a closed chain is curved — so the loader
computes the design pose and writes it to `m->q0` (quaternion + position for
the anchor shoe, a turn angle per pin), which consumers seed their initial q
from. Everything else in the format is translation-only at the design pose;
`q0` exists for the band.

Two things to know:

- **The belt is the convex hull of every roller circle.** The path is the taut
  band around all of them — gift-wrapped by outer tangents, so it comes out as a
  sequence of tangent runs and wrap arcs, with shoes laid evenly by arc length.
  Every roller shapes the path, not just the two `end` pulleys, which is why
  every one of them has to be locatable in the root frame at the design pose.

  This used to wrap the two `end` pulleys only, and the difference was not
  cosmetic. A two-pulley belt runs straight from sprocket to idler, so `m113`'s
  five road wheels poked up to **0.22 m** through their own track before
  anything had moved, and the perimeter came out 9.474 m against
  63 x 0.154 = 9.702 (2.3% short — almost exactly the road-wheel wrap the
  idealization left out). Wrapping the hull closes both gaps: the design pose
  now sits within **0.009 m** of equilibrium, and `linkLen` agrees with the
  declared pitch and with the sprocket's tooth pitch. For two rollers the hull
  walk reproduces the old two-pulley path exactly (9.474466 against 9.474470),
  which is how it was checked.

  The one subtlety is the tie-break: `m113`'s five road wheels are collinear and
  equal-radius, so the gift-wrap has to prefer the *farthest* candidate among
  those tied on turn angle, or it staples itself to the nearest wheel and skips
  the rest.
- **`pitch` is advisory.** The link length actually used is perimeter / N,
  reported as `linkLen`, since the shoes must close on the path exactly.
- The contact data (`pad`, `kn`, `kt`, `mu`, `teeth`) is parsed and kept for the
  contact force engines, which live outside this loader — see
  [Track contact](#track-contact) and [Ground contact](#ground-contact).
- **Shoe pitch has to stay well under the `end` rollers' diameter.** A sprocket
  or idler is carried by the pins either side of it, so once the pitch reaches
  `2R` it passes clean between two consecutive pins and can never carry load,
  however stiff the contact is. `modelLoad` warns.

  This applies to the **`end` rollers only**. Road wheels roll on the shoe's flat
  inner face, which cannot be fallen between at any pitch, so warning about them
  is warning about the contact primitive rather than about the model. The warning
  used to cover every roller, and `track_demo` was rebuilt from 12 shoes to 40 to
  satisfy a constraint its road wheels never had. 40 is still right, but for the
  other reason: 12 gives a 0.464 m pitch against 0.50 m end pulleys.

## Track contact

A tracked vehicle's only path to the ground is through the band, so until
there were contact forces the `m113` model could not move: the hull was in
free fall, the band hung free, and each idler tensioner — 0.27 m inside its
free length with no track tension to react it — shoved its carrier rearward
at ~2.7e5 N. The `track` block had been carrying `kn`/`kt`/`mu`, a ground
plane, and each roller's `radius`/`isDrive`/`teeth` since it landed, parsed
and validated and then read by nothing but the two `end` pulleys.

This stage implements shoe-against-roller contact for all seven rollers per
band. Shoe-against-ground came next and is written up under
[Ground contact](#ground-contact); the note on gravity below is why none of
this needed it in order to be tested.

**A caller-side force producer, not a change to the dynamics.** The step is

```
denseKinematics( mdl, q, qd, X, v );          // src/model/kinematics.h
trackContactFext( mdl, prm, X, v, fext );     // src/contact/track_contact.h
fwd_dyn_crba( q, qd, tau, fext, act, qdd );   // unchanged
```

**Why penalty contact and not exact unilateral constraints.** The reason is
structural rather than numerical, and it is forced by the architecture: the
loop machinery is baked at code-generation time. `nc` and `cidx[]` are frozen
in the model, `K` is emitted as straight-line code, and `gen_solve_qdd.h`
emits a solve sized for a fixed `A = K H⁻¹ Kᵀ`. A contact set that changes
every step — pins entering and leaving the sprocket wrap — would have to add
and remove constraint rows at runtime, which cannot be done without
abandoning code generation entirely. `fext` is a dense array read at fixed
offsets, so it is structurally invariant. It is also what the model file
already assumed: `m113.mech` says contact forces "have to be supplied from
outside: apply them through `fext`".

**Stateless, deliberately.** `trackContactFext` is a pure function of
`(q, qd)`: no contact anchors, no active set, no memory between calls. That
is what lets it be exercised at 1000 random, off-trajectory states by the
same machinery the rest of the project uses. It is also why friction is
regularized Coulomb rather than a bristle model, and why **`kt` remains
unused**: `kt` is dimensionally a tangential stiffness belonging to a
persistent-anchor stick model, which could not be evaluated at a random state
at all. The honest cost is that a loaded stationary vehicle creeps at up to
`vEps` ≈ 1 cm/s — measured at 0.88 cm/s on a settled `m113`. `vEps` used to set
the step size as well; see [Ground contact](#ground-contact) for the implicit
treatment that decoupled the two.

**Two primitives, not one, and which applies is set by `end`.** Pins engage the
sprocket and idler; road wheels roll on the flat inner face of the shoe between
them. Modelling a road wheel as pin-against-circle gets two things wrong, both
of which showed up:

- At rest the pin polygon cuts a sagitta `R − √(R² − ds²/4)` **into** the wheel
  — 9.4 mm on `m113` — so the wheel is embedded in the surface it is meant to
  be rolling on, and the error grows as the square of the pitch.
- A wheel resting on points can fall between two of them once the pitch reaches
  `2R`. A wheel resting on a face cannot, at any pitch.

The face law is worked entirely in the **shoe's** frame, where the running
surface is flat and axis-aligned: the plane `z = padZ + thickness`, spanning
`x ∈ [0, linkLen]`, with the rollers on its `+z` side because the band wraps the
outside. `thickness` is a new `track` key, and it is declared alongside `pad`
rather than independently precisely so the two cannot drift: the face is derived
from the pad, at `padZ + thickness`.

**Shoes tile the band end to end, so the face weight is a partition of unity.**
A hard `x ∈ [0, linkLen]` test hands a wheel's load from one shoe to the next
discontinuously; a symmetric fade leaves a support dip at every pin, since shoe
*i+1*'s face *begins* where shoe *i*'s ends. Instead each face's weight reaches
`fade` past both of its own ends with complementary smoothsteps, so
`w_i(x) + w_{i+1}(x − linkLen) ≡ 1` across the joint and both are exactly ½ at
the pin. Support transfers smoothly and is never doubled or dropped.

That weight makes this a **gated potential** `U = ½·w(x)·k·δ²`, so `−∇U` carries
a tangential term `½·k·δ²·w′(x)` beside the normal one — structurally the same
trap as the tooth law's missing `g′(δ)·U₀`, in a new place. It is included, and
the harness finite-differences the whole gradient rather than waiting for a
trajectory to notice.

**The two halves act at different points, so they take different velocities.**
The conservative part is a function of the axle's position in the shoe frame and
acts **at the axle**, conjugate to the axle's velocity relative to the shoe.
Friction is a surface effect and acts **at the foot** — the only place it can
spin the wheel — conjugate to the surface slip there. A rolling wheel drives the
slip to zero while its axle is still moving quite fast, so these are genuinely
two different vectors and the law needs both.

**A sign that cost a debugging session.** Both of the face law's velocities are
the *roller* against the shoe, where the pin law's `u` is the *shoe* against the
roller — while both report the force *on the shoe*. So friction's coefficient is
**positive** here and negative there. Copying the pin law's line verbatim turns
the friction into an anti-damper: with force `f = sc·u_s` on the shoe and `−f` on
the roller, the pair does work at `−f·u_s`, so a negative `sc` *injects*
`μ·f_N·|u_s|` watts. On `track_min` that drove the band round under its own
friction and grew total energy by 5e6 J over a four-second settle.

Every check in the harness passed while it did so, because all of them drove
`trkPairForce` — the pin path, which road wheels no longer take. The passivity
check in particular reported `worst power 0.000e+00 W`. There are now three
checks that drive `trkFaceForce` directly (partition of unity, `F = ∇U`, and
passivity + cone), and the passivity one catches the bug at a *state*, naming
the position and both velocities, rather than after a trajectory.

**One pin primitive for the `end` rollers.** The shoe frame's *origin is a pin*
(`model.h`: "the shoe frame is at its own pin, the next pin is at +linkLen on
x"), so walking shoes `shoe0 .. shoe0+N-1` visits all N pins exactly once
with no local-point transform, and a pin's world velocity is exactly
`Eᵀ·v[3..5]` with no `ω × ρ` term. The belt path was itself built tangent to
circles of the declared radii, so the query is just pin-centre against a
circle of `roller.radius`, in the roller's own rotating frame.

The force law, in that frame (band plane x–z, axle y):

- **Radial**, one-sided because the band wraps the outside:
  `f_n = kn·δ + c_n·ramp(δ)·δ̇`, clamped at zero. The clamp is what stops the
  damper pulling a separating shoe back, and it is exactly what the passivity
  check would catch if it were removed.
- **Tooth engagement** on the drive sprocket:
  `f = −g(δ)·[ F_skip·sin(N(θ−θ₀)) + c·u_t ]`. Because `N = teeth` is an
  **integer**, the sine is automatically periodic at the pocket pitch `2π/N`
  and automatically continuous across `θ = ±π` — no phase wrapping, no
  `round()`, no branch, which is most of why this form was chosen over an
  explicit nearest-pocket law. It vanishes at the pocket centre (stable) and
  at the tooth tip (unstable — correct: a pin balanced on a tip has no
  circumferential force), so the taper at half pitch is built in; `|f| ≤
  F_skip` makes tooth skip fall out of the law rather than being modelled
  separately; and it is the gradient of a real potential, which is what makes
  the energy accounting meaningful. Two-sided by construction, which is what
  a pin sitting between two flanks needs to carry drive *and* overrun.
- **Getting the gradient right, in both coordinates.** `g(δ)·U₀(θ)` is a
  potential in *two* variables, and two separate mistakes are easy here. Both
  were made, and the trajectory energy check found both.
  - `−∇U` has a **radial** part `g′(δ)·U₀` alongside the circumferential one.
    Applying only the circumferential part leaves the law non-conservative — a
    pin entering the wrap at one phase and leaving at another closes a loop
    with net work — and the term is not small: `g′` reaches 1.5/`gateFade` =
    50/m against a tooth potential reaching `2·F_skip·R/N` = 980 J, so the
    omitted force is up to ~5e4 N, the same order as contact itself.
  - In polar coordinates the conservative tangential force is
    `−(1/ρ)·∂U₀/∂θ`, so the moment arm is **ρ, not R**. Dividing by `R` leaves
    a residual of order `δ/ρ`, which reaches 14% at the edge of the gate band.
    Using `ρ` costs nothing and keeps the peak exactly `F_skip` *on* the pitch
    circle, where `ρ = R` — which is where the unit test measures, so it stays
    exact.

  Cumulatively, on the `m113` settle: worst single-step energy rise 13.04 J →
  0.22 J → **0.102 J**, and settled interpenetration 0.029 → 0.019 → **0.0078
  m**.
- **Lateral guide.** Not optional. The band is a floating planar chain
  parented to the **base**, closed only by a point-on-line loop constraining
  two in-plane translations, so its rigid-body y-translation, roll and yaw are
  completely unconstrained — and gravity is z-only. With purely in-plane
  contact those three modes would receive no force at all and the band would
  walk off its wheels. The guide horns are the real mechanism and `d[1]` is
  already in hand.
- **Friction** on the non-drive rollers only. On the sprocket the tooth term
  *is* the circumferential force, and Coulomb on top would double-count drive.
- **Damping** is derived from a critical-damping fraction and the declared
  shoe mass rather than new format keys — dimensionally correct and it
  rescales automatically with `kn`. For `m113` (`kn` 2e5, 18.02 kg shoes):
  ω = 105 rad/s, `c_crit` = 3797, so `c_n` = 759 at the running ζ = 0.2.

Reactions go on as `[p × F; F]` and its exact negation, with the moment
computed **once** from the same `p` and the same `F`. That makes each pair
cancel exactly; the *total* cancels only to accumulation roundoff (4.3e-16
relative on `m113`), because the two halves land in different running sums and
floating-point addition is not associative.

Sprockets (bodies 1/14) and idler wheels (3/16) gained `fext 1` — without it
their slots are never read, the reaction is silently dropped, and the vehicle
would gain momentum from nothing. Lumping the reaction into `tau` instead
would be wrong for the same reason: only the y-moment maps onto a joint
coordinate, so the radial bearing load carrying band tension into the hull
would vanish. `tau[6]`/`tau[19]` remain the drive input; the contact engine
adds the band's reaction on top.

That correctness is not free, and the cost is bigger than it looks: each
`fext 1` body emits its own `invertXform` plus a 6x6 force transform, and four
more of those inside `m113`'s 74k-line function cost **9.5%** (137.6 -> 150.7
us, interleaved A/B, min of four) — far more than the ~100 flops involved,
because what actually hurts is instruction-cache pressure. `trackContactCheckModel`
refuses to run on a model whose roller bodies lack the flag, since the failure
mode is otherwise invisible: the wrench still sums to zero in the `fext` array
while quietly never reaching the dynamics.

**Three things about the design pose.** All three were consequences of the belt
wrapping only the two `end` pulleys; the hull generalisation has since closed the
second and third, and they are kept here because they are what the design-pose
regression is shaped around:

- Pins do **not** sit exactly on the pitch circle. The loader takes each
  shoe's direction from the chord between pins but its length from the *arc*
  spacing `ds`, and on a wrap arc the true chord `2R·sin(ds/2R)` = 0.148040 m
  is shorter than `ds` = 0.150388 m, so every link overshoots by ~2.3 mm and
  the reconstructed pins drift a few millimetres *outside*. The radial force
  is zero there because the one-sided clamp fires, not because δ = 0 — so the
  regression asserts `δ ≤ 0`, never `|δ| ≈ 0`.
- The five road wheels per side used to poke **0.17–0.25 m** through the
  design-pose lower run. Wrapping the hull brings that to **0.0087 m**, which
  is what lets the gain ramp be short: it existed to relax a band jammed into
  its own rollers, and while the gains are low the whole vehicle free-falls, so
  0.2 s of ramp drops it 0.196 m — almost exactly the depth the wheels were then
  found buried at. The ramp is now 0.02 s and the drop is sub-millimetre.
- The sprocket's tooth arc pitch `2πR/N` = 0.153938 m is *exactly* the real
  0.154 m shoe pitch, and `linkLen` used to come out 2.36% short for the same
  wrap-shortfall reason. The hull belt makes the two agree. `modelLoad` warns
  when a drive sprocket's tooth pitch is more than 5% off the shoe pitch, and
  separately when the belt perimeter over `shoes` lands more than 5% off the
  declared `pitch` — that second warning only became meaningful once the belt
  wrapped everything, and it is what says the shoe count and the roller layout
  disagree about how long the track is.

**A stiffness consequence of the face primitive, worth stating separately.**
`kn` is per contact pair, and the two primitives share load very differently. A
pin wrap engages ~10 pins at once and they split the roller's load; the running
face is a partition of unity whose weights sum to exactly **1**, so one shoe's
stiffness carries a whole road wheel. `m113`'s `kn = 2e5` was sized for the pin
path and left the wheels deflecting 55 mm into the band once they moved to the
face — it is now `2.5e6`, giving 4.4 mm under a 1.11e4 N wheel share. `track_min`
and `track_demo`, both authored after the face law, sit at 3.5e6 already, which
is why neither showed the problem.

That also moves the **step-size limit**. As on the ground, it is the friction
regularisation and not `kn` that binds: `μ·f_N/vEps` against the shoe mass. A
tenfold rise in per-pair load is a tenfold rise in that slope, so
`trackContactStabilityStep` now exists beside `groundContactStabilityStep` and
`src/sim/main.c` takes the tighter of the two. Choosing from the ground alone
gave `m113` 6.2e-5 s and the run answered with 10469 single-step energy rises of
~0.05 J while its *net* energy still fell — the signature of a marginally
unstable explicit damper rather than of truncation.

**How the load is assumed to share is the whole content of that estimate, and
it is currently the weak part.** `f_N` is the vehicle weight over the roller
count. Dividing by road wheels only is right for `m113`, whose sprocket and
idler sit above the road wheels and carry no ground load; dividing by *all*
rollers is right for `track_demo`, where all four are tangent to the lower run.
The code does the latter, which halves `track_demo`'s settle (1.4M steps to
698k) and cuts `m113`'s by 28%, at the price of over-estimating `m113`'s step:

| `m113` dt | worst single-step rise | bound | |
| --- | --- | --- | --- |
| 2.02e-5 (road wheels) | 0.00029 J | 0.0088 J | pass, 30x margin |
| 2.83e-5 (all rollers) | 0.00161 J | 0.0094 J | pass, 5.8x margin |
| 6.22e-5 (ground only) | 0.05535 J | 0.0116 J | **fail** |

The rise goes as roughly `dt⁵` there, not the `dt²`–`dt³` of truncation, which
puts the real boundary near 4e-5 — *below* the 5.65e-5 the function returns.
The caller's 0.5 safety factor is what brings it back under, so that factor had
become load-bearing rather than belt-and-braces.

**That is the problem implicit friction solved**, and it is why the table above
is a historical record rather than current behaviour. With the damper solved
implicitly there is no friction stability boundary to estimate, both functions
report the normal spring instead, and 6.2e-5 passes with 48x margin. See
[Ground contact](#ground-contact) for the mechanism and the current sweep. Two
things from this episode still stand, though: the estimate's denominator is
still the roller count, which is right for `track_demo` and wrong for `m113`
(whose sprocket and idler carry no ground load), and the principled fix is still
to count the rollers actually near the lower run — `modelLoad` already has every
roller's centre and radius in the root frame from the hull walk, so "bottom
within a radius of the lowest bottom" is a few lines. It simply stopped being
urgent once the quantity it feeds stopped being the binding constraint.

**Verification, and why not against Chrono.** Chrono's `ChSprocketSinglePin`
runs its own narrowphase callback against a hand-authored tooth profile from
`M113_SprocketSinglePin_Left.json`, colliding against the pin *surface* — not
a pitch-circle abstraction. Reproducing it would mean reimplementing the
profile, at which point the check tests the reimplementation rather than the
physics. What is achievable later is aggregate agreement only: tractive force
against sprocket torque, steady track speed against `ω·R`.

So the checks are self-contained (`src/contact/main.c`, on `m113`):

| check | what it catches | result |
| --- | --- | --- |
| internal wrenches sum to zero over 1000 random states | any sign error, wrong body index, or dropped reaction | 4.3e-16 relative |
| determinism at a repeated state | hidden state creeping into the force law | bitwise |
| radial and circumferential sweeps | a discontinuity that would make an explicit integrator chatter | 1.1 N and 12.6 N per step |
| tooth law fixed points, slope signs, amplitude, periodicity, branch cut | the whole engagement model | exact |
| design pose on the path-defining pulleys | a band jammed into the geometry that defines it | no penetration, worst clearance 0.4 mm |
| sprocket axle torque, assembled from `fext` vs summed directly | a wrong moment reference point — which the balance check *cannot* see, since a wrong reference still cancels between the two bodies | 6.4e-15 relative |
| passivity of the velocity-dependent force | an energy-injecting damper | dissipates everywhere, worst power 0 W |
| face weights partition unity across a joint | support dropped or doubled at every pin, which an explicit integrator answers with chatter | 4.4e-16 on the sum, 3.5e-14 on the slopes |
| face elastic force vs `∇U` by finite difference against an independently written potential | the gated potential's missing tangential term `U₀·w′(x)` — the tooth law's trap in a new place | 1.3e-6 relative over 32004 loaded states |
| face passivity and friction cone | **the friction sign inversion**, which every check above missed because they all drive the *pin* path | dissipates everywhere, cone excess 0 |
| friction impulse bound over 2520 (load, slip, `v_eps`, `dt`, mass) combinations | an implicit-friction sign or reciprocal error, invisible until a trajectory ran long enough to notice | worst step removes `1 − 1.25e-8` of the slip, never all |

The last four postdate the road-wheel face primitive and semi-implicit friction.
The face ones exist because *every* check above them drives `trkPairForce`, which
road wheels stopped using — so all of them stayed clean, `passivity` included,
while the face law injected 5e6 J over a four-second settle. Passivity and the
cone are additionally re-run at `dt = 1e-4`, since the explicit and implicit
modes are different laws and the simulator uses the second one.

And the trajectory checks (`src/sim/main.c`): settle the band from `q0` with
ramped gains, then assert the band was never crushed onto its axles, that the
ground carries the vehicle's weight, that linear momentum tracks
`∫(M·g + F_gnd)dt`, and that the energy budget shrinks once the ramp is over.
On `m113`, at `dt` = 1.87e-4 over 21341 steps: worst pin penetration 0.0087 →
0.0049 m, settled ground load 0.38% off the vehicle weight, momentum matching to
1.2e-5 relative (10.3 N·s against an 8.8e5 N·s scale), creep 0.81 cm/s, and a net
energy change of −55.35 J with the largest single-step rise 0.008983 J against a
0.01354 J truncation bound.

Those are the assertions. The run itself can also just be looked at: `--trace
<file.csv>` writes it out as time histories — hull pose, the ground wrench,
both penetration metrics with their loaded counts, and the world position of
every roller and every shoe — and `tools/plot_trace.py` draws them.

```
build\sim_m113.exe ..\models\m113.mech --trace m113.csv
python tools\plot_trace.py build\m113.csv
```

The trace is off unless asked for, because `build.bat` runs this executable as
a gate and a gate has no business writing files. It is written on a stride
(about 2000 rows whatever `dt` turns out to be; `--stride N` overrides), since
m113 traces 140 bodies over 21341 steps and a full-rate dump is a few hundred
MB to look at four seconds of settling. The file carries its own constants —
weight, roller radii, the terrain plane — in its `#` preamble, so the plots
are drawn from the CSV alone and nothing has to be told which model it came
from.

The energy assertion is deliberately *not* strict monotonicity. A symplectic
scheme's energy error oscillates within a bounded band rather than decreasing,
so an explicit integrator on a stiff penalty law ripples upward occasionally no
matter how correct the forces are; what is asserted is that the net change is a
loss and that no single step adds a truncation-scale fraction of the kinetic
energy. The budget has to be complete for any of that to mean
anything — kinetic, gravitational, the model's own 18 springs, *and* the
contact springs. Omitting the force elements alone made the tensioners' stored
energy look like energy appearing from nowhere.

That last one deserves a note, because it is why **gravity does not need
turning off** even though nothing holds the vehicle up. Gravity is uniform, so
it contributes the same acceleration to every body and cancels out of every
*relative* position and velocity the contact law depends on: the band relaxes
onto its wheels in free fall exactly as it would in orbit, and only the
absolute altitude runs away. It also makes the momentum check exact rather
than approximate — with every other force internal, `dP/dt = M·g` holds
identically, whatever the contact model does. (The check starts after the
ramp, because the artificial viscosity that drains the relaxation transient
multiplies all of `qd`, rigid-body modes included, and is an external force
while it is on.)

That paragraph is now history rather than current state: with ground contact in
the loop, the vehicle *does* stand up, `dP/dt = M·g` no longer holds on its own,
and the momentum check had to be rebuilt around the ground impulse. See below.

## Ground contact

The M113 now stands on the ground under its own weight: 147 of its 504 shoe pads
settle onto flat terrain carrying **1.114e5 N against a 1.111e5 N vehicle
weight — 0.24% off**. That number is the headline check, and it is one an
internal-only band could never have offered.

**Flat ground and uneven ground are the same code path.** The force law reads a
height field `h(x,y)` and its two first derivatives and nothing else, so a
sloped plane, a sinusoid or a sampled grid are new cases in
`src/contact/terrain.h` alone — they change no force law, no wrench
bookkeeping, no energy budget and no test. Flat ground is not a special case at
all; it is `plane 0 0 1 z`.

**The vertical gap, not the perpendicular distance.** This is the one design
decision everything else follows from. Take

```
g  = h(p.x, p.y) - p.z          > 0 = penetrating
F  = ( kn·g + c·ramp(g)·ġ ) · ( -h_x, -h_y, 1 )        clamped at zero
U  = ½·kn·g²
```

with the direction deliberately **not** normalized. Then, because
`∇g = (h_x, h_y, -1)`:

- `-∇U = -kn·g·∇g = kn·g·(-h_x, -h_y, 1)`, which is the elastic force verbatim.
  The law is **exactly conservative on arbitrary terrain**, with no planarity
  assumption and using only *first* derivatives of `h`.
- `ġ = ∇g·u` exactly, so `power = F·u = -kn·g·ġ - c·ramp·ġ²`. The first term is
  `-dU/dt`; the second is a square. The law is **exactly passive**, whatever the
  terrain does, and the clamp can only ever zero it.

The obvious alternative — penetration as perpendicular distance, force along the
unit normal — is exact on a plane and *not the gradient of anything* on curved
terrain. That is the same class of error as the sprocket tooth law's missing
`g′(δ)·U₀` term, which cost two rounds of trajectory-energy debugging to find.
Here it costs nothing to get right, and a check confirms it: finite-differencing
`groundContactPotential` by translating the whole model matches the total elastic
force to **1.4e-9 relative over 153 loaded states**.

There is also **no engagement gate, no polar coordinates, no `atan2` and no
branch cut**. `trkGate` exists because the tooth and lateral terms are gated
potentials; nothing here is gated, and the one-sided `g > 0` test is the entire
activation logic. Ground contact is genuinely simpler than roller contact.

`ground` is a **top-level block**, not a `track` key — there is one ground and
two bands, and a per-band copy stops making sense the moment terrain is a field.
The old `ground <axis> <offset>` had been parsed and read by nothing since the
band landed; `modelLoad` now rejects it rather than ignoring it, telling the
author what to write instead. The shoe's own contact geometry is a `pad` line
giving the plate's shoe-local `z` and footprint, from which the engine takes four
corners per shoe — computed on demand, since materializing m113's 504 of them is
exactly the sort of new model array that overflows the dense validator's stack.

**Damping is derived against the supported mass, not the shoe mass.** Copying
`c = 2ζ√(k·m)` with the shoe's own `m` — correct for a shoe bouncing on a
roller — gives the *vehicle heave* mode ζ ≈ 0.06 and it rings for seconds. What
rests on the ground is the vehicle, so the share is `M_total/n` with `n` the
points expected to be carrying it. Only a band's lower run ever touches, so `n`
is taken as a quarter of the declared pads: the conservative end, since assuming
fewer points both stiffens the damping and shortens the step estimate. The
harness prints the realized ζ after settling so the assumption is checkable.

**The step size used to be set by friction, not by stiffness — and no longer is.**
Regularized Coulomb has tangential slope `μ·f_n/v_eps`, which on m113 is four
points of 4.3e4 N·s/m on an 18 kg shoe: `dt < 2m/c`, an order of magnitude
tighter than the normal damper. Since `dt ∝ v_eps` and `v_eps` is *also* the
creep rate of a parked vehicle, speed could only ever be bought with creep, and
overstepping did not fail loudly — it dripped single-step energy rises into a run
whose net energy still fell.

**Friction is now solved implicitly at each contact.** Writing the slope as
`k = μ·f_n/√(u² + v_eps²)`, both engines apply

```
k_implicit = k / ( 1 + k·dt/m_eff )
```

the backward-Euler solution of `m_eff·u̇ = −k·u` over the caller's own step,
divided back out to a force. The velocity change it implies is
`u·(k·dt/m_eff)/(1 + k·dt/m_eff)`, a factor strictly inside `[0,1)` — friction
can bring the slip to zero but never reverse it, at **any** `dt`. Three things
follow: `dt = 0` reproduces the explicit law bit for bit, so every fixed-state
test keeps testing what it tested before; the divisor exceeds 1, so `|f_t|` stays
strictly inside the cone and no cone bound had to be loosened; and it is smooth
in `u`, unlike the impulse cap `min(k·u, m_eff·u/dt)` which has the same
stability but a kink to chatter on. `src/contact/main.c` asserts all three over
2520 combinations of load, slip, `v_eps`, `dt` and mass.

`m_eff` is the mass the impulse actually accelerates — the shoe for a band
contact, a quarter of it per ground pad corner. **That is not the ground engine's
`mEff`**, which is `M_total/n_loaded`, the *supported* mass that sizes the normal
damper for the vehicle heave mode. The two are deliberately separate fields.

Measured on m113, sweeping `dt` with everything else fixed:

| `dt` | worst single-step rise / bound | equilibrium | creep | |
| --- | --- | --- | --- | --- |
| 2.83e-5 | 2.4e-5 / 1.12e-2 J | 0.36% off | 0.0106 m/s | pass |
| 6.22e-5 | 1.2e-3 / 1.18e-2 J | 0.42% off | 0.0095 m/s | pass |
| **1e-4** | **5.8e-3 / 1.26e-2 J** | **0.40% off** | **0.0088 m/s** | **pass** |
| 2e-4 | 1.5e-2 / 1.36e-2 J | 0.38% off | 0.0081 m/s | fail |
| 5e-4 | 1.9e-1 / 1.30e-2 J | 0.33% off | 0.0075 m/s | fail |
| 1e-3 | 2.2 / 1.21e-2 J | 0.29% off | 0.0061 m/s | fail |

At 6.22e-5 — the step that *failed* under explicit friction with a 0.05535 J
rise — the rise is now 1.157e-3 J, **48x smaller**. Nothing diverges even at
1e-3, 35x the old step: net energy still falls and the vehicle still stands to
0.29%. The failures from 2e-4 are the truncation gate (bound = 0.01% of the KE
scale), which is the *accuracy* limit reporting itself, and the right response is
a smaller step rather than a looser bound.

**Two things worth not glossing.** First, creep goes *down* with `dt`, which is
the opposite of the naive prediction that the effective regularization velocity
becomes `v_eps + μ·f_n·dt/m_eff`. That estimate assumes saturated friction; near
equilibrium the residual driving force is small, and more to the point the
implicit impulse tends to `m_eff·u` as `dt` grows, removing essentially the whole
slip each step. Large `dt` is therefore *stickier* than reality — ordinary
backward-Euler over-damping — so the low creep at 1e-3 is a numerical artifact,
not a better answer. Second, none of this makes true stick reachable: the law is
still stateless and still creeps, which is still the price of being testable at
random states.

Net effect on m113: `dt` goes 2.83e-5 → 1e-4 at a step cost of 175.4 µs, so the
real-time factor goes **0.16 → 0.57** — from 6.3x slower than real time to 1.75x.
The next limit is truncation, not stability, and that is what the section below
goes after.

**What the harnesses had to change, and what got stronger.** The ground is not a
body, so ground reactions are external and two long-standing checks stopped
meaning what they said:

| check | before | now |
| --- | --- | --- |
| Newton's third law | all contact internal, wrenches sum to zero | replaced: the wrench must match an independently written sum, **and land only on padded bodies** — with no reaction body, a wrench at the wrong index leaves no trace in any total |
| linear momentum | `ΔP = M·g·Δt` identically | `ΔP = ∫(M·g + F_gnd)dt`, measured against the wrench the engine reports applying |
| band interpenetration | must go down | a band that starts *tangent* has to compress to carry load; now "never crushed onto its axles", plus "if it started jammed it must have relaxed" |
| — | — | **new: static equilibrium.** Settled ground load = vehicle weight |

The momentum residual is asserted to be truncation rather than assumed to be:
halving `dt` twice from 6.22e-5 takes it **33.4 → 4.91 → 1.64 N·s** against a
scale that stays put at 8.25e5. A dropped reaction does not converge — it is
order 1e5. (Pass an explicit `dt` to `sim_<model>.exe` to repeat that.) The
normalization had to change too: dividing by the momentum present was right in
free fall and meaningless once the vehicle is standing still, where it measures
how quiet the settle happened to be rather than whether a force went missing.

Five new checks in `src/contact/main.c`, all inert on models with no `ground`
block, on top of the two above:

| check | result on m113 |
| --- | --- |
| wrench vs an independent sum; lands only on padded bodies; deterministic | 1.8e-15 relative |
| vertical sweep through the surface | 0.70 N step against a 2.24 N bound |
| horizontal sweep across terrain | 0 N over 20 m |
| **elastic force is exactly `-∇U`**, by finite-differencing the potential | 1.4e-9 relative, 153 loaded states |
| terrain gradient vs central differences of its own height | exact |
| passivity, friction cone `\|f_t\| ≤ μ·f_n`, and `f_t ⊥ n̂` | exact |

The horizontal sweep is aimed at terrain that does not exist yet. On a plane it
is trivially zero, which is the point of writing it now: when a sampled grid
arrives, bilinear interpolation gives C0 heights and therefore a normal
direction that snaps at every cell edge, and this is the check that will reject
it and force a bicubic.

**Not validated against Chrono, and not for the sprocket's reason.** Flat rigid
ground *is* reproducible there, but Chrono's contact is its own penalty model
with its own parameters, so accelerations would not agree to 1e-9 and forcing
them to would be tuning one model against another. The analytic gates are
strictly stronger and need no Chrono at all — which is also what the next stage
is built around: slope holding below `tan θ = μ`, slope release above it at
`g(sin θ - μ cos θ)`, sinkage `W/(n·kn)`, and tractive force `T/R` against
sprocket torque.

**What is deliberately not here.** Triangle meshes (C⁻¹ normals at edges, and
they need a broadphase); unilateral constraints, for the same frozen-`nc` reason
penalty contact was chosen for the band; and stateful stick friction, which would
forfeit the random-state testability the whole harness rests on. Sphere-against-
terrain — the honest replacement for `hmmwv`'s four shaker posts and `distance`
"tire" rods — reintroduces exactly the gradient subtlety the vertical-gap point
formulation avoids, so it waits.

The fidelity gap that stood here — the belt wrapping only the two `end` pulleys,
leaving m113's road wheels 0.22 m through the lower run for the settle to relax —
is closed. `trkBuildBelt` gift-wraps the convex hull of every roller circle, and
the road wheels roll on the shoe's flat inner face rather than on its pins; the
design pose now starts 0.0087 m from equilibrium and settles to 0.0049 m.

## Second-order stepping

Once friction was solved implicitly, nothing about *stability* bound the step any
more — both contact engines reported limits of 1e-3 to 5e-3 s, one to two orders
looser than the step in use. What bound it was **accuracy**: first-order
truncation, measured against a gate set at 0.01% of the kinetic energy scale. The
remaining 1.75x to real time was therefore an integrator problem, not a contact
problem.

`src/sim/integrate.h` now carries two schemes. Semi-implicit Euler is unchanged
and is kept as the control. The new one is **Störmer–Verlet**, composed by the
caller as drift / kick / drift:

```
q  = drift( q, qd, dt/2 )
qdd = f( q, qd_mid )                  <- ONE dynamics call, at the midpoint
qd += qdd * dt
q  = drift( q, qd, dt/2 )
```

**One force evaluation per step is the whole point.** `fwd_dyn_crba` is 82% of
the step budget, so a scheme needing two evaluations must more than double `dt`
just to break even. This needs one; the extra cost is a second `simIntegrateQ`,
measured by `bench_<model>.exe` at **292 ns of a 175 µs step — 0.2%**. The force
lands at the midpoint, which is why this form and not the algebraically
equivalent kick/drift/kick: every damper here is velocity-dependent, and
kick/drift/kick wants the force at the end state whose velocity is the unknown.

Order is asserted, not claimed. `energyDrift()` runs **both** schemes at `dt` and
`dt/2`, and each is gated on its own ratio — Euler ~2, Verlet ~4. Eight tree
models report 2.000 and 4.000 to three digits. Running both is what makes it a
measurement: if Verlet were quietly still first order, the two would agree, and
simply retargeting one band from 2 to 4 could not tell that apart from the real
thing. The absolute drift is the more striking half — `satellite` at `dt` = 1e-4
goes from 2.63e-1 J to 2.22e-6 J.

That accuracy also broke the gate's own noise floor, which was a flat 1e-12 J.
`gantry` lands at 4e-12 J and reports *more* drift at `dt/2` than at `dt` — the
signature of round-off accumulating over twice as many steps, not a scheme that
got worse. The floor is now relative, 1e-11 of the KE scale.

### The predictor, and why the track loop does not use one

Verlet is second order only if the force is evaluated at the midpoint in **both**
slots. The position comes free from the half-drift. The velocity needs
`qd + qdd_prev·dt/2`, and on the tree models that predictor is exactly what earns
the order: drop it and all eight fall from 4.000 to 2.000.

On the contact models it does the opposite, by a lot. A damper `-c·v` evaluated
at `v + δ` but applied to `v` gives `Δv = -(c·dt/m)(v + δ)`; the `δ` part is not
damping at all, it is forcing, and it does not care which way `v` points. The
implicit friction law's guarantee — an impulse that can zero the slip but never
reverse it — is stated about the velocity the impulse is *applied* to, so
predicting breaks it outright. Near equilibrium the true slip is on the order of
`v_eps` while `qdd_prev·dt/2` is set by the contact stiffness, so the predicted
slip can arrive with the wrong sign. Worst single-step energy rise on m113,
against a bound of ~0.013 J:

| `dt` | no predictor | predictor on both | predictor on `fwd_dyn_crba` only |
| --- | --- | --- | --- |
| 5e-5 | 0.000513 | 0.006219 | — |
| 1e-4 | 0.002437 | 0.034680 | 0.003175 |
| 2e-4 | 0.010160 | 0.270600 | 0.019170 |

Predicting got *relatively worse* as `dt` shrank, which is injection and not
truncation — a higher-order scheme cannot be uniformly worse and lose ground as
the step falls. The third column kills the tempting split, too: predict for the
dynamics, whose velocity dependence is conservative Coriolis, and not for
contact, whose velocity dependence is dissipative. It does not work because
m113's own ten torsion bars and six shocks are dampers as well, and they live
*inside* `fwd_dyn_crba`.

So the track loop keeps the second-order position and pays first order in the
velocity slot. That still beats semi-implicit Euler by **2.4x** at `dt` = 1e-4
and passes at 2e-4 where Euler does not. Recovering the full order here means
treating the contact damper and the model's own dampers implicitly, the way
friction already is — a separate change, and now one with a measured reason.

### The step is chosen by accuracy now

`src/sim/main.c` used to cap the step at a flat `1e-4` and let the stability
estimates lower it. With stability no longer binding, that constant *was* the
step for all three track models — and it only fit one of them. The cap is now a
fixed number of samples per contact **period**, `2π·sqrt(m/kn)/N`, from
`trackContactAccuracyStep` and its ground twin.

| model | contact period | largest passing `dt` | samples | `dt` chosen at N = 90 |
| --- | --- | --- | --- | --- |
| `m113` | 16.9 ms | 2.0e-4 (fails at 2.5e-4) | 84 | 1.87e-4 |
| `track_min` | 9.5 ms | 1.2e-4 (fails at 1.5e-4) | 79 | 1.06e-4 |
| `track_demo` | 7.5 ms | ≥ 2e-4 (never failed) | — | 8.34e-5 |

`track_demo` is the one the cap costs something. It settles at 2e-4 without a
single upward energy step and still gets 8.34e-5, because its band has the
highest natural frequency per shoe of the three and the cap only knows about
that. Conservative in the safe direction, and `argv[2]` is there for anyone who
wants to spend the margin.

Two models landing on the same number is what makes a period fraction the right
shape and one magic constant the wrong one. `track_min` has a single road wheel
carrying what m113 spreads over five, so it is not a scaled-down m113 — and
under semi-implicit Euler that never showed, because 1e-4 was inside every
model's limit. `N` = 90 is the default rather than the 77–84 the sweep marks as
the boundary: at 90, m113 uses 66% of its single-step energy budget, about the
margin Euler had at 1e-4, where 80 would use 83% and sit a few percent from the
cliff.

### Result

| | step | cost | real-time factor |
| --- | --- | --- | --- |
| explicit friction, Euler | 2.83e-5 | 178.6 µs | 0.16 (6.3x slower) |
| implicit friction, Euler | 1e-4 | 175.4 µs | 0.57 (1.75x slower) |
| **implicit friction, Verlet** | **1.87e-4** | **174.9 µs** | **1.07 (real time)** |

Those factors are `dt` divided by the benchmarked step cost, which is what the
earlier rows measured too, so the column is comparable down its length. It is
**not** the same as `sim_m113.exe`'s wall clock, and the difference is worth
stating rather than rounding away: the harness settles 4.00 s in 4.21–4.34 s
over five runs, **0.94x**, because it also pays for its own gates. The energy
budget runs on the last 35% of steps and costs an extra `denseKinematics`
(13.9 µs) plus `simEnergyK`, which calls `denseKinematics` a second time, plus
both contact potentials — about 0.23 s of that 4.24 s. Subtract just that and
the harness lands at 1.00x, where the step-cost figure says it should.

So m113 — 153 bodies, 882 shoe-roller pairs, 504 ground points, contact in the
loop — runs at real time within measurement noise: a bare step loop clears it,
the verification harness sits ~6% short. Run-to-run spread is ~3% on the wall
clock and ~4% on the benchmark, so "about 1.0" is as tight a claim as the
numbers support. `track_demo` runs at 2.3x real time and `track_min` at 10.2x.

The gates all of this passes are the ones that were already there, at the same
tolerances; none was loosened to get here.

## Force elements and points of interest

Built around two principles the project already followed: **declare what you
pay for** (nothing appears in the generated code that the `.mech` file didn't
ask for — the `fext 1` contract, extended) and **the dynamics API stays the
stable contract** while new machinery grows alongside it. `springpend.mech`
exercises every element code path in a loop-free model; `slidercrank.mech`
combines a passive element with a kinematic loop; `quartercar.mech` carries
the real HMMWV front coil spring and shock as two `tsda` elements (linearized
from the shipped nonlinear spring curve — this format's TSDA is linear-only —
at zero deflection; the shock's damping coefficient was already linear).

The design:

- **Points of interest (POIs)** — named points fixed to a body, declared in
  the `.mech` file (`poi <name> <body> <x y z>`); the generator emits their
  world kinematics. A TSDA anchor *is* a POI, so anchors and observation
  points share one declaration and one symbolic treatment.
- **Force elements are first-class** — `tsda <name>` blocks (translational
  spring-damper-actuator between two named POIs) and `rsda <name>` blocks
  (rotational — or, on a prismatic joint, translational; the name is carried
  over from vehicle-dynamics terminology — acting on a named tree joint's own
  coordinate). A TSDA is generated symbolically by the same machinery as the
  `distance` loop joint: the applied generalized force `Q = -(f/l) * row`
  reuses the path-walk row (`gen_elems.h`), so it touches only the joints on
  the path between the two anchors — no 6-vector spatial-force detour. An
  RSDA on a tree joint is a torque straight into that joint's own coordinate,
  `tau_e = -k*(q - q0) - c*qd + act[slot]` (RSDA across a loop-closing
  revolute is deferred; it would need an atan2 relative angle).
- **`fext` unchanged** — stays opt-in per body for genuinely external loads
  (tire contact, validation probes); elements never route through it.
- **Element inputs are opt-in** — every element carries a baked linear law
  (`k`, `c`, and a reference `l0`/`q0`) from the `.mech` file; elements
  additionally declared `input` get a slot in a compact `act[]` array, and
  the total force is `baked law + act[slot]`. Passive elements compile to
  pure expressions with no runtime input surface. The dynamics signature is
  `fwd_dyn_crba( q, qd, tau, fext, act, out_qdd )`, with `act` NULL-able
  exactly like `fext` (and `fwd_dyn_crba_nact` exported alongside `nq`/`nv`).
  Caller-side nonlinear laws (real damper curves) use the pattern
  `k 0 c 0 input`: compute the force outside from the element's observed
  state, feed it through `act`.
- **Name lookup via runtime tables** — the generated file exports
  `fwd_dyn_find( kind, name )` (`kind` one of `"poi"`, `"tsda"`, `"rsda"`,
  `"act"`) backed by const string tables (`gen_observe.h`), resolved once at
  caller init. Runtime lookup rather than a generated header because the
  eventual JIT consumer is a fixed host program that cannot recompile per
  model. Returns -1 for an unresolved kind or name.
- **Observation function** — generated `fwd_dyn_observe( q, qd, out_obs )`
  exposes POI world positions and velocities plus each element's state
  (`l`, `ldot` for TSDAs; angle, rate for RSDAs), indexed through the same
  lookup, at the layout `fwd_dyn_observe_len` gives the size of: all POIs (6
  doubles each), then all TSDAs (2 each), then all RSDAs (2 each). It runs
  the same per-body kinematics/velocity emission as the dynamics pass
  (`gen_kinematics.h`) in its own function scope, so it duplicates no
  dynamics (CRBA `H`, RNEA bias force) — and is skipped entirely, at zero
  runtime cost beyond the two exported constants, for a model with nothing
  declared.

Validation follows the same two-routes discipline as everything else in this
project: `src/validate/dense.h` applies each TSDA as **explicit
equal-and-opposite spatial forces at its two anchors in base coordinates** —
deliberately not the generator's path-walk row — so agreement is a real check
on that row's derivation (RBDA eq 8.18–8.20's epsilon structure), not a
shared-bug false positive; RSDA has no second formulation (there is no
alternative way to express "a torque on your own joint"), so both sides use
the same law there. `fwd_dyn_observe`'s TSDA `l`/`ldot` are cross-checked
against an independent recomputation from the same `Xa` transforms the
element-force block uses, and `fwd_dyn_find` is checked to round-trip every
declared name and reject a bogus one. All of it validates to 1e-9 on
`build.bat` and to ~1e-6 against Project Chrono on `build_chrono.bat` (see
`src/chrono` above for how — and how not — force elements map onto Chrono's
own link types).

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

`simIntegrateQ()` in `src/sim/integrate.h` writes the floating case **exactly**
rather than to first order, because a Verlet composition needs its half-steps to
be self-adjoint (see *Second-order stepping* above):

```
Q(tau) = Q0 x exp( 0.5 * w * tau )      exact for constant body-frame w
r(tau) = r0 + R( Q(tau/2) ) * v * tau   midpoint rule, self-adjoint
```

The exponential is a unit quaternion by construction, so the renormalisation the
truncated update needed is gone. Reversibility is the cheap way to see the
difference: drift by `+tau` then `-tau` and the exact map returns to the start to
1e-10, while the first-order one leaves a position residual of 9.4e-3 — and it
was the *position* half that was at fault, since the truncated quaternion update
is already self-adjoint on its own.

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

For every model in its list (all eighteen above), this builds
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

The generated file exports the dynamics function, the observation function,
name lookup, and their size constants:

```c
void fwd_dyn_crba( double *q, double *qd, double *tau, double *fext, double *act, double *out_qdd );
extern const int fwd_dyn_crba_nq;    // length of q
extern const int fwd_dyn_crba_nv;    // length of qd, tau, out_qdd
extern const int fwd_dyn_crba_nact;  // length of act (force-element input slots)

void fwd_dyn_observe( double *q, double *qd, double *out_obs );
extern const int fwd_dyn_observe_len; // length of out_obs

int fwd_dyn_find( const char *kind, const char *name ); // "poi"|"tsda"|"rsda"|"act"; -1 if not found
```

`fext` is `6 * nbodies` doubles (base-coordinate spatial force per body; see
[Conventions](#conventions)). `act` is `fwd_dyn_crba_nact` doubles, one per
force element declared `input` (in declaration order); NULL-able exactly like
`fext` when no element has one. `out_obs` is `fwd_dyn_observe_len` doubles;
see [Force elements and points of
interest](#force-elements-and-points-of-interest) for the layout and
`fwd_dyn_find`'s lookup. To use a mechanism of your own, drop a `.mech` file
in `models\` and add its name to the `MODELS` list at the top of `build.bat`
— it will be generated, validated, and benchmarked like the others.

### `validate_<model>.exe` — dense reference vs generated code

```
build\validate_threelink.exe models\threelink.mech
```

Runs 1000 random states (near-manifold states for loop models, random `act`
alongside `q`/`qd`/`tau`/`fext` when the model has force elements) through
both the compiled-in generated code and the dense 6x6 reference and compares
qdd to 1e-9, also checking rank(K) = nc for loop models at every trial and,
for models with elements, `fwd_dyn_observe`'s TSDA `l`/`ldot` against an
independent dense recomputation. `fwd_dyn_find` is checked once, up front, to
round-trip every declared name and reject a bogus one. Exit code 0 = PASS,
1 = FAIL:

```
validate: models\threelink.mech: 1000 trials, max scaled |qdd_generated - qdd_dense| = 1.699e-14 (worst at trial 969)
validate: PASS
```

```
validate: models\springpend.mech: 1000 trials, max scaled |qdd_generated - qdd_dense| = 8.105e-15 (worst at trial 888)
validate: models\springpend.mech: 1000 trials, max scaled |observe_generated - observe_dense| = 1.110e-16 (worst at trial 16)
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
acceleration cross-check (staging in random `act` alongside `q`/`qd`/`tau`/
`fext` for models with force elements), and the default simulation
comparison (T = 0.25 s, dt = 2e-5, repeated at dt/2 to show the gap is
integrator truncation; `act` held at zero there, so both sides agree by
construction):

```
chrono_check: models\threelink.mech: qd-vs-dq/dt consistency, max scaled error = 5.799e-10 (worst at trial 11)
chrono_check: models\threelink.mech: 200 trials, max scaled |qdd_chrono - qdd_dense| = 9.261e-10 (worst at trial 192)
chrono_check: models\threelink.mech: 200 trials, max scaled |qdd_chrono - qdd_generated| = 9.261e-10
chrono_check: PASS
chrono_check: models\threelink.mech: simulation loop, generated dynamics vs Chrono, T=0.25 dt=2e-05
    t= 0.031  gen b2( 2.8723  0.7004  0.0000)  chrono b2( 2.8723  0.7004 -0.0000)  max|gap|=1.07e-08 at b0
    ...
    t= 0.250  gen b2( 2.9017  0.6093  0.0000)  chrono b2( 2.9017  0.6093 -0.0000)  max|gap|=8.56e-08 at b0
chrono_check: max |gap| = 8.561e-08 at dt, 4.286e-08 at dt/2  (position in m or attitude in rad, whichever is worse)  (ratio 2.00 -- shrinks with the step, so the gap is integrator truncation)
chrono_check: SIM PASS (tolerance 1.0e-02)
```

`--sim` skips the first two checks and runs just the simulation comparison
at a chosen horizon and step, printing the trajectory table of both sides —
useful for watching a single model evolve:

```
build\chrono_check_elbow.exe models\elbow.mech --sim 0.25 2e-5
```