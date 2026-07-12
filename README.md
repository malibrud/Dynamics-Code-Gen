# dynamics

An implementation of Roy Featherstone's spatial-vector rigid-body dynamics
algorithms (RBDA), aimed at being fast enough for real-time use. Featherstone's
formulation uses 6x6 spatial transforms, but for single-DOF joints most of
that matrix is structurally zero or one — this project generates C code that
exploits that sparsity directly instead of doing dense 6x6 math at runtime.

The eventual goal is to JIT-compile the generated C with LLVM and call it
directly, rather than shelling out to a C compiler.

## Layout

- `src/codegen` — the active approach. A small C++ program that walks the
  mechanism definition and, instead of computing with numbers, computes with
  strings: each spatial-vector/transform component is a symbolic expression.
  `compose()`, `multFactors()`, `addTerms()` etc. simplify away structural
  zeros and ones as they go, then `emit()` prints each simplified expression
  as a named C temporary and rewires later composition to reference that name
  (so expressions never re-expand across pipeline stages). Running the
  generator produces `build/dynamics.c` — the specialized, sparsity-exploiting
  C implementation for the mechanism it was given.
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

This builds `build/codegen.exe` from `src/codegen`, runs it to emit
`build/dynamics.c`, and attempts to build `build/numeric.exe` from
`src/numeric` (currently expected to fail — see above).
