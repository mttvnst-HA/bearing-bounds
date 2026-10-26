# Clarabel in WebAssembly (vendored)

The interior-point conic solver Clarabel, compiled to WebAssembly, for the exact second-order-cone (SOC) finite element bounds. `scripts/build.py` inlines both files into the page, as it does for HiGHS; the tests load them with `tests/clarabel.js`.

1. `clarabel_bb_bg.wasm`: the module (about 300 kB; about 400 kB as base64 in the page).
2. `clarabel_bb.js`: the wasm-bindgen glue, target "no-modules" (defines the global `wasm_bindgen`; start it with `wasm_bindgen.initSync({ module: bytes })`).

## Source and versions

a. Built from `wasm/clarabel` in this repository: a minimal binding written for Bearing Bounds (not the cvxjs binding, which discards the solution unless Clarabel reports "Solved"). It returns the status, iteration count, primal x, dual z and slack s for every status, including "AlmostSolved"; the engine decides what to accept after its own admissibility check.

b. Clarabel.rs 0.11.1 (crates.io, pinned with `=0.11.1`), by Paul Goulart, Yuwen Chen and contributors, University of Oxford. Licence: Apache License 2.0, copied in `LICENSE-clarabel.md`. Source: https://github.com/oxfordcontrol/Clarabel.rs.

c. Other compiled-in crates (exact versions in `wasm/clarabel/Cargo.lock`): amd 0.2.2 (BSD-3-Clause, approximate minimum degree ordering), wasm-bindgen 0.2.100 and js-sys 0.3.77 (MIT or Apache-2.0), serde and serde_json (MIT or Apache-2.0), web-time 0.2.4 (MIT or Apache-2.0).

d. Toolchain: Docker image `rust:1.86-slim`, target `wasm32-unknown-unknown`, release profile with link-time optimisation, then `wasm-bindgen-cli` 0.2.100.

## Rebuilding

Only needed after a change to `wasm/clarabel`. With Docker Desktop running:

```bash
python scripts/build_clarabel.py
```

## Interface

`solve(n, m, q, A_colptr, A_rowidx, A_val, b, n_zero, n_nonneg, soc_dims, settings)` solves: minimise q·x subject to A x + s = b, with s in (zero cone of size n_zero) × (non-negative cone of size n_nonneg) × (second-order cones of the given sizes). A is in compressed sparse column form. It returns a Float64Array `[status, iterations, solve time (s), primal objective, dual objective, x (n), z (m), s (m)]`, or `[-1]` if Clarabel rejects the data. The status codes and the order of `settings` are listed in `wasm/clarabel/src/lib.rs`.
