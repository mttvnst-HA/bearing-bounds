# CLAUDE.md: Bearing Bounds

Single-page browser app for lower- and upper-bound bearing capacity of shallow foundations by plane-strain limit analysis. See README.md for features and docs/method.md for the mechanics and history.

## Who this is for and how to work

a. The owner is a geotechnical project manager (US units by default). Treat him as the project manager and yourself as the engineer: deliver complete, working changes, not snippets.

b. He works on Windows 11. Keep scripts OS-agnostic: Python with `pathlib`, Node with `path`; no bash-only constructs in npm scripts.

c. Define acronyms on first use in any user-facing text (UI copy, docs, messages).

d. Do not use bullet lists in docs or replies; use lettered or numbered lists.

e. When reviewing a document, cite a short searchable text snippet, never a paragraph number.

## Commands

```bash
npm install              # HiGHS solver + Playwright
npm run build            # python scripts/build.py -> dist/
npm run verify           # solver vs known solutions (Quick mesh, a few minutes); add -- --no-cone for HiGHS only
npm run verify:std       # Standard mesh, slower
npm run e2e              # Playwright test of dist/bearing-bounds-standalone.html (needs `npx playwright install chromium` once)
npm test                 # verify + build + e2e
python scripts/build_clarabel.py   # rebuild vendor/clarabel from wasm/clarabel (needs Docker; rarely needed)
```

Always run `npm run verify` after touching `src/engine.js`, and `npm run build && npm run e2e` after touching `src/app.html`. A change that breaks "LB <= reference <= UB" in verify is a correctness bug, not a tolerance to relax.

## Architecture

a. `src/engine.js` is a UMD-style module (`module.exports` in Node, `globalThis.FELA` in the browser). Main entry: `analyse(highs, model, opts)`. Model inputs are SI (m, kPa, kN/m³); `gammaW` = 9.81.

b. Coordinates: y up, footing base at y = 0, footing centre at x = 0, tension positive. Depths in the model input are measured below the top ground surface.

c. A symmetric half model is used only for a vertical, central load on level ground; otherwise the full section is solved.

d. Lower bound (`lowerBound`): linear stress triangles, edge discontinuities, inscribed yield surface, and an extension of the stress field beside and below the mesh. That extension is what makes the bound rigorous; never remove it or fit the domain tightly without it. The yield surface is a set of "cones" (`conesAt`: Mohr-Coulomb, tension cut-off, anisotropic undrained ellipse), each linearised by a `p`-sided polygon or, when `opts.nu` is set, the lifted Ben-Tal–Nemirovski polygon (`liftedCone`). Lower bounds must use the inscribed form (scale `cp`), upper bounds the circumscribed form or its LP dual (`addFlowUB`).

e. Upper bounds: `upperBoundDLO` (discontinuity layout optimisation, rigid footing with a no-tension interface, `columnWeight` for self-weight work), `upperBound` (linear finite element with jumps) and `upperBoundQ` (6-node quadratic, continuous; `opts.feQuad`). `dloAdaptive` grows the domain until the mechanism is clear of the edges.

f2. Loads are F = Ffix + alpha Fref with F = (V, H, M) (`footingFrame`, `footingPower`); every solver maximises or minimises alpha. `res.best` reports alpha and the collapse V, H, M. `analyse` runs `opts.adaptCycles` mesh cycles (indicators from the lower-bound row duals and upper-bound plastic strain) and keeps the best bound of each kind. `strengthFS` brackets the factor of safety on strength; each end must stay on its safe side.

f3. Axisymmetric (circular footing, `model.axisym`): `axiLowerBound` (quadratic r-weighted stresses, Bernstein yield checks, full 3D Mohr-Coulomb via `axiConesAt`) and `axiUpperBound` (quadratic velocities with quadratic jumps, Bernstein multipliers). The plane-strain DLO only sizes the domain and is never reported as a bound for the circle.

f. LP solves go through `runLP`. A bound is accepted only when HiGHS reports a feasible primal point (`primal_solution_status` = 2); otherwise it re-solves with crossover, then dual simplex. Do not report bounds from unconverged solves.

g. The finite element solvers (`lowerBound`, `upperBound`, `upperBoundQ`, `axiLowerBound`, `axiUpperBound`) are wrapped by `conicFirst`: they are first built with exact second-order cones (`opts._soc`, `lp.addSOC`) and solved by Clarabel (`runConic`, module registered with `setConic`), then rebuilt linearised for HiGHS if that fails. `runConic` accepts a point only after checking every row, bound and exact cone to 1e-6; never loosen that check or skip it. Stress points at the cone apex (traction-free, cohesionless) are fixed to zero in cone mode.

h. Meshing: `buildTriMesh` is a Bowyer-Watson Delaunay mesh with polar fans at the footing edges. It is triangulated band by band so layer boundaries, the water table, the slope face and toe are element edges. `buildQuadMesh` is a fallback for level ground.

i. Effective-stress basis throughout; pore pressure at the base is added back to report gross total pressure.

## UI (`src/app.html`)

a. The build replaces five markers: `/*__HIGHS__*/`, `/*__ENGINE__*/`, `__WASM_B64__`, `/*__CLARABEL__*/`, `__CLARABEL_B64__`. Keep them. Clarabel comes from `vendor/clarabel` (see its README for version and licence).

b. The solver runs in a Web Worker built from a Blob of the inlined HiGHS, Clarabel and engine sources. The worker is discarded after every run to release WebAssembly memory. Never run the solver on the main thread (it froze and crashed Claude Desktop).

c. State `S` stores every input in SI. `toSI`/`fromSI` convert at the field boundary. Strip loads are per unit length; other footings take totals. Switching footing type converts loads so the applied pressure stays the same.

d. Results are tied to `result.inputs` and `result.key` (from `runKey`). Changes to anything except the size of the load and the target factor of safety mark results stale.

e. The page is published as a Claude artifact. The artifact host wraps it in a document skeleton, so `src/app.html` has no `<html>`, `<head>` or `<body>` of its own. External scripts are not allowed except from approved CDNs, so everything is inlined. `dist/bearing-bounds.html` is the artifact file; the standalone file is for local use and tests.

f. Theme: colour tokens on `:root` with dark-mode overrides; keep every colour a token.

## Known limits and backlog

See the "Known limits" and "Backlog" sections of docs/method.md before starting new work.
