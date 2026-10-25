# Bearing Bounds

Bearing Bounds computes the collapse load of a shallow foundation as a **lower bound** and an **upper bound** from plane-strain limit analysis. It handles layered Mohr-Coulomb soil, inclined and eccentric loads, and footings near a slope. Everything runs in one HTML page in the browser; the linear programs are solved locally by the HiGHS solver compiled to WebAssembly.

It was built as an OPTUM-style finite element limit analysis (FELA) tool focused on foundations.

## What it does

a. **Foundations:** strip, square, rectangular, circular, mat and slab footings, with embedment and a rough or smooth base.

b. **Loads:** vertical load, horizontal load and eccentricity. All three are increased together until collapse.

c. **Ground and soil:** level ground or a slope on one side of the footing; up to six soil layers, each drained (c′, φ′) or undrained (sᵤ, optionally increasing with depth); and a water table.

d. **Results:** factor of safety on the lower bound, allowable and applied pressure, the collapse pressure range, and a classical Vesić/Meyerhof comparison with a reason for any difference.

e. **Views:** Model, Collapse animation, Mechanism (upper bound), Plastic zones + velocity, Mesh, and Dissipation (finite element upper bound).

f. **Units:** US (ft, psf, kips) or SI (m, kPa, kN). Inputs are stored in SI and converted only for display.

## Methods in one paragraph

The **lower bound** is FELA with linear stress triangles, stress discontinuities on every edge and an inscribed linearised Mohr-Coulomb yield surface. The stress field is extended to the infinite half-space, so the result is a rigorous lower bound at any mesh size. The **upper bound** comes from discontinuity layout optimisation (DLO). For layered soil, slopes and combined loads, a finite element upper bound with velocity discontinuities also runs, and the lower of the two is reported. The footing is a rigid body on a no-tension interface. Results are gross total pressures at the footing base. Non-strip shapes apply a classical shape-factor ratio to the plane-strain bounds, which is an estimate, not a rigorous bound. See [docs/method.md](docs/method.md) for detail, verification and known limits.

## Getting started (Windows, macOS or Linux)

You need [Node.js](https://nodejs.org) 18 or later and Python 3.9 or later.

```bash
npm install                      # gets the HiGHS solver (and Playwright for tests)
npm run build                    # writes dist/bearing-bounds.html and dist/bearing-bounds-standalone.html
```

Open `dist/bearing-bounds-standalone.html` in Chrome or Edge to use the app. `dist/bearing-bounds.html` is the version without a document skeleton, for publishing as a Claude artifact.

## Tests

```bash
npm run verify                   # solver against known solutions, Quick mesh (a few minutes)
npm run verify:std               # same with the Standard mesh (slower, tighter bounds)
npx playwright install chromium  # first time only
npm run e2e                      # browser test of the built page
npm test                         # all of the above, Quick mesh
```

`verify` checks that the bounds bracket exact solutions: Prandtl's Nc, Nq, Martin's Nγ, Green's inclined load, a crust over soft clay, an eccentric load and a footing at a slope crest.

## Project layout

| Path | What it is |
|---|---|
| `src/engine.js` | Solver: geometry, meshing, lower bound, upper bounds (DLO and finite element), driver. Runs in Node and in a browser worker. |
| `src/app.html` | User interface template (HTML, CSS and the page script). The build inlines the solver into it. |
| `scripts/build.py` | Build script; inlines HiGHS (JS and base64 WebAssembly) and the engine into one HTML file. |
| `tests/verify.js` | Solver verification against known solutions. |
| `tests/e2e.js` | Playwright browser test of the built page. |
| `docs/method.md` | Method, verification table, review history and known limits. |
| `CLAUDE.md` | Working notes for Claude Code sessions on this repo. |

## Status

This is an engineering prototype. It has had an independent solver review and an independent UX review, and the findings were fixed (see `docs/method.md`). Its changes since then for inclined loads, eccentric loads and slopes have been tested here but not independently reviewed. Check results against established methods before using them for design.
