# Bearing Bounds

Bearing Bounds computes the collapse load of a shallow foundation as a **lower bound** and an **upper bound** from limit analysis: plane strain, or axisymmetric for a circular footing. It handles layered Mohr-Coulomb soil, inclined and eccentric loads, and footings near a slope. Everything runs in one HTML page in the browser; the linear programs are solved locally by the HiGHS solver compiled to WebAssembly.

It was built as an OPTUM-style finite element limit analysis (FELA) tool focused on foundations.

## What it does

a. **Foundations:** strip, square, rectangular, circular, mat and slab footings, with embedment and a rough, smooth or partly rough base (interface factor R). A circular footing under a central vertical load on level ground is analysed axisymmetrically, which gives rigorous bounds for the circle itself.

b. **Loads:** vertical load, horizontal load and eccentricity, increased together until collapse; or V increased with H held (wind, seismic); or V held with H increased (sliding capacity). Optionally, the factor of safety on soil strength (strength reduction), for the applied load or for a slope alone.

c. **Ground and soil:** level ground or a slope on one side of the footing; up to six soil layers, each drained (c′, φ′, optional dilation angle ψ′ through Davis's reduced strength) or undrained (sᵤ, optionally increasing with depth and optionally anisotropic); an optional tension cut-off per layer; and a water table.

d. **Results:** factor of safety on the lower bound, allowable and applied pressure, the collapse pressure range with the mean of the bounds (an estimate, labelled as such), the strength factor of safety when asked for, and a classical Vesić/Meyerhof comparison with a reason for any difference.

e. **Views:** Model, Collapse animation, Mechanism (upper bound), Plastic zones + velocity, Mesh, and Dissipation (finite element upper bound).

f. **Units:** US (ft, psf, kips) or SI (m, kPa, kN). Inputs are stored in SI and converted only for display.

## Methods in one paragraph

The **lower bound** is FELA with linear stress triangles, stress discontinuities on every edge and an inscribed linearised Mohr-Coulomb yield surface (on the Standard and Fine meshes a 32- or 64-sided lifted polygon after Ben-Tal and Nemirovski, at about the cost of a 24-sided one). The stress field is extended to the infinite half-space, so the result is a rigorous lower bound at any mesh size. The mesh is refined in adaptive cycles from the plastic multipliers of the previous solve, after OptumG2. The **upper bound** comes from discontinuity layout optimisation (DLO) and, for layered or frictional soil, slopes and combined loads, a finite element upper bound (6-node quadratic velocity triangles on the Standard and Fine meshes); the lower is reported. The circular footing uses axisymmetric finite element bounds of its own. The footing is a rigid body on a no-tension interface. Results are gross total pressures at the footing base. Non-strip shapes apply a classical shape-factor ratio to the plane-strain bounds, which is an estimate, not a rigorous bound. See [docs/method.md](docs/method.md) for detail, verification and known limits.

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

`verify` checks that the bounds bracket exact solutions: Prandtl's Nc, Nq, Martin's Nγ, Green's inclined load (proportional, with H held, and sliding with V held), a crust over soft clay, an eccentric load, a footing at a slope crest, Davis's reduced strength, a tension cut-off, an anisotropic clay, a partly rough base, the strength factor of safety and a circular footing (Eason and Shield's Nc = 6.05).

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

This is an engineering prototype. It has had an independent solver review and an independent UX review, and the findings were fixed (see `docs/method.md`). Its changes since then (inclined and eccentric loads, slopes, and the 2026-10-25 additions from the OptumG2 manual review) have been tested here but not independently reviewed. Check results against established methods before using them for design.
