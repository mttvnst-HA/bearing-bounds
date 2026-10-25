# Method, verification and history

## Lower bound

Finite element limit analysis (FELA) with linear stress triangles, statically admissible stress discontinuities on every edge, and an inscribed linearised Mohr-Coulomb yield surface (after Sloan 1988). The stress field is extended to the infinite half-space beside and below the mesh. Beside the mesh it has zero shear and geostatic vertical stress. Below the mesh it has uniform horizontal stress plus equal geostatic increments. Yield is checked at the boundary nodes and at the top of every deeper layer. This makes the result a rigorous lower bound whatever the mesh size.

## Upper bounds

a. **Discontinuity layout optimisation (DLO)** (after Smith and Gilbert 2007). A translational mechanism is formed from candidate discontinuities joining every pair of nodes. Self-weight and surcharge work come from the weight vertically above each line. The footing is a rigid body (sliding, settling, rotating) on a no-tension interface. The domain grows until the mechanism is clear of its edges; then a finer grid is fitted to the mechanism.

b. **Finite element upper bound** (after Sloan and Kleeman 1995): linear velocity triangles with velocity discontinuities and a circumscribed yield surface. It runs automatically for layered profiles, and for slopes and combined loads on the Standard and Fine meshes. The lower of the two upper bounds is reported.

## Loading, ground and stresses

a. Vertical load V, horizontal load H and moment M = V·e are increased together until collapse. The factor of safety is the factor on that whole load.

b. Embedment is a surcharge equal to the overburden at the footing base; the strength of soil above the base is ignored (conservative).

c. The analysis is in effective stress, with buoyant unit weight below the water table for every layer. Gross total pressures add the pore pressure at the base.

d. A slope on one side is described by edge-to-crest distance, height and angle. The water table must be at or below the toe.

e. Non-strip footings scale the plane-strain bounds by the ratio of classical bearing capacity with and without Vesić/De Beer shape factors (with effective width B − 2e). This is an estimate, not a rigorous bound.

## Verification (Standard mesh, strip footing)

| Case | Reference | LB ÷ ref. | UB ÷ ref. |
|---|---|---|---|
| Undrained, weightless, rough (Nc = 5.14) | 5.14 sᵤ | 0.992 | 1.009 |
| c′ = 10, φ′ = 20°, weightless (Nc = 14.8) | 148.3 | 0.969 | 1.023 |
| Surcharge only, φ′ = 30°, smooth (Nq = 18.4) | 184.0 | 0.884 | 1.026 |
| Self-weight, φ′ = 30°, rough (Nγ = 14.75, Martin) | 265.5 | 0.910 | 1.189 |
| Self-weight, φ′ = 40°, rough (Nγ = 85.6, Martin) | 1541 | 0.833 | 1.266 |
| Crust: sᵤ 100 to 0.5 m over sᵤ 10, B = 2 m (true load ≤ 110.7) | ≤ 110.7 | 0.915 | 1.052 |
| Inclined, undrained, H/(B·sᵤ) = 0.5 (Green 1954) | 448.4 | 0.994 | 1.013 |
| Eccentric, undrained, e = 0.1B (Meyerhof effective width) | 411.4 | 1.006 | 1.062 |
| Footing at crest of a 30° undrained slope, 4 m high, weightless | no closed form | Nc 4.06 | Nc 4.13 |

The eccentric lower bound sits slightly above the effective-width estimate, consistent with published limit-analysis studies that find the effective-width method slightly conservative.

## Review history

a. **Independent solver review (2026-10-24).** Found that the lower bound was not rigorous. The domain was fitted to the mechanism with no stress-field extension; unconverged interior-point solutions were accepted; elements straddled layer boundaries; and the surcharge mixed total and effective stress. All were fixed: extension elements, domain growth, strict LP status checks, a layer-conforming mesh and an effective-stress basis. The review confirmed the core formulations (equilibrium, discontinuities, yield polygons, flow rule and DLO compatibility) were correct.

b. **Inclined, eccentric and slope loading added (2026-10-24).** Rigid footing with a no-tension interface in all three solvers; full-section models; slope geometry. During testing, two bugs were found and fixed: the self-weight integration next to the footing edge in the DLO, and a missing symmetry condition. These changes have been tested but not independently reviewed.

c. **Crash fix (2026-10-24).** The solver moved to a disposable Web Worker after it froze and crashed Claude Desktop.

d. **Independent UX review (2026-10-25).** Ten findings, all implemented:

1. results tied to their inputs and marked stale;
2. load conversion when the footing type changes;
3. SI storage so unit switching is lossless;
4. an answer-first summary;
5. a clearer bound chart with a classical-difference note;
6. plain-language run progress and field-level errors;
7. layer boundaries on every view and slope dimensions;
8. collapsible layers and method notes, and a sticky Run bar;
9. clearer layer cards;
10. accessibility fixes.

e. **Collapse animation (2026-10-25).** The load ramps up through the bounds, then the upper-bound mechanism plays out with exaggerated motion.

## Known limits

a. Plane strain only. 3D shapes use classical shape factors (an estimate).

b. Soil above the footing base is a surcharge with no strength. This is especially conservative near slopes.

c. Self-weight bounds in cohesionless soil converge slowly: about −10 to −17% / +19 to +27% at Standard.

d. H and e act across B only. No moment about the other axis.

e. Water table must be at or below the slope toe.

f. No settlement or load–displacement behaviour; limit analysis gives collapse only.

## Backlog (suggested next steps)

a. Independent review of the inclined, eccentric and slope code.

b. Model the soil above the footing base as real soil (footing in an excavation).

c. Tighter self-weight upper bound: more DLO nodes near the footing edge, or adaptive member addition.

d. Load–displacement (serviceability) analysis.
