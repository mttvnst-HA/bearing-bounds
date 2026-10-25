# Method, verification and history

## Lower bound

Finite element limit analysis (FELA) with linear stress triangles, statically admissible stress discontinuities on every edge, and an inscribed linearised yield surface (after Sloan 1988). The stress field is extended to the infinite half-space beside and below the mesh. Beside the mesh it has zero shear and geostatic vertical stress. Below the mesh it has uniform horizontal stress plus equal geostatic increments. Yield is checked exactly (zero shear needs no linearisation) at the boundary nodes and at the top of every deeper layer. This makes the result a rigorous lower bound whatever the mesh size.

a. **Yield surface.** Mohr-Coulomb in plane strain is a cone, ‖(σx − σy, 2τ)‖ ≤ 2c cos φ − (σx + σy) sin φ. The Quick mesh uses an inscribed polygon with p = 18 sides (strength factor cos(π/p) = 0.985). The Standard and Fine meshes use the lifted polygon of Ben-Tal and Nemirovski (2001): ν stages of auxiliary variables give a 2^(ν+1)-sided polygon with about 3ν rows per stress point (ν = 4: 32 sides, factor 0.995; ν = 5: 64 sides, factor 0.9988). The cone is scaled by that factor, which makes the lifted polygon lie inside Mohr-Coulomb. This was checked numerically: in every direction its support lies between 1 and 1/cos(π/2^(ν+1)). Before this change the linearisation was most of the undrained lower-bound error (0.992 against a factor of 0.991).

b. **Other strength models**, as further cones intersected with Mohr-Coulomb: a tension cut-off σ₁ ≤ kₜ, and anisotropic undrained strength (Davis and Christian 1971), an ellipse in ((σx − σy)/2, τ) with semi-axes (sᵤ,c + sᵤ,e)/2 and sᵤ and its centre shifted by (sᵤ,c − sᵤ,e)/2.

c. **Adaptive mesh** (after OptumG2: shear dissipation is "the most efficient and reliable" control, with 3 or 4 cycles). The first mesh is refined along the DLO mechanism. Later cycles refine on the plastic multipliers of the previous lower bound (its LP row duals, times the deviatoric stress: a shear-dissipation measure that also works for cohesionless soil) and on the plastic strain of the finite element upper bound. Standard runs 2 cycles; Fine runs 3. Every mesh gives valid bounds, so the best of each kind is kept.

## Upper bounds

a. **Discontinuity layout optimisation (DLO)** (after Smith and Gilbert 2007). A translational mechanism is formed from candidate discontinuities joining every pair of nodes. Self-weight and surcharge work come from the weight vertically above each line. The footing is a rigid body (sliding, settling, rotating) on a no-tension interface. The domain grows until the mechanism is clear of its edges; then a finer grid is fitted to the mechanism. Each line's strength is the maximum shear stress of the yield surface on that line (which differs for the two slip directions in anisotropic clay), plus an opening term at the tensile strength when there is a cut-off.

b. **Finite element upper bound, linear** (after Sloan and Kleeman 1995): linear velocity triangles with velocity discontinuities and a circumscribed yield surface (a polygon, or the LP dual of the lifted polygon, which gives its support function exactly).

c. **Finite element upper bound, quadratic** (Yu, Sloan and Kleeman 1994; OptumG2's "Upper" element): continuous 6-node triangles. Strain rates are linear, so the flow rule imposed at the vertices holds everywhere, and dissipation and self-weight work are integrated exactly. Used on the Standard and Fine meshes. It runs for layered profiles, slopes, combined loads and frictional soil, and the lowest valid upper bound is reported. It locks somewhat in undrained soil, where the DLO usually governs.

## Circular footings: axisymmetric bounds

A circular footing under a central vertical load on level ground is analysed as axisymmetric, which gives rigorous bounds for the circle itself instead of the shape-factor estimate. Mohr-Coulomb is used in three principal stresses: the hoop stress is principal, and the three pairs give three cones in the meridian plane. A test over 200 000 random stress states found no case where these cones and Mohr-Coulomb disagree.

a. **Lower bound.** The stress is the geostatic part plus a part whose r-weighted components r·(σr, σz, τ) are quadratic in each triangle (after Makrodimopoulos and Martin's use of r-weighted stresses). The vertical equilibrium equation then holds exactly, and the radial one defines the hoop stress, linear in the triangle. Each yield row multiplied by r is a quadratic and is imposed on its six Bernstein coefficients, which is sufficient for it to hold everywhere. The half-space extension is the axisymmetric form of the plane-strain one: beside the mesh σr = σθ, and below it uniform σr = σθ.

b. **Upper bound.** Quadratic velocities in each triangle, with quadratic velocity jumps between triangles; without the jumps the hoop strain locks undrained soil badly (bound 1.8 times the exact value). Plastic multipliers and slip rates are quadratics with non-negative Bernstein coefficients, and the flow rule is matched coefficient by coefficient. Dissipation and work are integrated exactly.

c. The plane-strain DLO still runs, but only to size the model; it is never reported as a bound for the circle. Tension cut-off and anisotropic strength are not available in the axisymmetric analysis. The solves are larger, so the meshes are coarser and the run is slower.

## Loading, ground and stresses

a. Loads are F = F_fixed + α F_ref with F = (V, H, M), and the solvers find the collapse multiplier α (after OptumG2's Fixed and Multiplier loads). There are three modes. "All loads together": V, H and M = V·e are increased in proportion, and the factor of safety is on the whole load. "Increase V, hold H": H is held, for example a wind or seismic load. "Hold V, increase H": V is held at eccentricity e, and the result is the sliding capacity.

b. **Base interface.** Smooth, or rough with an interface factor R: tan δ = R tan φ′ and adhesion = R c′ (or R sᵤ) under the footing. In the lower bound this is |τ| ≤ a − σ tan δ at the footing nodes, which is exact.

c. **Factor of safety on strength** (OptumG2's Strength Reduction analysis). c′, sᵤ and tan φ′ are divided by F with the applied loads held. Any F at which the lower-bound multiplier still reaches the applied load proves FS ≥ F, and any F at which the upper-bound multiplier falls below it proves FS ≤ F. Each bound is therefore the safe end of its own bracket. Secant steps on log F aim just inside the safe side, so a bound usually takes 2 to 4 extra solves. If the bracket is still wide (a slope with little footing load), the adaptive analysis is re-run at the estimated F. "Slope alone" uses an applied load of zero.

d. **Non-associated flow** (Davis 1968). Given a dilation angle ψ′ < φ′, the layer uses c* = β c′ and tan φ* = β tan φ′, with β = cos ψ′ cos φ′ / (1 − sin ψ′ sin φ′); checked against OptumG2's Table 8.2, φ′ = 40° and ψ′ = 10° give φ* = 35.47°. The bounds are rigorous for that reduced material. Its link to the real non-associated soil is an engineering estimate: OptumG2 notes that associated parameters overestimate the bearing capacity and Davis parameters usually underestimate it.

e. Embedment is a surcharge equal to the overburden at the footing base; the strength of soil above the base is ignored (conservative).

f. The analysis is in effective stress, with buoyant unit weight below the water table for every layer. Gross total pressures add the pore pressure at the base.

g. A slope on one side is described by edge-to-crest distance, height and angle. The water table must be at or below the toe. When the soil has weight, the mesh reaches 4 slope heights below the toe (see Known limits h). Optional mesh fans at the crest and toe (OptumG2 "Mesh Fans at singular points") are available; in tests they did not tighten the bounds, so they are off by default.

h. Non-strip footings other than the axisymmetric circle scale the plane-strain bounds by the ratio of classical bearing capacity with and without Vesić/De Beer shape factors (with effective width B − 2e). This is an estimate, not a rigorous bound.

i. The summary also shows the mean of the bounds with ± half the gap, labelled as an estimate (from OptumG2's verification and validation section). Where the upper bound is further from the true value than the lower bound, as for self-weight in cohesionless soil, the mean reads high.

j. LP status: a bound is accepted only from a feasible primal point. If the interior-point solver stops without proving optimality and gives a lower bound of zero or less, the LP is re-solved with crossover.

## Verification (Standard mesh)

Standard mesh: lifted 32-sided yield polygon, 2 adaptive cycles, 6-node finite element upper bound. The "before" columns are the 2026-10-24 values (24-sided polygon, one mesh).

| Case | Reference | LB ÷ ref. | UB ÷ ref. | Before: LB / UB |
|---|---|---|---|---|
| Undrained, weightless, rough (Nc = 5.14) | 5.14 sᵤ | 0.994 | 1.009 | 0.992 / 1.009 |
| c′ = 10, φ′ = 20°, weightless (Nc = 14.8) | 148.3 | 0.980 | 1.023 | 0.969 / 1.023 |
| Surcharge only, φ′ = 30°, smooth (Nq = 18.4) | 184.0 | 0.886 | 1.026 | 0.884 / 1.026 |
| Self-weight, φ′ = 30°, rough (Nγ = 14.75, Martin) | 265.5 | 0.927 | 1.108 | 0.910 / 1.189 |
| Self-weight, φ′ = 40°, rough (Nγ = 85.6, Martin) | 1541 | 0.887 | 1.171 | 0.833 / 1.266 |
| Crust: sᵤ 100 to 0.5 m over sᵤ 10, B = 2 m (true load ≤ 110.7) | ≤ 110.7 | 0.915 | 1.001 | 0.915 / 1.052 |
| Inclined, undrained, H/(B·sᵤ) = 0.5 (Green 1954) | 448.4 | 0.993 | 1.013 | 0.994 / 1.013 |
| Same, H = 100 held and V increased | 448.4 | 0.992 | 1.016 | new |
| Same, V held and H increased (sliding) | H = 100 | 0.947 | 1.085 | new |
| Eccentric, undrained, e = 0.1B (Meyerhof effective width) | 411.4 | 1.011 | 1.056 | 1.006 / 1.062 |
| Footing at crest of a 30° undrained slope, 4 m high, weightless | no closed form | Nc 4.06 | Nc 4.13 | 4.06 / 4.13 |
| Davis: c′ = 10, φ′ = 25°, ψ′ = 5°, weightless (c* Nc(φ*)) | 176.4 | 0.976 | 1.022 | new |
| Tension cut-off kₜ = 0, c′ = 10, φ′ = 20° (Prandtl is an upper limit) | ≤ 148.3 | 0.975 | 1.023 | new |
| Anisotropic clay, sᵤ,c = 1.6 sᵤ, sᵤ,e = 0.6 sᵤ (no closed form) | – | 546.8 | 555.2 | new |
| Inclined as Green, interface factor R = 0.6 (no closed form) | – | 445.4 | 454.3 | new |
| Strength factor of safety at half the Prandtl load (exact 2.0) | 2.0 | 0.990 | 1.013 | new |
| Circular, rough, undrained, axisymmetric (Nc = 6.05, Eason and Shield) | 605 | 0.968 | 1.048 | new |

Frictional circular footings have no exact reference here. Two consistency checks (axisymmetric, about 650 elements): a weightless smooth circle at φ′ = 30° with surcharge gives N_q from 25.8 to 34.1, which brackets the usual estimate of about 29 (strip 18.4 × Vesić s_q = 1.58); a rough circle on weighty sand at φ′ = 30° gives N_γ (on ½γB) from 12.1 to 20.1. The upper bound is loose for frictional self-weight. These frictional axisymmetric results need the independent review in the backlog before design use.

The eccentric lower bound sits slightly above the effective-width estimate, consistent with published limit-analysis studies that find the effective-width method slightly conservative. A Standard run of one of these cases now takes 2 to 7 minutes on a desktop (before: 1 to 2).

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

f. **OptumG2 manual review and changes (2026-10-25).** A review of the OptumG2 Analysis, Theory, Materials, Features and Results manuals led to twelve recommendations, all implemented: the lifted yield polygon, adaptive mesh cycles, the 6-node finite element upper bound, fixed and multiplied loads, the base interface factor, Davis parameters, the strength factor of safety, optional slope fans, a tension cut-off, the mean of the bounds, axisymmetric bounds for circular footings, and anisotropic undrained strength with a note that sᵤ is the plane-strain (direct simple shear) strength. One bug found during the work (misaligned auxiliary coefficients in the lifted polygon, which made the lower bound zero) was fixed before release. These changes have been verified against the cases below but not independently reviewed.

## Known limits

a. Plane strain, except the axisymmetric circular footing. Square, rectangular, mat and slab footings, and circles under H, e or near a slope, use classical shape factors (an estimate).

b. Soil above the footing base is a surcharge with no strength. This is especially conservative near slopes.

c. Self-weight bounds in cohesionless soil still converge slowly (see the verification table), though the lifted polygon, adaptive cycles and the 6-node upper bound have narrowed the gap.

d. H and e act across B only. No moment about the other axis.

e. Water table must be at or below the slope toe.

f. No settlement or load–displacement behaviour; limit analysis gives collapse only.

g. Associated flow. Davis parameters (when ψ′ < φ′) give an estimate of non-associated behaviour, not a bound on it.

h. **Lower bound beside slopes in clay.** The half-space extension has uniform horizontal stress at depth, while the far-field vertical stress differs by γ·H_slope between crest and toe. For an undrained slope this limits any lower bound to γ·H_slope ≤ 4 sᵤ, whatever the mesh. The slope's own strength factor of safety is therefore bracketed loosely, with the lower end at most 4 sᵤ / (γ H_slope). A rigorous extension whose stresses decay with depth would remove this limit.

i. Tension cut-off and anisotropic strength are not available in the axisymmetric analysis.

j. The axisymmetric and Fine runs are slow (minutes) because their LPs are large.

## Backlog (suggested next steps)

a. Independent review of the inclined, eccentric and slope code, and of the 2026-10-25 additions (axisymmetric bounds, lifted polygon, strength factor of safety, cones).

b. Model the soil above the footing base as real soil (footing in an excavation).

c. Tighter self-weight upper bound: more DLO nodes near the footing edge, or adaptive member addition.

d. Load–displacement (serviceability) analysis.

e. A decaying half-space extension for slopes (Known limits h).

f. Speed: the axisymmetric LPs, and a conic solver if one becomes available in WebAssembly.
