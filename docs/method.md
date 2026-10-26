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

j. LP status: a bound is accepted only from a feasible primal point. If the interior-point solver stops without proving optimality and gives a lower bound of zero or less, the LP is re-solved with crossover. The interior-point optimality tolerance is 1e-3: a bound needs only a feasible point, and tests showed the same feasibility (1e-7 or better) with objectives within 0.07% and solves 10 to 20% faster.

k. Speed: in the browser app, the lower and upper bound of each adaptive cycle are solved at the same time, each in its own child worker with its own WebAssembly copy of HiGHS. This roughly halves a Standard run (5:23 to 3:03 on a test case). If child workers are not available the solves run one after the other, and the solution detail says so. In Node (the tests) they always run one after the other.

## Exact cones and the Clarabel solver

The finite element bounds (plane-strain lower bound, linear and 6-node upper bounds, and both axisymmetric bounds) are first solved with the yield condition stated exactly, as second-order cones (SOC), by Clarabel (an interior-point conic solver from Oxford, compiled to WebAssembly; see vendor/clarabel/README.md for the version and licence). HiGHS with the linearised yield surface is the fallback. The discontinuity layout bound is a linear program (LP) by nature and always uses HiGHS.

a. **Lower bound cones.** Each stress point gets (2c cos φ − (σx + σy) sin φ, σx − σy, 2τ) in the cone; the other strength models are further cones, as in the Lower bound section. No strength factor is needed: the cone is the yield condition itself. In the axisymmetric lower bound, r times (t, x₁, x₂) of each of the three cones is a quadratic in the triangle, and each of its six Bernstein-coefficient triples is put in the cone. The quadratic at any point is a convex combination of those triples (the Bernstein basis is non-negative and sums to one), so it lies in the convex cone everywhere.

b. **Upper bound cones.** By conic duality, for any y = (y₀, y₁, y₂) in the cone with plastic strain rates (εx, εy, γxy) = (y₀ sin φ − y₁, y₀ sin φ + y₁, −2 y₂) (plus the hoop term in the axisymmetric case), the dissipation is at most 2c cos φ · y₀. Each polygon multiplier of the linearised upper bound becomes one such triple, so the flow rule and the dissipation integrals are unchanged in form, and in the axisymmetric upper bound each Bernstein coefficient gets its own triple.

c. **The cone apex.** A traction-free surface point of cohesionless soil with no surcharge can carry no stress at all: the exact cone leaves only its apex. The problem then has no strict interior, and Clarabel stalls there. In cone mode those points are fixed to zero before the solve, and so are the other stress points at the same surface nodes (for example the footing corner, where the optimum is at the apex too). Fixing stresses only restricts the field, so the bound stays rigorous; the cones that become constants are checked and dropped. In the axisymmetric bound the hoop stress at those points is also set to zero (the cones force that too). A zero tension cut-off at a traction-free point leaves its cone only a ray (σ along the surface ≤ 0); that linear condition is solved instead and the cone is still checked.

d. **Admissibility check.** Clarabel's point is accepted only after every constraint of the problem is checked at that point, with the absolute 1e-6 rule that the HiGHS path uses: every equality and inequality row (equilibrium, tractions, boundary and footing conditions, the extension beside and below the mesh, the flow rule and footing interface in the upper bound), every variable bound, and every cone, exactly. The bound is recomputed from that point. To make the check pass at points where the yield condition is active, Clarabel solves with every cone narrowed by a relative 1e-6; the check uses the exact cone. The point is accepted when Clarabel reports "Solved" or "AlmostSolved" and the check passes. Otherwise (any other status, a failed check, or no Clarabel module) the solver runs again with the linearised yield surface and HiGHS, and the solution detail says why.

e. **Numerical settings.** Gap tolerances 1e-7, feasibility tolerance 1e-9, at most 100 iterations. If the first solve gives no accepted point (a numerical error, or an "almost solved" stall whose point misses the 1e-6 check), it is repeated once with static regularisation 1e-7 instead of Clarabel's default 1e-8. The axisymmetric lower bound and some 6-node upper bounds (eccentric load, slope) needed that; with it they converge cleanly in 25 to 35 iterations.

f. **Effect.** See the Verification table. A cone solve takes 20 to 70 interior-point iterations and a few seconds where HiGHS took one to eight minutes on the Standard mesh, and the bounds are tighter because nothing is lost to the polygon: the Standard lower bounds rise by up to 2% (sliding 0.947 to 0.966, Nγ 0.926 to 0.942, circle 0.968 to 0.977).

g. **Remaining fallback.** A zero tension cut-off (kₜ = 0) in weightless soil leaves large zones at the apex of the cut-off cone (zero stress beside the Prandtl zone); Clarabel stalls there and the lower bound falls back to HiGHS. The result is still a valid bound; it just takes the HiGHS time plus the failed cone attempt.

## Verification (Standard mesh)

Standard mesh: exact cones (Clarabel) with the 32-sided lifted polygon and HiGHS as fallback, 2 adaptive cycles, 6-node finite element upper bound. The "HiGHS only" column is the same run with `--no-cone` (the linearised yield surface throughout); it agrees with the previous version of this table to within 0.002. Times are the whole case (all solves, both cycles) in seconds on the owner's laptop (Intel i5-1035G7, Node 24, one solve at a time), measured 2026-10-25.

| Case | Reference | LB ÷ ref. | UB ÷ ref. | HiGHS only: LB / UB | Time, cones (s) | Time, HiGHS only (s) |
|---|---|---|---|---|---|---|
| Undrained, weightless, rough (Nc = 5.14) | 5.14 sᵤ | 0.997 | 1.010 | 0.994 / 1.010 | 40 | 113 |
| c′ = 10, φ′ = 20°, weightless (Nc = 14.8) | 148.3 | 0.987 | 1.023 | 0.980 / 1.023 | 72 | 360 |
| Surcharge only, φ′ = 30°, smooth (Nq = 18.4) | 184.0 | 0.892 | 1.017 | 0.887 / 1.027 | 52 | 1708 |
| Self-weight, φ′ = 30°, rough (Nγ = 14.75, Martin) | 265.5 | 0.942 | 1.096 | 0.926 / 1.108 | 49 | 446 |
| Crust: sᵤ 100 to 0.5 m over sᵤ 10, B = 2 m (true load ≤ 110.7) | ≤ 110.7 | 0.917 | 0.999 | 0.915 / 1.002 | 18 | 239 |
| Inclined, undrained, H/(B·sᵤ) = 0.5 (Green 1954) | 448.4 | 0.996 | 1.014 | 0.993 / 1.014 | 22 | 422 |
| Same, H = 100 held and V increased | 448.4 | 0.995 | 1.016 | 0.992 / 1.016 | 19 | 295 |
| Same, V held and H increased (sliding) | H = 100 | 0.966 | 1.085 | 0.947 / 1.085 | 20 | 332 |
| Eccentric, undrained, e = 0.1B (Meyerhof effective width) | 411.4 | 1.014 | 1.053 | 1.010 / 1.057 | 43 | 1083 |
| Footing at crest of a 30° undrained slope, 4 m high, weightless | no closed form | 407.6 | 413.2 | 406.2 / 413.2 | 34 | 272 |
| Davis: c′ = 10, φ′ = 25°, ψ′ = 5°, weightless (c* Nc(φ*)) | 176.4 | 0.984 | 1.022 | 0.976 / 1.022 | 17 | 345 |
| Tension cut-off kₜ = 0, c′ = 10, φ′ = 20° (Prandtl is an upper limit) | ≤ 148.3 | 0.969 (HiGHS) | 1.023 | 0.975 / 1.023 | 405 | 765 |
| Anisotropic clay, sᵤ,c = 1.6 sᵤ, sᵤ,e = 0.6 sᵤ (no closed form) | – | 548.4 | 555.4 | 546.8 / 555.4 | 36 | 347 |
| Inclined as Green, interface factor R = 0.6 (no closed form) | – | 446.7 | 454.5 | 445.4 / 454.5 | 75 | 873 |
| Strength factor of safety at half the Prandtl load (exact 2.0) | 2.0 | 0.993 | 1.014 | 0.990 / 1.014 | 80 | 249 |
| Circular, rough, undrained, axisymmetric, FS at half capacity (exact 2.0) | 2.0 | 0.973 | 1.045 | 0.964 / 1.053 | 81 | 540 |
| Circular, rough, undrained, axisymmetric (Nc = 6.05, Eason and Shield) | 605 | 0.977 | 1.041 | 0.968 / 1.049 | 24 | 172 |

Every bound in the cone column came from Clarabel except the zero tension cut-off lower bound, which falls back to HiGHS (see Exact cones g). The whole Standard verification (17 analysed cases plus two fixed-mesh regression cases) took 18.8 minutes with cones against 143 minutes with HiGHS only (7.6 times faster; 11 times without the cut-off case). On the Quick mesh it took 5.0 minutes against 14.7. The φ′ = 40° self-weight case (run only on request) was not re-run; its last values were LB 0.887 and UB 1.171 (HiGHS). The "toe" regression case (a slope no deeper than the embedment) was added after these timings. Every lower bound is equal or higher with cones, except the cut-off case, whose HiGHS fallback ran on a slightly different adaptive mesh; every upper bound is equal or lower.

Frictional circular footings have no exact reference here. Two consistency checks (axisymmetric, about 650 elements): a weightless smooth circle at φ′ = 30° with surcharge gives N_q from 25.8 to 34.1, which brackets the usual estimate of about 29 (strip 18.4 × Vesić s_q = 1.58); a rough circle on weighty sand at φ′ = 30° gives N_γ (on ½γB) from 12.1 to 20.1. The upper bound is loose for frictional self-weight. These frictional axisymmetric results need the independent review in the backlog before design use.

The eccentric lower bound sits slightly above the effective-width estimate, consistent with published limit-analysis studies that find the effective-width method slightly conservative. With exact cones a Standard run of one of these cases takes about 15 to 80 seconds in Node (2 to 30 minutes with HiGHS only).

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

f. **OptumG2 manual review and changes (2026-10-25).** A review of the OptumG2 Analysis, Theory, Materials, Features and Results manuals led to twelve recommendations, all implemented: the lifted yield polygon, adaptive mesh cycles, the 6-node finite element upper bound, fixed and multiplied loads, the base interface factor, Davis parameters, the strength factor of safety, optional slope fans, a tension cut-off, the mean of the bounds, axisymmetric bounds for circular footings, and anisotropic undrained strength with a note that sᵤ is the plane-strain (direct simple shear) strength. One bug found during the work (misaligned auxiliary coefficients in the lifted polygon, which made the lower bound zero) was fixed before release. These changes were verified against the cases below and then reviewed independently (item g).

g. **Independent review of PR #1 (2026-10-25).** Three critical lower-bound bugs were found and fixed, each with a regression case that fails on the earlier code:

1. A layer whose top lay exactly at the mesh base was never yield-checked below the mesh (issue #2). In the test the LB was 389.8 kPa against a valid UB of 171.1.
2. The stress field beside the mesh assumed the full surcharge where a slope lowers the ground without reaching the base level. That inflated an LB by 18% in one test and made the LB problem infeasible in another.
3. The strength factor of safety of a circular (axisymmetric) footing used the plane-strain solvers and units. It reported 4.1 to 4.4 against an exact 2.0.

Three guards were also added: a negative sᵤ gradient is rejected, the axisymmetric analysis no longer falls back to the non-conforming quadtree mesh, and an "unbounded or infeasible" status no longer counts as proof of collapse.

A second review of these fixes found two more lower-bound gaps, also fixed:

1. Where a slope is no deeper than the embedment, the mesh could stop short of the slope toe, so the stress field beside it (which is uniform in x) assumed a surcharge the ground no longer carries beyond the edge. In a test with a hand-set domain the LB was 239.5 kN/m against an upper bound of 234.3. The domain now always reaches past the toe, and the lower bound refuses a domain that does not (regression case "toe").
2. An infeasible lower-bound problem was reported as a lower bound of zero. It proves nothing, so no lower bound is now reported.

Open minor points:

1. Layer bottoms and the water table are snapped by up to 0.01B to the base level or to each other.
2. Strength reduction leaves the tensile limit unreduced.

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

j. Fine runs are still slow (minutes): the discontinuity layout LP and any HiGHS fallback are large.

## Backlog (suggested next steps)

a. Independent review of the inclined, eccentric and slope code, and of the 2026-10-25 additions (axisymmetric bounds, lifted polygon, strength factor of safety, cones).

b. Model the soil above the footing base as real soil (footing in an excavation).

c. Tighter self-weight upper bound: more DLO nodes near the footing edge, or adaptive member addition.

d. Load–displacement (serviceability) analysis.

e. A decaying half-space extension for slopes (Known limits h).

f. Speed: the remaining HiGHS fallback for a zero tension cut-off (Exact cones g), and the discontinuity layout LP, which is now the slowest part of a Quick run.
