// Solver verification: runs the engine in Node against known solutions and checks that the
// bounds bracket them. Usage: node tests/verify.js [--std] [case,case,...]
// Quick mesh by default (a few minutes in total); --std uses the Standard mesh (slower, tighter).
const path = require("path");
const FELA = require(path.join(__dirname, "..", "src", "engine.js"));

const Nq = (p) => { const f = p * Math.PI / 180; return Math.exp(Math.PI * Math.tan(f)) * Math.tan(Math.PI / 4 + f / 2) ** 2; };
const L1 = (o) => [Object.assign({ bottom: 1e9, gammaSat: o.gamma }, o)];
const greenV = (h) => 1 + Math.PI / 2 + Math.acos(h) + Math.sqrt(1 - h * h); // V/(B su) for H/(B su) = h
const NcOf = (p) => (p < 1e-9 ? 2 + Math.PI : (Nq(p) - 1) / Math.tan(p * Math.PI / 180));
// Davis reduced strength for c' = 10, phi' = 25, psi' = 5 (weightless, so the exact load is c* Nc(phi*))
const davis = (() => { const f = 25 * Math.PI / 180, s = 5 * Math.PI / 180, b = Math.cos(s) * Math.cos(f) / (1 - Math.sin(s) * Math.sin(f)); return { c: 10 * b, phi: Math.atan(b * Math.tan(f)) * 180 / Math.PI }; })();
// rigorous: the reference is exact, so LB <= ref <= UB must hold
// lbOnly: the reference is itself an upper bound on the true load, so only LB <= ref must hold
// approx: the reference is an engineering estimate; only LB <= UB is checked
// q: "H" compares the horizontal load per unit length (sliding mode); fs: the strength factor of safety brackets ref
// axi: axisymmetric (circular footing); the reference is a pressure on the circle
const CASES = {
  tresca: { m: { B: 2, D: 0, rough: true, layers: L1({ c: 100, phi: 0, gamma: 0 }) }, ref: 514.2, kind: "rigorous", note: "Prandtl Nc = 5.14" },
  nc20: { m: { B: 2, D: 0, rough: true, layers: L1({ c: 10, phi: 20, gamma: 0 }) }, ref: 148.35, kind: "rigorous", note: "c' = 10, phi' = 20, weightless" },
  nq30: { m: { B: 2, D: 0, rough: false, surcharge: 10, layers: L1({ c: 0, phi: 30, gamma: 0 }) }, ref: 10 * Nq(30), kind: "rigorous", note: "Nq, phi' = 30, smooth" },
  ng30: { m: { B: 2, D: 0, rough: true, layers: L1({ c: 0, phi: 30, gamma: 18 }) }, ref: 18 * 14.75, kind: "rigorous", note: "Ngamma, phi' = 30, rough (Martin)" },
  ng40: { m: { B: 2, D: 0, rough: true, layers: L1({ c: 0, phi: 40, gamma: 18 }) }, ref: 18 * 85.6, kind: "rigorous", only: true, note: "Ngamma, phi' = 40, rough (Martin)" },
  crust: { m: { B: 2, D: 0, rough: true, layers: [{ bottom: 0.5, c: 100, phi: 0, gamma: 18, gammaSat: 18 }, { bottom: 1e9, c: 10, phi: 0, gamma: 18, gammaSat: 18 }] }, ref: 110.7, kind: "lbOnly", note: "stiff crust over soft clay" },
  green: { m: { B: 2, D: 0, rough: true, load: { tanA: 0.5 / greenV(0.5) }, layers: L1({ c: 100, phi: 0, gamma: 0 }) }, ref: 100 * greenV(0.5), kind: "rigorous", note: "inclined, Green (1954)" },
  ecc: { m: { B: 2, D: 0, rough: true, load: { e: 0.2 }, layers: L1({ c: 100, phi: 0, gamma: 0 }) }, ref: 514.2 * 0.8, kind: "approx", note: "eccentric e = 0.1B (Meyerhof)" },
  slope: { m: { B: 2, D: 0, rough: true, slope: { dist: 0, height: 4, angle: 30 }, layers: L1({ c: 100, phi: 0, gamma: 0 }) }, ref: NaN, kind: "approx", note: "footing at crest of 30 deg slope" },
  fixH: { m: { B: 2, D: 0, rough: true, load: { tanA: 0, Hfix: 100 }, layers: L1({ c: 100, phi: 0, gamma: 0 }) }, ref: 100 * greenV(0.5), kind: "rigorous", note: "H = 100 held, V increased (Green)" },
  slide: { m: { B: 2, D: 0, rough: true, load: { mode: "H", Vfix: 200 * greenV(0.5) }, layers: L1({ c: 100, phi: 0, gamma: 0 }) }, ref: 100, q: "H", kind: "rigorous", note: "V held, H increased (Green)" },
  davis: { m: { B: 2, D: 0, rough: true, layers: L1({ c: 10, phi: 25, psi: 5, gamma: 0 }) }, ref: davis.c * NcOf(davis.phi), kind: "rigorous", note: "Davis c*, phi* for phi' 25, psi' 5" },
  cutoff: { m: { B: 2, D: 0, rough: true, layers: L1({ c: 10, phi: 20, kt: 0, gamma: 0 }) }, ref: 148.35, kind: "lbOnly", note: "tension cut-off (Prandtl load is an upper limit)" },
  aniso: { m: { B: 2, D: 0, rough: true, layers: L1({ c: 100, phi: 0, gamma: 0, aniso: { rc: 1.6, re: 0.6 } }) }, ref: NaN, kind: "approx", note: "anisotropic su, compression 1.6, extension 0.6" },
  partial: { m: { B: 2, D: 0, rough: true, rInter: 0.6, load: { tanA: 0.5 / greenV(0.5) }, layers: L1({ c: 100, phi: 0, gamma: 0 }) }, ref: NaN, kind: "approx", note: "inclined, interface factor 0.6" },
  sfs: { m: { B: 2, D: 0, rough: true, layers: L1({ c: 100, phi: 0, gamma: 0 }) }, opts: { strength: { target: 514.2 } }, ref: 2.0, kind: "fs", note: "strength FS, half the Prandtl load (exact 2.0)" },
  circle: { m: { B: 2, D: 0, rough: true, axisym: true, layers: L1({ c: 100, phi: 0, gamma: 0 }) }, ref: 605, kind: "rigorous", axi: true, note: "circular, rough, undrained (Nc = 6.05)" },
};
// the same settings as the app (src/app.html QOPTS and QAXI)
const OPTS = {
  quick: { p: 18, dloNodes: 260, targetElems: 1100, nTheta: 12, fanR0: 0.08, feAutoCombined: false, ipmTol: 1e-5 },
  std: { p: 24, nu: 4, dloNodes: 420, targetElems: 1800, adaptCycles: 2, feQuad: true, feFrictional: true, ipmTol: 1e-5 },
};
const AXI = {
  quick: { p: 12, targetElems: 350, nTheta: 6, nf: 2, fanR0: 0.1, adaptCycles: 1, nu: 0 },
  std: { p: 16, targetElems: 650, nTheta: 8, nf: 3, adaptCycles: 1, nu: 0 },
};

(async () => {
  const args = process.argv.slice(2);
  const mode = args.includes("--std") ? "std" : "quick";
  const pick = args.find((a) => !a.startsWith("--"));
  const names = pick ? pick.split(",") : Object.keys(CASES).filter((k) => !CASES[k].only);   // "only" cases run when named
  const highs = await require("highs")();
  let fails = 0;
  for (const k of names) {
    const c = CASES[k];
    if (!c) { console.log(`unknown case ${k}`); fails++; continue; }
    const t = Date.now();
    const o = Object.assign({}, OPTS[mode], c.axi ? AXI[mode] : {}, c.opts || {});
    if (!o.nu) delete o.nu;
    const r = await FELA.analyse(highs, Object.assign({ gammaW: 9.81, waterDepth: null }, c.m), o);
    const ref = c.ref;
    const lb = c.kind === "fs" ? r.strength.fsLb : c.q === "H" ? r.best.Hlb : r.best.lb, ub = c.kind === "fs" ? r.strength.fsUb : c.q === "H" ? r.best.Hub : r.best.ub;
    const checks = [["LB <= UB", lb <= ub * (1 + 1e-6)]];
    if (c.kind === "rigorous" || c.kind === "fs") checks.push(["LB <= ref", lb <= ref * (1 + 1e-6)], ["UB >= ref", ub >= ref * (1 - 1e-6)]);
    if (c.kind === "lbOnly") checks.push(["LB <= ref", lb <= ref * (1 + 1e-6)]);
    const bad = checks.filter(([, ok]) => !ok).map(([n]) => n);
    if (bad.length) fails++;
    const ratio = (v) => (isFinite(ref) ? (v / ref).toFixed(3) : v.toFixed(1));
    console.log(`${bad.length ? "FAIL" : "ok  "} ${k.padEnd(7)} LB ${ratio(lb)}  UB ${ratio(ub)}  ${((Date.now() - t) / 1000).toFixed(0)} s  (${c.note})${bad.length ? "  failed: " + bad.join(", ") : ""}`);
  }
  console.log(fails ? `${fails} case(s) failed` : "All cases passed");
  process.exit(fails ? 1 : 0);
})();
