// Solver verification: runs the engine in Node against known solutions and checks that the
// bounds bracket them. Usage: node tests/verify.js [--std] [case,case,...]
// Quick mesh by default (a few minutes in total); --std uses the Standard mesh (slower, tighter).
const path = require("path");
const FELA = require(path.join(__dirname, "..", "src", "engine.js"));

const Nq = (p) => { const f = p * Math.PI / 180; return Math.exp(Math.PI * Math.tan(f)) * Math.tan(Math.PI / 4 + f / 2) ** 2; };
const L1 = (o) => [Object.assign({ bottom: 1e9, gammaSat: o.gamma }, o)];
const greenV = (h) => 1 + Math.PI / 2 + Math.acos(h) + Math.sqrt(1 - h * h); // V/(B su) for H/(B su) = h
// rigorous: the reference is exact, so LB <= ref <= UB must hold
// lbOnly: the reference is itself an upper bound on the true load, so only LB <= ref must hold
// approx: the reference is an engineering estimate; only LB <= UB is checked
const CASES = {
  tresca: { m: { B: 2, D: 0, rough: true, layers: L1({ c: 100, phi: 0, gamma: 0 }) }, ref: 514.2, kind: "rigorous", note: "Prandtl Nc = 5.14" },
  nc20: { m: { B: 2, D: 0, rough: true, layers: L1({ c: 10, phi: 20, gamma: 0 }) }, ref: 148.35, kind: "rigorous", note: "c' = 10, phi' = 20, weightless" },
  nq30: { m: { B: 2, D: 0, rough: false, surcharge: 10, layers: L1({ c: 0, phi: 30, gamma: 0 }) }, ref: 10 * Nq(30), kind: "rigorous", note: "Nq, phi' = 30, smooth" },
  ng30: { m: { B: 2, D: 0, rough: true, layers: L1({ c: 0, phi: 30, gamma: 18 }) }, ref: 18 * 14.75, kind: "rigorous", note: "Ngamma, phi' = 30, rough (Martin)" },
  crust: { m: { B: 2, D: 0, rough: true, layers: [{ bottom: 0.5, c: 100, phi: 0, gamma: 18, gammaSat: 18 }, { bottom: 1e9, c: 10, phi: 0, gamma: 18, gammaSat: 18 }] }, ref: 110.7, kind: "lbOnly", note: "stiff crust over soft clay" },
  green: { m: { B: 2, D: 0, rough: true, load: { tanA: 0.5 / greenV(0.5) }, layers: L1({ c: 100, phi: 0, gamma: 0 }) }, ref: 100 * greenV(0.5), kind: "rigorous", note: "inclined, Green (1954)" },
  ecc: { m: { B: 2, D: 0, rough: true, load: { e: 0.2 }, layers: L1({ c: 100, phi: 0, gamma: 0 }) }, ref: 514.2 * 0.8, kind: "approx", note: "eccentric e = 0.1B (Meyerhof)" },
  slope: { m: { B: 2, D: 0, rough: true, slope: { dist: 0, height: 4, angle: 30 }, layers: L1({ c: 100, phi: 0, gamma: 0 }) }, ref: NaN, kind: "approx", note: "footing at crest of 30 deg slope" },
};
const OPTS = {
  quick: { p: 18, dloNodes: 260, targetElems: 1100, nTheta: 12, fanR0: 0.08, feAutoCombined: false },
  std: { p: 24, dloNodes: 420, targetElems: 1800 },
};

(async () => {
  const args = process.argv.slice(2);
  const mode = args.includes("--std") ? "std" : "quick";
  const pick = args.find((a) => !a.startsWith("--"));
  const names = pick ? pick.split(",") : Object.keys(CASES);
  const highs = await require("highs")();
  let fails = 0;
  for (const k of names) {
    const c = CASES[k];
    if (!c) { console.log(`unknown case ${k}`); fails++; continue; }
    const t = Date.now();
    const r = await FELA.analyse(highs, Object.assign({ gammaW: 9.81, waterDepth: null }, c.m), Object.assign({}, OPTS[mode]));
    const lb = r.best.lb, ub = r.best.ub, ref = c.ref;
    const checks = [["LB <= UB", lb <= ub * (1 + 1e-6)]];
    if (c.kind === "rigorous") checks.push(["LB <= ref", lb <= ref * (1 + 1e-6)], ["UB >= ref", ub >= ref * (1 - 1e-6)]);
    if (c.kind === "lbOnly") checks.push(["LB <= ref", lb <= ref * (1 + 1e-6)]);
    const bad = checks.filter(([, ok]) => !ok).map(([n]) => n);
    if (bad.length) fails++;
    const ratio = (v) => (isFinite(ref) ? (v / ref).toFixed(3) : v.toFixed(1));
    console.log(`${bad.length ? "FAIL" : "ok  "} ${k.padEnd(7)} LB ${ratio(lb)}  UB ${ratio(ub)}  ${((Date.now() - t) / 1000).toFixed(0)} s  (${c.note})${bad.length ? "  failed: " + bad.join(", ") : ""}`);
  }
  console.log(fails ? `${fails} case(s) failed` : "All cases passed");
  process.exit(fails ? 1 : 0);
})();
