// Browser test of the built page (dist/bearing-bounds-standalone.html) with Playwright and Chromium.
// Usage: npm run build && node tests/e2e.js     (first time: npx playwright install chromium)
const { chromium } = require("playwright");
const path = require("path");
const fs = require("fs");
const { pathToFileURL } = require("url");

const page = pathToFileURL(path.join(__dirname, "..", "dist", "bearing-bounds-standalone.html")).href;
let fails = 0;
// a row of the solution detail table
const detail = (pg, label) => pg.$eval("#detail", (t, l) => { const r = [...t.rows].find((x) => x.cells[0].textContent === l); return r ? r.cells[1].textContent : ""; }, label);
const check = (name, ok, extra) => { if (!ok) fails++; console.log(`${ok ? "ok  " : "FAIL"} ${name}${extra ? "  (" + extra + ")" : ""}`); };

(async () => {
  const browser = await chromium.launch();
  const p = await browser.newPage({ viewport: { width: 1200, height: 1000 } });
  const errors = [];
  p.on("pageerror", (e) => errors.push(e.message));
  await p.goto(page);
  await p.waitForFunction(() => /Ready/.test(document.getElementById("status").textContent), null, { timeout: 60000 });
  const txt = (id) => p.$eval("#" + id, (e) => e.textContent.trim());
  const val = (id) => p.$eval("#" + id, (e) => e.value);

  // hidden fields really hidden
  check("Length L hidden for a square footing", (await p.$eval("#wrap-L", (e) => getComputedStyle(e).display)) === "none");
  // unit round trip leaves inputs unchanged
  const before = await Promise.all(["B", "D", "P", "wt"].map(val));
  await p.click("#u-si"); await p.click("#u-us");
  const after = await Promise.all(["B", "D", "P", "wt"].map(val));
  check("US -> SI -> US round trip", before.join() === after.join(), before.join() + " vs " + after.join());

  // run the default example
  await p.click("#go");
  await p.waitForFunction(() => /Finished|stopped|crashed/.test(document.getElementById("status").textContent), null, { timeout: 240000 });
  check("run finishes", /Finished/.test(await txt("status")), await txt("status"));
  const lb = parseFloat((await txt("fLB")).replace(/,/g, "")), ub = parseFloat((await txt("fUB")).replace(/,/g, ""));
  check("lower bound below upper bound", lb > 0 && lb <= ub, lb + " / " + ub);
  check("answer line shows the FS", /Factor of safety/.test(await txt("answer")));
  const lbSolver = await detail(p, "Lower-bound solver");
  check("lower bound solved with exact cones (Clarabel)", /^Clarabel/.test(lbSolver), lbSolver);
  check("animation tab opens after a run", (await p.$eval("[role=tab][aria-selected=true]", (e) => e.id)) === "tab-anim");

  // load-size change keeps results; geometry change makes them stale
  await p.fill("#P", "300");
  check("changing only the load keeps the results", !(await txt("verdict")).includes("re-run"));
  await p.fill("#B", "10");
  check("changing B marks results stale", (await txt("verdict")).includes("re-run"));
  await p.fill("#B", "6"); await p.fill("#P", "150");

  // footing type switch converts the load and keeps the applied pressure
  const qa = await txt("fApp");
  await p.click("#ftype button:nth-child(1)");
  check("strip switch keeps applied pressure", (await txt("fApp")) === qa, qa + " vs " + (await txt("fApp")));
  check("strip switch explains the conversion", (await p.$eval("#loadNote", (e) => !e.hidden)));
  await p.click("#ftype button:nth-child(2)");

  // inline validation
  await p.fill("#ec", "3"); await p.click("#go");
  check("invalid eccentricity is flagged on the field", (await p.$eval("#ec", (e) => e.getAttribute("aria-invalid"))) === "true");
  await p.fill("#ec", "0");

  // circular footing: axisymmetric option, and its conditions
  await p.click("#ftype button:nth-child(4)");
  check("axisymmetric option shown for a circular footing", !(await p.$eval("#wrap-axi", (e) => e.hidden)));
  await p.fill("#Hl", "10"); await p.click("#go");
  check("axisymmetric run with H is flagged", (await p.$eval("#axi", (e) => e.getAttribute("aria-invalid"))) === "true");
  await p.fill("#Hl", "0");
  await p.click("#ftype button:nth-child(1)");
  check("axisymmetric option hidden for a strip", await p.$eval("#wrap-axi", (e) => e.hidden));

  // sliding mode: V held, H increased
  await p.click("#lmode button:nth-child(3)");
  await p.click("#go");
  check("sliding without H is flagged", (await p.$eval("#Hl", (e) => e.getAttribute("aria-invalid"))) === "true");
  await p.fill("#Hl", "5"); await p.click("#go");
  await p.waitForFunction(() => /Finished|stopped|crashed/.test(document.getElementById("status").textContent), null, { timeout: 240000 });
  check("sliding run finishes", /Finished/.test(await txt("status")), await txt("status"));
  check("sliding answer is a factor of safety on H", /on H/.test(await txt("answer")), await txt("answer"));
  const hl = parseFloat((await txt("fLB")).replace(/,/g, "")), hu = parseFloat((await txt("fUB")).replace(/,/g, ""));
  check("sliding lower bound below upper bound", hl > 0 && hl <= hu, hl + " / " + hu);
  await p.click("#lmode button:nth-child(1)"); await p.fill("#Hl", "0");

  // the finite element upper bound, when asked for, also comes from Clarabel
  await p.$eval("#feub", (e) => { e.checked = true; e.dispatchEvent(new Event("change")); });
  await p.click("#go");
  await p.waitForFunction(() => /Finished|stopped|crashed/.test(document.getElementById("status").textContent), null, { timeout: 240000 });
  const ubSolver = await detail(p, "Finite element UB solver");
  check("finite element upper bound solved with exact cones (Clarabel)", /^Clarabel/.test(ubSolver), ubSolver);
  await p.$eval("#feub", (e) => { e.checked = false; e.dispatchEvent(new Event("change")); });

  // no horizontal scrolling at phone width
  await p.setViewportSize({ width: 390, height: 900 }); await p.waitForTimeout(300);
  check("no horizontal page scroll at 390 px", (await p.evaluate(() => document.documentElement.scrollWidth - document.documentElement.clientWidth)) <= 0);

  check("no script errors", errors.length === 0, errors.join("; "));

  // without the Clarabel module the page still works, with HiGHS alone
  const noCone = path.join(__dirname, "..", "dist", "e2e-no-clarabel.html");
  fs.writeFileSync(noCone, fs.readFileSync(path.join(__dirname, "..", "dist", "bearing-bounds-standalone.html"), "utf8")
    .replace(/(<script id="wasm-clarabel" type="application\/octet-stream">)[^<]*(<\/script>)/, "$1$2"));
  const p2 = await browser.newPage({ viewport: { width: 1200, height: 1000 } });
  const errors2 = [];
  p2.on("pageerror", (e) => errors2.push(e.message));
  await p2.goto(pathToFileURL(noCone).href);
  await p2.waitForFunction(() => /Ready/.test(document.getElementById("status").textContent), null, { timeout: 60000 });
  await p2.click("#go");
  await p2.waitForFunction(() => /Finished|stopped|crashed/.test(document.getElementById("status").textContent), null, { timeout: 240000 });
  const lb2 = await detail(p2, "Lower-bound solver");
  check("without Clarabel the run finishes with HiGHS", /Finished/.test(await p2.$eval("#status", (e) => e.textContent)) && /^HiGHS/.test(lb2), lb2);
  check("no script errors without Clarabel", errors2.length === 0, errors2.join("; "));
  fs.unlinkSync(noCone);
  await browser.close();
  console.log(fails ? `${fails} check(s) failed` : "All checks passed");
  process.exit(fails ? 1 : 0);
})();
