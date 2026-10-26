/* Plane-strain limit analysis of footings: rigorous lower and upper bounds.
 *
 * Lower bound (LB): finite element limit analysis (FELA) with linear stress triangles, stress
 *   discontinuities on every edge, an inscribed linearised Mohr-Coulomb yield surface, and an
 *   extension of the stress field to the infinite half-space beside and below the mesh.
 * Upper bound (UB): discontinuity layout optimisation (DLO), plus an optional finite element UB
 *   (linear velocity triangles, velocity discontinuities, circumscribed yield surface).
 * Footing: rigid body (sliding, settling, rotating) on a no-tension interface that is either rough
 *   (soil strength) or smooth. Loads follow a proportional path V, H = V tan(alpha), M = V e.
 * Ground: level, or a slope on the +x side. Layers are horizontal, measured from the top ground.
 * Coordinates: y up, footing base at y = 0, footing centre at x = 0; tension positive.
 * Effective-stress basis throughout; pore pressure under the footing is added to give total loads.
 * Units: any consistent set, e.g. (kPa, kN/m3, m) or (psf, pcf, ft).
 * When the load is vertical and central and the ground is level, a half model with symmetry is used.
 */
(function (root) {
  "use strict";

  /* ================= soil profile (depth below the top ground surface) ================= */
  function soilAt(model, depth) {
    const L = model.layers;
    for (let i = 0; i < L.length; i++) if (depth <= L[i].bottom + 1e-9) return L[i];
    return L[L.length - 1];
  }
  // effective-stress basis: buoyant weight below the water table for every layer. For undrained
  // (phi = 0) layers this only shifts normal stresses by the pore pressure, which leaves yield unchanged.
  function unitWeightAt(model, depth) {
    const s = soilAt(model, depth);
    const below = model.waterDepth != null && depth > model.waterDepth;
    return below ? s.gammaSat - model.gammaW : s.gamma;
  }
  function sigmaV(model, depth) {
    if (depth <= 0) return 0;
    const pts = [0, depth];
    model.layers.forEach((l) => { if (l.bottom > 0 && l.bottom < depth) pts.push(l.bottom); });
    if (model.waterDepth != null && model.waterDepth > 0 && model.waterDepth < depth) pts.push(model.waterDepth);
    pts.sort((a, b) => a - b);
    let s = 0;
    for (let i = 0; i < pts.length - 1; i++) s += unitWeightAt(model, 0.5 * (pts[i] + pts[i + 1])) * (pts[i + 1] - pts[i]);
    return s;
  }
  function porePressure(model, depth) {
    return model.waterDepth == null ? 0 : model.gammaW * Math.max(0, depth - model.waterDepth);
  }
  // strength of a layer at model y (depth below top ground = D - y)
  function cphiAt(model, layer, y) {
    const depth = model.D - y;
    return { c: layer.c + (layer.suGrad || 0) * Math.max(0, depth - layer.top), phi: layer.phi * Math.PI / 180 };
  }

  /* ================= geometry ================= */
  // Builds ground-surface functions for the current domain (xL..xR, depth H below the base).
  function geometry(m) {
    const b = m.B / 2, tol = 1e-9 * Math.max(1, m.B);
    const g = { b: b, sym: !!m.sym, xL: m.xL, xR: m.xR, H: m.H, slope: false, surfX: [] };
    const sl = m.slope;
    if (sl && sl.height > 1e-9 && sl.angle > 0) {
      const tb = Math.tan(sl.angle * Math.PI / 180);
      g.tb = tb; g.xcReal = b + sl.dist; g.HsReal = sl.height; g.xtReal = g.xcReal + sl.height / tb;
      if (sl.height > m.D + 1e-9) {        // the slope face cuts below the footing base
        g.slope = true;
        g.xc = g.xcReal + m.D / tb; g.Hs = sl.height - m.D; g.xt = g.xc + g.Hs / tb;
      }
    }
    // drop of the ground below the top surface at x
    g.drop = (x) => (g.tb === undefined ? 0 : Math.min(g.HsReal, Math.max(0, (x - g.xcReal) * g.tb)));
    g.groundY = (x) => -Math.max(0, g.drop(x) - m.D);
    g.slopeX = (y) => g.xc + (-y) / g.tb;               // face position at level y (slope case)
    g.rightX = (y) => (g.slope && y > -g.Hs + tol ? Math.min(m.xR, g.slopeX(y)) : m.xR);
    g.inside = (x, y) => x >= m.xL - tol && x <= m.xR + tol && y >= -m.H - tol && y <= g.groundY(x) + tol;
    // effective surcharge on the model surface (soil above the base level plus any user surcharge)
    g.qs = (x) => {
      if (Math.abs(x) <= b + tol) return 0;
      const dr = g.drop(x);
      if (dr >= m.D - tol && g.tb !== undefined && dr > 0) return 0;
      return sigmaV(m, m.D) - sigmaV(m, dr) + (dr <= tol ? (m.surcharge || 0) : 0);
    };
    // breakpoints along the surface where the surcharge changes slope
    if (g.tb !== undefined) {
      const xs = [g.xcReal];
      const ds = m.layers.map((l) => l.bottom).concat(m.waterDepth != null ? [m.waterDepth] : []).concat([m.D, g.HsReal]);
      ds.forEach((d) => { if (d > 0 && d <= Math.min(m.D, g.HsReal)) xs.push(g.xcReal + d / g.tb); });
      g.surfX = xs.filter((x) => x > b && x < m.xR);
    }
    // classify a boundary segment from its end points
    g.segType = (P, Q) => {
      const on = (f) => f(P) && f(Q);
      const top = (p) => Math.abs(p[1]) < tol && (!g.slope || p[0] <= g.xc + tol);
      const face = (p) => g.slope && p[1] <= tol && p[1] >= -g.Hs - tol && Math.abs(p[0] - g.slopeX(p[1])) < 1e-7 * Math.max(1, m.B);
      const toe = (p) => g.slope && Math.abs(p[1] + g.Hs) < tol && p[0] >= g.xt - tol;
      if (on(top)) return Math.max(Math.abs(P[0]), Math.abs(Q[0])) <= b + tol && Math.min(P[0], Q[0]) >= -b - tol ? "footing" : "surface";
      if (on(face) || on(toe)) return "surface";
      if (on((p) => Math.abs(p[0] - m.xL) < tol)) return g.sym ? "sym" : "left";
      if (on((p) => Math.abs(p[0] - m.xR) < tol)) return "right";
      if (on((p) => Math.abs(p[1] + m.H) < tol)) return "base";
      return "other";
    };
    return g;
  }

  /* ================= mesh ================= */
  function uniqSorted(arr, tol) {
    arr.sort((p, q) => p - q);
    const out = [];
    arr.forEach((v) => { if (!out.length || v - out[out.length - 1] > tol) out.push(v); });
    return out;
  }
  // horizontal lines that must be element edges: layer bottoms, water table, slope toe level
  function fixedDepths(model) {
    const extra = [];
    model.layers.forEach((l) => { const z = l.bottom - model.D; if (z > 1e-6 && z < model.H - 1e-6) extra.push(z); });
    if (model.waterDepth != null) { const z = model.waterDepth - model.D; if (z > 1e-6 && z < model.H - 1e-6) extra.push(z); }
    if (model.geo && model.geo.slope && model.geo.Hs < model.H - 1e-6) extra.push(model.geo.Hs);
    return uniqSorted(extra, 1e-6);
  }
  function grow(a, b, h0, r, hmax) {
    const pts = [a];
    let x = a, h = h0;
    while (x + h < b - 0.4 * h) { x += h; pts.push(x); h = Math.min(h * r, hmax); }
    pts.push(b);
    return pts;
  }
  function initialCells(model, nf) {
    const g = model.geo, b = g.b, h0 = b / nf, ext = model.ext;
    const side = (lim, reach) => {
      const hx = Math.max(h0, (reach - b) / 10);
      let xs = grow(b, Math.min(lim, Math.max(b + h0, reach * 1.05)), h0, 1.25, hx);
      if (xs[xs.length - 1] < lim - 1e-9) xs = xs.concat(grow(xs[xs.length - 1], lim, hx, 1.3, Infinity).slice(1));
      return xs;
    };
    let xs = [];
    for (let i = 0; i <= nf; i++) xs.push(b * i / nf);
    xs = xs.concat(side(model.xR, ext.xR).slice(1));
    if (!g.sym) {
      const left = side(-model.xL, -ext.xL).slice(1).map((x) => -x);
      for (let i = 1; i <= nf; i++) left.push(-b * i / nf);
      xs = xs.concat(left);
    }
    const key = [];
    if (g.slope) key.push(g.xc, g.xt);
    g.surfX.forEach((x) => key.push(x));
    xs = uniqSorted(xs.filter((x) => key.every((k) => Math.abs(x - k) > 0.3 * h0)).concat(key.filter((k) => k > model.xL && k < model.xR)), 1e-9);
    const hy = Math.max(h0, ext.depth / 7);
    let ds = grow(0, Math.min(model.H, ext.depth * 1.05), h0, 1.25, hy);
    if (ds[ds.length - 1] < model.H - 1e-9) ds = ds.concat(grow(ds[ds.length - 1], model.H, hy, 1.3, Infinity).slice(1));
    const extra = fixedDepths(model);
    if (extra.length) {
      ds = ds.filter((d) => d === 0 || d === model.H || extra.every((e) => Math.abs(e - d) > 0.3 * h0));
      ds = uniqSorted(ds.concat(extra), 1e-9);
    }
    const cells = [];
    for (let j = 0; j < ds.length - 1; j++) for (let i = 0; i < xs.length - 1; i++)
      cells.push({ x0: xs[i], x1: xs[i + 1], y0: -ds[j + 1], y1: -ds[j], lev: 0 });
    return cells;
  }
  function splitCell(c) {
    const xm = 0.5 * (c.x0 + c.x1), ym = 0.5 * (c.y0 + c.y1), l = c.lev + 1;
    return [
      { x0: c.x0, x1: xm, y0: c.y0, y1: ym }, { x0: xm, x1: c.x1, y0: c.y0, y1: ym },
      { x0: c.x0, x1: xm, y0: ym, y1: c.y1 }, { x0: xm, x1: c.x1, y0: ym, y1: c.y1 }].map((k) => Object.assign(k, { lev: l, ind: (c.ind || 0) / 4 }));
  }
  const key = (v) => Math.round(v * 1e7);
  function cellNeighbours(cells) {
    const V = new Map(), Hm = new Map();
    const push = (m, k, item) => { if (!m.has(k)) m.set(k, []); m.get(k).push(item); };
    cells.forEach((c, i) => {
      push(V, key(c.x0), { i: i, a: c.y0, b: c.y1 }); push(V, key(c.x1), { i: i, a: c.y0, b: c.y1 });
      push(Hm, key(c.y0), { i: i, a: c.x0, b: c.x1 }); push(Hm, key(c.y1), { i: i, a: c.x0, b: c.x1 });
    });
    const nb = cells.map(() => new Set());
    [V, Hm].forEach((m) => m.forEach((list) => {
      list.sort((p, q) => p.a - q.a);
      for (let s = 0; s < list.length; s++) for (let t = s + 1; t < list.length; t++) {
        if (list[t].a >= list[s].b - 1e-12) break;
        if (list[s].i === list[t].i) continue;
        const ov = Math.min(list[s].b, list[t].b) - Math.max(list[s].a, list[t].a);
        if (ov > 1e-12) { nb[list[s].i].add(list[t].i); nb[list[t].i].add(list[s].i); }
      }
    }));
    return nb;
  }
  function refineCells(cells, marked) {
    let cur = [];
    cells.forEach((c, i) => { if (marked.has(i)) cur.push(...splitCell(c)); else cur.push(c); });
    for (let it = 0; it < 30; it++) { // 2:1 balance
      const nb = cellNeighbours(cur);
      const bad = new Set();
      cur.forEach((c, i) => nb[i].forEach((j) => { if (cur[j].lev > c.lev + 1) bad.add(i); }));
      if (!bad.size) break;
      const nxt = [];
      cur.forEach((c, i) => { if (bad.has(i)) nxt.push(...splitCell(c)); else nxt.push(c); });
      cur = nxt;
    }
    return cur;
  }

  function delaunay(pts, areaScale) {
    let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
    pts.forEach(([x, y]) => { minx = Math.min(minx, x); miny = Math.min(miny, y); maxx = Math.max(maxx, x); maxy = Math.max(maxy, y); });
    const d = Math.max(maxx - minx, maxy - miny) * 2000, mx = (minx + maxx) / 2, my = (miny + maxy) / 2;
    const P = pts.concat([[mx - d, my - d], [mx + d, my - d], [mx, my + d]]);
    const n = pts.length;
    const mk = (a, b, c) => {
      const [ax, ay] = P[a], [bx, by] = P[b], [cx, cy] = P[c];
      const D = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
      const ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / D;
      const uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / D;
      return { v: [a, b, c], x: ux, y: uy, r2: (ax - ux) * (ax - ux) + (ay - uy) * (ay - uy) };
    };
    let tris = [mk(n, n + 1, n + 2)];
    const order = pts.map((_, i) => i).sort((i, j) => (pts[i][0] - pts[j][0]) || (pts[i][1] - pts[j][1]));
    for (const i of order) {
      const px = P[i][0], py = P[i][1];
      const keep = [], edges = new Map();
      for (const t of tris) {
        const dx = px - t.x, dy = py - t.y;
        if (dx * dx + dy * dy < t.r2 * (1 - 1e-11)) {
          for (let k = 0; k < 3; k++) {
            const a = t.v[k], b = t.v[(k + 1) % 3], kk = a < b ? a * 1e6 + b : b * 1e6 + a;
            const ex = edges.get(kk);
            if (ex) ex.c++; else edges.set(kk, { a: a, b: b, c: 1 });
          }
        } else keep.push(t);
      }
      tris = keep;
      edges.forEach((e) => { if (e.c === 1) tris.push(mk(e.a, e.b, i)); });
    }
    const out = [];
    tris.forEach((t) => {
      if (t.v[0] >= n || t.v[1] >= n || t.v[2] >= n) return;
      const [a, b, c] = t.v;
      const A2 = (P[b][0] - P[a][0]) * (P[c][1] - P[a][1]) - (P[c][0] - P[a][0]) * (P[b][1] - P[a][1]);
      if (Math.abs(A2) < 1e-12 * (areaScale || 1)) return;
      out.push(A2 > 0 ? [a, b, c] : [a, c, b]);
    });
    return out;
  }
  // polar fan of points around a singular point (cx, cy): rings growing geometrically from r0 to R, each
  // swept round the full circle and clipped to the domain, plus the points where each ring meets the slope
  // face or the toe level (so the fan conforms to those boundaries)
  function fanPoints(model, fan, f) {
    const g = model.geo, pts = [[f.cx, f.cy]];
    const nth = f.nTheta || fan.nTheta, dth = Math.PI / nth, R = f.R, tolI = 1e-9;
    let r = f.r0 || fan.r0;
    const radii = [];
    while (r < R) { radii.push(r); r *= 1 + dth; }
    radii.push(R);
    const ok = (x, y) => g.inside(x, y) && x > model.xL + tolI && x < model.xR - tolI && y > -model.H + tolI;
    radii.forEach((rr) => {
      for (let j = 0; j < 2 * nth; j++) {
        const th = j * dth, x = f.cx + rr * Math.cos(th), y = f.cy - rr * Math.sin(th);
        if (ok(x, y)) pts.push([x, y]);
      }
      if (g.slope && f.kind !== "footing") {
        const ca = Math.cos(Math.atan(g.tb)), sa = Math.sin(Math.atan(g.tb));
        [[f.cx + rr * ca, f.cy - rr * sa], [f.cx - rr * ca, f.cy + rr * sa], [f.cx + rr, f.cy], [f.cx - rr, f.cy]].forEach(([x, y]) => {
          const onF = y <= 1e-12 && y >= -g.Hs - 1e-12 && Math.abs(x - g.slopeX(y)) < 1e-9 * Math.max(1, model.B);
          const onLevel = (Math.abs(y) < 1e-12 && x <= g.xc + 1e-12) || (Math.abs(y + g.Hs) < 1e-12 && x >= g.xt - 1e-12);
          if ((onF || onLevel) && ok(x, y)) pts.push([x, y]);
        });
      }
      if (g.sym && f.kind === "footing" && rr > Math.abs(f.cx) + 1e-9) { const y = -Math.sqrt(rr * rr - f.cx * f.cx); if (y > -model.H + 1e-9) pts.push([model.xL, y]); }
    });
    return pts;
  }
  // the fans of a model: footing edges, and the slope crest and toe when they are clear of the footing
  function fansOf(model, opts) {
    const g = model.geo, b = g.b, nTheta = opts.nTheta || 16, r0 = b * (opts.fanR0 || 0.04);
    const gapC = g.slope ? g.xc - b : Infinity, crest = g.slope && gapC > 0.3 * b && opts.slopeFans === true;
    const R = g.slope ? Math.max(0.16 * b, Math.min(b, (crest ? 0.65 : 0.8) * gapC)) : b;
    const out = (g.sym ? [b] : [-b, b]).map((cx) => ({ kind: "footing", cx: cx, cy: 0, R: R, nTheta: nTheta, r0: r0 }));
    if (crest) {
      const Rc = Math.min(b, 0.3 * gapC, 0.45 * g.Hs);
      out.push({ kind: "crest", cx: g.xc, cy: 0, R: Rc, nTheta: Math.max(8, Math.round(nTheta * 0.75)), r0: Rc * 0.05 });
      const Rt = Math.min(b, 0.45 * g.Hs, 0.45 * (model.xR - g.xt));
      if (g.xt < model.xR - 0.2 * b && Rt > 0.05 * b) out.push({ kind: "toe", cx: g.xt, cy: -g.Hs, R: Rt, nTheta: Math.max(8, Math.round(nTheta * 0.75)), r0: Rt * 0.05 });
    }
    return out;
  }
  function geoOf(XY) {
    return XY.map((X) => {
      const A2 = (X[2] - X[0]) * (X[5] - X[1]) - (X[4] - X[0]) * (X[3] - X[1]);
      return { A: A2 / 2, b: [(X[3] - X[5]) / A2, (X[5] - X[1]) / A2, (X[1] - X[3]) / A2],
        c: [(X[4] - X[2]) / A2, (X[0] - X[4]) / A2, (X[2] - X[0]) / A2], cx: (X[0] + X[2] + X[4]) / 3, cy: (X[1] + X[3] + X[5]) / 3 };
    });
  }
  function domainArea(model) {
    const g = model.geo;
    let a = (model.xR - model.xL) * model.H;
    if (g.slope) {
      const x0 = Math.max(model.xL, g.xc), x1 = Math.min(model.xR, g.xt);
      if (x1 > x0) a -= ((-g.groundY(x0)) + (-g.groundY(x1))) / 2 * (x1 - x0);
      if (model.xR > g.xt) a -= g.Hs * (model.xR - Math.max(g.xt, model.xL));
    }
    return a;
  }

  // Unstructured mesh: polar fans at the footing edges, quadtree-sized background points, every layer
  // boundary, the water table, the slope face and the toe level as element edges. Each horizontal band is
  // convex (a trapezoid where it meets the slope face) and is triangulated on its own.
  function buildTriMesh(model, cells, fan) {
    const g = model.geo, b = g.b, W = model.xR - model.xL, H = model.H;
    const tol0 = 1e-9 * Math.max(W, H);
    const pts = [], seen = new Set();
    const add = (x, y) => { const k = key(x) + ":" + key(y); if (!seen.has(k)) { seen.add(k); pts.push([x, y]); } };
    const fanList = fan.list || (g.sym ? [b] : [-b, b]).map((cx) => ({ kind: "footing", cx: cx, cy: 0, R: fan.R }));
    fanList.forEach((f) => fanPoints(model, fan, f).forEach(([x, y]) => add(x, y)));
    const nearFan = (x, y, h) => fanList.some((f) => Math.hypot(x - f.cx, y - f.cy) < f.R + 0.6 * h);
    const onFace = (x, y) => g.slope && y > -g.Hs - tol0 && y < tol0 && Math.abs(x - g.slopeX(y)) < 1e-7 * Math.max(W, H);
    cells.forEach((c) => {
      const h = Math.max(c.x1 - c.x0, c.y1 - c.y0);
      [[c.x0, c.y0], [c.x1, c.y0], [c.x0, c.y1], [c.x1, c.y1]].forEach(([x, y]) => {
        if (nearFan(x, y, h) || !g.inside(x, y)) return;
        add(x, y);
      });
      // slope face points where cell edges meet the face
      if (g.slope) {
        [c.y0, c.y1].forEach((y) => { if (y <= tol0 && y >= -g.Hs - tol0) { const x = g.slopeX(y); if (x >= c.x0 - tol0 && x <= c.x1 + tol0 && x <= model.xR + tol0) add(x, y); } });
        [c.x0, c.x1].forEach((x) => { if (x > g.xc && x < g.xt && x <= model.xR) { const y = g.groundY(x); if (y >= c.y0 - tol0 && y <= c.y1 + tol0) add(x, y); } });
      }
    });
    g.surfX.forEach((x) => { if (x < model.xR) add(x, g.groundY(x)); });
    // band lines and their end points
    const lines = fixedDepths(model).map((z) => -z);
    [0].concat(lines).concat([-H]).forEach((y) => { add(model.xL, y); add(g.rightX(y), y); });
    if (g.slope && g.xc < model.xR) add(g.xc, 0);
    if (g.slope && g.xt < model.xR) add(g.xt, -g.Hs);
    if (g.slope && g.xt >= model.xR && g.xc < model.xR) add(model.xR, g.groundY(model.xR));
    // move points that crowd an interface line or the slope face onto it
    const hp = pts.map((p, i) => {
      let d = Infinity;
      for (let j = 0; j < pts.length; j++) if (j !== i) { const dd = Math.hypot(pts[j][0] - p[0], pts[j][1] - p[1]); if (dd < d) d = dd; }
      return d;
    });
    const keepPt = pts.map(() => true);
    const projections = lines.map((y0) => {
      const set = [];
      pts.forEach((p, i) => {
        const dy = Math.abs(p[1] - y0);
        if (dy > tol0 && dy < 0.55 * hp[i] && p[0] < g.rightX(y0) - tol0 && p[0] > model.xL + tol0) { keepPt[i] = false; set.push([p[0], hp[i]]); }
      });
      return { y0: y0, xs: set };
    });
    const faceCand = [];
    if (g.slope) {
      const t = Math.atan(g.tb), ux = Math.cos(t), uy = -Math.sin(t);   // unit vector down the face
      pts.forEach((p, i) => {
        if (!keepPt[i] || onFace(p[0], p[1]) || p[1] > -tol0 || p[1] < -g.Hs + tol0) return;
        const dx = g.slopeX(p[1]) - p[0];
        if (dx <= tol0) return;
        const dist = dx * Math.sin(t);                                    // perpendicular distance to the face
        if (dist < 0.55 * hp[i]) {
          keepPt[i] = false;
          const s = (p[0] - g.xc) * ux + p[1] * uy;                       // foot of the perpendicular
          faceCand.push([Math.max(-g.Hs, Math.min(0, s * uy)), hp[i]]);
        }
      });
    }
    const P2 = pts.filter((_, i) => keepPt[i]);
    const seen2 = new Set(P2.map((p) => key(p[0]) + ":" + key(p[1])));
    const add2 = (x, y) => { const k = key(x) + ":" + key(y); if (!seen2.has(k)) { seen2.add(k); P2.push([x, y]); } };
    const thin = (fixed, cand, scale) => {
      const all = fixed.map((v) => [v, Infinity, true]).concat(cand.map(([v, h]) => [v, h, false])).sort((p, q) => p[0] - q[0]);
      const out = [];
      all.forEach((e) => {
        if (!out.length) { out.push(e); return; }
        const last = out[out.length - 1];
        if ((e[0] - last[0]) * scale < 0.4 * Math.min(e[1], last[1])) { if (e[2] && !last[2]) out[out.length - 1] = e; return; }
        out.push(e);
      });
      return out.filter((e) => !e[2]).map((e) => e[0]);
    };
    projections.forEach(({ y0, xs }) => {
      const fixed = P2.filter((p) => Math.abs(p[1] - y0) < tol0).map((p) => p[0]);
      thin(fixed, xs, 1).forEach((x) => add2(x, y0));
    });
    if (g.slope) {
      const fixed = P2.filter((p) => onFace(p[0], p[1])).map((p) => p[1]);
      thin(fixed, faceCand, 1 / Math.sin(Math.atan(g.tb))).forEach((y) => add2(g.slopeX(y), y));
    }
    const bandsY = [0].concat(lines.slice().sort((a, c) => c - a)).concat([-H]);
    const tri = [];
    for (let k = 0; k < bandsY.length - 1; k++) {
      const yt = bandsY[k], yb = bandsY[k + 1];
      const ids = [];
      const aboveToe = g.slope && yb > -g.Hs - tol0;  // band bounded on the right by the slope face
      const bandRight = (y) => (aboveToe ? Math.min(model.xR, g.slopeX(y)) : model.xR);
      P2.forEach((p, i) => { if (p[1] <= yt + tol0 && p[1] >= yb - tol0 && p[0] <= bandRight(p[1]) + 1e-7 * W) ids.push(i); });
      delaunay(ids.map((i) => P2[i]), W * Math.abs(yt - yb)).forEach((tt) => tri.push(tt.map((j) => ids[j])));
    }
    const XY = tri.map(([a, c, d]) => [P2[a][0], P2[a][1], P2[c][0], P2[c][1], P2[d][0], P2[d][1]]);
    const geo = geoOf(XY);
    const edgeMap = new Map(), ifaces = [], bnd = [];
    tri.forEach((t, e) => { for (let k = 0; k < 3; k++) {
      const a = t[k], c = t[(k + 1) % 3], kk = a < c ? a + "," + c : c + "," + a;
      const ex = edgeMap.get(kk);
      if (ex) {
        ifaces.push({ P: P2[ex.ga], Q: P2[ex.gb], s1: { e: ex.e, i: ex.i, j: ex.j, tP: 0, tQ: 1 }, s2: { e: e, i: t.indexOf(ex.ga), j: t.indexOf(ex.gb), tP: 0, tQ: 1 } });
        edgeMap.delete(kk);
      } else edgeMap.set(kk, { e: e, i: k, j: (k + 1) % 3, ga: a, gb: c });
    } });
    let bad = 0;
    edgeMap.forEach((ex) => {
      const P = P2[ex.ga], Q = P2[ex.gb], type = g.segType(P, Q);
      if (type === "other") bad++;
      bnd.push({ e: ex.e, i: ex.i, j: ex.j, P: P, Q: Q, type: type });
    });
    orientIfaces(ifaces, geo);
    const cellOf = new Int32Array(XY.length).fill(-1);
    geo.forEach((gg, e) => {
      for (let ci = 0; ci < cells.length; ci++) {
        const c = cells[ci];
        if (gg.cx >= c.x0 && gg.cx <= c.x1 && gg.cy >= c.y0 && gg.cy <= c.y1) { cellOf[e] = ci; break; }
      }
    });
    let area = 0; geo.forEach((gg) => (area += gg.A));
    const target = domainArea(model);
    if (bad || Math.abs(area - target) > 1e-7 * target) throw new Error("Mesh generation failed (area " + area.toFixed(4) + " vs " + target.toFixed(4) + ", " + bad + " stray edges)");
    return { cells: cells, XY: XY, geo: geo, ifaces: ifaces, bnd: bnd, nE: XY.length, cellOf: cellOf, area: area, fan: fan, conforming: true };
  }
  function orientIfaces(ifaces, geo) {
    ifaces.forEach((f) => {
      const L = Math.hypot(f.Q[0] - f.P[0], f.Q[1] - f.P[1]);
      const tx = (f.Q[0] - f.P[0]) / L, ty = (f.Q[1] - f.P[1]) / L;
      let nx = ty, ny = -tx;
      const g1 = geo[f.s1.e], g2 = geo[f.s2.e];
      if ((g2.cx - g1.cx) * nx + (g2.cy - g1.cy) * ny < 0) { nx = -nx; ny = -ny; }
      f.L = L; f.t = [tx, ty]; f.n = [nx, ny];
    });
  }
  // Quadtree mesh (4 crossed triangles per cell, hanging nodes handled by the discontinuities).
  // Used only as a fallback for level ground when the unstructured mesher fails.
  function buildQuadMesh(model, cells) {
    const g = model.geo;
    if (g.slope) throw new Error("Mesh generation failed for this slope geometry");
    const XY = [], sides = [], inner = [];
    cells.forEach((c) => {
      const BL = [c.x0, c.y0], BR = [c.x1, c.y0], TR = [c.x1, c.y1], TL = [c.x0, c.y1];
      const M = [0.5 * (c.x0 + c.x1), 0.5 * (c.y0 + c.y1)];
      const e0 = XY.length;
      [[BL, BR], [BR, TR], [TR, TL], [TL, BL]].forEach(([P, Q]) => XY.push([P[0], P[1], Q[0], Q[1], M[0], M[1]]));
      sides.push({ e: e0, horiz: true, coord: c.y0, a: c.x0, b: c.x1, P: BL, Q: BR });
      sides.push({ e: e0 + 1, horiz: false, coord: c.x1, a: c.y0, b: c.y1, P: BR, Q: TR });
      sides.push({ e: e0 + 2, horiz: true, coord: c.y1, a: c.x0, b: c.x1, P: TR, Q: TL });
      sides.push({ e: e0 + 3, horiz: false, coord: c.x0, a: c.y0, b: c.y1, P: TL, Q: BL });
      for (let k = 0; k < 4; k++) {
        const ea = e0 + k, eb = e0 + (k + 1) % 4;
        inner.push({ P: XY[ea].slice(2, 4), Q: [M[0], M[1]], s1: { e: ea, i: 1, j: 2, tP: 0, tQ: 1 }, s2: { e: eb, i: 0, j: 2, tP: 0, tQ: 1 } });
      }
    });
    const geo = geoOf(XY);
    const lines = new Map();
    sides.forEach((s, k) => { const kk = (s.horiz ? "h" : "v") + key(s.coord); if (!lines.has(kk)) lines.set(kk, []); lines.get(kk).push(k); });
    const ifaces = inner.slice(), covered = sides.map(() => 0);
    const param = (s, v) => { const pa = s.horiz ? s.P[0] : s.P[1], qa = s.horiz ? s.Q[0] : s.Q[1]; return (v - pa) / (qa - pa); };
    lines.forEach((list) => {
      list.sort((p, q) => sides[p].a - sides[q].a);
      for (let u = 0; u < list.length; u++) for (let w = u + 1; w < list.length; w++) {
        const A = sides[list[u]], Bs = sides[list[w]];
        if (Bs.a >= A.b - 1e-12) break;
        const lo = Math.max(A.a, Bs.a), hi = Math.min(A.b, Bs.b);
        if (hi - lo <= 1e-12) continue;
        covered[list[u]] += hi - lo; covered[list[w]] += hi - lo;
        ifaces.push({ P: A.horiz ? [lo, A.coord] : [A.coord, lo], Q: A.horiz ? [hi, A.coord] : [A.coord, hi],
          s1: { e: A.e, i: 0, j: 1, tP: param(A, lo), tQ: param(A, hi) }, s2: { e: Bs.e, i: 0, j: 1, tP: param(Bs, lo), tQ: param(Bs, hi) } });
      }
    });
    const bnd = [];
    sides.forEach((s, k) => { if (covered[k] <= (s.b - s.a) * 0.5) bnd.push({ e: s.e, i: 0, j: 1, P: s.P, Q: s.Q, type: g.segType(s.P, s.Q) }); });
    orientIfaces(ifaces, geo);
    const cellOf = new Int32Array(XY.length); for (let e = 0; e < XY.length; e++) cellOf[e] = Math.floor(e / 4);
    return { cells: cells, XY: XY, geo: geo, ifaces: ifaces, bnd: bnd, nE: XY.length, cellOf: cellOf };
  }
  function assignMaterials(model, mesh) {
    return mesh.geo.map((g) => {
      const depth = model.D - g.cy;
      return { layer: soilAt(model, depth), gamma: unitWeightAt(model, depth) };
    });
  }

  /* ================= LP ================= */
  function LP() { this.nc = 0; this.cost = []; this.lo = []; this.up = []; this.ri = []; this.rv = []; this.rs = [0]; this.rlo = []; this.rup = []; this.soc = []; }
  // second-order cone: f0 >= ||(f1, f2)|| for affine forms f = [idx, val, constant] (cone mode only). A
  // check-only cone is not sent to the solver (an equivalent linear row is) but is checked like the others.
  LP.prototype.addSOC = function (forms, checkOnly) {
    this.soc.push(forms.map(([idx, val, c0]) => {
      const m = new Map();
      for (let i = 0; i < idx.length; i++) if (val[i] !== 0) m.set(idx[i], (m.get(idx[i]) || 0) + val[i]);
      const I = [], V = [];
      m.forEach((v, k) => { if (Math.abs(v) > 1e-14) { I.push(k); V.push(v); } });
      return { idx: I, val: V, c: c0 || 0 };
    }));
    if (checkOnly) (this.socCheckOnly || (this.socCheckOnly = new Set())).add(this.soc.length - 1);
    return this.soc.length - 1;
  };
  LP.prototype.addVar = function (lo, up, cost) { this.lo.push(lo); this.up.push(up); this.cost.push(cost || 0); return this.nc++; };
  LP.prototype.addRow = function (idx, val, lo, up) {
    const m = new Map();
    for (let i = 0; i < idx.length; i++) if (val[i] !== 0) m.set(idx[i], (m.get(idx[i]) || 0) + val[i]);
    m.forEach((v, k) => { if (Math.abs(v) > 1e-14) { this.ri.push(k); this.rv.push(v); } });
    this.rs.push(this.ri.length); this.rlo.push(lo); this.rup.push(up);
  };
  LP.prototype.toModel = function (highs, sense) {
    const inf = highs.infinity;
    const fix = (a) => Float64Array.from(a, (v) => (v === Infinity ? inf : v === -Infinity ? -inf : v));
    const nr = this.rlo.length;
    return {
      numCols: this.nc, numRows: nr,
      sense: sense === "max" ? highs.constants.objectiveSense.maximize : highs.constants.objectiveSense.minimize,
      colCost: Float64Array.from(this.cost), colLower: fix(this.lo), colUpper: fix(this.up),
      rowLower: fix(this.rlo), rowUpper: fix(this.rup),
      matrix: { format: "csr", numRows: nr, numCols: this.nc, starts: Int32Array.from(this.rs), indices: Int32Array.from(this.ri), values: Float64Array.from(this.rv) }
    };
  };
  // A bound needs a feasible point, not a proven optimum: any statically admissible stress field is a
  // lower bound and any admissible mechanism an upper bound. Accept a feasible primal solution; if the
  // fast interior-point run (no crossover) does not give one, re-solve with crossover.
  function runLP(highs, lp, sense, opts) {
    const data = lp.toModel(highs, sense);
    const attempt = (crossover, solver) => {
      const model = highs.createModel(data);
      try {
        const o = { log_to_console: !!globalThis.__LPLOG, output_flag: !!globalThis.__LPLOG, presolve: "on", solver: solver || opts.lpSolver || "ipm", run_crossover: crossover, time_limit: opts.timeLimit || 600 };
        if (opts.ipmTol) o.ipm_optimality_tolerance = opts.ipmTol;
        model.options.set(o);
        try { model.run(); } catch (e) { if (globalThis.__LPDEBUG) console.log("run error", solver, crossover, e.message); return { status: -1, pss: -1, pinf: NaN, obj: NaN, x: null }; }
        const status = model.getModelStatus();
        let pss = -1; try { pss = Number(model.info.get("primal_solution_status")); } catch (e) { /* not reported */ }
        let pinf = NaN; try { pinf = Number(model.info.get("max_primal_infeasibility")); } catch (e) { /* not reported */ }
        const sol = model.getSolution();
        return { status: status, pss: pss, pinf: pinf, obj: model.getObjectiveValue(), x: sol.colValue, y: sol.rowDual };
      } finally { model.dispose(); }
    };
    const good = (q) => q.pss === 2 && (q.status === 7 || q.status === 15) && !(q.pinf > 1e-6);
    // A feasible point that is not proven optimal is still a valid bound, but the interior-point solver can stop
    // early (status 15); its point is usually as good as the optimum, but a lower bound of zero or less is
    // suspect. Retry with crossover in that case (or when no valid point came back) and keep the better one.
    const better = (q, r0) => good(q) && (!good(r0) || (sense === "max" ? q.obj > r0.obj : q.obj < r0.obj));
    let r = attempt(opts.crossover || "off"), crossover = false;
    if (r.status !== 8 && r.status !== 9 && r.status !== 10 && (!good(r) || (r.status !== 7 && sense === "max" && !(r.obj > 0)))) {
      const r2 = attempt("on"); crossover = true;
      if (better(r2, r) || (!good(r) && r2.status !== -1)) r = r2;
    }
    // interior point can fail outright on badly scaled problems; dual simplex is slower but robust
    // (simplex only when interior point failed outright: on these large LPs it is very slow)
    if (r.status === -1) r = attempt("off", "simplex");
    return { status: r.status, optimal: r.status === 7, valid: good(r), crossover: crossover, pinf: r.pinf, obj: r.obj, x: good(r) ? r.x : null, y: good(r) ? r.y : null, nVars: data.numCols, nRows: data.numRows };
  }
  function statusText(s) { return { 7: "optimal", 8: "infeasible", 9: "unbounded or infeasible", 10: "unbounded", 13: "time limit", 14: "iteration limit", 15: "not proven optimal" }[s] || "status " + s; }

  /* ================= exact second-order cones: Clarabel ================= */
  // In cone mode the finite element solvers state the Mohr-Coulomb condition exactly as second-order cones
  // (lp.soc) and the problem goes to the Clarabel interior-point solver (WebAssembly, registered with
  // setConic). Clarabel's point is accepted only after a check of every constraint of the problem at that
  // point, with the 1e-6 absolute rule of runLP: equality and inequality rows, variable bounds and every
  // cone. The objective is recomputed from the point, so the bound is that of a field that passed the check.
  // Anything else falls back to the linearised problem and HiGHS (conicFirst).
  let CONIC = null;
  function setConic(c) { CONIC = c || null; }
  const CLARABEL_STATUS = ["unsolved", "solved", "primal infeasible", "dual infeasible", "almost solved", "almost primal infeasible",
    "almost dual infeasible", "iteration limit", "time limit", "numerical error", "insufficient progress", "stopped"];
  const FEAS_TOL = 1e-6;
  function runConic(lp, sense, opts) {
    const n0 = lp.nc, R = lp.rlo.length, INF = Infinity;
    const fail = (why, extra) => Object.assign({ status: -1, optimal: false, valid: false, obj: NaN, x: null, y: null, nVars: n0, nRows: R, nCones: lp.soc.length, solver: "Clarabel", reason: why }, extra || {});
    // eliminate fixed variables (lo = up), which also turns the cones at fixed stress points into constants
    const fixed = new Float64Array(n0), isFix = new Uint8Array(n0), col = new Int32Array(n0);
    let n = 0;
    for (let j = 0; j < n0; j++) {
      if (lp.lo[j] === lp.up[j]) { isFix[j] = 1; fixed[j] = lp.lo[j]; col[j] = -1; } else col[j] = n++;
    }
    // constraint rows of Clarabel's A x + s = b: zero cone, then non-negative cone, then the second-order cones
    const zr = [], nr = [];   // [kind (0 row, 1 variable), index, sign, rhs]
    for (let r = 0; r < R; r++) {
      let fc = 0, live = 0;
      for (let k = lp.rs[r]; k < lp.rs[r + 1]; k++) { const j = lp.ri[k]; if (isFix[j]) fc += lp.rv[k] * fixed[j]; else live++; }
      const lo = lp.rlo[r] - fc, up = lp.rup[r] - fc;
      if (!live) { if (lo > FEAS_TOL || up < -FEAS_TOL) return fail("a constraint on fixed values only is violated"); continue; }
      if (lp.rlo[r] === lp.rup[r]) zr.push([0, r, 1, up]);
      else { if (up < INF) nr.push([0, r, 1, up]); if (lo > -INF) nr.push([0, r, -1, -lo]); }
    }
    for (let j = 0; j < n0; j++) {
      if (isFix[j]) continue;
      if (lp.up[j] < INF) nr.push([1, j, 1, lp.up[j]]);
      if (lp.lo[j] > -INF) nr.push([1, j, -1, -lp.lo[j]]);
    }
    const socKeep = [];
    let constBad = false;
    lp.soc.forEach((F, s) => {
      if (lp.socCheckOnly && lp.socCheckOnly.has(s)) return;
      const f = F.map((fm) => { let c = fm.c, live = 0; fm.idx.forEach((j, k) => { if (isFix[j]) c += fm.val[k] * fixed[j]; else live++; }); return { c: c, live: live }; });
      // a cone on fixed values only (the apex at a fixed zero stress point, say) is checked and dropped
      if (f.every((q) => !q.live)) { if (Math.hypot(f[1].c, f[2].c) - f[0].c > FEAS_TOL) constBad = true; return; }
      socKeep.push(s);
    });
    if (constBad) return fail("a cone on fixed values only is violated");
    const m = zr.length + nr.length + 3 * socKeep.length;
    // Clarabel solves with every cone narrowed by a relative 1e-6 (f0 (1 - 1e-6) >= ||(f1, f2)||), so that its
    // point, which is only feasible to within its tolerance, lies inside the exact cone that is checked below
    const shrink = 1 - (opts.coneShrink == null ? 1e-6 : opts.coneShrink);
    // A in compressed sparse column form (rows emitted in order, so each column's rows are sorted)
    const cnt = new Int32Array(n + 1);
    const each = (fn) => {
      let i = 0;
      const rowOf = (q) => {
        if (q[0] === 0) { for (let k = lp.rs[q[1]]; k < lp.rs[q[1] + 1]; k++) { const j = lp.ri[k]; if (!isFix[j]) fn(i, col[j], q[2] * lp.rv[k]); } }
        else fn(i, col[q[1]], q[2]);
        i++;
      };
      zr.forEach(rowOf); nr.forEach(rowOf);
      socKeep.forEach((s) => lp.soc[s].forEach((fm, u) => { const f = u ? 1 : shrink; fm.idx.forEach((j, k) => { if (!isFix[j]) fn(i, col[j], -f * fm.val[k]); }); i++; }));
    };
    each((i, j) => cnt[j + 1]++);
    for (let j = 0; j < n; j++) cnt[j + 1] += cnt[j];
    const nnz = cnt[n], rowIdx = new Uint32Array(nnz), aVal = new Float64Array(nnz), pos = cnt.slice(0, n);
    each((i, j, v) => { rowIdx[pos[j]] = i; aVal[pos[j]++] = v; });
    const b = new Float64Array(m);
    let i0 = 0;
    zr.forEach((q) => (b[i0++] = q[3])); nr.forEach((q) => (b[i0++] = q[3]));
    socKeep.forEach((s) => lp.soc[s].forEach((fm, u) => { let c = fm.c; fm.idx.forEach((j, k) => { if (isFix[j]) c += fm.val[k] * fixed[j]; }); b[i0++] = (u ? 1 : shrink) * c; }));
    const q = new Float64Array(n), sg = sense === "max" ? -1 : 1;
    for (let j = 0; j < n0; j++) if (!isFix[j]) q[col[j]] = sg * lp.cost[j];
    const tol = opts.coneTol || 1e-7;
    const colPtr = Uint32Array.from(cnt), socDims = Uint32Array.from(socKeep, () => 3);
    let iters = 0;
    // one Clarabel solve and the check of its point
    const attempt = (reg) => {
      // settings: max_iter, time_limit, tol_gap_abs, tol_gap_rel, tol_feas, tol_ktratio, verbose, equilibrate,
      // presolve, static regularisation (see wasm/clarabel/src/lib.rs)
      const settings = Float64Array.from([opts.coneMaxIter || 100, opts.timeLimit || 600, tol, tol, opts.coneFeasTol || 1e-9, 1e-7, globalThis.__LPLOG ? 1 : 0, 1, 1, reg]);
      let out;
      try { out = CONIC.solve(n, m, q, colPtr, rowIdx, aVal, b, zr.length, nr.length, socDims, settings); }
      catch (e) { return fail("Clarabel failed: " + (e && e.message ? e.message : e)); }
      if (!out || out.length < 5) return fail("Clarabel rejected the problem data");
      const status = out[0];
      iters += out[1];
      const info = { coneStatus: CLARABEL_STATUS[status] || "status " + status, iterations: iters, nCones: lp.soc.length };
      if (status !== 1 && status !== 4) return fail("Clarabel stopped: " + info.coneStatus, info);
      // the point in the original variables
      const x = new Float64Array(n0);
      for (let j = 0; j < n0; j++) x[j] = isFix[j] ? fixed[j] : Math.min(lp.up[j], Math.max(lp.lo[j], out[5 + col[j]]));
      // dual of the first component of every cone (a plastic multiplier; 0 for dropped cones)
      const zs = new Float64Array(lp.soc.length), zBase = 5 + n + zr.length + nr.length;
      socKeep.forEach((s, k) => (zs[s] = out[zBase + 3 * k]));
      // check every constraint at x
      let rowViol = 0, coneViol = 0;
      for (let r = 0; r < R; r++) {
        let ax = 0; for (let k = lp.rs[r]; k < lp.rs[r + 1]; k++) ax += lp.rv[k] * x[lp.ri[k]];
        rowViol = Math.max(rowViol, lp.rlo[r] - ax, ax - lp.rup[r]);
      }
      const val = (fm) => { let v = fm.c; for (let k = 0; k < fm.idx.length; k++) v += fm.val[k] * x[fm.idx[k]]; return v; };
      lp.soc.forEach((F) => { coneViol = Math.max(coneViol, Math.hypot(val(F[1]), val(F[2])) - val(F[0])); });
      let obj = 0; for (let j = 0; j < n0; j++) obj += lp.cost[j] * x[j];
      Object.assign(info, { rowViol: rowViol, coneViol: coneViol });
      if (!(rowViol <= FEAS_TOL && coneViol <= FEAS_TOL)) return fail("the solution failed the admissibility check (rows " + rowViol.toExponential(1) + ", cones " + coneViol.toExponential(1) + ")", info);
      return Object.assign({ status: 7, optimal: status === 1, valid: true, crossover: false, pinf: Math.max(0, rowViol, coneViol), obj: obj, x: x, y: null, zs: zs, nVars: n0, nRows: R, solver: "Clarabel" }, info);
    };
    // Clarabel's default static regularisation first. When that gives no accepted point (a numerical error,
    // or an "almost solved" stall whose point misses the check, as on the axisymmetric lower bound and some
    // 6-node upper bounds), a larger regularisation usually converges cleanly.
    let r = attempt(1e-8);
    if (!r.valid && r.reason !== "Clarabel rejected the problem data") r = attempt(1e-7);
    return r;
  }
  // solve with the solver the problem was built for
  function solveLP(highs, lp, sense, opts) {
    if (opts._soc) return runConic(lp, sense, opts);
    const r = runLP(highs, lp, sense, opts); r.solver = "HiGHS"; return r;
  }
  // which solver produced a result, and how well its point satisfied the constraints
  function solverInfo(out, r) {
    out.solver = r.solver;
    ["nCones", "iterations", "coneStatus", "rowViol", "coneViol", "reason", "pinf"].forEach((k) => { if (r[k] != null) out[k] = r[k]; });
    return out;
  }
  // Run a finite element solver with exact cones and Clarabel when it is available; if that gives no
  // checked bound, run it again linearised with HiGHS. The result says which solver produced it.
  function conicFirst(fn) {
    return function (highs, model, mesh, mats, opts) {
      if (opts._soc != null) return fn(highs, model, mesh, mats, opts);
      if (CONIC && opts.cone !== false) {
        const t = Date.now();
        let r = null, why = "";
        try { r = fn(highs, model, mesh, mats, Object.assign({}, opts, { _soc: true })); } catch (e) { why = "Clarabel failed: " + (e && e.message ? e.message : e); }
        if (r && r.valid) return r;
        if (r === null && !why) return null;
        const r2 = fn(highs, model, mesh, mats, Object.assign({}, opts, { _soc: false }));
        if (r2) {
          r2.fallback = why || (r && r.reason ? r.reason + (r.coneStatus && !/stopped/.test(r.reason) ? " after " + r.iterations + " iterations, " + r.coneStatus : "") : "no valid cone solution");
          r2.coneTime = (Date.now() - t) / 1000;
        }
        return r2;
      }
      return fn(highs, model, mesh, mats, Object.assign({}, opts, { _soc: false }));
    };
  }

  // footing load geometry shared by the solvers
  // Loads on the whole footing are F = Ffix + alpha Fref with F = (V, H, M): V down, H along +x, M = V e
  // (positive for V to the right of the centre). Only alpha is increased to collapse. A half model carries
  // the fraction w = 1/2 of every load.
  function footingFrame(model) {
    const g = model.geo, b = g.b;
    const width = g.sym ? b : 2 * b;                     // modelled footing width
    const uFoot = porePressure(model, model.D) * width;  // water force on the footing base
    const ref = model.loadRef || [1, model.tanA || 0, model.ecc || 0];
    const fix = model.loadFix || [0, 0, 0];
    return { width: width, uFoot: uFoot, w: g.sym ? 0.5 : 1, ref: ref, fix: fix, tanA: g.sym ? 0 : (model.tanA || 0), e: g.sym ? 0 : (model.ecc || 0) };
  }
  // footing-soil interface: rough (soil strength times the factor rInter) or smooth
  function interfaceOf(model) {
    const base = cphiAt(model, soilAt(model, model.D + 1e-9), 0), r = model.rough ? (model.rInter == null ? 1 : model.rInter) : 0;
    return { tphi: r * Math.tan(base.phi), c: r * base.c, r: r };
  }

  /* ================= linearised yield surface ================= */
  // Plane-strain Mohr-Coulomb: ||(sigma_x - sigma_y, 2 tau)|| <= t = 2c cos(phi) - (sigma_x + sigma_y) sin(phi).
  // Two linearisations are available. (1) A p-sided polygon (Sloan 1988): p rows per point.
  // (2) The lifted polygon of Ben-Tal and Nemirovski (2001): nu stages of auxiliary variables give a
  //     2^(nu+1)-sided polygon P with disk(t) inside P inside disk(t / cos(pi / 2^(nu+1))), using only
  //     about 3 nu rows. The lower bound requires (x, cos(pi / 2^(nu+1)) t) in P, which implies the exact yield
  //     condition; the upper bound uses the support function of P (through LP duality), which is at least
  //     the exact dissipation. Both keep the bounds rigorous.
  const POLYS = new Map();
  function polygon(p) {
    let o = POLYS.get(p);
    if (!o) {
      const th = []; for (let k = 0; k < p; k++) th.push(2 * Math.PI * (k + 0.5) / p);
      o = { p: p, cp: Math.cos(Math.PI / p), th: th, co: th.map(Math.cos), si: th.map(Math.sin) };
      POLYS.set(p, o);
    }
    return o;
  }
  // rows of the lifted polygon in local symbols "x1", "x2", "t" and auxiliaries 0..nAux-1; each row is
  // { a: { symbol: coefficient }, lo, up } with the constant on the right
  const BTNS = new Map();
  function liftedCone(nu) {
    let o = BTNS.get(nu);
    if (o) return o;
    const xi = (j) => j, eta = (j) => nu + 1 + j, R = [];
    R.push({ a: { x1: 1, [xi(0)]: -1 }, lo: -Infinity, up: 0 }, { a: { x1: -1, [xi(0)]: -1 }, lo: -Infinity, up: 0 });
    R.push({ a: { x2: 1, [eta(0)]: -1 }, lo: -Infinity, up: 0 }, { a: { x2: -1, [eta(0)]: -1 }, lo: -Infinity, up: 0 });
    for (let j = 1; j <= nu; j++) {
      const an = Math.PI / Math.pow(2, j + 1), c = Math.cos(an), s = Math.sin(an);
      R.push({ a: { [xi(j)]: 1, [xi(j - 1)]: -c, [eta(j - 1)]: -s }, lo: 0, up: 0 });
      R.push({ a: { [xi(j - 1)]: -s, [eta(j - 1)]: c, [eta(j)]: -1 }, lo: -Infinity, up: 0 });
      R.push({ a: { [xi(j - 1)]: s, [eta(j - 1)]: -c, [eta(j)]: -1 }, lo: -Infinity, up: 0 });
    }
    const tRow = R.length;
    R.push({ a: { [xi(nu)]: 1, t: -1 }, lo: -Infinity, up: 0 });
    R.push({ a: { [eta(nu)]: 1, [xi(nu)]: -Math.tan(Math.PI / Math.pow(2, nu + 1)) }, lo: -Infinity, up: 0 });
    o = { nu: nu, rows: R, nAux: 2 * nu + 2, tRow: tRow, cp: Math.cos(Math.PI / Math.pow(2, nu + 1)) };
    BTNS.set(nu, o);
    return o;
  }
  const isAux = (key) => key !== "x1" && key !== "x2" && key !== "t";
  // A yield "cone" is the set ||(x1, x2)|| <= t with x1 = k1 (sx - sy) + c1, x2 = 2 tau, t = ct - st (sx + sy).
  //   Mohr-Coulomb:                k1 = 1, c1 = 0, ct = 2c cos(phi), st = sin(phi)
  //   tension cut-off (s1 <= kt):  k1 = 1, c1 = 0, ct = 2 kt, st = 1
  //   anisotropic undrained strength (Davis and Christian 1971): an ellipse in ((sx - sy)/2, tau) with centre
  //     X0 = su (rc - re)/2, semi-axes a = su (rc + re)/2 and su, where su is the plane-strain (simple shear)
  //     strength and rc, re are the compression and extension strengths divided by su:
  //     k1 = su / a, c1 = -2 X0 k1, ct = 2 su, st = 0
  // The soil at a point is the intersection of its cones.
  function conesAt(model, layer, y) {
    const st = cphiAt(model, layer, y), out = [];
    const an = layer.aniso;
    if (an && st.phi < 1e-12 && (Math.abs(an.rc - 1) > 1e-9 || Math.abs(an.re - 1) > 1e-9)) {
      const k1 = 2 / (an.rc + an.re);
      out.push({ k1: k1, c1: -st.c * (an.rc - an.re) * k1, ct: 2 * st.c, st: 0 });
    } else out.push({ k1: 1, c1: 0, ct: 2 * st.c * Math.cos(st.phi), st: Math.sin(st.phi) });
    if (layer.kt != null && (st.phi < 1e-12 || layer.kt < st.c / Math.tan(st.phi))) out.push({ k1: 1, c1: 0, ct: 2 * Math.max(0, layer.kt), st: 1 });
    return out;
  }
  // a lifted-cone row on the stress components with t scaled by f; returns [[c_sx, c_sy, c_tau], constant]
  // (axisymmetry adds the hoop stress: t = ct - st (sr + sz) - sth s_theta, and a fourth coefficient)
  function coneRowSigma(a, cone, f) {
    const x1 = a.x1 || 0, x2 = a.x2 || 0, t = a.t || 0;
    return [[x1 * cone.k1 - t * f * cone.st, -x1 * cone.k1 - t * f * cone.st, 2 * x2, -t * f * (cone.sth || 0)], x1 * cone.c1 + t * f * cone.ct];
  }
  // lower bound: inscribed yield conditions at one stress point (variables ix, iy, it); returns the factor on
  // the cone radius that the linearisation guarantees
  // (free: the point is on a traction-free surface, cone mode only)
  function addYieldLB(lp, ix, iy, it, cones, opts, free) {
    let cp = 1;
    cones.forEach((cn) => {
      if (opts._soc) {   // exact: (ct - st (sx + sy), k1 (sx - sy) + c1, 2 tau) in the second-order cone
        const forms = [[[ix, iy], [-cn.st, -cn.st], cn.ct], [[ix, iy], [cn.k1, -cn.k1], cn.c1], [[it], [2], 0]];
        // A zero tension cut-off at a traction-free point (stress sigma_t along the surface only) leaves the
        // cone no interior: there it is exactly sigma_t = sx + sy <= 0, which is solved instead (and the
        // cone still checked)
        if (free && cn.ct === 0 && cn.c1 === 0 && cn.st === 1 && cn.k1 === 1) { lp.addRow([ix, iy], [1, 1], -Infinity, 0); lp.addSOC(forms, true); return; }
        lp.addSOC(forms);
        return;
      }
      if (!opts.nu) {
        const pg = polygon(opts.p); cp = pg.cp;
        for (let k = 0; k < pg.p; k++) lp.addRow([ix, iy, it], [pg.co[k] * cn.k1 + pg.cp * cn.st, -pg.co[k] * cn.k1 + pg.cp * cn.st, 2 * pg.si[k]], -Infinity, pg.cp * cn.ct - pg.co[k] * cn.c1);
        return;
      }
      const cone = liftedCone(opts.nu), aux = []; cp = cone.cp;
      for (let k = 0; k < cone.nAux; k++) aux.push(lp.addVar(-Infinity, Infinity, 0));
      cone.rows.forEach((r) => {
        const [cs, k0] = coneRowSigma(r.a, cn, cone.cp), idx = [ix, iy, it], val = cs.slice(0, 3);
        Object.keys(r.a).forEach((key) => { if (isAux(key)) { idx.push(aux[+key]); val.push(r.a[key]); } });
        lp.addRow(idx, val, r.lo === -Infinity ? -Infinity : r.lo - k0, r.up - k0);
      });
    });
    return cp;
  }
  // exact conditions for a stress state with tau = 0: sx = (var ix or none) + cx, sy = (var iy or none) + cy
  function addYieldTau0(lp, ix, cx, iy, cy, cones) {
    cones.forEach((cn) => {
      [[cn.k1 + cn.st, -cn.k1 + cn.st, cn.ct - cn.c1], [-cn.k1 + cn.st, cn.k1 + cn.st, cn.ct + cn.c1]].forEach(([ax, ay, r]) => {
        const idx = [], val = [];
        if (ix != null) { idx.push(ix); val.push(ax); }
        if (iy != null) { idx.push(iy); val.push(ay); }
        lp.addRow(idx, val, -Infinity, r - ax * cx - ay * cy);
      });
    });
  }
  // upper bound: plastic strain rates eps = [[idx, val] for eps_x, eps_y, gamma_xy] must be a sum over the
  // cones of non-negative combinations of their outward normals (circumscribed polygon, or the dual of the
  // lifted polygon). Costs, already multiplied by `area`, give the dissipation. Returns { vars, lam } where
  // lam lists [variable, weight] for a plastic strain measure.
  function addFlowUB(lp, eps, cones, area, opts) {
    const nc = eps.length, vars = [], rowsS = [[], [], [], []], lam = [];
    cones.forEach((cn) => {
      if (opts._soc) {
        // exact: by conic duality, for any y = (y0, y1, y2) in the second-order cone with
        // eps = (st y0 - k1 y1, st y0 + k1 y1, -2 y2, sth y0), the dissipation is at most ct y0 + c1 y1
        const y = [0, 1, 2].map(() => lp.addVar(-Infinity, Infinity, 0));
        lp.cost[y[0]] = cn.ct * area; lp.cost[y[1]] = cn.c1 * area;
        rowsS[0].push([y[0], cn.st], [y[1], -cn.k1]); rowsS[1].push([y[0], cn.st], [y[1], cn.k1]); rowsS[2].push([y[2], -2]);
        if (nc === 4) rowsS[3].push([y[0], cn.sth || 0]);
        lp.addSOC(y.map((v) => [[v], [1], 0]));
        vars.push(...y); lam.push([y[0], 2]);
        return;
      }
      if (!opts.nu) {
        const pg = polygon(opts.p);
        for (let k = 0; k < pg.p; k++) {
          const v = lp.addVar(0, Infinity, (cn.ct - pg.co[k] * cn.c1) * area);
          rowsS[0].push([v, pg.co[k] * cn.k1 + cn.st]); rowsS[1].push([v, -pg.co[k] * cn.k1 + cn.st]); rowsS[2].push([v, 2 * pg.si[k]]);
          if (nc === 4) rowsS[3].push([v, cn.sth || 0]);
          vars.push(v); lam.push([v, 2]);
        }
        return;
      }
      // support function of {sigma : A sigma + G z <= b} = min b.y over y >= 0 (free for equality rows)
      // with A^T y = eps and G^T y = 0
      const cone = liftedCone(opts.nu), auxRows = [];
      for (let k = 0; k < cone.nAux; k++) auxRows.push([]);
      cone.rows.forEach((r, ri) => {
        const [cs, k0] = coneRowSigma(r.a, cn, 1);
        const v = lp.addVar(r.lo === r.up ? -Infinity : 0, Infinity, (r.up - k0) * area);
        for (let q = 0; q < nc; q++) if (cs[q] !== 0) rowsS[q].push([v, cs[q]]);
        Object.keys(r.a).forEach((key) => { if (isAux(key)) auxRows[+key].push([v, r.a[key]]); });
        vars.push(v);
        if (ri === cone.tRow) lam.push([v, 2]);
      });
      auxRows.forEach((ar) => lp.addRow(ar.map((q) => q[0]), ar.map((q) => q[1]), 0, 0));
    });
    for (let q = 0; q < nc; q++) lp.addRow(eps[q][0].concat(rowsS[q].map((z) => z[0])), eps[q][1].concat(rowsS[q].map((z) => -z[1])), 0, 0);
    return { vars: vars, lam: lam };
  }
  // dissipation per unit slip and unit length on a line with unit tangent t and normal n (velocity jump
  // measured from the side away from n to the side n points into), for slip along +t and -t, from the
  // strength of `layer` at level y: the maximum of +/- tau_nt over the yield surface (its c part; the
  // friction part is carried by the dilation). Also the tensile strength for opening, if any.
  function lineStrength(model, layer, y, t, n) {
    const st = cphiAt(model, layer, y), an = layer.aniso;
    let plus = st.c, minus = st.c;
    if (an && st.phi < 1e-12) {
      const lX = t[0] * n[0] - t[1] * n[1], lY = t[0] * n[1] + t[1] * n[0];
      const R = st.c * Math.hypot((an.rc + an.re) / 2 * lX, lY), m = st.c * (an.rc - an.re) / 2 * lX;
      plus = R + m; minus = R - m;
    }
    const kt = layer.kt != null && (st.phi < 1e-12 || layer.kt < st.c / Math.tan(st.phi)) ? Math.max(0, layer.kt) : null;
    return { plus: plus, minus: minus, tphi: Math.tan(st.phi), kt: kt };
  }

  /* ================= lower bound ================= */
  function lowerBound(highs, model, mesh, mats, opts) {
    const nE = mesh.nE, lp = new LP(), g = model.geo, fr = footingFrame(model);
    for (let k = 0; k < nE * 9; k++) lp.addVar(-Infinity, Infinity, 0);
    const sv = (e, i, k) => (e * 3 + i) * 3 + k;
    for (let e = 0; e < nE; e++) {
      const gg = mesh.geo[e];
      lp.addRow([sv(e, 0, 0), sv(e, 1, 0), sv(e, 2, 0), sv(e, 0, 2), sv(e, 1, 2), sv(e, 2, 2)], gg.b.concat(gg.c), 0, 0);
      lp.addRow([sv(e, 0, 2), sv(e, 1, 2), sv(e, 2, 2), sv(e, 0, 1), sv(e, 1, 1), sv(e, 2, 1)], gg.b.concat(gg.c), mats[e].gamma, mats[e].gamma);
    }
    mesh.ifaces.forEach((f) => {
      const [nx, ny] = f.n;
      const cn = [nx * nx, ny * ny, 2 * nx * ny], ct = [-nx * ny, nx * ny, nx * nx - ny * ny];
      ["tP", "tQ"].forEach((tk) => {
        [cn, ct].forEach((cf) => {
          const idx = [], val = [];
          [[f.s1, 1], [f.s2, -1]].forEach(([s, sg]) => {
            const t = s[tk];
            for (let k = 0; k < 3; k++) { idx.push(sv(s.e, s.i, k), sv(s.e, s.j, k)); val.push(sg * cf[k] * (1 - t), sg * cf[k] * t); }
          });
          lp.addRow(idx, val, 0, 0);
        });
      });
    });
    // footing: effective vertical force V' = -int(sigma_y), H = int(tau), M = -int(sigma_y x); total V = V' + uFoot.
    // The objective is the load multiplier alpha: V' + uFoot = w (Vfix + alpha Vref), and likewise for H and M.
    const Av = lp.addVar(-Infinity, Infinity, 1), itf = interfaceOf(model);
    const vI = [Av], vC = [-fr.w * fr.ref[0]], hI = [], hC = [], mI = [], mC = [];
    const tipAt = new Set(), freePt = new Set();   // traction-free surface nodes and stress points (cone mode)
    mesh.bnd.forEach((bd) => {
      const L = Math.hypot(bd.Q[0] - bd.P[0], bd.Q[1] - bd.P[1]);
      const ends = [[bd.i, bd.P], [bd.j, bd.Q]];
      if (bd.type === "surface") {
        const tx = (bd.Q[0] - bd.P[0]) / L, ty = (bd.Q[1] - bd.P[1]) / L, nx = ty, ny = -tx;
        ends.forEach(([i, X]) => {
          const q = g.qs(X[0]);
          lp.addRow([sv(bd.e, i, 0), sv(bd.e, i, 1), sv(bd.e, i, 2)], [nx * nx, ny * ny, 2 * nx * ny], -q, -q);
          lp.addRow([sv(bd.e, i, 0), sv(bd.e, i, 1), sv(bd.e, i, 2)], [-nx * ny, nx * ny, nx * nx - ny * ny], 0, 0);
          if (opts._soc && Math.abs(q) < 1e-12) { tipAt.add(key(X[0]) + ":" + key(X[1])); freePt.add(sv(bd.e, i, 0)); }
        });
      } else if (bd.type === "footing") {
        const x1 = bd.P[0], x2 = bd.Q[0];
        ends.forEach(([i], k) => {
          lp.addRow([sv(bd.e, i, 1)], [1], -Infinity, 0);            // no tension under the footing
          if (!model.rough) lp.addRow([sv(bd.e, i, 2)], [1], 0, 0);
          else if (itf.r < 1) {                                       // partial roughness: |tau| <= a - sigma_y tan(delta)
            lp.addRow([sv(bd.e, i, 2), sv(bd.e, i, 1)], [1, itf.tphi], -Infinity, itf.c);
            lp.addRow([sv(bd.e, i, 2), sv(bd.e, i, 1)], [-1, itf.tphi], -Infinity, itf.c);
          }
          vI.push(sv(bd.e, i, 1)); vC.push(-L / 2);
          hI.push(sv(bd.e, i, 2)); hC.push(L / 2);
          const w = k === 0 ? (2 * x1 + x2) : (x1 + 2 * x2);
          mI.push(sv(bd.e, i, 1)); mC.push(Math.abs(x2 - x1) / 6 * w);
        });
      } else if (bd.type === "sym") {
        [bd.i, bd.j].forEach((i) => lp.addRow([sv(bd.e, i, 2)], [1], 0, 0));
      }
    });
    const vRhs = fr.w * fr.fix[0] - fr.uFoot;
    lp.addRow(vI, vC, vRhs, vRhs);
    if (!g.sym) {
      // H = Hfix + alpha Href;  M = Mfix + alpha Mref, with M = -int(sigma_y x) (the uniform water force has no moment)
      lp.addRow(hI.concat([Av]), hC.concat([-fr.ref[1]]), fr.fix[1], fr.fix[1]);
      lp.addRow(mI.concat([Av]), mC.map((c) => -c).concat([-fr.ref[2]]), fr.fix[2], fr.fix[2]);
    }
    /* Extension of the stress field to the infinite half-space (makes the bound rigorous):
       beside the mesh the field is uniform in x with tau = 0 and geostatic sigma_y under the far ground;
       below the mesh tau = 0 and sigma_x = S (one value) plus equal geostatic increments in sigma_x and
       sigma_y, which never move a stress state closer to yield. Yield is checked at the boundary nodes
       and at the top of every deeper layer. */
    if (!opts.noExtension) {
      // the extension beside the mesh is uniform in x, so the ground must be level from xR outwards
      if (g.tb !== undefined && g.drop(model.xR) < g.HsReal - 1e-9 * Math.max(1, model.B)) throw new Error("Lower bound: the domain must extend past the toe of the slope");
      const sTop = sigmaV(model, model.D);
      // beside the mesh the vertical stress is geostatic under the ground there, with the surcharge g.qs that
      // actually acts on it (a slope lower than the embedment lowers the ground without reaching the base level)
      const geoL = (y) => g.qs(model.xL) + sigmaV(model, model.D - y) - sTop;
      const yR = g.groundY(model.xR);
      const geoR = (y) => (yR < -1e-12 ? sigmaV(model, model.D - y) - sigmaV(model, model.D - yR) : g.qs(model.xR) + sigmaV(model, model.D - y) - sTop);
      // every layer whose top is at or below the mesh base (same tolerance as fixedDepths, which does not
      // make a boundary within 1e-6 of the base an element edge)
      const deeper = model.layers.filter((l) => l.top > model.D + model.H - 1e-6);
      const H = model.H, inc = (L) => -(sigmaV(model, L.top) - sigmaV(model, model.D + H));
      const S = lp.addVar(-Infinity, Infinity, 0);
      const tied = { left: false, right: false };
      mesh.bnd.forEach((bd) => {
        if (bd.type !== "left" && bd.type !== "right" && bd.type !== "base") return;
        [[bd.i, bd.P], [bd.j, bd.Q]].forEach(([i, X]) => {
          const y = X[1], layer = mats[bd.e].layer;
          lp.addRow([sv(bd.e, i, 2)], [1], 0, 0);
          if (bd.type !== "base") {
            const geoF = bd.type === "left" ? geoL : geoR;
            addYieldTau0(lp, sv(bd.e, i, 0), 0, null, -geoF(y), conesAt(model, layer, y));
            if (Math.abs(y + H) < 1e-9 * Math.max(1, H) && !tied[bd.type]) {
              tied[bd.type] = true;
              lp.addRow([S, sv(bd.e, i, 0)], [1, -1], 0, 0);
              deeper.forEach((L) => { const Dl = inc(L); addYieldTau0(lp, S, Dl, null, -geoF(-H) + Dl, conesAt(model, L, model.D - L.top)); });
            }
          } else {
            addYieldTau0(lp, S, 0, sv(bd.e, i, 1), 0, conesAt(model, layer, y));
            deeper.forEach((L) => { const Dl = inc(L); addYieldTau0(lp, S, Dl, sv(bd.e, i, 1), Dl, conesAt(model, L, model.D - L.top)); });
          }
        });
      });
    }
    // A traction-free point of cohesionless soil can only carry zero stress (the apex of its cone), and at
    // the optimum the points that meet it at the same node (the footing corner, say) are at or near the apex
    // too. The exact cone has no interior there, which stalls an interior-point solver, so every stress
    // point at such a node is fixed to zero. Fixing stresses only restricts the field, so the bound stays
    // rigorous.
    if (tipAt.size) {
      for (let e = 0; e < nE; e++) {
        const X = mesh.XY[e];
        for (let i = 0; i < 3; i++) {
          if (!tipAt.has(key(X[2 * i]) + ":" + key(X[2 * i + 1]))) continue;
          if (!conesAt(model, mats[e].layer, X[2 * i + 1]).some((cn) => cn.ct === 0 && cn.c1 === 0 && cn.st < 1)) continue;
          for (let k = 0; k < 3; k++) { lp.lo[sv(e, i, k)] = 0; lp.up[sv(e, i, k)] = 0; }
        }
      }
    }
    // yield rows (linearised) or cones (exact) of each stress point, for the refinement indicator
    const yRow = new Int32Array(nE * 3 + 1), ySoc = new Int32Array(nE * 3 + 1);
    let cpY = 1;
    for (let e = 0; e < nE; e++) {
      const X = mesh.XY[e];
      for (let i = 0; i < 3; i++) {
        yRow[e * 3 + i] = lp.rlo.length; ySoc[e * 3 + i] = lp.soc.length;
        cpY = addYieldLB(lp, sv(e, i, 0), sv(e, i, 1), sv(e, i, 2), conesAt(model, mats[e].layer, X[2 * i + 1]), opts, freePt.has(sv(e, i, 0)));
      }
    }
    yRow[nE * 3] = lp.rlo.length; ySoc[nE * 3] = lp.soc.length;
    const r = solveLP(highs, lp, "max", opts);
    const V = fr.w * (fr.fix[0] + r.obj * fr.ref[0]);
    const out = { status: r.status, valid: r.valid, crossover: r.crossover, nVars: r.nVars, nRows: r.nRows, alpha: r.obj, V: V, qu: V / fr.width };
    solverInfo(out, r);
    if (r.status === 8 || r.status === 9) out.infeasible = true;
    if (r.x) {
      const util = new Float64Array(nE * 3), sig = new Float64Array(nE * 9), ind = new Float64Array(nE);
      for (let e = 0; e < nE; e++) for (let i = 0; i < 3; i++) {
        const sx = r.x[sv(e, i, 0)], sy = r.x[sv(e, i, 1)], t = r.x[sv(e, i, 2)];
        sig[sv(e, i, 0)] = sx; sig[sv(e, i, 1)] = sy; sig[sv(e, i, 2)] = t;
        const s = cphiAt(model, mats[e].layer, mesh.XY[e][2 * i + 1]);
        const R = 2 * s.c * Math.cos(s.phi) - (sx + sy) * Math.sin(s.phi), q = Math.hypot(sx - sy, 2 * t);
        util[e * 3 + i] = R > 1e-12 ? Math.min(1, q / (R * cpY)) : 1;
        // shear dissipation indicator: plastic multipliers (row duals) times the deviatoric stress
        if (r.y) { let lam = 0; for (let k = yRow[e * 3 + i]; k < yRow[e * 3 + i + 1]; k++) lam += Math.abs(r.y[k]); ind[e] += lam * q; }
        else if (r.zs) { let lam = 0; for (let k = ySoc[e * 3 + i]; k < ySoc[e * 3 + i + 1]; k++) lam += Math.abs(r.zs[k]); ind[e] += lam * q; }
      }
      out.util = util; out.sig = sig; if (r.y || r.zs) out.ind = ind;
    }
    return out;
  }

  // Rigid footing velocity (u0, v0, w0). The power of the multiplied (reference) load is normalised to 1, and
  // w is the modelled fraction, so the objective is w alpha = D - W_gravity - W_surcharge - uFoot v0 - w P_fix, where
  // P_fix = Vfix (-v0) + Hfix u0 - Mfix w0 is the power of the fixed loads. u0 is created by the caller.
  function footingPower(lp, fr, u0, sym) {
    const v0 = lp.addVar(-Infinity, Infinity, -fr.uFoot + fr.w * fr.fix[0]);
    const w0 = lp.addVar(sym ? 0 : -Infinity, sym ? 0 : Infinity, fr.w * fr.fix[2]);
    lp.addRow([v0, w0, u0], [-fr.ref[0], -fr.ref[2], fr.ref[1]], 1, 1);
    return [v0, w0];
  }

  /* ================= finite element upper bound ================= */
  function upperBound(highs, model, mesh, mats, opts) {
    const nE = mesh.nE, lp = new LP(), g = model.geo, fr = footingFrame(model);
    for (let k = 0; k < nE * 6; k++) lp.addVar(-Infinity, Infinity, 0);
    const uv = (e, i, k) => (e * 3 + i) * 2 + k;
    const flows = [];
    for (let e = 0; e < nE; e++) {
      const gg = mesh.geo[e];
      for (let i = 0; i < 3; i++) lp.cost[uv(e, i, 1)] += mats[e].gamma * gg.A / 3;
      const U = [uv(e, 0, 0), uv(e, 1, 0), uv(e, 2, 0)], V = [uv(e, 0, 1), uv(e, 1, 1), uv(e, 2, 1)];
      flows.push(addFlowUB(lp, [[U, gg.b.slice()], [V, gg.c.slice()], [U.concat(V), gg.c.concat(gg.b)]], conesAt(model, mats[e].layer, gg.cy), gg.A, opts));
    }
    const discVars = [];
    mesh.ifaces.forEach((f, fi) => {
      const [tx, ty] = f.t, [nx, ny] = f.n;
      const m1 = mats[f.s1.e].layer, m2 = mats[f.s2.e].layer;
      // a jump on a layer boundary may be taken in either layer: use the weaker one
      let lay = m1;
      if (m1 !== m2) {
        const sv0 = model.q0 + sigmaV(model, model.D - f.P[1]) - sigmaV(model, model.D);
        const tau = (L) => { const a1 = cphiAt(model, L, f.P[1]), a2 = cphiAt(model, L, f.Q[1]); return a1.c + a2.c + 2 * sv0 * Math.tan(a1.phi); };
        if (tau(m2) < tau(m1)) lay = m2;
      }
      const ls = [lineStrength(model, lay, f.P[1], f.t, f.n), lineStrength(model, lay, f.Q[1], f.t, f.n)];
      const tphi = ls[0].tphi, w6 = (k, j) => f.L / 6 * (j === 0 ? 2 * ls[0][k] + ls[1][k] : ls[0][k] + 2 * ls[1][k]);
      const dv = [];
      ["tP", "tQ"].forEach((tk, j) => {
        const dp = lp.addVar(0, Infinity, w6("plus", j)), dm = lp.addVar(0, Infinity, w6("minus", j));
        dv.push(dp, dm);
        // opening beyond dilatancy costs the tensile strength when there is a cut-off (otherwise none allowed)
        const so = ls[0].kt != null ? lp.addVar(0, Infinity, ls[0].kt * f.L / 2) : null;
        if (so != null) dv.push(so);
        const idx = [], vt = [], vn = [];
        [[f.s2, 1], [f.s1, -1]].forEach(([s, sg]) => {
          const t = s[tk];
          idx.push(uv(s.e, s.i, 0), uv(s.e, s.j, 0), uv(s.e, s.i, 1), uv(s.e, s.j, 1));
          vt.push(sg * tx * (1 - t), sg * tx * t, sg * ty * (1 - t), sg * ty * t);
          vn.push(sg * nx * (1 - t), sg * nx * t, sg * ny * (1 - t), sg * ny * t);
        });
        lp.addRow(idx.concat([dp, dm]), vt.concat([-1, 1]), 0, 0);
        lp.addRow(idx.concat([dp, dm], so != null ? [so] : []), vn.concat([-tphi, -tphi], so != null ? [-1] : []), 0, 0);
      });
      discVars.push([fi, dv]);
    });
    // rigid footing (u0, v0, w): velocity at base point x is (u0, v0 + w x)
    const u0 = lp.addVar(g.sym ? 0 : -Infinity, g.sym ? 0 : Infinity, -fr.w * fr.fix[1]);
    const [v0, w0] = footingPower(lp, fr, u0, g.sym);
    const fixV = (v, val) => { lp.lo[v] = val; lp.up[v] = val; };
    mesh.bnd.forEach((bd) => {
      const L = Math.hypot(bd.Q[0] - bd.P[0], bd.Q[1] - bd.P[1]);
      const ends = [[bd.i, bd.P], [bd.j, bd.Q]];
      if (bd.type === "surface") {
        const q1 = g.qs(bd.P[0]), q2 = g.qs(bd.Q[0]), dx = Math.abs(bd.Q[0] - bd.P[0]);
        lp.cost[uv(bd.e, bd.i, 1)] += dx / 6 * (2 * q1 + q2);   // surcharge acts per unit horizontal length
        lp.cost[uv(bd.e, bd.j, 1)] += dx / 6 * (q1 + 2 * q2);
      } else if (bd.type === "footing") {
        // no-tension interface: free opening s >= 0; slip costs c|slip| with dilation tan(phi) when rough
        const itf = interfaceOf(model), tphi = itf.tphi, cc = itf.c;
        ends.forEach(([i, X]) => {
          const dp = lp.addVar(0, Infinity, cc * L / 2), dm = lp.addVar(0, Infinity, cc * L / 2), so = lp.addVar(0, Infinity, 0);
          lp.addRow([u0, uv(bd.e, i, 0), dp, dm], [1, -1, -1, 1], 0, 0);
          lp.addRow([v0, w0, uv(bd.e, i, 1), dp, dm, so], [1, X[0], -1, -tphi, -tphi, -1], 0, 0);
        });
      } else if (bd.type === "sym") {
        [bd.i, bd.j].forEach((i) => fixV(uv(bd.e, i, 0), 0));
      } else {
        [bd.i, bd.j].forEach((i) => { fixV(uv(bd.e, i, 0), 0); fixV(uv(bd.e, i, 1), 0); });
      }
    });
    const r = solveLP(highs, lp, "min", opts);
    const V = fr.w * fr.fix[0] + r.obj * fr.ref[0];
    const out = solverInfo({ status: r.status, valid: r.valid, crossover: r.crossover, nVars: r.nVars, nRows: r.nRows, alpha: r.obj / fr.w, V: V, qu: V / fr.width }, r);
    if (r.x) {
      const vel = new Float64Array(nE * 6);
      for (let k = 0; k < nE * 6; k++) vel[k] = r.x[k];
      const diss = new Float64Array(nE), strain = new Float64Array(nE);
      for (let e = 0; e < nE; e++) {
        let s = 0, gg = 0;
        flows[e].vars.forEach((v) => (s += r.x[v] * lp.cost[v]));
        flows[e].lam.forEach(([v, w]) => (gg += w * Math.abs(r.x[v])));
        diss[e] = s; strain[e] = gg * mesh.geo[e].A;
      }
      discVars.forEach(([fi, dv]) => {
        const f = mesh.ifaces[fi];
        let s = 0, gg = 0; dv.forEach((v) => { s += r.x[v] * lp.cost[v]; gg += r.x[v] * f.L / 2; });
        diss[f.s1.e] += s / 2; diss[f.s2.e] += s / 2; strain[f.s1.e] += gg / 2; strain[f.s2.e] += gg / 2;
      });
      for (let e = 0; e < nE; e++) { diss[e] /= mesh.geo[e].A; strain[e] /= mesh.geo[e].A; }
      out.vel = vel; out.diss = diss; out.strain = strain;
      out.footing = { u: r.x[u0], v: r.x[v0], w: r.x[w0] };
    }
    return out;
  }

  /* ================= axisymmetric bounds (circular footing, central vertical load, level ground) ================= */
  // Coordinates: x is the radius r (axis at x = 0), y is up. Full Mohr-Coulomb in three principal stresses:
  // with s = (sr + sz)/2 and R the in-plane Mohr radius, the hoop stress s_theta is principal and the pairs give
  //   R <= c cos(phi) - s sin(phi)                                   (in-plane)
  //   R (1 - sin) <= 2c cos - s_theta (1 + sin) + s (1 - sin)        (s_theta against the minor in-plane stress)
  //   R (1 + sin) <= 2c cos - s (1 + sin) + s_theta (1 - sin)        (the major in-plane stress against s_theta)
  // each a cone ||(sr - sz, 2 tau)|| <= t, linearised as in plane strain.
  function axiConesAt(model, layer, y) {
    const st = cphiAt(model, layer, y), sn = Math.sin(st.phi), cc = st.c * Math.cos(st.phi);
    return [{ k1: 1, c1: 0, ct: 2 * cc, st: sn, sth: 0 },
      { k1: 1, c1: 0, ct: 4 * cc / (1 - sn), st: -1, sth: 2 * (1 + sn) / (1 - sn) },
      { k1: 1, c1: 0, ct: 4 * cc / (1 + sn), st: 1, sth: -2 * (1 - sn) / (1 + sn) }];
  }
  // geostatic effective vertical stress below the base level (tension positive), including the surcharge
  const sigGeo = (model, y) => -(model.q0 + sigmaV(model, model.D - y) - sigmaV(model, model.D));

  // Lower bound. The stress is split into the geostatic part (0, sz_g(y), 0, 0) and a part whose r-weighted
  // components S = r (sr, sz, tau) are quadratic in each triangle (six stress nodes per triangle, discontinuous
  // between triangles). The second equilibrium equation d(S_tau)/dr + d(S_z)/dz = 0 (a linear function) is
  // imposed at the three vertices, so it holds exactly; the first defines the hoop stress
  // s_theta = d(S_r)/dr + d(S_tau)/dz, linear in the triangle. Every yield row multiplied by r is then a
  // quadratic in the triangle; it is imposed on the six Bernstein coefficients of that quadratic, which is
  // sufficient for it to hold everywhere. Tractions are continuous across edges (quadratic, so at three
  // points). Beside the mesh the field is uniform in r (sr = s_theta, tau = 0, geostatic sz); below it
  // sr = s_theta = S0 plus geostatic increments, as in plane strain.
  const QPTS = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]];
  const QEDG = [[0, 1], [1, 2], [2, 0]];
  // Bernstein coefficients of a quadratic from its values at the six points, as linear forms [idx, val]
  function bernsteinOf(vals) {
    const out = vals.slice(0, 3);
    QEDG.forEach(([i, j], k) => {
      const m = vals[3 + k], A1 = vals[i], A2 = vals[j];
      out.push([m[0].concat(A1[0], A2[0]), m[1].map((v) => 2 * v).concat(A1[1].map((v) => -0.5 * v), A2[1].map((v) => -0.5 * v)), 2 * m[2] - 0.5 * (A1[2] + A2[2])]);
    });
    return out;
  }
  function axiLowerBound(highs, model, mesh, mats, opts) {
    const nE = mesh.nE, lp = new LP(), g = model.geo, b = g.b;
    for (let k = 0; k < nE * 18; k++) lp.addVar(-Infinity, Infinity, 0);
    const sv = (e, n, k) => (e * 6 + n) * 3 + k;   // n: 0..2 vertices, 3..5 mid-sides (0,1), (1,2), (2,0)
    const tolR = 1e-9 * Math.max(1, model.B);
    const ptXY = (e, n) => { const X = mesh.XY[e], L = QPTS[n]; return [L[0] * X[0] + L[1] * X[2] + L[2] * X[4], L[0] * X[1] + L[1] * X[3] + L[2] * X[5]]; };
    const dN = (e, L) => {
      const gg = mesh.geo[e], dx = [], dy = [];
      for (let i = 0; i < 3; i++) { dx.push((4 * L[i] - 1) * gg.b[i]); dy.push((4 * L[i] - 1) * gg.c[i]); }
      QEDG.forEach(([i, j]) => { dx.push(4 * (L[j] * gg.b[i] + L[i] * gg.b[j])); dy.push(4 * (L[j] * gg.c[i] + L[i] * gg.c[j])); });
      return { dx: dx, dy: dy };
    };
    // hoop stress at a point as a linear form
    const hoop = (e, L) => { const d = dN(e, L), idx = [], val = []; for (let n = 0; n < 6; n++) { idx.push(sv(e, n, 0), sv(e, n, 2)); val.push(d.dx[n], d.dy[n]); } return [idx, val]; };
    for (let e = 0; e < nE; e++) {
      for (let j = 0; j < 3; j++) {
        const d = dN(e, QPTS[j]), idx = [], val = [];
        for (let n = 0; n < 6; n++) { idx.push(sv(e, n, 2), sv(e, n, 1)); val.push(d.dx[n], d.dy[n]); }
        lp.addRow(idx, val, 0, 0);
      }
      for (let n = 0; n < 6; n++) if (Math.abs(ptXY(e, n)[0]) < tolR) for (let k = 0; k < 3; k++) { lp.lo[sv(e, n, k)] = 0; lp.up[sv(e, n, k)] = 0; }
    }
    const midOf = (sd) => 3 + ((sd.j === (sd.i + 1) % 3) ? sd.i : sd.j);
    mesh.ifaces.forEach((f) => {
      const [nx, ny] = f.n;
      const cn = [nx * nx, ny * ny, 2 * nx * ny], ct = [-nx * ny, nx * ny, nx * nx - ny * ny];
      [["i", "i"], ["m", "m"], ["j", "j"]].forEach(([which]) => {
        [cn, ct].forEach((cf) => {
          const idx = [], val = [];
          [[f.s1, 1], [f.s2, -1]].forEach(([sd, sg]) => {
            const n = which === "m" ? midOf(sd) : sd[which];
            for (let k = 0; k < 3; k++) { idx.push(sv(sd.e, n, k)); val.push(sg * cf[k]); }
          });
          lp.addRow(idx, val, 0, 0);
        });
      });
    });
    const q0 = model.q0, u = porePressure(model, model.D), itf = interfaceOf(model);
    const Av = lp.addVar(-Infinity, Infinity, 1);
    const vI = [Av], vC = [-1];
    // a quadratic along a boundary edge with values at (P, M, Q) given as linear forms: impose <= 0 on its
    // Bernstein coefficients (or = 0 on its values)
    const edgeLE = (fP, fM, fQ) => {
      const coef = [fP, [fM[0].concat(fP[0], fQ[0]), fM[1].map((v) => 2 * v).concat(fP[1].map((v) => -0.5 * v), fQ[1].map((v) => -0.5 * v)), 2 * fM[2] - 0.5 * (fP[2] + fQ[2])], fQ];
      coef.forEach(([idx, val, c0]) => lp.addRow(idx, val, -Infinity, -c0));
    };
    const deeper = model.layers.filter((l) => l.top > model.D + model.H - 1e-6);
    const H = model.H, inc = (L) => -(sigmaV(model, L.top) - sigmaV(model, model.D + H));
    let SR = null;
    const tauRows = (cones) => { const out = []; cones.forEach((cn) => out.push([cn.k1 + cn.st, -cn.k1 + cn.st, cn.ct - cn.c1], [-cn.k1 + cn.st, cn.k1 + cn.st, cn.ct + cn.c1])); return out; };
    mesh.bnd.forEach((bd) => {
      const pts = [[bd.e, bd.i], [bd.e, midOf(bd)], [bd.e, bd.j]];
      const xy = pts.map(([e, n]) => ptXY(e, n));
      if (bd.type === "surface") pts.forEach(([e, n]) => { lp.addRow([sv(e, n, 1)], [1], 0, 0); lp.addRow([sv(e, n, 2)], [1], 0, 0); });
      else if (bd.type === "footing") {
        // no tension: S_z - q0 r <= 0
        edgeLE(...pts.map(([e, n], k) => [[sv(e, n, 1)], [1], -q0 * xy[k][0]]));
        if (!model.rough) pts.forEach(([e, n]) => lp.addRow([sv(e, n, 2)], [1], 0, 0));
        else if (itf.r < 1) [1, -1].forEach((sg) => edgeLE(...pts.map(([e, n], k) => [[sv(e, n, 2), sv(e, n, 1)], [sg, itf.tphi], -(itf.c + q0 * itf.tphi) * xy[k][0]])));
        // V' = -2 pi int S_z dr + pi q0 b^2 (Simpson is exact for the quadratic)
        const L = Math.abs(bd.Q[0] - bd.P[0]);
        [[pts[0], 1], [pts[1], 4], [pts[2], 1]].forEach(([[e, n], w]) => { vI.push(sv(e, n, 1)); vC.push(-2 * Math.PI * L * w / 6); });
      } else if ((bd.type === "right" || bd.type === "base") && !opts.noExtension) {
        pts.forEach(([e, n]) => lp.addRow([sv(e, n, 2)], [1], 0, 0));
        const layer = mats[bd.e].layer;
        if (bd.type === "right") {
          const R = xy[0][0];
          // R (px sr + py sz - rr) = px S_r + R (py g(y) - rr(y)) <= 0, quadratic along the edge
          const rowsAt = xy.map((p2) => tauRows(conesAt(model, layer, p2[1])));
          for (let k = 0; k < rowsAt[0].length; k++) edgeLE(...pts.map(([e, n], m) => { const [px, py, rr] = rowsAt[m][k]; return [[sv(e, n, 0)], [px], R * (py * sigGeo(model, xy[m][1]) - rr)]; }));
          pts.forEach(([e, n], m) => {
            if (Math.abs(xy[m][1] + H) < 1e-9 * Math.max(1, H) && SR === null) {
              SR = lp.addVar(-Infinity, Infinity, 0);
              lp.addRow([SR, sv(e, n, 0)], [1, -1 / R], 0, 0);
              deeper.forEach((Ly) => { const Dl = inc(Ly); tauRows(conesAt(model, Ly, model.D - Ly.top)).forEach(([px, py, rr]) => lp.addRow([SR], [px], -Infinity, rr - px * Dl - py * (sigGeo(model, -H) + Dl))); });
            }
          });
        }
      }
    });
    if (!opts.noExtension) {
      if (SR === null) SR = lp.addVar(-Infinity, Infinity, 0);
      mesh.bnd.forEach((bd) => {
        if (bd.type !== "base") return;
        const pts = [[bd.e, bd.i], [bd.e, midOf(bd)], [bd.e, bd.j]], xy = pts.map(([e, n]) => ptXY(e, n)), layer = mats[bd.e].layer;
        const gz = sigGeo(model, -H);
        const sets = [[conesAt(model, layer, -H), 0]].concat(deeper.map((Ly) => [conesAt(model, Ly, model.D - Ly.top), inc(Ly)]));
        sets.forEach(([cones, D]) => tauRows(cones).forEach(([px, py, rr]) => {
          // r (px (S0 + D) + py (S_z / r + gz + D) - rr) <= 0
          edgeLE(...pts.map(([e, n], m) => [[SR, sv(e, n, 1)], [px * xy[m][0], py], xy[m][0] * (px * D + py * (gz + D) - rr)]));
        }));
      });
    }
    const Vrhs = -Math.PI * b * b * (q0 + u);
    lp.addRow(vI, vC, Vrhs, Vrhs);
    // Exact cones: points of a traction-free surface edge of cohesionless soil with no surcharge carry zero
    // stress (the apex of every cone, so no interior there): fix them, with a zero hoop stress, and leave out
    // the cones of the Bernstein coefficients that then vanish. Vertices of other triangles at those surface
    // nodes (the footing corner) are fixed too, as in plane strain; that only restricts the field.
    const tipPt = new Map();   // element -> six flags
    if (opts._soc && Math.abs(q0) < 1e-12) {
      const fix = (e, n) => {
        const fl = tipPt.get(e) || [0, 0, 0, 0, 0, 0];
        if (fl[n]) return;
        fl[n] = 1; tipPt.set(e, fl);
        for (let k = 0; k < 3; k++) { lp.lo[sv(e, n, k)] = 0; lp.up[sv(e, n, k)] = 0; }
        if (ptXY(e, n)[0] > tolR) { const [hi, hv] = hoop(e, QPTS[n]); lp.addRow(hi, hv, 0, 0); }
      };
      const nodes = new Set(), cohesionless = (e, n) => axiConesAt(model, mats[e].layer, ptXY(e, n)[1])[0].ct === 0;
      mesh.bnd.forEach((bd) => {
        if (bd.type !== "surface") return;
        const pts = [bd.i, midOf(bd), bd.j];
        if (!pts.every((n) => cohesionless(bd.e, n))) return;
        pts.forEach((n) => fix(bd.e, n));
        [bd.P, bd.Q].forEach((X) => nodes.add(key(X[0]) + ":" + key(X[1])));
      });
      if (nodes.size) for (let e = 0; e < nE; e++) for (let n = 0; n < 3; n++) {
        const X = ptXY(e, n);
        if (nodes.has(key(X[0]) + ":" + key(X[1])) && cohesionless(e, n)) fix(e, n);
      }
    }
    // yield: Bernstein coefficients of r * (row) <= 0 in each triangle
    const yRow = new Int32Array(nE + 1), ySoc = new Int32Array(nE + 1);
    let cpY = 1;
    for (let e = 0; e < nE; e++) {
      yRow[e] = lp.rlo.length; ySoc[e] = lp.soc.length;
      const layer = mats[e].layer, xy = [0, 1, 2, 3, 4, 5].map((n) => ptXY(e, n));
      const rP = xy.map((p2) => Math.max(0, p2[0])), gP = xy.map((p2) => sigGeo(model, p2[1]));
      const conesP = xy.map((p2) => axiConesAt(model, layer, p2[1]));
      const hoopP = QPTS.map((L) => hoop(e, L));
      if (opts._soc) {
        // r times (t, x1, x2) of each cone is quadratic in the triangle; each of its six Bernstein coefficient
        // triples in the cone puts every point in it (a convex combination of cone points, r >= 0)
        const tip = tipPt.get(e) || [0, 0, 0, 0, 0, 0];
        const coefTip = [tip[0], tip[1], tip[2]].concat(QEDG.map(([i, j], k) => tip[i] && tip[j] && tip[3 + k]));
        for (let ci = 0; ci < conesP[0].length; ci++) {
          const cn0 = conesP[0][ci];
          const comp = (which) => [0, 1, 2, 3, 4, 5].map((n) => {
            const cn = conesP[n][ci], idx = [], val = [];
            if (which === 0) {          // r t = r ct - st (S_r + S_z + r g) - sth r s_theta
              idx.push(sv(e, n, 0), sv(e, n, 1)); val.push(-cn0.st, -cn0.st);
              if (cn0.sth && rP[n] > 0) { const [hi, hv] = hoopP[n]; hi.forEach((v, k) => { idx.push(v); val.push(-cn0.sth * rP[n] * hv[k]); }); }
              return [idx, val, rP[n] * (cn.ct - cn0.st * gP[n])];
            }
            if (which === 1) {          // r x1 = k1 (S_r - S_z - r g) + r c1
              idx.push(sv(e, n, 0), sv(e, n, 1)); val.push(cn0.k1, -cn0.k1);
              return [idx, val, rP[n] * (cn.c1 - cn0.k1 * gP[n])];
            }
            idx.push(sv(e, n, 2)); val.push(2);   // r x2 = 2 S_tau
            return [idx, val, 0];
          });
          const B = [0, 1, 2].map((w) => bernsteinOf(comp(w)));
          for (let k = 0; k < 6; k++) if (!coefTip[k]) lp.addSOC([B[0][k], B[1][k], B[2][k]]);
        }
        continue;
      }
      for (let ci = 0; ci < conesP[0].length; ci++) {
        let tmpl;
        if (!opts.nu) {
          const pg = polygon(opts.p); cpY = pg.cp;
          tmpl = [];
          for (let k = 0; k < pg.p; k++) {
            const cn0 = conesP[0][ci];
            tmpl.push({ cs: [pg.co[k] * cn0.k1 + pg.cp * cn0.st, -pg.co[k] * cn0.k1 + pg.cp * cn0.st, 2 * pg.si[k], pg.cp * (cn0.sth || 0)], aux: {}, beta: conesP.map((c6) => pg.cp * c6[ci].ct - pg.co[k] * c6[ci].c1), eq: false });
          }
        } else {
          const cone = liftedCone(opts.nu); cpY = cone.cp;
          tmpl = cone.rows.map((r) => {
            const [cs] = coneRowSigma(r.a, conesP[0][ci], cone.cp), aux = {};
            Object.keys(r.a).forEach((key) => { if (isAux(key)) aux[key] = r.a[key]; });
            return { cs: cs, aux: aux, beta: conesP.map((c6) => r.up - coneRowSigma(r.a, c6[ci], cone.cp)[1]), eq: r.lo === r.up };
          });
        }
        const auxV = {};   // auxiliary quadratics (r times the lifted-cone variables): six values each
        tmpl.forEach((t) => {
          const vals = [0, 1, 2, 3, 4, 5].map((n) => {
            const idx = [sv(e, n, 0), sv(e, n, 1), sv(e, n, 2)], val = [t.cs[0], t.cs[1], t.cs[2]];
            if (t.cs[3] && rP[n] > 0) { const [hi, hv] = hoopP[n]; hi.forEach((v, k) => { idx.push(v); val.push(t.cs[3] * rP[n] * hv[k]); }); }
            Object.keys(t.aux).forEach((key) => {
              if (!auxV[key]) auxV[key] = [0, 1, 2, 3, 4, 5].map(() => lp.addVar(-Infinity, Infinity, 0));
              idx.push(auxV[key][n]); val.push(t.aux[key]);
            });
            return [idx, val, rP[n] * (t.cs[1] * gP[n] - t.beta[n])];
          });
          if (t.eq) vals.forEach(([idx, val, c0]) => lp.addRow(idx, val, -c0, -c0));
          else bernsteinOf(vals).forEach(([idx, val, c0]) => lp.addRow(idx, val, -Infinity, -c0));
        });
      }
    }
    yRow[nE] = lp.rlo.length; ySoc[nE] = lp.soc.length;
    const r = solveLP(highs, lp, "max", opts);
    const area = Math.PI * b * b;
    const out = solverInfo({ status: r.status, valid: r.valid, crossover: r.crossover, nVars: r.nVars, nRows: r.nRows, alpha: r.obj, V: r.obj, qu: r.obj / area, axisym: true }, r);
    if (r.status === 8 || r.status === 9) out.infeasible = true;
    if (r.x) {
      const util = new Float64Array(nE * 3), ind = new Float64Array(nE);
      for (let e = 0; e < nE; e++) {
        for (let i = 0; i < 3; i++) {
          const [rr, y] = ptXY(e, i);
          if (rr < tolR) continue;
          const sx = r.x[sv(e, i, 0)] / rr, sy = r.x[sv(e, i, 1)] / rr + sigGeo(model, y), t = r.x[sv(e, i, 2)] / rr;
          const st = cphiAt(model, mats[e].layer, y), R = 2 * st.c * Math.cos(st.phi) - (sx + sy) * Math.sin(st.phi), q = Math.hypot(sx - sy, 2 * t);
          util[e * 3 + i] = R > 1e-12 ? Math.min(1, q / (R * cpY)) : 1;
        }
        if (r.y) for (let k = yRow[e]; k < yRow[e + 1]; k++) ind[e] += Math.abs(r.y[k]);
        else if (r.zs) for (let k = ySoc[e]; k < ySoc[e + 1]; k++) ind[e] += Math.abs(r.zs[k]);
      }
      out.util = util; if (r.y || r.zs) out.ind = ind;
    }
    return out;
  }

  // Upper bound: quadratic velocities (u radial, v vertical) in each triangle, with quadratic velocity jumps
  // between triangles (which avoids the locking that the hoop strain causes in undrained soil). The strain rates multiplied by r
  // (r e_r, r e_z, r g_rz and r e_theta = u) are quadratic in each triangle. The plastic multipliers, also
  // multiplied by r, are quadratics written in Bernstein form; non-negative Bernstein coefficients make them
  // non-negative everywhere, and matching the flow rule coefficient by coefficient makes it hold everywhere.
  // Along each edge the slip is split into non-negative quadratics (Bernstein form) with the associated
  // dilation. The dissipation 2 pi int D r dA and the work terms are integrated exactly. Needs a conforming mesh.
  function axiUpperBound(highs, model, mesh, mats, opts) {
    if (!mesh.conforming) return null;
    const nE = mesh.nE, lp = new LP(), g = model.geo, b = g.b, TP = 2 * Math.PI;
    const tolR = 1e-9 * Math.max(1, model.B);
    const pos = [], E = [];
    for (let e = 0; e < nE; e++) {   // six nodes per triangle, not shared (jumps are allowed between triangles)
      const X = mesh.XY[e], cn = [0, 1, 2].map((i) => { pos.push([X[2 * i], X[2 * i + 1]]); return pos.length - 1; });
      const md = [0, 1, 2].map((i) => { const j = (i + 1) % 3; pos.push([(X[2 * i] + X[2 * j]) / 2, (X[2 * i + 1] + X[2 * j + 1]) / 2]); return pos.length - 1; });
      E.push({ cn: cn, md: md });
    }
    const N = pos.length;
    for (let k = 0; k < 2 * N; k++) lp.addVar(-Infinity, Infinity, 0);
    const U = (n) => 2 * n, Vv = (n) => 2 * n + 1;
    // Bernstein points of a quadratic on the triangle: vertices 0..2, then edges (0,1), (1,2), (2,0)
    const EDG = [[0, 1], [1, 2], [2, 0]];
    const strainOut = new Float64Array(nE), lamOf = [];
    for (let e = 0; e < nE; e++) {
      const gg = mesh.geo[e], el = E[e], bb = gg.b, cc = gg.c, X = mesh.XY[e];
      const rV = [X[0], X[2], X[4]].map((v) => Math.max(0, v));
      const lay = mats[e].layer, st = [0, 1, 2].map((j) => cphiAt(model, lay, X[2 * j + 1]));
      // self-weight: int N_n r dA for the six Lagrange shape functions
      const wV = (i, m) => (i === m ? gg.A / 30 : -gg.A / 60), wM = (i, j, m) => (m === i || m === j ? 2 * gg.A / 15 : gg.A / 15);
      for (let i = 0; i < 3; i++) { let w = 0; for (let m = 0; m < 3; m++) w += wV(i, m) * rV[m]; lp.cost[Vv(el.cn[i])] += mats[e].gamma * TP * w; }
      EDG.forEach(([i, j], k) => { let w = 0; for (let m = 0; m < 3; m++) w += wM(i, j, m) * rV[m]; lp.cost[Vv(el.md[k])] += mats[e].gamma * TP * w; });
      // derivatives of the six shape functions at a point with area coordinates L
      const dN = (L) => {
        const dx = [], dy = [], nodes = [];
        for (let i = 0; i < 3; i++) { nodes.push(el.cn[i]); dx.push((4 * L[i] - 1) * bb[i]); dy.push((4 * L[i] - 1) * cc[i]); }
        EDG.forEach(([i, j], k) => { nodes.push(el.md[k]); dx.push(4 * (L[j] * bb[i] + L[i] * bb[j])); dy.push(4 * (L[j] * cc[i] + L[i] * cc[j])); });
        return { nodes: nodes, dx: dx, dy: dy };
      };
      const Lpt = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]];
      const rAt = (L) => L[0] * rV[0] + L[1] * rV[1] + L[2] * rV[2];
      // r-weighted strain components at the six points (vertices, edge midpoints): [idx, val] each
      const comp = Lpt.map((L, k) => {
        const d = dN(L), r = rAt(L), uIdx = d.nodes.map(U), vIdx = d.nodes.map(Vv);
        const uHere = k < 3 ? [[U(el.cn[k])], [1]] : [[U(el.md[k - 3])], [1]];
        return [[uIdx, d.dx.map((v) => v * r)], [vIdx, d.dy.map((v) => v * r)], [uIdx.concat(vIdx), d.dy.concat(d.dx).map((v) => v * r)], uHere];
      });
      // Bernstein coefficients from point values: vertex = value; edge = 2 f(mid) - (f_i + f_j)/2
      const bern = (q) => {
        const out = [];
        for (let k = 0; k < 3; k++) out.push(comp[k][q]);
        EDG.forEach(([i, j], k) => {
          const m = comp[3 + k][q], A1 = comp[i][q], A2 = comp[j][q];
          out.push([m[0].concat(A1[0], A2[0]), m[1].map((v) => 2 * v).concat(A1[1].map((v) => -0.5 * v), A2[1].map((v) => -0.5 * v))]);
        });
        return out;
      };
      const B4 = [0, 1, 2, 3].map(bern);  // B4[q][k] = [idx, val] of Bernstein coefficient k of component q
      // effective cohesion for each Bernstein coefficient: (6/A) int B_k c dA
      const cE = [0, 1, 2].map((i) => 0.6 * st[i].c + 0.2 * (st[(i + 1) % 3].c + st[(i + 2) % 3].c))
        .concat(EDG.map(([i, j]) => 0.4 * (st[i].c + st[j].c) + 0.2 * st[3 - i - j].c));
      let lamW = [];
      for (let k = 0; k < 6; k++) {
        const lk = Object.assign({}, lay, { c: cE[k], suGrad: 0, top: 0 });
        const f = addFlowUB(lp, B4.map((q) => q[k]), axiConesAt(Object.assign({}, model, { D: 0 }), lk, 0), TP * gg.A / 6, opts);
        lamW = lamW.concat(f.lam);
      }
      lamOf.push(lamW);
    }
    const simpson = (L, cP, rP, cM, rM, cQ, rQ) => [L / 6 * (cP * rP + 2 * cM * rM), L / 6 * (2 * cM * rM + cQ * rQ)];
    const gp5 = [[0.1127016653792583, 0.2777777777777778], [0.5, 0.4444444444444444], [0.8872983346207417, 0.2777777777777778]];
    // velocity jumps on the inter-element edges: Bernstein coefficients (P, middle, Q) of the jump components
    mesh.ifaces.forEach((f) => {
      const [tx, ty] = f.t, [nx, ny] = f.n;
      const sides = [f.s1, f.s2].map((sd) => { const el = E[sd.e]; return { P: el.cn[sd.i], Q: el.cn[sd.j], M: el.md[(sd.j === (sd.i + 1) % 3) ? sd.i : sd.j] }; });
      const lay = mats[f.s1.e].layer === mats[f.s2.e].layer ? mats[f.s1.e].layer : [mats[f.s1.e].layer, mats[f.s2.e].layer].sort((A, B) => cphiAt(model, A, f.P[1]).c - cphiAt(model, B, f.P[1]).c)[0];
      const tphi = Math.tan(cphiAt(model, lay, f.P[1]).phi);
      // jump at a point = side 2 minus side 1; Bernstein coefficient of the middle: 2 f(M) - (f(P) + f(Q))/2
      const jumpCoef = (k, dir) => {
        const idx = [], val = [];
        [[sides[1], 1], [sides[0], -1]].forEach(([sd, sg]) => {
          const pts = k === 0 ? [[sd.P, 1]] : k === 2 ? [[sd.Q, 1]] : [[sd.M, 2], [sd.P, -0.5], [sd.Q, -0.5]];
          pts.forEach(([n, w]) => { idx.push(U(n), Vv(n)); val.push(sg * w * dir[0], sg * w * dir[1]); });
        });
        return [idx, val];
      };
      // dissipation weights: 2 pi int B_k(s) c(s) r(s) L ds
      const wB = [0, 0, 0];
      gp5.forEach(([sp, w]) => {
        const x = f.P[0] + sp * (f.Q[0] - f.P[0]), y = f.P[1] + sp * (f.Q[1] - f.P[1]), cr = cphiAt(model, lay, y).c * Math.max(0, x) * w * f.L;
        wB[0] += cr * (1 - sp) * (1 - sp); wB[1] += cr * 2 * sp * (1 - sp); wB[2] += cr * sp * sp;
      });
      for (let k = 0; k < 3; k++) {
        const dp = lp.addVar(0, Infinity, TP * wB[k]), dm = lp.addVar(0, Infinity, TP * wB[k]);
        const [it, vt] = jumpCoef(k, [tx, ty]), [inn, vn] = jumpCoef(k, [nx, ny]);
        lp.addRow(it.concat([dp, dm]), vt.concat([-1, 1]), 0, 0);
        lp.addRow(inn.concat([dp, dm]), vn.concat([-tphi, -tphi]), 0, 0);
      }
    });
    const u = porePressure(model, model.D), area = Math.PI * b * b;
    const v0 = lp.addVar(-Infinity, Infinity, -u * area);
    lp.addRow([v0], [-1], 1, 1);                       // the footing load does unit power
    const itf = interfaceOf(model), fixV = (v) => { lp.lo[v] = 0; lp.up[v] = 0; };
    const gp = [[0.1127016653792583, 0.2777777777777778], [0.5, 0.4444444444444444], [0.8872983346207417, 0.2777777777777778]];
    mesh.bnd.forEach((bd) => {
      const el = E[bd.e], a1 = el.cn[bd.i], q1 = el.cn[bd.j];
      const m1 = el.md[(bd.j === (bd.i + 1) % 3) ? bd.i : bd.j];
      if (bd.type === "surface") {
        const dxs = Math.abs(bd.Q[0] - bd.P[0]);
        gp.forEach(([sp, w]) => {
          const x = bd.P[0] + sp * (bd.Q[0] - bd.P[0]), qv = g.qs(x) * w * dxs * TP * Math.max(0, x);
          lp.cost[Vv(a1)] += qv * (1 - sp) * (1 - 2 * sp); lp.cost[Vv(m1)] += qv * 4 * sp * (1 - sp); lp.cost[Vv(q1)] += qv * sp * (2 * sp - 1);
        });
      } else if (bd.type === "footing") {
        lp.addRow([U(m1), U(a1), U(q1)], [2, -1, -1], 0, 0);
        lp.addRow([Vv(m1), Vv(a1), Vv(q1)], [2, -1, -1], 0, 0);
        const L = Math.abs(bd.Q[0] - bd.P[0]), w = simpson(L, itf.c, bd.P[0], itf.c, 0.5 * (bd.P[0] + bd.Q[0]), itf.c, bd.Q[0]);
        [a1, q1].forEach((n, k) => {
          const dp = lp.addVar(0, Infinity, TP * w[k]), dm = lp.addVar(0, Infinity, TP * w[k]), so = lp.addVar(0, Infinity, 0);
          lp.addRow([U(n), dp, dm], [-1, -1, 1], 0, 0);
          lp.addRow([v0, Vv(n), dp, dm, so], [1, -1, -itf.tphi, -itf.tphi, -1], 0, 0);
        });
      } else if (bd.type === "sym") [a1, m1, q1].forEach((n) => fixV(U(n)));
      else [a1, m1, q1].forEach((n) => { fixV(U(n)); fixV(Vv(n)); });
    });
    for (let n = 0; n < N; n++) if (Math.abs(pos[n][0]) < tolR) fixV(U(n));
    const r = solveLP(highs, lp, "min", opts);
    const out = solverInfo({ status: r.status, valid: r.valid, crossover: r.crossover, nVars: r.nVars, nRows: r.nRows, alpha: r.obj, V: r.obj, qu: r.obj / area, axisym: true, quad: true }, r);
    if (r.x) {
      for (let e = 0; e < nE; e++) { let gm = 0; lamOf[e].forEach(([v, w]) => (gm += w * Math.abs(r.x[v]))); strainOut[e] = gm / 6; }
      out.strain = strainOut; out.footing = { u: 0, v: r.x[v0], w: 0 };
    }
    return out;
  }

  /* ================= quadratic (6-node) finite element upper bound ================= */
  // Continuous quadratic velocity triangles (Yu, Sloan and Kleeman 1994; Optum's "Upper" element). Strain
  // rates are linear, so the flow rule imposed at the three vertices with plastic multipliers interpolated
  // linearly holds everywhere in the element, and the dissipation and self-weight work are integrated
  // exactly. There are no velocity jumps between elements; the footing interface is kept linear along each
  // footing edge so that its conditions at the edge ends hold along the whole edge. Needs a conforming mesh.
  function upperBoundQ(highs, model, mesh, mats, opts) {
    if (!mesh.conforming) return null;
    const nE = mesh.nE, lp = new LP(), g = model.geo, fr = footingFrame(model);
    const nodeId = new Map(), mid = new Map(), pos = [];
    const corner = (x, y) => { const k = key(x) + ":" + key(y); let id = nodeId.get(k); if (id === undefined) { id = pos.length; nodeId.set(k, id); pos.push([x, y]); } return id; };
    const E = [];
    for (let e = 0; e < nE; e++) {
      const X = mesh.XY[e], cn = [corner(X[0], X[1]), corner(X[2], X[3]), corner(X[4], X[5])];
      const md = [0, 1, 2].map((i) => {
        const a = cn[i], b = cn[(i + 1) % 3], kk = a < b ? a + "," + b : b + "," + a;
        let id = mid.get(kk);
        if (id === undefined) { id = pos.length; mid.set(kk, id); pos.push([(pos[a][0] + pos[b][0]) / 2, (pos[a][1] + pos[b][1]) / 2]); }
        return id;
      });
      E.push({ cn: cn, md: md }); // md[i] lies between corners i and i+1
    }
    const N = pos.length;
    for (let k = 0; k < 2 * N; k++) lp.addVar(-Infinity, Infinity, 0);
    const U = (n) => 2 * n, Vv = (n) => 2 * n + 1;
    const flows = [];
    for (let e = 0; e < nE; e++) {
      const gg = mesh.geo[e], el = E[e], b = gg.b, c = gg.c;
      el.md.forEach((m) => (lp.cost[Vv(m)] += mats[e].gamma * gg.A / 3));    // self-weight acts on the mid-side nodes
      const lay = mats[e].layer, st = [0, 1, 2].map((j) => cphiAt(model, lay, mesh.XY[e][2 * j + 1]));
      const csum = st[0].c + st[1].c + st[2].c;
      // cone at vertex j with the cohesion replaced by (c_j + sum c)/4 (exact integral of linear times linear)
      const conesJ = (j) => conesAt(Object.assign({}, model, { D: 0 }), Object.assign({}, lay, { c: (st[j].c + csum) / 4, suGrad: 0, top: 0 }), 0);
      const fl = [];
      for (let j = 0; j < 3; j++) {
        // derivatives of the six shape functions at vertex j
        const nodes = [], dx = [], dy = [];
        for (let i = 0; i < 3; i++) { nodes.push(el.cn[i]); dx.push((i === j ? 3 : -1) * b[i]); dy.push((i === j ? 3 : -1) * c[i]); }
        for (let i = 0; i < 3; i++) {
          const k = (i + 1) % 3; // mid-side node between corners i and k: N = 4 L_i L_k
          nodes.push(el.md[i]); dx.push(4 * ((k === j ? b[i] : 0) + (i === j ? b[k] : 0))); dy.push(4 * ((k === j ? c[i] : 0) + (i === j ? c[k] : 0)));
        }
        const eps = [[nodes.map(U), dx], [nodes.map(Vv), dy], [nodes.map(U).concat(nodes.map(Vv)), dy.concat(dx)]];
        // exact integral of (linear multiplier) x (linear cohesion): weight (c_j + sum c) / 4 on area A / 3
        fl.push(addFlowUB(lp, eps, conesJ(j), gg.A / 3, opts));
      }
      flows.push(fl);
    }
    const u0 = lp.addVar(g.sym ? 0 : -Infinity, g.sym ? 0 : Infinity, -fr.w * fr.fix[1]);
    const [v0, w0] = footingPower(lp, fr, u0, g.sym);
    const fixV = (v) => { lp.lo[v] = 0; lp.up[v] = 0; };
    const gp = [[0.1127016653792583, 0.2777777777777778], [0.5, 0.4444444444444444], [0.8872983346207417, 0.2777777777777778]];
    const itf = interfaceOf(model);
    mesh.bnd.forEach((bd) => {
      const el = E[bd.e], a = el.cn[bd.i], q = el.cn[bd.j];
      const m = el.md[(bd.j === (bd.i + 1) % 3) ? bd.i : bd.j];
      if (bd.type === "surface") {
        const dxs = Math.abs(bd.Q[0] - bd.P[0]);
        gp.forEach(([s, w]) => {
          const x = bd.P[0] + s * (bd.Q[0] - bd.P[0]), qv = g.qs(x) * w * dxs;
          lp.cost[Vv(a)] += qv * (1 - s) * (1 - 2 * s); lp.cost[Vv(m)] += qv * 4 * s * (1 - s); lp.cost[Vv(q)] += qv * s * (2 * s - 1);
        });
      } else if (bd.type === "footing") {
        // soil velocity linear along the footing edge, then the no-tension interface at both ends
        lp.addRow([U(m), U(a), U(q)], [2, -1, -1], 0, 0);
        lp.addRow([Vv(m), Vv(a), Vv(q)], [2, -1, -1], 0, 0);
        const L = Math.hypot(bd.Q[0] - bd.P[0], bd.Q[1] - bd.P[1]);
        [[a, bd.P], [q, bd.Q]].forEach(([n, X]) => {
          const dp = lp.addVar(0, Infinity, itf.c * L / 2), dm = lp.addVar(0, Infinity, itf.c * L / 2), so = lp.addVar(0, Infinity, 0);
          lp.addRow([u0, U(n), dp, dm], [1, -1, -1, 1], 0, 0);
          lp.addRow([v0, w0, Vv(n), dp, dm, so], [1, X[0], -1, -itf.tphi, -itf.tphi, -1], 0, 0);
        });
      } else if (bd.type === "sym") {
        [a, m, q].forEach((n) => fixV(U(n)));
      } else {
        [a, m, q].forEach((n) => { fixV(U(n)); fixV(Vv(n)); });
      }
    });
    const r = solveLP(highs, lp, "min", opts);
    const V = fr.w * fr.fix[0] + r.obj * fr.ref[0];
    const out = solverInfo({ status: r.status, valid: r.valid, crossover: r.crossover, nVars: r.nVars, nRows: r.nRows, alpha: r.obj / fr.w, V: V, qu: V / fr.width, quad: true }, r);
    if (r.x) {
      const diss = new Float64Array(nE), strain = new Float64Array(nE);
      for (let e = 0; e < nE; e++) {
        let s = 0, gm = 0;
        flows[e].forEach((f) => { f.vars.forEach((v) => (s += r.x[v] * lp.cost[v])); f.lam.forEach(([v, w]) => (gm += w * Math.abs(r.x[v]) / 3)); });
        diss[e] = s / mesh.geo[e].A; strain[e] = gm;
      }
      out.diss = diss; out.strain = strain;
      out.footing = { u: r.x[u0], v: r.x[v0], w: r.x[w0] };
    }
    return out;
  }

  /* ================= discontinuity layout optimisation ================= */
  // weight of material (plus surcharge) vertically above the line between its end points
  function columnWeight(model, x1, y1, x2, y2) {
    if (Math.abs(x2 - x1) < 1e-12) return 0;
    if (x2 < x1) { [x1, x2] = [x2, x1]; [y1, y2] = [y2, y1]; }
    const g = model.geo, b = g.b;
    const yl = (x) => y1 + (y2 - y1) * (x - x1) / (x2 - x1);
    const cuts = [x1, x2];
    [b, -b].concat(g.slope ? [g.xc, g.xt] : []).concat(g.surfX).concat(g.tb !== undefined ? [g.xcReal] : []).forEach((x) => { if (x > x1 && x < x2) cuts.push(x); });
    const depths = model.layers.map((l) => l.bottom).concat(model.waterDepth != null ? [model.waterDepth] : []);
    depths.forEach((d) => {
      const y = model.D - d;
      if ((y1 - y) * (y2 - y) < 0) cuts.push(x1 + (x2 - x1) * (y - y1) / (y2 - y1));
      if (g.tb !== undefined) { // where the ground surface crosses this depth on the slope
        const xf = g.xcReal + d / g.tb;
        if (d > 0 && d < g.HsReal && xf > x1 && xf < x2) cuts.push(xf);
      }
    });
    cuts.sort((a, c) => a - c);
    const f = (x) => g.qs(x) + sigmaV(model, model.D - yl(x)) - sigmaV(model, model.D - g.groundY(x));
    let W = 0;
    for (let k = 0; k < cuts.length - 1; k++) {
      const a = cuts[k], c = cuts[k + 1];
      if (c - a < 1e-14) continue;
      W += (c - a) * f(0.5 * (a + c)); // f is linear within each piece, so the midpoint rule is exact
    }
    return W;
  }
  function dloAxis(model, h) {
    const g = model.geo, b = g.b, ext = model.ext;
    const nf = Math.max(2, Math.round(b / h)), hf = b / nf;
    const side = (lim, reach) => {
      const xs = [];
      let x = b, step = hf;
      while (x + step < lim - 0.4 * step) { x += step; xs.push(x); if (x > Math.min(lim, reach)) step *= 1.3; else step = Math.min(step * 1.18, 4 * hf); }
      xs.push(lim);
      return xs;
    };
    let xs = [];
    for (let i = 0; i <= nf; i++) xs.push(b * i / nf);
    xs = xs.concat(side(model.xR, ext.xR));
    if (!g.sym) { for (let i = 1; i <= nf; i++) xs.push(-b * i / nf); xs = xs.concat(side(-model.xL, -ext.xL).map((x) => -x)); }
    if (g.slope) [g.xc, g.xt].forEach((k) => { if (k < model.xR) xs = xs.filter((x) => Math.abs(x - k) > 0.3 * hf).concat([k]); });
    return uniqSorted(xs, 1e-9);
  }
  function dloDepth(model, h) {
    let ds = [0], d = 0, step = h;
    const H = model.H, dext = model.ext.depth;
    while (d + step < H - 0.4 * step) { d += step; ds.push(d); if (d > dext) step *= 1.3; else step = Math.min(step * 1.18, 4 * h); }
    ds.push(H);
    const extra = fixedDepths(model);
    if (extra.length) { ds = ds.filter((v) => v === 0 || v === H || extra.every((e) => Math.abs(e - v) > 0.3 * h)); ds = uniqSorted(ds.concat(extra), 1e-9); }
    return ds;
  }
  function dloNodes(model, h) {
    const g = model.geo, xs = dloAxis(model, h), ds = dloDepth(model, h), tol = 1e-9 * Math.max(1, model.B);
    const nodes = [], seen = new Set();
    const add = (x, y) => { const k = key(x) + ":" + key(y); if (!seen.has(k) && g.inside(x, y)) { seen.add(k); nodes.push([x, y]); } };
    ds.forEach((d) => xs.forEach((x) => {
      // drop grid nodes that would crowd the slope face; face nodes are added below
      if (g.slope && -d > -g.Hs + tol && -d < -tol && x < g.slopeX(-d) && g.slopeX(-d) - x < 0.3 * h) return;
      add(x, -d);
    }));
    if (g.slope) {
      ds.forEach((d) => { if (d <= g.Hs + tol) { const x = g.slopeX(-d); if (x <= model.xR + tol) add(x, -d); } });
      xs.forEach((x) => { if (x > g.xc && x < g.xt) add(x, g.groundY(x)); });
    }
    return nodes;
  }
  function upperBoundDLO(highs, model, opts) {
    const g = model.geo, b = g.b, tol = 1e-9 * Math.max(1, model.B), fr = footingFrame(model);
    let lo = b / 40, hi = b, nodes = dloNodes(model, hi);
    for (let it = 0; it < 30; it++) { const mid = Math.sqrt(lo * hi), n = dloNodes(model, mid); if (n.length > (opts.dloNodes || 400)) lo = mid; else { hi = mid; nodes = n; } }
    const N = nodes.length;
    const onTop = (p) => Math.abs(p[1]) < tol && (!g.slope || p[0] <= g.xc + tol);
    const onFace = (p) => g.slope && p[1] <= tol && p[1] >= -g.Hs - tol && Math.abs(p[0] - g.slopeX(p[1])) < 1e-7 * Math.max(1, model.B);
    const onToe = (p) => g.slope && Math.abs(p[1] + g.Hs) < tol && p[0] >= g.xt - tol;
    const onL = (p) => Math.abs(p[0] - model.xL) < tol, onR = (p) => Math.abs(p[0] - model.xR) < tol, onBase = (p) => Math.abs(p[1] + model.H) < tol;
    const onGround = (p) => onTop(p) || onFace(p) || onToe(p) || (onR(p) && Math.abs(p[1] - g.groundY(p[0])) < tol);
    const underFooting = (p) => onTop(p) && Math.abs(p[0]) <= b + tol;
    // a line stays inside the (non-convex) domain if it stays left of the slope face
    const inDomain = (A, C) => {
      if (!g.slope) return true;
      const yTop = Math.min(0, Math.max(A[1], C[1])), yBot = Math.max(-g.Hs, Math.min(A[1], C[1]));
      if (yTop < yBot) return true;
      const sl = 1e-7 * model.B;
      if (Math.abs(C[1] - A[1]) < 1e-14) return Math.max(A[0], C[0]) <= g.slopeX(A[1]) + sl;
      const xAt = (y) => A[0] + (C[0] - A[0]) * (y - A[1]) / (C[1] - A[1]);
      return [yTop, yBot].every((y) => xAt(y) <= g.slopeX(y) + sl);
    };
    const lp = new LP(), lines = [];
    const layerRange = model.layers.map((l) => [l.top - model.D, l.bottom - model.D]);
    for (let a = 0; a < N; a++) for (let c = a + 1; c < N; c++) {
      const A = nodes[a], C = nodes[c];
      const same = (f) => f(A) && f(C);
      if (same(onTop) || same(onFace) || same(onToe) || same(onL) || same(onR) || same(onBase)) continue;
      if (!inDomain(A, C)) continue;
      const dA = -A[1], dC = -C[1], dm = 0.5 * (dA + dC);
      let li = 0; while (li < layerRange.length - 1 && dm > layerRange[li][1] + 1e-12) li++;
      const lr = layerRange[li];
      if (Math.min(dA, dC) < lr[0] - 1e-9 || Math.max(dA, dC) > lr[1] + 1e-9) continue; // crosses a layer boundary
      let layer = model.layers[li];
      if (Math.abs(dA - dC) < 1e-9 && li < layerRange.length - 1 && Math.abs(dA - lr[1]) < 1e-9) {
        const L2 = model.layers[li + 1], sv0 = model.q0 + sigmaV(model, model.D + dA) - sigmaV(model, model.D);
        const tauOf = (ly) => { const st = cphiAt(model, ly, -dA); return st.c + sv0 * Math.tan(st.phi); };
        if (tauOf(L2) < tauOf(layer)) layer = L2;
      }
      const L = Math.hypot(C[0] - A[0], C[1] - A[1]);
      const tx = (C[0] - A[0]) / L, ty = (C[1] - A[1]) / L, nx = -ty, ny = tx;
      // slip along +t and -t with dilation tan(phi), and an opening that costs the tensile strength when the
      // soil has a cut-off; the jump is the velocity on the n side minus the other side
      const ls = lineStrength(model, layer, 0.5 * (A[1] + C[1]), [tx, ty], [nx, ny]), tphi = ls.tphi;
      const vs = [lp.addVar(0, Infinity, ls.plus * L), lp.addVar(0, Infinity, ls.minus * L)];
      const jl = [[tx + tphi * nx, ty + tphi * ny], [-tx + tphi * nx, -ty + tphi * ny]];
      if (ls.kt != null) { vs.push(lp.addVar(0, Infinity, ls.kt * L)); jl.push([nx, ny]); }
      const up = ny > 1e-12 ? 1 : (ny < -1e-12 ? -1 : 0);
      const Wt = up !== 0 ? columnWeight(model, A[0], A[1], C[0], C[1]) : 0;
      vs.forEach((v, k) => (lp.cost[v] += Wt * up * jl[k][1]));
      lines.push({ a, c, vs, jl, up, x1: Math.min(A[0], C[0]), x2: Math.max(A[0], C[0]) });
    }
    const rows = nodes.map(() => [[], [], [], []]);
    lines.forEach((ln) => {
      [[ln.a, 1], [ln.c, -1]].forEach(([nd, sg]) => {
        const r = rows[nd];
        ln.vs.forEach((v, k) => { r[0].push(v); r[1].push(sg * ln.jl[k][0]); r[2].push(v); r[3].push(sg * ln.jl[k][1]); });
      });
    });
    // compatibility: closed loops at interior, side and base nodes; none on the ground or under the
    // footing (the footing interface carries those); x only on the symmetry axis
    nodes.forEach((p, k) => {
      if (onGround(p) || underFooting(p)) return;
      if (g.sym && onL(p) && !onBase(p)) { lp.addRow(rows[k][0], rows[k][1], 0, 0); return; }
      lp.addRow(rows[k][0], rows[k][1], 0, 0); lp.addRow(rows[k][2], rows[k][3], 0, 0);
    });
    // rigid footing and its no-tension interface with the soil under each footing segment
    const u0 = lp.addVar(g.sym ? 0 : -Infinity, g.sym ? 0 : Infinity, -fr.w * fr.fix[1]);
    const [v0, w0] = footingPower(lp, fr, u0, g.sym);
    const itf = interfaceOf(model), tphiI = itf.tphi, cI = itf.c;
    // Soil velocity under each footing segment: one column sum anchors the first segment; the others follow
    // from closing the loop at each footing node (keeps the constraint matrix sparse).
    const fx = uniqSorted(nodes.filter(underFooting).map((p) => p[0]), 1e-9);
    const nodeAt = new Map(); nodes.forEach((p, k) => { if (underFooting(p)) nodeAt.set(key(p[0]), k); });
    const us = [], vs = [];
    for (let i = 0; i < fx.length - 1; i++) { us.push(lp.addVar(-Infinity, Infinity, 0)); vs.push(lp.addVar(-Infinity, Infinity, 0)); }
    {
      const xm = 0.5 * (fx[0] + fx[1]);
      const iu = [us[0]], cu = [-1], iv = [vs[0]], cv = [-1];
      lines.forEach((ln) => {
        if (ln.up === 0 || xm <= ln.x1 || xm >= ln.x2) return;
        ln.vs.forEach((v, k) => { iu.push(v); cu.push(ln.up * ln.jl[k][0]); iv.push(v); cv.push(ln.up * ln.jl[k][1]); });
      });
      lp.addRow(iu, cu, 0, 0); lp.addRow(iv, cv, 0, 0);
    }
    if (g.sym) {
      // corner of the axis and the footing: soil next to the axis has u = 0, so the loop from the axis
      // region to the first footing segment gives u under that segment
      const k0 = nodeAt.get(key(fx[0]));
      lp.addRow(rows[k0][0].concat([us[0]]), rows[k0][1].concat([-1]), 0, 0);
    }
    for (let i = 1; i < fx.length - 1; i++) {
      // loop sum around the footing node = (soil right) - (soil left)
      const k = nodeAt.get(key(fx[i])), rw = rows[k];
      lp.addRow(rw[0].concat([us[i], us[i - 1]]), rw[1].concat([-1, 1]), 0, 0);
      lp.addRow(rw[2].concat([vs[i], vs[i - 1]]), rw[3].concat([-1, 1]), 0, 0);
    }
    for (let i = 0; i < fx.length - 1; i++) {
      const Lseg = fx[i + 1] - fx[i];
      [fx[i], fx[i + 1]].forEach((xe) => {
        const dp = lp.addVar(0, Infinity, cI * Lseg / 2), dm = lp.addVar(0, Infinity, cI * Lseg / 2), so = lp.addVar(0, Infinity, 0);
        lp.addRow([u0, us[i], dp, dm], [1, -1, -1, 1], 0, 0);
        lp.addRow([v0, w0, vs[i], dp, dm, so], [1, xe, -1, -tphiI, -tphiI, -1], 0, 0);
      });
    }
    const r = runLP(highs, lp, "min", opts);
    const V = fr.w * fr.fix[0] + r.obj * fr.ref[0];
    const out = { status: r.status, valid: r.valid, crossover: r.crossover, nVars: r.nVars, nRows: r.nRows, alpha: r.obj / fr.w, V: V, qu: V / fr.width, nodes: nodes, nLines: lines.length };
    if (r.x) {
      const act = [];
      let jmax = 0;
      lines.forEach((ln) => {
        let jx = 0, jy = 0;
        ln.vs.forEach((v, k) => { jx += ln.jl[k][0] * r.x[v]; jy += ln.jl[k][1] * r.x[v]; });
        const m = Math.hypot(jx, jy);
        if (m > 1e-6) { act.push({ A: nodes[ln.a], C: nodes[ln.c], jx, jy, m, up: ln.up, x1: ln.x1, x2: ln.x2 }); jmax = Math.max(jmax, m); }
      });
      out.lines = act; out.jmax = jmax;
      out.footing = { u: r.x[u0], v: r.x[v0], w: r.x[w0] };
      out.footSoil = fx.slice(0, -1).map((x, i) => ({ xm: 0.5 * (x + fx[i + 1]), u: r.x[us[i]], v: r.x[vs[i]] }));
    }
    return out;
  }
  function dloVelocity(ub, x, y) {
    let u = 0, v = 0;
    for (const l of ub.lines) {
      if (l.up === 0 || x <= l.x1 || x >= l.x2) continue;
      const t = (x - l.A[0]) / (l.C[0] - l.A[0]);
      const yl = l.A[1] + t * (l.C[1] - l.A[1]);
      if (yl < y) { u += l.up * l.jx; v += l.up * l.jy; }
    }
    return [u, v];
  }
  function mechanismExtent(ub, frac) {
    let xmax = 0, xmin = 0, dm = 0;
    ub.lines.forEach((l) => { if (l.m > frac * ub.jmax) { xmax = Math.max(xmax, l.A[0], l.C[0]); xmin = Math.min(xmin, l.A[0], l.C[0]); dm = Math.max(dm, -l.A[1], -l.C[1]); } });
    return { xR: xmax, xL: xmin, depth: dm };
  }
  // copy of the model on a new domain (geometry rebuilt)
  function onDomain(model, xL, xR, H, ext) {
    const m = Object.assign({}, model, { xL: xL, xR: xR, H: H });
    m.geo = geometry(m);
    m.ext = ext || { xL: Math.max(xL, model.ext.xL), xR: Math.min(xR, model.ext.xR), depth: Math.min(H, model.ext.depth) };
    return m;
  }
  function touches(ub, m) {
    const e = mechanismExtent(ub, 0.01), span = m.xR - m.xL;
    return { left: !m.geo.sym && e.xL < m.xL + 0.08 * span, right: e.xR > m.xR - 0.08 * span, base: e.depth > 0.9 * m.H, e: e };
  }
  // DLO in stages: coarse solves on a domain grown until the mechanism is clear of its edges, then a
  // finer solve on a box fitted to the mechanism (grown again if the finer mechanism reaches it)
  function dloAdaptive(highs, model, opts) {
    const n = opts.dloNodes || 420;
    const coarse = Object.assign({}, opts, { dloNodes: Math.round(n * 0.45) });
    let m1 = model, r1 = upperBoundDLO(highs, m1, coarse), grows = 0;
    for (; r1.lines && grows < 4; grows++) {
      const t = touches(r1, m1);
      if (!t.left && !t.right && !t.base) break;
      const step = 0.7 * Math.max((m1.xR - m1.xL) / 2, model.B);
      m1 = onDomain(model, t.left ? m1.xL - step : m1.xL, t.right ? m1.xR + step : m1.xR, t.base ? m1.H * 1.7 : m1.H, model.ext);
      r1 = upperBoundDLO(highs, m1, coarse);
    }
    if (!r1.lines) return r1;
    r1.domain = { xL: m1.xL, xR: m1.xR, H: m1.H };
    const e = mechanismExtent(r1, 0.01), B = model.B, g = model.geo;
    const fitR = Math.min(m1.xR, Math.max(1.35 * e.xR, 1.2 * B));
    const fitL = g.sym ? 0 : Math.max(m1.xL, Math.min(1.35 * e.xL, -1.2 * B));
    const fitH = Math.min(m1.H, Math.max(1.4 * e.depth, 0.6 * B));
    let m2 = onDomain(model, fitL, fitR, fitH, { xL: fitL, xR: fitR, depth: fitH });
    let r2 = upperBoundDLO(highs, m2, opts);
    for (let k = 0; k < 3 && r2.lines; k++) {
      const t = touches(r2, m2);
      if (!t.left && !t.right && !t.base) break;
      const nL = t.left ? Math.max(m1.xL, m2.xL * 1.5) : m2.xL, nR = t.right ? Math.min(m1.xR, m2.xR * 1.5) : m2.xR, nH = t.base ? Math.min(m1.H, m2.H * 1.5) : m2.H;
      if (nL === m2.xL && nR === m2.xR && nH === m2.H) break;
      m2 = onDomain(model, nL, nR, nH, { xL: nL, xR: nR, depth: nH });
      r2 = upperBoundDLO(highs, m2, opts);
    }
    if (r2.lines) r2.domain = { xL: m2.xL, xR: m2.xR, H: m2.H };
    const best = r2.lines && r2.valid && (!r1.valid || r2.V < r1.V) ? r2 : r1;
    best.coarse = r1.V; best.extent = mechanismExtent(best, 0.01); best.grows = grows;
    best.fullDomain = { xL: m1.xL, xR: m1.xR, H: m1.H };
    best.touches = touches(r1, m1);
    return best;
  }

  /* ================= adaptivity ================= */
  const findCell = (cells, x, y) => {
    for (let i = 0; i < cells.length; i++) { const c = cells[i]; if (x >= c.x0 && x <= c.x1 && y >= c.y0 && y <= c.y1) return i; }
    return -1;
  };
  const inFan = (c, fans) => fans.some((f) => { const dx = Math.max(c.x0 - f.cx, 0, f.cx - c.x1), dy = Math.max(c.y0 - f.cy, 0, f.cy - c.y1); return Math.hypot(dx, dy) + Math.max(c.x1 - c.x0, c.y1 - c.y0) < f.R; });
  // refinement indicator from the DLO mechanism: jump magnitude times length of the active lines in each cell
  function dloIndicator(cells, dlo, fans) {
    const ind = new Float64Array(cells.length);
    let hmin = Infinity; cells.forEach((c) => (hmin = Math.min(hmin, c.x1 - c.x0, c.y1 - c.y0)));
    dlo.lines.forEach((l) => {
      if (l.m < 0.02 * dlo.jmax) return;
      const L = Math.hypot(l.C[0] - l.A[0], l.C[1] - l.A[1]), n = Math.max(2, Math.ceil(L / hmin));
      for (let k = 0; k < n; k++) {
        const t = (k + 0.5) / n, ci = findCell(cells, l.A[0] + t * (l.C[0] - l.A[0]), l.A[1] + t * (l.C[1] - l.A[1]));
        if (ci >= 0) ind[ci] += l.m * L / n;
      }
    });
    cells.forEach((c, i) => { ind[i] *= Math.sqrt((c.x1 - c.x0) * (c.y1 - c.y0)); if (inFan(c, fans)) ind[i] = 0; });
    return ind;
  }
  // indicator from a finite element solution: per-element values summed into their cells
  function meshIndicator(cells, mesh, perElem, fans) {
    const ind = new Float64Array(cells.length);
    for (let e = 0; e < mesh.nE; e++) { const ci = mesh.cellOf[e]; if (ci >= 0 && perElem[e] > 0) ind[ci] += perElem[e]; }
    cells.forEach((c, i) => { if (inFan(c, fans)) ind[i] = 0; });
    return ind;
  }
  const normalise = (a) => { let t = 0; a.forEach((v) => (t += v)); return t > 0 ? a.map((v) => v / t) : a; };
  function markByIndicator(ind, frac) {
    const order = Array.from(ind.keys()).sort((p, q) => ind[q] - ind[p]);
    const tot = ind.reduce((p, q) => p + q, 0), marked = new Set();
    let cum = 0;
    for (const i of order) { if (cum >= 0.97 * tot || marked.size >= frac * ind.length || ind[i] <= 0) break; marked.add(i); cum += ind[i]; }
    return marked;
  }
  function markCellsDLO(cells, dlo, frac, fans) { return markByIndicator(dloIndicator(cells, dlo, fans), frac); }
  /* ================= classical bearing capacity (comparison) ================= */
  // Vesic N factors, De Beer/Vesic shape factors, Meyerhof effective width and inclination factors.
  // Returns average ultimate pressures on the full width B (effective-width pressure x B'/B).
  function classical(phiDeg, c, gamma, q, B, L, shape, opt) {
    opt = opt || {};
    const f = phiDeg * Math.PI / 180;
    let Nq, Nc, Ng;
    if (phiDeg < 1e-6) { Nq = 1; Nc = Math.PI + 2; Ng = 0; }
    else { Nq = Math.exp(Math.PI * Math.tan(f)) * Math.pow(Math.tan(Math.PI / 4 + f / 2), 2); Nc = (Nq - 1) / Math.tan(f); Ng = 2 * (Nq + 1) * Math.tan(f); }
    const e = Math.abs(opt.e || 0), Bp = Math.max(1e-9, B - 2 * e);
    const aDeg = Math.abs(opt.alphaDeg || 0);
    const ic = Math.pow(Math.max(0, 1 - aDeg / 90), 2), iq = ic, ig = phiDeg > 1e-6 ? Math.pow(Math.max(0, 1 - aDeg / phiDeg), 2) : 1;
    let sc = 1, sq = 1, sg = 1;
    if (shape !== "strip") {
      const BL = shape === "circle" ? 1 : Math.min(1, Bp / L);
      sc = 1 + BL * Nq / Nc; sq = 1 + BL * Math.tan(f); sg = Math.max(0.6, 1 - 0.4 * BL);
    }
    const term = (s1, s2, s3) => c * Nc * s1 * ic + q * Nq * s2 * iq + 0.5 * gamma * Bp * Ng * s3 * ig;
    return { Nc, Nq, Ng, sc, sq, sg, ic, iq, ig, Bp, qult: term(sc, sq, sg) * Bp / B, qultStrip: term(1, 1, 1) * Bp / B };
  }

  /* ================= driver ================= */
  function prandtlExtent(b, phiDeg) {
    const f = phiDeg * Math.PI / 180, psi = Math.PI / 4 + f / 2;
    const r0 = b / Math.cos(psi), a0 = Math.PI - psi, a1 = Math.PI / 4 - f / 2;
    let depth = b * Math.tan(psi), r1 = r0;
    for (let k = 0; k <= 60; k++) {
      const a = a0 + (a1 - a0) * k / 60, r = r0 * Math.exp((a0 - a) * Math.tan(f));
      depth = Math.max(depth, r * Math.sin(a)); r1 = r;
    }
    return { x: b + 2 * r1 * Math.cos(a1), depth: depth, r0: r0 };
  }
  // model0: {B, D, rough, gammaW, waterDepth, surcharge, layers, load: {tanA, e}, slope: {dist, height, angle} | null}
  function prepare(model0) {
    const m = Object.assign({}, model0);
    m.input = model0;
    let top = 0;
    m.layers = model0.layers.map((l) => { const o = Object.assign({}, l, { top: top }); top = l.bottom; return o; });
    m.layers[m.layers.length - 1].bottom = 1e12;
    const snapTol = 0.01 * m.B;
    m.layers.forEach((l, i) => { if (i < m.layers.length - 1 && Math.abs(l.bottom - m.D) < snapTol) l.bottom = m.D; });
    if (m.waterDepth != null) {
      if (Math.abs(m.waterDepth - m.D) < snapTol) m.waterDepth = m.D;
      for (const l of m.layers) if (Math.abs(m.waterDepth - l.bottom) < snapTol) m.waterDepth = l.bottom;
    }
    const ld = model0.load || {};
    m.tanA = ld.tanA || 0; m.ecc = ld.e || 0;
    // load vectors (V, H, M) on the whole footing, per unit length. Mode "V": V is increased to collapse with
    // H = V tan(alpha) + Hfix and M = V e + Mfix. Mode "H": V = Vfix (at eccentricity e) is held and H is
    // increased (sliding capacity).
    if (ld.mode === "H") { m.loadRef = [0, 1, 0]; m.loadFix = [ld.Vfix || 0, 0, (ld.Vfix || 0) * m.ecc + (ld.Mfix || 0)]; }
    else { m.loadRef = [1, m.tanA, m.ecc]; m.loadFix = [0, ld.Hfix || 0, ld.Mfix || 0]; }
    m.loadMode = ld.mode === "H" ? "H" : "V";
    m.slope = model0.slope && model0.slope.height > 0 ? Object.assign({}, model0.slope) : null;
    m.sym = !m.slope && m.loadMode === "V" && Math.abs(m.tanA) < 1e-12 && Math.abs(m.ecc) < 1e-12 && !(ld.Hfix) && !(ld.Mfix);
    // axisymmetric analysis of a circular footing of diameter B (central vertical load on level ground only)
    m.axisym = !!model0.axisym;
    // non-associated flow (Davis 1968): with dilation angle psi < phi, use c* = beta c, tan(phi*) = beta tan(phi),
    // beta = cos(psi) cos(phi) / (1 - sin(psi) sin(phi)). The bounds are rigorous for the reduced material; the
    // link to the real non-associated soil is approximate.
    m.layers.forEach((l) => {
      if (l.psi == null || !(l.phi > 0) || l.psi >= l.phi) return;
      const f = l.phi * Math.PI / 180, ps = Math.max(0, l.psi) * Math.PI / 180;
      const beta = Math.cos(ps) * Math.cos(f) / (1 - Math.sin(ps) * Math.sin(f));
      l.davis = { beta: beta, c: l.c, phi: l.phi };
      l.c = beta * l.c; l.suGrad = beta * (l.suGrad || 0); l.phi = Math.atan(beta * Math.tan(f)) * 180 / Math.PI;
    });
    m.q0 = sigmaV(m, m.D) + (m.surcharge || 0);
    let phi = 0;
    m.layers.forEach((l) => { if (l.top - m.D < 2 * m.B) phi = Math.max(phi, l.phi); });
    const ext = prandtlExtent(m.B / 2, phi), sc = m.domainScale || 1;
    m.xR = m.xR || Math.max(1.15 * ext.x, 1.6 * m.B) * sc;
    m.xL = m.sym ? 0 : (m.xL != null ? m.xL : -m.xR);
    m.H = m.H || Math.max(1.25 * ext.depth, 1.0 * m.B) * sc;
    m.geo = geometry(m);
    if (m.geo.slope) {
      m.xR = Math.max(m.xR, m.geo.xt + 1.0 * m.B);
      m.H = Math.max(m.H, m.geo.Hs + 1.0 * m.B);
      m.geo = geometry(m);
    } else if (m.geo.tb !== undefined && m.xR < m.geo.xtReal + 0.5 * m.B) {
      // a slope no deeper than the embedment: the ground must be level again inside the mesh (see lowerBound)
      m.xR = m.geo.xtReal + 0.5 * m.B;
      m.geo = geometry(m);
    }
    m.ext = { xL: m.sym ? 0 : Math.max(m.xL, -ext.x * 1.1), xR: Math.min(m.xR, Math.max(ext.x * 1.1, m.geo.slope ? m.geo.xt : 0)), depth: Math.min(m.H, Math.max(ext.depth * 1.1, m.geo.slope ? m.geo.Hs : 0)) };
    m.phiDomain = phi;
    return m;
  }
  function checkModel(m) {
    // the extensions check yield only at the top of each layer, which needs strength that does not fall with depth
    if (m.layers.some((l) => (l.suGrad || 0) < 0)) throw new Error("The undrained strength cannot decrease with depth; use a gradient of zero or more.");
    if (m.axisym && !m.sym) throw new Error("The axisymmetric (circular) analysis needs a central vertical load on level ground.");
    if (m.axisym && m.layers.some((l) => l.kt != null || (l.aniso && (l.aniso.rc !== 1 || l.aniso.re !== 1)))) throw new Error("Tension cut-offs and anisotropic strength are not available in the axisymmetric analysis.");
    if (Math.abs(m.ecc) >= 0.95 * m.B / 2) throw new Error("The eccentricity must be less than 0.475 B.");
    if (!m.rough && (Math.abs(m.loadRef[1]) > 1e-12 || Math.abs(m.loadFix[1]) > 1e-12)) throw new Error("A smooth base cannot carry horizontal load; tick Rough base or set the horizontal load to zero.");
    if (m.rInter != null && !(m.rInter > 0 && m.rInter <= 1)) throw new Error("The interface strength factor must be greater than 0 and at most 1.");
    if (m.slope && m.waterDepth != null && m.waterDepth < m.slope.height - 1e-9) throw new Error("A water table above the slope toe is not supported; set it at or below the toe depth.");
  }
  const tick = () => new Promise((r) => setTimeout(r, 0));

  /* opts: {p, nf, targetElems, adaptSteps, dloNodes, lower, feUpper, onProgress(msg, frac)} */
  async function analyse(highs, model0, opts) {
    const t0 = Date.now();
    const say = opts.onProgress || (() => {});
    const steps = opts.adaptSteps == null ? 3 : opts.adaptSteps;
    let model = prepare(model0);
    checkModel(model);
    const prepared = model;
    const total = 3 + (opts.lower !== false ? 3 : 0) + 2;
    let done = 0;
    const res = { warnings: [] };
    say("Upper bound: discontinuity layout optimisation", 0); await tick();
    let t1 = Date.now();
    res.dlo = dloAdaptive(highs, model, opts); res.dlo.time = (Date.now() - t1) / 1000; done += 3;
    if (!res.dlo.lines) {
      if ((res.dlo.status === 10 || res.dlo.status === 9) && !model.axisym) throw new Error(model.loadFix.some((v) => v) ? "The fixed loads alone (with the ground's own weight) exceed the capacity: no multiplied load can be added. Reduce the fixed loads." : "The ground fails under its own weight before any footing load is applied. Check the slope angle against the soil strength.");
      throw new Error("Upper bound solve failed (" + statusText(res.dlo.status) + ")");
    }
    // finite element domain sized from the mechanism, with room for the stress field to spread.
    // The lower bound stays rigorous at any domain size because the field is extended to infinity.
    const ex = res.dlo.extent, B = model.B, g0 = model.geo;
    let xR = Math.max(2.2 * ex.xR, 2.5 * B);
    if (g0.slope) xR = Math.max(xR, g0.xt + 0.75 * B);
    else if (g0.tb !== undefined) xR = Math.max(xR, g0.xtReal + 0.5 * B);   // level ground beyond the mesh (see lowerBound)
    const xL = g0.sym ? 0 : Math.min(2.2 * ex.xL, -2.5 * B);
    let H = Math.max(2.8 * ex.depth, 1.5 * B);
    // below a slope the base must be deep enough for the topographic stress difference to spread by shear
    // before the extension (which has uniform horizontal stress) takes over
    const weighty = model.layers.some((l) => l.gamma > 0 || l.gammaSat > model.gammaW);
    if (g0.slope) H = Math.max(H, g0.Hs + Math.max(1.0 * B, (opts.slopeDepth == null ? (weighty ? 4 : 0) : opts.slopeDepth) * g0.Hs));
    model = onDomain(model, xL, xR, H, { xL: Math.min(-0.5 * B, 1.1 * ex.xL), xR: Math.min(xR, Math.max(1.1 * ex.xR, g0.slope ? g0.xt : 0)), depth: Math.min(H, 1.1 * ex.depth) });
    if (g0.sym) model.ext.xL = 0;
    const b = B / 2;
    let fans = fansOf(model, opts);
    const fan = { nTheta: opts.nTheta || 16, r0: b * (opts.fanR0 || 0.04), R: fans[0].R, list: fans };
    let cells = initialCells(model, opts.nf || 4);
    const target = (opts.targetElems || 1800) * (g0.sym ? 1 : 1.3);
    const mk = (cs) => {
      try { return buildTriMesh(model, cs, fan); }
      catch (e) {
        if (fan.list.length > fans.filter((f) => f.kind === "footing").length) {  // drop the slope fans and retry
          fans = fans.filter((f) => f.kind === "footing"); fan.list = fans; res.slopeFansDropped = true;
          return mk(cs);
        }
        // the axisymmetric solvers need a conforming mesh, so there is no quadtree fallback for them
        if (model.geo.slope || model.axisym) throw e; res.meshFallback = true; return buildQuadMesh(model, cs);
      }
    };
    // Adaptive cycles (after Optum: start coarse, refine on the shear dissipation, 3 or 4 cycles). The first
    // mesh is refined from the DLO mechanism; later meshes from the plastic multipliers of the lower bound
    // (row duals) and the plastic strain of the finite element upper bound. Every mesh gives valid bounds,
    // so the best of each kind is kept.
    const cycles = Math.max(1, opts.adaptCycles || 1);
    let n0 = 0;
    const targetAt = (k) => (k === 0 ? (cycles === 1 ? target : target * (opts.startFrac || 0.5))
      : Math.max(1.15 * n0 * Math.pow(1.15, k - 1), n0 * Math.pow(target / n0, k / (cycles - 1))));
    const refineTo = (tgt, indicatorOf, steps) => {
      let m = mk(cells);
      for (let s = 0; s < steps && m.nE < tgt * 0.9; s++) {
        const growth = Math.pow(Math.max(1.05, tgt / m.nE), 1 / (steps - s));
        const ind = indicatorOf(cells, m);
        cells.forEach((c, i) => (c.ind = ind[i]));
        cells = refineCells(cells, markByIndicator(ind, Math.min(0.6, Math.max(0.05, (growth - 1) / 3))));
        m = mk(cells);
      }
      return m;
    };
    const dloInd = (cs) => dloIndicator(cs, res.dlo, fans);
    let mesh = refineTo(targetAt(0), dloInd, steps);
    n0 = mesh.nE;
    const wantLB = opts.lower !== false;
    res.layered = fixedDepths(model).some((z) => model.layers.some((l) => Math.abs(l.bottom - model.D - z) < 1e-9));
    const axi = !!model.axisym;
    const wantFE = axi || opts.feUpper === true || (opts.feUpper !== false && (res.layered || (!g0.sym && opts.feAutoCombined !== false) || (opts.feQuad && opts.feFrictional && model.layers.some((l) => l.phi > 0 && l.top < model.D + 3 * model.B))));
    res.cycles = [];
    for (let cyc = 0; cyc < cycles; cyc++) {
      const mats = assignMaterials(model, mesh);
      const tag = cycles > 1 ? ` – cycle ${cyc + 1} of ${cycles}` : "";
      const lbName = axi ? "axiLowerBound" : "lowerBound";
      const ubName = axi ? "axiUpperBound" : opts.feQuad && mesh.conforming ? "upperBoundQ" : "upperBound";
      let lb = null, ub = null;
      // The lower and upper bound on one mesh are independent, so when the host can run solves in parallel
      // (the browser app gives each its own worker and WebAssembly instance) they run together.
      if (wantLB && wantFE && opts.parallel) {
        say(`Lower bound and upper bound in parallel: finite elements (${mesh.nE} elements)${tag}`, (done + 3 * cyc / cycles) / total); await tick();
        try {
          const plainModel = Object.assign({}, model); delete plainModel.geo;
          const plainOpts = JSON.parse(JSON.stringify(Object.assign({}, opts, { onProgress: null, parallel: null, strength: null })));
          [lb, ub] = await opts.parallel([lbName, ubName].map((fn) => ({ fn: fn, model: plainModel, mesh: mesh, opts: plainOpts })));
          res.parallelUsed = true;
        } catch (e) { lb = ub = null; res.parallelFailed = e && e.message ? e.message : String(e); }
      }
      if (wantLB && !lb) {
        say(`Lower bound: finite element limit analysis (${mesh.nE} elements)${tag}`, (done + 3 * cyc / cycles) / total); await tick();
        t1 = Date.now();
        lb = api[lbName](highs, model, mesh, mats, opts); lb.time = (Date.now() - t1) / 1000;
      }
      if (lb) {
        lb.mesh = mesh;
        if (!res.lb || (lb.valid && (!res.lb.valid || lb.V > res.lb.V)) || (!res.lb.valid && !lb.valid)) res.lb = lb;
      }
      if (wantFE && !ub) {
        say(`Upper bound: finite elements (${mesh.nE} elements)${tag}`, (done + 3 * (cyc + 0.5) / cycles) / total); await tick();
        t1 = Date.now();
        ub = api[ubName](highs, model, mesh, mats, opts);
        if (ub) ub.time = (Date.now() - t1) / 1000;
      }
      if (ub) {
        ub.mesh = mesh;
        if (!res.ub || (ub.valid && (!res.ub.valid || ub.V < res.ub.V)) || (!res.ub.valid && !ub.valid)) res.ub = ub;
        res.feReason = opts.feUpper === true ? "requested" : (!g0.sym ? "combined load or slope" : "layered profile");
      }
      res.cycles.push({ nE: mesh.nE, lb: lb && lb.valid ? lb.V : NaN, ub: ub && ub.valid ? ub.V : NaN });
      if (cyc === cycles - 1) break;
      // indicators on the current cells, then refine toward the next target
      const parts = [normalise(dloInd(cells))];
      if (lb && lb.ind) parts.push(normalise(meshIndicator(cells, mesh, lb.ind, fans)));
      if (ub && ub.strain) parts.push(normalise(meshIndicator(cells, mesh, Array.from(ub.strain, (v, e) => v * mesh.geo[e].A), fans)));
      const comb = new Float64Array(cells.length);
      parts.forEach((a) => a.forEach((v, i) => (comb[i] += v / parts.length)));
      let first = true;
      mesh = refineTo(targetAt(cyc + 1), (cs) => { if (first) { first = false; return comb; } return Float64Array.from(cs, (c) => c.ind || 0); }, 2);
    }
    done += 3;
    if (res.lb) { res.mesh = res.lb.mesh; } else res.mesh = mesh;
    res.model = model;
    const ubs = [];
    // (the plane-strain DLO mechanism is not an upper bound for a circular footing; it only sized the domain)
    if (res.dlo.valid && !axi) ubs.push(["Discontinuity layout optimisation", res.dlo.alpha, res.dlo]);
    if (res.ub && res.ub.valid) ubs.push(["Finite element", res.ub.alpha, res.ub]);
    ubs.sort((p, q) => p[1] - q[1]);
    const fr = footingFrame(model);
    res.uBase = porePressure(model, model.D);
    res.sym = g0.sym;
    // collapse loads are totals per unit length for the whole footing: F = Ffix + alpha Fref
    const aLb = res.lb && res.lb.valid ? res.lb.alpha : NaN;   // an infeasible lower-bound problem proves nothing
    const aUb = ubs.length ? ubs[0][1] : NaN;
    const at = (a, k) => fr.fix[k] + a * fr.ref[k];
    // pressures: per unit length / B (plane strain), or force / (pi B^2 / 4) (axisymmetric)
    const per = axi ? Math.PI * B * B / 4 : B;
    res.axisym = axi;
    res.best = { mode: model.loadMode, alphaLb: aLb, alphaUb: aUb, lb: at(aLb, 0) / per, ub: at(aUb, 0) / per, ubMethod: ubs.length ? ubs[0][0] : "", axisym: axi,
      Vlb: at(aLb, 0), Vub: at(aUb, 0), Hlb: at(aLb, 1), Hub: at(aUb, 1), Mlb: at(aLb, 2), Mub: at(aUb, 2),
      ref: fr.ref, fix: fr.fix, tanA: fr.tanA, e: fr.e, mech: ubs.length ? ubs[0][2].footing : null };
    if (res.lb && !res.lb.valid) res.warnings.push(res.lb.infeasible ? "No statically admissible stress field was found on this mesh, so no lower bound is reported." : "The lower-bound solve did not reach a feasible solution (" + statusText(res.lb.status) + "), so no lower bound is reported.");
    if (!ubs.length) res.warnings.push("No upper-bound solve reached a feasible solution.");
    if (res.dlo.touches && (res.dlo.touches.left || res.dlo.touches.right || res.dlo.touches.base)) res.warnings.push("The mechanism still reaches the edge of the largest domain tried; the upper bound may be loose.");
    if (aUb <= 0) res.warnings.push("The upper bound is zero or negative: the ground fails under its own weight" + (fr.fix.some((v) => v) ? " and the fixed loads." : "."));
    if (opts.strength && res.lb) {
      res.strength = await strengthFS(highs, prepared, model, res, opts, say);
      if (res.strength.note) res.warnings.push(res.strength.note);
    }
    res.time = (Date.now() - t0) / 1000;
    say("Done", 1);
    return res;
  }

  /* ================= factor of safety on strength ================= */
  // Strength reduction (as in Optum's Strength Reduction analysis): c, su and tan(phi) are divided by F and
  // the applied loads are held. With alpha(F) the collapse multiplier of the reduced soil and t the applied
  // multiplier, the factor of safety is the F at which alpha(F) = t. Any F with alpha_LB(F) >= t proves
  // FS >= F, and any F with alpha_UB(F) <= t proves FS <= F, so each bound is the safe end of its own bracket.
  function reducedModel(model, F) {
    const layers = model.layers.map((l) => Object.assign({}, l, { c: l.c / F, suGrad: (l.suGrad || 0) / F, phi: Math.atan(Math.tan(l.phi * Math.PI / 180) / F) * 180 / Math.PI }));
    return Object.assign({}, model, { layers: layers });
  }
  async function strengthFS(highs, prepared, feModel, res, opts, say) {
    const t = opts.strength.target, tol = opts.strength.tol || 0.01, maxEval = opts.strength.maxEval || 9;
    const mesh = res.mesh, dloOpts = Object.assign({}, opts, { dloNodes: Math.round((opts.dloNodes || 420) * 0.8) }), axi = !!feModel.axisym;
    const evals = [];
    const lbAt = async (F) => {
      say(`Factor of safety on strength: lower bound at F = ${F.toFixed(3)}`, 0.9); await tick();
      const m = reducedModel(feModel, F), r = (axi ? axiLowerBound : lowerBound)(highs, m, mesh, assignMaterials(m, mesh), opts);
      const a = r.valid ? r.alpha : -Infinity;
      evals.push({ kind: "LB", F: F, alpha: a });
      return a;
    };
    // (status 10, unbounded, proves collapse; status 9, unbounded or infeasible, proves nothing)
    const ubAt = async (F) => {
      say(`Factor of safety on strength: upper bound at F = ${F.toFixed(3)}`, 0.95); await tick();
      let a;
      if (axi) {
        // the plane-strain DLO is not an upper bound for the circle: use the axisymmetric solver on its mesh
        const um = (res.ub && res.ub.mesh) || mesh, m = reducedModel(feModel, F), r = axiUpperBound(highs, m, um, assignMaterials(m, um), opts);
        a = r && r.valid ? r.alpha : (r && r.status === 10 ? -Infinity : NaN);
      } else {
        const r = dloAdaptive(highs, reducedModel(prepared, F), dloOpts);
        a = r.lines && r.valid ? r.alpha : (r.status === 10 ? -Infinity : NaN);
      }
      evals.push({ kind: "UB", F: F, alpha: a });
      return a;
    };
    // Secant steps on (log F, log alpha). Each step aims just past the estimated root on whichever side is
    // still missing from the bracket (first the proving side), with geometric bisection as the fallback.
    // proven(alpha) says whether an evaluation proves the bound at that F; proofs hold for F below the FS
    // (lower bound, provenBelow) or above it (upper bound).
    const bracket = async (fn, F0, a0, proven, provenBelow) => {
      let pf = null, nf = null, n = 0;
      const pts = [[F0, a0]];
      if (proven(a0)) pf = F0; else nf = F0;
      const near = (F) => pts.some((q) => Math.abs(Math.log(q[0] / F)) < 0.2 * tol);
      while (n < maxEval && !(pf !== null && nf !== null && Math.abs(Math.log(nf / pf)) <= Math.log(1 + tol))) {
        const good = pts.filter((q) => q[1] > 0 && isFinite(q[1])).sort((u, v) => Math.abs(Math.log(u[1] / t)) - Math.abs(Math.log(v[1] / t))).slice(0, 2);
        let est = NaN;
        if (good.length === 2 && Math.abs(Math.log(good[1][0] / good[0][0])) > 1e-9) {
          const k = -Math.log(good[1][1] / good[0][1]) / Math.log(good[1][0] / good[0][0]);
          if (k > 0.05) est = good[0][0] * Math.pow(good[0][1] / t, 1 / k);
        } else if (good.length >= 1 && t > 0) est = good[0][0] * good[0][1] / t;   // alpha roughly proportional to 1/F
        // which side to aim at: the proving side first, then the other side, then the side whose bracket end
        // is further from the estimate
        const wantProof = pf === null ? true : nf === null ? false : !(Math.abs(Math.log(est / pf)) <= Math.abs(Math.log(est / nf)));
        const up = provenBelow ? !wantProof : wantProof;    // aim above the root?
        let F = est * (up ? 1 + 0.4 * tol : 1 - 0.4 * tol);
        const lo = provenBelow ? pf : nf, hi = provenBelow ? nf : pf;
        if (!(F > 0) || !isFinite(F)) F = (pf !== null ? pf : nf) * (up ? 1.6 : 1 / 1.6);
        if (lo !== null && F <= lo) F = hi !== null ? Math.sqrt(lo * hi) : lo * (1 + tol);
        if (hi !== null && F >= hi) F = lo !== null ? Math.sqrt(lo * hi) : hi / (1 + tol);
        if (near(F)) F = lo !== null && hi !== null ? Math.sqrt(lo * hi) : F * (up ? 1 + tol : 1 / (1 + tol));
        F = Math.min(100, Math.max(0.01, F));
        const a = await fn(F); n++;
        if (Number.isNaN(a)) break;
        pts.push([F, a]);
        if (proven(a)) pf = F; else nf = F;
        if ((F === 100 || F === 0.01) && (pf === null || nf === null)) break;
      }
      return { F: pf, other: nf };
    };
    const aLb1 = res.lb.valid ? res.lb.alpha : -Infinity;
    const aUb1 = Math.min(res.dlo.valid && !axi ? res.dlo.alpha : Infinity, res.ub && res.ub.valid ? res.ub.alpha : Infinity);
    // lower bound on FS: largest F with alpha_LB(F) >= t (safe side is "alpha above t")
    let lo = await bracket(lbAt, 1, aLb1, (a) => a >= t, true);
    // upper bound on FS: smallest F with alpha_UB(F) <= t
    let hi = await bracket(ubAt, 1, aUb1, (a) => a <= t, false);
    // If the bracket is still wide (typically a slope with little or no footing load, whose mechanism differs
    // from the one the mesh was adapted to), re-run the adaptive analysis at the estimated factor of safety
    // and bracket the lower bound again on that mesh.
    if (opts.strength.remesh !== false && lo.F && hi.F && hi.F / lo.F > 1 + 3 * tol) {
      const Fm = Math.sqrt(lo.F * hi.F);
      say(`Factor of safety on strength: new mesh at F = ${Fm.toFixed(3)}`, 0.9); await tick();
      const layers = prepared.layers.map((l) => { const o = Object.assign({}, l); delete o.psi; return o; });
      const sub = Object.assign({}, opts, { strength: null, onProgress: null });
      const r2 = await analyse(highs, Object.assign({}, prepared.input, { layers: reducedModel({ layers: layers }, Fm).layers }), sub);
      const a2 = r2.lb && r2.lb.valid ? r2.lb.alpha : -Infinity, u2 = r2.best.alphaUb;
      evals.push({ kind: "LB", F: Fm, alpha: a2 }, { kind: "UB", F: Fm, alpha: u2 });
      if (u2 <= t && Fm < hi.F) hi = { F: Fm };
      const mesh2 = r2.mesh, fe2 = r2.model;
      const lb2 = await bracket(async (F) => {
        say(`Factor of safety on strength: lower bound at F = ${F.toFixed(3)}`, 0.95); await tick();
        const m = reducedModel(fe2, F / Fm), r = (axi ? axiLowerBound : lowerBound)(highs, m, mesh2, assignMaterials(m, mesh2), opts);
        const a = r.valid ? r.alpha : -Infinity;
        evals.push({ kind: "LB", F: F, alpha: a });
        return a;
      }, Fm, a2, (a) => a >= t, true);
      if (lb2.F && lb2.F > lo.F) lo = lb2;
    }
    const out = { fsLb: lo.F == null ? 0 : lo.F, fsUb: hi.F == null ? NaN : hi.F, target: t, evals: evals };
    if (hi.F == null) out.note = "The upper bound on the factor of safety on strength could not be bracketed.";
    return out;
  }

  // every caller of the finite element solvers gets the exact-cone solve first, with the HiGHS fallback
  lowerBound = conicFirst(lowerBound); upperBound = conicFirst(upperBound); upperBoundQ = conicFirst(upperBoundQ);
  axiLowerBound = conicFirst(axiLowerBound); axiUpperBound = conicFirst(axiUpperBound);

  // run one finite element solve from plain (cloned) data, as sent to a parallel worker
  function runJob(highs, job) {
    const m = job.model; m.geo = geometry(m);
    const t = Date.now(), r = api[job.fn](highs, m, job.mesh, assignMaterials(m, job.mesh), job.opts);
    if (r) { r.time = (Date.now() - t) / 1000; delete r.mesh; }
    return r;
  }
  const api = { analyse, runJob, axiLowerBound, axiUpperBound, prepare, geometry, onDomain, buildTriMesh, buildQuadMesh, initialCells, refineCells, assignMaterials,
    lowerBound, upperBound, upperBoundQ, upperBoundDLO, dloAdaptive, dloVelocity, prandtlExtent, classical, sigmaV, unitWeightAt, soilAt,
    porePressure, statusText, fixedDepths, columnWeight, setConic };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.FELA = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
