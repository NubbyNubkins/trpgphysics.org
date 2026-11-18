/* trpgphysics.org validation checks. Each check recomputes a published result from its inputs.
   The same calculations, in more detail, are in the downloadable Python scripts (MIT License). */
"use strict";
var TRPG = (function () {
  var PI = Math.PI;
  // Constants (CODATA 2018 / IAU) and inputs (PDG 2024, Planck 2018, NuFIT 6.0)
  var G = 6.67430e-11, c = 299792458, GMsun = 1.32712440018e20, Rsun = 6.957e8, Msun = 1.98847e30;
  var hc = 197.3269804, hbarGeVs = 6.582119569e-25;
  var me = 0.51099895, mmu = 105.6583755, mtau = 1776.93;

  function fmt(x, d) { return Number(x).toFixed(d); }

  // ---------- 1. Light deflection ----------
  function deflection() {
    var gammas = { tensor: 1, timeOnly: 0 };
    var rad = 4 * GMsun / (Rsun * c * c);
    var arc = rad * 180 / PI * 3600;
    return { value: arc, timeOnly: arc / 2 };
  }
  // ---------- 2. Shapiro delay (one way, Earth-Mars, ray grazing the Sun) ----------
  function shapiro() {
    var rE = 1.496e11, rM = 2.279e11, b = Rsun;
    return 2 * GMsun / Math.pow(c, 3) * Math.log(4 * rE * rM / (b * b)) * 1e6; // microseconds
  }
  // ---------- 3. Age of the universe ----------
  function age() {
    var H0 = 67.4, Om = 0.315, Or = 9.1e-5, OL = 1 - Om - Or;
    var H0s = H0 * 1000 / 3.0856775814913673e22; // s^-1
    // t = int_0^1 da/(a H) = int dx with a = x^2 substitution handled by fine Simpson in ln a
    var N = 20000, lo = Math.log(1e-12), hi = 0, h = (hi - lo) / N, s = 0;
    for (var i = 0; i <= N; i++) {
      var la = lo + i * h, a = Math.exp(la);
      var f = 1 / (H0s * Math.sqrt(Om / (a * a * a) + Or / (a * a * a * a) + OL));
      s += f * ((i === 0 || i === N) ? 1 : (i % 2 ? 4 : 2));
    }
    return s * h / 3 / (3.15576e7 * 1e9);
  }
  // ---------- 4-6. Charged leptons ----------
  function koide() {
    var s = Math.sqrt(me) + Math.sqrt(mmu) + Math.sqrt(mtau);
    return (me + mmu + mtau) / (s * s);
  }
  function twist() {
    var M0 = (Math.sqrt(me) + Math.sqrt(mmu) + Math.sqrt(mtau)) / 3;
    return Math.acos((Math.sqrt(mtau) / M0 - 1) / Math.SQRT2);
  }
  function leptonsFromTau() {
    var d = 2 / 9;
    var M0 = Math.sqrt(mtau) / (1 + Math.SQRT2 * Math.cos(d));
    var m = [0, 1, 2].map(function (k) { var r = M0 * (1 + Math.SQRT2 * Math.cos(d + 2 * PI * k / 3)); return r * r; });
    // k=0 tau, k=1 electron, k=2 muon
    return { e: m[1], mu: m[2], scale3: 3 * M0 * M0 };
  }
  // ---------- 7. Neutrinos ----------
  function neutrinos() {
    var d21 = 7.49e-5, d31 = 2.513e-3; // eV^2, normal ordering
    function g(m1) {
      var m2 = Math.sqrt(m1 * m1 + d21), m3 = Math.sqrt(m1 * m1 + d31);
      var s = -Math.sqrt(m1) + Math.sqrt(m2) + Math.sqrt(m3);
      return (m1 + m2 + m3) - (2 / 3) * s * s;
    }
    var a = 1e-7, b = 0.02, ga = g(a);
    for (var i = 0; i < 200; i++) { var mid = 0.5 * (a + b), gm = g(mid); if ((gm > 0) === (ga > 0)) { a = mid; ga = gm; } else b = mid; }
    var m1 = 0.5 * (a + b), m2 = Math.sqrt(m1 * m1 + d21), m3 = Math.sqrt(m1 * m1 + d31);
    return { m1: m1 * 1e3, m2: m2 * 1e3, m3: m3 * 1e3, sum: m1 + m2 + m3 };
  }
  // ---------- 8. Cabibbo ----------
  function cabibbo() { return Math.sqrt(4.70 / 93.5); }
  // ---------- 9. Lifetimes ----------
  function lifetimes() {
    var GF = 1.1663788e-5, m = mmu / 1000;
    var tmu = 192 * Math.pow(PI, 3) * hbarGeVs / (GF * GF * Math.pow(m, 5));
    var S = 0.5 * Math.log(mtau / mmu);
    var ttau = tmu * 0.1782 * Math.exp(-10 * S);
    return { tmu: tmu, ttau: ttau * 1e15 };
  }
  // ---------- 10. Charges from anomaly cancellation ----------
  function anomalies(N) {
    var YQ = 1 / (2 * N), Yu = YQ + 0.5, Yd = YQ - 0.5, YL = -0.5, Ye = -1;
    var c1 = 2 * YQ - Yu - Yd, c2 = N * YQ + YL, c3 = N * (2 * YQ - Yu - Yd) + 2 * YL - Ye;
    var c4 = N * (2 * Math.pow(YQ, 3) - Math.pow(Yu, 3) - Math.pow(Yd, 3)) + 2 * Math.pow(YL, 3) - Math.pow(Ye, 3);
    return { Qu: Yu, Qd: Yd, maxResidual: Math.max(Math.abs(c1), Math.abs(c2), Math.abs(c3), Math.abs(c4)) };
  }
  // ---------- 11. Crossing-wave geometry (Gauss-Bonnet) ----------
  function crossingWave() {
    var N = 600, h = 2 * PI / N, tot = 0, pos = 0;
    for (var i = 0; i < N; i++) for (var j = 0; j < N; j++) {
      var x = (i + 0.5) * h, y = (j + 0.5) * h;
      var zx = -Math.sin(x), zy = -Math.sin(y), zxx = -Math.cos(x), zyy = -Math.cos(y);
      var W = 1 + zx * zx + zy * zy;
      var K = zxx * zyy / (W * W);
      tot += K * Math.sqrt(W) * h * h;
      if (K > 0) pos++;
    }
    return { total: tot, positive: pos / (N * N) };
  }
  // ---------- 12. Containment potential ----------
  function containment(B, mG) {
    var B4 = B * Math.pow(hc, 3);
    var lam = Math.pow(mG, 4) / (16 * B4);
    return { lam: lam, v: mG / Math.sqrt(2 * lam), wall: hc / mG };
  }
  // ---------- 13. Ringdown ----------
  function ringdown() {
    var Msun_s = 4.925491e-6, M = 62 * 1.09, a = 0.67;
    var Mw = 1.5251 - 1.1568 * Math.pow(1 - a, 0.1292), Q = 0.7 + 1.4187 * Math.pow(1 - a, -0.499);
    var w = Mw / (M * Msun_s);
    return { f: w / (2 * PI), tau: 2 * Q / w * 1e3 };
  }

  // ---------- 14. Neutron stars: TOV + tidal Love number ----------
  var CONV = 1.3234e-6, MSUNKM = 1.4766;
  function k2of(C, y) {
    var t1 = 8 * Math.pow(C, 5) / 5 * Math.pow(1 - 2 * C, 2) * (2 + 2 * C * (y - 1) - y);
    var t2 = 2 * C * (6 - 3 * y + 3 * C * (5 * y - 8)) + 4 * Math.pow(C, 3) * (13 - 11 * y + C * (3 * y - 2) + 2 * C * C * (1 + y))
      + 3 * Math.pow(1 - 2 * C, 2) * (2 - y + 2 * C * (y - 1)) * Math.log(1 - 2 * C);
    return t1 / t2;
  }
  // EOS given as monotone tables P (MeV/fm^3) vs E (MeV/fm^3)
  function makeEos(P, E) {
    function interp(xs, ys, x) {
      if (x <= xs[0]) return ys[0];
      var n = xs.length; if (x >= xs[n - 1]) return ys[n - 1];
      var lo = 0, hi = n - 1;
      while (hi - lo > 1) { var m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
      var t = (x - xs[lo]) / (xs[hi] - xs[lo]); return ys[lo] + t * (ys[hi] - ys[lo]);
    }
    var es = E[0];
    return {
      es: es, Emax: E[E.length - 1],
      p: function (e) { return Math.max(interp(E, P, e), 0); },
      e: function (p) { return p > 0 ? interp(P, E, p) : es; },
      dpde: function (e) { var h = 0.5, a = Math.max(e - h, es); return Math.max((this.p(e + h) - this.p(a)) / (e + h - a), 1e-4); }
    };
  }
  function star(eos, ec) {
    var psurf = 1e-9, r = 1e-4, pc = eos.p(ec);
    var Y = [pc, 4 / 3 * PI * r * r * r * ec * CONV, 2.0];
    function rhs(r, Y) {
      var p = Math.max(Y[0], psurf), m = Y[1], y = Y[2];
      var eraw = eos.e(p), e = eraw * CONV, P = p * CONV, f = 1 - 2 * m / r, cs2 = eos.dpde(eraw);
      var A = m + 4 * PI * r * r * r * P;
      var dp = -(e + P) * A / (r * r * f);
      var Q = 4 * PI / f * (5 * e + 9 * P + (e + P) / cs2) - 6 / (r * r * f) - Math.pow(2 * A / (r * r * f), 2);
      var F = (1 - 4 * PI * r * r * (e - P)) / f;
      return [dp / CONV, 4 * PI * r * r * e, -(y * y + y * F + r * r * Q) / r];
    }
    var h = 0.004;
    for (var it = 0; it < 100000; it++) {
      var k1 = rhs(r, Y);
      var Y2 = Y.map(function (v, i) { return v + h / 2 * k1[i]; }), k2 = rhs(r + h / 2, Y2);
      var Y3 = Y.map(function (v, i) { return v + h / 2 * k2[i]; }), k3 = rhs(r + h / 2, Y3);
      var Y4 = Y.map(function (v, i) { return v + h * k3[i]; }), k4 = rhs(r + h, Y4);
      var Yn = Y.map(function (v, i) { return v + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]); });
      if (Yn[0] <= psurf) { // linear step to the surface
        var t = (Y[0] - psurf) / (Y[0] - Yn[0]);
        Y = Y.map(function (v, i) { return v + t * (Yn[i] - v); }); r += t * h; break;
      }
      Y = Yn; r += h;
    }
    var R = r, M = Y[1], y = Y[2] - 4 * PI * R * R * R * eos.es * CONV / M;
    var C = M / R, k2 = k2of(C, y);
    return { M: M / MSUNKM, R: R, L: 2 / 3 * k2 / Math.pow(C, 5) };
  }
  function curve(eos, n) {
    n = n || 60; var rows = [], e0 = eos.es * 1.01, e1 = Math.min(eos.Emax * 0.95, 3000);
    for (var i = 0; i < n; i++) rows.push(star(eos, e0 * Math.pow(e1 / e0, i / (n - 1))));
    return rows;
  }
  function summarize(rows) {
    var imax = 0; rows.forEach(function (s, i) { if (s.M > rows[imax].M) imax = i; });
    var L14 = NaN, R14 = NaN;
    for (var j = 1; j <= imax; j++) if (rows[j - 1].M <= 1.4 && rows[j].M >= 1.4) {
      var t = (1.4 - rows[j - 1].M) / (rows[j].M - rows[j - 1].M);
      L14 = Math.exp(Math.log(rows[j - 1].L) + t * (Math.log(rows[j].L) - Math.log(rows[j - 1].L)));
      R14 = rows[j - 1].R + t * (rows[j].R - rows[j - 1].R);
    }
    return { Mmax: rows[imax].M, Rmax: rows[imax].R, L14: L14, R14: R14, imax: imax };
  }
  function sequence(eos, n) { return summarize(curve(eos, n)); }
  function bagEos(B) {
    var P = [], E = [];
    for (var i = 0; i <= 4000; i++) { var e = 4 * B + i * (3000 - 4 * B) / 4000; E.push(e); P.push((e - 4 * B) / 3); }
    return makeEos(P, E);
  }
  function alphaS(Q) {
    function run(a, Q0, Q1, nf) {
      var b0 = (33 - 2 * nf) / (12 * PI), b1 = (153 - 19 * nf) / (24 * PI * PI);
      var t0 = Math.log(Q0 * Q0), t1 = Math.log(Q1 * Q1), n = 400, h = (t1 - t0) / n;
      function f(a) { return -(b0 * a * a + b1 * a * a * a); }
      for (var i = 0; i < n; i++) { var k1 = f(a), k2 = f(a + h * k1 / 2), k3 = f(a + h * k2 / 2), k4 = f(a + h * k3); a += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6; }
      return a;
    }
    var a = run(0.1180, 91.1876, 4.18, 5); a = run(a, 4.18, 1.27, 4); return run(a, 1.27, Q, 3);
  }
  function quarkEos(B, X, gap) {
    var ms = 100, n = 3000, mu = [], p = [];
    for (var i = 0; i < n; i++) {
      var m = 560 / X + i * (1200 - 560 / X) / (n - 1); mu.push(m);
      var a4 = 1 - 2 * alphaS(X * m / 1000) / PI;
      p.push(3 / (4 * PI * PI) * (a4 * Math.pow(m, 4) - (ms * ms - 4 * gap * gap) * m * m) / Math.pow(hc, 3) - B);
    }
    var nn = p.map(function (_, i) {
      if (i === 0) return (p[1] - p[0]) / (mu[1] - mu[0]);
      if (i === n - 1) return (p[n - 1] - p[n - 2]) / (mu[n - 1] - mu[n - 2]);
      return (p[i + 1] - p[i - 1]) / (mu[i + 1] - mu[i - 1]);
    });
    var e = mu.map(function (m, i) { return m * nn[i] - p[i]; });
    var i0 = 0; while (p[i0] <= 0) i0++;
    var es = e[i0 - 1] + (0 - p[i0 - 1]) * (e[i0] - e[i0 - 1]) / (p[i0] - p[i0 - 1]);
    var P = [0], E = [es];
    for (var j = i0; j < n; j++) { P.push(p[j]); E.push(e[j]); }
    return makeEos(P, E);
  }

  return {
    fmt: fmt, deflection: deflection, shapiro: shapiro, age: age, koide: koide, twist: twist,
    leptonsFromTau: leptonsFromTau, neutrinos: neutrinos, cabibbo: cabibbo, lifetimes: lifetimes,
    anomalies: anomalies, crossingWave: crossingWave, containment: containment, ringdown: ringdown,
    sequence: sequence, curve: curve, summarize: summarize, bagEos: bagEos, quarkEos: quarkEos
  };
})();
if (typeof module !== "undefined") module.exports = TRPG;
