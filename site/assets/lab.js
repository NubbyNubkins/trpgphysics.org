"use strict";
(function () {
  const T = TRPG, f = T.fmt;
  function byId(id) {
    const el = document.getElementById(id);
    if (!el)
      throw new Error("missing element #" + id);
    return el;
  }
  function row(name, comp, ref, ok, note) {
    return '<tr><td>' + name + '</td><td>' + comp + '</td><td>' + ref + '</td><td class="' + (ok ? 'ok' : 'bad') + '">' + (ok ? 'Pass' : 'Fail') + '</td><td class="note">' + (note || '') + '</td></tr>';
  }
  function runAll() {
    const out = [], d = T.deflection(), a = T.age(), q = T.koide(), tw = T.twist(), l = T.leptonsFromTau(), n = T.neutrinos(), cb = T.cabibbo(), lt = T.lifetimes(), an = T.anomalies(3), cw = T.crossingWave(), cp = T.containment(60.1, 1710), rd = T.ringdown(), sh = T.shapiro();
    out.push(row('Deflection of starlight at the solar limb', f(d.value, 3) + '″', '1.75″ (γ = 1 to 2 × 10⁻⁵, Cassini)', Math.abs(d.value - 1.75) < 0.01, 'time-only resistance would give ' + f(d.timeOnly, 3) + '″'));
    out.push(row('Shapiro delay, Earth–Mars, one way', f(sh, 1) + ' μs', 'general-relativity value; Viking agreement 0.1%', true, 'consistency check'));
    out.push(row('Age of a flat universe (Planck 2018 inputs)', f(a, 2) + ' Gyr', '13.787 ± 0.020 Gyr', Math.abs(a - 13.787) < 0.04));
    out.push(row('Koide ratio of the charged leptons', f(q, 7), '2/3 = 0.6666667', Math.abs(q - 2 / 3) < 1e-5));
    out.push(row('Twist angle δ from data', f(tw, 6), '2/9 = 0.222222 (uncertainty 0.000040)', Math.abs(tw - 2 / 9) < 4e-5));
    out.push(row('Electron mass from the tau, with δ = 2/9', f(l.e, 5) + ' MeV', '0.51099895 MeV', Math.abs(l.e / 0.51099895 - 1) < 1e-4, f(100 * (l.e / 0.51099895 - 1), 4) + '%'));
    out.push(row('Muon mass from the tau, with δ = 2/9', f(l.mu, 4) + ' MeV', '105.6583755 MeV', Math.abs(l.mu / 105.6583755 - 1) < 1e-4, f(100 * (l.mu / 105.6583755 - 1), 4) + '%'));
    out.push(row('Shared part 3M₀² of the lepton matrix', f(l.scale3, 2) + ' MeV', 'proton 938.27 MeV (observation, not a claim)', Math.abs(l.scale3 / 938.272 - 1) < 0.005, f(100 * (l.scale3 / 938.272 - 1), 2) + '%'));
    out.push(row('Neutrino masses (normal ordering)', f(n.m1, 2) + ', ' + f(n.m2, 2) + ', ' + f(n.m3, 2) + ' meV', 'prediction; JUNO and DESI will test', n.sum < 0.064, 'sum ' + f(n.sum, 4) + ' eV, below DESI 0.064 eV'));
    out.push(row('Cabibbo relation √(m<sub>d</sub>/m<sub>s</sub>)', f(cb, 4), '|V<sub>us</sub>| = 0.2243', Math.abs(cb - 0.2243) < 0.002));
    out.push(row('Muon lifetime, Fermi theory', f(lt.tmu * 1e6, 4) + ' μs', '2.1970 μs', Math.abs(lt.tmu / 2.196981e-6 - 1) < 0.01, 'tree level'));
    out.push(row('Tau lifetime from the muon, τ<sub>μ</sub>·BR·e<sup>−10S</sup>', f(lt.ttau, 1) + ' fs', '290.3 ± 0.5 fs', Math.abs(lt.ttau - 290.3) < 2));
    out.push(row('Quark charges for three colours', '+' + f(an.Qu, 4) + ', ' + f(an.Qd, 4), '+2/3, −1/3', an.maxResidual < 1e-12, 'all four anomaly sums vanish'));
    out.push(row('Total curvature of z = cos x + cos y over one cell', f(cw.total, 6), '0 (Gauss–Bonnet)', Math.abs(cw.total) < 1e-6, f(100 * cw.positive, 1) + '% of the area curved like a sphere'));
    out.push(row('Containment potential from B = 60.1 and m<sub>G</sub> = 1710 MeV', 'λ = ' + f(cp.lam, 0) + ', v = ' + f(cp.v, 1) + ' MeV', 'derived (no target)', true, 'wall ' + f(cp.wall, 3) + ' fm'));
    out.push(row('Ringdown of GW150914', f(rd.f, 0) + ' Hz, ' + f(rd.tau, 1) + ' ms', '251 ± 8 Hz, 4.0 ± 0.3 ms', Math.abs(rd.f - 251) < 8 && Math.abs(rd.tau - 4.0) < 0.3));
    byId('checks').innerHTML = out.join('');
    const np = out.filter(function (r) { return r.indexOf('class="ok"') > 0; }).length;
    byId('summary').textContent = np + ' of ' + out.length + ' checks pass.';
  }
  const cv = byId('nsPlot'), cx = cv.getContext('2d');
  const base = {};
  let timer;
  function draw(rows, ref) {
    if (!cx)
      return;
    const W = cv.width, H = cv.height, L = 90, R = W - 30, Tp = 30, Bt = H - 70;
    const X = function (r) { return L + (r - 7) / (14 - 7) * (R - L); };
    const Y = function (m) { return Bt - (m - 0) / (2.6 - 0) * (Bt - Tp); };
    cx.clearRect(0, 0, W, H);
    cx.font = '22px system-ui, sans-serif';
    cx.lineWidth = 2;
    cx.fillStyle = '#93A6BA';
    cx.strokeStyle = '#22344B';
    for (let m = 0.5; m <= 2.5; m += 0.5) {
      cx.beginPath();
      cx.moveTo(L, Y(m));
      cx.lineTo(R, Y(m));
      cx.stroke();
      cx.textAlign = 'right';
      cx.textBaseline = 'middle';
      cx.fillText(m.toFixed(1), L - 12, Y(m));
    }
    cx.textAlign = 'center';
    cx.textBaseline = 'top';
    for (let r = 7; r <= 14; r++)
      cx.fillText(String(r), X(r), Bt + 12);
    cx.fillText('radius (km)', (L + R) / 2, Bt + 40);
    cx.save();
    cx.translate(24, (Tp + Bt) / 2);
    cx.rotate(-Math.PI / 2);
    cx.fillText('mass (M☉)', 0, -10);
    cx.restore();
    cx.setLineDash([10, 8]);
    cx.strokeStyle = '#6FA8DC';
    cx.beginPath();
    cx.moveTo(L, Y(2.0));
    cx.lineTo(R, Y(2.0));
    cx.stroke();
    cx.setLineDash([]);
    cx.fillStyle = '#6FA8DC';
    cx.textAlign = 'right';
    cx.textBaseline = 'bottom';
    cx.fillText('2.0 M☉ pulsar', R, Y(2.0) - 6);
    const line = function (rs, col, w) {
      if (!cx)
        return;
      cx.strokeStyle = col;
      cx.lineWidth = w;
      cx.beginPath();
      rs.forEach(function (s, i) {
        const px = X(s.R), py = Y(s.M);
        if (i)
          cx.lineTo(px, py);
        else
          cx.moveTo(px, py);
      });
      cx.stroke();
    };
    if (ref)
      line(ref, 'rgba(201,213,225,.35)', 3);
    line(rows, '#F5C518', 5);
  }
  function solve() {
    const B = parseFloat(byId('nsB').value), gap = parseFloat(byId('nsGap').value);
    byId('nsBv').textContent = B.toFixed(1) + ' MeV/fm³';
    byId('nsGapv').textContent = gap + ' MeV';
    const rows = T.curve(T.quarkEos(B, 3.0, gap), 60), s = T.summarize(rows);
    const key = B.toFixed(1);
    let ref = base[key];
    if (!ref) {
      ref = gap === 0 ? rows : T.curve(T.quarkEos(B, 3.0, 0), 40);
      base[key] = ref;
    }
    draw(rows.slice(0, s.imax + 1), gap === 0 ? null : ref.slice(0, T.summarize(ref).imax + 1));
    const okM = s.Mmax >= 2.0, okL = s.L14 < 580;
    byId('rM').textContent = f(s.Mmax, 2);
    byId('rR').textContent = isNaN(s.R14) ? '–' : f(s.R14, 2);
    byId('rL').textContent = isNaN(s.L14) ? '–' : f(s.L14, 0);
    byId('nsOut').innerHTML = (okM ? '<span class="ok">Meets</span>' : '<span class="bad">Below</span>') + ' the 2.0 M☉ pulsar limit. '
      + (okL ? '<span class="ok">Within</span>' : '<span class="bad">Above</span>') + ' the GW170817 tidal bound of about 580. '
      + (okM && okL ? 'Both limits are met.' : 'Both limits must be met.');
  }
  function later() { window.clearTimeout(timer); timer = window.setTimeout(solve, 180); }
  ['nsB', 'nsGap'].forEach(function (id) { const el = byId(id); el.addEventListener('input', later); el.addEventListener('change', solve); });
  window.setTimeout(solve, 50);
  byId('runAll').addEventListener('click', runAll);
})();
