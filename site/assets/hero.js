"use strict";
(function () {
  const canvas = document.getElementById('field');
  if (!(canvas instanceof HTMLCanvasElement))
    return;
  const g = canvas.getContext('2d');
  if (!g)
    return;
  const c = canvas, ctx = g;
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const still = !!window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  let W = 0, H = 0;
  function size() {
    const r = c.getBoundingClientRect();
    W = r.width;
    H = r.height;
    c.width = W * dpr;
    c.height = H * dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }
  const N = 34, span = 4 * Math.PI;
  function frame(t) {
    ctx.clearRect(0, 0, W, H);
    const ph = t * 0.00018, cy = Math.cos(0.95), sy = Math.sin(0.95), s = Math.min(W, H) * 0.074;
    function P(i, j) {
      const x = -span / 2 + span * i / N, y = -span / 2 + span * j / N;
      const z = Math.cos(x + ph) + Math.cos(y - 0.6 * ph);
      const a = 0.62, xr = x * Math.cos(a) - y * Math.sin(a), yr = x * Math.sin(a) + y * Math.cos(a);
      return [W / 2 + xr * s, H * 0.52 + (yr * cy - z * 1.15 * sy) * s * 0.62, Math.cos(x + ph) * Math.cos(y - 0.6 * ph)];
    }
    const pts = [];
    for (let i = 0; i <= N; i++) {
      const row = [];
      for (let j = 0; j <= N; j++)
        row.push(P(i, j));
      pts.push(row);
    }
    ctx.lineWidth = 1.1;
    for (let i = 0; i <= N; i++)
      for (let j = 0; j <= N; j++) {
        const p = pts[i][j], k = p[2];
        const col = k > 0 ? 'rgba(245,197,24,' : 'rgba(111,168,220,', al = (0.18 + 0.6 * Math.abs(k)).toFixed(2);
        ctx.strokeStyle = col + al + ')';
        if (i < N) {
          const q = pts[i + 1][j];
          ctx.beginPath();
          ctx.moveTo(p[0], p[1]);
          ctx.lineTo(q[0], q[1]);
          ctx.stroke();
        }
        if (j < N) {
          const q2 = pts[i][j + 1];
          ctx.beginPath();
          ctx.moveTo(p[0], p[1]);
          ctx.lineTo(q2[0], q2[1]);
          ctx.stroke();
        }
      }
  }
  let run = true;
  function loop(t) {
    if (!run)
      return;
    frame(t);
    window.requestAnimationFrame(loop);
  }
  size();
  window.addEventListener('resize', function () {
    size();
    if (still)
      frame(0);
  });
  if (still)
    frame(0);
  else
    window.requestAnimationFrame(loop);
  document.addEventListener('visibilitychange', function () {
    run = !document.hidden;
    if (run && !still)
      window.requestAnimationFrame(loop);
  });
})();
