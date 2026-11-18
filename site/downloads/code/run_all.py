#!/usr/bin/env python3
"""Run every script in this folder and save each output as <script>.txt in results/.
Usage: python run_all.py            (run all)
       python run_all.py ringdown   (run scripts whose name contains 'ringdown')
Requires Python 3 with numpy, scipy, sympy, pandas and openpyxl."""
import os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
FIRST = ["build_relational_matrix.py", "build_antiparticle_and_decays.py"]
def main():
    res = os.path.join(HERE, "results"); os.makedirs(res, exist_ok=True)
    scripts = FIRST + sorted(f for f in os.listdir(HERE) if f.endswith(".py") and f not in FIRST + ["run_all.py"])
    if len(sys.argv) > 1:
        scripts = [s for s in scripts if sys.argv[1] in s]
    failed = []
    for s in scripts:
        t = time.time()
        r = subprocess.run([sys.executable, "-I", os.path.join(HERE, s)], cwd=HERE, capture_output=True, text=True)
        with open(os.path.join(res, s[:-3] + ".txt"), "w", encoding="utf-8") as f:
            f.write(r.stdout + (("\n--- errors ---\n" + r.stderr) if r.returncode else ""))
        if r.returncode: failed.append(s)
        print(f"{'ok' if r.returncode == 0 else 'FAILED':6} {time.time()-t:6.1f}s  {s}", flush=True)
    print(f"\n{len(scripts)-len(failed)} of {len(scripts)} scripts ran. Outputs are in results/.")
    sys.exit(1 if failed else 0)
if __name__ == "__main__":
    main()
