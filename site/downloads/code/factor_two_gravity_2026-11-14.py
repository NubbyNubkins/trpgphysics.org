"""Where the '2' in resistance n = 1 + 2GM/(r c^2) comes from (14 Nov 2026).
Light deflection at the Sun's edge = (1+gamma)/2 * 4GM/(R c^2); gamma = how much mass-energy stretches SPACE relative to how much it slows TIME (PPN; Eddington 1922, Will 1971)."""
import math
G,c,M,R=6.674e-11,2.998e8,1.989e30,6.957e8
full=4*G*M/(R*c*c)*206265
cases=[("scalar resistance (one number, couples to trace of energy; Nordstrom 1913)",None),
       ("time-only resistance (Newton/Soldner 1801)",0.0),
       ("time + space resistance, equal shares (tensor; Einstein 1915)",1.0)]
for name,g in cases:
    d=0.0 if g is None else (1+g)/2*full
    print("%-75s %.3f arcsec"%(name,d))
print("measured: 1.75 arcsec; Cassini (2003) gamma - 1 = (2.1 +/- 2.3) x 10^-5")
# split of the 2: time part and space part each contribute GM/(r c^2)
r=R; print("at the Sun's edge: time part %.2e + space part %.2e = %.2e = n - 1"%(G*M/(r*c*c),G*M/(r*c*c),2*G*M/(r*c*c)))
