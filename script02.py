import math
SIGMA = 5.670374e-8
Cc, Cs, G, eps, Ar, Poth = 300.0, 6000.0, 4.0, 0.85, 0.5, 30.0
Pidle, Qsun, Qir = 20.0, 70.0, 25.0
RE, MU, h0 = 6378.137, 398600.4418, 500.0
TORB = 2*math.pi*math.sqrt((RE+h0)**3/MU)
k = eps*SIGMA*Ar
cond = (Pidle/Cc, (Qsun+Qir+Poth)/Cs)            # 0.0667 > 0.0208 -> invariancia ok
Tinv = ((Qir+Poth)/k)**0.25                       # 218.57 K
r = 4*k*Tinv**3                                   # 1.0065 W/K
lam = (2*G + r - math.sqrt(r*r + 4*G*G))/2        # 0.4717 W/K
print(cond, cond[0] > cond[1], Tinv, r, lam, lam/Cs, math.exp(-lam/Cs*TORB))
