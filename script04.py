import math, json, random, bisect

SIG = 5.670374e-8
RE = 6378.137
MU = 398600.4418
J2 = 1.08263e-3
EPS_OBL = math.radians(23.4393)
PAR = dict(Cc=300.0, Cs=6000.0, G=4.0, emi=0.85, Ar=0.5, Poth=30.0, Pidle=20.0, Pact=100.0,
           duty=0.35, S0=1361.0, Qir=237.0, alpha=0.15, ac=2e-6, asub=16e-6,
           F=0.90, L=19.80, hs=0.18, ef2=0.65)
H = 500.0
A0 = RE + H
TORB = 2*math.pi*math.sqrt(A0**3/MU)
RHO = math.asin(RE/A0)
F_SIDE = (RHO - math.sin(RHO)*math.cos(RHO))/math.pi
OPD = 86400.0/TORB
GRID = [0, 5, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35, 40, 45, 50, 55, 58, 60, 62, 64, 66,
        67, 67.5, 68, 68.5, 69, 70, 72, 75, 80, 85, 90]
ATTS = ['zenith', 'lateral']
POLS = ['const', 'noon', 'mid']

def efrac(bdeg):
    cb = math.cos(math.radians(bdeg))
    if cb <= 1e-9:
        return 0.0
    x = math.sqrt(H*H + 2*RE*H)/(A0*cb)
    return 0.0 if x >= 1.0 else math.acos(x)/math.pi

def sun_vec(d):
    M = 2*math.pi*(d-3.0)/365.2422
    lam = 2*math.pi*(d-79.5)/365.2422 + 0.0334*math.sin(M)
    return (math.cos(lam), math.sin(lam)*math.cos(EPS_OBL), math.sin(lam)*math.sin(EPS_OBL))

def sun_ra(d):
    s = sun_vec(d)
    return math.atan2(s[1], s[0])

def sso_incl(a):
    n = math.sqrt(MU/a**3)
    rate = 2*math.pi/(365.2422*86400.0)
    return math.acos(-rate/(1.5*J2*(RE/a)**2*n))

def raan_rate(a, inc):
    n = math.sqrt(MU/a**3)
    return -1.5*J2*(RE/a)**2*n*math.cos(inc)

def beta_series(inc, ltan_h=None, om0=0.0):
    rate = raan_rate(A0, inc)*86400.0
    if ltan_h is not None:
        om0 = sun_ra(80.0) + (ltan_h-12.0)*math.radians(15.0)
    out = []
    for d in range(365):
        om = om0 + rate*(d-80.0)
        s = sun_vec(d)
        nv = (math.sin(inc)*math.sin(om), -math.sin(inc)*math.cos(om), math.cos(inc))
        sb = s[0]*nv[0] + s[1]*nv[1] + s[2]*nv[2]
        out.append(math.degrees(math.asin(max(-1.0, min(1.0, sb)))))
    return out

def orbitas():
    inc = sso_incl(A0)
    return {'SSO_6h': (inc, dict(ltan_h=6.0)), 'SSO_10h30': (inc, dict(ltan_h=10.5)),
            'i51.6': (math.radians(51.6), dict(om0=0.0))}

def fdyn(q, Tc, Ts, Pw, Qw):
    return ((Pw - q['G']*(Tc-Ts))/q['Cc'],
            (q['G']*(Tc-Ts) + Qw + q['Poth'] - q['emi']*SIG*q['Ar']*Ts**4)/q['Cs'])

def rk4(q, Tc, Ts, Pw, Qw, dt):
    k1 = fdyn(q, Tc, Ts, Pw, Qw)
    k2 = fdyn(q, Tc+0.5*dt*k1[0], Ts+0.5*dt*k1[1], Pw, Qw)
    k3 = fdyn(q, Tc+0.5*dt*k2[0], Ts+0.5*dt*k2[1], Pw, Qw)
    k4 = fdyn(q, Tc+dt*k3[0], Ts+dt*k3[1], Pw, Qw)
    return (Tc+dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]), Ts+dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]))

def wdist(a, b):
    d = abs(a-b) % 1.0
    return min(d, 1.0-d)

def load(q, pol, ph):
    if pol == 'const':
        return q['duty']*q['Pact'] + (1-q['duty'])*q['Pidle']
    c = 0.0 if pol == 'noon' else 0.5
    return q['Pact'] if wdist(ph, c) < q['duty']/2 else q['Pidle']

def qext(q, att, bdeg, ph, ecl):
    cb = math.cos(math.radians(bdeg))
    sb = math.sin(math.radians(bdeg))
    phi = 2*math.pi*ph
    ae = 0.0 if ecl else 1.0
    if att == 'zenith':
        return ae*q['alpha']*q['S0']*q['Ar']*max(0.0, cb*math.cos(phi))
    return ae*q['alpha']*q['S0']*q['Ar']*abs(sb) + q['emi']*F_SIDE*q['Qir']*q['Ar']

def turning(v):
    pts = [v[0]]
    for i in range(1, len(v)-1):
        if (v[i]-v[i-1])*(v[i+1]-v[i]) < 0:
            pts.append(v[i])
    pts.append(v[-1])
    return pts

def rainflow(seq):
    st = []
    cy = []
    for pt in seq:
        st.append(pt)
        while len(st) >= 3:
            xr = abs(st[-1]-st[-2])
            yr = abs(st[-2]-st[-3])
            if xr < yr:
                break
            if len(st) == 3:
                cy.append((yr, 0.5))
                st.pop(0)
            else:
                cy.append((yr, 1.0))
                del st[-3:-1]
    for i in range(len(st)-1):
        cy.append((abs(st[i+1]-st[i]), 0.5))
    return cy

def dmg(TC, TS, keep, q):
    X = [q['ac']*u - q['asub']*w for u, w in zip(TC, TS)]
    cyc = rainflow(turning(X))
    Tsj = sum((u+w)/2 for u, w in zip(TC, TS))/len(TC) - 273.15
    c = -0.442 - 6e-4*Tsj + 0.0174*math.log(1+OPD)
    D = 0.0
    for dX, n in cyc:
        dg = q['F']*q['L']/q['hs']*dX
        if dg < 1e-6:
            continue
        lnN = math.log(0.5) + (1.0/c)*math.log(dg/q['ef2'])
        D += n*math.exp(-max(min(lnN, 700.0), -700.0))
    return D/keep, Tsj, sum(X)/len(X)

def sim_beta(q, bdeg, att, pol, norb=8, dt=20.0, keep=3):
    n = int(round(norb*TORB/dt))
    dtt = norb*TORB/n
    spo = n/norb
    kst = n - int(keep*spo)
    fe = efrac(bdeg)
    Tc = 285.0
    Ts = 285.0
    TC = []
    TS = []
    for k in range(n):
        ph = ((k*dtt) % TORB)/TORB
        ecl = abs(ph-0.5) < fe/2
        Tc, Ts = rk4(q, Tc, Ts, load(q, pol, ph), qext(q, att, bdeg, ph, ecl), dtt)
        if k >= kst:
            TC.append(Tc)
            TS.append(Ts)
    D, Tsj, Xm = dmg(TC, TS, keep, q)
    return D, Tsj, Xm, max(TS)-min(TS), max(TC)-min(TC)

def interp(grid, vals, b):
    b = abs(b)
    if b <= grid[0]:
        return vals[0]
    if b >= grid[-1]:
        return vals[-1]
    i = bisect.bisect_right(grid, b) - 1
    t = (b-grid[i])/(grid[i+1]-grid[i])
    return vals[i]*(1-t) + vals[i+1]*t

# ---------------- script04: robustez e validacao de grade ----------------
NMC = 16
SIGLN = dict(Cc=0.15, Cs=0.15, G=0.25, emi=0.05, Ar=0.05, alpha=0.30, S0=0.02, Qir=0.20, Poth=0.20)
GRID_MC = {'SSO_10h30': [10, 15, 20, 25, 30], 'SSO_6h': [58, 64, 68, 72, 80, 90]}

def abs_beta(on):
    ic, kw = orbitas()[on]
    return [abs(b) for b in beta_series(ic, **kw)]

def annual_life(grid, Dl, ab):
    Dy = sum(interp(grid, Dl, b) for b in ab)*OPD
    return 1.0/Dy if Dy > 1e-30 else 1e30

def table_sub(q, grid, combos, norb=8, dt=20.0):
    return {c: [sim_beta(q, b, c[0], c[1], norb, dt, 3)[0] for b in grid] for c in combos}

def validar():
    q = dict(PAR)
    out = {}
    combos = [('lateral', 'const'), ('zenith', 'mid')]
    for on in ('SSO_10h30', 'SSO_6h'):
        ab = abs_beta(on)
        lo, hi = (10, 30) if on == 'SSO_10h30' else (58, 90)
        gfine = [b for b in GRID if lo <= b <= hi]
        tf = table_sub(q, gfine, combos)
        tc = table_sub(q, GRID_MC[on], combos, norb=8, dt=30.0)
        for c in combos:
            Lf = annual_life(gfine, tf[c], ab)
            Lc = annual_life(GRID_MC[on], tc[c], ab)
            out[on + '/' + c[0] + '/' + c[1]] = dict(
                pontos_fina=len(gfine), vida_grade_fina=Lf, vida_grade_MC=Lc, dif_MC_vs_fina=Lc/Lf-1.0)
    return out

def mc_orbit(on, n, seed):
    rng = random.Random(seed)
    grid = GRID_MC[on]
    ab = abs_beta(on)
    combos = [(a, p) for a in ATTS for p in POLS]
    lives = []
    for i in range(n):
        q = dict(PAR)
        for k, s in SIGLN.items():
            q[k] = PAR[k]*math.exp(rng.gauss(0.0, s))
        tab = table_sub(q, grid, combos, norb=8, dt=30.0)
        lives.append({c: annual_life(grid, tab[c], ab) for c in combos})
    res = dict(N=n)
    for p in POLS:
        res['zenith_melhor_que_lateral/' + p] = sum(1 for L in lives if L[('zenith', p)] > L[('lateral', p)])
    for a in ATTS:
        res['mid_melhor_que_noon/' + a] = sum(1 for L in lives if L[(a, 'mid')] > L[(a, 'noon')])
        rr = sorted(L[(a, 'mid')]/L[(a, 'noon')] for L in lives)
        res['razao_mid_noon_min_med_max/' + a] = [round(rr[0], 3), round(rr[n//2], 3), round(rr[-1], 3)]
    best = {}
    for L in lives:
        cb = max(combos, key=lambda c: L[c])
        key = cb[0] + '/' + cb[1]
        best[key] = best.get(key, 0) + 1
    res['melhor_combinacao_contagem'] = best
    return res

out04 = dict(validacao=validar(), mc={'SSO_10h30': mc_orbit('SSO_10h30', NMC, 2026),
                                      'SSO_6h': mc_orbit('SSO_6h', NMC, 2027)})
result = {"json": json.dumps(out04, ensure_ascii=False)}
