import math, json

SIG = 5.670374e-8
RE = 6378.137
MU = 398600.4418
J2 = 1.08263e-3
EPS_OBL = math.radians(23.4393)

# ---- parametros termicos (ilustrativos) ----
Cc = 300.0; Cs = 6000.0; G = 4.0
EMI = 0.85; AR = 0.5; POTH = 30.0
PIDLE = 20.0; PACT = 100.0; DUTY = 0.35
# ---- ambiente/atitude (ASSUMIDOS) ----
S0 = 1361.0
QIR = 237.0
ALPHA = 0.15
H = 500.0
A0 = RE + H
TORB = 2*math.pi*math.sqrt(A0**3/MU)
RHO = math.asin(RE/A0)
BSTAR = math.degrees(RHO)
F_NADIR = math.sin(RHO)**2
F_SIDE = (RHO - math.sin(RHO)*math.cos(RHO))/math.pi
OPD = 86400.0/TORB

def efrac(bdeg):
    cb = math.cos(math.radians(bdeg))
    if cb <= 1e-9: return 0.0
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
    ci = -rate/(1.5*J2*(RE/a)**2*n)
    return math.acos(ci)

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
        nvec = (math.sin(inc)*math.sin(om), -math.sin(inc)*math.cos(om), math.cos(inc))
        sb = s[0]*nvec[0] + s[1]*nvec[1] + s[2]*nvec[2]
        out.append(math.degrees(math.asin(max(-1.0, min(1.0, sb)))))
    return out

def fdyn(Tc, Ts, P, Q):
    return (P - G*(Tc-Ts))/Cc, (G*(Tc-Ts) + Q + POTH - EMI*SIG*AR*Ts**4)/Cs

def rk4(Tc, Ts, P, Q, dt):
    k1 = fdyn(Tc, Ts, P, Q)
    k2 = fdyn(Tc+0.5*dt*k1[0], Ts+0.5*dt*k1[1], P, Q)
    k3 = fdyn(Tc+0.5*dt*k2[0], Ts+0.5*dt*k2[1], P, Q)
    k4 = fdyn(Tc+dt*k3[0], Ts+dt*k3[1], P, Q)
    return (Tc+dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]), Ts+dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]))

def qext(att, bdeg, ph, ecl):
    cb = math.cos(math.radians(bdeg)); sb = math.sin(math.radians(bdeg))
    phi = 2*math.pi*ph
    if att == 'zenith':
        sun = 0.0 if ecl else ALPHA*S0*AR*max(0.0, cb*math.cos(phi))
        ir = 0.0
    elif att == 'side':
        sun = 0.0 if ecl else ALPHA*S0*AR*abs(sb)
        ir = EMI*F_SIDE*QIR*AR
    else:
        sun = 0.0
        ir = EMI*F_SIDE*QIR*AR
    return sun + ir

def wdist(a, b):
    d = abs(a-b) % 1.0
    return min(d, 1.0-d)

def load(pol, ph):
    if pol == 'const':
        return DUTY*PACT + (1-DUTY)*PIDLE
    c = 0.0 if pol == 'noon' else 0.5
    return PACT if wdist(ph, c) < DUTY/2 else PIDLE

AC = 2e-6; AS = 16e-6; FF = 0.90; LL = 19.80; HS = 0.18; EF2 = 0.65

def turning(v):
    pts = [v[0]]
    for i in range(1, len(v)-1):
        if (v[i]-v[i-1])*(v[i+1]-v[i]) < 0: pts.append(v[i])
    pts.append(v[-1]); return pts

def rainflow(seq):
    st = []; cy = []
    for pt in seq:
        st.append(pt)
        while len(st) >= 3:
            xr = abs(st[-1]-st[-2]); yr = abs(st[-2]-st[-3])
            if xr < yr: break
            if len(st) == 3:
                cy.append((yr, 0.5)); st.pop(0)
            else:
                cy.append((yr, 1.0)); del st[-3:-1]
    for i in range(len(st)-1): cy.append((abs(st[i+1]-st[i]), 0.5))
    return cy

def dmg_per_orbit(TC, TS, keep):
    X = [AC*u - AS*w for u, w in zip(TC, TS)]
    cyc = rainflow(turning(X))
    Tsj = sum((u+w)/2 for u, w in zip(TC, TS))/len(TC) - 273.15
    c = -0.442 - 6e-4*Tsj + 0.0174*math.log(1+OPD)
    D = 0.0
    for dX, n in cyc:
        dg = FF*LL/HS*dX
        if dg < 1e-6: continue
        lnN = math.log(0.5) + (1.0/c)*math.log(dg/EF2)
        D += n*math.exp(-max(min(lnN, 700.0), -700.0))
    return D/keep, Tsj

def sim_beta(bdeg, att, pol, norb=8, dt=20.0, keep=3):
    n = int(round(norb*TORB/dt)); dtt = norb*TORB/n
    fe = efrac(bdeg)
    Tc = 285.0; Ts = 285.0
    TC = []; TS = []; spo = n/norb
    for k in range(n):
        ph = ((k*dtt) % TORB)/TORB
        ecl = abs(ph-0.5) < fe/2
        Tc, Ts = rk4(Tc, Ts, load(pol, ph), qext(att, bdeg, ph, ecl), dtt)
        if k >= n - int(keep*spo):
            TC.append(Tc); TS.append(Ts)
    D, Tsj = dmg_per_orbit(TC, TS, keep)
    return D, dict(dTc=max(TC)-min(TC), dTs=max(TS)-min(TS), Tsj=Tsj, Tsmax=max(TS), Tcmax=max(TC))

GRID = [0, 15, 30, 45, 55, 60, 64, 66, 68, 70, 72, 75, 80, 85, 90]
ATTS = ['zenith', 'side', 'avoid']
POLS = ['const', 'noon', 'mid']

TAB = {}
EXTRA = {}
for att in ATTS:
    for pol in POLS:
        Dl = []
        for b in GRID:
            D, ex = sim_beta(b, att, pol)
            Dl.append(D)
            EXTRA[(att, pol, b)] = ex
        TAB[(att, pol)] = Dl

def interp(Dl, b):
    b = abs(b)
    for i in range(len(GRID)-1):
        if GRID[i] <= b <= GRID[i+1]:
            t = (b-GRID[i])/(GRID[i+1]-GRID[i])
            return Dl[i]*(1-t) + Dl[i+1]*t
    return Dl[-1]

inc_sso = sso_incl(A0)
ORBS = {
    'SSO_6h_dawn_dusk': (inc_sso, dict(ltan_h=6.0)),
    'SSO_10h30': (inc_sso, dict(ltan_h=10.5)),
    'i51.6_precessando': (math.radians(51.6), dict(om0=0.0)),
}
out = dict(T_orb_min=TORB/60, beta_star_deg=BSTAR, F_nadir=F_NADIR, F_side=F_SIDE,
           inc_sso_deg=math.degrees(inc_sso),
           raan_rate_sso_deg_dia=math.degrees(raan_rate(A0, inc_sso))*86400.0,
           raan_rate_51_6_deg_dia=math.degrees(raan_rate(A0, math.radians(51.6)))*86400.0,
           Qir_side_W=EMI*F_SIDE*QIR*AR, Qsun_max_W=ALPHA*S0*AR)
res = {}
for on, (inc, kw) in ORBS.items():
    bs = beta_series(inc, **kw)
    ab = [abs(b) for b in bs]
    ecl_days = sum(1 for b in ab if b < BSTAR)
    emax = max(efrac(b)*TORB/60 for b in ab)
    res[on] = dict(beta_min=min(ab), beta_max=max(ab), dias_com_eclipse=ecl_days,
                   eclipse_max_min=emax, beta_dia0=bs[0], beta_dia80=bs[80],
                   beta_dia172=bs[172], beta_dia355=bs[355])
    life = {}
    for att in ATTS:
        for pol in POLS:
            Dl = TAB[(att, pol)]
            Dy = sum(interp(Dl, b) for b in bs)*OPD
            life[att+'/'+pol] = (1.0/Dy if Dy > 1e-30 else 1e30)
    res[on]['vida_anos'] = life
out['orbitas'] = res
tabout = {}
for att in ATTS:
    for pol in POLS:
        tabout[att+'/'+pol] = {str(b): TAB[(att, pol)][i] for i, b in enumerate(GRID) if b in (0, 30, 60, 68, 75, 90)}
out['D_orbita'] = tabout
tt = {}
for key in [('zenith','const',0),('zenith','const',60),('side','const',60),('side','const',90),
            ('avoid','const',60),('zenith','noon',60),('zenith','mid',60)]:
    e = EXTRA[key]
    tt[str(key)] = dict(dTc=round(e['dTc'],2), dTs=round(e['dTs'],2), Tsj=round(e['Tsj'],1), Tsmax=round(e['Tsmax'],1))
out['temps'] = tt
cond = {}
k = EMI*SIG*AR
for att in ATTS:
    if att == 'zenith':
        qmin, qmax = 0.0, ALPHA*S0*AR
    elif att == 'side':
        qmin, qmax = 0.0, ALPHA*S0*AR + EMI*F_SIDE*QIR*AR
    else:
        qmin, qmax = EMI*F_SIDE*QIR*AR, EMI*F_SIDE*QIR*AR
    cond[att] = dict(Pidle_Cc=PIDLE/Cc, rhs=(qmax+POTH)/Cs, ok=PIDLE/Cc > (qmax+POTH)/Cs,
                     Tinv_K=((qmin+POTH)/k)**0.25)
out['invariancia'] = cond
result = {"json": json.dumps(out, ensure_ascii=False)}

if __name__ == "__main__":
    print(json.dumps(out, ensure_ascii=False, indent=1))
