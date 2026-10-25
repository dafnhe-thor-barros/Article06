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
BSTAR = math.degrees(RHO)
F_NADIR = math.sin(RHO)**2
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

def table_grid(q, grid, atts, pols, norb=8, dt=20.0):
    T = {}
    for att in atts:
        for pol in pols:
            Dl = []
            Tl = []
            Xl = []
            for b in grid:
                D, Tsj, Xm, dTs, dTc = sim_beta(q, b, att, pol, norb, dt)
                Dl.append(D)
                Tl.append(Tsj)
                Xl.append(Xm)
            T[(att, pol)] = (Dl, Tl, Xl)
    return T

def slow_damage(xs, tsj, q):
    i0 = 0
    for i in range(len(xs)):
        if xs[i] > xs[i0]:
            i0 = i
    seq = xs[i0:] + xs[:i0] + [xs[i0]]
    cyc = rainflow(turning(seq))
    c = -0.442 - 6e-4*tsj + 0.0174*math.log(1+1.0/365.25)
    D = 0.0
    for dX, n in cyc:
        dg = q['F']*q['L']/q['hs']*dX
        if dg < 1e-9:
            continue
        lnN = math.log(0.5) + (1.0/c)*math.log(dg/q['ef2'])
        D += n*math.exp(-max(min(lnN, 700.0), -700.0))
    return D

def wn(q, a, b):
    return math.sqrt(q['Cc']*a*a + q['Cs']*b*b)

def contr_obs(q, att, pol, bdeg, norb=7, dt=20.0):
    n = int(round(norb*TORB/dt))
    dtt = norb*TORB/n
    spo = int(round(n/norb))
    fe = efrac(bdeg)
    xa = [285.0, 285.0]
    xb = [320.0, 255.0]
    dprev = wn(q, xa[0]-xb[0], xa[1]-xb[1])
    rat = []
    for k in range(n):
        ph = ((k*dtt) % TORB)/TORB
        ecl = abs(ph-0.5) < fe/2
        Pw = load(q, pol, ph)
        Qw = qext(q, att, bdeg, ph, ecl)
        xa = list(rk4(q, xa[0], xa[1], Pw, Qw, dtt))
        xb = list(rk4(q, xb[0], xb[1], Pw, Qw, dtt))
        if (k+1) % spo == 0:
            d = wn(q, xa[0]-xb[0], xa[1]-xb[1])
            if dprev > 1e-9:
                rat.append(d/dprev)
            dprev = d
    return rat

def vf_num(rho, mode, n=300):
    s = 0.0
    dth = rho/n
    if mode == 'nadir':
        for i in range(n):
            th = (i+0.5)*dth
            s += math.cos(th)*math.sin(th)*dth
        return 2*math.pi*s/math.pi
    dph = math.pi/n
    for i in range(n):
        th = (i+0.5)*dth
        for j in range(n):
            ph = -math.pi/2 + (j+0.5)*dph
            s += math.sin(th)*(math.sin(th)*math.cos(ph))*dth*dph
    return s/math.pi

def parte_a():
    q = dict(PAR)
    T = table_grid(q, GRID, ATTS, POLS)
    inc = sso_incl(A0)
    out = dict(T_orb_min=TORB/60, beta_star=BSTAR, F_nadir=F_NADIR, F_side=F_SIDE,
               inc_sso=math.degrees(inc), opd=OPD,
               raan_sso=math.degrees(raan_rate(A0, inc))*86400.0,
               F_nadir_num=vf_num(RHO, 'nadir'), F_side_num=vf_num(RHO, 'lateral'))
    orb = {}
    for on, (ic, kw) in orbitas().items():
        bs = beta_series(ic, **kw)
        ab = [abs(b) for b in bs]
        life = {}
        slow = {}
        tsj = {}
        for att in ATTS:
            for pol in POLS:
                Dl, Tl, Xl = T[(att, pol)]
                Dy = sum(interp(GRID, Dl, b) for b in ab)*OPD
                key = att + '/' + pol
                life[key] = round(1.0/Dy, 4) if Dy > 1e-30 else 1e30
                xs = [interp(GRID, Xl, b) for b in ab]
                ts = sum(interp(GRID, Tl, b) for b in ab)/len(ab)
                sd = slow_damage(xs, ts, q)
                slow[key] = round(sd/Dy, 6) if Dy > 1e-30 else 0.0
                tsj[key] = round(ts, 1)
        orb[on] = dict(beta_min=round(min(ab), 2), beta_max=round(max(ab), 2),
                       dias_eclipse=sum(1 for b in ab if b < BSTAR),
                       eclipse_max_min=round(max(efrac(b)*TORB/60 for b in ab), 2),
                       vida_anos=life, razao_dano_lento_sobre_orbital=slow, Tsj_medio=tsj)
    out['orbitas'] = orb
    dz = {}
    for pol in POLS:
        Dl = T[('zenith', pol)][0]
        dz[pol] = {str(b): Dl[GRID.index(b)] for b in (0, 30, 60, 68, 75, 90)}
    out['D_zenith'] = dz
    co = {}
    for (att, pol, b) in (('zenith', 'mid', 20), ('zenith', 'noon', 20), ('lateral', 'mid', 20),
                          ('lateral', 'noon', 20), ('zenith', 'const', 70), ('lateral', 'mid', 70)):
        co[att + '/' + pol + '/b' + str(b)] = [round(v, 4) for v in contr_obs(q, att, pol, b)]
    out['contracao_obs'] = co
    return out, T

def parte_b(T):
    q = dict(PAR)
    res = {}
    casos = (('SSO_6h', 'zenith', 'mid'), ('SSO_6h', 'lateral', 'const'),
             ('SSO_10h30', 'lateral', 'const'), ('SSO_10h30', 'zenith', 'mid'),
             ('i51.6', 'zenith', 'const'), ('i51.6', 'lateral', 'mid'))
    for (on, att, pol) in casos:
        ic, kw = orbitas()[on]
        ab = [abs(b) for b in beta_series(ic, **kw)]
        Dl = T[(att, pol)][0]
        days = list(range(0, 365, 10))
        esc = 365.0/len(days)
        direto = sum(sim_beta(q, ab[d], att, pol)[0] for d in days)*esc*OPD
        i37 = sum(interp(GRID, Dl, ab[d]) for d in days)*esc*OPD
        i365 = sum(interp(GRID, Dl, b) for b in ab)*OPD
        res[on + '/' + att + '/' + pol] = dict(
            direto=direto, interp_mesmos_dias=i37, interp_365_dias=i365,
            erro_mesmos_dias=(i37/direto-1.0) if direto > 1e-30 else None,
            erro_365_dias=(i365/direto-1.0) if direto > 1e-30 else None)
    return res

def run_ekf(truthq, bias, sc, qc, ss, qs, seed, att, pol, bdeg, dtf=10.0, norb=7, keep=3):
    rng = random.Random(seed)
    n = int(round(norb*TORB/dtf))
    dtt = norb*TORB/n
    kst = n - int(keep*n/norb)
    fe = efrac(bdeg)
    q0 = PAR
    Tc = 285.0
    Ts = 285.0
    xh = [290.0, 280.0]
    Pm = [[25.0, 0.0], [0.0, 25.0]]
    TCt = []
    TSt = []
    TCe = []
    TSe = []
    e2 = 0.0
    cnt = 0
    for k in range(n):
        ph = ((k*dtt) % TORB)/TORB
        ecl = abs(ph-0.5) < fe/2
        Tc, Ts = rk4(truthq, Tc, Ts, load(truthq, pol, ph), qext(truthq, att, bdeg, ph, ecl), dtt)
        d = fdyn(q0, xh[0], xh[1], load(q0, pol, ph), qext(q0, att, bdeg, ph, ecl)*(1+bias))
        xh = [xh[0]+dtt*d[0], xh[1]+dtt*d[1]]
        r = 4*q0['emi']*SIG*q0['Ar']*xh[1]**3
        Jm = [[-q0['G']/q0['Cc'], q0['G']/q0['Cc']],
              [q0['G']/q0['Cs'], (-q0['G']-r)/q0['Cs']]]
        Fm = [[1+dtt*Jm[0][0], dtt*Jm[0][1]], [dtt*Jm[1][0], 1+dtt*Jm[1][1]]]
        Qn = [[1e-4*dtt, 0.0], [0.0, 4e-4*dtt]]
        FP = [[sum(Fm[i][u]*Pm[u][j] for u in range(2)) for j in range(2)] for i in range(2)]
        Pm = [[sum(FP[i][u]*Fm[j][u] for u in range(2)) + Qn[i][j] for j in range(2)] for i in range(2)]
        meas = [(0, round((Tc+rng.gauss(0, sc))/qc)*qc, sc**2+qc**2/12)]
        if ss is not None:
            meas.append((1, round((Ts+rng.gauss(0, ss))/qs)*qs, ss**2+qs**2/12))
        for idx, z, Rv in meas:
            S = Pm[idx][idx] + Rv
            K = [Pm[0][idx]/S, Pm[1][idx]/S]
            inn = z - xh[idx]
            xh = [xh[0]+K[0]*inn, xh[1]+K[1]*inn]
            Pm = [[Pm[i][j]-K[i]*Pm[idx][j] for j in range(2)] for i in range(2)]
        if k >= kst:
            TCt.append(Tc)
            TSt.append(Ts)
            TCe.append(xh[0])
            TSe.append(xh[1])
            e2 += (xh[1]-Ts)**2
            cnt += 1
    Dt = dmg(TCt, TSt, keep, q0)[0]
    De = dmg(TCe, TSe, keep, q0)[0]
    return Dt, De, math.sqrt(e2/cnt)

def parte_d():
    pert = dict(PAR)
    pert['G'] = PAR['G']*1.2
    pert['Cs'] = PAR['Cs']*0.9
    pert['emi'] = PAR['emi']*1.05
    casos = {"A_ruido_planta_nominal": (PAR, 0.0, 0.5, 0.5, 1.0, 1.0),
             "B_erro_modelo_sem_ruido": (pert, 0.10, 1e-6, 1e-6, 1e-6, 1e-6),
             "C_ruido_e_erro_modelo": (pert, 0.10, 0.5, 0.5, 1.0, 1.0),
             "D_so_diodo": (pert, 0.10, 0.5, 0.5, None, None)}
    res = {}
    for att, pol in (('zenith', 'noon'), ('zenith', 'mid'), ('lateral', 'mid')):
        for cn, (tq, b, sc, qc, ss, qs) in casos.items():
            rs = [run_ekf(tq, b, sc, qc, ss, qs, sd, att, pol, 20.0) for sd in (1, 2, 3, 4, 5)]
            rat = [r[0]/r[1] for r in rs if r[1] > 1e-30]
            res[att + '/' + pol + '/' + cn] = dict(
                razao_media=round(sum(rat)/len(rat), 3), razao_min=round(min(rat), 3),
                razao_max=round(max(rat), 3), rmse_Ts=round(sum(r[2] for r in rs)/5, 3))
    return res

outA, TT = parte_a()
outB = parte_b(TT)
outD = parte_d()
result = {"json": json.dumps(dict(A=outA, B=outB, D=outD), ensure_ascii=False)}
