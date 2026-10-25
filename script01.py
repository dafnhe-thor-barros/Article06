import math, random, json

SIGMA = 5.670374e-8
RE = 6378.137
MU = 398600.4418
P0 = dict(h0=500.0, Cc=300.0, Cs=6000.0, G=4.0, eps=0.85, Ar=0.5, Poth=30.0,
          Pidle=20.0, Pact=100.0, duty=0.35, Qsun=70.0, Qir=25.0,
          ac=2e-6, asub=16e-6, F=0.90, L=19.80, hs=0.18, ef2=0.65)
A0 = RE + P0['h0']
TORB = 2*math.pi*math.sqrt(A0**3/MU)
FE = math.acos(math.sqrt(P0['h0']**2 + 2*RE*P0['h0'])/A0)/math.pi
QAVG = (1-FE)*(P0['Qsun']+P0['Qir']) + FE*P0['Qir']
S_SUN = 0.10
S_ECL = (1-FE+0.01) % 1.0

def load(ph, start, q, mode):
    if mode == 'const':
        return q['duty']*q['Pact'] + (1-q['duty'])*q['Pidle']
    return q['Pact'] if ((ph-start) % 1.0) < q['duty'] else q['Pidle']

def qext(ph, q, scen):
    if scen == 'sso':
        return QAVG
    return (q['Qsun']+q['Qir']) if ph < (1-FE) else q['Qir']

def fdyn(Tc, Ts, P, Q, q):
    dTc = (P - q['G']*(Tc-Ts))/q['Cc']
    dTs = (q['G']*(Tc-Ts) + Q + q['Poth'] - q['eps']*SIGMA*q['Ar']*Ts**4)/q['Cs']
    return dTc, dTs

def rk4(Tc, Ts, P, Q, dt, q):
    k1 = fdyn(Tc, Ts, P, Q, q)
    k2 = fdyn(Tc+0.5*dt*k1[0], Ts+0.5*dt*k1[1], P, Q, q)
    k3 = fdyn(Tc+0.5*dt*k2[0], Ts+0.5*dt*k2[1], P, Q, q)
    k4 = fdyn(Tc+dt*k3[0], Ts+dt*k3[1], P, Q, q)
    return (Tc+dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),
            Ts+dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]))

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

def life(TC, TS, keep, q, gate=0.0):
    X = [q['ac']*u - q['asub']*w for u, w in zip(TC, TS)]
    cyc = rainflow(turning(X))
    Tsj = sum((u+w)/2 for u, w in zip(TC, TS))/len(TC) - 273.15
    c = -0.442 - 6e-4*Tsj + 0.0174*math.log(1+86400.0/TORB)
    D = 0.0
    gmax = 0.0
    for dX, n in cyc:
        if dX < gate:
            continue
        dg = q['F']*q['L']/q['hs']*dX
        if dg < 1e-6:
            continue
        lnN = math.log(0.5) + (1.0/c)*math.log(dg/q['ef2'])
        D += n*math.exp(-max(min(lnN, 700.0), -700.0))
        gmax = max(gmax, dg)
    D = D/keep
    yrs = (1.0/D*TORB/86400.0/365.25) if D > 1e-30 else 1e30
    return yrs, c, Tsj, gmax

def simulate(q, start, mode, scen, norb=7, dt=20.0, keep=3, x0=(285.0, 285.0)):
    n = int(round(norb*TORB/dt))
    dtt = norb*TORB/n
    Tc, Ts = x0
    TC = []
    TS = []
    orbend = []
    spo = n/norb
    for k in range(n):
        ph = ((k*dtt) % TORB)/TORB
        Tc, Ts = rk4(Tc, Ts, load(ph, start, q, mode), qext(ph, q, scen), dtt, q)
        if k >= n - int(keep*spo):
            TC.append(Tc)
            TS.append(Ts)
        if (k+1) % int(round(spo)) == 0:
            orbend.append((Tc, Ts))
    return TC, TS, orbend

def wnorm(a, b, q):
    return math.sqrt(q['Cc']*a*a + q['Cs']*b*b)

def rank(v):
    idx = sorted(range(len(v)), key=lambda i: v[i])
    r = [0]*len(v)
    for k, i in enumerate(idx):
        r[i] = k
    return r

def spearman(x, y):
    rx = rank(x)
    ry = rank(y)
    n = len(x)
    mx = sum(rx)/n
    my = sum(ry)/n
    cxy = sum((a-mx)*(b-my) for a, b in zip(rx, ry))
    sx = math.sqrt(sum((a-mx)**2 for a in rx))
    sy = math.sqrt(sum((b-my)**2 for b in ry))
    return cxy/(sx*sy)

def run_ekf(truthq, bias, sc, qc, ss, qs, seed, gate, dtf=10.0, norb=7, keep=3):
    rng = random.Random(seed)
    n = int(round(norb*TORB/dtf))
    dtt = norb*TORB/n
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
        Tc, Ts = rk4(Tc, Ts, load(ph, S_ECL, truthq, 'pulse'), qext(ph, truthq, 'leo'), dtt, truthq)
        d = fdyn(xh[0], xh[1], load(ph, S_ECL, P0, 'pulse'), qext(ph, P0, 'leo')*(1+bias), P0)
        xh = [xh[0]+dtt*d[0], xh[1]+dtt*d[1]]
        r = 4*P0['eps']*SIGMA*P0['Ar']*xh[1]**3
        Jm = [[-P0['G']/P0['Cc'], P0['G']/P0['Cc']],
              [P0['G']/P0['Cs'], (-P0['G']-r)/P0['Cs']]]
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
        if k >= n - int(keep*n/norb):
            TCt.append(Tc)
            TSt.append(Ts)
            TCe.append(xh[0])
            TSe.append(xh[1])
            e2 += (xh[1]-Ts)**2
            cnt += 1
    Lt = life(TCt, TSt, keep, P0, gate)[0]
    Le = life(TCe, TSe, keep, P0, gate)[0]
    return Lt, Le, math.sqrt(e2/cnt)

def main():
    out = {}
    out['orbita'] = dict(T_min=TORB/60, fe=FE, eclipse_min=FE*TORB/60,
                         beta_star_deg=math.degrees(math.asin(RE/A0)),
                         orbitas_por_dia=86400/TORB, Qavg_W=QAVG)
    # A) reproducao aritmetica das Tabelas 7 e 8 do manuscrito original
    k = P0['F']*(14e-6)*P0['L']/P0['hs']
    c = -0.442 - 6e-4*25 + 0.0174*math.log(1+15.2)
    rep = {}
    for nm, dT, e in (("LEO_maior", 65, 2.887), ("LEO_menor", 65, 0.65),
                      ("LEO_TEMS", 20, 2.887), ("SSO_TEMS", 7.5, 2.887),
                      ("SSO_menor", 7.5, 0.65)):
        rep[nm] = (k*dT, 0.5*(k*dT/e)**(1/c))
    out['reproducao_tab7'] = dict(k=k, c=c, casos=rep,
                                  fator_m20=0.8**(1/c), fator_p20=1.2**(1/c),
                                  razao_2887_065=(2.887/0.65)**(-1/c))
    # B) cenarios
    cen = {}
    for nm, (mode, scen, st) in {"LEO_sol": ('pulse', 'leo', S_SUN),
                                 "LEO_eclipse": ('pulse', 'leo', S_ECL),
                                 "LEO_const": ('const', 'leo', 0.0),
                                 "SSO_pulsado": ('pulse', 'sso', S_SUN),
                                 "SSO_const": ('const', 'sso', 0.0)}.items():
        TC, TS, oe = simulate(P0, st, mode, scen, dt=10.0)
        yrs, cc, Tsj, gm = life(TC, TS, 3, P0)
        cen[nm] = dict(dTc=max(TC)-min(TC), dTs=max(TS)-min(TS), dgam_pct=100*gm,
                       Tsj=Tsj, c=cc, vida_anos=yrs,
                       erro_periodicidade=max(abs(oe[-1][0]-oe[-2][0]), abs(oe[-1][1]-oe[-2][1])))
    out['cenarios'] = cen
    # C) contracao
    n = int(7*TORB/10.0)
    dtt = 7*TORB/n
    xa = [285.0, 285.0]
    xb = [320.0, 255.0]
    d0 = wnorm(xa[0]-xb[0], xa[1]-xb[1], P0)
    dec = []
    for kk in range(n):
        ph = ((kk*dtt) % TORB)/TORB
        P = load(ph, S_ECL, P0, 'pulse')
        Q = qext(ph, P0, 'leo')
        xa = list(rk4(xa[0], xa[1], P, Q, dtt, P0))
        xb = list(rk4(xb[0], xb[1], P, Q, dtt, P0))
        if (kk+1) % int(round(n/7)) == 0:
            dec.append(wnorm(xa[0]-xb[0], xa[1]-xb[1], P0)/d0)
    cont = dict(razao_norma_por_orbita=dec)
    for Tmin in (230.0, 274.0):
        r = 4*P0['eps']*SIGMA*P0['Ar']*Tmin**3
        lam = (2*P0['G']+r-math.sqrt(r*r+4*P0['G']**2))/2
        cont['Tmin_%d' % int(Tmin)] = dict(lam_W_K=lam, taxa_1_s=lam/P0['Cs'],
                                           fator_por_orbita=math.exp(-lam/P0['Cs']*TORB))
    out['contracao'] = cont
    # D) fase da carga
    ph = {}
    for dm in (-60, -45, -30, -20, -10, -5, 0, 5, 10, 20, 30, 45, 60):
        st0 = (S_ECL + dm*60.0/TORB) % 1.0
        TC, TS, oe = simulate(P0, st0, 'pulse', 'leo', dt=10.0)
        ph[str(dm)] = life(TC, TS, 3, P0)[0]
    out['fase_min_vida_anos'] = ph
    # E) Monte Carlo
    random.seed(2026)
    sig = dict(Cs=0.15, G=0.25, eps=0.05, Qsun=0.06, Qir=0.20, Poth=0.20,
               ef2=0.30, hs=0.08, F=0.05)
    NMC = 150
    mc = {}
    for nm, (mode, scen, st) in {"LEO_sol": ('pulse', 'leo', S_SUN),
                                 "LEO_eclipse": ('pulse', 'leo', S_ECL),
                                 "SSO_pulsado": ('pulse', 'sso', S_SUN)}.items():
        samples = {kk: [] for kk in sig}
        lv = []
        for i in range(NMC):
            q = dict(P0)
            for kk, s in sig.items():
                z = random.gauss(0, s)
                q[kk] = P0[kk]*math.exp(z)
                samples[kk].append(z)
            TC, TS, oe = simulate(q, st, mode, scen, norb=7, dt=25.0)
            lv.append(life(TC, TS, 3, q)[0])
        sl = sorted(lv)
        ly = [math.log(v) for v in lv]
        mc[nm] = dict(P5=sl[int(0.05*NMC)], P50=sl[NMC//2], P95=sl[int(0.95*NMC)-1],
                      spearman={kk: spearman(samples[kk], ly) for kk in sig})
    out['monte_carlo'] = mc
    # F) EKF: decomposicao ruido x erro de modelo
    pert = dict(P0)
    pert['G'] = P0['G']*1.2
    pert['Cs'] = P0['Cs']*0.9
    pert['eps'] = P0['eps']*1.05
    casos = {"A_ruido_sem_erro_modelo": (P0, 0.0, 0.5, 0.5, 1.0, 1.0),
             "C_ruido_e_erro_modelo": (pert, 0.10, 0.5, 0.5, 1.0, 1.0),
             "D_so_diodo_ruido_e_erro": (pert, 0.10, 0.5, 0.5, None, None)}
    # caso B (erro de modelo sem ruido): adicionar (pert, 0.10, 1e-6, 1e-6, 1e-6, 1e-6)
    ekf = {}
    for cn, (tq, b, sc, qc, ss, qs) in casos.items():
        rs = [run_ekf(tq, b, sc, qc, ss, qs, sd, 0.0) for sd in (1, 2, 3, 4, 5)]
        rat = [r[1]/r[0] for r in rs]
        ekf[cn] = dict(vida_real=sum(r[0] for r in rs)/5, razao_media=sum(rat)/5,
                       razao_min=min(rat), razao_max=max(rat), rmse_Ts=sum(r[2] for r in rs)/5)
    out['ekf'] = ekf
    # G) faixa de Ts observada
    mins = {}
    for nm, q in (("nominal", P0), ("perturbado", pert)):
        Tc = 285.0
        Ts = 285.0
        m = 1e9
        M = -1e9
        n2 = int(round(7*TORB/10.0))
        dt2 = 7*TORB/n2
        for kk in range(n2):
            phh = ((kk*dt2) % TORB)/TORB
            Tc, Ts = rk4(Tc, Ts, load(phh, S_ECL, q, 'pulse'), qext(phh, q, 'leo'), dt2, q)
            if kk > n2//2:
                m = min(m, Ts)
                M = max(M, Ts)
        mins[nm] = (m, M)
    out['Ts_min_max_K'] = mins
    # H) sensibilidade a erro de medicao de DeltaT
    out['sens_medicao'] = {str(dT): {str(d): ((dT+d)/dT)**(1/c) for d in (-2, -1, 1, 2)}
                           for dT in (7.5, 20, 65)}
    return out

result = {"json": json.dumps(main(), ensure_ascii=False)}
