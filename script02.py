import math, json

SIGMA = 5.670374e-8
Cc, Cs, G = 300.0, 6000.0, 4.0
EPS, AR, POTH, PIDLE = 0.85, 0.5, 30.0, 20.0
RE, MU, H0 = 6378.137, 398600.4418, 500.0
S0, QIR_FLUX, ALPHA = 1361.0, 237.0, 0.15      # premissas (sem fonte verificada)
A0 = RE + H0
TORB = 2*math.pi*math.sqrt(A0**3/MU)
RHO = math.asin(RE/A0)
F_SIDE = (RHO - math.sin(RHO)*math.cos(RHO))/math.pi
K = EPS*SIGMA*AR

def contracao(qmin, qmax):
    """Proposicoes 1 e 2: invariancia de {Tc>=Ts>=Tinv} e taxa de contracao."""
    tinv = ((qmin + POTH)/K)**0.25
    r = 4*K*tinv**3
    lam = (2*G + r - math.sqrt(r*r + 4*G*G))/2
    lhs = PIDLE/Cc
    rhs = (qmax + POTH)/Cs
    return dict(qmin_W=qmin, qmax_W=qmax, cond_lhs=lhs, cond_rhs=rhs,
                cond_ok=lhs > rhs, Tinv_K=tinv, r_min_W_K=r, lambda_W_K=lam,
                taxa_1_s=lam/Cs, fator_por_orbita=math.exp(-lam/Cs*TORB))

q_sun_max = ALPHA*S0*AR
q_ir_side = EPS*F_SIDE*QIR_FLUX*AR
out = {
    "referencia_beta0": contracao(25.0, 95.0),          # Qir=25 W, Qsun+Qir=95 W
    "estendido_zenith": contracao(0.0, q_sun_max),
    "estendido_side": contracao(0.0, q_sun_max + q_ir_side),
    "TORB_s": TORB,
}
result = {"json": json.dumps(out, ensure_ascii=False)}

if __name__ == "__main__":
    print(json.dumps(out, ensure_ascii=False, indent=1))
