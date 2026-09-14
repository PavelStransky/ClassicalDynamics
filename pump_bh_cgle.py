"""
Incoherently pumped Bose-Hubbard with two-body loss  ->  discrete CGLE.
Third companion to driven_bh.py / undriven_sitedep_bh.py.  Same conventions:
psi_i = b_i/sqrt(N), hbar_eff = 1/N, g = UN.

Lindblad
    H   = -J sum<ij>(b_i^d b_j + h.c.) - mu sum n_i + (U/2) sum n_i(n_i-1)
    L1  = sqrt(gamma_p) b_i^d      L2 = sqrt(gamma_l) b_i
    L3  = sqrt(Gamma_2) b_i^2                     (two-body loss)
    L4  = sqrt(Gamma_b) (b_i - b_j)               (bond-correlated loss)

Classical limit (U ~ g/N, Gamma_2 ~ Gam/2N, gamma's and Gamma_b O(1)):

    dpsi_j/dt = (d + iJ) lap_j + i Dt psi_j + (P/2) psi_j
                - (ig + Gam/2) |psi_j|^2 psi_j

    P   = gamma_p - gamma_l      net linear gain   (= -kappa)
    Gam = 2 N Gamma_2            two-body loss
    d   = Gamma_b / 2            real diffusion (bond loss only)
    Dt  = mu + Jz                pure gauge here (U(1) unbroken) -> 0

Canonical coordinates  psi_j = (q_j + i p_j)/sqrt2 :

    H = -J sum<jk>(q_j q_k + p_j p_k) - (Dt/2) sum (q_j^2+p_j^2)
        + (g/8) sum (q_j^2+p_j^2)^2

    q_j' = +dH/dp_j + [P/2 - (Gam/4) r_j^2] q_j + d (lap q)_j
    p_j' = -dH/dq_j + [P/2 - (Gam/4) r_j^2] p_j + d (lap p)_j     r_j^2 = q_j^2+p_j^2

i.e. a lattice of Stuart-Landau (Hopf) oscillators with conservative coupling.

Exact results checked below
---------------------------
  plane wave  psi_j = A e^{i(Qj - wt)} :   |A|^2 = (P - 2 d D_Q)/Gam,  D_Q = 2(1-cos Q)
  trapping annulus (d=0):  P/Gam <= S <= P L/Gam,  origin repels at rate P/2
  trace rule:   sum_k lam_k = (P - 2 d z) L - 2 Gam <S>     (z = coordination = 2)
  continuum BFN:  Re[(d+iJ)(Gam/2 - ig)] < 0   <=>   J g < -d Gam / 2
  ON A LATTICE that is necessary but NOT sufficient: waves with |Q| > pi/2 have
  effective mass J cos Q of the opposite sign and are stable at d = 0.  A real
  d is what closes that escape (it makes |A|^2 fall with Q).  At J=P=Gam=1,
  g=-1 no plane wave is stable above d* ~ 0.094, and only there is the attractor
  strange.
"""
import numpy as np

# ----------------------------------------------------------------- flow -----

def make(L, J, g, P, Gam, d=0.0, Dt=0.0):
    W = d + 1j*J
    def f(p):
        n = np.abs(p)**2
        return (W*(np.roll(p, -1) + np.roll(p, 1) - 2*p) + 1j*Dt*p + 0.5*P*p
                - (1j*g + 0.5*Gam)*n*p)
    def jac(p, D):
        pp = p[:, None] if D.ndim == 2 else p
        n = np.abs(pp)**2
        return (W*(np.roll(D, -1, axis=0) + np.roll(D, 1, axis=0) - 2*D)
                + 1j*Dt*D + 0.5*P*D
                - (1j*g + 0.5*Gam)*(2*n*D + pp**2*np.conj(D)))
    return f, jac


def rk4(f, p, dt):
    k1 = f(p); k2 = f(p+.5*dt*k1); k3 = f(p+.5*dt*k2); k4 = f(p+dt*k3)
    return p + dt/6*(k1 + 2*k2 + 2*k3 + k4)


def rk4_tan(f, jc, p, D, dt):
    k1 = f(p);           K1 = jc(p, D)
    k2 = f(p+.5*dt*k1);  K2 = jc(p+.5*dt*k1, D+.5*dt*K1)
    k3 = f(p+.5*dt*k2);  K3 = jc(p+.5*dt*k2, D+.5*dt*K2)
    k4 = f(p+dt*k3);     K4 = jc(p+dt*k3,    D+dt*K3)
    return p+dt/6*(k1+2*k2+2*k3+k4), D+dt/6*(K1+2*K2+2*K3+K4)


c2r = lambda D: np.vstack([D.real, D.imag])
r2c = lambda V: V[:V.shape[0]//2] + 1j*V[V.shape[0]//2:]

def to_qp(psi):  return np.sqrt(2)*psi.real, np.sqrt(2)*psi.imag
def to_psi(q, p): return (q + 1j*p)/np.sqrt(2)

def hamiltonian(q, p, J, g, Dt=0.0):
    r2 = q**2 + p**2
    return (-J*np.sum(q*np.roll(q, -1) + p*np.roll(p, -1))
            - 0.5*Dt*np.sum(r2) + (g/8)*np.sum(r2**2))

def kydim(l):
    l = np.sort(l)[::-1]; c = np.cumsum(l)
    if c[0] < 0: return 0.0
    k = int(np.max(np.where(c >= 0)[0]))
    return float(len(l)) if k == len(l)-1 else (k+1) + c[k]/abs(l[k+1])


# ------------------------------------------------ plane waves & stability ---

def pw_amp2(P, Gam, d, Q):
    return (P - 2*d*2*(1-np.cos(Q)))/Gam        # <= 0  =>  wave does not exist

def pw_growth(L, J, g, P, Gam, d, Q):
    """max_k Re lambda for the plane wave of wavenumber Q (analytic 2x2 blocks)"""
    A2 = pw_amp2(P, Gam, d, Q)
    if A2 <= 0: return None
    W = d + 1j*J; c = 1j*g*A2 + 0.5*Gam*A2
    k = 2*np.pi*np.arange(L)/L
    DQ = 2*(1-np.cos(Q))
    A11 = W*(DQ - 2*(1-np.cos(Q+k))) - c
    A22 = np.conj(W)*(DQ - 2*(1-np.cos(Q-k))) - np.conj(c)
    tr = A11 + A22; det = A11*A22 - c*np.conj(c)
    r = np.sqrt(tr**2 - 4*det + 0j)
    return max(np.max(((tr+r)/2).real), np.max(((tr-r)/2).real))

def pw_growth_jac(L, J, g, P, Gam, d, Q):
    A2 = pw_amp2(P, Gam, d, Q)
    if A2 <= 0: return None
    A = np.sqrt(A2)
    _, jc = make(L, J, g, P, Gam, d, Dt=g*A2 + J*2*(1-np.cos(Q)))
    p = A*np.exp(1j*Q*np.arange(L))
    return np.linalg.eigvals(c2r(jc(p, r2c(np.eye(2*L))))).real.max()

def d_star(J, g, P, Gam, L=256, hi=3.0):
    """smallest d above which NO plane wave is linearly stable"""
    lo = 0.0
    for _ in range(40):
        mid = 0.5*(lo+hi)
        stable = any((x is not None and x < 1e-10)
                     for x in (pw_growth(L, J, g, P, Gam, mid, 2*np.pi*m/L)
                               for m in range(L)))
        lo, hi = (mid, hi) if stable else (lo, mid)
    return hi


# ------------------------------------------------------------- Lyapunov ----

def lyap(L, J, g, P, Gam, d, dt=0.01, ttrans=400.0, trun=2500.0, seed=5):
    f, jc = make(L, J, g, P, Gam, d)
    rng = np.random.default_rng(seed)
    p = np.sqrt(P/Gam)*(rng.normal(size=L) + 1j*rng.normal(size=L))
    for _ in range(int(ttrans/dt)):
        p = rk4(f, p, dt)
    m = 2*L
    D = r2c(np.linalg.qr(rng.normal(size=(2*L, m)))[0])
    acc = np.zeros(m); n = int(trun/dt); Ssum = 0.0
    for i in range(n):
        p, D = rk4_tan(f, jc, p, D, dt)
        Ssum += np.sum(np.abs(p)**2)
        if (i+1) % 10 == 0:
            Q, R = np.linalg.qr(c2r(D))
            s = np.sign(np.diag(R)); Q, R = Q*s, (R.T*s).T
            acc += np.log(np.abs(np.diag(R))); D = r2c(Q)
    return np.sort(acc/(n*dt))[::-1], Ssum/n, p


# ------------------------------------------------------------------ main ----

if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)
    J, g, P, Gam, z = 1.0, -1.0, 1.0, 1.0, 2

    print("A. (q,p) form reproduces the complex form")
    rng = np.random.default_rng(0); L = 8
    psi = rng.normal(size=L) + 1j*rng.normal(size=L)
    q, pp = to_qp(psi)
    eps = 1e-6
    Dt = -2*J          # make() writes the hopping as a laplacian: Dt_H = Dt_make - Jz
    dHdq = np.array([(hamiltonian(q+eps*np.eye(L)[i], pp, J, g, Dt)
                      - hamiltonian(q-eps*np.eye(L)[i], pp, J, g, Dt))/(2*eps)
                     for i in range(L)])
    dHdp = np.array([(hamiltonian(q, pp+eps*np.eye(L)[i], J, g, Dt)
                      - hamiltonian(q, pp-eps*np.eye(L)[i], J, g, Dt))/(2*eps)
                     for i in range(L)])
    r2 = q**2 + pp**2; gl = 0.5*P - 0.25*Gam*r2
    qdot = dHdp + gl*q
    pdot = -dHdq + gl*pp
    f, _ = make(L, J, g, P, Gam)
    print(f"   max |psi-form - (q,p)-form| = "
          f"{np.max(np.abs(f(psi) - to_psi(qdot, pdot))):.2e}")

    print("\nB. general-Q Bogoliubov formula vs full 2Lx2L Jacobian (L=16, d=0.1)")
    for m in range(0, 8, 2):
        Q = 2*np.pi*m/16
        print(f"   Q={Q:6.3f}  |A|^2={pw_amp2(P,Gam,0.1,Q):+.4f}"
              f"  analytic {pw_growth(16,J,g,P,Gam,0.1,Q):+.8f}"
              f"   Jacobian {pw_growth_jac(16,J,g,P,Gam,0.1,Q):+.8f}")

    print("\nC. the lattice escape: how many plane waves are STABLE?")
    for d in (0.0, 0.05, 0.094, 0.1, 0.25):
        L = 64
        gr = [pw_growth(L, J, g, P, Gam, d, 2*np.pi*m/L) for m in range(L)]
        ex = [x for x in gr if x is not None]
        print(f"   d={d:<6} exist {len(ex):3d}/{L}   stable "
              f"{sum(1 for x in ex if x < 1e-9):3d}")
    print(f"   d* (no stable plane wave above this) = {d_star(J,g,P,Gam):.5f}")

    print("\nD. Lyapunov spectra, L=32")
    for d in (0.0, 0.1, 0.25):
        lam, Sb, pf = lyap(32, J, g, P, Gam, d)
        L = 32
        print(f"\n   d={d}   <S>={Sb:.4f}  (uniform value P L/Gam = {P*L/Gam:.1f})")
        print(f"     lam[:8]  = {np.round(lam[:8], 4)}")
        print(f"     n_pos = {int((lam>1e-4).sum()):2d}   sum lam = {lam.sum():+.5f}")
        print(f"     (P - 2 d z) L - 2 Gam <S> = {(P-2*d*z)*L - 2*Gam*Sb:+.5f}"
              f"   <-- trace rule")
        print(f"     D_KY = {kydim(lam):.4f}  of {2*L}")