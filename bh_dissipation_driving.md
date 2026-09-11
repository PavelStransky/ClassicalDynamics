# Coherently driven dissipative Bose–Hubbard model

Companion note to `dissipative_bh_summary.pdf`. Same conventions throughout:
ψ_i = b_i/√N, ħ_eff = 1/N, g = UN fixed, κ = γ_l − γ_p, γ_d = dephasing rate.
Reference implementation: `driven_bh.py` (tested); figure `driven_bh.png`.

## 1. Model

Lab frame: on-site frequency ω_0, drive at ω_p, amplitude F.
Rotating frame at ω_p, with Δ = ω_p − ω_0:

    H = −J Σ_<ij> (b_i† b_j + h.c.) − Δ Σ_i n_i + (U/2) Σ_i n_i(n_i−1)
        + Σ_i (F b_i† + F* b_i)

    dρ/dt = −i[H,ρ] + Σ_i ( γ_p D[b_i†] + γ_l D[b_i] + γ_d D[n_i] ) ρ

D[A]ρ = AρA† − ½{A†A, ρ}. The dissipators are unchanged by the rotating-frame
transformation, since D[b] is invariant under b → e^{−iω_p t} b.

## 2. Scaling of the drive — F ∝ √N

F b_i† = F √N ψ_i*, and H_cl = H/N, so the drive term in H_cl is (F/√N)ψ_i* + c.c.
Therefore the correct fixed quantity is

    f ≡ F / √N          (NOT F/N^{3/2})

i.e. drive amplitude ∝ √N, drive *power* ∝ N — extensive, matching the O(N) loss.
Same statement as in the Dicke model, where (λ/√N)(aJ_+ + h.c.) ~ λN.

## 3. Classical equation of motion

With z the coordination number, Δ̃ ≡ Δ + Jz the detuning measured from the bottom
of the band (z = 2, Δ̃ = Δ + 2J for a periodic chain), and
lap_j = ψ_{j+1} + ψ_{j−1} − 2ψ_j:

    dψ_j = [ iJ lap_j + iΔ̃ ψ_j − ig|ψ_j|²ψ_j − i f − (κ/2)ψ_j ] dt
           − i √γ_d ψ_j ∘ dW_j

This is the **discrete Lugiato–Lefever equation** (driven-dissipative Kerr lattice).
f real positive WLOG. The drive **breaks U(1)** — the structural change that makes
everything else possible.

### Norm equation no longer closes

    d/dt Σ_i |ψ_i|² = −κ Σ_i |ψ_i|² − 2f Σ_i Im ψ_i

Cauchy–Schwarz gives dS/dt ≤ −κS + 2f√(LS), i.e.

    **absorbing ball:  Σ_i|ψ_i|² ≤ 4f²L/κ²**

a compact absorbing set of *fixed* radius, versus the shrinking ball e^{−κt}S(0) of
the undriven model. This is the single structural fact that permits a nontrivial
attractor.

## 4. Uniform fixed points — optical bistability

ψ_j = ψ (lap = 0), n = |ψ|²:

    ψ = f / [ (Δ̃ − g n) + i κ/2 ]
    **n [ (Δ̃ − g n)² + κ²/4 ] = f²**
    ⇔ g²n³ − 2gΔ̃ n² + (Δ̃² + κ²/4) n − f² = 0

Turning points from d(f²)/dn = 3g²n² − 4gΔ̃n + Δ̃² + κ²/4 = 0:

    n_± = [ 2Δ̃ ± √(Δ̃² − ¾κ²) ] / (3g)

    **Bistability ⇔ gΔ̃ > 0 and |Δ̃| > (√3/2) κ**

(the textbook Kerr condition Δ > √3 γ with γ = κ/2 the amplitude decay rate).

## 5. Bogoliubov / modulational stability

ψ_j = ψ_0 + δψ_j, Fourier with ε̄_k = 2J(1 − cos k), A_k ≡ ε̄_k − Δ̃ + 2gn:

    **λ_k^± = −κ/2 ∓ i √(A_k² − g²n²)**

- A_k² > g²n²: complex pair, Re λ = −κ/2 → damped Bogoliubov mode.
- A_k² < g²n²: **real** pair, Re λ = −κ/2 ± √(g²n² − A_k²).

Instability ⇔ **A_k² < g²n² − κ²/4**, which needs |g|n > κ/2.

### Two consequences

**(a) k = 0 instability ⇔ negative slope.** At k=0, A_0² − g²n² + κ²/4 =
3g²n² − 4gΔ̃n + Δ̃² + κ²/4 = d(f²)/dn *identically*. The k=0-unstable set is exactly
the middle (negative-slope) branch — a saddle, as required.

**(b) No limit cycles in the uniform sector.** Whenever λ_0 is complex,
Re λ_0 = −κ/2 exactly, so a Hopf bifurcation cannot occur at k = 0. Stronger: the
single-site (L=1) flow has divergence −κ everywhere, so by **Bendixson–Dulac** it
has *no* periodic orbits at all. **The single driven Kerr resonator is purely
bistable — chaos requires the lattice.**

### Where MI lives — a band-sliding picture (CORRECTED)

Only the combination ε̄_k − Δ̃ enters A_k. So picture a **fixed** instability window
|A| < σ ≡ √(g²n² − κ²/4), and a band occupying A ∈ [A_0, A_0 + 4J] with
A_0 = 2gn − Δ̃. **Δ̃ translates the band rigidly across that window.**

- k = 0 sits at A_0; it is in the window ⟺ |2gn − Δ̃| < σ ⟺ the saddle branch.
- k ≠ 0 in the window while k = 0 is out ⟺ MI. For J > 0 (A_k increasing in k):

      2gn + σ  <  Δ̃  <  2gn + σ + 4J

  and for J < 0:   2gn − σ − 4|J|  <  Δ̃  <  2gn − σ.   Plus |g|n > κ/2 for σ real.
  On a small lattice an *allowed* ε̄_k = 2J(1 − cos 2πm/L) must land in the window.

**Correction to an earlier claim in this note:** MI does NOT require Jg < 0. That
argument only covered the upper branch (where gn ≳ Δ̃). With J = g = +1, κ = 1 there
are many MI points on the *lower* branch at large Δ̃ (e.g. Δ̃ = 3, f = 1.95,
n = 0.669, k* = π/2; strongest found Δ̃ = 17.5, f = 28.2, n = 5.73, Re λ = +4.86).

What *does* hold empirically across every scan run here: with **Jg > 0 the MI
saturates into a stable stationary pattern** (λ_max = −κ/2, D_KY = 0) even at growth
rate 4.9, whereas **chaos was found only for Jg < 0**. So the sign governs the
*nonlinear outcome*, not linear stability. (Optics: Jg < 0 is the
anomalous-dispersion/LLE condition.) The exact symmetry ψ → ψ*,
(Δ̃, g, J, f) → −(Δ̃, g, J, f) exchanges the two Jg < 0 realisations.

## 5b. The meaning of Δ, and the Δ̃ = 0 limit

**Δ = ω_p − ω_0** is the drive detuning from the bare site frequency; it enters as
−Δ Σ n_i, i.e. formally a chemical potential μ = Δ. In the *undriven* model it is
**pure gauge** — removable by ψ → e^{iΔt}ψ — which is exactly why the summary could
set μ = 0 without loss. The drive breaks U(1) and thereby **promotes Δ from gauge to
physics**: the rotation would turn −if into −i f e^{−iΔt}.

Δ̃ = Δ + Jz is the detuning from the k = 0 Bloch mode (band ε_k = −2J cos k, bottom
ε_0 = −Jz), the only mode a uniform drive couples to. Δ̃ − gn is the Kerr-shifted
detuning; the response peaks at Δ̃ = gn.

**Δ̃ = 0** (drive resonant with the k = 0 mode, i.e. Δ = −Jz):
- state equation g²n³ + (κ²/4)n = f², with d(f²)/dn = 3g²n² + κ²/4 > 0 always
  ⇒ **single-valued, monostable, no folds, no hysteresis** (consistent with the
  threshold |Δ̃| > (√3/2)κ — Δ̃ = 0 is on the monostable side of the cusp);
- k = 0 can never be unstable (no middle branch exists);
- small f: n → 4f²/κ², the **peak of the Lorentzian** — maximal linear response;
- large f: n → (f/|g|)^{2/3}. The Kerr shift detunes the mode from its own drive as
  it fills, so the response saturates to a power law instead of growing as f².
- MI at k ≠ 0 survives. Numerically (J=1, g=−1, κ=1, L=8): stable to f ≈ 0.9, MI
  beyond; λ_max ≈ +0.07…+0.17 near f ≈ 1.5–2, D_KY ≈ 4.6–5.1 of 16 — chaos, but far
  weaker than at Δ̃ = −4 (λ_max ≈ 0.74, D_KY ≈ 9.4). With g = +1: stable at all f.

**Δ = 0** (drive on the bare-site resonance) is *not* the resonant case: Δ̃ = Jz = 2J.
At J = κ = 1 (Δ̃ = 2): repulsive g = +1 gives gΔ̃ > 0 and 2 > 0.866 ⇒ **bistable**,
folds n = (0.732, 1.934), f ∈ (0.701, 1.166) — hysteresis, but only fixed-point
attractors. Attractive g = −1 gives gΔ̃ < 0 ⇒ monostable, MI from f ≈ 2, settling to
a stationary pattern. So at Δ = 0 the lattice itself supplies the detuning Jz, and
the sign of g decides whether the system is bistable.

As f → 0 at fixed Δ, n → f²/(Δ̃² + κ²/4) → 0, the attractor collapses back to ψ = 0,
λ → −κ/2, and Δ becomes gauge again: the driven model connects continuously to the
undriven one, with f the parameter that switches Δ on.

## 6. What survives from the undriven analysis

- The drive is a Hamiltonian term ⇒ the flow is **still conformally symplectic**,
  M^T Ω M = e^{−κt} Ω. Hence **λ_k + λ_{2L+1−k} = −κ** and **Σλ_k = −κL** exactly.
  These replace the norm law (16) as the numerical validation invariant.
- Pairing ⇒ at most L positive exponents.
- U(1) is broken, so the undriven model's extra pair of exponents at −κ/2
  disappears; the only vanishing exponent is the flow direction.
- Dephasing γ_d acts exactly as before (exact random rotations in the splitting);
  it smears the attractor into a random attractor with a smooth density.

## 7. Numerical results (J = 1, g = −1, κ = 1, Δ̃ = −4, periodic)

Folds at n = (1.365, 3.968), f ∈ (0.998, 3.134) bistable. Most unstable mode k = π.

Route as f increases: stable uniform (low branch) → non-uniform fixed point →
limit cycle (f ≈ 2.5–2.9) → **chaos (f ≳ 3.0)**, with strong multistability.

| L | f | λ_max | Σλ | pairing dev | D_KY | 2L |
|---|---|-------|-----|-------------|------|-----|
| 8 | 3.2 | +0.74 | −7.999996 | 3.7e−3 | 9.41 | 16 |
| 16 | 3.2 | +0.78 | −15.999992 | 3.9e−3 | 18.03 | 32 |

**Fractional D_KY ⇒ a genuine fractal strange attractor** — precisely what the
undriven model cannot have, where the only attractor is ψ = 0 and λ → −κ/2.

Checks passed: fold formula vs numerical turning points; analytic λ_k reproduces the
full 2L×2L Jacobian spectrum to 1e−13; k=0 instability ⇔ d(f²)/dn < 0 at every
sampled point; Σλ = −κL to 1e−5; pairing to 4e−3; Δt-halving consistent.