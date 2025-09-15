import math, random
import numpy as np
#import WINNER_ALL_v3 as W

HBAR = 1.054_571_817e-34
H = 2.0 * math.pi * HBAR
Q_E = 1.602_176_634e-19  # electron charge

# CNF klauzule jako [(var_idx, is_positive), (..), (..)]
# u_c(sigma) = 1 (false), jinak 0
def clause_unsat(clause, sigma):
    sat = False
    for (i, is_pos) in clause:
        lit = (sigma[i] == +1) if is_pos else (sigma[i] == -1)
        sat = sat or lit
        if sat: break
    return 0 if sat else 1

def unsat_count(cnf, sigma):
    return sum(clause_unsat(c, sigma) for c in cnf)

# Fázová chyba pro dané sigma (bez auto-locku, čistá geometrie toku)
def delta_phi(cnf, sigma, Phi_base, Phi_unit, x, t, omega, q=Q_E):
    u = unsat_count(cnf, sigma)
    #print("UNSAT clauses u =", u)
    Phi0 = (2*math.pi*HBAR)/q
    Phi  = Phi_base + u*Phi_unit
    dphi_geo = 2*math.pi*((Phi / Phi0) % 1.0)
    total = math.pi*x + omega*t + dphi_geo
    n_star = int(round(total/(2.0*math.pi)))
    return ((total - 2.0*math.pi*n_star + math.pi) % (2.0*math.pi)) - math.pi

# Jednoduchý GSAT: hledá sigma s |δφ| <= tol
def solve_sat_by_resonance(cnf, n, tol=0.02, steps=2000):
    # nastavení rezonance: chci aby ideál (u=0) dával δφ=0 ⇒ nastav Phi_base podle x,t,omega
    x, t, omega = 0.25, 1.0, 2*math.pi
    Phi_unit = (2*math.pi*HBAR)/Q_E / 4.0  # Φ0/4 ⇒ π/2 na špatnou klauzuli
    # vyrob Phi_base, aby pro u=0 vycházelo δφ=0
    # řeš: πx + ωt + 2π frac(Phi_base/Φ0) = 2πk ⇒ zvol k = round((πx+ωt)/2π)
    Phi0 = (2*math.pi*HBAR)/Q_E
    k0 = int(round((math.pi*x + omega*t)/(2.0*math.pi)))
    frac_base = (k0 - (math.pi*x + omega*t)/(2.0*math.pi)) % 1.0
    Phi_base = frac_base * Phi0

    sigma = [random.choice([-1,+1]) for _ in range(n)]
    best = delta_phi(cnf, sigma, Phi_base, Phi_unit, x, t, omega)
    for step in range(steps):
        if abs(best) <= tol:
            return True, step, sigma, best
        # vyzkoušej lokálně nejlepší flip
        j_star, val_star = None, best

        #for j in range(n):
        #    sigma[j] *= -1
        #    val = delta_phi(cnf, sigma, Phi_base, Phi_unit, x, t, omega)
        #    if abs(val) < abs(val_star):
        #        j_star, val_star = j, val
        #    sigma[j] *= -1

        val = np.array([delta_phi(cnf, sigma[:j] + [-sigma[j]] + sigma[j + 1:], Phi_base, Phi_unit, x, t, omega)
                                       for j in range(n)])

        if j_star is None:
            # náhodný krok (WalkSAT prvek)
            j_star = random.randrange(n)
            sigma[j_star] *= -1
            best = delta_phi(cnf, sigma, Phi_base, Phi_unit, x, t, omega)
        else:
            sigma[j_star] *= -1
            best = val_star
    return False, steps, sigma, best

def random_3sat(n: int, m: int):
    cnf = []
    for _ in range(m):
        lits = []
        vars_ = random.sample(range(n), 3) if n >= 3 else [0,1,2][:n]
        for v in vars_:
            lits.append((v, random.choice([True, False])))
        while len(lits) < 3:
            lits.append((random.randrange(max(1,n)), True))
        cnf.append(tuple(lits[:3]))
    return cnf

if __name__=="__main__":
    formula = random_3sat(n=10000, m=42000)
    print(formula)
    out = solve_sat_by_resonance(cnf=formula, n=10000)
    print(out)