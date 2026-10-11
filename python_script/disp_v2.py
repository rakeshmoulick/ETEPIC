# IAW Dispersion solver considering ion inertia

import numpy as np
import matplotlib.pyplot as plt
# from scipy.constants import e, me, eps0
from scipy.constants import value as constants
import sys
import configparser 
#--------------------------------------------
# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')
#--------------------------------------------
config_path = sys.argv[1]
# Read the configuration file
config = configparser.ConfigParser()
config.read(config_path)

# ---------------IMPORT REQUIRED DATA FROM INPUT.INP-----------
NUM_TS = config.getint('time','NUM_TS')
DT_coeff = config.getfloat('time','DT_coeff')
NC = config.getint('grid','NC')
write_interval = config.getint('diagnostics','write_interval')
n0 = config.getfloat('population','density')
alp = config.getfloat('population','alp')
mi = config.getfloat('population','massI')
mi = mi*AMU
nParticlesE = config.getint('population','nParticlesE')
nParticlesI = config.getint('population','nParticlesI')
tempEcold = config.getfloat('population','tempEcold') 
tempEhot = config.getfloat('population','tempEhot')
tempI = config.getfloat('population','tempI')

# =============================================================================
# 1. PARAMETERS
# =============================================================================
ni0 = n0
nec0 = n0/(1+alp)
neh0 = alp*nec0
Tec = tempEcold*e
Teh = tempEhot*e
Ti = tempI*e

# Derived plasma frequencies and thermal speeds
omega_pi = np.sqrt((ni0 * e**2) / (eps0 * mi))
omega_pe = np.sqrt((ni0 * e**2) / (eps0 * me))
v_ti = np.sqrt(Ti / mi)
v_tc = np.sqrt(Tec / me)
v_th = np.sqrt(Teh / me)

# Effective Debye length for normalization
L_eff = np.sqrt(eps0 / (e**2 * (nec0/Tec + neh0/Teh)))

# Collisional frequency (nu)
nu = 0.01 * omega_pi

# Define k range (normalized by L_eff)
k_norm = np.logspace(-2, 2, 500)
k_values = k_norm / L_eff

# =============================================================================
# 2. DISPERSION SOLVER
# =============================================================================
def solve_IAW_physics(k_list):
    omega_sat = []   # Model A: Simplified (Saturates)
    omega_drop = []  # Model B: Inertia corrected (Drops)
    gamma = []


    for k in k_list:
        # --- MODEL A: THE SIMPLIFIED QUADRATIC (Saturates) ---
        # (1 + k^2L^2)w^2 + i*nu*(1 + k^2L^2)w - wp_i^2*k^2L^2 = 0
        y = (k * L_eff)**2
        roots_a = np.roots([(1 + y), 1j * nu * (1 + y), -(omega_pi**2 * y)])
        omega_sat.append(np.max(np.real(roots_a)))

        # --- MODEL B: TInertia corrected (Drops)) ---
        # Based on: 1 - chi_i + chi_ec + chi_eh = 0
        # Incorporates ion thermal speed (v_ti) which causes the breakdown.
        chi_e = ( (nec0 * e**2 / (eps0 * me)) / (k**2 * v_tc**2) +
                  (neh0 * e**2 / (eps0 * me)) / (k**2 * v_th**2) )

        # The effective restoring force term (RHS)
        # Note: - k^2 v_ti^2 is the thermal pressure term that leads to the drop
        RHS = (omega_pi**2 / (1 + chi_e)) - (k**2 * v_ti**2)

        # Solving the quadratic: w^2 + i*nu*w - RHS = 0
        roots_b = np.roots([1.0, 1j * nu, -RHS])

        re_vals = np.real(roots_b[np.argmax(np.real(roots_b))])
        im_vals = np.imag(roots_b)

        # If the real part is negligible, the wave is effectively overdamped (frequency = 0)
        omega_drop.append(re_vals if re_vals > 1e-3 else 0.0)
        gamma.append(im_vals)

    return np.array(omega_sat), np.array(omega_drop), np.array(gamma)

# Execute Solver
w_sat, w_drop, w_gamma = solve_IAW_physics(k_values)

# =============================================================================
# 3. PLOTTING
# =============================================================================
plt.figure(1,figsize=(10, 6))

plt.plot(k_norm, w_sat / omega_pi, 'b-',  linewidth=2.5, label='Model 1')
plt.plot(k_norm, w_drop / omega_pi, 'r--', linewidth=2.5, label='Model 2')

plt.xscale('log')
plt.axvline(1.0, color='gray', linestyle=':', label=r'Debye Scale ($k\lambda_D=1$)')
plt.tick_params(axis='both', which='major', labelsize=24)
plt.xlabel(r'$k\lambda_D$', fontsize=24)
plt.ylabel(r'$\omega_r / \omega_{pi}$', fontsize=24)
# plt.title('IAW Dispersion: Comparing Saturation vs. Thermal Breakdown', fontsize=14)
plt.legend(fontsize=20)
plt.grid(True, which='both', alpha=0.3)
plt.tight_layout()
# ------------------------------------------------------------------------------------
plt.figure(2,figsize=(10, 6))
plt.plot(k_norm, w_gamma[:,0] / omega_pi, 'r-', linewidth=2.5, label='Model 2')
plt.plot(k_norm, w_gamma[:,1] / omega_pi, 'r-', linewidth=2.5, label='_Hidden')
# plt.ylim(-5,5)
plt.xscale('log')
plt.axvline(1.0, color='gray', linestyle=':', label=r'Debye Scale ($k\lambda_D=1$)')
plt.tick_params(axis='both', which='major', labelsize=24)
plt.xlabel(r'$k\lambda_D$', fontsize=24)
plt.ylabel(r'$\gamma / \omega_{pi}$', fontsize=24)

plt.legend(fontsize=20)
plt.grid(True, which='both', alpha=0.3)
plt.tight_layout()
plt.show()
