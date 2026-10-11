# ============================================================
# IAW dispersion relation used in the PIC dispersion plot
# ============================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
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
file_name = 'processed_results_all.npz'
#path = './'
path = '../data/files/'
# -------------------------------------------
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

wpi = np.sqrt(ni0*e**2/(eps0*mi))
wp  = np.sqrt(n0*e**2/(eps0*me)) # Total Plasma Frequency
wpec = np.sqrt((nec0*e**2)/(eps0*me))
wpeh = np.sqrt(neh0*e**2/(eps0*me))

# Derived plasma frequencies and thermal speeds
omega_pi = np.sqrt((ni0 * e**2) / (eps0 * mi))
omega_pe = np.sqrt((ni0 * e**2) / (eps0 * me))
v_ti = np.sqrt(Ti / mi)
v_tc = np.sqrt(Tec / me)
v_th = np.sqrt(Teh / me)

# Use the same Debye length as dispersion_npz_EF.py
LD = np.sqrt(eps0 * Tec / (n0 * e**2))
LDC = np.sqrt(eps0 * Tec / (nec0 * e**2))
DATA_TS = int(NUM_TS/write_interval) + 1
DT = DT_coeff*(1.0/wp)

# =============================================================================
# NUMERICAL DISPERSION SOLVER FROM disp_v2.py
# =============================================================================
# Define the collision frequency based on the PIC ion plasma frequency
nu = 0.01 * wpi 

# Calculate the effective multi-species Debye length used by the solver
L_eff = np.sqrt(eps0 / (e**2 * (nec0/Teh + neh0/Teh))) # Note: Ensure nec0/Tec + neh0/Teh matches your physical model

def solve_IAW_numerical(k_list):
    omega_model1 = []   # Model 1: Simplified (Saturates)
    omega_model2 = []   # Model 2: Inertia corrected (Drops)

    # Derived thermal speeds needed for Model 2
    v_ti = np.sqrt(Ti / mi)
    v_tc = np.sqrt(Tec / me)
    v_th = np.sqrt(Teh / me)

    for k in k_list:
        if k == 0:
            omega_model1.append(0.0)
            omega_model2.append(0.0)
            continue

        # --- MODEL 1: Numerical Quadratic (Saturates) ---
        y = (k * L_eff)**2
        roots_a = np.roots([(1 + y), 1j * nu * (1 + y), -(wpi**2 * y)])
        omega_model1.append(np.max(np.real(roots_a)))

        # --- MODEL 2: Numerical Inertia Corrected (Drops) ---
        chi_e = ( (nec0 * e**2 / (eps0 * me)) / (k**2 * v_tc**2) +
                  (neh0 * e**2 / (eps0 * me)) / (k**2 * v_th**2) )

        RHS = (wpi**2 / (1 + chi_e)) - (k**2 * v_ti**2)
        roots_b = np.roots([1.0, 1j * nu, -RHS])

        re_vals = np.real(roots_b[np.argmax(np.real(roots_b))])
        omega_model2.append(re_vals if re_vals > 1e-3 else 0.0)

    return np.array(omega_model1), np.array(omega_model2)


#++++++++++++++++++++ UNPACK ++++++++++++++++++++++++++++++++++++++++++++++++++
if os.path.exists(pjoin(path,file_name)):
    data = np.load(pjoin(path,file_name))
    x = data['x']
    ndec = data['ndec']
    ndeh = data['ndeh']
    rho = data['rho']
    phi = data['phi']
    EF = data['EF']
else:
    print('No data')
    exit()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
x = x.reshape(DATA_TS,(NC+1))
x = x[0,:] # SPATIAL GRID: Select only a column, all the columns have the same value
dx = x[1]-x[0] # Normalized StepSize

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
"""
In each row of EF, there are NC+1 = 1025 columns. Data is put in the first row first
and then it goes into the second row. Therefore, all the values of the EF for the
first time step is written in first row in 1025 columns. The data of next time step
is in the second row and so on. Similar is the case for writing and reshaping x.
w_pe*t = NUM_TS*DT
"""
EF = EF.reshape(DATA_TS,(NC+1))
print("The shape of EF is: ", EF.shape)

wpet_1 = 0
wpet_2 = NUM_TS*DT_coeff
y1 = wpet_1/(DT_coeff*write_interval)
y2 = wpet_2/(DT_coeff*write_interval)
E = EF[int(y1):int(y2),:]
print("The shape of E (reduced EF) is: ", E.shape)
#+++++++++++++++++++++++++ FFT of E-field data ++++++++++++++++++++++++++++++++
F = np.fft.fftn(E, norm='ortho')
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define Omega and K
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NUM_TS1 = wpet_1/DT_coeff
NUM_TS2 = wpet_2/DT_coeff
actual_sim_time = (NUM_TS2 - NUM_TS1)*(DT) # Previously it was (NUM_TS*DT)
# -----------------------------------------------------------------------------
omega   = 2*np.pi*np.arange(NUM_TS)/(actual_sim_time) #(DT*NUM_TS) #(actual_sim_time) # Unnormalized Omega
k  = 2*np.pi*np.arange(NC+1)/(NC*dx*LD)          # un-normalized k

# --------------------------------------------------------------------
# Run the numerical solver using the actual spatial frequencies of your PIC grid
w_num1, w_num2 = solve_IAW_numerical(k)

# Convert from absolute rad/s to the PIC plot scale (divided by electron plasma frequency)
num_model1_normalized = w_num1 / wp
num_model2_normalized = w_num2 / wp
# --------------------------------------------------------------------


# +++++++++++++++++++++ Plotting Numerical Models ++++++++++++++++++++++++++++++
# Plotting the numerical solutions on top of the linear PIC background map
X_axis_normalized = k * LD

# Model 1 (Blue dashed line to match disp_v2 styling)
plt.plot(X_axis_normalized, num_model1_normalized, color='black', 
        linestyle='--', lw=1.5, label='Numerical Model 1')

# Model 2 (Red solid line to match disp_v2 styling)
# plt.plot(X_axis_normalized, num_model2_normalized, color='red', 
#         linestyle='-', lw=1.5, label='Numerical Model 2')
plt.legend(loc=0)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

plt.show()