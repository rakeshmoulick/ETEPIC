"""
Plots the dispersion graph of from the electric field data
Use %reset -f to clear the python workspace
Data File Invoked: processed_results_all.npz
Run as: 
"""
# TO RUN: python3 dispersion_npz.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import configparser 

file_name = 'processed_results_all.npz'
#path = './'
path = '../data/files/'
# ------------------ Comments -------------------------------------------------
# input parameters specific file
# path to the data folder is ../data/data002_vd_20/files for vd=20
# path to the data folder is ../data/data001_vd_80/files for vd=80
#------------------------------------------------------------------------------
# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
#--------------------------------------------------------------
Tec = tempEcold*e
Teh = tempEhot*e
Ti = tempI*e
# ----------------------------------------------------
Time = 0
#-----------------------------------------------------
ni0 = n0
nec0 = n0/(1+alp)
neh0 = alp*nec0
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DATA_TS = int(NUM_TS/write_interval) + 1

LDH = np.sqrt(eps0*Teh/(neh0*e**2)) # Hot electron Debye Length
LDC = np.sqrt(eps0*Tec/(nec0*e**2)) # Cold electron Debye Length
LD = np.sqrt(eps0*Tec/(n0*e**2)) # Characteristic Debye length

wpi = np.sqrt(ni0*e**2/(eps0*mi))
wp  = np.sqrt(n0*e**2/(eps0*me)) # Total Plasma Frequency
wpec = np.sqrt((nec0*e**2)/(eps0*me))
wpeh = np.sqrt(neh0*e**2/(eps0*me))


print('LDC:',LDC)
print('LDH:',LDH)
print('LD:',LD)

print('wp:%g'%(wp))
print('wpc:%g'%(wpec))
print('wph:%g'%(wpeh))

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
print("The shape of ndeh is: ", ndeh.shape)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
omega = 2*np.pi*np.arange(NUM_TS)/(actual_sim_time) #(DT*NUM_TS) #(actual_sim_time) # Unnormalized Omega
k     = 2*np.pi*np.arange(NC+1)/(NC*dx*LD)          # un-normalized k
kappa = k[1:]*LD
print('Length of k: ',len(k))
print('Max of k: ',np.max(k))

Omega, K = np.meshgrid(omega, k, indexing='ij')
print('Shape of Omega: ',Omega.shape)
#------------------------------------------------------------------------------
halflen = np.array(F.shape, dtype=int)//2
Omega = Omega[:halflen[0],:halflen[1]]
K = K[:halflen[0],:halflen[1]]
F = F[:halflen[0],:halflen[1]]
# Omega /=wp  # Normalized Omega in terms of Total Electron Plasma Frequency
Omega /=wpi # Normalized Omega in terms of Ion Plasma Frequency
K = K*LD  # Normalized K

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++ Raw Analytic Part Calculation ++++++++++++++++++++++++++++++
raw_analytic = True
if raw_analytic:
    wea = np.sqrt((wpec**2)*( (1+3*(k*LDC)**2)/(1 + (1/(k*LDH)**2)) ))
    wla = np.sqrt((wpec**2*(1 + 3*k**2*LDC**2)) + (wpeh**2*(1 + 3*k**2*LDH**2)))
    ea = wea/wpi    # Electron Acoustic Wave
    ep = wla/wpi    # Electron Plasma Wave (Langmuir Wave)
        
    # CF = 1.0/(1.0 + k[1:]**2*LDC**2) # Correction Factor due to electron
    # CS = np.sqrt((Tec*CF + Ti)/mi)
    # wia = CS*k[1:] # ion acoustic mode
    # ia = wia/wp

    
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# For plotting we need to take the absolute value of the electric field. 
# The electric field contains complex data. We take the logarithm of the 
# absolute magnitude of the FFT of the electric field. Because the logarithm 
# of a small positive number is negative, we receive negative FFT amplitude in the plots. 
F = np.abs(F)**2
Z = np.log(np.abs(F))
#print(F)
# ----------------------------------------------------------------------------
# Run the numerical solver using the actual spatial frequencies of your PIC grid
w_num1, w_num2 = solve_IAW_numerical(kappa/LD) 
# kappa is normalized wave number starting from k[1:] instead of k[0:]. 
# Here we need need un-normalized k for w_num1 and w_num2, hence we divided with LD

# Convert w_num1 and w_num2 from absolute rad/s to the PIC plot scale (divided by ion plasma frequency)
ia_num_model1 = w_num1 / wpi
ia_num_model2 = w_num2 / wpi

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([80,80/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 1200                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fig,ax = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)

# plt.pcolor(K, Omega, Z,cmap='rainbow',shading='auto',vmin=np.min(Z),vmax=np.max(Z))
# plt.pcolor(K, Omega, Z,cmap='rainbow',shading='auto',vmin=-3,vmax=2)
# vmin = -3.0
# vmax = None
# extent = [np.min(K), np.max(K), np.min(Omega), np.max(Omega)]
# ax.imshow(Z, extent=extent, interpolation='nearest', aspect='auto', origin='lower', cmap='rainbow', vmin=vmin, vmax=vmax)
# cbar = fig.colorbar(mp.cm.ScalarMappable(cmap='rainbow'), cax=None, ax=ax, orientation='vertical')
# cbar.set_label(r'$\zeta$')

vmin = -3.0
vmax =  1.0
extent = [np.min(K), np.max(K), np.min(Omega), np.max(Omega)]
ax.imshow(Z, extent=extent, interpolation='nearest', aspect='auto', 
          origin='lower', cmap='plasma', vmin=vmin, vmax=vmax)
# Update the colorbar mapping to match the new vibrant colormap
cbar = fig.colorbar(mp.cm.ScalarMappable(norm=mp.colors.Normalize(vmin=vmin, vmax=vmax), cmap='turbo'), 
                    cax=None, ax=ax, orientation='vertical')
cbar.set_label(r'$\zeta$')


if raw_analytic:
    # plt.plot(k*LD, ep, color='k', linestyle='-.', lw = 1.0, alpha = 0.8, label='$EPW$')
    # plt.plot(k*LD, ea, color='k', linestyle='--', lw = 1.0, alpha = 0.8, label='$EAW$')
    plt.plot(kappa, ia_num_model1, color='lime', linestyle='--',  lw = 1.0, alpha = 0.8, label='$IAW (M1)$')
    plt.plot(kappa, ia_num_model2, color='white', linestyle='-.', lw = 1.0, alpha = 0.8, label='$IAW (M2)$')    


ax.set_xlabel(r'$k \lambda_{D}$')
ax.set_ylabel(r'$\omega/\omega_{pi}$')
#ax.set_xlim([0, np.max(K)])
ax.set_xlim([0, 3.0])
ax.set_ylim([0, 0.05* (wp / wpi)])
leg = ax.legend(loc='upper right',framealpha=0.5)
plt.savefig(pjoin('./','dispersion.jpg'),dpi=dpi)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plt.show()
