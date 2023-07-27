# This file contains the parameters that are used in the simulations.
# In some situations the parameters are overwritten in the simulation scripts.
# This script is run in the beginning of all simulation scripts, that 
# load the parameters from this file.

# Load the functions that are used to generate the parameters
include("fermentation_utils_functions.jl")

### Fixed parameters
## Molar masses
MW_p = 150 # g/mol just a generic product with 8 carbon atoms
MW_s = 180.156 # g/mol
MW_co2 = 44.01 #g/mol
MW_x = 23.58 # gDW/cmol (for s. cerevisiae from equation 3.19 in Bioreaction Engineering Principples)

# physical
tspan = (0., 60.)
x0 = 0.5 # g/L
p0 = 0. # g/L
co2_0 = 0. # g/L
V0 = 1000. # µL 
save_ode_timesteps = LinRange(0,60, 1000)

# biological parameters
Yxs = 1.85 # gGLC/gDW from Table 7.3 in Bioreaction Engineering principples
product_yield_frac = 0.4
mu0 = 0.1 # 1/h
mu_max = 0.3 # 1/h
Kc_s = 150 / 1000 # mg / ml from table 7.1 from Bioreaction Engineering Principples 
s_f = 100 # gGLC/L = µgGLC/µL
n_samples = 8
sample_volume = 170 # µL

## Generated parameters
# To avoid an initial adaption phase, the initial substrate concentration is
# calculated to be at steady state at the target growth rate
s0 = monod_ss_substrate_conc(mu0, mu_max, Kc_s) # g/L, to ensure 

# The yield coefficients for CO2 and product are calculatedto ensure stoichiometric
# consistency.
Yxs_cmoles = gg_to_cmolcmol(Yxs, MW_x, MW_s, 1, 6)
Yxp_cmoles = Yxs_cmoles * product_yield_frac # Yxp is a fraction of the Yxs
Yxco2_cmoles = Yxs_cmoles - Yxp_cmoles # Yxco2 is the remianing carbon
Yxp = cmolcmol_to_gg(Yxp_cmoles, MW_x, MW_p, 1, 8)
Yxco2 = cmolcmol_to_gg(Yxco2_cmoles, MW_x, MW_co2, 1, 1)

# Calculating feeding profile and sampling times
F0 = calc_F0(Yxs, x0, V0, mu0, s_f, s0)
sampling_times = LinRange(tspan[1]+10, tspan[2], Int(round(n_samples)))

