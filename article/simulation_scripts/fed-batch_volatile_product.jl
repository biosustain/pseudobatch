# This file setups, and runs the ODE systems for the standard exponential fedbatch 
# process with product inhibition and a volatile product simulation. Finally, some 
# simple data processing is carried out to format the data into a self containing 
# .csv file.

using OrdinaryDiffEq, DiffEqCallbacks, CSV, DataFrames
include("fermentation_simulator_functions.jl")
include("fermentation_utils_functions.jl")
include("standard_parameters.jl")

# Override saveat to save more data points
save_ode_timesteps = LinRange(0,60, 2000)

# Product inhibition parameter
MW_o2 = 32 # g/mol
Ki_p = 2
evaporation_rate_product = 0.1 # 1/h
O2_flow_rate_in = 10 # g/h 
Yxo2 = 0.01 # gO2 / gDW

## follow numbers are obtain from example 10.2 in Bioprocess Engineering Principles,
## Villadsen et al.
kla = 748 # 1/h 
c_o2_sat = (0.265*10e-3) * MW_o2 # g/L
c_o2_init = c_o2_sat # g/L

# # calculate realistic Yxo2
# MW_co2 = 44.01 # g/mol
# Yco_moles = 1 # RQ in moles
# Yco_g = Yco_moles * (MW_co2/MW_o2) # RQ in g/g
# Yxo2 = Yxco2/Yco_g # (gCO2 / dDW ) * (gCO2 / gO2)^-1 = gO2 / gDW

#
state_variable_names = [:m_Glucose, :m_Biomass, :m_Product, :m_CO2, :v_Volume, :v_Feed_accum, :m_CO2_gas, :m_product_gas, :m_O2_in, :m_O2, :m_O2_gas]
init_cond = [s0*V0, x0*V0, p0*V0, co2_0*V0, V0, 0., 0., 0., 0., c_o2_init*V0, 0.] # input for the model is masses not concentrations, accum feed at time zero is 0
sample_volume_dict = Dict(zip(sampling_times, repeat([sample_volume], n_samples)))
ode_input_p = [Kc_s, mu_max, Yxs, Yxp, Yxco2, F0, mu0, s_f, Ki_p, evaporation_rate_product, Yxo2, O2_flow_rate_in, kla, c_o2_sat]

ode_func = ODEFunction(fedbatch_prod_inhib_volatile!, syms=state_variable_names)
prob = ODEProblem(ode_func, init_cond, tspan, ode_input_p)
cb = PresetTimeCallback(sampling_times, affect_sample_volatile_compounds!, filter_tstops = true) # creates new sampling callback
sol = solve(prob, Tsit5(), callback=cb, saveat=save_ode_timesteps, abstol = 1e-12) # solve ode system with new sampling callback

df = DataFrame(sol)
for state_variable in state_variable_names
    # It is not meaningfull to calculate the concentration of the feed_accum or volume
    if state_variable in [:v_Feed_accum, :v_Volume, :m_CO2_gas, :m_product_gas, :m_O2_in, :m_O2_gas]
        continue
    end
    metabolite = split(string(state_variable), "_")[2]
    new_colname = string("c_", metabolite)
    df[!, new_colname] = df[!, state_variable] ./ df[!, :v_Volume] 
end
transform!(df, [:c_Glucose, :c_Product] => ByRow((c_s, c_p) -> monod_with_product_inhibition(c_s, c_p, mu_max, Kc_s, Ki_p)) => :mu_true)
transform!(df, [:c_O2] => ByRow((c_O2) -> volumetric_gas_to_liquid_transfer_rate(c_o2_sat, c_O2, kla)) => :qt_O2)

# adding sampling information
insertcols!(df, 1, "sample_volume" => NaN)
df[[x in sampling_times for x in df.timestamp], "sample_volume"] .= sample_volume

# adding the parameter values to the dataframe
output_header = [:Kc_s, :mu_max, :Yxs, :Yxp, :Yxco2, :F0, :mu0, :s_f, :Ki_p, :evaporation_rate_product, :Yxo2, :O2_flow_rate_in, :kla, :c_o2_sat] 
insertcols!(df, 1, (output_header .=> ode_input_p)...)

CSV.write(string("data/volatile_product.csv"), df)

