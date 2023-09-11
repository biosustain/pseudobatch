# This file setups, and runs the ODE systems for the standard exponential fedbatch 
# process simulation. Finally, some simple data processing is carried out to format 
# the data into a self containing .csv file.
using OrdinaryDiffEq, DiffEqCallbacks, CSV, DataFrames
include("fermentation_simulator_functions.jl")
include("fermentation_utils_functions.jl")
include("standard_parameters.jl")


# Setting up ODE input
state_variable_names = [:m_Glucose, :m_Biomass, :m_Product, :m_CO2, :v_Volume, :v_Feed_accum, :m_CO2_gas]
init_cond = [s0*V0, x0*V0, p0*V0, co2_0*V0, V0, 0., 0.] # input for the model is masses not concentrations, accum feed at time zero is 0
sample_volume_dict = Dict(zip(sampling_times, repeat([sample_volume], n_samples))) # defines timepoints and sample volumes
ode_input_p = [Kc_s, mu_max, Yxs, Yxp, Yxco2, F0, mu0, s_f, evaporation_k]

ode_func = ODEFunction(fedbatch!, syms=state_variable_names)
prob = ODEProblem(ode_func, init_cond, tspan, ode_input_p)
cb = PresetTimeCallback(sampling_times, affect_sample!, filter_tstops = true) # creates new sampling callback
sol = solve(prob, Tsit5(), callback=cb, saveat=save_ode_timesteps, abstol = 1e-12) # solve ode system with new sampling callback

# Processing simulated data 
df = DataFrame(sol)
for state_variable in state_variable_names
    # It is not meaningfull to calculate the concentration of the feed_accum or volume
    if state_variable in [:v_Feed_accum, :v_Volume, :m_CO2_gas]
        continue
    end
    metabolite = split(string(state_variable), "_")[2]
    new_colname = string("c_", metabolite)
    df[!, new_colname] = df[!, state_variable] ./ df[!, :v_Volume] 
end

# adding the true growth rate
transform!(df, [:c_Glucose] => ByRow((c_s) -> monod_kinetics(c_s, mu_max, Kc_s)) => :mu_true)

# adding sampling information
insertcols!(df, 1, "sample_volume" => NaN)
df[[x in sampling_times for x in df.timestamp], "sample_volume"] .= sample_volume

# adding the parameter values to the dataframe
output_header = [:Kc_s, :mu_max, :Yxs, :Yxp, :Yxco2, :F0, :mu0, :s_f, :evaporation_k] 
insertcols!(df, 1, (output_header .=> ode_input_p)...)

CSV.write(string("data/standard_fed-batch_process.csv"), df)
