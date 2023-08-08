# This file setups, and runs the ODE systems for the CHO cell like fedbatch 
# process simulation. Finally, some simple data processing is carried out to format 
# the data into a self containing .csv file.

using OrdinaryDiffEq, DiffEqCallbacks, CSV, DataFrames
include("fermentation_simulator_functions.jl")
include("fermentation_utils_functions.jl")
include("standard_parameters.jl")

# Additianal parameters
Kc_s2 = 1 # mg / ml from table 7.1 from Bioreaction Engineering Principples
mu_max2 = 0.8
Yxs2 = 0.02 # gSubstrate2/gDW 
tspan = (0., 80.)

# Override standard sampling times
sampling_times = collect(tspan[1]+10 : 6 : tspan[2]) # sample every 12 hour starting from hour 10

# Initial concentrations
init_glucose = 50 # mg / ml
init_glutamate = 1 # mg / ml

s_f1 = [100, 0, 0, 0, 0, 0, 0, 1, 0] # Feed medium 1, indexed the same as state variables, thus there are 100 mg of Glucose and 50 mg of Glutamine
s_f2 = [0, 0, 0, 0, 0, 0, 0, 10, 0] # Feed medium 2, indexed the same as state variables
feeding_times = sampling_times .+ 0.1 # Add 0.1 to avoid feeding at the same time as sampling
feeding_volumes1 = repeat([100], length(feeding_times))
feeding_volumes2 = repeat([5], length(feeding_times))
feeding_dict1 = Dict(zip(feeding_times, feeding_volumes1))
feeding_dict2 = Dict(zip(feeding_times, feeding_volumes2))

"""
Event handler function that simulates step feed addtion to the ODE system.
"""
function affect_add_feed!(integrator, feed_state_variable_idx::Int, feeding_dict::Dict, feed_concentrations::Array)
    volume = feeding_dict[integrator.t]
    
    integrator.u[5] += volume # add fed volume to the volume state variable
    integrator.u[feed_state_variable_idx] += volume # add fed volume to the volume state variable

    for idx in eachindex(feed_concentrations)
        integrator.u[idx] += volume * feed_concentrations[idx]
    end
    return 
end

# Create affect functions for each feed
affect_feed1!(integrator) = affect_add_feed!(integrator, 6, feeding_dict1, s_f1)
affect_feed2!(integrator) = affect_add_feed!(integrator, 7, feeding_dict2, s_f2)

# Setup callback, when the execute the affects
cb_samples = PresetTimeCallback(sampling_times, affect_sample_multiple_impulse_feeds!, filter_tstops = true) # creates new sampling callback
cb_feed1 = PresetTimeCallback(feeding_times, affect_feed1!, filter_tstops = true) # creates new sampling callback
cb_feed2 = PresetTimeCallback(feeding_times, affect_feed2!, filter_tstops = true) # creates new sampling callback

cbs = CallbackSet(cb_samples, cb_feed1, cb_feed2)

# Setting up ODE input
# make sure that the sampling and feeding times are included in the timesteps
save_ode_timesteps = unique(sort([feeding_times; sampling_times; LinRange(tspan[1], tspan[2], 1000)]))
state_variable_names = [:m_Glucose, :m_Biomass, :m_Product, :m_CO2, :v_Volume, :v_Feed_accum1, :v_Feed_accum2, :m_Glutamine, :m_CO2_gas]
init_cond = [init_glucose*V0, x0*V0, p0*V0, co2_0*V0, V0, 0, 0., init_glutamate*V0, 0.] # input for the model is masses not concentrations, accum feed at time zero is 0
sample_volume_dict = Dict(zip(sampling_times, repeat([sample_volume], length(sampling_times)))) # defines timepoints and sample volumes
ode_input_p = [Kc_s, Kc_s2, mu_max, mu_max2, Yxs, Yxs2, Yxp, Yxco2]

ode_func = ODEFunction(fedbatch_multiple_step_feeds!, syms=state_variable_names)
prob = ODEProblem(ode_func, init_cond, tspan, ode_input_p)
sol = solve(prob, Tsit5(), callback=cbs, saveat=save_ode_timesteps, abstol = 1e-12) # solve ode system with new sampling callback

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

# Adding the feed concentrations only for species which are fed
for idx in eachindex(state_variable_names)
    state_variable = state_variable_names[idx]
    metabolite = split(string(state_variable), "_")[2]

    if s_f1[idx] > 0 
        new_colname = string("c_", metabolite, "_feed1")
        insertcols!(df, new_colname => s_f1[idx])
    end

    if s_f2[idx] > 0 
        new_colname = string("c_", metabolite, "_feed2")
        insertcols!(df, new_colname => s_f2[idx])
    end
end
# adding the true growth rate
transform!(df, [:c_Glucose, :c_Glutamine] => ByRow((c_s1, c_s2) -> growth_rate_two_substrates(c_s1, c_s2, mu_max, mu_max2, Kc_s, Kc_s2)) => :mu_true)

# adding sampling information
insertcols!(df, 1, "sample_volume" => NaN)
df[[x in sampling_times for x in df.timestamp], "sample_volume"] .= sample_volume

# adding the parameter values to the dataframe
output_header = [:Kc_s, :Kc_s2, :mu_max, :mu_max2, :Yxs, :Yxs2, :Yxp, :Yxco2] 
insertcols!(df, 1, (output_header .=> ode_input_p)...)

# Substrate concentrations in feed

CSV.write(string("data/multiple_impulse_feed_process.csv"), df)