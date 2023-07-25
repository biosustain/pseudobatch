############################## UNITILITY FUNCTIONS ###########################
function calc_F0(Yxs, x0, V0, mu0, s_f, s0)
    F0 = (Yxs * x0 * V0 * mu0 / (s_f - s0))
    return F0
end

v_feed(t, F0, mu0) = F0 * exp(mu0*t)
integral_of_v_feed(t, F0, mu0) = F0 * exp(mu0*t) / mu0

V(t, V0, F0, mu0) = V0 + integral_of_v_feed(t, F0, mu0) - integral_of_v_feed(0, F0, mu0)

function monod_kinetics(c_s, mu_max, Kc_s)
    if c_s < 0
        return 0
    end 
    return mu_max * c_s / (c_s + Kc_s)
end

function monod_ss_substrate_conc(mu, mu_max, Kc_s)
    # Calcualtes the substrate concentration required to achieve a certain growth rate
    return (mu * Kc_s) / (mu_max - mu)
end

monod_with_product_inhibition = function(c_s, c_p, mu_max, Kc_s, Ki_p)
    return mu_max * (c_s / (c_s + Kc_s)) * (1 + (c_p/Ki_p))^(-1)
end

function growth_rate_two_substrates(c_s1, c_s2, mu_max1, mu_max2, Kc_s1, Kc_s2)
    """
    Implements the growth rate as funciton of two substrates concentrations. This 
    equation is obtained eq. 7.22 from Bioreaction Engineering Principples by Villadsen, et. al.
    """
    if c_s1 < 0 || c_s2 < 0
        return 0
    end 
    return (mu_max1 * mu_max2 * c_s1 * c_s2) / ((c_s1 + Kc_s1) * (c_s2 + Kc_s2))
end

############################ ODE SYSTEMS #######################################
function fedbatch!(dudt, u, p, t)
    Kc_s, mu_max, Yxs, Yxp, Yxco2, F0, mu0, s_f = p      

    # Refacotring ideas 
    ### reorganise dudt, so generic ones are first, e.g. volume and feed

    # sV = u[1] moles susbtrate
    # xV = u[2]
    # pV = u[3]

    # Growth kinetics consider moving this as a seperate function
    c_s = u[1]/u[5]  #V(t, V0, F0, mu0) # calculating concentration
    mu = monod_kinetics(c_s, mu_max, Kc_s)

    dudt[1] = -Yxs * mu * u[2] + v_feed(t, F0, mu0) * s_f
    dudt[2] = mu * u[2]
    dudt[3] = Yxp * mu * u[2]
    dudt[4] = Yxco2 * mu * u[2] - Yxco2 * mu * u[2] # co2 production - co2 evaporation
    dudt[5] = v_feed(t, F0, mu0) # volume
    dudt[6] = v_feed(t, F0, mu0) # feed used to integrate feed_accum
    dudt[7] = Yxco2 * mu * u[2] # all co2 evaporates immidately into off-gas analyzer


    return dudt
end

function fedbatch_prod_inhib!(dudt, u, p, t)
    Kc_s, mu_max, Yxs, Yxp, Yxco2, F0, mu0, s_f, Ki_p = p      

    # Refacotring ideas 
    ### reorganise dudt, so generic ones are first, e.g. volume and feed

    # sV = u[1] moles susbtrate
    # xV = u[2]
    # pV = u[3]

    # Growth kinetics consider moving this as a seperate function
    c_s = u[1]/u[5]  #V(t, V0, F0, mu0) # calculating concentration
    c_p = u[3]/u[5]
    mu = monod_with_product_inhibition(c_s, c_p, mu_max, Kc_s, Ki_p)

    dudt[1] = -Yxs * mu * u[2] + v_feed(t, F0, mu0) * s_f
    dudt[2] = mu * u[2]
    dudt[3] = Yxp * mu * u[2]
    dudt[4] = Yxco2 * mu * u[2] - Yxco2 * mu * u[2] # co2 production - co2 evaporation
    dudt[5] = v_feed(t, F0, mu0) # volume
    dudt[6] = v_feed(t, F0, mu0) # feed used to integrate feed_accum
    dudt[7] = Yxco2 * mu * u[2] # all co2 evaporates immidately into off-gas analyzer


    return dudt
end

function fedbatch_multiple_step_feeds!(dudt, u, p, t)
    Kc_s1, Kc_s2, mu_max1, mu_max2, Yxs1, Yxs2, Yxp, Yxco2 = p      

    # Refacotring ideas 
    ### reorganise dudt, so generic ones are first, e.g. volume and feed

    # sV = u[1] moles susbtrate
    # xV = u[2]
    # pV = u[3]

    # Growth kinetics consider moving this as a seperate function
    c_s1 = u[1]/u[5]  #V(t, V0, F0, mu0) # calculating concentration
    c_s2 = u[8]/u[5]
    mu = growth_rate_two_substrates(c_s1, c_s2, mu_max1, mu_max2, Kc_s1, Kc_s2)

    dudt[1] = -Yxs1 * mu * u[2] 
    dudt[2] = mu * u[2]
    dudt[3] = Yxp * mu * u[2]
    dudt[4] = Yxco2 * mu * u[2] - Yxco2 * mu * u[2] # co2 production - co2 evaporation
    dudt[5] = 0 # volume is only changed through sampling or step feeds
    dudt[6] = 0 # Feed1 only changed through feeding
    dudt[7] = 0 # Feed2 only changed through feeding
    dudt[8] = -Yxs2 * mu * u[2] 
    dudt[9] = Yxco2 * mu * u[2] # all co2 evaporates immidately into off-gas analyzer


    return dudt
end

############################# EVENT HANDLING / CALLBACKS #####################
function condition_sampling_time(sampling_times)
    # function that sets conditions which triggers the event
    return condition(u,t,integrator) = t ∈ sampling_times
end

function remove_mass_through_sampling(mass_state_variable, volume_state_variable, sample_volume)
    return mass_state_variable - (mass_state_variable / volume_state_variable) * sample_volume
end

function sample_volume_from_volume_fraction(volume_state_variable, sample_volume_fraction)
    return volume_state_variable * sample_volume_fraction
end

function affect_sample!(integrator)
    sample_vol = sample_volume_dict[integrator.t]
    integrator.u[1] = remove_mass_through_sampling(integrator.u[1], integrator.u[5], sample_vol)
    integrator.u[2] = remove_mass_through_sampling(integrator.u[2], integrator.u[5], sample_vol)
    integrator.u[3] = remove_mass_through_sampling(integrator.u[3], integrator.u[5], sample_vol)
    integrator.u[4] = remove_mass_through_sampling(integrator.u[4], integrator.u[5], sample_vol)
    integrator.u[5] -= sample_vol

    integrator.p[6] *= integrator.u[5]/(integrator.u[5]+sample_vol) # adjusting feed to account for removed volume
    
end

function affect_sample_volume_fraction!(integrator)
    sample_vol = sample_volume_from_volume_fraction(integrator.u[5], sample_volume_fraction) # calculates the sample volume to constantly sample the same volume fraction
    integrator.u[1] = remove_mass_through_sampling(integrator.u[1], integrator.u[5], sample_vol)
    integrator.u[2] = remove_mass_through_sampling(integrator.u[2], integrator.u[5], sample_vol)
    integrator.u[3] = remove_mass_through_sampling(integrator.u[3], integrator.u[5], sample_vol)
    integrator.u[4] = remove_mass_through_sampling(integrator.u[4], integrator.u[5], sample_vol)
    integrator.u[5] -= sample_vol

    integrator.p[6] *= integrator.u[5]/(integrator.u[5]+sample_vol) # adjusting feed to account for removed volume

    sample_volume_dict[integrator.t] = sample_vol # logs the sample volume and time
end

function affect_sample_volume_fraction_no_log!(integrator, sample_volume_fraction)
    sample_vol = sample_volume_from_volume_fraction(integrator.u[5], sample_volume_fraction) # calculates the sample volume to constantly sample the same volume fraction
    integrator.u[1] = remove_mass_through_sampling(integrator.u[1], integrator.u[5], sample_vol)
    integrator.u[2] = remove_mass_through_sampling(integrator.u[2], integrator.u[5], sample_vol)
    integrator.u[3] = remove_mass_through_sampling(integrator.u[3], integrator.u[5], sample_vol)
    integrator.u[4] = remove_mass_through_sampling(integrator.u[4], integrator.u[5], sample_vol)
    integrator.u[5] -= sample_vol

    integrator.p[6] *= integrator.u[5]/(integrator.u[5]+sample_vol) # adjusting feed to account for removed volume

end

function affect_sample_multiple_impulse_feeds!(integrator)
    sample_vol = sample_volume_dict[integrator.t]
    integrator.u[1] = remove_mass_through_sampling(integrator.u[1], integrator.u[5], sample_vol)
    integrator.u[2] = remove_mass_through_sampling(integrator.u[2], integrator.u[5], sample_vol)
    integrator.u[3] = remove_mass_through_sampling(integrator.u[3], integrator.u[5], sample_vol)
    integrator.u[4] = remove_mass_through_sampling(integrator.u[4], integrator.u[5], sample_vol)
    integrator.u[5] -= sample_vol
    integrator.u[8] = remove_mass_through_sampling(integrator.u[8], integrator.u[5], sample_vol)
    
end