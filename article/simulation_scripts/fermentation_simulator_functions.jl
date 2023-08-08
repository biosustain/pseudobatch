# This file contains all the functions used to simulate the fed-batch fermentations.
# Thus here the ODE systems, utility functions and event handling functions are defined.

############################## UNITILITY FUNCTIONS ###########################

"""
Calcualtes the initial feeding rate to achive a constants growth rate.
This equation is obtained from eq. 9.53 from Bioreaction Engineering Principples by Villadsen, et. al.
"""
function calc_F0(Yxs, x0, V0, mu0, s_f, s0)
    F0 = (Yxs * x0 * V0 * mu0 / (s_f - s0))
    return F0
end

"""
Calculates feeding rate at any give time point of a exponential feed profile.
This equation is obtained from eq. 9.53 from Bioreaction Engineering Principples by Villadsen, et. al.
"""
function v_feed(t, F0, mu0) 
    return F0 * exp(mu0*t)
end


"""
Monods equation for growth rate as a function of substrate concentration.
This equation is obtained from eq. 7.16 from Bioreaction Engineering Principples by Villadsen, et. al.
"""
function monod_kinetics(c_s, mu_max, Kc_s)
    if c_s < 0
        return 0
    end 
    return mu_max * c_s / (c_s + Kc_s)
end


"""
Calculates the substrate concentration required to achieve a certain growth rate.
The equation is obtained via isolating c_s from the monod kinetics equation 
(eq. 7.16 in Bioreaction Engineering Principples by Villadsen, et. al.).
"""
function monod_ss_substrate_conc(mu, mu_max, Kc_s)
    # Calcualtes the substrate concentration required to achieve a certain growth rate
    return (mu * Kc_s) / (mu_max - mu)
end


"""
Calculated the growth rate with product inhibition. This equation is obtained
from eq. 7.19 from Bioreaction Engineering Principples by Villadsen, et. al.
"""
function monod_with_product_inhibition(c_s, c_p, mu_max, Kc_s, Ki_p)
    return mu_max * (c_s / (c_s + Kc_s)) * (1 + (c_p/Ki_p))^(-1)
end


"""
Implements the growth rate as funciton of two substrates concentrations. This 
equation is obtained eq. 7.22 from Bioreaction Engineering Principples by Villadsen, et. al.
"""
function growth_rate_two_substrates(c_s1, c_s2, mu_max1, mu_max2, Kc_s1, Kc_s2)
    if c_s1 < 0 || c_s2 < 0
        return 0
    end 
    return (mu_max1 * mu_max2 * c_s1 * c_s2) / ((c_s1 + Kc_s1) * (c_s2 + Kc_s2))
end

"""
Simple calculation the evaporation rate based only based a scalar and the concentration
of the component in the liquid. return 0 if the concentration is negative.
"""
function evaporation_rate(c, k)
    if c < 0
        return 0
    end
    return c * k
end

"""
Calculates the volumetric gas to liquid transfer rate based on the saturation
concentration of the component, liquid concentration and the mass transfer coefficient (kla).
The equation is obtained from eq. 10.1 from Bioreaction Engineering Principples by 
Villadsen, et. al.
"""
function volumetric_gas_to_liquid_transfer_rate(c_sat, c_l, kla)
    return kla * (c_sat - c_l)
end


############################ ODE SYSTEMS #######################################
"""
Defines the ODE system for a fed-batch fermentation operated with an exponential
feed profile. The growth of the microorganism is modelled through monods equation.

The system holds 7 state variables with are defined as follows:
    dudt[1] : mass of substrate in the liquid
    dudt[2] : mass of biomass in the liquid
    dudt[3] : mass of product in the liquid
    dudt[4] : mass of co2 in liquid
    dudt[5] : volume of liquid
    dudt[6] : feeding rate
    dudt[7] : mass of co2 in off-gas

The systems is build under the assumption that the CO2 evaporates instantly
after production.
"""
function fedbatch!(dudt, u, p, t)
    Kc_s, mu_max, Yxs, Yxp, Yxco2, F0, mu0, s_f = p      

    # Growth kinetics consider moving this as a seperate function
    c_s = u[1]/u[5] 
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

"""
Defines the ODE system for a fed-batch fermentation operated with an exponential
feed profile. The growth of the microorganism is modelled through monods equation 
with product inhibition.

The system holds 7 state variables with are defined as follows:
    dudt[1] : mass of substrate in the liquid
    dudt[2] : mass of biomass in the liquid
    dudt[3] : mass of product in the liquid
    dudt[4] : mass of co2 in liquid
    dudt[5] : volume of liquid
    dudt[6] : feeding rate
    dudt[7] : mass of co2 in off-gas

The systems is build under the assumption that the CO2 evaporates instantly
after production.
"""
function fedbatch_prod_inhib!(dudt, u, p, t)
    Kc_s, mu_max, Yxs, Yxp, Yxco2, F0, mu0, s_f, Ki_p = p      
    
    # Growth kinetics consider moving this as a seperate function
    c_s = u[1]/u[5]  
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

"""
Defines the ODE system for a fed-batch fermentation operated step feed (e.g. manually 
pipeting a certain volume into the reactor). The growth of the microorganism is modelled 
through monods equation with two limiting substrates. This system models a fermentation 
of e.g. Chinese hamster ovary cells (CHO cells) which are grown on a mixture of glucose
and glutamine. Furthermore the system utilises two different feed mediums which can be 
added independently of each other.

The system holds 9 state variables with are defined as follows:
    dudt[1] : mass of substrate in the liquid
    dudt[2] : mass of biomass in the liquid
    dudt[3] : mass of product in the liquid
    dudt[4] : mass of co2 in liquid
    dudt[5] : volume of liquid
    dudt[6] : feeding rate of feed1
    dudt[7] : feeding rate feed2
    dudt[8] : mass of substrate2 in the liquid
    dudt[9] : mass of co2 in off-gas

The systems is build under the assumption that the CO2 evaporates instantly
after production.
"""
function fedbatch_multiple_step_feeds!(dudt, u, p, t)
    Kc_s1, Kc_s2, mu_max1, mu_max2, Yxs1, Yxs2, Yxp, Yxco2 = p      

    # Growth kinetics consider moving this as a seperate function
    c_s1 = u[1]/u[5]  # calculating concentration
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

"""
Defines the ODE system for a fed-batch fermentation operated with an exponential
feed profile. The growth of the microorganism is modelled through monods equation 
with product inhibition. The product is volatile and evaporates from the reactor.
The product evaporation is modelled through a first order reaction.

The system holds 7 state variables with are defined as follows:
    dudt[1] : mass of substrate in the liquid
    dudt[2] : mass of biomass in the liquid
    dudt[3] : mass of product in the liquid
    dudt[4] : mass of co2 in liquid
    dudt[5] : volume of liquid
    dudt[6] : feeding rate
    dudt[7] : mass of co2 in off-gas
    dudt[8] : mass of product in off-gas

The systems is build under the assumption that the CO2 evaporates instantly
after production.
"""
function fedbatch_prod_inhib_volatile!(dudt, u, p, t)
    Kc_s, mu_max, Yxs, Yxp, Yxco2, F0, mu0, s_f, Ki_p, evaporation_rate_product, Yxo2, O2_flow_rate_in, kla, c_o2_sat = p      
    
    # Growth kinetics consider moving this as a seperate function
    c_s = u[1]/u[5]  
    c_p = u[3]/u[5]
    c_o2 = u[10]/u[5]
    mu = monod_with_product_inhibition(c_s, c_p, mu_max, Kc_s, Ki_p)
    qt_o2 = volumetric_gas_to_liquid_transfer_rate(c_o2_sat, c_o2, kla)

    dudt[1] = -Yxs * mu * u[2] + v_feed(t, F0, mu0) * s_f
    dudt[2] = mu * u[2]
    dudt[3] = Yxp * mu * u[2] - evaporation_rate(u[3], evaporation_rate_product)
    dudt[4] = Yxco2 * mu * u[2] -  Yxco2 * mu * u[2] # co2 production - co2 evaporation
    dudt[5] = v_feed(t, F0, mu0) # volume
    dudt[6] = v_feed(t, F0, mu0) # feed used to integrate feed_accum
    dudt[7] =  Yxco2 * mu * u[2] # all co2 evaporates immidately into off-gas analyzer
    dudt[8] = evaporation_rate(u[3], evaporation_rate_product) # product evaporation
    dudt[9] = O2_flow_rate_in # O2 flow rate into the reactor
    dudt[10] = u[5] * volumetric_gas_to_liquid_transfer_rate(c_o2_sat, c_o2, kla) - Yxo2 * mu * u[2] # O2 liquid concentration
    dudt[11] = u[5] * volumetric_gas_to_liquid_transfer_rate(c_o2_sat, c_o2, kla) # O2 that flows into the off-gas analyzer


    return dudt
end

############################# EVENT HANDLING / CALLBACKS #####################
"""
Calculates the amount of mass removed when a give volume is sampled from the
fermenter.
"""
function remove_mass_through_sampling(mass_state_variable, volume_state_variable, sample_volume)
    return mass_state_variable - (mass_state_variable / volume_state_variable) * sample_volume
end

"""
Defines an event handling function which that modifies the state variables when the event is 
triggered. This function is used to simulate the sampling of the fermenter. Removes mass from
the fermenter and adjusts the feed rate to account for the removed biomass.

The adjustment of the feed ensures that growth rate is kept constant throughout the fermentation.
"""
function affect_sample!(integrator)
    sample_vol = sample_volume_dict[integrator.t]
    integrator.u[1] = remove_mass_through_sampling(integrator.u[1], integrator.u[5], sample_vol)
    integrator.u[2] = remove_mass_through_sampling(integrator.u[2], integrator.u[5], sample_vol)
    integrator.u[3] = remove_mass_through_sampling(integrator.u[3], integrator.u[5], sample_vol)
    integrator.u[4] = remove_mass_through_sampling(integrator.u[4], integrator.u[5], sample_vol)
    integrator.u[5] -= sample_vol

    integrator.p[6] *= integrator.u[5]/(integrator.u[5]+sample_vol) # adjusting feed to account for removed volume
    
end

"""
Defines an event handling function which that modifies the state variables when the event is
triggered. This function is used to simulate the sampling of the fermenter. Removes mass from
the fermenter. This does not adjust the feeding rate because the feeding is predetermined 
volume additions at certain time points and not attempting a reach a specific growth rate.
"""
function affect_sample_multiple_impulse_feeds!(integrator)
    sample_vol = sample_volume_dict[integrator.t]
    integrator.u[1] = remove_mass_through_sampling(integrator.u[1], integrator.u[5], sample_vol)
    integrator.u[2] = remove_mass_through_sampling(integrator.u[2], integrator.u[5], sample_vol)
    integrator.u[3] = remove_mass_through_sampling(integrator.u[3], integrator.u[5], sample_vol)
    integrator.u[4] = remove_mass_through_sampling(integrator.u[4], integrator.u[5], sample_vol)
    integrator.u[5] -= sample_vol
    integrator.u[8] = remove_mass_through_sampling(integrator.u[8], integrator.u[5], sample_vol)
    
end

function affect_sample_volatile_compounds!(integrator)
    sample_vol = sample_volume_dict[integrator.t]
    integrator.u[1] = remove_mass_through_sampling(integrator.u[1], integrator.u[5], sample_vol)
    integrator.u[2] = remove_mass_through_sampling(integrator.u[2], integrator.u[5], sample_vol)
    integrator.u[3] = remove_mass_through_sampling(integrator.u[3], integrator.u[5], sample_vol)
    integrator.u[4] = remove_mass_through_sampling(integrator.u[4], integrator.u[5], sample_vol)
    integrator.u[5] -= sample_vol
    integrator.u[10] = remove_mass_through_sampling(integrator.u[10], integrator.u[5], sample_vol)

    integrator.p[6] *= integrator.u[5]/(integrator.u[5]+sample_vol) # adjusting feed to account for removed volume
    
end