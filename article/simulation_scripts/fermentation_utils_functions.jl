# This file contains utility functions that is used during the simulation 
# of fedbatch processes.
"""
Calculates the yield coefficient of CO2 assuming the carbon which isn't used for either product or biomass is turned into CO2. 
Furthermore the number of carbons in the molecules (biomass, substrate, product and CO2) is currently hard coded.
"""
function calculate_co2_yield_gg(Yxp, Yxs, MW_x, MW_s, MW_p, MW_co2)
    # solve carbon balance to set theoretical CO2 yield
    ## we need the yield coef's in substrate basis
    Ysp = Yxp/Yxs
    Ysx = Yxs^(-1)

    ## Calculate the C-molarmasses
    cM_s = MW_s / 6 # gGLC / cmol
    cM_p = MW_p / 8 # gProduct / cmol
    cM_x = MW_x / 1 
    cM_co2 = MW_co2 / 1 # gCO2 / cmol

    # convert yield from g/g to cmol/cmol
    Ysx_cmol = Ysx * (cM_s/cM_x)
    Ysp_cmol = Ysp * (cM_s/cM_p)

    c_balance = -1 + Ysx_cmol + Ysp_cmol # does the carbon balance

    Ysco2_cmol = -1 * c_balance # sets co2 yield to the account for the remaining carbon
    Yxco2_cmol = Ysco2_cmol / Ysx_cmol # (cmol_co2 / cmol_s) / (cmol_dw / cmol_s)
    Yxco2 = Yxco2_cmol * (cM_co2/cM_x) 

    return Yxco2
end

"""
Convert a yield coefficient of the form Y_{g2/g1} to the equivalent yield coefficient in cmol/cmol
"""
function gg_to_cmolcmol(yield_coef, MW1, MW2, n_c1, n_c2)
    return yield_coef * (MW1/MW2) * (n_c1/n_c2)
end

"""
Convert a yield coefficient of the form Y_{cmol2/cmol1} to the equivalent yield coefficient in g/g
"""
function cmolcmol_to_gg(yield_coef, MW1, MW2, n_c1, n_c2)
    return yield_coef * (MW2/MW1) * (n_c2/n_c1)
end

