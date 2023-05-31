# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 23:33:34 2023

@author: uqtmccub
"""

# To do - attempt to follow Viktor's tutorials using the excel output in order to create a new, how to use the excel template tutorial

# To do - check excel tempalte is properly formatted to calculate DFs for first few rows. Is there an essential requirement for first row to not be considered? Seems to be hard-coded that not used?
# Check what happens if want to start from 0 but multi-matches for first two samples?



#First lets populate an excel with the same data as used by the multiple feed example

# Then run through all tutes generating figures, calcualting rates, volumetric productivities, etc while utilising the extra features of my approach, such as restriction to times for the fitting of rates
# 

fedbatch_df_measurement, phase_bounds, pseudo_cpds, fed_cpds = process_template('Pseudo_batch_template_test.xlsx', phase="1", outlier_removal = False)

#Typically we would then be interested in calculating rates. Defining the generic function to fit rates and yields:

def fit_ols_model(formula_like: str, data: pd.DataFrame) -> sm.regression.linear_model.RegressionResultsWrapper:
    y, X = dmatrices(formula_like, data)
    model = sm.OLS(endog=y, exog=X)
    res = model.fit()
    return res

# In previous tutorials, plotting the data and fits has been an important step. Here we assume that this has already occurred using Excel
# and we are instead concentrating on pulling the data into Python to perform more advanced data manipulation or statistics

#First let's add a log-trasnformed pseudo-biomass column to calculate growth rate

#We assume here that all time points have a corresponding biomass measurement.
log_bio_name = "np.log(" + pseudo_cpds[0] + ")"
fedbatch_df_measurement[log_bio_name] = np.log(fedbatch_df_measurement[pseudo_cpds[0]])


#Next let's reduce the fedbatch measurement dataframe to consider only the phase of interest if we are considering phases. We will do this by censoring values that fall outside of the phase range
fedbatch_df_measurement_phase = fedbatch_df_measurement
#Find matching rows to Sample ID
if isinstance(phase_bounds, pd.DataFrame):
    phase_bounds = phase_bounds.rename(columns = {pseudo_cpds[0]: log_bio_name}) #update name in phase_bounds
    for col in phase_bounds:
        row_locs = fedbatch_df_measurement_phase.index[fedbatch_df_measurement_phase['Sample ID'].isin(phase_bounds[col])].to_list()
        valid_indices = fedbatch_df_measurement_phase.loc[row_locs[0]:row_locs[1]].index.to_list() 
        fedbatch_df_measurement_phase.loc[~fedbatch_df_measurement_phase.index.isin(valid_indices),col] = np.nan
    
# First, let's calculate growth rate. First let's add log-transformed data to the fedbatch_df_measurement dataframe:
# To keep generic column name structure and avoid having to sacrifice human readable names for Python friendly names we will use patsy
# Noting biomass is always first in pseudo_cpds:

model_fit = "Q('" + log_bio_name + "')" + " ~ " + "Q('Time (h)')" #using patsy to handle disgusting variable names and hard coding the time column name in line with current spreadsheet template
mask = ~pd.isnull(fedbatch_df_measurement_phase[log_bio_name]) #create a mask to select for rows containing non NaN values
res_mu_hat_pseudo = fit_ols_model(model_fit, fedbatch_df_measurement_phase)    
growth_rate = res_mu_hat_pseudo.params[1]    

#We can also calculate yields. Store results iteratively. Note that here we will demonstrate a slightly different way of calculating yields
# to the previous tutorials (see Tutorial 2). Here we will assume that there may be a large number of compounds in the spreadsheet. When feeding amino acids
# it could be the case that some of the amino acids are produced transiently (or entirely) across the culture. Rather than needing to check for
# net consumption/production and doing any additional data permutaiton, we will simply use the direct pseudo concentrations (which can give negative yields)
# and use this information to correclty assign rate directions and then adjsut the yield sign as required.

parameter_df = pd.DataFrame(columns=['Compound','Yield','Rate'])

#First add in growth rate:
parameter_df.loc[0] = [pseudo_cpds[0],1,growth_rate]

for cpd in pseudo_cpds[1:]:
    mask = ~pd.isnull(fedbatch_df_measurement_phase[cpd]) #create a mask to select for rows containing non NaN values
    model_fit = "Q('" + cpd + "')" + " ~ " + "Q('" + pseudo_cpds[0] +"')" #using patsy to handle disgusting variable names and hard coding the time column name in line with current spreadsheet template
    res_yield = fit_ols_model(model_fit, fedbatch_df_measurement_phase[mask])
    #Append to DF. We will apply the convention of yields always being positive here
    parameter_df.loc[len(parameter_df)] = [cpd,abs(res_yield.params[1]),res_yield.params[1]*growth_rate]
    
#While much of the above could be done in the Excel template already, more advanced data processing is 
# more readily achieved in Python. For example, implementation of the finite differences method.
#Here we will perform this on the full data set - simply because boundary points will have insufficient information


#No need - the mass calculations and the consumed compound calculations are a lot of mess around for nothing. Just rely on pseudo concentrations
#Then we are really interested in volumetric rates and yields and automating the propagation of this across the dataframe

# Need to 1) automate rates and yields calculations
# Need to 2) automate specific consumption and production rate calculations

#Viktor feedback - why muck around with mass/pseudomass and consumed values etc? Can be somewhat confusing. Only thing you avoid is the negative yield coefficient...
#Also note no volumetric producitivity implementation



res_yxp_corrected = fit_ols_model(formula_like = "Q('c_Glucose_pseudo') ~ Q('Biomass (OD)_pseudo')", data= fedbatch_df_measurement_phase)

import numpy as np
from patsy import dmatrices
import statsmodels.api as sm



def finite_difference_derivative(df: pd.DataFrame, x_colname: str, y_colname: str)->np.ndarray:
    x = df[x_colname].to_numpy()
    y = df[y_colname].to_numpy()
    return np.gradient(y, x)


#Write a code to iterate through and calcualte volumetric productivities, rates etc. based on phase information
import matplotlib.pyplot as plt


mu_hat = finite_difference_derivative(fedbatch_df_measurement, "Time (h)", "Biomass (OD)_pseudo") / fedbatch_df_measurement["Biomass (OD)_pseudo"]
growth_rate_fig, growth_rate_ax = plt.subplots(figsize=(6, 5))
#growth_rate_ax.plot(fedbatch_df_measurement["Time (h)"], fedbatch_df_measurement.mu_true, label="Growth rate (simulation)", color = "grey")
growth_rate_ax.scatter(fedbatch_df_measurement["Time (h)"], mu_hat, label="Estimated growth rate (pseudo batch)", color = "C1")
growth_rate_ax.set_title("Growth rate")
growth_rate_ax.set_xlabel("Feeding time [h]")
growth_rate_ax.set_ylabel("Growth rate [1/h]")
growth_rate_ax.legend()