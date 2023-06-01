# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 23:03:57 2023

@author: uqtmccub
"""

import pandas as pd
import numpy as np
from .data_correction import pseudobatch_transform_pandas
from typing import Union, Literal, Tuple, List
import pathlib

def process_excel_template(
        file_name: Union[str, pathlib.Path], 
        phase: Literal['all', "1", "2"] = 'all', 
        outlier_removal: bool = False
    )-> Tuple[pd.DataFrame, Union[pd.DataFrame, pd.Series], List[str], List[str]]:
    """
    Reads in data from excel template and perform pseudo-batch correction. Excel/csv-type 
    files are typical ways to store bioreactor data on and offline.
    Measurements are populated into the excel file either programmatically or
    manually from a mix of on and off-line data. 

    Data can be quality controlled by the user with guidance from the spreadsheet.

    Parameters
    ----------
    file_name : str or pathlib.Path
        Filename and path to the file if outside of the directory.
    phase : {'all', '1', '2'}, optional
        Specifies the boundaries of rate calculations and the source of outlier information.
        Practically, this picks where the outlier information is taken from and records 
        phase bounds which can be considered later, e.g., restricting rate calculations 
        to a predefined phase. Default is 'all'.
    outlier_removal : bool, optional
        Flag to indicate whether to remove data points based on the Excel document. 
        Default is False.

    Usage
    -----
    The idea is to read in individual phases into their own data frame.
    Note: if outlier removal is set to True and phase is set to 'all', we will assume 
    that only values in Phase 1 on the template are to be considered outliers.

    May need to consider whether it is meaningful/worthwhile modifying the existing code 
    to only return average values across a specified phase, for example, the finite 
    differences method.

    Notes
    -----
    Created on Tue Mar 7 15:01:51 2023
    Author: uqtmccub
    """

    #Check function inputs
    assert phase in ["1","2","all", 1, 2], "phase is not set to an allowable value of '1', '2' or 'all'"
    assert type(outlier_removal) is bool, "outlier_removal may only be True or False"

    #Bioreactor data#
    
    df = pd.read_excel(file_name, 'Bioreactor_data',header=1)
    
    #Identify the biomass column
    bio_name = [col for col in df.columns if 'Biomass (' in col]
    assert len(bio_name) == 1, "The raw biomass column can't be identified. Only one column name beginning with 'Biomass (' is permitted"
    
    #Keep useful columns only and move biomass to the end
    col_list = ['Sample ID','Time (h)','Volume sampled (mL)','Cumulative Feed 1 (mL)','Cumulative Feed 2 (mL)','Cumulative Feed 3 (mL)','Volume at sample (mL)',bio_name[0]]
    
    #Before dropping all columns, we need to keep some if we want to retain phase and outlier removal information
    ###NOTE these lines of code depend on the excel columns remaining constant - a little lazy coding here.

    #Create a mapper column to translate times to sample IDs
    mapper = dict(zip(df.iloc[:,4],df.iloc[:,1]))
    if phase in ['all','1']:  #Note: we assume the phase definition points are always set in AF4 and AF5 in Bioreactor_data
        phase_bounds = pd.Series(df.iloc[[1,2],31], name = bio_name[0])
    else:
        phase_bounds = pd.Series(df.iloc[[1,2],39], name = bio_name[0])
    # Map time to sample ID                 
    phase_bounds = phase_bounds.map(mapper)   
        
    #If considering outliers or phases, retain additional information
    if outlier_removal == True or phase != 'all':
        mapping_df = df.iloc[:,[1,4,17]]

    #Drop other columns
    df = df[col_list]
    
    #Feed media concentations
    feeds = pd.read_excel(file_name, 'Feed definitions',header=0,nrows=3)
    #Clean up df to remove columns not needed
    feeds = feeds.drop(columns = feeds.columns[0])
    feeds = feeds.loc[:, ~feeds.columns.str.contains('^Unnamed')]
    
    #Raw concentration data
    measurements = pd.read_excel(file_name, 'Raw concentration data',header=0)
    #Remove any columns with no measurements
    measurements = measurements.dropna(axis=1, how='all')
    
    #We need to ensure that the column names for measurements are consistent with the feed column names
    
    #If measurement not in feed, we assume feed is 0 and provide warning.
    measured_cpds = measurements.columns[1:].to_list()
    #Add in biomass
    measured_cpds.insert(0,bio_name[0])
    
    #Generate feed list
    feed_cpds = feeds.columns[1:].to_list()
    
    #Calculate mismatches
    measure_mm = list(set(measured_cpds).difference(feed_cpds))
    
    if len(measure_mm) > 0:
        print("The following compounds are not defined in feeds but are measured. They are assumed to not be fed. Compounds: "+", ".join(measure_mm))
        for cpd in measure_mm:
            feeds[cpd] = 0    
    
    #If feed compound is not measured, possibly an error. We will assume an error as there is no need to define if not measured
    
    feed_mm = list(set(feed_cpds).difference(measured_cpds))
    assert len(feed_mm) == 0, "A substrate in the feed is not measured. Check if this is in error; if not then remove substrate from feed."
    
    #Measurement and feed definition lists are now known to contain the same headers, we need to ensure that the order is conserved between the two.
    feeds = feeds[measured_cpds]
    
    #Finally, we need to define the list of feed concentrations
    feed_conc = feeds.T.values.tolist()
    
    #As a further convenience, let's retain a list of whether a compound is fed or not. If fed, we will have to treat it different for rate calculations
    fed_cpds = feeds.any() #Returns a Boolean series if compound is in any of the feeds.
    fed_cpds = fed_cpds.index[fed_cpds.values == True].tolist()
    #For extra convenience we will take on a pseudo:
    fed_cpds = [cpd + "_pseudo" for cpd in fed_cpds]    
    
    #Merge df and measurement dataframes
    
    fedbatch_df_measurement = pd.merge(df, measurements, on ='Sample ID')
 
    #Prior to performing pseudobatch conversion, drop datapoints that have been flagged as outliers
    #We need a little bit of wrangling here - outlier values are based on corrected values in excel sheet. Will need to keep columns from bioreactor data and search across dataframes
    if outlier_removal == True or phase != 'all':
        mapper = dict(zip(mapping_df.iloc[:,2],mapping_df.iloc[:,0]))
        phase_data = pd.read_excel(file_name, 'Visualisations',header=0)
    if phase in ['all','1']:
        if outlier_removal == True:
            outliers = phase_data.iloc[40:43,1:]
        if phase == "1":
            bounds = phase_data.iloc[38:40,1:]    
      
    if phase == "2":
        if outlier_removal == True:
            outliers = phase_data.iloc[90:93,1:]
        bounds = phase_data.iloc[88:90,1:]  
    
    if outlier_removal == True:
        for col in outliers.columns:
            outliers[col]=outliers[col].map(mapper)
        #Censor flagged data - extract row info first
        fedbatch_df_measurement.loc[fedbatch_df_measurement['Sample ID'].isin(outliers.loc[:,col].values), col] = np.nan
    if phase != 'all':
        for col in bounds.columns:
           bounds[col]=bounds[col].map(mapper)

            
    #Perform pseudobatch calculation
    
    #Define new columns for pseudo concentrations   
    pseudo_cpds = [cpd + "_pseudo" for cpd in measured_cpds] 
    
    #Is it too restrictive to automatically create the pseudo conc. columns on import? For example, let's imagine we have
    # more than 3 feeds we need captured (unlikely) so we want control to do extra processing within Python
    #Run function   
    fedbatch_df_measurement[pseudo_cpds] = pseudobatch_transform_pandas(
        df=fedbatch_df_measurement,
        measured_concentration_colnames= measured_cpds,
        reactor_volume_colname='Volume at sample (mL)',
        accumulated_feed_colname=['Cumulative Feed 1 (mL)', 'Cumulative Feed 2 (mL)', 'Cumulative Feed 3 (mL)'],
        concentration_in_feed= feed_conc,
        sample_volume_colname='Volume sampled (mL)'
    )
    
    #As a final step, return all phase bounds
    #Combine the biomass and compound specific bounds and rename to match the pseudo annotation
    if phase != 'all':
        phase_bounds = pd.concat([phase_bounds.reset_index(drop=True),bounds.reset_index(drop=True)],axis=1).add_suffix('_pseudo')
        phase_bounds = phase_bounds[pseudo_cpds]

    #Now remove any columns not associated with measured compounds
    
    return fedbatch_df_measurement, phase_bounds, pseudo_cpds, fed_cpds