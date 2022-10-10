*Env_corr_analysis.R* - code that analyzes correlation patterns between different environmental variables and CEC concentrations. 
*Env_data_analysis.R* - code that processes environemntal variables and calculates different averages to approximate immediate and longer term climate conditions, as
CECs remain in the system for multiple days.
*CEC_aver_dataset.csv* - contaminant concentration values (will be added after this work is published)

# Environmental Data

We used hydrological observations to (1) characterize similarities and differences in hydrologic functioning across study sites and (2) to evaluate instantaneous and 
long-term hydrologic conditions during sampling periods. Hydrological conditions associated with each sampling date were used to evaluate how hydrological signatures 
differed in each site as well as between sampling dates. Daily streamflow values were obtained from USGS for 2018 and 2019 and normalized to mm per day. Precipitation 
(mm) and air temperature (C) data was acquired from DAYMET (Daily Surface Weather and Climatological Summaries) (Thornton et al., 1997; 2021). 

Hydrological characterizations were derived for annual and sub-annual timescales to identify short-term and long-term changes at each sampling point. Since certain CEC 
contributions and transport are often linked to runoff (Kolpin et al., 2004; Gray et al., 2017; Fairbairn et al., 2018), we performed baseflow separation on daily 
streamflow values using the EcoHydRology package (Fuka et al., 2018) to analyze daily baseflow and runoff contributions to the stream. As CECs can survive in the 
environment for multiple days before degradation (Ashraf, 2017), we analyzed streamflow and precipitation contributions over different time periods. We calculated a 
set of hydrologic signatures that capture hydrologic conditions at each sampling location to analyze how hydrologic variables impact observed CEC dynamics. Using daily 
precipitation and streamflow observations, we extracted same-day, 4-day average, 7-day average, and 14-day average streamflow and precipitation values, which were 
calculated to reflect the hydrologic conditions prior to sampling event. 

We correlated each of these hydrologic variables to CEC concentrations of the most commonly detected compounds to identify which environmental variables may empirically 
explain temporal CEC concentration changes.  The goal of this exercise was to better understand the impact hydroclimate has on CEC dynamics in surface waters. 
Relationships were quantified with the spearman rank correlation coefficient.
 
\
References:\
Thornton, P.E., Running, S.W., White, M.A. 1997. Generating surfaces of daily meteorological variables over large regions of complex terrain. Journal of Hydrology 190: 
214 - 251. https://doi.org/10.1016/S0022-1694(96)03128-9 \
Thornton, P. E., R. Shrestha, M. Thornton, S.-C. Kao, Y. Wei, and B. E. Wilson. 2021. Gridded daily weather data for North America with comprehensive uncertainty 
quantification. Scientific Data 8. https://doi.org/10.1038/s41597-021-00973-0 \
Kolpin, D.W., Skopec, M., Meyer, M.T., Furlong, E.T., Zaugg, S.D. (2004). Urban contribution of pharmaceuticals and other organic wastewater contaminants to streams 
during differing flow conditions. Sci. Total Environ. 328, 119–130. https://doi.org/10.1016/j.scitotenv.2004.01.015 \
Gray, J.L., Borch, T., Furlong, E.T., Davis, J., Yager, T.J., Yang, Y.-Y., Kolpin, D.W. (2017). Rainfall-runoff of anthropogenic waste indicators from agricultural 
fields applied with municipal biosolids. Sci. Total Environ. 580, 83–89. https://doi.org/10.1016/j.scitotenv.2016.03.033 \
Fairbairn, D.J., Elliott, S.M., Kiesling, R.L., Schoenfuss, H.L., Ferrey, M.L., Westerhoff, B.M. (2018). Contaminants of emerging concern in urban stormwater: 
spatiotemporal patterns and removal by iron-enhanced sand filters (IESFs). Water Res. 145, 332–345. https://doi.org/10.1016/j.watres.2018.08.020 \
Fuka, D. T., Walter, M. T., Archibald, J. A., Steenhuis, T. S., Easton, Z. M. (2018). A Community Modeling Foundation for Eco-Hydrology. R package version 0.4.12 \
Ashraf, M.A. (2017). Persistent organic pollutants (POPs): a global issue, a global challenge. Environ Sci Pollut Res 24, 4223–4227.
https://doi.org/10.1007/s11356-015-5225-9 
