# Spatial Analysis

Spatial analysis can be an essential tool to evaluate non-point sources of CECs, as CECs are often applied with unknown rates, timing, and frequency, and transported 
through watersheds via various, often unobservable pathways. In this sense, spatial characteristics (e.g. agricultural and urban land cover in watershed) can be reliable
proxies for estimating observed CEC dynamics (Kiesling et al., 2019). To examine spatial relationships between land use and CECs, we first delineated watershed 
boundaries with StreamStats (USGS, 2016) using the sampling location as the watershed outlet. Within each watershed, we calculated several spatial characteristics. For 
land cover, we extracted the percentage of agricultural (summing hay/pasture and cultivated crops), urban (summing open space, low, medium, high developed cover), and 
forested (summing deciduous, evergreen, and mixed forests) land cover using the 2016 dataset from USGS National Land Cover Database (NLCD) (Jin et al., 2019; Homer et 
al., 2020; MRLCC, 2020). We also report septic tank use in each watershed as septic systems can leach CECs to the environment (Schaider et al., 2017), which can then be 
diverted to surface streams via baseflow. In addition, as proximity to Waste Water Treatment Plants (WWTPs) have been proven to exert influence on detection of CECs (Kolpin et al., 2002; Kolpin et 
al., 2004; Veach and Bernot, 2011; Vidal-Dorsch et al., 2012), we report number of WWTPs in each study watershed as well as sampling location proximity to the closest 
WWTP facility (i.e. hydrologic distance along the stream). 

We also explored empirical predictors beyond watershed aggregated land cover, as some CECs are less persistent in the environment and their contributing area might be 
smaller than the watershed scale. We focused on urban, pasture, and cultivated crops land covers as they are likely the main land use proxies associated with sources of 
pharmaceuticals, personal care products, and pesticides. We used a distance weighting model (Van Sickle and Burch Johnson, 2008), converting each watershed NLCD raster
cell to point data and calculating the distance between each raster point and sampling location for each watershed. Then, the distribution of land cover in each 
watershed as a function of proximity to the sampling location was analyzed for all three land cover types.  

Using these land cover distributions, we analyzed if land cover proximity to the sampling location impacts the ability to predict observed CEC numbers by type. We 
created an Ordinary Least Square (OLS) model for each land cover type to determine relationship between detected total number of CECs at each watershed and percent of 
land cover (cultivated crops, pasture, and urban) at a certain distance, with r2 values were used to report model accuracy. Land cover (cultivated crops, pasture, and 
urban) percentages at a given distance (from 0.3 to 9.7 km at 0.1 km distance increments from the sampling location) in all study watersheds were used as independent 
variables to predict number of detected pharmaceuticals, personal care products, pesticides and their transformation products in each sampled location. We also used the 
same procedure to assess empirical relationships between incremented land cover percentages and concentrations of the most frequently detected CECs. Only CEC models 
that exhibited r2 values exceeding 0.3 were reported in this study, to focus on the compounds that show the strongest empirical relationships with watershed land cover.

References: \

