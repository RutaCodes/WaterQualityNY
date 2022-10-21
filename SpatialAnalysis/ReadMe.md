*LC_developed_allWS_by_dist.csv; LC_crop_allWS_by_dist.csv; LC_pasture_allWS_by_dist.csv; LC_agr_allWS_by_dist.csv* - datasets with developed, agricultural, pasture, and cultivated crops land cover percentages at different distances from sampling locations. \
\
*LC_models_pred_conc.R* - develops models to predict CECs concentrations using developed, pasture, and cultivated crops land cover. Separate models are developed for 7 compounds.\
\
*LC_models_pred_numb.R* - develops models to predict number of detected CECs by type using developed, pasture, and cultivated crops land cover. Two methods are compared: 1) using whole watershed scale land cover percentages and 2) developing models using land cover percentages at varying distances from sampling location a.k.a. distance weighted models. \
\
*Land_cov_distrib.R* - calculates land cover percentage distributions at different distances from sampling locationthat is later used to develop distance weighted 
models to predict detected number of compounds.\
\
*NR_vs_landcover.R* - land cover data and number of compounds found by type in each study site.

*Total_CEC_Nrs.csv* - number of compounds found in each study site.

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

\
**References:** \
-Kiesling, R. L., Elliott, S. M., Kammel, L. E., Choy, S. J., Hummel, S. L. (2019). Predicting the occurrence of chemicals of emerging concern in surface water and sediment across the U.S. portion of the Great Lakes Basin. Sci. Total Environ. 651, 838 – 850. https://doi.org/10.1016/j.scitotenv.2018.09.201 \
-U.S. Geological Survey, (2016). The StreamStats program, online at http://streamstats.usgs.gov, accessed on September 27, 2020 \
-Jin, S., Homer, C.G., Yang, L., Danielson, P., Dewitz, J., Li, C., Zhu, Z., Xian, G., Howard, D. (2019). Overall methodology design for the United States National Land Cover Database 2016 products. Remote Sensing, 11(24); https://doi.org/10.3390/rs11242971 \
-Homer, C., Dewitz, J., Jin, S., Xian, G., Costello, C., Danielson, P., Gass, L., Funk, M., Wickham, J., Stehman, S., Auch, R., Riitters, K. (2020). Conterminous United States land cover change patterns 2001–2016 from the 2016 National Land Cover Database. ISPRS Journal of Photogrammetry and Remote Sensing, v. 162, p. 184–199. https://doi.org/10.1016/j.isprsjprs.2020.02.019 \
-Multi-Resolution Land Characteristics Consortium (MRLCC). National Land Cover Database (NLCD) 2016. Retrieved on May 1 2020 from https://www.mrlc.gov/national-land-cover-database-nlcd-2016 \
-Schaider, L. A., Rodgers, K. M., Rudel, R. A. (2017). Review of Organic Wastewater Compound Concentrations and Removal in Onsite Wastewater Treatment Systems. Environ. Sci. Technol., 51, 13, 7304–7317 DOI: 10.1021/acs.est.6b04778 \
-Kolpin, D.W., Furlong, E.T., Meyer, M.T., Thurman, E.M., Zaugg, S.D., Barber, L.B., Buxton, H.T. (2002). Pharmaceuticals, hormones, and other organic wastewater contaminants in U.S. streams, 1999–2000: a national reconnaissance. Environ. Sci. Technol. 36, 1202–1211. https://pubs.acs.org/doi/abs/10.1021/es011055j \
-Kolpin, D.W., Skopec, M., Meyer, M.T., Furlong, E.T., Zaugg, S.D. (2004). Urban contribution of pharmaceuticals and other organic wastewater contaminants to streams during differing flow conditions. Sci. Total Environ. 328, 119–130. https://doi.org/10.1016/j.scitotenv.2004.01.015\
-Veach, A.M., Bernot, M.J. (2011). Temporal variation of pharmaceuticals in an urban and agriculturally influenced stream. Sci. Total Environ. 409 (21), 4553–4563. DOI: 10.1016/j.scitotenv.2011.07.022 \
-Vidal-Dorsch, D.E., Bay, S. M., Maruya, K., Snyder, S. A., Trenholm, R. A., Vanderford, B. J. (2012). Contaminants of emerging concern in municipal wastewater effluents and marine receiving water. Environ. Toxicol. Chem. 31 (12), 2674–2682. DOI: 10.1002/etc.2004 \
-Van Sickle, J., Burch Johnson, C. (2008). Parametric distance weighting of landscape influence on streams. Landscape Ecol 23, 427–438. https://doi.org/10.1007/s10980-008-9200-4 
