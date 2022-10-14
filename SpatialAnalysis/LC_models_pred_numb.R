#Ruta Basijokaite
#----------------
#This code develops models to predict number of detected CECs by type using developed, pasture, and cultivated crops land cover.
#Two methods are compared: 1) using whole watershed scale land cover percentages and 2) developing models using land cover
#percentages at varying distances from sampling location a.k.a. distance weighted models
#----------------

### LOADING DATA

#Developed, agricultural, pasture, and cultivated crops land cover percentages at different distances from sampling locations
Perc_devel_allWS = read.csv(file="LC_developed_allWS_by_dist.csv",sep=",",header=T) 
Perc_crop_allWS = read.csv(file="LC_crop_allWS_by_dist.csv",sep=",",header=T) 
Perc_pasture_allWS = read.csv(file="LC_pasture_allWS_by_dist.csv",sep=",",header=T) 
Perc_agr_allWS = read.csv(file="LC_agr_allWS_by_dist.csv",sep=",",header=T) 

#Loading number of CECs found by type
Total_CEC_nrs = read.csv(file="Total_CEC_Nrs.csv",sep=",",header=T) 

#Preparing variables
Tot_CEC_nr_count = cbind(Total_CEC_nrs[,2],Total_CEC_nrs[,2]+Total_CEC_nrs[,3],
                         Total_CEC_nrs[,4],Total_CEC_nrs[,4]+Total_CEC_nrs[,5],
                         Total_CEC_nrs[,6],Total_CEC_nrs[,6]+Total_CEC_nrs[,7],
                         Total_CEC_nrs[,8])
colnames(Tot_CEC_nr_count) = c('PHARonly','PHARandTP','PESTonly','PESTandTP',
                               'PCHConly','PCHCand TP','TOTAL')

CEC_nr_agr_r2 = CEC_nr_devel_r2 = CEC_nr_crop_r2 = CEC_nr_pasture_r2 = CEC_nr_allLC_r2 = 
  matrix(0,length(seq(from=0.3, to=9.7, by=0.1)),dim(Tot_CEC_nr_count)[2]) 
colnames(CEC_nr_agr_r2) = colnames(CEC_nr_devel_r2) = colnames(CEC_nr_crop_r2) = colnames(CEC_nr_pasture_r2) = colnames(CEC_nr_allLC_r2) = 
  c('PHARonly','PHARandTP','PESTonly','PESTandTP','PCHConly','PCHCand TP','TOTAL')
rownames(CEC_nr_agr_r2) = rownames(CEC_nr_devel_r2) = rownames(CEC_nr_crop_r2) = rownames(CEC_nr_pasture_r2) = rownames(CEC_nr_allLC_r2) =
  seq(from=0.3, to=9.7, by=0.1)

#--------------------------------------------------------------
#### DEVELOPING DISTANCE WEIGHTED MODELS TO PREDICT NUMBER OF CECs

#Pulling R2 value from OLS model created to predict total number of CECs using land cover percentages as explanatory variable
for (i in seq(1:7)){ #Predicting each type of CEC separately
  for (j in seq(1:95)){
    #Predicting number of CECs using agricultural land cover (agricultural = pasture + cultivated crops)
    CEC_nr_agr_r2[j,i]=summary(lm(Tot_CEC_nr_count[,i] ~ as.numeric(Perc_agr_allWS[j,])))$r.squared
    #Predicting number of CECs using developed land cover
    CEC_nr_devel_r2[j,i]=summary(lm(Tot_CEC_nr_count[,i] ~ as.numeric(Perc_devel_allWS[j,])))$r.squared
    #Predicting number of CECs using cultivated crops land cover
    CEC_nr_crop_r2[j,i]=summary(lm(Tot_CEC_nr_count[,i] ~ as.numeric(Perc_crop_allWS[j,])))$r.squared
    #Predicting number of CECs using pasture land cover
    CEC_nr_pasture_r2[j,i]=summary(lm(Tot_CEC_nr_count[,i] ~ as.numeric(Perc_pasture_allWS[j,])))$r.squared
    #Predicting number of CECs using developed and agricultural land cover
    CEC_nr_allLC_r2[j,i]=summary(lm(Tot_CEC_nr_count[,i] ~ as.numeric(Perc_devel_allWS[j,]) + as.numeric(Perc_crop_allWS[j,]) + as.numeric(Perc_pasture_allWS[j,])))$r.squared
  }
}


#### ANALYZING MODEL RESULTS FOR EACH CEC TYPE SEPARATELY #####
#--------------------------------------------------------------
# PREDICTING NUMBER OF PESTICIDES 

Dist_seq=seq(from=0.3, to=9.7, by=0.1)
#Predicting # of PEST - R2 changes with distance from sampling location 
plot(Dist_seq,CEC_nr_allLC_r2[,3],pch=19,xlab = 'Distance from sampling location',ylab = 'R2',ylim=c(0,0.76), main='# PEST ~ % Crop., Pasture, Developed land cover')
#Predicting # of PEST + TP
plot(Dist_seq,CEC_nr_allLC_r2[,4],pch=19, xlab = 'Distance from sampling location',ylab = 'R2', ylim=c(0,0.76),main='# PEST + TP ~ % Crop., Pasture, Developed land cover')
#Model summary
summary(lm(Tot_CEC_nr_count[,4] ~ as.numeric(Perc_devel_allWS[1,]) + as.numeric(Perc_crop_allWS[1,]) + as.numeric(Perc_pasture_allWS[1,])))
#Predicted vs observed number of PEST at a distance with highest R2 - Max r2 value with PEST is closest to sampling location
PEST_mod_multi=lm(Tot_CEC_nr_count[,3] ~ as.numeric(Perc_devel_allWS[1,]) + as.numeric(Perc_crop_allWS[1,]) + as.numeric(Perc_pasture_allWS[1,]))
plot(Tot_CEC_nr_count[,3],PEST_mod_multi$fitted.values,xlab='Observed # PEST',ylab='Predicted # PEST')

#--------------------------------------------------------------
### PREDICTING NUMBER OF PHARMACEUTICALS

#Predicting # of PHAR - R2 changes with distance from sampling location
plot(Dist_seq,CEC_nr_allLC_r2[,1],pch=19,xlab = 'Distance from sampling location',ylab = 'R2',ylim=c(0,0.76), main='# PHAR ~ % Crop., Pasture, Developed land cover')
#Predicting # of PHAR + TP 
plot(Dist_seq,CEC_nr_allLC_r2[,2],pch=19,xlab = 'Distance from sampling location',ylab = 'R2',ylim=c(0,0.76), main='# PHAR + TP ~ % Crop., Pasture, Developed land cover')
#At what distance R2 is the highest?
nr_max=which(CEC_nr_allLC_r2[,1]==max(CEC_nr_allLC_r2[,1])) 
Dist_seq[nr_max] #distance with highest R2
#Model summary
summary(lm(Tot_CEC_nr_count[,1] ~ as.numeric(Perc_devel_allWS[nr_max,]) + as.numeric(Perc_crop_allWS[1,]) + as.numeric(Perc_pasture_allWS[nr_max,])))
#Predicted vs observed number of PHAR at a distance with highest R2
PHAR_mod_multi=lm(Tot_CEC_nr_count[,1] ~ as.numeric(Perc_devel_allWS[nr_max,]) + as.numeric(Perc_crop_allWS[nr_max,]) + as.numeric(Perc_pasture_allWS[nr_max,]))
plot(Tot_CEC_nr_count[,1],PHAR_mod_multi$fitted.values,xlab='Observed # PHAR',ylab='Predicted # PHAR')

#--------------------------------------------------------------
### PREDICTING NUMBER OF PERSONAL CARE PRODUCTS

#Predicting # of PCHC - R2 changes with distance from sampling location
plot(Dist_seq,CEC_nr_allLC_r2[,5],pch=19,xlab = 'Distance from sampling location',ylab = 'R2',ylim=c(0,0.76), main='# PCHC ~ % Crop., Pasture, Developed land cover')
#Predicting # of PHAR + TP
plot(Dist_seq,CEC_nr_allLC_r2[,6],pch=19,xlab = 'Distance from sampling location',ylab = 'R2',ylim=c(0,0.76), main='# PCHC + TP ~ % Crop., Pasture, Developed land cover')
#At what distance R2 is the highest?
nr_max2=which(CEC_nr_allLC_r2[,5]==max(CEC_nr_allLC_r2[,5])); Dist_seq[nr_max2]
#Model summary
summary(lm(Tot_CEC_nr_count[,5] ~ as.numeric(Perc_devel_allWS[nr_max2,]) + as.numeric(Perc_crop_allWS[nr_max2,]) + as.numeric(Perc_pasture_allWS[nr_max2,])))
#Predicted vs observed number of PCHC at a distance with highest R2
PCHC_mod_multi=lm(Tot_CEC_nr_count[,5] ~ as.numeric(Perc_devel_allWS[nr_max2,]) + as.numeric(Perc_crop_allWS[nr_max2,]) + as.numeric(Perc_pasture_allWS[nr_max2,]))
plot(Tot_CEC_nr_count[,5],PCHC_mod_multi$fitted.values,xlab='Observed # PCHC',ylab='Predicted # PCHC')

#--------------------------------------------------------------
### PREDICTING TOTAL NUMBER OF COMPOUNDS

plot(Dist_seq,CEC_nr_allLC_r2[,7],pch=19,xlab = 'Distance from sampling location',ylab = 'R2', main='# TOT ~ % Crop., Pasture, Developed land cover')
#At what distance R2 is the highest?
nr_max3=which(CEC_nr_allLC_r2[,7]==max(CEC_nr_allLC_r2[,7])); Dist_seq[nr_max3]
#Model summary
summary(lm(Tot_CEC_nr_count[,7] ~ as.numeric(Perc_devel_allWS[nr_max3,]) + as.numeric(Perc_crop_allWS[nr_max3,]) + as.numeric(Perc_pasture_allWS[nr_max3,])))
#Predicted vs observed number of CECs at a distance with highest R2
TOT_mod_multi=lm(Tot_CEC_nr_count[,7] ~ as.numeric(Perc_devel_allWS[nr_max3,]) + as.numeric(Perc_crop_allWS[nr_max3,]) + as.numeric(Perc_pasture_allWS[nr_max3,]))
plot(Tot_CEC_nr_count[,7],TOT_mod_multi$fitted.values,xlab='Observed # TOT',ylab='Predicted # TOT')
#It is closest to PHAR model, since in most of the samples PHAR compounds dominate

#--------------------------------------------------------------
### INDIVIDUAL LAND COVER CONTRIBUTIONS

#correlation with inditidual land cover layers
plot(Dist_seq,CEC_nr_agr_r2[,3],xlab = 'Distance from sampling location',ylab = 'R2', main='# PEST ~ % Agr. land cover')
plot(Dist_seq,CEC_nr_agr_r2[,4],xlab = 'Distance from sampling location',ylab = 'R2', main='# PEST + TP ~ % Agr. land cover')
plot(Dist_seq,rowMeans(Perc_agr_allWS[1:95,]),xlab = 'Distance from sampling location',ylab = 'Average % Agr. land cover',main='Aver. Agr. land cover distribution')

plot(Dist_seq,CEC_nr_crop_r2[,3],xlab = 'Distance from sampling location',ylab = 'R2', main='# PEST ~ % Crop land cover')
plot(Dist_seq,CEC_nr_crop_r2[,4],xlab = 'Distance from sampling location',ylab = 'R2', main='# PEST + TP ~ % Crop land cover')
plot(Dist_seq,rowMeans(Perc_crop_allWS[1:95,]),xlab = 'Distance from sampling location',ylab = 'Average % Crop. land cover')

plot(Dist_seq,CEC_nr_pasture_r2[,1],xlab = 'Distance from sampling location',ylab = 'R2', main='# PHAR ~ % Pasture land cover')
plot(Dist_seq,CEC_nr_pasture_r2[,2],xlab = 'Distance from sampling location',ylab = 'R2', main='# PHAR + TP ~ % Pasture land cover')
plot(Dist_seq,rowMeans(Perc_pasture_allWS[1:95,]),xlab = 'Distance from sampling location',ylab = 'Average % Pasture land cover')

plot(Dist_seq,CEC_nr_devel_r2[,5],xlab = 'Distance from sampling location',ylab = 'R2', main='# PCHC ~ % Developed land cover')
plot(Dist_seq,CEC_nr_devel_r2[,6],xlab = 'Distance from sampling location',ylab = 'R2', main='# PCHC + TP ~ % Developed land cover')


#How does land cover distribution looks like as a function of distance from sampling location?
plot(Dist_seq,rowMeans(Perc_devel_allWS[1:95,]),xlab = 'Distance from sampling location',ylab = 'Average % Developed land cover')

#Cumulative % of agricultural and developed land cover 
Cum_LC=rowMeans(Perc_agr_allWS[1:95,]+Perc_devel_allWS[1:95,])
plot(Dist_seq,Cum_LC,xlab = 'Distance from sampling location',ylab = 'Aver % Developed + Agr land cover')

plot(Dist_seq,rowMeans(Perc_crop_allWS[1:95,]),xlab = 'Distance from sampling location',ylab = 'Aver % Crop land cover') #averaged cultivated crop land cover distribution with distance
plot(Dist_seq,rowMeans(Perc_pasture_allWS[1:95,]),xlab = 'Distance from sampling location',ylab = 'Aver % Pasture land cover') #averaged pasture land cover distribution with distance

#Sample watershed 
plot(Dist_seq,Perc_crop_allWS[1:95,4],xlab = 'Distance from sampling location',ylab = '% Crop land cover',main='Chenango') #col 4 - sample watershed


###########################################################
### PREDICTING NUMBER OF CECs USING WHOLE WATERSHED LAND COVER 

#Loading whole watershed land cover percentages
LC_CEC_nrs = read.csv(file="NR_vs_landcover.csv",sep=",",header=T)[,-8] #removing PCHC TP count from dataset, as there is only one compound present in all the samples 

#Pulling R2 value from OLS model created to predict number of CECs using whole watershed land cover %
CEC_LC_r2 = CEC_LC_r2_agr = CEC_LC_r2_urb = matrix(0,8,1)
for (i in seq(from=4, to=11)){
  CEC_LC_r2[i-3]=summary(lm(LC_CEC_nrs[,i] ~ LC_CEC_nrs[,2]+LC_CEC_nrs[,3]))$r.squared
  CEC_LC_r2_agr[i-3]=summary(lm(LC_CEC_nrs[,i] ~ LC_CEC_nrs[,2]))$r.squared
  CEC_LC_r2_urb[i-3]=summary(lm(LC_CEC_nrs[,i] ~ LC_CEC_nrs[,3]))$r.squared
}

#What is the highest R2 value?
max(CEC_LC_r2)
max(CEC_LC_r2_agr)
max(CEC_LC_r2_urb)

#--------------
### CONCLUSION:
#Using distance weighting model produces better model accuracy compared to models developed using land cover 
#percentages from the whole watershed
