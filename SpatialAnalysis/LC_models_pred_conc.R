#Ruta Basijokaite
#----------------
#This code develops models to predict CECs concentrations using developed, pasture, and cultivated crops land cover.
#Separate models are developed for the compounds that occured most frequently.
#----------------
########################
### UPLOADING DATA

#Uploading land cover distributions
Perc_devel_allWS = read.csv(file="LC_developed_allWS_by_dist.csv",sep=",",header=T) 
Perc_crop_allWS = read.csv(file="LC_crop_allWS_by_dist.csv",sep=",",header=T) 
Perc_pasture_allWS = read.csv(file="LC_pasture_allWS_by_dist.csv",sep=",",header=T) 
Perc_agr_allWS = read.csv(file="LC_agr_allWS_by_dist.csv",sep=",",header=T) 

#Uploading CEC concnetrations (all detected compounds)
CEC_aver_conc = read.csv(file="CEC_aver_dataset.csv",sep=",",header=T) 
#Creating duplicate variable, since one variable will have 0s instead of NAs,
#and another will have NAs. Both variables will be used to test if model results 
#depend on how non-detects are treated
CEC_aver_conc_NA = CEC_aver_conc #Variable with NAs
#Site names
Names_u = as.character(unique(CEC_aver_conc[,1]))
#Replacing NAs with 0s in order to calculate averages compound concentrations at a site
CEC_aver_conc[is.na(CEC_aver_conc)] = 0 

########################
### CALCULATING AVERAGES 

Site_conc_aver = Site_conc_aver_NA = matrix(NA, length(Names_u),dim(CEC_aver_conc)[2]-1) #don't need dates
colnames(Site_conc_aver) = colnames(Site_conc_aver_NA) = colnames(CEC_aver_conc)[-2] #removing 'Date'
Site_conc_aver[,1] = Site_conc_aver_NA[,1] = Names_u
for (i in seq(1, length(Names_u))){
  R = which(CEC_aver_conc[,1] == Names_u[i])
  #Calculate CEC average at a site (NAs replaced by 0s)
  Site_conc_aver[i,2:dim(Site_conc_aver)[2]] = as.numeric(lapply(CEC_aver_conc[R,3:dim(CEC_aver_conc)[2]],mean))
  #Calculate CEC averages at a site (keeping NAs)
  Site_conc_aver_NA[i,2:dim(Site_conc_aver_NA)[2]] = as.numeric(lapply(CEC_aver_conc_NA[R,3:dim(CEC_aver_conc_NA)[2]],mean,na.rm=TRUE))
  #After using lapply NA values turns into NaN
}
#Removing NaN and replacing them with NAs
Site_conc_aver_NA[Site_conc_aver_NA =='NaN'] = NA

########################
### CHANGING SITE NAMES 

#Before modeling, change column (site name) names
#Need to replace column names and shorten them, so that they could be matched in land cover 
# and average concentration variables

#Selecting 1st character of a string and capitalize that letter
for (i in seq(1,length(colnames(Perc_devel_allWS)))){
  substr(colnames(Perc_devel_allWS)[i], start = 1, stop = 1) =
    toupper(substr(colnames(Perc_devel_allWS)[i], start = 1, stop = 1))
}
#Check if all site names start with capital letter
colnames(Perc_devel_allWS)
#Find name with the lowest number of characters
N_min = min(nchar((colnames(Perc_devel_allWS))))
#In order to match names, shorten site names to N_min
Site_conc_aver_names_sh = substring(Site_conc_aver[,1],1,N_min)
Perc_devel_allWS_names_sh = substring(colnames(Perc_devel_allWS),1,N_min)
#Manually correct 'Hono' to 'Hone'(misspelled)
Perc_devel_allWS_names_sh[Perc_devel_allWS_names_sh == 'Hono'] = 'Hone'
#Editing Ninemile Creek Site names separately
Perc_devel_allWS_names_sh[Perc_devel_allWS_names_sh == 'Nine'] = c('NineLK', 'NineMAR')
Site_conc_aver_names_sh[Site_conc_aver_names_sh == 'Nine'] = c('NineLK', 'NineMAR')

#Changing column names in Perc_devel_allWS, Perc_crop_allWS, Perc_pasture_allWS, Perc_agr_allWS
colnames(Perc_devel_allWS) = colnames(Perc_crop_allWS) = colnames(Perc_pasture_allWS) =
                             colnames(Perc_agr_allWS) = Perc_devel_allWS_names_sh

#Assign row names to Site_conc_aver, Site_conc_aver_NA as Site_conc_aver_names_sh
rownames(Site_conc_aver) = rownames(Site_conc_aver_NA) = Site_conc_aver_names_sh
#Site_conc_aver and Site_conc_aver_NA have character data type - change that to numeric before using values to model

########################
### PICKING COMPOUNDS TO MODEL

#Most frequently occuring compounds
Ocur_nr = 0
for (i in seq(from=2, to=dim(Site_conc_aver)[2])){
  Ocur_nr[i] = length(which(Site_conc_aver[,i]>0))
} 
#Checking occurance numbers
sort(Ocur_nr, decreasing = TRUE)
#Compounds that occured in 75% of the sites at least once  
colnames(Site_conc_aver)[which(Ocur_nr > 14)]
Freq_CECs = data.frame(colnames(Site_conc_aver)[which(Ocur_nr > 14)],(Ocur_nr[which(Ocur_nr > 14)]/20*100))
Comp_to_mod = as.character(Freq_CECs[,1]) #all compounds names that were in 75% of the sites
#Distance intervals
Dist_seq = seq(from=0.3, to=9.7, by=0.1)

#Changing site order in Site_conc_aver and Site_conc_aver_NA to match site order in land cover variables
for(i in seq(1,length(Perc_devel_allWS_names_sh))) {Ord[i] = which(Perc_devel_allWS_names_sh[i] == rownames(Site_conc_aver))}
#Changing order
Site_conc_aver_ord = Site_conc_aver[Ord,]
Site_conc_aver_NA_ord = Site_conc_aver_NA[Ord,]

########################
### MODELING

#Variables to store r2 values from OLS models
Mod_r2_agr = Mod_NA_r2_agr = matrix(0,length(Dist_seq),length(Comp_to_mod))
Mod_r2_urb = Mod_NA_r2_urb = matrix(0,length(Dist_seq),length(Comp_to_mod))
Mod_r2_agr_urb = Mod_NA_r2_agr_urb = matrix(0,length(Dist_seq),length(Comp_to_mod))
colnames(Mod_r2_agr) = colnames(Mod_NA_r2_agr) = Comp_to_mod
#Testing agricultural (cultivated crops + pasture) land cover for predicting average compound concentrations
for (i in seq(1,length(Comp_to_mod))){
  for (j in seq(1,length(Dist_seq))){
    Col = which(colnames(Site_conc_aver) == Comp_to_mod[i])
    Mod_r2_agr[j,i] = summary(lm(Site_conc_aver_ord[,Col] ~ as.numeric(Perc_agr_allWS[j,])))$r.squared
    Mod_NA_r2_agr[j,i] = summary(lm(Site_conc_aver_NA_ord[,Col] ~ as.numeric(Perc_agr_allWS[j,])))$r.squared
    Mod_r2_urb[j,i] = summary(lm(Site_conc_aver_ord[,Col] ~ as.numeric(Perc_devel_allWS[j,])))$r.squared
    Mod_NA_r2_urb[j,i] = summary(lm(Site_conc_aver_NA_ord[,Col] ~ as.numeric(Perc_devel_allWS[j,])))$r.squared
    Mod_r2_agr_urb[j,i] = summary(lm(Site_conc_aver_ord[,Col] ~ as.numeric(Perc_agr_allWS[j,]) + as.numeric(Perc_devel_allWS[j,])))$r.squared
    Mod_NA_r2_agr_urb[j,i] = summary(lm(Site_conc_aver_NA_ord[,Col] ~ as.numeric(Perc_agr_allWS[j,]) + as.numeric(Perc_devel_allWS[j,])))$r.squared
  }
}

########################
### PLOTTING MODEL RESULTS
#Can land cover predict average compound concentrations?

###PEST - pesticides
#Atrazine
plot(Dist_seq,Mod_r2_agr_urb[,7],pch=19,xlab = 'Distance from sampling location (km)',ylab = 'R2',ylim=c(0,0.63), main='Atrazine average conc. ~ % Agricultural + Developed LC')
plot(Dist_seq,Mod_NA_r2_agr_urb[,7],pch=19,xlab = 'Distance from sampling location (km)',ylab = 'R2',ylim=c(0,0.63), main='Atrazine average conc. ~ % Agricultural + Developed LC')
#Metolachlor
plot(Dist_seq,Mod_r2_agr_urb[,8],pch=19,xlab = 'Distance from sampling location (km)',ylab = 'R2',ylim=c(0,0.63), main='Metolachlor average conc. ~ % Agricultural + Developed LC')
plot(Dist_seq,Mod_NA_r2_agr_urb[,8],pch=19,xlab = 'Distance from sampling location (km)',ylab = 'R2',ylim=c(0,0.63), main='Metolachlor average conc. ~ % Agricultural + Developed LC')
#Atrazine 2H (TP)
plot(Dist_seq,Mod_r2_agr_urb[,9],pch=19,xlab = 'Distance from sampling location (km)',ylab = 'R2',ylim=c(0,0.63), main='Atrazine 2H (TP) average conc. ~ % Agricultural + Developed LC')
plot(Dist_seq,Mod_NA_r2_agr_urb[,9],pch=19,xlab = 'Distance from sampling location (km)',ylab = 'R2',ylim=c(0,0.63), main='Atrazine 2H (TP) average conc. ~ % Agricultural + Developed LC')
#Almost all the contribution comes from agricultural cover
plot(Dist_seq,Mod_r2_agr[,9],pch=19,xlab = 'Distance from sampling location (km)',ylab = 'R2',ylim=c(0,0.63), main='Atrazine 2H (TP) average conc. ~ % Agricultural LC')
#Atrazine des (TP)
#Averages that were calculated for Atrazine des (TP) assuming non-detects were 0, produced better R2 
plot(Dist_seq,Mod_r2_agr_urb[,10],pch=19,xlab = 'Distance from sampling location (km)',ylab = 'R2',ylim=c(0,0.7), main='Atrazine des (TP) average conc. ~ % Agricultural + Developed LC')
plot(Dist_seq,Mod_NA_r2_agr_urb[,10],pch=19,xlab = 'Distance from sampling location (km)',ylab = 'R2',ylim=c(0,0.63), main='Atrazine des (TP) average conc. ~ % Agricultural + Developed LC')
#Metolachlor ESA (TP)
plot(Dist_seq,Mod_r2_agr_urb[,11],pch=19,xlab = 'Distance from sampling location (km)',ylab = 'R2',ylim=c(0,0.63), main='Metolachlor ESA (TP) average conc. ~ % Agricultural + Developed LC')
plot(Dist_seq,Mod_NA_r2_agr_urb[,11],pch=19,xlab = 'Distance from sampling location (km)',ylab = 'R2',ylim=c(0,0.63), main='Metolachlor ESA (TP) average conc. ~ % Agricultural + Developed LC')

###PHAR - pharmaceuticals 
#Lamotrigine
plot(Dist_seq,Mod_r2_agr_urb[,1],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='Lamotrigine average conc. ~ % Developed LC')
plot(Dist_seq,Mod_NA_r2_agr_urb[,1],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='Lamotrigine average conc. ~ % Developed LC')
#Lidocaine
plot(Dist_seq,Mod_r2_agr_urb[,2],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='Lidocaine average conc. ~ % Developed LC')
plot(Dist_seq,Mod_NA_r2_agr_urb[,2],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='Lidocaine average conc. ~ % Developed LC')
#Metformin
#Averages that were calculated for Metformin assuming non-detects were 0, produced better R2 
plot(Dist_seq,Mod_r2_agr_urb[,3],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='Metformin average conc. ~ % Developed LC')
max(Mod_r2_agr_urb[,3]) #= 0.3208141
plot(Dist_seq,Mod_NA_r2_agr_urb[,3],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='Metformin average conc. ~ % Developed LC')
max(Mod_NA_r2_agr_urb[,3]) #= 0.2534893
#Metoprolol
plot(Dist_seq,Mod_r2_agr_urb[,4],pch=19,xlab = 'Distance from sampling location (km)',ylab = 'R2', main='Metoprolol average conc. ~ % Developed LC')
max(Mod_r2_agr_urb[,4]) #= 0.4983866
plot(Dist_seq,Mod_NA_r2_agr_urb[,4],pch=19,xlab = 'Distance from sampling location (km)',ylab = 'R2', main='Metoprolol average conc. ~ % Developed LC')
max(Mod_NA_r2_agr_urb[,4]) #= 0.5048412

###PCHC - personal care products
#5Methyl
#Averages that were calculated for 5Methyl assuming non-detects were 0, produced better R2 
plot(Dist_seq,Mod_r2_agr_urb[,12],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='5Methyl average conc. ~ % Developed LC')
plot(Dist_seq,Mod_NA_r2_agr_urb[,12],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='5Methyl average conc. ~ % Developed LC')
#Benzophenone
#Averages that were calculated for Benzophenone assuming non-detects were NA, produced better R2
plot(Dist_seq,Mod_r2_agr_urb[,13],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='Benzophenone average conc. ~ % Developed LC')
plot(Dist_seq,Mod_NA_r2_agr_urb[,13],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='Benzophenone average conc. ~ % Developed LC')
max(Mod_NA_r2_agr_urb[,13]) #= 0.3235227
#Benzothiazole
#Averages that were calculated for Benzothiazole assuming non-detects were 0, produced better R2 
plot(Dist_seq,Mod_r2_agr_urb[,14],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='Benzothiazole average conc. ~ % Developed LC')
max(Mod_r2_agr_urb[,14]) #= 0.351153
plot(Dist_seq,Mod_NA_r2_agr_urb[,14],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='Benzothiazole average conc. ~ % Developed LC')
#DEET
plot(Dist_seq,Mod_r2_agr_urb[,15],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='DEET average conc. ~ % Developed LC')
plot(Dist_seq,Mod_NA_r2_agr_urb[,15],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='DEET average conc. ~ % Developed LC')
#Sucralose
plot(Dist_seq,Mod_r2_agr_urb[,16],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='Sucralose average conc. ~ % Developed LC')
plot(Dist_seq,Mod_NA_r2_agr_urb[,16],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='Sucralose average conc. ~ % Developed LC')
#Galaxolidone
plot(Dist_seq,Mod_r2_agr_urb[,17],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='Galaxolidone average conc. ~ % Developed LC')
plot(Dist_seq,Mod_NA_r2_agr_urb[,17],xlab = 'Distance from sampling location (km)',ylab = 'R2', main='Galaxolidone average conc. ~ % Developed LC')
max(Mod_r2_agr_urb[,17]) #= 0.2859576
