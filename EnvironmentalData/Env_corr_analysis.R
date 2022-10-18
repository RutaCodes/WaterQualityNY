#Ruta Basijokaite
#----------------
#This code analyzes correlation between detected compound concentrations and environmental conditions
#It also analyzes correlations among sites 
#----------------

#Loading libraries
library(corrplot)
library(Hmisc) 

#Uploading data
CEC_aver = read.csv(file="CEC_aver_dataset.csv",sep=",",header=T) 
Climate_var = read.csv(file="Climate_var.csv",sep=",",header=T)
CEC_use_info = read.csv(file="CEC_use_info_dataset.csv",sep=",",header=T) #info about CEC use, type and detection frequency

#Make sure that column names do not have extra symbols (after upload spaces become .)
colnames(CEC_use_info)
colnames(CEC_use_info)[3] = 'Det.Freq.'
colnames(CEC_use_info)

#-------------
### Preparing CEC data

#Have row names as sites - otherwise in corrplot labels will just be numbers; Remove dates
CEC_aver_comp = as.matrix(CEC_aver[,3:dim(CEC_aver)[2]])
#Keeping only 10 characters in CEC names as column names will be used as labels
colnames(CEC_aver_comp) = substring(colnames(CEC_aver_comp),1,10)

#-------------
### IDENTIFYING COMPOUNDS

#Finding compounds that have higher that 45% detection frequency to do further analysis

#Finding compounds that have detection frequency > 45 
Which_high = which(CEC_use_info$Det.Freq. > 45)
Which_comp_high = CEC_use_info$Name[Which_high]

#Concentrations of those compounds
Conc_high = CEC_aver_comp[,Which_high]

#---------------
### Correlations between environmental variables and compound concentrations 

#gives correlation between matrix itself and other variable (interested only in upper right quadrant)
rcor_vals_clim=rcorr(as.matrix(cbind(Conc_high,Climate_var[,3:14])))
corrplot(rcor_vals_clim$r,tl.col = 'black')
#Therefore, it is more convenient to use 'cor' function
cor_mat_clim=cor(Conc_high,Climate_var[,3:14],use="complete.obs")
corrplot(cor_mat_clim,tl.col = 'black')
#When we use complete.obs, we discard entire row if an NA is present 
#pairwise.complete.obs uses the non-NA values when calculating the correlation between variables 
#In order to preserve as much data as possible to test correlation patterns between variables, 
#pairwise.complete.obs is used. If complete.obs were used, correlation patterns would be determined based
#on samples that have all highly detected compounds, instead of using all available data
cor_mat_clim_pair=cor(Conc_high,Climate_var[,3:14],use="pairwise.complete.obs")
corrplot(cor_mat_clim_pair,tl.col = 'black')
#As correlation coefficient values are low, it would be easier to visually assess values if color scale is adjusted
#Changing color scale values along with color distribution
corrplot(cor_mat_clim_pair,tl.col = 'black',is.corr = FALSE,col.lim = c(-0.5,0.5))

#If ‘use’ is ‘"everything"’, ‘NA’s will propagate conceptually, i.e., a resulting value will be ‘NA’ whenever 
#one of its contributing observations is ‘NA’.

#If ‘use’ is ‘"all.obs"’, then the presence of missing observations will produce an error. If ‘use’ is 
#‘"complete.obs"’ then missing values are handled by casewise deletion (and if there are no complete cases, 
#that gives an error).

#Correlation results in this case highly depend on 'use method chosen', since dataset has a lot of NA values

#-----------------------------------------
#Could add correlation between sites (not samples) that come from 17 highest detected compounds

#First, fix site names
CEC_aver[37:38,1] = CEC_aver[39,1] #Fixing dual names for Ninemile creek Lakeland
CEC_aver[41:42,1] = CEC_aver[43,1] #Fixing dual names for Ninemile creek Marieta
CEC_aver[48:49,1] = CEC_aver[50,1] #Fixing misspelling in Oatka
CEC_aver[25:26,1] = CEC_aver[27,1] #Fixing misspelling in Ganargua

#Update CEC_aver variable -  save with corrected site names
write.table(CEC_aver,file="CEC_aver_dataset.csv",sep=",",row.names=F,col.names=T,append=T)

WS_names = as.character(unique(CEC_aver$Site))
Site_aver = matrix(NA, length(WS_names),dim(CEC_aver)[2]-2) #no need for dates
#Calculating average concentration of each compound
for (i in seq(from=1,to=length(WS_names))){
  rows = which(CEC_aver$Site == WS_names[i])
  #When lapply is used, it returns list format. To assign results from lapply
  #to a variable, first list needs to be converted to numeric data type using as.numeric
  Site_aver[i,] = as.numeric(lapply(CEC_aver[rows,-c(1,2)],mean,na.rm=T))
}
row.names(Site_aver) = substring(WS_names,1,5) #shorten site names as they will be used as labels for correlation plot

#Correlation between sites - using average compound concentrations - see how this changes, if NAs are counted as 0s
CEC_sites_trasf=t(Site_aver)
rcor_CEC_sites_spear=rcorr(as.matrix(CEC_sites_trasf),type = c("spearman"))
corrplot(rcor_CEC_sites_spear$r,tl.col = 'black')

#Replacing NAs with 0s and comparing averages
Site_aver_DF = Site_aver
#After converting from list to numeric, NA values become NaN
Site_aver_DF[Site_aver_DF == "NaN"] = 0

CEC_sites_DF_trasf=t(Site_aver_DF)
rcor_CEC_sites_DF_spear=rcorr(as.matrix(CEC_sites_DF_trasf),type = c("spearman"))
corrplot(rcor_CEC_sites_DF_spear$r,tl.col = 'black')
#As correlation values among sites are similar, adjust color scale for better visualization
corrplot(rcor_CEC_sites_DF_spear$r,tl.col = 'black',is.corr = TRUE,col.lim = c(0,1))

#Big difference in results. Results from correlation analysis using 0s instead of NA values are more realistic, as 
#NA indicates that compound was not present in the sample (or at concentrations that were not detectable) 
