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
###Preparing CEC data

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
###Correlations between environmental variables and compound concentrations 

#gives correlation between matrix itself and other variable (interested only in upper right quadrant)
rcor_vals_clim=rcorr(as.matrix(cbind(Conc_high,Climate_var[,3:14])))
corrplot(rcor_vals_clim$r,tl.col = 'black')
#Therefore, it is more convenient to use 'cor' function
cor_mat_clim=cor(Conc_high,Climate_var[,3:14],use="complete.obs")
corrplot(cor_mat_clim,tl.col = 'black')

#-----------------------------------------
#Could add correlation between sites (not samples) that come from 17 highest detected compounds

#First, fix site names
CEC_aver[37:38,1] = CEC_aver[39,1] #Fixing dual names for Ninemile creek Lakeland
CEC_aver[41:42,1] = CEC_aver[43,1] #Fixing dual names for Ninemile creek Marieta
CEC_aver[48:49,1] = CEC_aver[50,1] #Fixing misspelling in Oatka
CEC_aver[25:26,1] = CEC_aver[27,1] #Fixing misspelling in Ganargua

WS_names = as.character(unique(CEC_aver$Site))
Site_aver = matrix(NA, length(WS_names),dim(CEC_aver)[2]-2) # no need for dates
#colnames(Site_aver) = colnames(CEC_aver)
#Site_aver[,1] = WS_names #make sure this does not change matrix type to character
#calculating average concnetration of each compound
for (i in seq(from=1,to=length(WS_names))){
  rows = which(CEC_aver$Site == WS_names[i])
  Site_aver[i,] = as.numeric(lapply(CEC_aver[rows,-c(1,2)],mean,na.rm=T))
}
row.names(Site_aver) = substring(WS_names,1,5) #shorten site names as they will be used as labels for correlation plot

#Correlation between sites - using average compound concentrations - see how this changes, if NAs are counted as 0s
CEC_sites_trasf=t(Site_aver)
rcor_CEC_sites_spear=rcorr(as.matrix(CEC_sites_trasf),type = c("spearman"))
corrplot(rcor_CEC_sites_spear$r,tl.col = 'black')

#Try adding 0s and compare results
Site_aver_DF = Site_aver
Site_aver_DF[Site_aver_DF=="NaN"] = 0

CEC_sites_DF_trasf=t(Site_aver_DF)
rcor_CEC_sites_DF_spear=rcorr(as.matrix(CEC_sites_DF_trasf),type = c("spearman"))
corrplot(rcor_CEC_sites_DF_spear$r,tl.col = 'black')

#Big difference

