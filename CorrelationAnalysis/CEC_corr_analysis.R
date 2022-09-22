#Ruta Basijokaite
#-----------------
#This code explores correlation pattenrs among samples as well as frequently found compounds 
#-----------------

#Loading libraries
library(corrplot) #corrplot plot function
library(Hmisc) #import rcorr function

#Uploading data 
CEC_aver = read.csv(file="CEC_aver_dataset.csv",sep=",",header=T) #array with CEC concentrations, site names and sample dates
CEC_use_info = read.csv(file="CEC_use_info_dataset.csv",sep=",",header=T) #info about CEC use, type and detection frequency

#Make sure that column names do not have extra symbols (after upload spaces become .)
colnames(CEC_use_info)
colnames(CEC_use_info)[3] = 'Det.Freq.'

#Have row names as sites - otherwise in corrplot labels will just be numbers; Remove dates
CEC_aver_comp = as.matrix(CEC_aver[,3:dim(CEC_aver)[2]])
#Keeping only 7 characters in CEC and site names as row and column names will be used as labels
row.names(CEC_aver_comp) = substring(as.character(CEC_aver[,1]),1,7)
colnames(CEC_aver_comp) = substring(colnames(CEC_aver_comp),1,7)

#-----------------
#### Correlations among samples

#Correlation between samples
CEC_aver_conc_trasf=t(CEC_aver_comp)
rcor_CEC_samp_spear=rcorr(as.matrix(CEC_aver_conc_trasf),type = c("spearman"))
corrplot(rcor_CEC_samp_spear$r,tl.col = 'black') 

#-----------------
#### Correlations among highest occuring compounds

#Since some compounds had very low detection frequency, correlation coefficient could no be calculated for all compounds
#Instead, compounds with higher than 45% detection frequency were selected to study co-occurance patterns 

#Finding compounds that have detection frequency > 45 
Which_high = which(CEC_use_info$Det.Freq. > 45)
Which_comp_high = CEC_use_info$Name[Which_high]
#Most frequently detected compound info
CEC_type_high = CEC_use_info[Which_high,]

#Concentrations of those compounds
Conc_high = CEC_aver_comp[,Which_high]

#Correlation between frequently occuring compounds
rcor_CEC_high_df_spear=rcorr(as.matrix(Conc_high),type = c("spearman")) 
corrplot(rcor_CEC_high_df_spear$r,tl.col = 'black')

#If using 'hclust', corrplot() can draw rectangles around the plot of correlation matrix based on the 
#results of hierarchical clustering
corrplot(rcor_CEC_high_df_spear$r, order = 'hclust', addrect = 2,tl.col = 'black')

#Arrange by detection frequency within each CEC type
CEC_type_high_ord = CEC_type_high %>% group_by(Type) %>% arrange(Type,desc(Det.Freq.))

