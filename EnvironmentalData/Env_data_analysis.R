#Ruta Basijokaite
#----------------
#This code:
#1) Corrects missing date values in streamflow dataset
#2) Converts streamflow values from cfs to mm/day
#3) Performs baseflow separation
#4) Finds 4,7,14 day streamflow averages and cumulative precipitation values
#5) Saves this variable as csv 
#----------------

#### UPLOADING STREAMFLOW DATA

#Read in all daily streamflow values from a folder
library(readr)
files_ext <- list.files(path = "~/Desktop/Stream_Flow_data/Q_extended", pattern = "*.csv", full.names = T)
for (i in 1:length(files_ext)) assign(files_ext[i], read.csv(files_ext[i]))

#----------------
#### CLEANING DATASET

#Filling missing dates 
nr_sites=length(files_ext)
site_numb_list=0
#Creating complete datasets without missing dates
date_file= as.POSIXlt(seq(as.Date("2000/01/01"), as.Date("2019/12/31"), by="day"))
date_file_col=data.frame(date_file$year+1900, date_file$mon+1, date_file$mday)
names(date_file_col)=c("year","month","day")
Q_ext_WS=matrix(0,dim(date_file_col)[1],(nr_sites+3)) #adding year, month, day
Q_ext_WS[,1]=date_file_col[,1] #year
Q_ext_WS[,2]=date_file_col[,2] #month
Q_ext_WS[,3]=date_file_col[,3] #day
#library(dplyr)
for (i in 1:nr_sites){
  data_ext=read.csv(files_ext[i])
  dataF_ext=data.frame(data_ext)
  names(dataF_ext)=c("site","year","month","day","Q")
  site_nr=dataF_ext[1,1]
  site_numb_list[i]=site_nr
  data_cor=left_join(date_file_col,dataF_ext[,2:5], on=c("year","month","day"))
  Q_ext_WS[,i+3]=data_cor[,4]
}
Q_ex_WS_DF=data.frame(Q_ext_WS)
names(Q_ex_WS_DF)=c("year","month","day",site_numb_list)
Q_ex_WS_DF_names=Q_ex_WS_DF
names(Q_ex_WS_DF_names)=c("year","month","day","Oriskany","Sauquoit","Chenango","Otselic","Conesus","Honeoye","Oatka","Sixmile","Salmon","Ganargua","Flint","Skaneateles","Harbor","Ninemile MAR","Ninemile LK","Oneida","Cowaselon","Chittenango","Scriba")

#Creating new variable that will be have streamflow in mm/day
Q_val_mm=Q_ex_WS_DF_names[,-17] #excluding Ninemile MAR since there is no record

#----------------
### CONVERTING UNITS AND PERFORMING BASEFLOW SEPARATION

source("BaseflowSeparation.R")
nr_sites=dim(Q_val_mm)[2]-3 #3 columns for year, month, day
data_length=dim(Q_val_mm)[1]
#Manually entering watershed drainage area values
WS_size=c(373.74,155.14,677.42,380.63,187.58,502.86,530.64,125.80,229.61,296.83,262.79,193.38,31.20,
          298.13,295.14,107.29,175.37,106.43)
Baseflow_c = Runoff_c = data.frame(matrix(NA,dim(Q_val_mm)[1],dim(Q_val_mm)[2])) 
colnames(Baseflow_c) = colnames(Runoff_c) = colnames(Q_val_mm)
Baseflow_c[,1:3] = Runoff_c[,1:3] = Q_val_mm[,1:3]
for (i in 1:nr_sites){
    Q_val=Q_val_mm[,(i+3)]
    Q_val_mm[,i+3]=(Q_val*0.0283168/(WS_size[i]*10^6))*1000*3600*24
    
    #Finding NA values - baseflow separation function cannot deal with NAs
    NonNA_ind = which(!is.na(Q_val_mm[,i+3]))
    Basef_run = data.frame(matrix(NA, length(Q_val_mm[,i]),2))
    #Baseflow separation
    Basef_run[NonNA_ind,] = BaseflowSeparation(Q_val_mm[NonNA_ind,i+3], filter_parameter = 0.925, passes = 3)  
    Baseflow_c[,i+3] = Basef_run[,1]
    Runoff_c[,i+3] = Basef_run[,2]
}

#----------------
### FINDING 4,7,14 DAY AVERAGES 

#Uploading precipitation data and sample info
precip_ed = read.csv(file="USDA_precip_ed.csv",sep=",",header=T) #in mm/day
CEC_aver = read.csv(file="CEC_aver_dataset.csv",sep=",",header=T) 
samp_date = as.POSIXlt(CEC_aver$Date, format= "%m/%d/%y")

#Removing watersheds that do not have Q values
WS_col = c(1:16,18:21,23)
precip_samp_WS = precip_ed[,WS_col]

#Finding rows in Q and P datasets that store climate values for sampling days
#However, precipitation dataset does not have leap days, therefore row numbers for the same date will 
#vary between these two datasets
N=length(samp_date)
row = rowP = 0
for (i in seq(1:N)){
  row[i]=which(samp_date[i]$year+1900==Q_val_mm$year & samp_date[i]$mon+1==Q_val_mm$month & samp_date[i]$mday==Q_val_mm$day)
  rowP[i]=which(samp_date[i]$year+1900==precip_samp_WS$year & samp_date[i]$mon+1==precip_samp_WS$month & samp_date[i]$mday==precip_samp_WS$day)
}

#Shorten site names, so that sites could be matched by names. Full names vary between datasets
Aver_n = substring(CEC_aver[,1],1,5)
Precip_nm = substring(colnames(precip_ed),1,5) 
Q_nm = substring(colnames(Q_val_mm),1,5) 

#Butternut creek and Ninemile don't have Q data to calculate averages
View(CEC_aver)
#Butternut creek - rows 1:4
#Ninemile at Marieta - rows 41:47

#Creating data frame with climate variable averages
Par = 13 #number of columns in climate matrix - first two columns: site, date, others are averages 
Clim_var = data.frame(matrix(NA, dim(CEC_aver)[1], Par))
colnames(Clim_var) = c('Site','Date','P','P_4d_cum','P_7d_cum','P_14d_cum','Q','Q_4d_aver','Q_7d_aver','Q_14d_aver','Baseflow',
                       'Runoff','RR_7d')
Clim_var[,1:2] = CEC_aver[,1:2] #First two columns - site and date

Q_7d_cum = 0
for (i in seq(1,90)[-c(1:4)]){ #Butternut creek is removed
  if (i > 40 & i < 47){ #for Ninemile MAR use Q from Ninemile Lakeland, as that site is within the watershed
    Ind_P = 17 
  }else{
    Ind_P = which(Aver_n[i] == Precip_nm)
  }
  Ind_Q = which(Aver_n[i] == Q_nm)
  #Precipitation and streamflow values during the day of sampling
  Clim_var$P[i] = precip_ed[rowP[i],Ind_P]
  Clim_var$Q[i] = Q_val_mm[row[i],Ind_Q]
  #Baseflow and Runoff values during the day of sampling
  Clim_var$Baseflow[i] = Baseflow_c[row[i],Ind_Q]
  Clim_var$Runoff[i] = Runoff_c[row[i],Ind_Q]
  #4, 7, 14 day precipitation and streamflow values 
  Clim_var$P_4d_cum[i] = sum(precip_ed[(rowP[i] - 3) : rowP[i],Ind_P])
  Clim_var$Q_4d_aver[i] = mean(Q_val_mm[(row[i] - 3) : row[i],Ind_Q])
  Clim_var$P_7d_cum[i] = sum(precip_ed[(rowP[i] - 6) : rowP[i],Ind_P])
  Clim_var$Q_7d_aver[i] = mean(Q_val_mm[(row[i] - 6) : row[i],Ind_Q])
  Q_7d_cum  = sum(Q_val_mm[(row[i] - 6) : row[i],Ind_Q])
  Clim_var$P_14d_cum[i] = sum(precip_ed[(rowP[i] - 13) : rowP[i],Ind_P])
  Clim_var$Q_14d_aver[i] = mean(Q_val_mm[(row[i] - 13) : row[i],Ind_Q])
  Clim_var$RR_7d[i] = Q_7d_cum/Clim_var$P_7d_cum[i]
}

#Calculating P averages for Butternut Creek, Q values are not available
for (i in seq(1,4)){
  Ind_P = which(Aver_n[i] == Precip_nm)
  Clim_var$P[i] = precip_ed[rowP[i],Ind_P]
  Clim_var$P_4d_cum[i] = sum(precip_ed[(rowP[i] - 3) : rowP[i],Ind_P])
  Clim_var$P_7d_cum[i] = sum(precip_ed[(rowP[i] - 6) : rowP[i],Ind_P])
  Clim_var$P_14d_cum[i] = sum(precip_ed[(rowP[i] - 13) : rowP[i],Ind_P])
}

#Saving calculated averages as csv file
write.table(Clim_var,file="CEC_clim_aver.csv",sep=",",row.names=F,col.names=T,append=T)


#From stackoverflow about converting from list to numeric:
#DF[4, ] gives a one row data frame which is a list, while matrix is an atomic vector which can hold only 
#one data type. You need to unlist the data frame row and convert it to an atomic vector before assigning it to the matrix
#A[,1] = unlist(DF[4,])

#Converting data frame to numeric matrix, so that matrix could be used with 'cor' function as it requires numeric matrix
Clim_var_num = as.matrix(sapply(Clim_var[,3:13],as.numeric))

#Checking correlation results, as 4,7,14 day Q averages are slightly different than in saved Climate_var.csv
cor_mat_clim_pair_c = cor(Conc_high,Clim_var_num,use="pairwise.complete.obs")
corrplot(cor_mat_clim_pair_c,tl.col = 'black')
corrplot(cor_mat_clim_pair_c,tl.col = 'black',is.corr = FALSE,col.lim = c(-0.5,0.5))

