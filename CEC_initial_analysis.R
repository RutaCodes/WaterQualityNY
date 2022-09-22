#Ruta Basijokaite
#----------------
#This code analyzes CEC concentrations from 90 samples taken over two separate sampling periods from
#22 separate sampling sites
#---------------
library(dplyr)

#Uploading data 
CEC_conc_not_aver = read.csv(file="CEC_info.csv",sep=",",header=T) 
#This dataset contains sampling site names, sample dates, and duplicate compound concentrations.
#Each sample was analyzed twice producing two sample estimates. 

#Getting names of detected compounds
Name_index=seq(from=3, to=length(colnames(CEC_conc_not_aver)), by=2)
CEC_detected_names=colnames(CEC_conc_not_aver)[Name_index]
#Removing . from names that got added after data that was imported from csv (names had spaces)
CEC_detected_names=gsub('.','',CEC_detected_names,fixed=TRUE)
#Removing X at the start of some CEC names
CEC_detected_names=gsub('X','',CEC_detected_names)

#----------------
#### DETERMINING COMPOUND CONCENTRATIONS

#Overall compound concentration is determined by averaging concentrations of two sample 
#estimates and are reported in ng/L. 
CEC_detected_aver=matrix(0,nrow=dim(CEC_conc_not_aver)[1],ncol=length(CEC_detected_names))
colnames(CEC_detected_aver)=CEC_detected_names
for (i in seq(1:length(CEC_detected_names))){
  for (j in seq(1:dim(CEC_conc_not_aver)[1])){ #90 samples
    CEC_detected_aver[j,i]=(CEC_conc_not_aver[j,Name_index[i]]+CEC_conc_not_aver[j,(Name_index[i]+1)])/2
  }
}

paste("Total number of detected compounds:",length(CEC_detected_names))

#Creating data frame combining sample site, sample date, and concentrations
CEC_aver=data.frame(CEC_conc_not_aver[,1:2],CEC_detected_aver)
colnames(CEC_aver)=c("Site","Date",CEC_detected_names)

#----------------
#### SUMMARIZING COMPOUNDS BY TYPE

#Uploading CEC use data
CEC_use = read.csv(file="CEC_use.csv",sep=",",header=T) 

#CEC names should match colnames in CEC_detected_names
CEC_use_names = as.character(CEC_use$Name)
ch = c("-",",","[()]"," ") 
for (i in seq(1,length(ch))){
  CEC_use_names = gsub(ch[i],'',CEC_use_names)
}

#Making sure that all names match in two variables 
length(CEC_detected_names %in% CEC_use_names) == length(CEC_detected_names) #if all names match, this should be TRUE

#Replacing names
CEC_use[,1] = CEC_use_names

#How many CECs of each type were detected?
paste("Number of Personal Care Products (PCHC) detected:", length(which(CEC_use[,2]=='PCHC')))
paste("Number of Personal Care Product transformation products (PCHC TP) detected:", length(which(CEC_use[,2]=='PCHC TP')))
paste("Number of Pesticides (PEST) detected:", length(which(CEC_use[,2]=='PEST')))
paste("Number of Pesticide transformation products (PEST) detected:", length(which(CEC_use[,2]=='PEST TP')))
paste("Number of Pharmaceuticals (PHAR) detected:", length(which(CEC_use[,2]=='PHAR')))
paste("Number of Pharmaceutical transformation products (PHAR TP) detected:", length(which(CEC_use[,2]=='PHAR TP')))

#----------------
#### ANALYZING COMPOUND CONCENTRATIONS

#Calculating min, max, average, and cv for each compound 
#Compounds that were found in only one sample will return NA value for some of these metrics
Min_v = apply(CEC_detected_aver,MARGIN=2,FUN=min,na.rm=T)
Max_v = apply(CEC_detected_aver,MARGIN=2,FUN=max,na.rm=T)
Median_v = apply(CEC_detected_aver,MARGIN=2,FUN=median,na.rm=T)
Mean_v = apply(CEC_detected_aver,MARGIN=2,FUN=mean,na.rm=T)
Stdev_v = apply(CEC_detected_aver,MARGIN=2,FUN=sd,na.rm=T)
cv_v = Stdev_v/Mean_v*100

#Calculate detection frequency for each compound
Det_freq=0 
for (i in seq(from=1,to=dim(CEC_detected_aver)[2])){
  #detection frequency=number of times compound was detected/total number of samples*100
  Det_freq[i]=length(which(!is.na(CEC_detected_aver[,i])))/dim(CEC_detected_aver)[1]*100
}

#Summarizing compound occurance
#Table with all stats
CEC_stats = data.frame(round(Det_freq,digits=1),round(Median_v,digits=1),round(Mean_v,digits=1),
                       round(Stdev_v,digits=1),round(cv_v,digits=1),round(Min_v,digits=1),round(Max_v,digits=1))
colnames(CEC_stats) = c("Detec.Freq.","Median","Mean","St.Dev.","CV","Min","Max")

#Table with info about compound use (following the same order as in CEC_stats and CEC_detected_aver)
CEC_order = CEC_use[match(CEC_detected_names,CEC_use_names),]
CEC_use_info = data.frame(CEC_order[,1:2],round(Det_freq,digits=1),CEC_order[,5])
colnames(CEC_use_info) = c("Name","Type","Detec. Freq.","Use")

#----------------
#### TOTAL NUMBER OF CECs AND CUMULATIVE SAMPLE CONCENTRATIONS

#Calculating number of compounds and cumulative concentration in each sample
Nr_per_sampl = Cum_conc = 0
Com_types_array = matrix(0,dim(CEC_detected_aver)[1],length(unique(CEC_use_info[,2]))) #6 CEC types 
colnames(Com_types_array) = c("PCHC","PCHC TP","PEST","PEST TP","PHAR","PHAR TP")
for(i in seq(from=1,to=dim(CEC_detected_aver)[1])){
  Nr_per_sampl[i] = dim(CEC_detected_aver)[2] - length(which(is.na(CEC_detected_aver[i,]))) 
  #How many compounds of each type?
  Com_types_array[i,] = table(CEC_use_info[which(is.na(CEC_detected_aver[i,])==F),2])
  Cum_conc[i] = sum(CEC_detected_aver[i,],na.rm=T) 
}

#Data frame with site name, sample date, and number of compounds detected in that sample as well as cumulative concnetrations
Sample_nrs_cum_conc = data.frame(CEC_aver[,1:2],Nr_per_sampl,Com_types_array,Cum_conc)
colnames(Sample_nrs_cum_conc) = c("Site","Date","Total Nr",paste(colnames(Com_types_array)),"Cum. Conc.")

#----------------
#### SUMMARIZING SITES

#Fixing names
Sample_nrs_cum_conc[37:38,1] = Sample_nrs_cum_conc[39,1]
Sample_nrs_cum_conc[41:42,1] = Sample_nrs_cum_conc[43,1]

#Summarizing sites
Site_sum_gr = Sample_nrs_cum_conc %>% group_by(Site) 
Site_sum = Site_sum_gr %>%  summarise(`Aver.Total Nr` = round(mean(`Total Nr`)), `Aver.PCHC` = round(mean(`PCHC`)),
                                      `Aver.PCHC TP` = round(mean(`PCHC TP`)), `Aver.PEST` = round(mean(`PEST`)),
                                      `Aver.PHAR` = round(mean(`PHAR`)), `Aver.PHAR TP` = round(mean(`PHAR TP`)),
                                      `Aver.Cum. Conc.` = round(mean(`Cum. Conc.`),digits=2))

#----------------
#### COMPOUNDS DETECTED IN ALL SAMPLES

#Coumpounds (N=4) that were found in each sample
Which_highest = which(Det_freq==100)
Which_comp_highest = colnames(CEC_detected_aver)[Which_highest]
#Concentrations
Aver_highest = CEC_detected_aver[,Which_highest]
#Double check names, as there are two atrazine abbreviations
CEC_highest_names = CEC_detected_names[Which_highest]

#Creating data frame 
CECs_in_all = data.frame(Sample_nrs_cum_conc[,1],CEC_aver[,2],Aver_highest)
colnames(CECs_in_all) = c("Site","Date",paste(CEC_highest_names))

#Summarizing
Site_highestCEC_gr = CECs_in_all %>% group_by(Site) 
Site_highestCEC = Site_highestCEC_gr %>% summarise(`Aver.Atrazine` = round(mean(`Atrazine`),digits=2),
                                                   `St.D.Atrazine` = round(sd(`Atrazine`),digits=2),
                                                   `Aver.Metolachlor` = round(mean(`Metolachlor`),digits=2),
                                                   `St.D.Metolachlor` = round(sd(`Metolachlor`),digits=2),
                                                   `Aver.Atrazine2H` = round(mean(`Atrazine2hydroxy`),digits=2),
                                                   `St.D.Atrazine2H` = round(sd(`Atrazine2hydroxy`),digits=2),
                                                   `Aver.Galaxolidone` = round(mean(`Galaxolidone`),digits=2),
                                                   `St.D.Galaxolidone` = round(sd(`Galaxolidone`),digits=2))
View(Site_highestCEC)

#----------------
#### SUMMARIZING METRICS FOR COMPOUNDS DETECTED IN ALL SAMPLES

#ATRAZINE
Site_highestCEC_atrazine = Site_highestCEC_gr %>% summarise(`Aver. Conc.` = round(mean(`Atrazine`),digits=2),
                                                            `St. D.` = round(sd(`Atrazine`),digits=2))
#Calculating CV 
CV_atrazine = round(Site_highestCEC_atrazine[,3]/Site_highestCEC_atrazine[,2]*100,digits=2)
Site_highestCEC_atrazine2 = data.frame(Site_highestCEC_atrazine,CV_atrazine)
colnames(Site_highestCEC_atrazine2) = c("Site","Aver.Conc","St.Dev.","CV")

#METOLACHLOR
Site_highestCEC_metolachlor = Site_highestCEC_gr %>% summarise(`Aver. Conc.` = round(mean(`Metolachlor`),digits=2),
                                                               `St. D.` = round(sd(`Metolachlor`),digits=2))
CV_metolachlor = round(Site_highestCEC_metolachlor[,3]/Site_highestCEC_metolachlor[,2]*100,digits=2)
Site_highestCEC_metolachlor2 = data.frame(Site_highestCEC_metolachlor,CV_metolachlor)
colnames(Site_highestCEC_metolachlor2) = c("Site","Aver.Conc","St.Dev.","CV")

#ATRAZINE-2-HYDROXY
Site_highestCEC_atrazine2H = Site_highestCEC_gr %>% summarise(`Aver. Conc.` = round(mean(`Atrazine2hydroxy`),digits=2),
                                                              `St. D.` = round(sd(`Atrazine2hydroxy`),digits=2))
CV_atrazine2H = round(Site_highestCEC_atrazine2H[,3]/Site_highestCEC_atrazine2H[,2]*100,digits=2)
Site_highestCEC_atrazine2H2 = data.frame(Site_highestCEC_atrazine2H,CV_atrazine2H)
colnames(Site_highestCEC_atrazine2H2) = c("Site","Aver.Conc","St.Dev.","CV")

#GALAXOLIDONE
Site_highestCEC_galaxolidone = Site_highestCEC_gr %>% summarise(`Aver. Conc.` = round(mean(`Galaxolidone`),digits=2),
                                                                `St. D.` = round(sd(`Galaxolidone`),digits=2))
CV_galaxolidone = round(Site_highestCEC_galaxolidone[,3]/Site_highestCEC_galaxolidone[,2]*100,digits=2)
Site_highestCEC_galaxolidone2 = data.frame(Site_highestCEC_galaxolidone,CV_galaxolidone)
colnames(Site_highestCEC_galaxolidone2) = c("Site","Aver.Conc","St.Dev.","CV")

#Saving data
write.table(CEC_aver,file="CEC_aver_dataset.csv",sep=",",row.names=F,col.names=T,append=T)
write.table(CEC_use_info,file="CEC_use_info_dataset.csv",sep=",",row.names=F,col.names=T,append=T)
