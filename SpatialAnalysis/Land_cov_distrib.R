#Ruta Basijokaite
#----------------
#This code calculates land cover percentage distributions at different distances from sampling location
#that is later used to develop distance weighted models to predict detected number of compounds 
#----------------

#Create matrix that will have all percentages recorded at varying distances from sampling locations
#Number of rows is random number, as it will modified based on number of available cells inside for loop
Perc_devel_allWS = Perc_crop_allWS = Perc_pasture_allWS = Perc_agr_allWS = matrix(0,10,20)

#Create char vector with site names
WS_names = c('conesus','chitt','ganarg','chenan','cowas','flint','butter',
             'oatka','honoeye','nineLK','nineMAR','harbor','sauq','oneida',
             'otselic','salmon','orisk','scriba','sixmile','skaneat')
colnames(Perc_devel_allWS) = colnames(Perc_crop_allWS) = colnames(Perc_pasture_allWS) = colnames(Perc_agr_allWS) = WS_names

#----------------
#### Calculating land cover percentages at different distances from sampling location

#Loop all sites
for (j in seq(1:20)){
  #Upload landcover distances - each watershed has a separate file
  WS_lc <- read.delim(file=paste(WS_names[j],"_dist_lc.txt",sep=""),sep=",") #Large files
  
  #Sorting cells based on distance from sampling location
  WS_lc_sort = sort(WS_lc$Distance, decreasing = FALSE,index.return = TRUE)
  WS_lc_sort_mtrx = WS_lc[WS_lc_sort$ix,]
  
  cnt = 0
  perc_devel = perc_crop = perc_pasture = perc_agr = 0
  Dist = WS_lc_sort_mtrx$Distance/1000 #to km
  for (i in seq(from=0.3, to=max(Dist), by=0.1)){
    #Cell index that is closest to 'i' distance from sampling location 
    ind = which(abs(Dist-i)==min(abs(Dist-i)))
    #Matrix with cells within 'i' radius from sampling location
    LC_at_dist = WS_lc[1:ind,]
    cnt=cnt+1
    #Total number of cells within 'i' radius 
    N_lines=dim(LC_at_dist)[1]
    #Calculating land cover percentage: number of cells that match land cover type divided by total number of cells multiplied by 100 
    perc_devel[cnt]=length(which(LC_at_dist$grid_code==21 | LC_at_dist$grid_code==22 | LC_at_dist$grid_code==23 | LC_at_dist$grid_code==24))/N_lines*100
    perc_crop[cnt]=length(which(LC_at_dist$grid_code==82))/N_lines*100
    perc_pasture[cnt]=length(which(LC_at_dist$grid_code==81))/N_lines*100
    perc_agr[cnt]=length(which(LC_at_dist$grid_code==82 | LC_at_dist$grid_code==81))/N_lines*100 #combining cult. crop and pasture
  }
  
  #Combining variables
  vl=seq(from=0.3, to=max(Dist), by=0.1)
  LC_perc_values=cbind(vl,perc_pasture,perc_crop,perc_agr,perc_devel)
  
  #As matrix is created to store percentages of land cover distance 'i' away from sampling location for every site,
  #Matrix needs to be adjusted as each site varies in size, hence number of rows will differ 
  #If site is bigger (i.e. max radius is higher from previous site), add number of rows with NA values
  if (length(perc_crop)>dim(Perc_crop_allWS)[1]){
    Perc_devel_allWS=rbind(Perc_devel_allWS,matrix(NA,length(perc_devel)-dim(Perc_devel_allWS)[1],20))
    Perc_crop_allWS=rbind(Perc_crop_allWS,matrix(NA,length(perc_crop)-dim(Perc_crop_allWS)[1],20))
    Perc_pasture_allWS=rbind(Perc_pasture_allWS,matrix(NA,length(perc_pasture)-dim(Perc_pasture_allWS)[1],20))
    Perc_agr_allWS=rbind(Perc_agr_allWS,matrix(NA,length(perc_agr)-dim(Perc_agr_allWS)[1],20))
  }
  Perc_devel_allWS[1:length(perc_devel),j]=perc_devel
  Perc_crop_allWS[1:length(perc_crop),j]=perc_crop
  Perc_pasture_allWS[1:length(perc_pasture),j]=perc_pasture
  Perc_agr_allWS[1:length(perc_agr),j]=perc_agr
}

#Saving variables
write.table(Perc_devel_allWS,file=paste( "LC_developed_allWS_by_dist.csv"),sep=",",row.names=F,col.names=T,append=T)
write.table(Perc_crop_allWS,file=paste( "LC_crop_allWS_by_dist.csv"),sep=",",row.names=F,col.names=T,append=T)
write.table(Perc_pasture_allWS,file=paste( "LC_pasture_allWS_by_dist.csv"),sep=",",row.names=F,col.names=T,append=T)
write.table(Perc_agr_allWS,file=paste( "LC_agr_allWS_by_dist.csv"),sep=",",row.names=F,col.names=T,append=T)

#Calculate average land cover percentage distributions 
Perc_devel_allWS_mean = apply(Perc_devel_allWS,FUN=mean,na.rm=TRUE, MARGIN=1)
Perc_crop_allWS_mean = apply(Perc_crop_allWS,FUN=mean,na.rm=TRUE, MARGIN=1)
Perc_pasture_allWS_mean = apply(Perc_pasture_allWS,FUN=mean,na.rm=TRUE, MARGIN=1)
Perc_agr_allWS_mean = apply(Perc_agr_allWS,FUN=mean,na.rm=TRUE, MARGIN=1)
Perc_agr_devel_allWS_mean = Perc_devel_allWS_mean + Perc_agr_allWS_mean

#Plot averaged land cover percentage distribution
plot(Perc_devel_allWS_mean)
plot(Perc_crop_allWS_mean)
plot(Perc_pasture_allWS_mean)
plot(Perc_agr_allWS_mean)
plot(Perc_agr_devel_allWS_mean)

