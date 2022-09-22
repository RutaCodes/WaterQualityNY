#Ruta Basijokaite
#----------------
#This code creates matrix similar to correlation matrix that shows how many compouds match
#i.e. how many commmon compounds were detected among samples
#----------------

#Upload data
CEC_aver = read.csv(file="CEC_aver_dataset.csv",sep=",",header=T) #array with CEC concentrations, site names and sample dates

#Initiating 
Nr_all_match=matrix(0,dim(CEC_aver)[1],dim(CEC_aver)[1])
#Labeling 
colnames(Nr_all_match) = row.names(Nr_all_match) = substring(as.character(CEC_aver[,1]),1,7)

#Counting number of common compounds (total)
for (a in seq(1:89)){
  for (i in seq(a+1,90)){
    count=0
    for (j in seq(3,81)) {
      if (is.na(CEC_aver[a,j]) || is.na(CEC_aver[i,j])){}
      else{count=count+1}
    } 
    Nr_all_match[a,i]=count
    Nr_all_match[i,a]=count #mirror matrix
  }
}
Nr_all_match_diag_rem=Nr_all_match
diag(Nr_all_match_diag_rem)=NA

#Visualizing results
library('plot.matrix')
plot(Nr_all_match_diag_rem,breaks = c(1,5,10,15,20,50),main="Number of common compounds",axis.col=NULL, axis.row=NULL,xlab='', ylab='')
axis(1, at = seq(1,dim(Nr_all_match_diag_rem)[1], by=4), labels = 
       colnames(Nr_all_match_diag_rem)[seq(1,dim(Nr_all_match_diag_rem)[1], by=4)],cex.axis=0.7, las = 2)
axis(2, at = seq(1,dim(Nr_all_match_diag_rem)[1], by=4), labels = 
       row.names(Nr_all_match_diag_rem)[seq(dim(Nr_all_match_diag_rem)[1],1, by=-4)],cex.axis=0.7, las = 1)

