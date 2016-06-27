rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}
Placeholder2<-FACStShort
a<-c(9:33,37:61,65:90)
for (i in 1:length(a)){ Placeholder2[,a[i]]<-rz.transform(as.numeric(Placeholder2[,a[i]]))}

rankzANOVA<-c()
for (i in 1:length(c(9:33,37:61,65:89))) {
  g<-c(2,3,a[i],90)
  pls<-na.omit(Placeholder2[,g])
  m1<-lm(pls[,4]~pls[,1]+pls[,2]+pls[,3])#gender+generation+variable
  m2<-lm(pls[,4]~pls[,1]+pls[,2])#gender+generation
  b<-anova(m1,m2)[2,6]
  if (b<.05){
    rankzANOVA<-rbind(rankzANOVA,c(colnames(Placeholder2)[a[i]],b))
  }
  else {}
}
rm(i,m1,m2,g,b,pls)

#ExtractedData<-as.numeric(Placeholder2[,a[1]])
#for (i in 2:length(a)){
#  rrg<-as.numeric(Placeholder2[,a[i]])
#  ExtractedData<-cbind(ExtractedData,rrg)
#}
#rm(i,rrg)
#
#library(qtlcharts)
#colnames(ExtractedData)<-colnames(Placeholder2[a])
#iplotCorr(ExtractedData,reorder=TRUE)

##### Garbage Test #####
Placeholder3<-FACStShort
for (i in 1:length(a)){
  m1<-resid(lm(Placeholder3[,a[i]]~Placeholder3[,2]+Placeholder3[,3]))
  b<-c()
  for (j in 1:length(Placeholder3[,a[i]])){
    c<-unname(m1[paste(j)])
    b<-c(b,c)
  }
  Placeholder3[,a[i]]<-b
}
rm(b,c,i,j,m1)
#Actually works?
for (i in 1:length(a)){ Placeholder3[,a[i]]<-rz.transform(as.numeric(Placeholder3[,a[i]]))}
rm(i,rz.transform)

rankzANOVA2<-c()
for (i in 1:length(c(9:33,37:61,65:89))) {
  g<-c(2,3,a[i],90)
  pls<-na.omit(Placeholder3[,g])
  m1<-lm(pls[,4]~pls[,1]+pls[,2]+pls[,3])#gender+generation+variable
  m2<-lm(pls[,4]~pls[,1]+pls[,2])#gender+generation
  b<-anova(m1,m2)[2,6]
  if (b<.05){
    rankzANOVA2<-rbind(rankzANOVA2,c(colnames(Placeholder3)[a[i]],b))
  }
  else {}
}
rm(i,m1,m2,g,b,pls)

#ExtractedData<-as.numeric(Placeholder3[,a[1]])
#for (i in 2:length(a)){
#  rrg<-as.numeric(Placeholder3[,a[i]])
#  ExtractedData<-cbind(ExtractedData,rrg)
#}
#rm(i,rrg)
#
#library(qtlcharts)
#colnames(ExtractedData)<-colnames(Placeholder3[a])
#iplotCorr(ExtractedData,reorder=TRUE)

Placeholder4<-FACStShort
for (i in 1:length(a)){
  Placeholder4[,a[i]]<-log(as.numeric(FACStShort[,a[i]]),base=2)
}

logANOVA<-c()
for (i in 1:length(c(9:33,37:61,65:89))) {
  g<-c(2,3,a[i],90)
  pls<-na.omit(Placeholder4[,g])
  pls[!is.finite(pls[,3]),3]<-0
  pls[!is.finite(pls[,4]),4]<-0
  m1<-lm(pls[,4]~pls[,1]+pls[,2]+pls[,3])#gender+generation+variable
  m2<-lm(pls[,4]~pls[,1]+pls[,2])#gender+generation
  b<-anova(m1,m2)[2,6]
  if (b<.05){
    logANOVA<-rbind(logANOVA,c(colnames(Placeholder4)[a[i]],b))
  }
  else {}
}
rm(i,m1,m2,g,b,pls)

