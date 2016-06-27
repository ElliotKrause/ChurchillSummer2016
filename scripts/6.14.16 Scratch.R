#This is a way of determining the ability of lifespan in to predict each variable
anovatest1<-c()
a<-c(9:33,37:61,65:89)
for (i in 1:length(a)){
m1<-lm(FACStShort[,a[i]]~FACStShort[,2]+FACStShort[,3]+FACStShort[,90])#gender+generation+lifespan
m2<-lm(FACStShort[,a[i]]~FACStShort[,2]+FACStShort[,3])#gender+generation
b<-anova(m1,m2)[2,6]
if (b<.05){
  anovatest1<-rbind(anovatest1,c(colnames(FACStShort)[a[i]],b))
}
}
rm(a,i,m1,m2,b)

#This is a way of determining the ability of a variable to predict lifespan
anovatest2<-c()
a<-c(9:33,37:61,65:89)
for (i in 1:length(a)) {
g<-c(2,3,a[i],90)
pls<-na.omit(FACStShort[,g])
m1<-lm(pls[,4]~pls[,1]+pls[,2]+pls[,3])#gender+generation+variable
m2<-lm(pls[,4]~pls[,1]+pls[,2])#gender+generation
b<-anova(m1,m2)[2,6]
if (b<.05){
  anovatest2<-rbind(anovatest2,c(colnames(FACStShort)[a[i]],b))
}
else {}
}
rm(a,i,m1,m2,g,b,pls)

#This looks at the ability of gender in predicting the different variables
anovatest1sex<-c()
a<-c(9:33,37:61,65:89)
for (i in 1:length(a)){
  m1<-lm(FACStShort[,a[i]]~FACStShort[,2]+FACStShort[,3]+FACStShort[,90])#gender+generation+lifespan
  m2<-lm(FACStShort[,a[i]]~FACStShort[,90]+FACStShort[,3])#generation+lifespan
  b<-anova(m1,m2)[2,6]
  if (b<.05){
    anovatest1sex<-rbind(anovatest1sex,c(colnames(FACStShort)[a[i]],b))
  }
}
rm(a,i,m1,m2,b)

anovatest1gen<-c()
a<-c(9:33,37:61,65:89)
for (i in 1:length(a)){
  m1<-lm(FACStShort[,a[i]]~FACStShort[,2]+FACStShort[,3]+FACStShort[,90])#gender+generation+lifespan
  m2<-lm(FACStShort[,a[i]]~FACStShort[,2]+FACStShort[,90])#gender+generation
  b<-anova(m1,m2)[2,6]
  if (b<.05){
    anovatest1gen<-rbind(anovatest1gen,c(colnames(FACStShort)[a[i]],b))
  }
}
rm(a,i,m1,m2,b)

#for (i in 1:nrow(anovatest1gen)){
#  a<-match(anovatest1gen[i],colnames(FACStShort))
#  plots<-ggplot(data=FACStShort,aes(x=(Generation),y=as.numeric(FACStShort[,a])))+geom_boxplot(notch=TRUE,varwidth = TRUE)+theme_bw()
#  ggsave(plots,filename=paste("GenBoxPlot",anovatest1gen[i,1],".png"))
#}
#rm(plots)
#
#for (i in 1:nrow(anovatest1sex)){
#  a<-match(anovatest1sex[i],colnames(FACStShort))
#  plots<-ggplot(data=FACStShort,aes(x=(Sex),y=as.numeric(FACStShort[,a])))+geom_boxplot(notch=TRUE,varwidth = TRUE)+theme_bw()
#  ggsave(plots,filename=paste("SexBoxPlot",anovatest1gen[i,1],".png"))
#}
#rm(plots,a,i)

# This is a way of reducing the number of variables
FACStShortMale<-FACStShort[FACStShort$Generation=="G7" | FACStShort$Generation=="G8" | FACStShort$Generation=="G9" | FACStShort$Generation=="G11",]
ExtractedDataMale<-as.numeric(FACStShortMale[,a[1]])
for (i in 2:length(a)){
  rrg<-as.numeric(FACStShortMale[,a[i]])
  ExtractedDataMale<-cbind(ExtractedDataMale,rrg)
}
rm(i,rrg)

library(qtlcharts)
colnames(ExtractedDataMale)<-colnames(FACStShortMale[a])
iplotCorr(ExtractedDataMale,reorder=TRUE)


