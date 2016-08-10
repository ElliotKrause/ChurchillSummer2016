## Correlation and Discrete Plots ##

#### User Notes ####

# Before running any of the other sections, run 'Loading and Formatting of data'
  # This will give you the formatted files necessary for other sections
    # One thing to note, within that section all mice that die before 365 are removed
  # You can comment out that line if you'd rather see all mice
  # Output:
    # FACSt: the original data file, each row represents individual experiment
    # FACStLong: the data has been put in the long format
    # FACStShort: short format, 1 row per mouse
    # Lifespan: lifespan of each mouse. Lifespan has already been added to FACStShort

# 'Correlation Plots' use Karl Bromen's package 'qtlcharts'
  # First graph has all mice, and its sorted by correlation
  # Second graph has only male, not sorted
  # Third graph has only femail, not sorted
  # Be careful, first graph you see in viewer will be third graph
  # Output
    # 3 correlation plots in viewer

# 'Discrete Plots' is a good way to view how a phenotype changes with time and differs by gender
  # Data is separated into 6, 12, and 18 months, then split by gender.
  # Then mean and standard error are calculated and plotted. Same gender values are connected with a line.
  # Currently looped to save pdf folder with graphs to somewhere in your directory, however
  # you can plot specific graphs by changing 'COI' (originally columns of interest) to a single value
  # Output:
    # Discrete plots in pdf

#### Loading and Formatting of Data ####

FACSt <- read.csv("/Users/s-krauss/Desktop/data2/shock_FACS_T_long_JCMS.csv")
Lifespan <- read.csv("/Users/s-krauss/Desktop/data2/shock_lifespan_long.csv")
FACSt[626,7]<-229

# Long format

FACStLong <- reshape(FACSt,
                     varying = c("Age.at.Exp.Date", "Single.Cells", "Viable.Cells", "All.Lymphocytes", "B.Cells","B.Cells....of.viable.","CD62L..B.Cells","CD4..T.Cells","CD4..T.Cells....of.viable.","CD4..Effector.Memory.T.Cells","CD4..Central.Memory.T.Cells","CD4..Central.Memory.T.Cells","CD4..Naive.T.Cells","CD4..Effect.T.Cells","CD8..T.Cells","CD8..T.Cells....of.viable.","CD8..Effector.Memory.T.Cells","CD8..Central.Memory.T.Cells","CD8..Naive.T.Cells","CD8..Effect.T.Cells","NK.Cells","NK.Cells....viable.","NKG2A.C.E..T.Cell","NKG2A.C.E.....of.viable.","NKG2A.C.E..CD4..T.Cells","NKG2A.C.E..CD8..T.Cells","NKG2A.C.E..DN.T.Cells"), 
                     v.names = "number",
                     timevar = "cell.type", 
                     times =   c("Age.at.Exp.Date", "Single.Cells", "Viable.Cells", "All.Lymphocytes", "B.Cells","B.Cells....of.viable.","CD62L..B.Cells","CD4..T.Cells","CD4..T.Cells....of.viable.","CD4..Effector.Memory.T.Cells","CD4..Central.Memory.T.Cells","CD4..Central.Memory.T.Cells","CD4..Naive.T.Cells","CD4..Effect.T.Cells","CD8..T.Cells","CD8..T.Cells....of.viable.","CD8..Effector.Memory.T.Cells","CD8..Central.Memory.T.Cells","CD8..Naive.T.Cells","CD8..Effect.T.Cells","NK.Cells","NK.Cells....viable.","NKG2A.C.E..T.Cell","NKG2A.C.E.....of.viable.","NKG2A.C.E..CD4..T.Cells","NKG2A.C.E..CD8..T.Cells","NKG2A.C.E..DN.T.Cells"), 
                     new.row.names = 1:(42980),
                     direction = "long")



FACStLong$Mouse.ID<-as.character(FACStLong$Mouse.ID)
FACStLong$Sex<-as.character(FACStLong$Sex)
FACStLong$Generation<-as.character(FACStLong$Generation)
FACStLong$penName<-as.character(FACStLong$penName)
ShortRowLabels <- unique(FACSt[,1],incomparables = FALSE)
FACStShort <-data.frame()
Placeholder <-data.frame()

for (i in 1:length(ShortRowLabels)) {
  for(j in 1:nrow(FACSt)) {
    if (as.character(ShortRowLabels[i]) == as.character(FACSt[j,1]) )
      Placeholder <- rbind(Placeholder,FACSt[j,])
  }
  Placeholder$Mouse.ID<-as.character(Placeholder$Mouse.ID)
  Placeholder$Sex<-as.character(Placeholder$Sex)
  Placeholder$Generation<-as.character(Placeholder$Generation)
  Placeholder$penName<-as.character(Placeholder$penName)
  nnn<-rep(NA,length=89)
  nnn<-data.frame(t(data.frame(nnn)))
  nnn[1,c(1:5)]<-Placeholder[1,c(1:5)]
  for(k in 1:nrow(Placeholder)) {
    if (Placeholder[k,8]=="FCM (T-B-NK) 33wks"){
      nnn[1,c(6:33)]<-as.character(Placeholder[k,c(6:33)])
    } else {
      if (Placeholder[k,8]=="FCM (T-B-NK) 59wks"){
        nnn[1,c(34:61)]<-as.character(Placeholder[k,c(6:33)])
      } else {
        if (Placeholder[k,8]=="FCM (T-B-NK) 85wks"){
          nnn[1,c(62:89)]<-as.character(Placeholder[k,c(6:33)])
        }}}
  }
  FACStShort <-rbind(FACStShort,nnn)
  Placeholder <-data.frame()
}
rm(i,j,k,ShortRowLabels)

# Add lifespan data

k=1;
for (i in 1:nrow(FACStShort)){
  for (j in 1:nrow(Lifespan)){
    if (FACStShort[i,1]==Lifespan[j,1]){
      Placeholder[k,1]<-toString(Lifespan[j,43])
      k=k+1
    }
  }
}
FACStShort<-cbind(FACStShort,Placeholder)
rm(i,j,k,nnn,Placeholder)

#colnames(FACStShort)<-c("Mouse.ID","Sex","Generation","PenName","DOB","Exp.Date.33w","Age.at.Exp.Date.33w", "Test.Name.33w", "Single.Cells.33w", "Viable.Cells.33w", "All.Lymphocytes.33w", "B.Cells.33w","B.Cells....of.viable..33w","CD62L..B.Cells.33w","CD4..T.Cells.33w","CD4..T.Cells....of.viable..33w","CD4..Effector.Memory.T.Cells.33w","CD4..Central.Memory.T.Cells.33w","CD4..Naive.T.Cells.33w","CD4..Effect.T.Cells.33w","CD8..T.Cells.33w","CD8..T.Cells....of.viable..33w","CD8..Effector.Memory.T.Cells.33w","CD8..Central.Memory.T.Cells.33w","CD8..Naive.T.Cells.33w","CD8..Effect.T.Cells.33w","NK.Cells.33w","NK.Cells....viable..33w","NKG2A.C.E..T.Cell.33w","NKG2A.C.E.....of.viable..33w","NKG2A.C.E..CD4..T.Cells.33w","NKG2A.C.E..CD8..T.Cells.33w","NKG2A.C.E..DN.T.Cells.33w",
#                        "Exp.Date.59w","Age.at.Exp.Date.59w", "Test.Name.59w", "Single.Cells.59w", "Viable.Cells.59w", "All.Lymphocytes.59w", "B.Cells.59w","B.Cells....of.viable..59w","CD62L..B.Cells.59w","CD4..T.Cells.59w","CD4..T.Cells....of.viable..59w","CD4..Effector.Memory.T.Cells.59w","CD4..Central.Memory.T.Cells.59w","CD4..Naive.T.Cells.59w","CD4..Effect.T.Cells.59w","CD8..T.Cells.59w","CD8..T.Cells....of.viable..59w","CD8..Effector.Memory.T.Cells.59w","CD8..Central.Memory.T.Cells.59w","CD8..Naive.T.Cells.59w","CD8..Effect.T.Cells.59w","NK.Cells.59w","NK.Cells....viable..59w","NKG2A.C.E..T.Cell.59w","NKG2A.C.E.....of.viable..59w","NKG2A.C.E..CD4..T.Cells.59w","NKG2A.C.E..CD8..T.Cells.59w","NKG2A.C.E..DN.T.Cells.59w",
#                        "Exp.Date.85w","Age.at.Exp.Date.85w", "Test.Name.85w", "Single.Cells.85w", "Viable.Cells.85w", "All.Lymphocytes.85w", "B.Cells.85w","B.Cells....of.viable..85w","CD62L..B.Cells.85w","CD4..T.Cells.85w","CD4..T.Cells....of.viable..85w","CD4..Effector.Memory.T.Cells.85w","CD4..Central.Memory.T.Cells.85w","CD4..Naive.T.Cells.85w","CD4..Effect.T.Cells.85w","CD8..T.Cells.85w","CD8..T.Cells....of.viable..85w","CD8..Effector.Memory.T.Cells.85w","CD8..Central.Memory.T.Cells.85w","CD8..Naive.T.Cells.85w","CD8..Effect.T.Cells.85w","NK.Cells.85w","NK.Cells....viable..85w","NKG2A.C.E..T.Cell.85w","NKG2A.C.E.....of.viable..85w","NKG2A.C.E..CD4..T.Cells.85w","NKG2A.C.E..CD8..T.Cells.85w","NKG2A.C.E..DN.T.Cells.85w",
#                        "Lifespan")

colnames(FACStShort)<-c("Mouse.ID","Sex","Generation","PenName","DOB","Exp.Date.33w","Age.at.Exp.Date.33w", "Test.Name.33w", "Single Cells 6mon", "Viable Cells 6mon", "All Lymphocytes 6mon", "B Cells 6mon","Viable B Cells 6mon","CD62L B Cells 6mon","CD4 T Cells 6mon","Viable CD4 T Cells 6mon","CD4 Effector Memory T Cells 6mon","CD4 Central Memory T Cells 6mon","CD4 Naive T Cells 6mon","CD4 Effect T Cells 6mon","CD8 T Cells 6mon","Viable CD8 T Cells 6mon","CD8 Effector Memory T Cells 6mon","CD8 Central Memory T Cells 6mon","CD8 Naive T Cells 6mon","CD8 Effect T Cells 6mon","NK Cells 6mon","Viable NK Cells 6mon","NKG2ACE T Cell 6mon","Viable NKG2ACE 6mon","NKG2ACE CD4 T Cells 6mon","NKG2ACE CD8 T Cells 6mon","NKG2ACE DN T Cells 6mon",
                        "Exp Date 12m","Age at Exp Date 12mon", "Test Name 12mon", "Single Cells 12mon", "Viable Cells 12mon", "All Lymphocytes 12mon", "B Cells 12mon","Viable B Cells 12mon","CD62L B Cells 12mon","CD4 T Cells 12mon","Viable CD4 T Cells 12mon","CD4 Effector Memory T Cells 12mon","CD4 Central Memory T Cells 12mon","CD4 Naive T Cells 12mon","CD4 Effect T Cells 12mon","CD8 T Cells 12mon","Viable CD8 T Cells 12mon","CD8 Effector Memory T Cells 12mon","CD8 Central Memory T Cells 12mon","CD8 Naive T Cells 12mon","CD8 Effect T Cells 12mon","NK Cells 12mon","Viable NK Cells 12mon","NKG2ACE T Cell 12mon","Viable NKG2ACE 12mon","NKG2ACE CD4 T Cells 12mon","NKG2ACE CD8 T Cells 12mon","NKG2ACE DN T Cells 12mon",
                        "Exp Date 18m","Age at Exp Date 18mon", "Test Name 18mon", "Single Cells 18mon", "Viable Cells 18mon", "All Lymphocytes 18mon", "B Cells 18mon","Viable B Cells 18mon","CD62L B Cells 18mon","CD4 T Cells 18mon","Viable CD4 T Cells 18mon","CD4 Effector Memory T Cells 18mon","CD4 Central Memory T Cells 18mon","CD4 Naive T Cells 18mon","CD4 Effect T Cells 18mon","CD8 T Cells 18mon","Viable CD8 T Cells 18mon","CD8 Effector Memory T Cells 18mon","CD8 Central Memory T Cells 18mon","CD8 Naive T Cells 18mon","CD8 Effect T Cells 18mon","NK Cells 18mon","Viable NK Cells 18mon","NKG2ACE T Cell 18mon","Viable NKG2ACE 18mon","NKG2ACE CD4 T Cells 18mon","NKG2ACE CD8 T Cells 18mon","NKG2ACE DN T Cells 18mon",
                        "Lifespan")

# Event data

FACStShort$status<-rep(1, times = nrow(FACStShort))
for (i in 1:length(FACStShort[,90])){
  if(is.na(as.numeric(FACStShort[i,90]))==TRUE){
    FACStShort[i,91]<-2
  }
} 
rm(i)
b<-c()
for (i in 1:length(FACStShort[,90])){
  b<-c(b,substr(FACStShort[i,90],1,nchar(FACStShort[,90])))
}
FACStShort[,90]<-b
rm(b,i)

# Removing all mice dead before 1 year(564 mice to 539)

FACStShort<-FACStShort[as.numeric(FACStShort$Lifespan)>365,]

# Adding mitochondrial data

Mito <- read.csv("/Users/s-krauss/Desktop/Data/Raw Data/shock_DO_inventory2.csv")
Placeholder<-data.frame()
k=1;
for (i in 1:nrow(FACStShort)){
  for (j in 1:nrow(Mito)){
    if (FACStShort[i,1]==Mito[j,2]){
      Placeholder[k,1]<-toString(Mito[j,43])
      k=k+1
    }
  }
}
FACStShort<-cbind(FACStShort,Placeholder)
rm(i,j,k,Placeholder,Mito)
colnames(FACStShort)[92]<-c("Mitochondrial")

#### Correlation Plots ####

# Creating a matrix for data
FACStShort1<-FACStShort[FACStShort$Sex == "M",]
FACStShort2<-FACStShort[FACStShort$Sex == "F",]

a<-c(7,9:33,35,37:61,63,65:90)
ExtractedData<-as.numeric(FACStShort[,a[1]])
ExtractedData1<-as.numeric(FACStShort1[,a[1]])
ExtractedData2<-as.numeric(FACStShort2[,a[1]])


for (i in 2:length(a)){
  rrg<-as.numeric(FACStShort[,a[i]])
  ExtractedData<-cbind(ExtractedData,rrg)
  rrg1<-as.numeric(FACStShort1[,a[i]])
  ExtractedData1<-cbind(ExtractedData1,rrg1)
  rrg2<-as.numeric(FACStShort2[,a[i]])
  ExtractedData2<-cbind(ExtractedData2,rrg2)
}
rm(i,rrg,a,rrg1,rrg2)

# Plotting data
a<-c(7,9:33,35,37:61,63,65:90)
library(qtlcharts)
colnames(ExtractedData)<-colnames(FACStShort[a])
colnames(ExtractedData1)<-colnames(FACStShort1[a])
colnames(ExtractedData2)<-colnames(FACStShort2[a])

iplotCorr(ExtractedData,group=FACStShort$Sex,reorder=TRUE)
iplotCorr(ExtractedData1,reorder=FALSE)
iplotCorr(ExtractedData2,reorder=FALSE)
rm(ExtractedData,a,ExtractedData1,ExtractedData2,FACStShort1,FACStShort2)

#### Discrete Plots ####

# Manipulating data for discrete plots

b<-c(1:33,90)
m6<-data.frame(FACStShort[,b])
b<-c(-6:-33,-62:-89,-91:-92)
m12<-data.frame(FACStShort[,b])
b<-c(-6:-61,-91:-92)
m18<-data.frame(FACStShort[,b])
m6M<-data.frame()
m12M<-data.frame()
m18M<-data.frame()
m6F<-data.frame()
m12F<-data.frame()
m18F<-data.frame()
for (i in 1:nrow(FACStShort)){
  if (m6[i,2]=="M"){
    m6M<-rbind(m6M,m6[i,])
  }
  else {
    m6F<-rbind(m6F,m6[i,])
  }
}
for (i in 1:nrow(FACStShort)){
  if (m12[i,2]=="M"){
    m12M<-rbind(m12M,m12[i,])
  }
  else {
    m12F<-rbind(m12F,m12[i,])
  }
}
for (i in 1:nrow(FACStShort)){
  if (m18[i,2]=="M"){
    m18M<-rbind(m18M,m18[i,])
  }
  else {
    m18F<-rbind(m18F,m18[i,])
  }
}
rm(m12, m18,m6,b)

# Caluculating mean and standard error for Discrete Plots
plc<-c("Mouse.ID","Sex","Generation","PenName","DOB","Exp.Date.33w","Age.at.Exp.Date.33w", "Test.Name.33w", "Single Cells", "Viable Cells", "All Lymphocytes", "B Cells","Viable B Cells","CD62L B Cells","CD4 T Cells","Viable CD4 T Cells","CD4 Effector Memory T Cells","CD4 Central Memory T Cells","CD4 Naive T Cells","CD4 Effect T Cells","CD8 T Cells","Viable CD8 T Cells","CD8 Effector Memory T Cells","CD8 Central Memory T Cells","CD8 Naive T Cells","CD8 Effect T Cells","NK Cells","Viable NK Cells","NKG2ACE T Cell","Viable NKG2ACE","NKG2ACE CD4 T Cells","NKG2ACE CD8 T Cells","NKG2ACE DN T Cells")
library(ggplot2)
months<-c(6,12,18)
COI<-c(9:33)#Valid Values: 9:33
pdf(file="Discrete Plots 2")
for (i in 1:length(COI)){
  means<-c(mean(as.numeric(na.omit(m6M[,COI[i]]))),mean(as.numeric(na.omit(m12M[,COI[i]]))),mean(as.numeric(na.omit(m18M[,COI[i]]))))
  std<-c(sd(as.numeric(na.omit(unlist(m6M[,COI[i]]))))/(length(as.numeric(na.omit(unlist(m6M[,COI[i]])))))^.5,sd(as.numeric(na.omit(unlist(m12M[,COI[i]]))))/(length(as.numeric(na.omit(unlist(m12M[,COI[i]])))))^.5,sd(as.numeric(na.omit(unlist(m18M[,COI[i]]))))/(length(as.numeric(na.omit(unlist(m18M[,COI[i]])))))^.5)
  AvgValues<-data.frame(months,means)
  sex<-c("M","M","M")
  AvgValuesM<-cbind(AvgValues,sex,std)
  
  means<-c(mean(as.numeric(na.omit(m6F[,COI[i]]))),mean(as.numeric(na.omit(m12F[,COI[i]]))),mean(as.numeric(na.omit(m18F[,COI[i]]))))
  std<-c(sd(as.numeric(na.omit(unlist(m6F[,COI[i]]))))/(length(as.numeric(na.omit(unlist(m6F[,COI[i]])))))^.5,sd(as.numeric(na.omit(unlist(m12F[,COI[i]]))))/(length(as.numeric(na.omit(unlist(m12F[,COI[i]])))))^.5,sd(as.numeric(na.omit(unlist(m18F[,COI[i]]))))/(length(as.numeric(na.omit(unlist(m18F[,COI[i]])))))^.5)
  AvgValues<-data.frame(months,means)
  sex<-c("F","F","F")
  AvgValuesF<-cbind(AvgValues,sex,std)
  
  # Plotting discrete plots
  
  AvgValues<-data.frame()
  AvgValues<-rbind(AvgValuesM,AvgValuesF)
  print(ggplot(AvgValues,aes(x=months,y=means,colour=sex,group=sex))+geom_line()+
          geom_errorbar(aes(ymin=means-std,ymax=means+std),position=position_dodge(0.9),width=1)+
          labs(x="Months", y=plc[COI[i]])+ggtitle("Change in Phenotype with Age")+
          theme_bw()+theme(axis.text=element_text(size=15,face="bold"),axis.title=element_text(size=15,face="bold")))
  
}
dev.off()
rm(AvgValues,AvgValuesF,AvgValuesM,m12M,m12F,m18F,m18M,m6F,m6M,COI,i,means,months,sex,std,plc)
