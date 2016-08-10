## ANOVA and Outlier Detection ##

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

# 'ANOVA Test 1' runs an ANOVA test to see how significant each variable is in predicting lifespan
  # This data has yet to be transformed at all
  # Column titles are Variable and P value
  # Output:
    # ANOVAtest1: significance of variables at predicting lifespan. Untransformed data

# 'RankZ Transformation of Data' does a normal scores transformation of the data
  # normal scores distribution forces the data into bell curve format
  # Output:
    # FACStShortRankZ: rankZ transformed data

# 'ANOVA Test 2' runs an ANOVA test on the rankz transformed data
  # Prerequisite:
    # Must run 'RankZ Transformation of Data' for the FACStShortRankZ
  # Output:
    # ANOVA2: significance of variables at predicting lifespan. Rankz transformed data

# 'ANOVA Test 3' runs an ANOVA test on the RankZ of residuals
  # Creates a linear model using sex and generation, then calculates residuals
  # Residuals are then RankZ transformed and run through ANOVA test
  # Output:
    # ANOVA3: significance of variables at predicting lifespan. RankZ of residuals

# 'ANOVA Test 4' runs an ANOVA test on log base 2 transformed data
  # Transforms data with log base 2
  # Calculates ANOVA
  # Output:
    # ANOVA4: significance of variables at predicting lifespan. Log 2 transform

# 'Whisker Plot Outliers' finds outliers on whisker plots
  # Prerequisite:
    # Must run 'ANOVA Test 4' for FACStShortLog
  # Runs whisker plots for each variable and finds all the outlier mice and stores
  # the number of times each mouse is an outlier in a file called 'WhiskerOutliers'
  # Also saves each whisker plot to a pdf file
  # Output:
    # WhiskerOutliers: the number of times each mouse is an outlier
    # Whisker Plots saved in pdf

# 'Outliers using Prediction Interval' finds outliers using a prediction interval
  # Prerequisite:
    # Must run 'ANOVA Test 4' for FACStShortLog
  # Note: This probably isn't the best nor most efficient way to find outliers, but it looks cool
  # Way this works is for variables with a high correlation it builds a prediction model that can
  # be used to find a range for acceptable values within certain parameters(I used 99%), and then
  # it makes a list of the mice that are outside of the predicted interval. 
  # First thing it does is pull out all the relationships between variables that have a correlation >.50
    # if there isn't a high correlation, predict function won't be of much use
  # Then it goes through each of the variable pairs with high corr values and creates 2 linear models
    # one linear model with x~y, other with y~x
    # predict function finds outliers of the 99% interval in both models, then the outliers are overlapped
      # overlapping mice are stored because its more likely there is something weird with the mouse
  # 'PredictOutlierList' contains a list of all mice that were considered outliers after this test
    # Also contains the column values of the two variables compared in the test
  # 'PredictOutlierFreq' contains number of times each mouse is considered an outlier
    # 1 mouse was an outlier in 29 comparisons!!!
  # Output:
    # PredictOutlierList: List of the mice and the 2 variables(given as column index) where it was an outlier
    # PredictOutlierFreq: Number of times each mouse was an outlier

#### Loading and Formatting of Data ####

FACSt <- read.csv("/Users/s-krauss/Projects/ChurchillSummer2016/data/shock_FACS_T_long_JCMS.csv")
Lifespan <- read.csv("/Users/s-krauss/Projects/ChurchillSummer2016/data/shock_lifespan_long.csv")
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

colnames(FACStShort)<-c("Mouse.ID","Sex","Generation","PenName","DOB","Exp.Date.33w","Age.at.Exp.Date.33w", "Test.Name.33w", "Single.Cells.33w", "Viable.Cells.33w", "All.Lymphocytes.33w", "B.Cells.33w","B.Cells....of.viable..33w","CD62L..B.Cells.33w","CD4..T.Cells.33w","CD4..T.Cells....of.viable..33w","CD4..Effector.Memory.T.Cells.33w","CD4..Central.Memory.T.Cells.33w","CD4..Naive.T.Cells.33w","CD4..Effect.T.Cells.33w","CD8..T.Cells.33w","CD8..T.Cells....of.viable..33w","CD8..Effector.Memory.T.Cells.33w","CD8..Central.Memory.T.Cells.33w","CD8..Naive.T.Cells.33w","CD8..Effect.T.Cells.33w","NK.Cells.33w","NK.Cells....viable..33w","NKG2A.C.E..T.Cell.33w","NKG2A.C.E.....of.viable..33w","NKG2A.C.E..CD4..T.Cells.33w","NKG2A.C.E..CD8..T.Cells.33w","NKG2A.C.E..DN.T.Cells.33w",
                        "Exp.Date.59w","Age.at.Exp.Date.59w", "Test.Name.59w", "Single.Cells.59w", "Viable.Cells.59w", "All.Lymphocytes.59w", "B.Cells.59w","B.Cells....of.viable..59w","CD62L..B.Cells.59w","CD4..T.Cells.59w","CD4..T.Cells....of.viable..59w","CD4..Effector.Memory.T.Cells.59w","CD4..Central.Memory.T.Cells.59w","CD4..Naive.T.Cells.59w","CD4..Effect.T.Cells.59w","CD8..T.Cells.59w","CD8..T.Cells....of.viable..59w","CD8..Effector.Memory.T.Cells.59w","CD8..Central.Memory.T.Cells.59w","CD8..Naive.T.Cells.59w","CD8..Effect.T.Cells.59w","NK.Cells.59w","NK.Cells....viable..59w","NKG2A.C.E..T.Cell.59w","NKG2A.C.E.....of.viable..59w","NKG2A.C.E..CD4..T.Cells.59w","NKG2A.C.E..CD8..T.Cells.59w","NKG2A.C.E..DN.T.Cells.59w",
                        "Exp.Date.85w","Age.at.Exp.Date.85w", "Test.Name.85w", "Single.Cells.85w", "Viable.Cells.85w", "All.Lymphocytes.85w", "B.Cells.85w","B.Cells....of.viable..85w","CD62L..B.Cells.85w","CD4..T.Cells.85w","CD4..T.Cells....of.viable..85w","CD4..Effector.Memory.T.Cells.85w","CD4..Central.Memory.T.Cells.85w","CD4..Naive.T.Cells.85w","CD4..Effect.T.Cells.85w","CD8..T.Cells.85w","CD8..T.Cells....of.viable..85w","CD8..Effector.Memory.T.Cells.85w","CD8..Central.Memory.T.Cells.85w","CD8..Naive.T.Cells.85w","CD8..Effect.T.Cells.85w","NK.Cells.85w","NK.Cells....viable..85w","NKG2A.C.E..T.Cell.85w","NKG2A.C.E.....of.viable..85w","NKG2A.C.E..CD4..T.Cells.85w","NKG2A.C.E..CD8..T.Cells.85w","NKG2A.C.E..DN.T.Cells.85w",
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

#### ANOVA Test 1####

ANOVAtest1<-c()
a<-c(9:33,37:61,65:89)
for (i in 1:length(a)) {
  g<-c(2,3,a[i],90)
  pls<-na.omit(FACStShort[,g])
  m1<-lm(pls[,4]~pls[,1]+pls[,2]+pls[,3])#gender+generation+mdNA+Lifespan
  m2<-lm(pls[,4]~pls[,1]+pls[,2])#gender+generation+Lifespan
  b<-anova(m1,m2)[2,6]
  ANOVAtest1<-rbind(ANOVAtest1,c(colnames(FACStShort)[a[i]],b))
}
colnames(ANOVAtest1)<-c("Variable","P Value")
rm(a,i,m1,m2,g,b,pls)

#### RankZ Transformation of Data ####

# Rank-Z form

rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}

# Rank-Z of original data

FACStShortRankZ<-FACStShort
a<-c(9:33,37:61,65:90)
for (i in 1:length(a)){ FACStShortRankZ[,a[i]]<-rz.transform(as.numeric(FACStShortRankZ[,a[i]]))}
rm(a,i,rz.transform)

#### ANOVA Test 2 ####

ANOVA2<-c()
a<-c(9:33,37:61,65:90)
for (i in 1:length(c(9:33,37:61,65:89))) {
  g<-c(2,3,a[i],90)
  pls<-na.omit(FACStShortRankZ[,g])
  m1<-lm(pls[,4]~pls[,1]+pls[,2]+pls[,3])#gender+generation+variable
  m2<-lm(pls[,4]~pls[,1]+pls[,2])#gender+generation
  b<-anova(m1,m2)[2,6]
  ANOVA2<-rbind(ANOVA2,c(colnames(FACStShortRankZ)[a[i]],b))
}
rm(i,m1,m2,g,b,pls,a)
#rm(FACStShortRankZ)

#### ANOVA Test 3 ####

# Finding residuals of original data

a<-c(9:33,37:61,65:90)
FACStShortRankZRes<-FACStShort
for (i in 1:length(a)){
  m1<-resid(lm(FACStShortRankZRes[,a[i]]~FACStShortRankZRes[,2]+FACStShortRankZRes[,3]))
  b<-c()
  for (j in 1:length(FACStShortRankZRes[,a[i]])){
    c<-unname(m1[paste(j)])
    b<-c(b,c)
  }
  FACStShortRankZRes[,a[i]]<-b
}
rm(b,c,i,j,m1)

# Rank-Z form

rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}

# Rank-Z of residual data

for (i in 1:length(a)){ FACStShortRankZRes[,a[i]]<-rz.transform(as.numeric(FACStShortRankZRes[,a[i]]))}
rm(i,rz.transform)

# ANOVA of residual

ANOVA3<-c()
for (i in 1:length(c(9:33,37:61,65:89))) {
  g<-c(2,3,a[i],90)
  pls<-na.omit(FACStShortRankZRes[,g])
  m1<-lm(pls[,4]~pls[,1]+pls[,2]+pls[,3])#gender+generation+variable
  m2<-lm(pls[,4]~pls[,1]+pls[,2])#gender+generation
  b<-anova(m1,m2)[2,6]
  ANOVA3<-rbind(ANOVA3,c(colnames(FACStShortRankZRes)[a[i]],b))
}
rm(i,m1,m2,g,b,pls,a)
rm(FACStShortRankZRes)

#### ANOVA Test 4 ####

# Log2 Transforming Data

a<-c(9:33,37:61,65:90)
FACStShortLog<-FACStShort
for (i in 1:length(a)){
  FACStShortLog[,a[i]]<-log(as.numeric(FACStShort[,a[i]]),base=2)
}

# ANOVA of log2 transformed data

ANOVA4<-c()
for (i in 1:length(c(9:33,37:61,65:89))) {
  g<-c(2,3,a[i],90)
  pls<-na.omit(FACStShortLog[,g])
  pls[!is.finite(pls[,3]),3]<-0
  pls[!is.finite(pls[,4]),4]<-0
  m1<-lm(pls[,4]~pls[,1]+pls[,2]+pls[,3])#gender+generation+variable
  m2<-lm(pls[,4]~pls[,1]+pls[,2])#gender+generation
  b<-anova(m1,m2)[2,6]
  
  ANOVA4<-rbind(ANOVA4,c(colnames(FACStShortLog)[a[i]],b))
  
}
rm(a,i,b,g,m1,m2,pls)
#rm(FACStShortLog)

#### Whisker Plot Outliers ####

pdf(file='Whisker Plots')
a<-c(9:33,37:61,65:90)
b<-FACStShortLog[,1]
WhiskerOutliers<-cbind(b,rep(0,times=length(FACStShortLog[,1])))
for (i in 1:length(a)){
  b<-boxplot(FACStShortLog[,a[i]]~FACStShortLog[,2]+FACStShortLog[,3])
  c<-match(b$out,FACStShortLog[,a[i]])
  for (j in 1:length(c)){
    WhiskerOutliers[c[j],2]<-as.numeric(WhiskerOutliers[c[j],2])+1
  }
}
dev.off()
rm(b,c,i,j,a)

#### Outliers using Prediction Interval ####

# Correlation pairs with >.50 correlations

a<-c(9:33,37:61,65:90)
stor1<-c()
stor2<-c()
stor3<-c()
for (i in 1:length(a)){
  for (j in 1:length(a)){
    b<-as.numeric(FACStShortLog[,a[i]])
    c<-as.numeric(FACStShortLog[,a[j]])
    b[!is.finite(b)] <- 0
    c[!is.finite(c)] <- 0
    d<-cor(b,c,use="pairwise.complete.obs")
    if ((i>j) & (abs(d)>.5)){
      stor1<-c(stor1,a[i])
      stor2<-c(stor2,a[j])
      stor3<-c(stor3,d)
    }
  }
}
rm(a,b,c,d,i,j)

# Predict plot

PredictOutlierList<-c()
for (i in 1:length(stor1)){
  test1<-c()
  test2<-c()
  Placeholder5<-na.omit(FACStShortLog[,c(1,stor2[i],stor1[i])])
  Placeholder5[,2][is.infinite(Placeholder5[,2])]<-0
  Placeholder5[,3][is.infinite(Placeholder5[,3])]<-0
  model1<-lm(Placeholder5[,3]~Placeholder5[,2],data=Placeholder5)
  model2<-lm(Placeholder5[,2]~Placeholder5[,3],data=Placeholder5)
  pred1<-predict(model1,interval="prediction",level=.99)
  pred2<-predict(model2,interval="prediction",level=.99)
  combine1<-cbind(Placeholder5[,c(1,2,3)],pred1)
  combine2<-cbind(Placeholder5[,c(1,2,3)],pred2)
  for (j in 1:nrow(combine1)){
    if ((combine1[j,3]>=combine1[j,6])|(combine1[j,3]<=combine1[j,5])){
      test1<-rbind(test1,c(combine1[j,1],stor2[i],stor1[i]))
      
    }
    else{}
    if ((combine2[j,2]>=combine2[j,6])|(combine2[j,2]<=combine2[j,5])){
      test1<-rbind(test1,c(combine2[j,1],stor1[i],stor2[i]))
      
    }
    else{}
  }
  b<-test1[duplicated(test1[,1]),1]
  if (length(b)>0){
    for (j in 1:length(b)){
      PredictOutlierList<-rbind(PredictOutlierList,c(test1[duplicated(test1[,1]),1][j],test1[duplicated(test1[,1]),2][j],test1[duplicated(test1[,1]),3][j]))
      
    }
  }
}
rm(stor1,stor2,stor3,test2,b,i,j,model1,model2,combine1,combine2,pred1,pred2,test1,Placeholder5)

# Number of times as an outlier

PredictOutlierFreq<-c()
for (i in 1:length(unique(PredictOutlierList[,1]))){
  PredictOutlierFreq<-rbind(PredictOutlierFreq,c(unique(PredictOutlierList[,1])[i],as.numeric(sum(PredictOutlierList[,1]==unique(PredictOutlierList[,1])[i]))))
}
rm(i)
