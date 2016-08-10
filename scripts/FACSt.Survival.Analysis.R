## FACSt Survival Analysis ##

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

# 'Survival Graphs' gives you 4 survival graphs
  # All mice included, split by gender, generation, and mitochondrial DNA
  # Output:
    # Four survival graphs in plots

# 'Cox Proportional Hazards Test' as a means of analyzing survival data
  # Cox PH is a good way to analyze survival models
  # The P values of the likelihood a variable is correlated to lifespan is stored in theCoxPHValues
  # Done with untransformed data
  # Output:
    # CoxPHValues: Variable and the P value of it being a factor


#### Loading and Formatting of Data ####

# Loading the data

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

#### Survival Graphs ####

# Creation of graph

library(survival)
FACStShort$SurvObj<-with(FACStShort,Surv(as.numeric(Lifespan),as.numeric(status)==1))
SurvivalFit.One<-survfit(SurvObj~1,data=FACStShort,status==1,conf.type="log-log")
SurvivalFit.Sex<-survfit(SurvObj~Sex,data=FACStShort,status==1,conf.type="log-log")
SurvivalFit.Generation<-survfit(SurvObj~Generation,data=FACStShort,status==1,conf.type="log-log")
SurvivalFit.Mito<-survfit(SurvObj~Mitochondrial,data=FACStShort,status==1,conf.type="log-log")
plot(SurvivalFit.One, xlab = "Time", ylab="Survival")
plot(SurvivalFit.Sex, xlab = "Time", ylab="Survival",col=c("red", "blue"))
legend("topright", c("Female", "Male"), lty=1:2,col=c("red", "blue"))
plot(SurvivalFit.Generation, xlab = "Time", ylab="Survival",col=c("black", "red","blue","green","purple"))
legend("topright",c("G10","G11","G7","G8","G9"),lty=1:2,col=c("black","red","blue","green","purple"))
plot(SurvivalFit.Mito, xlab = "Time", ylab="Survival",col=c("black", "red","blue","green","purple","goldenrod"))
legend("topright",c("ABCD","E","F","G","H","NA"),lty=1:2,col=c("black","red","blue","green","purple","goldenrod"))
rm(SurvivalFit.Sex,SurvivalFit.One,SurvivalFit.Generation,SurvivalFit.Mito)

#### Cox Proportional Hazards Test ####

# Running the Cox PH test using sex, generation, and variable as covariates

a<-c(9:33,37:61,65:89)
CoxPHValues<-c()
for (i in 1:length(a)){
  b<-coxph(Surv(as.numeric(Lifespan),status==1)~Sex+Generation+as.numeric(FACStShort[,a[i]]),data=FACStShort)
  CoxPHValues<-rbind(CoxPHValues,c(colnames(FACStShort)[a[i]],summary(b)$coefficients[6,5]))
}
rm(a,b,i)