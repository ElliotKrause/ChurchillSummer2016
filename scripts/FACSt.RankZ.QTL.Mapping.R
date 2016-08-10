## RankZ QTL Mapping Setup ##

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

# 'RankZ Transformation of Data' does a normal scores transformation of the data
  # normal scores distribution forces the data into bell curve format
  # Output:
    # FACStShortRankZ: rankZ transformed data


# 'Pre-QTL' loads probability data, SNP data, and creates Kinship matrix
  # Probability data is a 3D array containing the proportion of each founder haplotype at each marker
  # for each DO sample. The samples are in the first dimension, 8 founder strains in second
  # dimension, and markers along mouse genome are in third dimension. 
  # SNP data comes from the Muga array, called MM_snps
  # Kinship array accounts for kinship relationships between mice
  # Output:
    # probs: contains probability data(warning: Large!)
    # MM_snps: contains SNP data
    # K: contains kinship data(warning: takes time to run, ~20 min)

# 'RankZ QTL Mapping' does the QTL mapping of the RankZ data
  # Prerequisite:
    # Must run 'Pre-QTL' to load SNPs, probs data, and create kinship matrix
    # Must also run 'RankZ Transformation of Data' beforehand
  # f3 has been matched to the probs data and all mice without a sequencing have been removed
  # e3 is the probs data after all mice not in f1 have been removed
  # addcovar3 has sex and generation as predictors
  # each data column is transformed into numeric(was character type beforehand)
  # a QTL test for each of the variables is run. completed QTL tests are stored as 'QTL3[insert variable name]'
  # then all QTLs are plotted
  # Output:
    # f3: matched FACStShortRankZ data
    # e1: matched probs data
    # addcovar1: covariate for QTL with sex and generation as predictors
    # all the QTL variables(warning: takes time to run, ~6 hours)
      # will be dispalyed as QTL3Lifespan, QTL3CD8..T.Cells.85w, etc (for all 76 variables!)


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

#### RankZ Transformation of Data ####

# Rank-Z form

rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}

#### pre-QTL ####

# Loading the data

library(DOQTL)
load("/Users/s-krauss/Desktop/Data/probs.Rdata")
load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))

# Creating Kinship matrix

K = kinship.probs(probs, snps = MM_snps, bychr = TRUE)

#### RankZ QTL mapping ####

# Matching up mouse IDs

a3<-rownames(probs)
b3<-FACStShortRankZ[,1]
c3<-which(a3 %in% b3)
d3<-which(b3 %in% a3)
e3<-probs[c3,,]
f3<-FACStShortRankZ[d3,]
rm(a3,b3,c3,d3)

# Creating covar, not looped

rownames(FACStShortRankZ)<-FACStShortRankZ[,1]
rownames(f3)<-f3[,1]
addcovar3=model.matrix(~Sex+Generation,data=FACStShortRankZ)[,-1]
colnames(addcovar3)[1]="Sex"

# Transforming data into numeric, not looped

a<-c(9:33,37:61,65:90)
for (i in 1:length(a)){
  f3[a[i]] <- sapply( f3[a[i]], as.numeric )
}

# Running QTL mapping, looped

for (i in 1:length(a)){
  qtl = scanone(pheno = f3, pheno.col = a[i], probs = e3, K = K, 
                addcovar = addcovar3, snps = MM_snps)
  assign(paste0("QTL3",colnames(f3)[a[i]]), qtl)
}

# Loop that plots all QTLs

pdf(file = "RankZ QTL Plots")
for (i in 1:length(a)){
  plot(eval(parse(text=paste0("QTL3",colnames(f3)[a[i]]))),main=paste0("QTL",colnames(f3)[a[i]]))
}
dev.off()
rm(i,a)
