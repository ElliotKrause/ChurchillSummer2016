## QTL Mapping of untransformed data ##

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

# 'Basic QTL Mapping' does the QTL mapping of the untransformed data
  # Prerequisite:
    # Must run 'Pre-QTL' to load SNPs, probs data, and create kinship matrix
  # f1 has been matched to the probs data and all mice without a sequencing have been removed
  # e1 is the probs data after all mice not in f1 have been removed
  # addcovar1 has sex and generation as predictors
  # each data column is transformed into numeric(was character type beforehand)
  # a QTL test for each of the variables is run. completed QTL tests are stored as 'QTL1[insert variable name]'
  # then all QTLs are plotted
  # Output:
    # f1: matched FACStShort data
    # e1: matched probs data
    # addcovar1: covariate for QTL with sex and generation as predictors
    # all the QTL variables(warning: takes time to run, ~6 hours)
      # will be dispalyed as QTL1Lifespan, QTL1CD8..T.Cells.85w, etc (for all 76 variables!)

# 'Permutations for QTL Graph' run the permutations that help determine significant peaks
  # Prerequisite:
    # 'Basic QTL Mapping'
      # Which requires 'Pre-QTL' run before
  # Permutations are useful for determining what QTL peaks are significant
  # They are computationaly intensive though and take a while to run
  # Output:
    # perms: contains permutation for Lifespan
     # number of times its run now is 100, recommended number is 1000
    # thr: contains thresholds for p<.05, <.1, and <.63
    # plots the permutation on top of QTLLifespan plot
  # Note: You should loop the permutations and plot them onto respective QTLs and save to pdf

# 'Coef Plots' plots the allele effects at a specific chromosome
  # Prerequisite:
    # 'Basic QTL Mapping'
      # Which requires 'Pre-QTL' run before
  # Coef plots allow you to view allele effects on each chromosome
  # They way I wrote it, is you specify the chromosome and it makes plots for all the variables
  # Output:
    # pdf of coef plots at 1 specific chromosome

# 'pxg plot' is a cool way of viewing the effects of different allele combinations of a gene
  # Prerequisite:
    # 'Basic QTL Mapping'
      # Which requires 'Pre-QTL' run before
    # Recommended: Having 'Coef Plots' done will allow you to view what chromosomes of a trait have more distinct intervals
  # pxg plots allow you to view the distribution of a phenotype across the 36 possible DO genotypes
  # the interval finds the most significant chromosome in a coefplot(what we run above) and looks
  # at the 36 founder strain combinations of that allele and displays the phenotype of that combination
  # Output:
    # pxg plot(can be looped)

# 'GWAS scan' is genome wide association mapping
  # Prerequisite:
    # 'Basic QTL Mapping'
      # Which requires 'Pre-QTL' run before
    # Also requires a 'DO_Sanger_SDPs.txt.bgz' file to be loaded
  # Really cool, I just never was able to get it to work
  # Output:
    # None, however if you figure out whats wrong with my code, it'll output a gwas scan

# 'Association Mapping' lets you see the significance of the genes in an interval(from coef scan)
  # Prerequisite:
    # 'Basic QTL Mapping'
      # Which requires 'Pre-QTL' run before
    # Recommended: Having 'Coef Plots' done will allow you to see where a good interval to study is
  # The way I have it written now is you first select the QTL map you want to focus in on by changing
  # the qtl in the interval and then the pheno.col
  # While it can only look at 1 QTL at a time, by changing the chr vector, you can look at multiple
  # chromosomes on a QTL at once
  # By changing your thr(threshold) value, you can change the cutoff for significant genes
  # Output:
    # association maps saved in a pdf

# 'Save assoc map genes to a text file' saves all genes in interval to text file
  # Prerequisite:
    # 'Basic QTL Mapping'
      # Which requires 'Pre-QTL' run before
  # By saving the genes to a text file, it makes it easier to search them in a browser later
  # Really smart to run this after running association mapping so you get a list of genes that are on assoc map
  # Output:
    # Text file of genes in an interval on a specified chromosome

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

#### pre-QTL ####

# Loading the data

library(DOQTL)
load("/Users/s-krauss/Desktop/Data/probs.Rdata")
load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))

# Creating Kinship matrix

K = kinship.probs(probs, snps = MM_snps, bychr = TRUE)

#### Basic QTL mapping ####

# Matching up mouse IDs

a1<-rownames(probs)
b1<-FACStShort[,1]
c1<-which(a1 %in% b1)
d1<-which(b1 %in% a1)
e1<-probs[c1,,]
f1<-FACStShort[d1,]
rm(a1,b1,c1,d1)

# Creating covar, not looped

rownames(FACStShort)<-FACStShort[,1]
rownames(f1)<-f1[,1]
addcovar1=model.matrix(~Sex+Generation,data=FACStShort)[,-1]
colnames(addcovar1)[1]="Sex"

# Transforming data into numeric, not looped

a<-c(9:33,37:61,65:90)
for (i in 1:length(a)){
  f1[a[i]] <- sapply( f1[a[i]], as.numeric )
}

# Running QTL mapping, looped

for (i in 1:length(a)){
  qtl = scanone(pheno = f1, pheno.col = a[i], probs = e1, K = K, 
                addcovar = addcovar1, snps = MM_snps)
  assign(paste0("QTL1",colnames(f1)[a[i]]), qtl)
}

# Loop that plots all QTLs

pdf(file = "Basic QTL Plots")
for (i in 1:length(a)){
  plot(eval(parse(text=paste0("QTL1",colnames(f1)[a[i]]))),main=paste0("QTL",colnames(f1)[a[i]]))
}
dev.off()
rm(i,a)

#### Permutations for QTL graph ####

# Create the permutation

perms = scanone.perm(pheno = f1, pheno.col = "Lifespan", probs = e1,
                     addcovar = addcovar1, snps = MM_snps, nperm = 100)

# Plot permutation on graph

thr = get.sig.thr(perms, alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)#perms[,,1]
plot(QTL1Lifespan, sig.thr = thr, sig.col = c("red", "orange", "goldenrod"), main = "Lifespan")

#### Coef Plots #### 

a<-c(9:33,37:61,65:90)
pdf(file = "Coef Plots Chr 18")
for (i in 1:length(a)){
  coefplot(eval(parse(text=paste0("QTL1",colnames(f1)[a[i]]))),chr=18,main=paste0("QTL",colnames(f1)[a[i]]))
}
dev.off()
rm(a)

#### pxg plot ####

# Used for picking out interval of interest, looped, not now

interval = bayesint(QTL1Lifespan, chr = 18)
knitr::kable(interval)

# Used for lookin at more closely one of the markers that the interval & knit above creates

pxg.plot(pheno = f1, pheno.col = "Lifespan", probs = e1,
         snp.id = interval[2,1], snps = MM_snps)
rm(interval)

#### GWAS scan ####

#gwas = scanone.assoc(pheno = f1, pheno.col = "Lifespan", probs = e1, K = K, 
#                     addcovar = addcovar1, markers = MM_snps,  sdp.file = "/Users/s-krauss/Desktop/DO_Sanger_SDPs.txt.bgz", ncl = 3)

#### Association Mapping ####

# Runs the Association Mapping

chr<-c(6:9)
pdf(file = "Assoc Mapping")
for (i in 1:length(chr)){
  interval = bayesint(QTL1Lifespan, chr = chr[i])
  knitr::kable(interval)
  
  assoc = assoc.map(pheno = f1, pheno.col ="Lifespan", probs = e1, K = K[[chr[i]]],
                    addcovar = addcovar1, snps = MM_snps, chr = chr[i], start = interval[1,3],
                    end = interval[3,3], output = "p-value")
  
  tmp = assoc.plot(assoc, thr = 4, show.sdps = TRUE)
}
dev.off()
rm(chr,interval,assoc,tmp)

#### Save assoc map genes to a text file ####

interval = bayesint(QTL1Lifespan, chr = 6)
knitr::kable(interval)

plcholder = get.mgi.features(chr = interval[1,2], start = interval[1,3],
                             end = interval[3,3], type = "gene", source = "MGI")

plc <- plcholder[,10]

write(plc, file = "Lifespan Genes Chr 6")

rm(interval,plcholder,plc)
