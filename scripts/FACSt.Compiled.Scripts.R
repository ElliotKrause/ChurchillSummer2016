#### Finalized Project ####

## Author: Samuel Elliot Krause
# Institution: University of North Carolina at Chapel Hill
# Project Date: 6/8/16 through 8/12/16

# Purpose: To identify relationships between lifespan and healthy aging.

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
      # FACStShortLog currently being removed at end of 'ANOVA Test 4' with rm(), you must manually
      # comment this part out in order to get this output before running 'Whisker Plot Outliers'
  # Runs whisker plots for each variable and finds all the outlier mice and stores
  # the number of times each mouse is an outlier in a file called 'WhiskerOutliers'
  # Also saves each whisker plot to a pdf file
  # Output:
    # WhiskerOutliers: the number of times each mouse is an outlier
    # Whisker Plots saved in pdf

# 'Outliers using Prediction Interval' finds outliers using a prediction interval
  # Prerequisite:
    # Must run 'ANOVA Test 4' for FACStShortLog
      # FACStShortLog currently being removed at end of 'ANOVA Test 4' with rm(), you must manually
      # comment this part out in order to get this output before running 'Whisker Plot Outliers'
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

# 'Cox Proportional Hazards Test' as a means of analyzing survival data
  # Cox PH is a good way to analyze survival models
  # The P values of the likelihood a variable is correlated to lifespan is stored in theCoxPHValues
  # Done with untransformed data
  # Output:
    # CoxPHValues: Variable and the P value of it being a factor

# 'Comparison between mice that died before and after 365 days' uses a T test
  # Prerequisite:
    # So mice that die are automatically removed in the 'Loading and Formatting of Data' section
      # Comment out the section 'Removing all mice dead before 1 year(564 mice to 539)' to keep these mice
  # This uses a T test to compare the mice that died before and after 365 days
  # P value: likelihood of significance
  # Mean before & after: before is mean of the mice that died before a year, after is after
  # SE: Standard Error before and after
  # Output:
    # TTest365: gives ttest comparing mice with death before and after 1 year

# 'Loading SNP Data' loads different SNP data from Annotation Hub
  # So this loads data from annotation hub
    # Pulls out the sets Ensembl, Mus Musculus , and gtf
  # you also only look at the most recent set, AH47076
  # you also only look at coding strands(CDS)
  # Outputs:
    # CDS: the coding strands

# 'Overlapping SNPs' pulls out the SNPs that are both above the threshold and intersecting with coding strands
  # Prerequisite:
    # 'Basic QTL Mapping'
      # Which requires 'Pre-QTL' run before
    # 'Loading SNP Data'
  # So what this does is it looks at 1 chromosome from 1 qtl map and
  # pulls out all genes above a certain threshhold. And with those genes
  # it intersects them with the CDS snps. The CDS SNPs consists of exons and
  # the other coding parts, however doesn't include introns. This reduces the
  # candidate list of genes significantly.
  # Output:
    # writes the significant genes to a text file

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

library(ggplot2)
months<-c(6,12,18)
COI<-c(9:33)#Valid Values: 9:33
pdf(file="Discrete Plots")
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
  labs(x="Months", y=colnames(m12M)[COI[i]])+
  theme_bw())

}
dev.off()
rm(AvgValues,AvgValuesF,AvgValuesM,m12M,m12F,m18F,m18M,m6F,m6M,COI,i,means,months,sex,std)

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
rm(FACStShortLog)

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

#### Cox Proportional Hazards Test ####

# Running the Cox PH test using sex, generation, and variable as covariates

a<-c(9:33,37:61,65:89)
CoxPHValues<-c()
for (i in 1:length(a)){
  b<-coxph(Surv(as.numeric(Lifespan),status==1)~Sex+Generation+as.numeric(FACStShort[,a[i]]),data=FACStShort)
  CoxPHValues<-rbind(CoxPHValues,c(colnames(FACStShort)[a[i]],summary(b)$coefficients[6,5]))
}
rm(a,b,i)

#### Comparison between mice that died before and after 365 days ####

# Remove all mice that died before 365 days

FACStBefore365<-c()
FACStAfter365<-c()
for (i in 1:nrow(FACStShort)){
  if (as.numeric(FACStShort[i,90])<365){
    FACStBefore365<-rbind(FACStBefore365,FACStShort[i,])
  }
  else{
    FACStAfter365<-rbind(FACStAfter365,FACStShort[i,])
  }
}
rm(i)

# T test between variables of pre and post 365 mice

a<-c(9:33)
TTest365<-c()
for (i in 1:length(a)){
  TTest365<-rbind(TTest365,c(colnames(FACStShort)[a[i]],t.test(as.numeric(FACStBefore365[,a[i]]),as.numeric(FACStAfter365[,a[i]]))$p.value,mean(as.numeric(na.omit(FACStBefore365[,a[i]]))),mean(as.numeric(na.omit(FACStAfter365[,a[i]]))),sd(as.numeric(na.omit(FACStBefore365[,a[i]]))/(length(as.numeric(na.omit(FACStBefore365[,a[i]]))))^.5),sd(as.numeric(na.omit(FACStAfter365[,a[i]]))/(length(as.numeric(na.omit(FACStAfter365[,a[i]]))))^.5)))
}
colnames(TTest365)<-c("variables","p-values","mean before","mean after","se before","se after")
rm(a,i,FACStAfter365,FACStBefore365)

#### Loading SNP data ####

#source("http://bioconductor.org/biocLite.R")
#biocLite("AnnotationHub")
library(AnnotationHub)

#Just want coding sequences so you can intersect with QTL data
hub=AnnotationHub()
hub=query(hub,c("Ensembl","Mus Musculus","gtf"))
ensembl=hub[["AH47076"]]

#Want coding sequences
CDS=ensembl[ensembl$type=="CDS"]
rm(ensembl,hub)

#### Overlapping SNPs ####

# CD62L B Cell 33 Chr 15
c<-15 #chr
t<-3.2 #thr

interval = bayesint(QTL1CD62L..B.Cells.33w, chr = c)
knitr::kable(interval)

PLSWORK = get.mgi.features(chr = interval[1,2], start = interval[1,3],
                           end = interval[3,3], type = "gene", source = "MGI")


assoc = assoc.map(pheno = f1, pheno.col ="CD62L..B.Cells.33w", probs = e1, K = K[[c]],
                  addcovar = addcovar1, snps = MM_snps, chr = c, start = interval[1,3],
                  end = interval[3,3], output = "p-value")
assocplot <- assoc.plot(assoc, thr = t, show.sdps = TRUE)
assocplot2<-assocplot
assocplot2[,2]<-assocplot[,2]*1e6

aa<-GRanges(assocplot2$chr,range=IRanges(start=assocplot2[,2],width=1))

bb<-subsetByOverlaps(aa,CDS)

bbs<-start(bb)
bbe<-end(bb)

SigGene<-c()
for (i in 1:length(bbs))
{
  for (j in 1:nrow(PLSWORK))
  {
    if ((bbs[i]>=PLSWORK[j,4])&(bbs[i]<=PLSWORK[j,5]))
    {
      SigGene<-c(SigGene,PLSWORK[j,10])
    }
  }
}
SigGene<-unique(SigGene)

write(SigGene, file = "Significant Genes CD62L B Cells 33w Chr 15")
rm(SigGene,bbe,bbs,aa,bb,PLSWORK,interval, assoc,assocplot,assocplot2,c,t)




