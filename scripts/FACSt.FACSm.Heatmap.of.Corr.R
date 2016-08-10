## FACSt and FACSm heatmaps for Corr with lifespan ##

# Runs independently from other code
# Will need to manually change location of files for each graph though

#### FACSt Both Genders ####

library(RColorBrewer)
options(stringsAsFactors = FALSE)

FACSt <- read.csv("/Users/s-krauss/Desktop/data2/shock_FACS_T_long_JCMS.csv")
Lifespan <- read.csv("/Users/s-krauss/Desktop/data2/shock_lifespan_long.csv")
lifespan_data<-Lifespan
last.FACSt<-FACSt
colnames(last.FACSt)[9:33]<-c("Single Cells", "Viable Cells", "All Lymphocytes","B Cells", "Viable B Cells", "CD62L B Cells", "CD4 T Cells", "Viable CD4 T Cells", "CD4 Effector Memory T Cells", "CD4 Central Memory T Cells", "CD4 Naive T Cells", "CD4 Effect  T Cells", "CD8 T Cells","Viable CD8 T Cells", "CD8 Effector Memory T Cells","CD8 Central Memory T Cells", "CD8 Naive T Cells","CD8 Effect T Cells", "NK Cells", "Viable NK Cells", "NKG2 T Cells","Viable NKG2 T Cells","NKG2 CD4 T Cells","CD8 NKG2 T Cells","NKG2 DN T Cells")

#mouse.id.ear.marks <- str_split_fixed(lifespan_data$Mouse.ID, " ", n = 2)
#lifespan_data <- cbind(mouse.id.ear.marks, lifespan_data)
colnames(lifespan_data)[1] <- "Mouse.ID"
#lifespan_data <- lifespan_data[,-c(2:4)]

a <- last.FACSt$Mouse.ID
b <- lifespan_data$Mouse.ID
c <- which(a %in% b)
d <- which(b %in% a)
last.FACSt <- last.FACSt[c,]
lifespan_data <- lifespan_data[d,]
rm(a,b,c,d)


m = merge(last.FACSt, lifespan_data, by = "Mouse.ID", all = FALSE)

#9:33

m<-m[,c(-4:-6,-34:-74)]
colnames(m)[31]<-"Lifespan"

m$Lifespan = as.numeric(m$Lifespan)
#m$BW = as.numeric(m$BW)
#m$Glucose = as.numeric(m$Glucose)
#m$IGF1..ng.mL. = as.numeric(m$IGF1..ng.mL.)
#colnames(m)[11] <- "IGF1"



# Keep mice that lived > 365 days.
m = m[m$Lifespan >= 365,]

rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}

# should I transform lifespan too?
a<-c(6:31)
for (i in 1:length(a)){ m[,a[i]]<-rz.transform(as.numeric(m[,a[i]]))}
rm(i,rz.transform,a)


#as.numeric(gsub("[^0-9]", "", m$Test.Name[i]))
pheno = colnames(m)[6:30]

cor = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
               list(c(6,12,18), pheno))
pv = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
              list(c(6,12,18), pheno))

m<-m[!is.na(m$Mouse.ID),]

for (i in 1:nrow(m))
{
  if (as.numeric(gsub("[^0-9]", "", m$Test.Name[i]))==85)
  {
    m$Age.at.Exp.Date[i]<-18
  }
  else if (as.numeric(gsub("[^0-9]", "", m$Test.Name[i]))==59)
  {
    m$Age.at.Exp.Date[i]<-12
  }
  else if (as.numeric(gsub("[^0-9]", "", m$Test.Name[i]))==33)
  {
    m$Age.at.Exp.Date[i]<-6
  }
}

time = unique(m$Age.at.exp.date)

layout(matrix(1:length(cor), nrow(cor), ncol(cor)))

par(plt = c(0,1,0,1))
for(i in 1:length(pheno)) {
  
  spl = split(m[,c("Mouse.ID", "Sex.x", "Generation.x", "Lifespan", pheno[i])], m$Test.Name)
  for(j in 1:length(spl)) 
  {
    
    mod1 = lm(Lifespan ~ Sex.x + Generation.x, data = spl[[j]],#Sex.x + Generation.x
              na.action = na.exclude)
    spl[[j]][,4] = residuals(mod1)
    mod2 = lm(spl[[j]][,5] ~ Sex.x + Generation.x, data = spl[[j]],#Sex.x + Generation.x
              na.action = na.exclude)
    spl[[j]][,5] = residuals(mod2)
    
    t = cor.test(spl[[j]][,4], spl[[j]][,5])
    cor[j,i] = t$estimate
    pv[j,i]  = t$p.value
    plot(spl[[j]][,5], spl[[j]][,4], pch = 16)
    
  }
  
}

colnames(cor) = colnames(m)[c(6:30)]
colnames(cor) = c("Single Cells", "Viable Cells", "All Lymphocytes", "B Cells","Viable B Cells","CD62L B Cells","CD4 T Cells","Viable CD4 T Cells","CD4 Effector Memory T Cells","CD4 Central Memory T Cells","CD4 Naive T Cells","CD4 Effect T Cells","CD8 T Cells","Viable CD8 T Cells","CD8 Effector Memory T Cells","CD8 Central Memory T Cells","CD8 Naive T Cells","CD8 Effect T Cells","NK Cells","Viable NK Cells","NKG2ACE T Cell","Viable NKG2ACE","NKG2ACE CD4 T Cells","NKG2ACE CD8 T Cells","NKG2ACE DN T Cells")


ord = order(cor[3,])
ord1<-ord
cor = cor[,ord]
pv  = pv[,ord]

png("FACStCorrelationLifespanBothSexes.png", width = 1500, height = 1000,#1500
    res = 164)
layout(matrix(1:2, 1, 2), width = c(0.85, 0.15))
brks = -400:400/1000
col = rev(colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(brks) - 1))
par(plt = c(0.15, 0.95, 0.05, 0.92))
image(1:nrow(cor), 1:ncol(cor), cor, breaks = brks, col = col, axes = F, 
      ann = F)
mtext(side = 3, line = 0.25, at = 1:3, text = paste(c(6, 12, 18), "Mo"))
mtext(side = 2, at = 1:ncol(cor), text = colnames(cor), las = 1,cex=.6)
for(i in 1:ncol(cor)) {
  text(x = 1:3, y = rep(i, 3), labels = format(pv[,i], digits = 3))
} # for(i)

par(plt = c(0, 0.5, 0.05, 0.92))
image(matrix(brks, nrow = 1), breaks = brks, col = col, axes = F,
      ann = F)
at = seq(1,length(brks), 50)
mtext(side = 4, line = 0.25, at = at / length(brks), text = brks[at], las = 1)
mtext(side = 3, line = 0.25, text = "Cor")
dev.off()
rm(cor,FACSt,last.FACSt,Lifespan,lifespan_data,m,pv,at,brks,col,i,j,mod1,mod2,ord,pheno,spl,t,time)


#### FACSt Male ####

library(RColorBrewer)
options(stringsAsFactors = FALSE)

FACSt <- read.csv("/Users/s-krauss/Desktop/data2/shock_FACS_T_long_JCMS.csv")
Lifespan <- read.csv("/Users/s-krauss/Desktop/data2/shock_lifespan_long.csv")
lifespan_data<-Lifespan
last.FACSt<-FACSt
colnames(last.FACSt)[9:33]<-c("Single Cells", "Viable Cells", "All Lymphocytes","B Cells", "Viable B Cells", "CD62L B Cells", "CD4 T Cells", "Viable CD4 T Cells", "CD4 Effector Memory T Cells", "CD4 Central Memory T Cells", "CD4 Naive T Cells", "CD4 Effect  T Cells", "CD8 T Cells","Viable CD8 T Cells", "CD8 Effector Memory T Cells","CD8 Central Memory T Cells", "CD8 Naive T Cells","CD8 Effect T Cells", "NK Cells", "Viable NK Cells", "NKG2 T Cells","Viable NKG2 T Cells","NKG2 CD4 T Cells","CD8 NKG2 T Cells","NKG2 DN T Cells")

#mouse.id.ear.marks <- str_split_fixed(lifespan_data$Mouse.ID, " ", n = 2)
#lifespan_data <- cbind(mouse.id.ear.marks, lifespan_data)
colnames(lifespan_data)[1] <- "Mouse.ID"
#lifespan_data <- lifespan_data[,-c(2:4)]

a <- last.FACSt$Mouse.ID
b <- lifespan_data$Mouse.ID
c <- which(a %in% b)
d <- which(b %in% a)
last.FACSt <- last.FACSt[c,]
lifespan_data <- lifespan_data[d,]
rm(a,b,c,d)


m = merge(last.FACSt, lifespan_data, by = "Mouse.ID", all = FALSE)

#9:33

m<-m[,c(-4:-6,-34:-74)]
colnames(m)[31]<-"Lifespan"

m$Lifespan = as.numeric(m$Lifespan)
#m$BW = as.numeric(m$BW)
#m$Glucose = as.numeric(m$Glucose)
#m$IGF1..ng.mL. = as.numeric(m$IGF1..ng.mL.)
#colnames(m)[11] <- "IGF1"



# Keep mice that lived > 365 days.
m = m[m$Lifespan >= 365,]
m = m[m$Sex.x =="M",]

rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}

# should I transform lifespan too?
a<-c(6:31)
for (i in 1:length(a)){ m[,a[i]]<-rz.transform(as.numeric(m[,a[i]]))}
rm(i,rz.transform,a)


#as.numeric(gsub("[^0-9]", "", m$Test.Name[i]))
pheno = colnames(m)[6:30]

cor = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
               list(c(6,12,18), pheno))
pv = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
              list(c(6,12,18), pheno))

m<-m[!is.na(m$Mouse.ID),]

for (i in 1:nrow(m))
{
  if (as.numeric(gsub("[^0-9]", "", m$Test.Name[i]))==85)
  {
    m$Age.at.Exp.Date[i]<-18
  }
  else if (as.numeric(gsub("[^0-9]", "", m$Test.Name[i]))==59)
  {
    m$Age.at.Exp.Date[i]<-12
  }
  else if (as.numeric(gsub("[^0-9]", "", m$Test.Name[i]))==33)
  {
    m$Age.at.Exp.Date[i]<-6
  }
}

time = unique(m$Age.at.exp.date)

layout(matrix(1:length(cor), nrow(cor), ncol(cor)))

par(plt = c(0,1,0,1))
for(i in 1:length(pheno)) {
  
  spl = split(m[,c("Mouse.ID", "Sex.x", "Generation.x", "Lifespan", pheno[i])], m$Test.Name)
  for(j in 1:length(spl)) 
  {
    
    mod1 = lm(Lifespan ~Generation.x, data = spl[[j]],#Sex.x + Generation.x
              na.action = na.exclude)
    spl[[j]][,4] = residuals(mod1)
    mod2 = lm(spl[[j]][,5] ~Generation.x, data = spl[[j]],#Sex.x + Generation.x
              na.action = na.exclude)
    spl[[j]][,5] = residuals(mod2)
    
    t = cor.test(spl[[j]][,4], spl[[j]][,5])
    cor[j,i] = t$estimate
    pv[j,i]  = t$p.value
    plot(spl[[j]][,5], spl[[j]][,4], pch = 16)
    
  }
  
}

colnames(cor) = c("Single Cells", "Viable Cells", "All Lymphocytes", "B Cells","Viable B Cells","CD62L B Cells","CD4 T Cells","Viable CD4 T Cells","CD4 Effector Memory T Cells","CD4 Central Memory T Cells","CD4 Naive T Cells","CD4 Effect T Cells","CD8 T Cells","Viable CD8 T Cells","CD8 Effector Memory T Cells","CD8 Central Memory T Cells","CD8 Naive T Cells","CD8 Effect T Cells","NK Cells","Viable NK Cells","NKG2ACE T Cell","Viable NKG2ACE","NKG2ACE CD4 T Cells","NKG2ACE CD8 T Cells","NKG2ACE DN T Cells")

ord <- ord1
cor = cor[,ord]
pv  = pv[,ord]

png("FACStCorrelationLifespanMale.png", width = 1500, height = 1000,#1500
    res = 164)
layout(matrix(1:2, 1, 2), width = c(0.85, 0.15))
brks = -400:400/1000
col = rev(colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(brks) - 1))
par(plt = c(0.15, 0.95, 0.05, 0.92))
image(1:nrow(cor), 1:ncol(cor), cor, breaks = brks, col = col, axes = F, 
      ann = F)
mtext(side = 3, line = 0.25, at = 1:3, text = paste(c(6, 12, 18), "Mo"))
mtext(side = 2, at = 1:ncol(cor), text = colnames(cor), las = 1,cex=.6)
for(i in 1:ncol(cor)) {
  text(x = 1:3, y = rep(i, 3), labels = format(pv[,i], digits = 3))
} # for(i)

par(plt = c(0, 0.5, 0.05, 0.92))
image(matrix(brks, nrow = 1), breaks = brks, col = col, axes = F,
      ann = F)
at = seq(1,length(brks), 50)
mtext(side = 4, line = 0.25, at = at / length(brks), text = brks[at], las = 1)
mtext(side = 3, line = 0.25, text = "Cor")
dev.off()
rm(cor,FACSt,last.FACSt,Lifespan,lifespan_data,m,pv,at,brks,col,i,j,mod1,mod2,ord,pheno,spl,t,time)

#### FACSt Female ####

library(RColorBrewer)
options(stringsAsFactors = FALSE)

FACSt <- read.csv("/Users/s-krauss/Desktop/data2/shock_FACS_T_long_JCMS.csv")
Lifespan <- read.csv("/Users/s-krauss/Desktop/data2/shock_lifespan_long.csv")
lifespan_data<-Lifespan
last.FACSt<-FACSt
colnames(last.FACSt)[9:33]<-c("Single Cells", "Viable Cells", "All Lymphocytes","B Cells", "Viable B Cells", "CD62L B Cells", "CD4 T Cells", "Viable CD4 T Cells", "CD4 Effector Memory T Cells", "CD4 Central Memory T Cells", "CD4 Naive T Cells", "CD4 Effect  T Cells", "CD8 T Cells","Viable CD8 T Cells", "CD8 Effector Memory T Cells","CD8 Central Memory T Cells", "CD8 Naive T Cells","CD8 Effect T Cells", "NK Cells", "Viable NK Cells", "NKG2 T Cells","Viable NKG2 T Cells","NKG2 CD4 T Cells","CD8 NKG2 T Cells","NKG2 DN T Cells")

#mouse.id.ear.marks <- str_split_fixed(lifespan_data$Mouse.ID, " ", n = 2)
#lifespan_data <- cbind(mouse.id.ear.marks, lifespan_data)
colnames(lifespan_data)[1] <- "Mouse.ID"
#lifespan_data <- lifespan_data[,-c(2:4)]

a <- last.FACSt$Mouse.ID
b <- lifespan_data$Mouse.ID
c <- which(a %in% b)
d <- which(b %in% a)
last.FACSt <- last.FACSt[c,]
lifespan_data <- lifespan_data[d,]
rm(a,b,c,d)


m = merge(last.FACSt, lifespan_data, by = "Mouse.ID", all = FALSE)

#9:33

m<-m[,c(-4:-6,-34:-74)]
colnames(m)[31]<-"Lifespan"

m$Lifespan = as.numeric(m$Lifespan)
#m$BW = as.numeric(m$BW)
#m$Glucose = as.numeric(m$Glucose)
#m$IGF1..ng.mL. = as.numeric(m$IGF1..ng.mL.)
#colnames(m)[11] <- "IGF1"



# Keep mice that lived > 365 days.
m = m[m$Lifespan >= 365,]
m = m[m$Sex.x =="F",]

rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}

# should I transform lifespan too?
a<-c(6:31)
for (i in 1:length(a)){ m[,a[i]]<-rz.transform(as.numeric(m[,a[i]]))}
rm(i,rz.transform,a)


#as.numeric(gsub("[^0-9]", "", m$Test.Name[i]))
pheno = colnames(m)[6:30]

cor = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
               list(c(6,12,18), pheno))
pv = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
              list(c(6,12,18), pheno))

m<-m[!is.na(m$Mouse.ID),]

for (i in 1:nrow(m))
{
  if (as.numeric(gsub("[^0-9]", "", m$Test.Name[i]))==85)
  {
    m$Age.at.Exp.Date[i]<-18
  }
  else if (as.numeric(gsub("[^0-9]", "", m$Test.Name[i]))==59)
  {
    m$Age.at.Exp.Date[i]<-12
  }
  else if (as.numeric(gsub("[^0-9]", "", m$Test.Name[i]))==33)
  {
    m$Age.at.Exp.Date[i]<-6
  }
}

time = unique(m$Age.at.exp.date)

layout(matrix(1:length(cor), nrow(cor), ncol(cor)))

par(plt = c(0,1,0,1))
for(i in 1:length(pheno)) {
  
  spl = split(m[,c("Mouse.ID", "Sex.x", "Generation.x", "Lifespan", pheno[i])], m$Test.Name)
  for(j in 1:length(spl)) 
  {
    
    mod1 = lm(Lifespan ~Generation.x, data = spl[[j]],#Sex.x + Generation.x
              na.action = na.exclude)
    spl[[j]][,4] = residuals(mod1)
    mod2 = lm(spl[[j]][,5] ~Generation.x, data = spl[[j]],#Sex.x + Generation.x
              na.action = na.exclude)
    spl[[j]][,5] = residuals(mod2)
    
    t = cor.test(spl[[j]][,4], spl[[j]][,5])
    cor[j,i] = t$estimate
    pv[j,i]  = t$p.value
    plot(spl[[j]][,5], spl[[j]][,4], pch = 16)
    
  }
  
}

colnames(cor) = c("Single Cells", "Viable Cells", "All Lymphocytes", "B Cells","Viable B Cells","CD62L B Cells","CD4 T Cells","Viable CD4 T Cells","CD4 Effector Memory T Cells","CD4 Central Memory T Cells","CD4 Naive T Cells","CD4 Effect T Cells","CD8 T Cells","Viable CD8 T Cells","CD8 Effector Memory T Cells","CD8 Central Memory T Cells","CD8 Naive T Cells","CD8 Effect T Cells","NK Cells","Viable NK Cells","NKG2ACE T Cell","Viable NKG2ACE","NKG2ACE CD4 T Cells","NKG2ACE CD8 T Cells","NKG2ACE DN T Cells")

ord <- ord1
cor = cor[,ord]
pv  = pv[,ord]

png("FACStCorrelationLifespanFemale.png", width = 1500, height = 1000,#1500
    res = 164)
layout(matrix(1:2, 1, 2), width = c(0.85, 0.15))
brks = -400:400/1000
col = rev(colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(brks) - 1))
par(plt = c(0.15, 0.95, 0.05, 0.92))
image(1:nrow(cor), 1:ncol(cor), cor, breaks = brks, col = col, axes = F, 
      ann = F)
mtext(side = 3, line = 0.25, at = 1:3, text = paste(c(6, 12, 18), "Mo"))
mtext(side = 2, at = 1:ncol(cor), text = colnames(cor), las = 1,cex=.6)
for(i in 1:ncol(cor)) {
  text(x = 1:3, y = rep(i, 3), labels = format(pv[,i], digits = 3))
} # for(i)

par(plt = c(0, 0.5, 0.05, 0.92))
image(matrix(brks, nrow = 1), breaks = brks, col = col, axes = F,
      ann = F)
at = seq(1,length(brks), 50)
mtext(side = 4, line = 0.25, at = at / length(brks), text = brks[at], las = 1)
mtext(side = 3, line = 0.25, text = "Cor")
dev.off()
rm(cor,FACSt,last.FACSt,Lifespan,lifespan_data,m,pv,at,brks,col,i,j,mod1,mod2,ord,pheno,spl,t,time)

#### FACSm Both Genders ####

library(RColorBrewer)
options(stringsAsFactors = FALSE)

FACSm <- read.csv("/Users/s-krauss/Desktop/data2/shock_FACS_M_long_JCMS.csv")
Lifespan <- read.csv("/Users/s-krauss/Desktop/data2/shock_lifespan_long.csv")
lifespan_data<-Lifespan
last.FACSm<-FACSm
colnames(last.FACSm)[9:20]<-c("Single Cells", "Viable Cells", "CD11b", "Eosinophils", "Viable Eosinophils", "Granulocytes", "Viable Granulocytes", "Monocytes", "Viable Monocytes", "Inflammatory Monocytes", "Other Monocytes", "Resident Monocytes")

#mouse.id.ear.marks <- str_split_fixed(lifespan_data$Mouse.ID, " ", n = 2)
#lifespan_data <- cbind(mouse.id.ear.marks, lifespan_data)
colnames(lifespan_data)[1] <- "Mouse.ID"
#lifespan_data <- lifespan_data[,-c(2:4)]

a <- last.FACSm$Mouse.ID
b <- lifespan_data$Mouse.ID
c <- which(a %in% b)
d <- which(b %in% a)
last.FACSm <- last.FACSm[c,]
lifespan_data <- lifespan_data[d,]
rm(a,b,c,d)


m = merge(last.FACSm, lifespan_data, by = "Mouse.ID", all = FALSE)

#9:33

m<-m[,c(-4:-6,-21:-61)]
colnames(m)[18]<-"Lifespan"

m$Lifespan = as.numeric(m$Lifespan)
#m$BW = as.numeric(m$BW)
#m$Glucose = as.numeric(m$Glucose)
#m$IGF1..ng.mL. = as.numeric(m$IGF1..ng.mL.)
#colnames(m)[11] <- "IGF1"



# Keep mice that lived > 365 days.
m = m[m$Lifespan >= 365,]

rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}

# should I transform lifespan too?
a<-c(6:18)
for (i in 1:length(a)){ m[,a[i]]<-rz.transform(as.numeric(m[,a[i]]))}
rm(i,rz.transform,a)



pheno = colnames(m)[6:17]

cor = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
               list(c(6,12,18), pheno))
pv = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
              list(c(6,12,18), pheno))
m<-m[!is.na(m$Mouse.ID),]

for (i in 1:nrow(m))
{
  if (m$Test.Name[i]=="FCM (M) 85wks")
  {
    m$Age.at.Exp.Date[i]<-18
  }
  else if (m$Test.Name[i]=="FCM (M) 59wks")
  {
    m$Age.at.Exp.Date[i]<-12
  }
  else if (m$Test.Name[i]=="FCM (M) 33wks")
  {
    m$Age.at.Exp.Date[i]<-6
  }
}

time = unique(m$Age.at.exp.date)

layout(matrix(1:length(cor), nrow(cor), ncol(cor)))

par(plt = c(0,1,0,1))
for(i in 1:length(pheno)) {
  
  spl = split(m[,c("Mouse.ID", "Sex.x", "Generation.x", "Lifespan", pheno[i])], m$Test.Name)
  for(j in 1:length(spl)) 
  {
    
    mod1 = lm(Lifespan ~ Sex.x + Generation.x, data = spl[[j]],#Sex.x + Generation.x
              na.action = na.exclude)
    spl[[j]][,4] = residuals(mod1)
    mod2 = lm(spl[[j]][,5] ~ Sex.x + Generation.x, data = spl[[j]],#Sex.x + Generation.x
              na.action = na.exclude)
    spl[[j]][,5] = residuals(mod2)
    
    t = cor.test(spl[[j]][,4], spl[[j]][,5])
    cor[j,i] = t$estimate
    pv[j,i]  = t$p.value
    plot(spl[[j]][,5], spl[[j]][,4], pch = 16)
    
  }
  
}


colnames(cor) = colnames(m)[c(6:17)]
cor1<-c("Single Cells", "Viable Cells", "All Lymphocytes", "B Cells","Viable B Cells","CD62L B Cells","CD4 T Cells","Viable CD4 T Cells","CD4 Effector Memory T Cells","CD4 Central Memory T Cells","CD4 Naive T Cells","CD4 Effect T Cells","CD8 T Cells","Viable CD8 T Cells","CD8 Effector Memory T Cells","CD8 Central Memory T Cells","CD8 Naive T Cells","CD8 Effect T Cells","NK Cells","Viable NK Cells","NKG2ACE T Cell","Viable NKG2ACE","NKG2ACE CD4 T Cells","NKG2ACE CD8 T Cells","NKG2ACE DN T Cells")

ord = order(cor[3,])
ord1<-ord
cor = cor[,ord]
pv  = pv[,ord]

png("FACSmCorrelationLifespanBothSexes.png", width = 1500, height = 1000,#1500
    res = 164)
layout(matrix(1:2, 1, 2), width = c(0.85, 0.15))
brks = -350:350/1000
col = rev(colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(brks) - 1))
par(plt = c(0.15, 0.95, 0.05, 0.92))
image(1:nrow(cor), 1:ncol(cor), cor, breaks = brks, col = col, axes = F, 
      ann = F)
mtext(side = 3, line = 0.25, at = 1:3, text = paste(c(6, 12, 18), "Mo"))
mtext(side = 2, at = 1:ncol(cor), text = cor1, las = 1,cex=.6)
for(i in 1:ncol(cor)) {
  text(x = 1:3, y = rep(i, 3), labels = format(pv[,i], digits = 3))
} # for(i)

par(plt = c(0, 0.5, 0.05, 0.92))
image(matrix(brks, nrow = 1), breaks = brks, col = col, axes = F,
      ann = F)
at = seq(1,length(brks), 50)
mtext(side = 4, line = 0.25, at = at / length(brks), text = brks[at], las = 1)
mtext(side = 3, line = 0.25, text = "Cor")
dev.off()
rm(cor, FACSm,last.FACSm,Lifespan,lifespan_data,m,pv,at,brks,col,i,j,mod1,mod2,ord,pheno,spl,t,time)

#### FACSm Males ####

library(RColorBrewer)
options(stringsAsFactors = FALSE)

FACSm <- read.csv("/Users/s-krauss/Desktop/data2/shock_FACS_M_long_JCMS.csv")
Lifespan <- read.csv("/Users/s-krauss/Desktop/data2/shock_lifespan_long.csv")
lifespan_data<-Lifespan
last.FACSm<-FACSm
colnames(last.FACSm)[9:20]<-c("Single Cells", "Viable Cells", "CD11b", "Eosinophils", "Viable Eosinophils", "Granulocytes", "Viable Granulocytes", "Monocytes", "Viable Monocytes", "Inflammatory Monocytes", "Other Monocytes", "Resident Monocytes")

#mouse.id.ear.marks <- str_split_fixed(lifespan_data$Mouse.ID, " ", n = 2)
#lifespan_data <- cbind(mouse.id.ear.marks, lifespan_data)
colnames(lifespan_data)[1] <- "Mouse.ID"
#lifespan_data <- lifespan_data[,-c(2:4)]

a <- last.FACSm$Mouse.ID
b <- lifespan_data$Mouse.ID
c <- which(a %in% b)
d <- which(b %in% a)
last.FACSm <- last.FACSm[c,]
lifespan_data <- lifespan_data[d,]
rm(a,b,c,d)


m = merge(last.FACSm, lifespan_data, by = "Mouse.ID", all = FALSE)

#9:33

m<-m[,c(-4:-6,-21:-61)]
colnames(m)[18]<-"Lifespan"

m$Lifespan = as.numeric(m$Lifespan)
#m$BW = as.numeric(m$BW)
#m$Glucose = as.numeric(m$Glucose)
#m$IGF1..ng.mL. = as.numeric(m$IGF1..ng.mL.)
#colnames(m)[11] <- "IGF1"



# Keep mice that lived > 365 days.
m = m[m$Lifespan >= 365,]
m<-m[m$Sex.x == "M",]

rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}

# should I transform lifespan too?
a<-c(6:18)
for (i in 1:length(a)){ m[,a[i]]<-rz.transform(as.numeric(m[,a[i]]))}
rm(i,rz.transform,a)



pheno = colnames(m)[6:17]

cor = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
               list(c(6,12,18), pheno))
pv = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
              list(c(6,12,18), pheno))
m<-m[!is.na(m$Mouse.ID),]

for (i in 1:nrow(m))
{
  if (m$Test.Name[i]=="FCM (M) 85wks")
  {
    m$Age.at.Exp.Date[i]<-18
  }
  else if (m$Test.Name[i]=="FCM (M) 59wks")
  {
    m$Age.at.Exp.Date[i]<-12
  }
  else if (m$Test.Name[i]=="FCM (M) 33wks")
  {
    m$Age.at.Exp.Date[i]<-6
  }
}

time = unique(m$Age.at.exp.date)

layout(matrix(1:length(cor), nrow(cor), ncol(cor)))

par(plt = c(0,1,0,1))
for(i in 1:length(pheno)) {
  
  spl = split(m[,c("Mouse.ID", "Sex.x", "Generation.x", "Lifespan", pheno[i])], m$Test.Name)
  for(j in 1:length(spl)) 
  {
    
    mod1 = lm(Lifespan ~ Generation.x, data = spl[[j]],#Sex.x + Generation.x
              na.action = na.exclude)
    spl[[j]][,4] = residuals(mod1)
    mod2 = lm(spl[[j]][,5] ~ Generation.x, data = spl[[j]],#Sex.x + Generation.x
              na.action = na.exclude)
    spl[[j]][,5] = residuals(mod2)
    
    t = cor.test(spl[[j]][,4], spl[[j]][,5])
    cor[j,i] = t$estimate
    pv[j,i]  = t$p.value
    plot(spl[[j]][,5], spl[[j]][,4], pch = 16)
    
  }
  
}


colnames(cor) = colnames(m)[c(6:17)]
cor1<-c("Single Cells", "Viable Cells", "All Lymphocytes", "B Cells","Viable B Cells","CD62L B Cells","CD4 T Cells","Viable CD4 T Cells","CD4 Effector Memory T Cells","CD4 Central Memory T Cells","CD4 Naive T Cells","CD4 Effect T Cells","CD8 T Cells","Viable CD8 T Cells","CD8 Effector Memory T Cells","CD8 Central Memory T Cells","CD8 Naive T Cells","CD8 Effect T Cells","NK Cells","Viable NK Cells","NKG2ACE T Cell","Viable NKG2ACE","NKG2ACE CD4 T Cells","NKG2ACE CD8 T Cells","NKG2ACE DN T Cells")

ord <-ord1
cor = cor[,ord]
pv  = pv[,ord]

png("FACSmCorrelationLifespanMales.png", width = 1500, height = 1000,#1500
    res = 164)
layout(matrix(1:2, 1, 2), width = c(0.85, 0.15))
brks = -350:350/1000
col = rev(colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(brks) - 1))
par(plt = c(0.15, 0.95, 0.05, 0.92))
image(1:nrow(cor), 1:ncol(cor), cor, breaks = brks, col = col, axes = F, 
      ann = F)
mtext(side = 3, line = 0.25, at = 1:3, text = paste(c(6, 12, 18), "Mo"))
mtext(side = 2, at = 1:ncol(cor), text = cor1, las = 1,cex=.6)
for(i in 1:ncol(cor)) {
  text(x = 1:3, y = rep(i, 3), labels = format(pv[,i], digits = 3))
} # for(i)

par(plt = c(0, 0.5, 0.05, 0.92))
image(matrix(brks, nrow = 1), breaks = brks, col = col, axes = F,
      ann = F)
at = seq(1,length(brks), 50)
mtext(side = 4, line = 0.25, at = at / length(brks), text = brks[at], las = 1)
mtext(side = 3, line = 0.25, text = "Cor")
dev.off()
rm(cor, FACSm,last.FACSm,Lifespan,lifespan_data,m,pv,at,brks,col,i,j,mod1,mod2,ord,pheno,spl,t,time)


#### FACSm Females ####

library(RColorBrewer)
options(stringsAsFactors = FALSE)

FACSm <- read.csv("/Users/s-krauss/Desktop/data2/shock_FACS_M_long_JCMS.csv")
Lifespan <- read.csv("/Users/s-krauss/Desktop/data2/shock_lifespan_long.csv")
lifespan_data<-Lifespan
last.FACSm<-FACSm
colnames(last.FACSm)[9:20]<-c("Single Cells", "Viable Cells", "CD11b", "Eosinophils", "Viable Eosinophils", "Granulocytes", "Viable Granulocytes", "Monocytes", "Viable Monocytes", "Inflammatory Monocytes", "Other Monocytes", "Resident Monocytes")

#mouse.id.ear.marks <- str_split_fixed(lifespan_data$Mouse.ID, " ", n = 2)
#lifespan_data <- cbind(mouse.id.ear.marks, lifespan_data)
colnames(lifespan_data)[1] <- "Mouse.ID"
#lifespan_data <- lifespan_data[,-c(2:4)]

a <- last.FACSm$Mouse.ID
b <- lifespan_data$Mouse.ID
c <- which(a %in% b)
d <- which(b %in% a)
last.FACSm <- last.FACSm[c,]
lifespan_data <- lifespan_data[d,]
rm(a,b,c,d)


m = merge(last.FACSm, lifespan_data, by = "Mouse.ID", all = FALSE)

#9:33

m<-m[,c(-4:-6,-21:-61)]
colnames(m)[18]<-"Lifespan"

m$Lifespan = as.numeric(m$Lifespan)
#m$BW = as.numeric(m$BW)
#m$Glucose = as.numeric(m$Glucose)
#m$IGF1..ng.mL. = as.numeric(m$IGF1..ng.mL.)
#colnames(m)[11] <- "IGF1"



# Keep mice that lived > 365 days.
m = m[m$Lifespan >= 365,]
m<-m[m$Sex.x == "F",]

rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}

# should I transform lifespan too?
a<-c(6:18)
for (i in 1:length(a)){ m[,a[i]]<-rz.transform(as.numeric(m[,a[i]]))}
rm(i,rz.transform,a)



pheno = colnames(m)[6:17]

cor = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
               list(c(6,12,18), pheno))
pv = matrix(0, nrow = 3, ncol = length(pheno), dimnames = 
              list(c(6,12,18), pheno))
m<-m[!is.na(m$Mouse.ID),]

for (i in 1:nrow(m))
{
  if (m$Test.Name[i]=="FCM (M) 85wks")
  {
    m$Age.at.Exp.Date[i]<-18
  }
  else if (m$Test.Name[i]=="FCM (M) 59wks")
  {
    m$Age.at.Exp.Date[i]<-12
  }
  else if (m$Test.Name[i]=="FCM (M) 33wks")
  {
    m$Age.at.Exp.Date[i]<-6
  }
}

time = unique(m$Age.at.exp.date)

layout(matrix(1:length(cor), nrow(cor), ncol(cor)))

par(plt = c(0,1,0,1))
for(i in 1:length(pheno)) {
  
  spl = split(m[,c("Mouse.ID", "Sex.x", "Generation.x", "Lifespan", pheno[i])], m$Test.Name)
  for(j in 1:length(spl)) 
  {
    
    mod1 = lm(Lifespan ~ Generation.x, data = spl[[j]],#Sex.x + Generation.x
              na.action = na.exclude)
    spl[[j]][,4] = residuals(mod1)
    mod2 = lm(spl[[j]][,5] ~ Generation.x, data = spl[[j]],#Sex.x + Generation.x
              na.action = na.exclude)
    spl[[j]][,5] = residuals(mod2)
    
    t = cor.test(spl[[j]][,4], spl[[j]][,5])
    cor[j,i] = t$estimate
    pv[j,i]  = t$p.value
    plot(spl[[j]][,5], spl[[j]][,4], pch = 16)
    
  }
  
}


colnames(cor) = colnames(m)[c(6:17)]
cor1<-c("Single Cells", "Viable Cells", "All Lymphocytes", "B Cells","Viable B Cells","CD62L B Cells","CD4 T Cells","Viable CD4 T Cells","CD4 Effector Memory T Cells","CD4 Central Memory T Cells","CD4 Naive T Cells","CD4 Effect T Cells","CD8 T Cells","Viable CD8 T Cells","CD8 Effector Memory T Cells","CD8 Central Memory T Cells","CD8 Naive T Cells","CD8 Effect T Cells","NK Cells","Viable NK Cells","NKG2ACE T Cell","Viable NKG2ACE","NKG2ACE CD4 T Cells","NKG2ACE CD8 T Cells","NKG2ACE DN T Cells")

ord <- ord1
cor = cor[,ord]
pv  = pv[,ord]

png("FACSmCorrelationLifespanFemales.png", width = 1500, height = 1000,#1500
    res = 164)
layout(matrix(1:2, 1, 2), width = c(0.85, 0.15))
brks = -350:350/1000
col = rev(colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(brks) - 1))
par(plt = c(0.15, 0.95, 0.05, 0.92))
image(1:nrow(cor), 1:ncol(cor), cor, breaks = brks, col = col, axes = F, 
      ann = F)
mtext(side = 3, line = 0.25, at = 1:3, text = paste(c(6, 12, 18), "Mo"))
mtext(side = 2, at = 1:ncol(cor), text = cor1, las = 1,cex=.6)
for(i in 1:ncol(cor)) {
  text(x = 1:3, y = rep(i, 3), labels = format(pv[,i], digits = 3))
} # for(i)

par(plt = c(0, 0.5, 0.05, 0.92))
image(matrix(brks, nrow = 1), breaks = brks, col = col, axes = F,
      ann = F)
at = seq(1,length(brks), 50)
mtext(side = 4, line = 0.25, at = at / length(brks), text = brks[at], las = 1)
mtext(side = 3, line = 0.25, text = "Cor")
dev.off()
rm(cor, FACSm,last.FACSm,Lifespan,lifespan_data,m,pv,at,brks,col,i,j,mod1,mod2,ord,pheno,spl,t,time)
