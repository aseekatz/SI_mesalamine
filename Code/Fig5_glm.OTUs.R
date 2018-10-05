##### Fig 5: glm, identifying OTUs that correlate with pH or mesalamine
###### Anna M. Seekatz
###### 09.18.18



####### Fig. 5: files used:
	- mesal_otufrac.w.meta.txt
	- mesal_duod_filteredOTUfrac.txt (created)
	- mesal_stomach_filteredOTUfrac.txt (created)
	- mesal.duodenum_gee.results.txt (created)
	- mesal.stomach_gee.results.txt (created)
	
####### Step 1: generalized mixed model of OTUs correlated with pH or mesalamine
	- for duodenum and stomach
	- for each, apply to full data set
	
```
library(TeachingDemos)
library(dplyr)
library(tidyr)
library(geepack)
library(pROC)
library(QICpack)

library(tidyverse)
library(lme4)
library(car)
library(lmerTest) #this overloads the lmer function to give P values
library(ggalluvial)

options(scipen=999)  # turn off scientific notation like 1e+06

### read in OTUs (rel abundance), meta:
data<-read.table(file="mesal_otufrac.w.meta.txt", header=TRUE)
data<-data[data$seq_status %in% c("OK"), ]
summary(data$site)		# need to combine stomach/antrum
data$site[data$site %in% c("Antrum")]<-"Stomach"
	# note: this already has crappy samples removed

### will do this for stomach and duodenum

### first let's look at the distribution of mesalamine and pH
hist(data$pH)
summary(data$pH)
hist(data$Mesalamineconcentration)
summary(data$Mesalamineconcentration)

# looks like mesalamine is definitely not a normal distribution--for that, might want to dichotomize into low/high
# for pH, can probably just use normal range

# plan:
#	- filter out OTUs that are prevalent in at least half of the samples taken from an individual
#	- in duodenum and stomach each, use glm to identify correlated OTUs (with either pH or mesalamine) --> these will be output to a table
#	- then, from each of these results, graph the RA vs. pH or mesalamine for ALL samples to see if it holds against all data

### go ahead and dichotomize Mesalamine concentration now:
data$mesal[data$Mesalamineconcentration < median(data$Mesalamineconcentration, na.rm=T)]<-"low"
data$mesal[data$Mesalamineconcentration > median(data$Mesalamineconcentration, na.rm=T)]<-"high"

data<-cbind(data[, 1:13], data[, 604], data[, 14:603])
names(data)[14]<-"mesal"


###------------
### duodenum first

# eliminate all OTUs with 0 total
rownames(data)<-data$seqID
dfs<-droplevels(data[data$site %in% c("Duodenum"), c(16:ncol(data)-1)])

### now, let's filter out OTUs that occur in at least 50% of the samples BY PATIENT

#tapply(d$seqID, d$pID, length)		# we have 5-8 samples per person

## define data:
df.frac<-dfs[,  which((colSums(dfs)) > 0)]
df<-df.frac
	#convert to binary:
df[df > 0]<-1

df<-merge(data[,1:14], df, by.x="seqID", by.y="row.names")
df<-arrange(df, pID, Hour)
rownames(df)<-df$seqID

# let's set a cutoff in that we want an OTU to appear in at least 50% of the samples available per patient (so by prevalence within a patient)
# function:
listFunction <- function(n) {
	otus <- df[df$pID==(n),  15:length(df)]					# this defines your dataframe for each patient
	otus <- otus[!is.na(otus$Otu0001), ]					# this eliminates samples that are NA
	minnum<-as.numeric(nrow(otus)/2)									# this sets the minimum number of samples you want that OTU to be observed in
	otus_list<-names(which(colSums(otus == 1) >= minnum))	# this gets a list of the OTUs that appear in at least half of the patient samples
	}
# get list of patients and get lists of OTUs that appear in at least 50% of the samples:
patient.list<-as.character(unique(df$pID))
min.otu.list<-lapply(patient.list, FUN=listFunction)
names(min.otu.list)<-patient.list
lapply(min.otu.list, FUN=length)

# now, collapse this list--THIS will be the list we will use to filter out our data set:
otulist<-as.character(unlist(min.otu.list, recursive = FALSE))
otulist<-unique(otulist)

# ok, let's use this list to get the OTUs we want
df.filtered<-df.frac[, colnames(df.frac) %in% otulist]
df.filtered<-merge(data[,1:14], df.filtered, by.x="seqID", by.y="row.names")
df<-arrange(df.filtered, pID, Hour)
#write.table(df, file="mesal_duod_filteredOTUfrac.txt", , sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

### let's test out some modeling on singular OTUs, first:

#### now, use updated code:
### this is basic code that can be ran:
# first testing for 1 OTU:

# testing out intra-individual variation? this is null model (I think)
OTU_lmer<-lmer(formula = Otu0001~ (1|pID), data=df, REML = F)

# add predictor (pH)
OTU.pH_lmer<-lmer(formula = Otu0001~pH + (1|pID), data=df, REML = F)
summary(OTU_lmer)
VarCorr(OTU_lmer)

# see if difference from null model is significantly different
anova(OTU_lmer, OTU.pH_lmer)

# add in sample number (a random variable?)
OTU.pH.hour_lmer<-lmer(formula = Otu0001~pH + (Hour|pID), data=df, REML = F)
summary(OTU.pH.hour_lmer)
VarCorr(OTU.pH.hour_lmer)
	# since the AIC also drops here, can check if statistically different
anova(OTU.pH.hour_lmer, OTU_lmer)		# also significant
AIC(OTU.pH.hour_lmer,OTU_lmer)

# could also run it as a fixed effect?
OTU.pH.hour_lmer<-lmer(formula = Otu0001~pH + Hour + (Hour|pID), data=df, REML = F)
summary(OTU.pH.hour_lmer)
VarCorr(OTU.pH.hour_lmer)		# lookas about same?

# test the differences between these two:
OTU.pH.hour_lmer<-lmer(formula = Otu0001~pH + Hour + (Hour|pID), data=df, REML = F)
temp<-lmer(formula = Otu0001~pH + (Hour|pID), data=df, REML = F)
anova(OTU.pH.hour_lmer, temp)
remove(temp)			# no difference here

### tested out a couple OTUs from GEE modeling:
OTU.pH.hour_lmer<-lmer(formula = Otu0009~ + (Hour|pID), data=df, REML = F)
temp<-lmer(formula = Otu0009~pH + (Hour|pID), data=df, REML = F)
anova(OTU.pH.hour_lmer, temp)
remove(temp)
	# overall, looks pretty good!
	# filtered out a LOT of OTUs that somehow passed the GEE test
	
### now, let's run for all OTUs:
#df<-read.table(file="mesal_duod_filteredOTUfrac.txt", header=T)

start<-match("Otu0001", names(df))
end<-(length(names(df)))

#now let's create the absence / presence of OTUs
catotus<-dplyr::select(df, contains("Otu"))
catotus<-catotus>0
temp<-dplyr::select(df, -contains("Otu"))
temp<-cbind(temp, catotus)
cattest<-temp
remove(temp)
remove(catotus)

#finally, let's sum the rows to get the total number of OTUs present
df$totalotus <- rowSums(cattest[start:end])
summary(df$totalotus)

# run code to model against all OTUs:
results<-matrix(data=NA, nrow = end-start+1, ncol = 3)
colnames(results)<-c("variable","P_value","coefficient")

# glmer(formula = post~I(as.factor(phylopam3_cluster)) + (sample_number|id), family=binomial, data=test, REML = F)

testotu<-df
for(j in start:end)
{
  tryCatch({
    variable<-names(testotu)[j]
    results[j-start+1, 1]<-variable
    results[j-start+1, 2]<-summary(lmer(as.formula(paste(variable, "~", "pH + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Pr(>|t|)"]
    results[j-start+1, 3]<-summary(lmer(as.formula(paste(variable, "~", "pH + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Estimate"] 
  }, error=function(e){}
  )
}

### without trycatch
testotu<-df
for(j in start:end)
{
    variable<-names(testotu)[j]
    results[j-start+1, 1]<-variable
    results[j-start+1, 2]<-summary(lmer(as.formula(paste(variable, "~", "pH + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Pr(>|t|)"]
    results[j-start+1, 3]<-summary(lmer(as.formula(paste(variable, "~", "pH + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Estimate"] 
}
#####
# both of these work!!!!!
# reformat and filter out only significantly correlated results:

temp<-as.data.frame(results)
temp$P_value<-as.numeric(as.character(temp$P_value))
temp$coefficient<-as.numeric(as.character(temp$coefficient))

results<-temp
remove(temp)

taxonomy<-read.table(file="mesalamine.taxnames03.txt", header=TRUE)

sigresults<-as.data.frame(results)
sigresults$P_value<-as.numeric(as.character(sigresults$P_value))
sigresults<-sigresults[which(sigresults$P_value<.001), ]
topresults<-sigresults
topresults$abscoef<-abs(as.numeric(as.character(topresults$coefficient)))
topresults1<-arrange(topresults, desc(abscoef))
topresults1$metabolite<-"pH"

### then, run for mesalamine:

# let's use the dichotomized value
df$mesal <- as.factor(ifelse(df$mesal=="high", 1, 0))
df$mesal<-as.numeric(as.character(df$mesal))
testotu<-df
for(j in start:end)
{
    variable<-names(testotu)[j]
    results[j-start+1, 1]<-variable
    results[j-start+1, 2]<-summary(lmer(as.formula(paste(variable, "~", "mesal + (Hour|pID)")), REML = F, data=testotu, family="binomial"))$coefficients[2,"Pr(>|t|)"]
    results[j-start+1, 3]<-summary(lmer(as.formula(paste(variable, "~", "mesal + (Hour|pID)")), REML = F, data=testotu, family="binomial"))$coefficients[2,"Estimate"] 
}
#Error in eval(expr, envir, enclos) : y values must be 0 <= y <= 1
#In addition: Warning message:
#extra argument(s) ‘REML’ disregarded 

# ok, can't get this to work for now
# will just use the concentration as is...
testotu<-df
for(j in start:end)
{
    variable<-names(testotu)[j]
    results[j-start+1, 1]<-variable
    results[j-start+1, 2]<-summary(lmer(as.formula(paste(variable, "~", "Mesalamineconcentration + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Pr(>|t|)"]
    results[j-start+1, 3]<-summary(lmer(as.formula(paste(variable, "~", "Mesalamineconcentration + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Estimate"] 
}

temp<-as.data.frame(results)
temp$P_value<-as.numeric(as.character(temp$P_value))
temp$coefficient<-as.numeric(as.character(temp$coefficient))

results<-temp
remove(temp)

sigresults<-as.data.frame(results)
sigresults$P_value<-as.numeric(as.character(sigresults$P_value))
sigresults<-sigresults[which(sigresults$P_value<.001), ]
topresults<-sigresults
topresults$abscoef<-abs(as.numeric(as.character(topresults$coefficient)))
topresults2<-arrange(topresults, desc(abscoef))
topresults2$metabolite<-"mesalamine"

### combine:
geedf<-rbind(topresults1, topresults2)
geedf$site<-"Duodenum"

# read in taxonomic file to merge with this:
taxa<-read.table(file="mesalamine.taxnames03.txt", header=TRUE)
# add family names:
tax <- taxa$Taxonomy
tax <- gsub("\\(\\d*\\)", "", tax)
tax <- gsub(";unclassified", "", tax)
tax <- gsub("_1", "", tax)
tax <- gsub(";$", "", tax)
tax <- gsub("/.*", "", tax)
tax <- gsub(".*les;", "", tax)
tax <- gsub(";.*", "", tax)
taxa$family<-tax
ftaxa<-taxa[, c("OTU", "taxname", "phylum", "family", "Size")]
taxanames<-taxa[, c("OTU", "taxname", "phylum", "family", "Size")]
geedf_classified<-merge(geedf, taxanames, by.x="variable", by.y="OTU")
#write.table(geedf_classified, file="GEE/mesal.duodenum_gee.results.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

### Finally, let's see how relatable these are to the rest of the data (all 157 samples)
# filter out only these OTUs from the data file that has all samples:
gee_otus<-as.character(geedf_classified$variable)
datadf<-data[, 15:length(data)]
data_filtered<-datadf[, colnames(datadf) %in% gee_otus]
data_filtered<-cbind(data[1:14], data_filtered)

# now run this through the model:
df<-data_filtered
start<-match("Otu0001", names(df))
end<-(length(names(df)))

catotus<-dplyr::select(df, contains("Otu"))
catotus<-catotus>0
temp<-dplyr::select(df, -contains("Otu"))
temp<-cbind(temp, catotus)
cattest<-temp
remove(temp)
remove(catotus)

testotu<-df
for(j in start:end)
{
    variable<-names(testotu)[j]
    results[j-start+1, 1]<-variable
    results[j-start+1, 2]<-summary(lmer(as.formula(paste(variable, "~", "pH + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Pr(>|t|)"]
    results[j-start+1, 3]<-summary(lmer(as.formula(paste(variable, "~", "pH + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Estimate"] 
}

temp<-as.data.frame(results)
temp$P_value<-as.numeric(as.character(temp$P_value))
temp$coefficient<-as.numeric(as.character(temp$coefficient))

results<-temp
remove(temp)

# also with mesalamine (since Otu0079 was for that):
df<-data_filtered
start<-match("Otu0001", names(df))
end<-(length(names(df)))

catotus<-dplyr::select(df, contains("Otu"))
catotus<-catotus>0
temp<-dplyr::select(df, -contains("Otu"))
temp<-cbind(temp, catotus)
cattest<-temp
remove(temp)
remove(catotus)

testotu<-df
for(j in start:end)
{
    variable<-names(testotu)[j]
    results[j-start+1, 1]<-variable
    results[j-start+1, 2]<-summary(lmer(as.formula(paste(variable, "~", "Mesalamineconcentration + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Pr(>|t|)"]
    results[j-start+1, 3]<-summary(lmer(as.formula(paste(variable, "~", "Mesalamineconcentration + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Estimate"] 
}

temp<-as.data.frame(results)
temp$P_value<-as.numeric(as.character(temp$P_value))
temp$coefficient<-as.numeric(as.character(temp$coefficient))

results<-temp
remove(temp)

### overall, all of these were significant in ALL samples, even for Otu0079
#	- these will be graphed in next section



###------------
### then for stomach samples:

# eliminate all OTUs with 0 total
rownames(data)<-data$seqID
dfs<-droplevels(data[data$site %in% c("Stomach"), c(16:ncol(data)-1)])

### now, let's filter out OTUs that occur in at least 50% of the samples BY PATIENT

## define data:
df.frac<-dfs[,  which((colSums(dfs)) > 0)]
df<-df.frac
	#convert to binary:
df[df > 0]<-1

df<-merge(data[,1:14], df, by.x="seqID", by.y="row.names")
df<-arrange(df, pID, Hour)
rownames(df)<-df$seqID

# let's set a cutoff in that we want an OTU to appear in at least 50% of the samples available per patient (so by prevalence within a patient)
# function:
listFunction <- function(n) {
	otus <- df[df$pID==(n),  13:length(df)]					# this defines your dataframe for each patient
	otus <- otus[!is.na(otus$Otu0001), ]					# this eliminates samples that are NA
	minnum<-as.numeric(nrow(otus)/2)									# this sets the minimum number of samples you want that OTU to be observed in
	otus_list<-names(which(colSums(otus == 1) >= minnum))	# this gets a list of the OTUs that appear in at least half of the patient samples
	}
# get list of patients and get lists of OTUs that appear in at least 50% of the samples:
patient.list<-as.character(unique(df$pID))
min.otu.list<-lapply(patient.list, FUN=listFunction)
names(min.otu.list)<-patient.list
lapply(min.otu.list, FUN=length)

# now, collapse this list--THIS will be the list we will use to filter out our data set:
otulist<-as.character(unlist(min.otu.list, recursive = FALSE))
otulist<-unique(otulist)

# ok, let's use this list to get the OTUs we want
df.filtered<-df.frac[, colnames(df.frac) %in% otulist]
df.filtered<-merge(data[,1:14], df.filtered, by.x="seqID", by.y="row.names")
df<-arrange(df.filtered, pID, Hour)
#write.table(df, file="mesal_stomach_filteredOTUfrac.txt", , sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

### run model for pH:
#df<-read.table(file="mesal_stomach_filteredOTUfrac.txt", header=T)

start<-match("Otu0001", names(df))
end<-(length(names(df)))

#now let's create the absence / presence of OTUs
catotus<-dplyr::select(df, contains("Otu"))
catotus<-catotus>0
temp<-dplyr::select(df, -contains("Otu"))
temp<-cbind(temp, catotus)
cattest<-temp
remove(temp)
remove(catotus)

#finally, let's sum the rows to get the total number of OTUs present
df$totalotus <- rowSums(cattest[start:end])
summary(df$totalotus)

# run code to model against all OTUs:
results<-matrix(data=NA, nrow = end-start+1, ncol = 3)
colnames(results)<-c("variable","P_value","coefficient")

# glmer(formula = post~I(as.factor(phylopam3_cluster)) + (sample_number|id), family=binomial, data=test, REML = F)

testotu<-df
for(j in start:end)
{
  tryCatch({
    variable<-names(testotu)[j]
    results[j-start+1, 1]<-variable
    results[j-start+1, 2]<-summary(lmer(as.formula(paste(variable, "~", "pH + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Pr(>|t|)"]
    results[j-start+1, 3]<-summary(lmer(as.formula(paste(variable, "~", "pH + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Estimate"] 
  }, error=function(e){}
  )
}

temp<-as.data.frame(results)
temp$P_value<-as.numeric(as.character(temp$P_value))
temp$coefficient<-as.numeric(as.character(temp$coefficient))

results<-temp
remove(temp)

taxonomy<-read.table(file="mesalamine.taxnames03.txt", header=TRUE)

sigresults<-as.data.frame(results)
sigresults$P_value<-as.numeric(as.character(sigresults$P_value))
sigresults<-sigresults[which(sigresults$P_value<.001), ]
topresults<-sigresults
topresults$abscoef<-abs(as.numeric(as.character(topresults$coefficient)))
topresults1<-arrange(topresults, desc(abscoef))
topresults1$metabolite<-"pH"

### then, run for mesalamine:
testotu<-df
for(j in start:end)
{
    variable<-names(testotu)[j]
    results[j-start+1, 1]<-variable
    results[j-start+1, 2]<-summary(lmer(as.formula(paste(variable, "~", "Mesalamineconcentration + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Pr(>|t|)"]
    results[j-start+1, 3]<-summary(lmer(as.formula(paste(variable, "~", "Mesalamineconcentration + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Estimate"] 
}

temp<-as.data.frame(results)
temp$P_value<-as.numeric(as.character(temp$P_value))
temp$coefficient<-as.numeric(as.character(temp$coefficient))

results<-temp
remove(temp)

sigresults<-as.data.frame(results)
sigresults$P_value<-as.numeric(as.character(sigresults$P_value))
sigresults<-sigresults[which(sigresults$P_value<.001), ]
topresults<-sigresults
topresults$abscoef<-abs(as.numeric(as.character(topresults$coefficient)))
topresults2<-arrange(topresults, desc(abscoef))
topresults2$metabolite<-"mesalamine"

### combine:
geedf<-rbind(topresults1, topresults2)
geedf$site<-"Stomach"

# read in taxonomic file to merge with this:
taxa<-read.table(file="mesalamine.taxnames03.txt", header=TRUE)
# add family names:
tax <- taxa$Taxonomy
tax <- gsub("\\(\\d*\\)", "", tax)
tax <- gsub(";unclassified", "", tax)
tax <- gsub("_1", "", tax)
tax <- gsub(";$", "", tax)
tax <- gsub("/.*", "", tax)
tax <- gsub(".*les;", "", tax)
tax <- gsub(";.*", "", tax)
taxa$family<-tax
ftaxa<-taxa[, c("OTU", "taxname", "phylum", "family", "Size")]
taxanames<-taxa[, c("OTU", "taxname", "phylum", "family", "Size")]
geedf_classified<-merge(geedf, taxanames, by.x="variable", by.y="OTU")
#write.table(geedf_classified, file="GEE/mesal.stomach_gee.results.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


### so, in the stomach glm does identify some correlated OTUs
# but, let's check again if these are transferable to all samples
gee_otus<-as.character(unique(geedf_classified$variable))
datadf<-data[, 15:length(data)]
data_filtered<-datadf[, colnames(datadf) %in% gee_otus]
data_filtered<-cbind(data[1:14], data_filtered)

# now run this through the model:
df<-data_filtered
start<-match("Otu0014", names(df))
end<-(length(names(df)))

catotus<-dplyr::select(df, contains("Otu"))
catotus<-catotus>0
temp<-dplyr::select(df, -contains("Otu"))
temp<-cbind(temp, catotus)
cattest<-temp
remove(temp)
remove(catotus)

testotu<-df
for(j in start:end)
{
    variable<-names(testotu)[j]
    results[j-start+1, 1]<-variable
    results[j-start+1, 2]<-summary(lmer(as.formula(paste(variable, "~", "pH + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Pr(>|t|)"]
    results[j-start+1, 3]<-summary(lmer(as.formula(paste(variable, "~", "pH + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Estimate"] 
}

temp<-as.data.frame(results)
temp$P_value<-as.numeric(as.character(temp$P_value))
temp$coefficient<-as.numeric(as.character(temp$coefficient))

results<-temp
remove(temp)
results1<-results[1:13,]
#### for this, get a new error:
#Warning messages:
#1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#  unable to evaluate scaled gradient
#2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#  Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
#3: Model failed to converge with 1 negative eigenvalue: -1.1e+02

# now for mesalamine OTUs
df<-data_filtered
start<-match("Otu0014", names(df))
end<-(length(names(df)))

catotus<-dplyr::select(df, contains("Otu"))
catotus<-catotus>0
temp<-dplyr::select(df, -contains("Otu"))
temp<-cbind(temp, catotus)
cattest<-temp
remove(temp)
remove(catotus)

testotu<-df
for(j in start:end)
{
    variable<-names(testotu)[j]
    results[j-start+1, 1]<-variable
    results[j-start+1, 2]<-summary(lmer(as.formula(paste(variable, "~", "Mesalamineconcentration + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Pr(>|t|)"]
    results[j-start+1, 3]<-summary(lmer(as.formula(paste(variable, "~", "Mesalamineconcentration + (Hour|pID)")), REML = F, data=testotu))$coefficients[2,"Estimate"] 
}

temp<-as.data.frame(results)
temp$P_value<-as.numeric(as.character(temp$P_value))
temp$coefficient<-as.numeric(as.character(temp$coefficient))

results<-temp
remove(temp)
results2<-results
#### did not get same warning as before
# BUT, none of the OTUs were significant when all patient-sites are considered
# perhaps just mention these in a supplemental table, and then mention in text?


### Table S2: combined the stomach/duodenum lists together


```

####### Step 2A: visualizing / graphing results
	- Figure S5:  relative abundance vs. pH or mesalamine for all significantly correlated OTUs, duodenum
	- Figure S6:  relative abundance vs. pH for all significantly correlated OTUs, stomach
	- Figure S7:  relative abundance vs. mesalamine for all significantly correlated OTUs, stomach
	
```

###-----------------###
#### First, let's graph the RA vs pH or mesalamine for all OTUs (in duodenum and stomach)

### Figure S5:

# read in data:
gee<-read.table(file="GEE/mesal.duodenum_gee.results.txt", header=T)
data<-read.table(file="mesal_otufrac.w.meta.txt", header=TRUE)

# get only duodenum samples:
data<-data[data$seq_status %in% c("OK"), ]
summary(data$site)		# need to combine stomach/antrum
data$site[data$site %in% c("Antrum")]<-"Stomach"
	# note: this already has crappy samples removed
rownames(data)<-data$seqID
dfs<-droplevels(data[data$site %in% c("Duodenum"), ])

# let's rename the Otu column names so they list the taxonomy:
taxa<-read.table(file="mesalamine.taxnames03.txt", header=TRUE)
otus<-colnames(dfs[, 15:length(dfs)-1])
taxanames<-as.character(taxa[taxa$OTU %in% otus, c("taxname")])
names(dfs)[15:length(dfs)-1]<-taxanames

# filter out OTUs of interest:
# g1: pH correlations
gee<-gee[order(gee$metabolite, -gee$abscoef),]
ph.otus<-gee[gee$metabolite %in% c("pH") & gee$site %in% c("Duodenum"), c("variable")]
ph.otus.name<-as.character(gee[gee$metabolite %in% c("pH") & gee$site %in% c("Duodenum"), c("taxname")])

### let's plot all the Otus that are relevant:

data_shared<-dfs
plotCorrs <- function(y, x) {
	plot(data_shared[,x], data_shared[,y], xlab=names(data_shared[x]), ylab=names(data_shared[y]), tck=-0.03, bty="n",
		xaxt="n", yaxt="n"#, ylim=c(0,100)
		)
	axis(side=2,tck=-0.03, line=0, mgp=c(2,0.5,0))
	axis(side=1,tck=-0.03, line=0, mgp=c(2,0.5,0))
	corrtest<-cor.test(data_shared[,x], data_shared[,y])
	print(corrtest)
	#mtext(ph.otus.name[y], 2, line=0)
	text(max(data_shared[,x])*.4, max(data_shared[,y]), labels=paste0("p = ", signif(corrtest$p.value, digits=3)), cex=0.6)
	text(max(data_shared[,x])*.4, max(data_shared[,y])*.93, labels=paste0("corr =  ", signif(corrtest$estimate, digits=2)), cex=0.6)
	abline(lm(data_shared[,y]~data_shared[,x]), col="red")
}

#par(mfrow=c(round(nrow(sig.corrs)/2), 2))
#otulist<-as.character(pos.corrs$var2)
# OR
#par(mfrow=c(round(nrow(sig.corrs)/3), 3), mar=c(3,4,1,0))
par(mfrow=c(4,4), mar=c(3,4,2,0), mgp=c(2,1,0.5))
lapply(ph.otus.name, FUN=plotCorrs, x=c("pH"))


### also graph for mesalamine:
mes.otus<-as.character(gee[gee$metabolite %in% c("mesalamine") & gee$site %in% c("Duodenum"), c("taxname")])
lapply(mes.otus, FUN=plotCorrs, x=c("Mesalamineconcentration"))

### Figure S6:

# read in data:
gee<-read.table(file="GEE/mesal.stomach_gee.results.txt", header=T)
data<-read.table(file="mesal_otufrac.w.meta.txt", header=TRUE)

# get only duodenum samples:
data<-data[data$seq_status %in% c("OK"), ]
summary(data$site)		# need to combine stomach/antrum
data$site[data$site %in% c("Antrum")]<-"Stomach"
	# note: this already has crappy samples removed
rownames(data)<-data$seqID
dfs<-droplevels(data[data$site %in% c("Stomach"), ])

# let's rename the Otu column names so they list the taxonomy:
taxa<-read.table(file="mesalamine.taxnames03.txt", header=TRUE)
otus<-colnames(dfs[, 15:length(dfs)-1])
taxanames<-as.character(taxa[taxa$OTU %in% otus, c("taxname")])
names(dfs)[15:length(dfs)-1]<-taxanames

# filter out OTUs of interest:
# g1: pH correlations
gee<-gee[order(gee$metabolite, -gee$abscoef),]
ph.otus<-gee[gee$metabolite %in% c("pH") & gee$site %in% c("Stomach"), c("variable")]
ph.otus.name<-as.character(gee[gee$metabolite %in% c("pH") & gee$site %in% c("Stomach"), c("taxname")])

### let's plot all the Otus that are relevant:

data_shared<-dfs
plotCorrs <- function(y, x) {
	plot(data_shared[,x], data_shared[,y], xlab=names(data_shared[x]), ylab=names(data_shared[y]), tck=-0.03, bty="n",
		xaxt="n", yaxt="n"#, ylim=c(0,100)
		)
	axis(side=2,tck=-0.03, line=0, mgp=c(2,0.5,0))
	axis(side=1,tck=-0.03, line=0, mgp=c(2,0.5,0))
	corrtest<-cor.test(data_shared[,x], data_shared[,y])
	print(corrtest)
	#mtext(ph.otus.name[y], 2, line=0)
	text(max(data_shared[,x])*.4, max(data_shared[,y]), labels=paste0("p = ", signif(corrtest$p.value, digits=3)), cex=0.6)
	text(max(data_shared[,x])*.4, max(data_shared[,y])*.93, labels=paste0("corr =  ", signif(corrtest$estimate, digits=2)), cex=0.6)
	abline(lm(data_shared[,y]~data_shared[,x]), col="red")
}

#par(mfrow=c(round(nrow(sig.corrs)/2), 2))
#otulist<-as.character(pos.corrs$var2)
# OR
#par(mfrow=c(round(nrow(sig.corrs)/3), 3), mar=c(3,4,1,0))
par(mfrow=c(3,3), mar=c(3,4,2,0), mgp=c(2,1,0.5))
lapply(ph.otus.name, FUN=plotCorrs, x=c("pH"))

### also graph for mesalamine:

# Figure S7:

mes.otus<-as.character(gee[gee$metabolite %in% c("mesalamine") & gee$site %in% c("Stomach"), c("taxname")])
par(mfrow=c(3,3), mar=c(3,4,2,0), mgp=c(2,1,0.5))
lapply(mes.otus, FUN=plotCorrs, x=c("Mesalamineconcentration"))


```

####### Step 2B: visualizing / graphing results
	- Figure 5:  pH / mesalamine for duodenal samples
	- Figure S8:  pH / mesalamine for stomach samples
	
```

###-----------------###
#### let's make a heatmap of the correlated OTUs

### Figure 5: OTUs correlated with pH in duodenum

library(RColorBrewer)
library(gplots)
library(vegan)
library(plyr)

# data:
gee<-read.table(file="GEE/mesal.duodenum_gee.results.txt", header=T)
data<-read.table(file="mesal_otufrac.w.meta.txt", header=TRUE)

# get only duodenal samples:
data<-data[data$seq_status %in% c("OK"), ]
summary(data$site)		# need to combine stomach/antrum
data$site[data$site %in% c("Antrum")]<-"Stomach"
	# note: this already has crappy samples removed
rownames(data)<-data$seqID
dfs<-droplevels(data[data$site %in% c("Duodenum"), ])

# order by pH:
rownames(dfs)<-dfs$seqID
df<-dfs[order(dfs$pH), ]

# now, get only OTUs that are in your GEE list:
rownames(df)<-df$seqID
df.matrix<-df[, 13:ncol(df)]
df.matrix<-df.matrix[, colnames(df.matrix) %in% as.character(gee$variable)]

# also convert to log10(1+RA)
df.matrix<-log10(df.matrix+1)

# set color:
pH.cols<-colorRampPalette(c("yellow","orange","red", "darkred"))(n=nrow(df))

hm.cols <- colorRampPalette(c("white", "lightblue", "dodgerblue1", "dodgerblue3", "blue", "darkblue"))

heatmap.2(as.matrix(df.matrix), 
			col=hm.cols(15), 
			#breaks=my.breaks, 
			trace="none", 
			density.info="none", 
			Colv=F, 			
			Rowv=F,
			dendrogram="none",
			RowSideColors = pH.cols
			)
# doesn't look bad... let's also sort by OTUs
gee<-gee[order(gee$Size, gee$phylum, gee$family), ]

tax.col<-gee[,c("variable", "phylum", "Size", "taxname")] #take the column with the group
summary(tax.col$phylum)			# see what levels are available
tax.col$phylum <- factor(tax.col$phylum, levels = c("Bacteroidetes","Firmicutes","Proteobacteria", "Actinobacteria", "Fusobacteria"))
tax.col$color<-mapvalues(tax.col$phylum, from = c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Fusobacteria", "Proteobacteria"), to = c("red4", "green4", "skyblue", "red", "yellow2"))
tax.col<- tax.col[order(tax.col$phylum, -tax.col$Size) , ]

# order heatmap according to topotus, then bind to column variable:
col.order<-as.character(tax.col[, 1])								#convert the order of the OTUs into a list
df.matrix<-df.matrix[,col.order]				#order your matrix by the ordered list
rbind(colnames(df.matrix), tax.col)									#check order of the matrix and list to ensure that colors will be correct (they should match!)

# check if this works:

heatmap.2(as.matrix(df.matrix), 
			col=hm.cols(15), 
			#breaks=my.breaks, 
			trace="none", 
			density.info="none", 
			Colv=F, 			
			Rowv=F,
			dendrogram="none",
			RowSideColors = pH.cols,
			ColSideColors = as.character(tax.col$color)
			)
			
# final heatmap for Figure 5:

heatmap.2(as.matrix(df.matrix),
                    notecol="black",     
                    density.info="none",
                    key.xlab="",
                    key.title="",  
                    trace="none",         
                    margins =c(2,4),    
                    col=hm.cols(15),        
                    #breaks=my.breaks,    
                    key=TRUE,  
                  	dendrogram="none",    
                    Colv=F,
                    Rowv=F,
                    srtCol=45,
                    cexRow= 0.7,
                    cexCol = 0.8,
                    keysize=1,
                    lwid = c(3,5),
                    lhei = c(2,5),
                    labCol=gsub("_", " ", as.character(tax.col$taxname)),
                    labRow=as.character(df$pH),
					RowSideColors = pH.cols,
					ColSideColors = as.character(tax.col$color)
					)
					
```
					
####### Step 2C: visualizing / graphing results
	- Figure 5 (other option):  pH / mesalamine for duodenal samples
	
```

###-----------------###
#### let's try to visualize the OTU RAs vs. pH in ONE graph

```
### read in data:
gee<-read.table(file="GEE/mesal.duodenum_gee.results.txt", header=T)
data<-read.table(file="mesal_otufrac.w.meta.txt", header=TRUE)

# get only duodenum samples:
data<-data[data$seq_status %in% c("OK"), ]
summary(data$site)		# need to combine stomach/antrum
data$site[data$site %in% c("Antrum")]<-"Stomach"
	# note: this already has crappy samples removed
rownames(data)<-data$seqID
dfs<-droplevels(data[data$site %in% c("Duodenum"), ])

# get only RA's that are in gee list:
rownames(dfs)<-dfs$seqID
df.matrix<-dfs[, 14:ncol(dfs)]
df.matrix<-df.matrix[, colnames(df.matrix) %in% as.character(gee$variable)]

# also convert to log10(1+RA)
df.matrix<-log10(df.matrix+1)

# combine with meta, again:
dfs<-cbind(dfs[, 1:13], df.matrix)

# let's rename the Otu column names so they list the taxonomy:
names(dfs)[14:length(dfs)]<-as.character(gee$taxname)

# let's also assign a color according to taxonomy:
tax.col<-gee[,c("variable", "phylum", "Size", "taxname", "coefficient")] #take the column with the group
summary(tax.col$phylum)			# see what levels are available
tax.col$phylum <- factor(tax.col$phylum, levels = c("Bacteroidetes","Firmicutes","Proteobacteria", "Actinobacteria", "Fusobacteria"))
tax.col$color<-mapvalues(tax.col$phylum, from = c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Fusobacteria", "Proteobacteria"), to = c("red4", "green4", "skyblue", "red", "yellow2"))
tax.col<- tax.col[order(tax.col$phylum, -tax.col$Size) , ]

### first panel: all Firmicutes:

# get list of OTUs that fall under Firmicutes:
flist<-tax.col[tax.col$phylum=="Firmicutes", c("taxname")]

# assign a new color:
#firm<-colorRampPalette(c("midnightblue","mediumblue","blue3","blue","dodgerblue4","dodgerblue1","deepskyblue4","deepskyblue1","skyblue3","skyblue","steelblue4","steelblue1","royalblue4","royalblue1","slateblue4","purple3","orchid3","plum4","plum1","pink3","pink","lightpink1","lightpink3","palevioletred4","palevioletred1","magenta4","deeppink4","mediumvioletred","magenta3","magenta1","thistle"))(n=25)
firm.col<-colorRampPalette(c("midnightblue","dodgerblue4","deepskyblue1","orchid3", "pink"))(n=length(flist))

# then, plot onto the same graph:

par(mfrow=c(1,3))

# singular:
data_shared<-dfs
plot(type='n', data_shared$pH, data_shared$Otu0001_Streptococcus, xlab="pH", ylab="relative abundance (%)", tck=-0.03, bty="n", ylim=c(0,2), xlim=c(1,8)
		)
points(data_shared$pH, data_shared$Otu0001_Streptococcus, col=firm.col[1], bg=adjustcolor(firm.col[1], alpha=0.5), pch=21)
abline(lm(data_shared$Otu0001_Streptococcus~data_shared$pH), col=firm.col[1], lwd=3)
points(data_shared$pH, data_shared$Otu0005_Gemella, col=firm.col[2], bg=adjustcolor(firm.col[2], alpha=0.5), pch=21)
abline(lm(data_shared$Otu0005_Gemella~data_shared$pH), col=firm.col[2], lwd=3)
points(data_shared$pH, data_shared$Otu0007_Streptococcus, col=firm.col[3], bg=adjustcolor(firm.col[3], alpha=0.5), pch=21)
abline(lm(data_shared$Otu0007_Streptococcus~data_shared$pH), col=firm.col[3], lwd=3)
points(data_shared$pH, data_shared$Otu0008_Lactobacillales_unclassified, col=firm.col[4], bg=adjustcolor(firm.col[4], alpha=0.5), pch=21)
abline(lm(data_shared$Otu0008_Lactobacillales_unclassified~data_shared$pH), col=firm.col[4], lwd=3)
points(data_shared$pH, data_shared$Otu0194_Streptococcus, col=firm.col[5], bg=adjustcolor(firm.col[5], alpha=0.5), pch=21)
abline(lm(data_shared$Otu0194_Streptococcus~data_shared$pH), col=firm.col[5], lwd=3)
legend("top", legend=as.character(gsub("_", " ", flist)), col=firm.col, pt.bg=adjustcolor(firm.col, alpha=0.5), pch=21, cex=0.7)

# all Bacteroidetes:
blist<-tax.col[tax.col$phylum=="Bacteroidetes", c("taxname")]

#bac<-colorRampPalette(c("darkgreen","green3","lightgreen","seagreen", "lightgreen"))(n=10)
bac.col<-colorRampPalette(c("darkgreen","green3","lightgreen","seagreen", "lightgreen"))(n=length(blist))

data_shared<-dfs
plot(type='n', data_shared$pH, data_shared$Otu0001_Streptococcus, xlab="pH", ylab="relative abundance (%)", tck=-0.03, bty="n", ylim=c(0,2), xlim=c(1,8)
		)
points(data_shared$pH, data_shared$Otu0011_Prevotellaceae_unclassified, col=bac.col[1], bg=adjustcolor(bac.col[1], alpha=0.5), pch=21)
abline(lm(data_shared$Otu0011_Prevotellaceae_unclassified~data_shared$pH), col=bac.col[1], lwd=3)
points(data_shared$pH, data_shared$Otu0022_Prevotella, col=bac.col[2], bg=adjustcolor(bac.col[2], alpha=0.5), pch=21)
abline(lm(data_shared$Otu0022_Prevotella~data_shared$pH), col=bac.col[2], lwd=3)
points(data_shared$pH, data_shared$Otu0025_Porphyromonas, col=bac.col[3], bg=adjustcolor(bac.col[3], alpha=0.5), pch=21)
abline(lm(data_shared$Otu0025_Porphyromonas~data_shared$pH), col=bac.col[3], lwd=3)
points(data_shared$pH, data_shared$Otu0030_Prevotella, col=bac.col[4], bg=adjustcolor(bac.col[4], alpha=0.5), pch=21)
abline(lm(data_shared$Otu0030_Prevotella~data_shared$pH), col=bac.col[4], lwd=3)
points(data_shared$pH, data_shared$Otu0033_Prevotella, col=bac.col[5], bg=adjustcolor(bac.col[5], alpha=0.5), pch=21)
abline(lm(data_shared$Otu0033_Prevotella~data_shared$pH), col=bac.col[5], lwd=3)
points(data_shared$pH, data_shared$Otu0082_Flavobacteriaceae_unclassified, col=bac.col[6], bg=adjustcolor(bac.col[6], alpha=0.5), pch=21)
abline(lm(data_shared$Otu0082_Flavobacteriaceae_unclassified~data_shared$pH), col=bac.col[6], lwd=3)
legend("top", legend=as.character(gsub("_", " ", blist)), col=bac.col, pt.bg=adjustcolor(bac.col, alpha=0.5), pch=21, cex=0.7)

# all other:
olist<-tax.col[tax.col$phylum %in% c("Proteobacteria", "Actinobacteria", "Fusobacteria"), c("taxname")]

#colorRampPalette(c("yellow2","darkgoldenrod3","goldenrod2","orange2","yellow4"))(n=12)
pro.col<-colorRampPalette(c("yellow2","darkgoldenrod3","goldenrod2","orange2","yellow4"))(n=(3))
other.col<-c(pro.col, "tan", "red")

data_shared<-dfs
plot(type='n', data_shared$pH, data_shared$Otu0001_Streptococcus, xlab="pH", ylab="relative abundance (%)", tck=-0.03, bty="n", ylim=c(0,2), xlim=c(1,8)
		)
points(data_shared$pH, data_shared$Otu0003_Pasteurellaceae_unclassified, col=other.col[1], bg=adjustcolor(other.col[1], alpha=0.5), pch=21)
abline(lm(data_shared$Otu0003_Pasteurellaceae_unclassified~data_shared$pH), col=other.col[1], lwd=3)
points(data_shared$pH, data_shared$Otu0009_Haemophilus, col=other.col[2], bg=adjustcolor(other.col[2], alpha=0.5), pch=21)
abline(lm(data_shared$Otu0009_Haemophilus~data_shared$pH), col=other.col[2], lwd=3)
points(data_shared$pH, data_shared$Otu0054_Pasteurellaceae_unclassified, col=other.col[3], bg=adjustcolor(other.col[3], alpha=0.5), pch=21)
abline(lm(data_shared$Otu0054_Pasteurellaceae_unclassified~data_shared$pH), col=other.col[3], lwd=3)
#points(data_shared$pH, data_shared$Otu0010_Actinomyces, col=other.col[4], bg=adjustcolor(other.col[4], alpha=0.5), pch=21)
#abline(lm(data_shared$Otu0010_Actinomyces~data_shared$pH), col=other.col[4], lwd=3)
#points(data_shared$pH, data_shared$Otu0079_Leptotrichiaceae_unclassified, col=other.col[5], bg=adjustcolor(other.col[5], alpha=0.5), pch=21)
#abline(lm(data_shared$Otu0079_Leptotrichiaceae_unclassified~data_shared$pH), col=other.col[5], lwd=3)
legend("top", legend=as.character(gsub("_", " ", olist[1:3])), col=other.col[1:3], pt.bg=adjustcolor(other.col[1:3], alpha=0.5), pch=21, cex=0.7)


# if separating actinomyces:
data_shared<-dfs
plot(type='n', data_shared$pH, data_shared$Otu0001_Streptococcus, xlab="pH", ylab="relative abundance (%)", tck=-0.03, bty="n", ylim=c(0,2), xlim=c(1,8)
		)
points(data_shared$pH, data_shared$Otu0010_Actinomyces, col=other.col[4], bg=adjustcolor(other.col[4], alpha=0.5), pch=21)
abline(lm(data_shared$Otu0010_Actinomyces~data_shared$pH), col=other.col[4], lwd=3)
legend("top", legend=as.character(gsub("_", " ", olist[4])), col=other.col[4], pt.bg=adjustcolor(other.col[4], alpha=0.5), pch=21, cex=0.7)




# now, do for all firmicutes:
# could not quite get it to work...for now
data_shared<-dfs
plotCorrs.SamePlot <- function(y, x, z) {
	#plot(type='n', data_shared[,x], data_shared[,y], xlab=names(data_shared[x]), 
	#	ylab="relative abundance (%)", tck=-0.03, bty="n", ylim=c(0,100)
	#	)
	points(data_shared[,x], data_shared[,y], col=firm.col[z], bg=adjustcolor(firm.col[z], alpha=0.5), pch=21)
	#axis(side=2,tck=-0.03, line=0, mgp=c(2,0.5,0))
	#axis(side=1,tck=-0.03, line=0, mgp=c(2,0.5,0))
	#corrtest<-cor.test(data_shared[,x], data_shared[,y])
	#print(corrtest)
	#mtext(ph.otus.name[y], 2, line=0)
	#text(max(data_shared[,x])*.4, max(data_shared[,y]), labels=paste0("p = ", signif(corrtest$p.value, digits=3)), cex=0.6)
	#text(max(data_shared[,x])*.4, max(data_shared[,y])*.93, labels=paste0("corr =  ", signif(corrtest$estimate, digits=2)), cex=0.6)
	abline(lm(data_shared[,y]~data_shared[,x]), col=firm.col[z], lwd=3)
}
lapply(flist, FUN=plotCorrs.SamePlot, x=c("pH"), z=firm.col)


```


