#### Mesalamine Study Manuscript:  Figure 3, similarity and "core" OTUs
##### Anna M. Seekatz
##### 08.22.18

###### Figures created:
	- Fig. 4: shared OTU % at each site by patient (based on presence / absence)
		- A: % OTUs at a given subject-site out of total detected OTUs in that individual's site
		- B: % OTUs of that were detected in all subject-site samples ("core" OTUs) out of total OTUs at a given time point
		- C: % relative abundance that is occupied by the detected "core" OTUs


####### Files used:
	- mesal_subs.otucounts.w.meta.txt
	- mesal_otufrac.w.meta.txt


```
### get samples:
rawotus<-read.table(file="mesal_subs.otucounts.w.meta.txt", header=T, row.names=1)
	# note: this already has samples with < 4000 removed
	
# let's also add a new column that specifies the organ AND patient
rawotus$pID_site<-paste(rawotus$pID, rawotus$location, sep="_")
summary(as.factor(rawotus$pID_site))		
	# some sites only have 1 sample collected--these are not useful to us for this analysis
	# let's remove those samples
multiples.list<-names(which(summary(as.factor(rawotus$pID_site)) > 1))
data<-droplevels(rawotus[rawotus$pID_site %in% multiples.list, ])
	
# get only OTUs and sampleID, patientID:
df<-data[, c(1, 603, 13:length(data)-2)]		#for now, removing the "rareOTUs" column (will need to add this to final count eventually...
rownames(df)<-data$seqID
df$pID_site<-as.factor(df$pID_site)

# convert to 1/0 (presence/absence)
df[df > 0]<-1	

### Let's calculate how many OTUs are:
	- detected in all samples taken from a specific site per subject (our denominator)
	- detected at a given timepoint out of all the OTUs possible for that subject-site
	- detected in all samples across time within a subject site ("core" OTUs)
### also, we will calculate the % of relative abundance represented by the core OTUs
	

## 1) otu lists by patient-site:
# function:
listFunction <- function(n) {
	otus <- df[df$pID_site==(n),  5:length(df)]
	otus <- otus[!is.na(otus$Otu0001), ]
	min <- nrow(otus)
	otus_list<-names(which(colSums(otus == 1) > 0))
	}
# get list of patients:
patient.list<-as.character(unique(df$pID_site))
p.otu.list<-lapply(patient.list, FUN=listFunction)
names(p.otu.list)<-patient.list
otu.totals.across.sites<-lapply(p.otu.list, FUN=length)
	# this gives you the # of unique OTUs that ever show up in a patient (note: excluding the "rare" OTUs column)
	# might need this, but not for the next section
### mean number of OTUs ever detected in a subject-site:
mean(unlist(otu.totals.across.sites))	# 135.9
min(unlist(otu.totals.across.sites))		# 78
max(unlist(otu.totals.across.sites))		# 212

# could also calculate which OTUs are "rare" within a patient (never show up in anyone else):
unique.p.otus<-lapply(1:length(p.otu.list), function(n) setdiff(p.otu.list[[n]], unlist(p.otu.list[-n])))
	# might need this, but not for the next section


## 2) Also need to calculate which OTUs are unique across the samples:
# the previous lists specify the potential OTUs found in each patient
# need to also see which OTUs are found across the conditions
# so first, we will create a function that subsets the list of OTUs into a list of the appropriate OTUs

# to do this per patient, need to loop the following:
test<-df[, colnames(df) %in% p.otu.list[[1]]]

# loop this to get a list of specific dataframes:
df.list<-lapply(1:length(patient.list), function (n) df[, colnames(df) %in% p.otu.list[[n]]])

# let's also loop to get ONLY the samples that we want per pID, and apply this to the df as well:
sample.listFunction <- function(n) {
	otus <- df[df$pID_site==(n),  c("seqID")]
	}
sample.list<-lapply(patient.list, FUN=sample.listFunction)

# now, pull only those samples from the lists:
df.list2<-lapply(1:length(patient.list), function (n) df[rownames(df) %in% sample.list[[n]], colnames(df) %in% p.otu.list[[n]]])
	# this is now a list of all the samples found in each patient
	
	
### 3) now, we need to calculate % of the "core" mb found in each patient's dataframe
# you can do this by counting the number observed (OTU at == 1) divided by the total number of columns in that dataframe
# this will give you the % stable across patient samples

# example:
rowSums(df.list2[[2]])		#gives you the total # of OTUs found in that sample
rowSums(df.list2[[2]]/length(df.list2[[2]]))		# the percent of OTUs ever observed in that patient observed at that moment

percent.otus.observed<-lapply(1:length(patient.list), function(n) rowSums(df.list2[[n]])/length(df.list2[[n]])*100)
	# this is a list of the percent of 'total observed' OTUs that are observed across the patient at any given point (by sample)
mean(unlist(percent.otus.observed))	# 36.63
min(unlist(percent.otus.observed))		# 4.48
max(unlist(percent.otus.observed))		# 67.62

######## this is a smaller number than I thought, showing that overall, there is a lot of flux going on
# what about PERCENT of the community that these OTUs describe at any given moment? should we look at that?
# or, what about OTUs that are ALWAYS present?
# how about we look at the 'core' OTUs (so which OTUS are present all the time), then see overall how much of the community this core OTU represents

# let's create a list of the OTUs per patient-site that are present in ALL of the samples for that patient-site:
# this creates a dataframe of those OTUs
core.otu.df<-lapply(1:length(patient.list), function(n) df.list2[[n]][ , which(colSums(df.list2[[n]]) == nrow(df.list2[[n]]))])
	# however, this is a dataframe--we just need a list:
core.otu.list<-lapply(1:length(patient.list), function(n) names(df.list2[[n]][ , which(colSums(df.list2[[n]]) == nrow(df.list2[[n]]))]))

## for number of core OTUs per patient-site:
core.totals.across.sites<-lapply(core.otu.list, FUN=length)
mean(unlist(core.totals.across.sites))	# 14.14
min(unlist(core.totals.across.sites))		# 2
max(unlist(core.totals.across.sites))	# 45

percent.core.otus.total<-lapply(1:length(patient.list), function(n) length(df.list2[[n]][ , which(colSums(df.list2[[n]]) == nrow(df.list2[[n]]))]) / length(df.list2[[n]]) *100)
names(percent.core.otus.total)<-patient.list
	# this reflects the % of total OTUs that are observed ALL the time within a patient-site
	# note: there is only 1 number for this list per patient-site (this is not by sample...)
percent.core.otus<-lapply(1:length(patient.list), function(n) length(df.list2[[n]][ , which(colSums(df.list2[[n]]) == nrow(df.list2[[n]]))]) / rowSums(df.list2[[n]]) *100)
	# this is a BY SAMPLE list: it reflects the % of each sample that the core is observed in (still at only presence/absence)
## what is mean across the board:
mean(unlist(percent.core.otus.total))		# 10.15%
min(unlist(percent.core.otus.total))		# 1.49%
max(unlist(percent.core.otus.total))		# 25%
## what about mean total # of core OTUs across the board:
core.otu.totals<-lapply(1:length(patient.list), function(n) length(df.list2[[n]][ , which(colSums(df.list2[[n]]) == nrow(df.list2[[n]]))]))
mean(unlist(core.otu.totals))		# 14.14
min(unlist(core.otu.totals))		# 2
max(unlist(core.otu.totals))		# 45


######## finally, we can use the list of core OTUs to identify how much ABUNDANCE is taken up by the core (OTUs that we observe all the time)
# what about PERCENT of the community that these OTUs describe at any given moment? should we look at that?

# read in relabund file and format as before:
relabund<-read.table(file="mesal_otufrac.w.meta.txt", header=T, row.names=1)
relabund$site[relabund$site=="Antrum"]<-"Stomach"
relabund<-droplevels(relabund)
relabund$site <- factor(relabund$site, levels = c("Stomach", "Duodenum", "Proximal_Jejunum", "Mid_Jejunum", "Distal_Jejunum", "Stool"))
relabund$pID_site<-paste(relabund$pID, rawotus$location, sep="_")
summary(as.factor(relabund$pID_site))		
multiples.list<-names(which(summary(as.factor(relabund$pID_site)) > 1))
ra.data<-droplevels(relabund[relabund$pID_site %in% multiples.list, ])
	
# get only OTUs and sampleID, patientID:
ra.df<-ra.data[, c(1, 604, 14:length(ra.data)-2)]		#for now, removing the "rareOTUs" column (will need to add this to final count eventually...
rownames(ra.df)<-ra.data$seqID
ra.df$pID_site<-as.factor(ra.df$pID_site)

# the relabund file has been calculated as relative abundance of ALL otus
# we want to get relative abundance of ONLY core otus...in a list form, again.
# we will use the core.otu.list file we just created, and filter our ra.df multiple times to get it:
ra.df.list<-lapply(1:length(patient.list), function (n) ra.df[rownames(ra.df) %in% sample.list[[n]], colnames(ra.df) %in% core.otu.list[[n]]])
names(ra.df.list)<-patient.list
	# this is a list of the % relative abundance of ONLY the core OTUs within an individual
# you can also get a list of all the OTUs, separated by patient-site
all.ra.df.list<-lapply(1:length(patient.list), function (n) ra.df[rownames(ra.df) %in% sample.list[[n]], ])
names(all.ra.df.list)<-patient.list

# now, we can calculate the totals per sample that make up the relative abundance of each sample:
percent.ra.core.otus<-lapply(1:length(patient.list), function(n) rowSums(ra.df.list[[n]]) )
names(percent.ra.core.otus)<-patient.list
### what is the mean % relative abundance occupied by these:
mean(unlist(percent.ra.core.otus))		# 71.995%
min(unlist(percent.ra.core.otus))		# 0.0667%
max(unlist(percent.ra.core.otus))		# 99.199%

#### we have 2 values now:
	- percent.otus.observed:  % of total OTUs ever observed in that person per site within a particular sample
	- percent.core.otus:  % of the OTUs in that sample that are core OTUs (for that patient-site)
	- percent.ra.core.otus:  % of relative abundance that is covered by CORE OTUs (OTUs that are always present in patient-site)
	
# combine these:
alldf<-do.call(rbind, Map(data.frame, p_otus_observed=percent.otus.observed, pn_core_otus=percent.core.otus, pra_core_otus=percent.ra.core.otus))

# combine with metadata, and graph:
f<-merge(relabund[,c(1:13, 604)], alldf, by.x="seqID", by.y="row.names")
#write.table(f, file="mesal_OTU.percent.totals", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

div.boxplot<- function (n, m) {
				plot(data[,n] ~ as.factor(data[,m]), data = data, boxwex=0.7, ylab=names(data[n]), 
					xlab="n", xaxt="n", outline=FALSE, mpg=c(2,1,0), ylim=c(0,100), cex=1, cex.lab=1, cex.axis=1)
				points(data[,n] ~ jitter(as.numeric(as.factor(data[,m]), factor=0)), data = data, col=adjustcolor("grey45", alpha=0.5), pch=21, bg="grey80", cex=0.7)
				print(as.character(unique(data[,m])))
				names<-as.character(levels(data[,m]))
				text(x =  seq(1,length(unique(data[,m])),by=1), y = par("usr")[3]-max(data[n], na.rm=TRUE)/12, srt = 45, adj = 1, labels = paste0(names, "\n n = ", as.numeric(summary(data$site))), xpd = TRUE, cex=0.8)
				stat<-kruskal.test(data[,n]~data[,m], data)
				print(stat$p.value)
				text(length(unique(data[,m]))/2, max(data[,n], na.rm=TRUE), labels=paste0("Kruskal-Wallis test \np < ", signif(stat$p.value, digits=3)), cex=0.5)
				#mtext(names(data[m]), 1, line=2)
				}
par(mar=c(5,3,2,5))
data<-droplevels(f)
par(mfrow=c(3,1))
lapply(c("p_otus_observed"), FUN=div.boxplot, m=c("site"))
title('% n OTUs detected from total OTUs detected at patient site')
lapply(c("pn_core_otus"), FUN=div.boxplot, m=c("site"))
title('% n core OTUs detected at patient site\n(core OTU = OTU always detected at patient site)')
lapply(c("pra_core_otus"), FUN=div.boxplot, m=c("site"))
title('% relative abundance of core OTUs detected at patient site')


```

### also did the same calculations to individual who was admitted 3X:

```
### subset samples only from this individual:
rawotus<-read.table(file="mesal_subs.otucounts.w.meta.txt", header=T, row.names=1)
	# note: this already has samples with < 4000 removed
	
# let's also add a new column that specifies the organ AND patient
# use subject instead of pID, which collates all samples to that individual
rawotus$sub_site<-paste(rawotus$subject, rawotus$location, sep="_")
summary(as.factor(rawotus$sub_site))		
	# now, only 4 totals for this individual:
	# M046_Duodenum  M046_Jejunum  M046_Stomach    M046_Stool
# subset to only M046:
rawotus<-rawotus[rawotus$subject %in% c("M046"), ]
multiples.list<-names(which(summary(as.factor(rawotus$sub_site)) > 2))		# this will get rid of stool...
data<-droplevels(rawotus[rawotus$sub_site %in% multiples.list, ])
	
# get only OTUs and sampleID, patientID:
df<-data[, c(1, 603, 13:length(data)-2)]		#for now, removing the "rareOTUs" column (will need to add this to final count eventually...
rownames(df)<-data$seqID
df$sub_site<-as.factor(df$sub_site)

# convert to 1/0 (presence/absence)
df[df > 0]<-1	

### now, calculate which OTUs are present at each site, again:
## 1) otu lists by patient-site:
# function:
listFunction <- function(n) {
	otus <- df[df$sub_site==(n),  5:length(df)]
	otus <- otus[!is.na(otus$Otu0001), ]
	min <- nrow(otus)
	otus_list<-names(which(colSums(otus == 1) > 0))
	}
# get list of patients:
patient.list<-as.character(unique(df$sub_site))
p.otu.list<-lapply(patient.list, FUN=listFunction)
names(p.otu.list)<-patient.list
otu.totals.across.sites<-lapply(p.otu.list, FUN=length)
# $M046_Duodenum [1] 28
# $M046_Jejunum [1] 29
# $M046_Stomach [1] 78

## 2) Also need to calculate which OTUs are unique across the samples:
# the previous lists specify the potential OTUs found in each patient
# need to also see which OTUs are found across the conditions
# so first, we will create a function that subsets the list of OTUs into a list of the appropriate OTUs

# to do this per patient, need to loop the following:
test<-df[, colnames(df) %in% p.otu.list[[1]]]

# loop this to get a list of specific dataframes:
df.list<-lapply(1:length(patient.list), function (n) df[, colnames(df) %in% p.otu.list[[n]]])

# let's also loop to get ONLY the samples that we want per pID, and apply this to the df as well:
sample.listFunction <- function(n) {
	otus <- df[df$sub_site==(n),  c("seqID")]
	}
sample.list<-lapply(patient.list, FUN=sample.listFunction)

# now, pull only those samples from the lists:
df.list2<-lapply(1:length(patient.list), function (n) df[rownames(df) %in% sample.list[[n]], colnames(df) %in% p.otu.list[[n]]])
	# this is now a list of all the OTUs found in all the samples per patient-site


### 3) now, we need to calculate % of the "core" mb found in this patient's dataframe

# example:
rowSums(df.list2[[2]])		#gives you the total # of OTUs found in that sample
rowSums(df.list2[[2]]/length(df.list2[[2]]))		# the percent of OTUs ever observed in that patient observed at that moment

percent.otus.observed<-lapply(1:length(patient.list), function(n) rowSums(df.list2[[n]])/length(df.list2[[n]])*100)
	# this is a list of the percent of 'total observed' OTUs that are observed across the patient at any given point (by sample)
mean(unlist(percent.otus.observed))		# 22.4%
min(unlist(percent.otus.observed))		# 9.25%
max(unlist(percent.otus.observed))		# 50.9%
	# this is across all bodysites--meaning that overall, the patient's totals are more conserved than compared to the full population


# let's create a list of the OTUs per patient-site that are present in ALL of the samples for that patient-site:
# this creates a dataframe of those OTUs
core.otu.df<-lapply(1:length(patient.list), function(n) df.list2[[n]][ , which(colSums(df.list2[[n]]) == nrow(df.list2[[n]]))])
	# however, this is a dataframe--we just need a list:
core.otu.list<-lapply(1:length(patient.list), function(n) names(df.list2[[n]][ , which(colSums(df.list2[[n]]) == nrow(df.list2[[n]]))]))

percent.core.otus.total<-lapply(1:length(patient.list), function(n) length(df.list2[[n]][ , which(colSums(df.list2[[n]]) == nrow(df.list2[[n]]))]) / length(df.list2[[n]]) *100)
names(percent.core.otus.total)<-patient.list
	# this reflects the % of total OTUs that are observed ALL the time within a patient-site
	# note: there is only 1 number for this list per patient-site (this is not by sample...)
percent.core.otus<-lapply(1:length(patient.list), function(n) length(df.list2[[n]][ , which(colSums(df.list2[[n]]) == nrow(df.list2[[n]]))]) / rowSums(df.list2[[n]]) *100)
	# this is a BY SAMPLE list: it reflects the % of each sample that the core is observed in (still at only presence/absence)
## what is mean across the board:
mean(unlist(percent.core.otus.total))		# 7.78%
min(unlist(percent.core.otus.total))		# 2.59%
max(unlist(percent.core.otus.total))		#17.89%
## what about mean total # of core OTUs across the board:
core.otu.totals<-lapply(1:length(patient.list), function(n) length(df.list2[[n]][ , which(colSums(df.list2[[n]]) == nrow(df.list2[[n]]))]))
names(core.otu.totals)<-patient.list
mean(unlist(core.otu.totals))		# 17
min(unlist(core.otu.totals))		# 6
max(unlist(core.otu.totals))		# 39


######## finally, we can use the list of core OTUs to identify how much ABUNDANCE is taken up by the core (OTUs that we observe all the time)
# what about PERCENT of the community that these OTUs describe at any given moment? should we look at that?

# read in relabund file and format as before:
relabund<-read.table(file="mesal_otufrac.w.meta.txt", header=T, row.names=1)
relabund$site[relabund$site=="Antrum"]<-"Stomach"
relabund<-droplevels(relabund)
relabund$site <- factor(relabund$site, levels = c("Stomach", "Duodenum", "Proximal_Jejunum", "Mid_Jejunum", "Distal_Jejunum", "Stool"))
relabund$sub_site<-paste(relabund$subject, rawotus$location, sep="_")
summary(as.factor(relabund$sub_site))		
	# use multiples.list from above
ra.data<-droplevels(relabund[relabund$sub_site %in% multiples.list, ])
	
# get only OTUs and sampleID, patientID:
ra.df<-ra.data[, c(1, 604, 14:length(ra.data)-2)]		#for now, removing the "rareOTUs" column (will need to add this to final count eventually...
rownames(ra.df)<-ra.data$seqID
ra.df$sub_site<-as.factor(ra.df$sub_site)

# the relabund file has been calculated as relative abundance of ALL otus
# we want to get relative abundance of ONLY core otus...in a list form, again.
# we will use the core.otu.list file we just created, and filter our ra.df multiple times to get it:
ra.df.list<-lapply(1:length(patient.list), function (n) ra.df[rownames(ra.df) %in% sample.list[[n]], colnames(ra.df) %in% core.otu.list[[n]]])
names(ra.df.list)<-patient.list
	# this is a list of the % relative abundance of ONLY the core OTUs within an individual
# you can also get a list of all the OTUs, separated by patient-site
all.ra.df.list<-lapply(1:length(patient.list), function (n) ra.df[rownames(ra.df) %in% sample.list[[n]], ])
names(all.ra.df.list)<-patient.list

# now, we can calculate the totals per sample that make up the relative abundance of each sample:
percent.ra.core.otus<-lapply(1:length(patient.list), function(n) rowSums(ra.df.list[[n]]) )
names(percent.ra.core.otus)<-patient.list
### what is the mean % relative abundance occupied by these:
mean(unlist(percent.ra.core.otus))		# 74.22%
min(unlist(percent.ra.core.otus))		# 10.56%
max(unlist(percent.ra.core.otus))		# 95.596%

#### we have 2 values now:
	- percent.otus.observed:  % of total OTUs ever observed in that person per site within a particular sample
	- percent.core.otus:  % of the OTUs in that sample that are core OTUs (for that patient-site)
	- percent.ra.core.otus:  % of relative abundance that is covered by CORE OTUs (OTUs that are always present in patient-site)
	
# combine these:
alldf<-do.call(rbind, Map(data.frame, p_otus_observed=percent.otus.observed, pn_core_otus=percent.core.otus, pra_core_otus=percent.ra.core.otus))

# combine with metadata, and graph:
f<-merge(relabund[,c(1:13, 604)], alldf, by.x="seqID", by.y="row.names")
#write.table(f, file="mesal.M046_OTU.percent.totals.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

### graph for this individual:
div.boxplot<- function (n, m) {
				plot(data[,n] ~ as.factor(data[,m]), data = data, boxwex=0.7, ylab=names(data[n]), 
					xlab="n", xaxt="n", outline=FALSE, mpg=c(2,1,0), ylim=c(0,100), cex=1, cex.lab=1, cex.axis=1)
				points(data[,n] ~ jitter(as.numeric(as.factor(data[,m]), factor=0)), data = data, col=adjustcolor("grey45", alpha=0.5), pch=21, bg="grey80", cex=0.7)
				print(as.character(unique(data[,m])))
				names<-as.character(levels(data[,m]))
				text(x =  seq(1,length(unique(data[,m])),by=1), y = par("usr")[3]-max(data[n], na.rm=TRUE)/12, srt = 45, adj = 1, labels = paste0(names, "\n n = ", as.numeric(summary(data$site))), xpd = TRUE, cex=0.8)
				stat<-kruskal.test(data[,n]~data[,m], data)
				print(stat$p.value)
				text(length(unique(data[,m]))/2, max(data[,n], na.rm=TRUE), labels=paste0("Kruskal-Wallis test \np < ", signif(stat$p.value, digits=3)), cex=0.5)
				#mtext(names(data[m]), 1, line=2)
				}
par(mar=c(5,3,2,5))
data<-droplevels(f)
par(mfrow=c(3,1))
lapply(c("p_otus_observed"), FUN=div.boxplot, m=c("site"))
title('% n OTUs detected from total OTUs detected at patient site')
lapply(c("pn_core_otus"), FUN=div.boxplot, m=c("site"))
title('% n core OTUs detected at patient site\n(core OTU = OTU always detected at patient site)')
lapply(c("pra_core_otus"), FUN=div.boxplot, m=c("site"))
title('% relative abundance of core OTUs detected at patient site')

Compared to treating these as different admissions:
allp<-read.table(file="mesal_OTU.percent.totals", header=T)
allp$site <- factor(allp$site, levels = c("Stomach", "Duodenum", "Proximal_Jejunum", "Mid_Jejunum", "Distal_Jejunum", "Stool"))

p<-allp[allp$subject %in% c("M046"), ]

par(mar=c(5,3,2,5))
data<-droplevels(p)
par(mfrow=c(3,1))
lapply(c("p_otus_observed"), FUN=div.boxplot, m=c("site"))
title('% n OTUs detected from total OTUs detected at patient site')
lapply(c("pn_core_otus"), FUN=div.boxplot, m=c("site"))
title('% n core OTUs detected at patient site\n(core OTU = OTU always detected at patient site)')
lapply(c("pra_core_otus"), FUN=div.boxplot, m=c("site"))
title('% relative abundance of core OTUs detected at patient site')

### seems to be about the same rate...
### let's graph this for:
	- only the duodenum and mid-jejunum, which were the only "repeat" sites
	- graph it for the individual treated as 3 separate visits and 1 visit at these sites:
	
### add data for m46 as sep entity to the one where they are the same:
allp<-read.table(file="mesal_OTU.percent.totals", header=T)
allp$site <- factor(allp$site, levels = c("Stomach", "Duodenum", "Proximal_Jejunum", "Mid_Jejunum", "Distal_Jejunum", "Stool"))
p<-allp[allp$subject %in% c("M046"), ]

# get only duodenum / mid-jejunum for each dataset:
f<-read.table(file="mesal.M046_OTU.percent.totals.txt", header=T)
f<-droplevels(f[f$site %in% c("Mid_Jejunum", "Duodenum"), ])
p<-droplevels(p[p$site %in% c("Mid_Jejunum", "Duodenum"), ])

# change the 'site' in one of thedata sets to reflect a "new" category (so we can graph on the same panel)
p$site<-as.character(p$site)
p$site[p$site=="Duodenum"]<-"A_Duodenum"
p$site[p$site=="Mid_Jejunum"]<-"A_Mid_Jejunum"
p$site<-as.factor(p$site)

# now, combine and reset levels for desired order:
colnames(f)[14]<-"pID_site"
d<-rbind(p, f)
d$site <- factor(d$site, levels = c("A_Duodenum", "A_Mid_Jejunum", "Duodenum", "Mid_Jejunum"))

now, graph this:
data<-droplevels(d)
lapply(c("p_otus_observed"), FUN=div.boxplot, m=c("site"))
title('% n OTUs detected from total OTUs detected at patient site')
lapply(c("pn_core_otus"), FUN=div.boxplot, m=c("site"))
title('% n core OTUs detected at patient site\n(core OTU = OTU always detected at patient site)')
lapply(c("pra_core_otus"), FUN=div.boxplot, m=c("site"))
title('% relative abundance of core OTUs detected at patient site')


```


##### Extra comparisons found in manuscript:

```
allp<-read.table(file="mesal_OTU.percent.totals", header=T)

# mean %OTUs detected at a given time point across a given patient-site:
mean(allp$p_otus_observed)

# mean % of OTUs that were core:
mean(allp$pn_core_otus)

# mean % RA explained by core OTUs:
mean(allp$pra_core_otus)

# mean percent in only jejunal sites:
jejunum<-allp[allp$site %in% c("Mid_Jejunum", "Proximal_Jejunum", "Distal_Jejunum"), ]
mean(jejunum$pra_core_otus)

# by site:
tapply(allp$p_otus_observed, allp$site, FUN=mean)
tapply(allp$pn_core_otus, allp$site, FUN=mean)
tapply(allp$pra_core_otus, allp$site, FUN=mean)

# same means, but across same subject:
f<-read.table(file="mesal.M046_OTU.percent.totals.txt", header=T)
mean(f$p_otus_observed)
mean(f$pn_core_otus)
mean(f$pra_core_otus)

# by site:
tapply(f$p_otus_observed, f$site, FUN=mean)
tapply(f$pn_core_otus, f$site, FUN=mean)
tapply(f$pra_core_otus, f$site, FUN=mean)


## compare only mid-jejunum, duodenum in this individual since this was the only part sampled multiple times:
# across ALL patients:
tapply(allp$pra_core_otus[allp$site %in% c("Mid_Jejunum", "Duodenum")], allp$site[allp$site %in% c("Mid_Jejunum", "Duodenum")], FUN=mean)
tapply(allp$pn_core_otus[allp$site %in% c("Mid_Jejunum", "Duodenum")], allp$site[allp$site %in% c("Mid_Jejunum", "Duodenum")], FUN=mean)
tapply(allp$p_otus_observed[allp$site %in% c("Mid_Jejunum", "Duodenum")], allp$site[allp$site %in% c("Mid_Jejunum", "Duodenum")], FUN=mean)

# across M46, considering all visits together:
tapply(f$pra_core_otus[f$site %in% c("Mid_Jejunum", "Duodenum")], f$site[f$site %in% c("Mid_Jejunum", "Duodenum")], FUN=mean)
tapply(f$pn_core_otus[f$site %in% c("Mid_Jejunum", "Duodenum")], f$site[f$site %in% c("Mid_Jejunum", "Duodenum")], FUN=mean)
tapply(f$p_otus_observed[f$site %in% c("Mid_Jejunum", "Duodenum")], f$site[f$site %in% c("Mid_Jejunum", "Duodenum")], FUN=mean)

# finally, in M46, but considering each visit separately:
m<-allp[allp$subject %in% c("M046"), ]
tapply(m$pra_core_otus[m$site %in% c("Mid_Jejunum", "Duodenum")], m$site[m$site %in% c("Mid_Jejunum", "Duodenum")], FUN=mean)
tapply(m$pn_core_otus[m$site %in% c("Mid_Jejunum", "Duodenum")], m$site[m$site %in% c("Mid_Jejunum", "Duodenum")], FUN=mean)
tapply(m$p_otus_observed[m$site %in% c("Mid_Jejunum", "Duodenum")], m$site[m$site %in% c("Mid_Jejunum", "Duodenum")], FUN=mean)


```

