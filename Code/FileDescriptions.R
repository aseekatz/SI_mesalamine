#### File Descriptions
##### Anna M. Seekatz
##### 08.10.18

###### Raw data file used in analyses:
	- sample information and de-identified metadata: mesalamine_metadata.txt
	- OTU counts (filtered in mothur): mesalamine.final.shared
	- OTU classifications: ibumes.final.0.03.cons.taxonomy
		- note: this includes all OTUs processed together with a larger study 
	- summary files calculated in mothur (based on thetayc):
		- mesalamine.final.thetayc.0.03.lt.nmds.axes
		- mesalamine.final.thetayc.0.03.lt.pcoa.axes
		- mesalamine.final.groups.summary
	- phylotype file (up to genus-level): ibumes.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.tax.summary
		- note: this included counts for a larger data set, and has been filtered
	- pairwise calculations of different beta diversity measures: mesalamine.final.summary
		- this was combined with metadata 
		
###### Data files merged with meta data or manipulated otherwise:
	- mesal_subs.otucounts.w.meta.txt, mesal_otufrac.w.meta.txt --> relative abundance of OTUs
	- mesal_summary.txt --> metadata with summary scores
	- mesalamine.taxnames03.txt --> specific OTU classifications
	- mesal_genfrac2p.meta.txt, mesal_genfrac2p.txt --> phylotype files 
	- mesal_alldist.txt --> beta diversity distance calculations
	
	
##### Generating files used for analysis of data (from raw mothur data files)

```
#-------------------# # diversity and PCOA summary calculations 
### merging alpha diversity, PCOA/NMDS calculations:
meta<-read.table(file="mesalamine_metadata.txt", header=TRUE)
	# added nseqs from previous mothurfile, which had all the samples
allsums<-read.table(file="../ibumes//mothurfiles/ibumes.final.groups.summary", header=TRUE)
nseqs<-allsums[, c("group", "nseqs")]
meta<-merge(meta, nseqs, by.x="sampleID", by.y="group")

pcoa3<-read.table(file="mothurfiles/mesalamine.final.thetayc.0.03.lt.pcoa.axes", header=TRUE)
	pcoa3<-pcoa3[,1:4]
	colnames(pcoa3)[2:4] <- paste("pcoa03", colnames(pcoa3)[2:4], sep = "_")
	colnames(pcoa3)[1]<-"sampleID"
nmds3<-read.table(file="mothurfiles/mesalamine.final.thetayc.0.03.lt.nmds.axes", header=TRUE)
	nmds3<-nmds3[1:4]
	colnames(nmds3)[2:4] <- paste("nmds03", colnames(nmds3)[2:4], sep = "_")
	colnames(nmds3)[1]<-"sampleID"
sum3<-read.table(file="mothurfiles/mesalamine.final.groups.summary", header=TRUE)
	sum3<-subset(sum3, select=-c(label, nseqs))
	colnames(sum3)[2:16] <- paste(colnames(sum3)[2:16], "03", sep = "_")
	colnames(sum3)[1]<-"sampleID"

combined.pcoa<-merge(meta, pcoa3, by=c("sampleID"), all.x=TRUE)
combined.nmds<-merge(combined.pcoa, nmds3, by=c("sampleID"), all.x=TRUE)
combined.sum<-merge(combined.nmds, sum3, by.x=c("sampleID"), all.x=TRUE)

# note: some samples were missing in pcoa/nmds/sums file
# checked on these in the larger (ibumes) groups.summary file--they had failed
# these were the failed ones, or the contaminated one:
combined.sum[is.na(combined.sum$sobs_03), c("sampleID", "sobs_03", "nseqs")]
#write.table(combined.sum, 'mesal_summary.txt', quote=FALSE,sep="\t", col.names=TRUE, row.names=FALSE)

# also added a 'failed' status to the meta, and re-wrote it:
meta$seq_status[meta$nseqs<3000]<-"failed"
meta$seq_status[meta$nseqs>3000]<-"OK"
meta$seq_status[meta$sampleID=="M046_IIC_1h_1A_and_M046_IIC_1h_2A"]<-"contam"
meta<-subset(meta, select=-c(nseqs))

# let's add some more info about the stomach location to the meta file:
meta$site<-as.character(meta$location)
meta$location<-as.character(meta$location)
meta$location[meta$site %in% c("Distal_Jejunum", "Mid_Jejunum", "Proximal_Jejunum")]<-"Jejunum"
meta$site[meta$location %in% c("Stomach") & meta$pH < 3]<-"Antrum"
meta$location<-as.factor(meta$location)
meta$site<-as.factor(meta$site)
	# note: the stomach designation in this could be either fundus or body
# let's also add another subject ID to reflect the other subject
meta$pID<-as.character(meta$subject)
meta$pID[meta$subject %in% c("M046") & meta$Med %in% c("A")]<-"M046-A"
meta$pID[meta$subject %in% c("M046") & meta$Med %in% c("B")]<-"M046-B"
meta$pID[meta$subject %in% c("M046") & meta$Med %in% c("C")]<-"M046-C"
meta$pID<-as.factor(meta$pID)
#write.table(meta, 'mesalamine_metadata.txt', quote=FALSE,sep="\t", col.names=TRUE, row.names=FALSE)
	# remember to rewrite the meta to the summary above


#-------------------# # relative abundance of OTUs (w/ subsampling and filtering)
### Getting relative abundance of OTUs (using ALL OTUs in mesalamine data set):
	# note: if you want to skip subsampling in R, you could also use the filtered file created by mothur (this has only 151 OTUs)
	#shared<-read.table(file="mothurfiles/mesalamine.final.0.03.filter.shared", header=TRUE, row.names=2)
library(vegan)
rawshared<-read.table(file="mothurfiles/mesalamine.final.shared", header=TRUE, row.names=2)
rawshared$numOtus <- NULL
rawshared$label <- NULL

# if you want to subsample, follow through with the following:
	# if wanting to subsample by quantile:
#sub_size <- floor(as.numeric(quantile(rowSums(shared), probs=0.09)))
	# if want to subsample by minimum:
sub_size <- as.numeric(min(rowSums(rawshared)))


# Rarefy shared file to appropriate min (subsample):
shared <- as.data.frame(t(rawshared))
shared <- shared[, colSums(shared) >= sub_size]	#if you are filtering by quantile, you will lose more samples
set.seed(1984)	# to get consistent results
for (index in 1:ncol(shared)){
  shared[,index] <- t(rrarefy(shared[,index], sample=sub_size))}
rm(index, sub_size)

# remove OTUs that are no longer existent due to rarification (0 totals):
#colSums(shared)		# if you want to check that everything was rarefied
shared<-shared[rowSums(shared) > 0,]
dim(shared)
	# this reduced OTUs to 1234

# last step, if you want, is to filter by prevalence
# Filter out columns that have values in at least 2 samples (ignores first column if needed)
filter_table <- function(data) {
  drop <- c()
  if (class(data[,1]) != 'character') {
    if (sum(data[,1] != 0) < 2) {
      drop <- c(drop, colnames(data)[1])
    }
  }
  for (index in 2:ncol(data)) {
    if (sum(data[,index] != 0) < 2) {
      drop <- c(drop, colnames(data)[index])
    }
  }
  filtered_data <- data[,!(colnames(data) %in% drop)]
  #filtered_data$rareOTUs <- rowSums(data[,(colnames(data) %in% drop)])
  #return(filtered_data)
}
filtered.shared <- as.data.frame(t(shared))
filtered.shared <- filter_table(filtered.shared) # OTU set reduced to 561--this is probably more appropriate than the method used in mothur, since all the samples are so different from each other
dim(filtered.shared)
rowSums(filtered.shared)
	# there appears to be some weirdness in some samples, where ~50-75% of the OTUs are 1 OTU that never appears in another: 
	#M048_IIB_GI_2h_2A, M063_IIC_1h_1A, M061_IIA_4h_2A, M046_IIA_4h_2A
	#track back to see what is going on by repeating above steps WITHOUT filtering steps, then:
filtered.otus<-colnames(filtered.shared)
shared2<-as.data.frame(t(shared))
test<-shared2[,!(colnames(shared2) %in% filtered.otus)]
test$remainder<-rowSums(shared2[,(colnames(shared2) %in% filtered.otus)])
relabund<-as.matrix(test/rowSums(test))
relabund<-as.data.frame(relabund)
relabund<-relabund[order(relabund$remainder),]
relabund<-subset(relabund, select=-c(remainder))
max.relabund<-relabund[,colSums(relabund) > 0.01]
weirdOTUs<-as.character(colnames(max.relabund))
	# in looking at these, we can also see that some OTUs (these weird ones) are also being removed despite showing up in more than 2 samples
	# not sure if this code is correct...
	# let's look at what kind of OTUs they are?
# note: fixed this by adding 'rareOTUs' column in filter_table script...
# however, if we want to include these 'weird otus', let's add those, and delete the rest
	# get OTU numbers from filtered list, as well as these weird OTUs
#setdiff(weirdOTUs, colnames(filtered.shared))
included.otus<-unique(c(weirdOTUs, colnames(filtered.shared)))
shared2<-as.data.frame(t(shared))
filtered.otutable<-shared2[, colnames(shared2) %in% included.otus]
filtered.otutable$rareOTUs <- rowSums(shared2[,!(colnames(shared2) %in% included.otus)])
	# I can live with this amount...
	# way to phrase this filtering is 'present in at least 2 samples across the data set OR present at a relative abundance of > 0.01% in a single sample

# convert to relative abundance:
otufrac<-as.data.frame(as.matrix(filtered.otutable)/rowSums(as.matrix(filtered.otutable))*100)
otufrac$seqID<-rownames(otufrac)
otufrac.meta<-merge(meta, otufrac, by="seqID")
otucounts.meta<-merge(meta, filtered.otutable, by.x="seqID", by.y='row.names')
#write.table(otufrac.meta, file="mesal_otufrac.w.meta.txt", sep="\t", quote=FALSE, col.names=NA)
#write.table(otucounts.meta, file="mesal_subs.otucounts.w.meta.txt", sep="\t", quote=FALSE, col.names=NA)

#-------------------# # OTU assignment names
### renaming taxonomy file:
taxonomy_file<-read.table(file="../ibumes/mothurfiles/ibumes.final.0.03.cons.taxonomy", header=TRUE)
tax <- taxonomy_file$Taxonomy
tax <- gsub("\\(\\d*\\)", "", tax)
tax <- gsub(";unclassified", "", tax)
tax <- gsub("_1", "", tax)
tax <- gsub(";$", "", tax)
tax <- gsub("/.*", "", tax)
tax <- gsub(".*;", "", tax)
tax.names <-paste(taxonomy_file$OTU, tax, sep="_")
#tax.names <-gsub("000", "", tax.names)
taxonomy_file$taxname<-tax.names
phylum <- taxonomy_file$Taxonomy
phylum <- gsub("\\(\\d*\\)", "", phylum)
phylum <- gsub("Bacteria;", "", phylum)
phylum <- gsub(";$", "", phylum)
phylum <- gsub(";.*", "", phylum)
taxonomy_file$phylum<-phylum
tax<-taxonomy_file
tax<-tax[order(tax$OTU),]

# then, filter to get OTUs identified in specific CX study set (after filtering)
shared<-read.table(file="mothurfiles/mesalamine.final.shared", header=TRUE)
#tax03<-read.table(file="../rushfinal.0.03.taxonomy.names.txt", header=TRUE)
tax03<-tax
OTUs<-as.character(colnames(shared)[4:ncol(shared)])
filtered.tax03<-tax03[tax03$OTU %in% OTUs, ]
#write.table(filtered.tax03, file="mesalamine.taxnames03.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


#-------------------# # phylotype (genus-level) relative abundance: 
### collating genus-level phylotype file:
taxo<-read.table(file="mothurfiles/ibumes.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.tax.summary", header=TRUE)
taxo.samples<-taxo[, 6:ncol(taxo)]
cx<-read.table(file="mesal_otufrac.w.meta.txt", header=TRUE, sep="\t")
names<-as.character(cx$seqID)
taxo.names<-taxo.samples[, colnames(taxo.samples) %in% names]	#get only mesalamine samples
taxo.filtered<-cbind(taxo[1:5], taxo.names)
taxo.present<-taxo.filtered[rowSums(taxo.filtered[6:ncol(taxo.filtered)]) > 0, ]
#write.table(taxo.present, file="mothurfiles/mesal.only_wang.tax.summary.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)	

# filter out genus-level assignments only, and assign a phylum:
tax<-taxo.present
tax2<-tax[which(tax$taxlevel==2), ]
tax2[, c("rankID", "taxon")]
tax6<-tax[which(tax$taxlevel==6), ]
tax6$rankID<-gsub("^0.1.1.*", "20_Euryarchaeota", tax6$rankID)
tax6$rankID<-gsub("^0.2.1\\..*", "10_Acidobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.2\\..*", "04_Actinobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.4\\..*", "11_Bacteria_unclassified", tax6$rankID)
tax6$rankID<-gsub("^0.2.5\\..*", "01_Bacteroidetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.6\\..*", "20_Candidatus_Saccharibacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.8\\..*", "20_Chloroflexi", tax6$rankID)
tax6$rankID<-gsub("^0.2.10..*", "20_Cyanobacteria/Chloroplast", tax6$rankID)
tax6$rankID<-gsub("^0.2.11..*", "20_Deinococcus-Thermus", tax6$rankID)
tax6$rankID<-gsub("^0.2.12..*", "02_Firmicutes", tax6$rankID)
tax6$rankID<-gsub("^0.2.13..*", "06_Fusobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.20..*", "20_Parcubacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.21..*", "20_Planctomycetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.22..*", "03_Proteobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.23..*", "20_SR1", tax6$rankID)
tax6$rankID<-gsub("^0.2.24..*", "09_Spirochaetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.25..*", "08_Synergistetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.26..*", "07_Tenericutes", tax6$rankID)
tax6$rankID<-gsub("^0.2.28..*", "05_Verrucomicrobia", tax6$rankID)
tax6$rankID<-gsub("^0.3.1..*", "11_unknown_unclassified", tax6$rankID)
colnames(tax6)[2]<-"phylum"

# filter some columns, turn into a matrix and check for duplication:
subtax6<-subset(tax6, select=-c(taxlevel, daughterlevels))
subtax6<-subtax6[order(subtax6$phylum, -subtax6$total), ]
taxmatrix<-subtax6[, c(4:ncol(subtax6))]
which(duplicated(subtax6$taxon))
	# none are duplicated, but follow trhough with remainder if they were:
	# fix the duplicated row:
	#subtax6$taxon<-as.character(subtax6$taxon)
	#subtax6$taxon[268]<-"Actinobacteria_unclassified2"
	#subtax6$taxon<-as.factor(subtax6$taxon)
rownames(taxmatrix)<-subtax6$taxon
#genera<- taxmatrix[, colSums(taxmatrix)>5000,]		# should already be filtered
	# get rel. abund fraction:
genera<-taxmatrix
genmatrix<-as.data.frame(t(genera))
genera.fr<-genmatrix/rowSums(genmatrix)*100
genus.fr<-t(genera.fr)
all.genera<-cbind(subtax6[1:3], genus.fr)
#write.table(all.genera, file="mesal_all.genera.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
	# get top 1%:
phyla<-subtax6[1:3]
genus1<- genus.fr[rowSums(genus.fr>=1)>=1,]
namelist<-as.character(rownames(genus1))
phyla1p<-phyla[phyla$taxon %in% namelist, ]
genera1<-cbind(phyla1p, genus1)
	# get top 2%
genus2<- genus.fr[rowSums(genus.fr>=2)>=2,]
namelist<-as.character(rownames(genus2))
phyla2p<-phyla[phyla$taxon %in% namelist, ]
genera2<-cbind(phyla2p, genus2)
#write.table(genera2, file="mesal_genfrac2p.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


# let's also add metadata and colors for later:
summary(as.factor(genera2$phylum))
	bac<-colorRampPalette(c("darkgreen","green3","lightgreen","seagreen", "lightgreen"))(n=10)
	firm<-colorRampPalette(c("midnightblue","mediumblue","blue3","blue","dodgerblue4","dodgerblue1","deepskyblue4","deepskyblue1","skyblue3","skyblue","steelblue4","steelblue1","royalblue4","royalblue1","slateblue4","purple3","orchid3","plum4","plum1","pink3","pink","lightpink1","lightpink3","palevioletred4","palevioletred1","magenta4","deeppink4","mediumvioletred","magenta3","magenta1","thistle"))(n=25)
	pro<-colorRampPalette(c("yellow2","darkgoldenrod3","goldenrod2","orange2","yellow4"))(n=12)
	actino<-colorRampPalette(c("tan", "brown", "darkred"))(n=5)
	verruco<-c("hotpink")
	fuso<-colorRampPalette(c("tan", "red", "darkred"))(n=3)
	uncl<-c("black")
	other<-colorRampPalette(c("cyan", "darkcyan"))(n=2)
	#other<-colorRampPalette(c("grey50"))(n=1)
color<-c(bac, firm, pro, actino, verruco, fuso, other, uncl)
genera2<-cbind(genera2[1:3], color, genera2[4:ncol(genera2)])
#write.table(genera2, file="mesal_genfrac2p.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# read in file and combine with meta:
#genera2<-read.table(file="mesal_genfrac2p.txt", header=TRUE, sep="\t", comment.char = "")
genbar<-genera2
rownames(genbar)<-genbar$taxon
	rm_g<-subset(genbar, select =-c(phylum, taxon, color, total) )
	barg<-as.data.frame(t(rm_g))
	barg$other<-100-rowSums(barg)
	others<-100-colSums(barg)
	barg$sampleID<-rownames(barg)
	col.gen<-c(as.character(genbar$color), "grey47")
	#barg$sampleID<-gsub("X19_", "19_", barg$sampleID)
bar<-merge(meta, barg, by.x=c("seqID"), by.y=c("sampleID"), all.x=TRUE)
#write.table(bar, 'mesal_genfrac2p.meta.txt',quote=FALSE,sep="\t", col.names=TRUE, row.names=FALSE)


#-------------------# # pairwise distances of samples, with metadata: 
# let's combine our meta data with the beta diversity pairwise distances calculated in mothur:
allsums<-read.table(file="mothurfiles/mesalamine.final.summary", header=TRUE)
meta<-read.table(file="mesalamine_metadata.txt", header=TRUE)
meta<-meta[meta$seq_status %in% c("OK"), ]
	# if combining antrum etc, after all:
meta$site[meta$site %in% c("Antrum")]<-"Stomach"			#note: not distinguishing between antrum/stomach in final analysis
meta<-droplevels(meta)
meta$site<-factor(meta$site, levels = c("Stomach", "Duodenum", "Proximal_Jejunum", "Mid_Jejunum", "Distal_Jejunum", "Stool"))

var<-meta[, 1:12]

# combine info per sample-sample comparison:
#add sample information for 'sample1':
m<-merge(var, allsums, by.x=c("seqID"), by.y=c("label"))	
s1names<-paste("s1", colnames(m[1:12]), sep = "_")
allnames<-c(s1names, colnames(m[13:21]))			
colnames(m) <- allnames

#add sample information for 2nd sample: ('sample2'):
m2<-merge(var, m, by.x=c("seqID"), by.y=c("comparison")) 			
s2names<-paste("s2", colnames(m2[1:12]), sep = "_")
allnames<-c(s2names, colnames(m2[13:32]))			
colnames(m2) <- allnames
#write.table(m2, file="mesal_alldist.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# this file represents all defined distance measures between ALL samples

```

##### Also made a table S1 that includes OTUfrac (with names), some summary info, metadata
	- note: this file includes negative controls as well
	- indicated in a column whether or not this was included in the study analysis


```
#-------------------# # Combine metadata, OTUfrac, and summary file for Table S1
### note: including negative controls here, too
### these were any samples included on Run51, Run197 from the full sampleinfo sheet

sums<-read.table(file="mesal_summary.txt", header=T)		#includes samples we did NOT include, as well
meta<-read.table(file="../ibumes/ibumes_temp_data.txt", header=T)		# full run information
shared<-read.table(file="../ibumes/mothurfiles/ibumes.final.shared", header=T, row.names=2)	#all OTU counts for all samples

# filter out samples we want (166 total, plus controls)
# subset samples that were used in mesalamine project
mes.samples<-meta[meta$Run %in% c("Run51", "Run197"), ]		#184 samples total
samples<-as.character(mes.samples$seqID)
shared<-shared[rownames(shared) %in% samples, ]

# let's work with raw data first
### now, process shared in same way that it was processed above to produce the OTUfrac files
shared$numOtus <- NULL
shared$label <- NULL
	# use the otufrac file we previously made to get only OTUs used in this file for the extra (not used) samples
mesal_otus<-read.table(file="mesal_otufrac.w.meta.txt", header=T)
mesal_names<-as.character(mesal_otus$seqID)
mesal_otu_names<-as.character(colnames(mesal_otus[14:length(mesal_otus)-1]))

xtra_names<-samples[!samples %in% mesal_names]		#these are the samples that were not previously excluded
# get otu totals for only these samples
xtra_otus<-shared[rownames(shared) %in% xtra_names, colnames(shared) %in% mesal_otu_names]
# but, also need to total for rareOTUs!
	# so, get totals of OTUs NOT in the official list
rare<-shared[rownames(shared) %in% xtra_names, !colnames(shared) %in% mesal_otu_names]
# now, just total this list to get a total for them
xtra_otus$rareOTUs<-rowSums(rare)
# finally, convert this to a percentage
xtra_otufrac<-as.data.frame(as.matrix(xtra_otus)/rowSums(as.matrix(xtra_otus))*100)

### now, add all together:
rownames(mesal_otus)<-mesal_otus$seqID
otus<-rbind(xtra_otufrac, mesal_otus[14:length(mesal_otus)])
	# also, lets add the taxnames to this
tax<-read.table(file="mesalamine.taxnames03.txt", header=T)
taxnames<-as.character(tax[tax$OTU %in% colnames(otus), c("taxname")])
names(otus)[1:length(otus)-1]<-taxnames

# merge with summary info
ts1<-merge(sums, otus, by.x="seqID", by.y="row.names", all.y=T)
#write.table(ts1, file="TableS1_sampleinfo.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# will also need to add nseqs to controls:
rowSums(shared)		# just added these to the table itself



```