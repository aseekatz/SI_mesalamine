#### Mesalamine Study Manuscript:  Figure 2, similarity and stability (intra- and inter-individuality)
##### Anna M. Seekatz
##### 08.16.18

###### Figures created:
	- Fig. 2A, B:  inter- and intra-individual distances--comparing all sites
	- Fig. 2C, D:  inter- and intra-individual distances--comparing stomach, duodenum, jejunum to themselves


####### Files used:
	- mesalamine.final.summary
	- mesalamine_metadata.txt
	- mesal_alldist.txt
	
```
library(reshape2)
library(plyr)
library(gplots)
library(RColorBrewer)

# get file with distances:
rawdata<-read.table(file="mesal_alldist.txt", header=T)

# note: this file already filtered out unsuccessful samples
# it also has antrum/stomach --> stomach
# let's try to get the squares in order by giving the body parts numbers:
rawdata$s1_site<-as.character(rawdata$s1_site)
rawdata$s2_site<-as.character(rawdata$s2_site)
rawdata$s1_site[rawdata$s1_site =="Duodenum"]<-"2_Duodenum"
rawdata$s1_site[rawdata$s1_site =="Stomach"]<-"1_Stomach"
rawdata$s1_site[rawdata$s1_site =="Proximal_Jejunum"]<-"3_Proximal_Jejunum"
rawdata$s1_site[rawdata$s1_site =="Mid_Jejunum"]<-"4_Mid_Jejunum"
rawdata$s1_site[rawdata$s1_site =="Distal_Jejunum"]<-"5_Distal_Jejunum"
rawdata$s1_site[rawdata$s1_site =="Stool"]<-"6_Stool"
rawdata$s2_site[rawdata$s2_site =="Duodenum"]<-"2_Duodenum"
rawdata$s2_site[rawdata$s2_site =="Stomach"]<-"1_Stomach"
rawdata$s2_site[rawdata$s2_site =="Proximal_Jejunum"]<-"3_Proximal_Jejunum"
rawdata$s2_site[rawdata$s2_site =="Mid_Jejunum"]<-"4_Mid_Jejunum"
rawdata$s2_site[rawdata$s2_site =="Distal_Jejunum"]<-"5_Distal_Jejunum"
rawdata$s2_site[rawdata$s2_site =="Stool"]<-"6_Stool"
rawdata$s1_site<-as.factor(rawdata$s1_site)
rawdata$s2_site<-as.factor(rawdata$s2_site)

# these are pairwise distances
# sometimes, sample types will not all be in same column
# thus, need to redefine a variable comparison that collapses these before doing any kind of mean comparison for a heatmap
d<-rawdata
all.comps<-paste0(d$s1_site, " ", d$s2_site)
split.comps<-strsplit(all.comps, split=" ")
ordered.comps<-sapply(split.comps, function(x) sort(x), simplify=F)
d$comp1<-sapply(ordered.comps, function(x) paste0(x[1]))
d$comp2<-sapply(ordered.comps, function(x) paste0(x[2]))
	# now, you have these distances, ready to be collapsed into averages per bodysite comparison
	
# let's also define which ones are inter, and which ones are intra-comparisons:
d$type[d$s1_pID==d$s2_pID]<-"intra"
d$type[!d$s1_pID==d$s2_pID]<-"inter"
#write.table(d.intra, 'mesal_test.txt', quote=FALSE,sep="\t", col.names=TRUE, row.names=FALSE)

# let's average out distance of choice:
#sums$comp1<-factor(sums$comp1, levels = c("Stomach", "Duodenum", "Proximal_Jejunum", "Mid_Jejunum", "Distal_Jejunum", "Stool"))
#sums$comp2<-factor(sums$comp2, levels = c("Stomach", "Duodenum", "Proximal_Jejunum", "Mid_Jejunum", "Distal_Jejunum", "Stool"))

# calculate differences for intra:
d.intra<-d[d$type %in% c("intra"), ]
sums <- ddply(d.intra, c("comp2", "comp1"), summarise,
               N    = length(thetayc),
               mean = mean(thetayc),
               sd   = sd(thetayc),
               se   = sd / sqrt(N) )
#sums$comp1<-factor(sums$comp1, levels = c("Stomach", "Duodenum", "Proximal_Jejunum", "Mid_Jejunum", "Distal_Jejunum", "Stool"))
#sums$comp2<-factor(sums$comp2, levels = c("Stomach", "Duodenum", "Proximal_Jejunum", "Mid_Jejunum", "Distal_Jejunum", "Stool"))
sums<-sums[order(sums$comp1), ]
sums.matrix<-acast(sums, comp1~comp2, value.var="mean")
sums.theta.intra<-sums.matrix

### calculate differences for inter:
d.inter<-d[d$type %in% c("inter"), ]
sums <- ddply(d.inter, c("comp2", "comp1"), summarise,
               N    = length(thetayc),
               mean = mean(thetayc),
               sd   = sd(thetayc),
               se   = sd / sqrt(N) )
#sums$comp1<-factor(sums$comp1, levels = c("Stomach", "Duodenum", "Proximal_Jejunum", "Mid_Jejunum", "Distal_Jejunum", "Stool"))
#sums$comp2<-factor(sums$comp2, levels = c("Stomach", "Duodenum", "Proximal_Jejunum", "Mid_Jejunum", "Distal_Jejunum", "Stool"))
sums<-sums[order(sums$comp1), ]
sums.matrix<-acast(sums, comp1~comp2, value.var="mean")
sums.theta.inter<-sums.matrix

# make a heatmap out of this:
# graph--subset samples:
# breaks need to be +1 more than colors
my.breaks = c(seq(0,0.30,length=1),  								#just remember: your breaks must always = 1 more than the number of colors you are defining!
               seq(0.30,0.50,length=10),
               seq(0.50,0.70,length=10),
               seq(0.70,0.90,length=10),
               seq(0.90,1,length=5))
# does R 3.3.0 make these differently?!
edited.breaks<-my.breaks[!duplicated(my.breaks)]
# colors
value.cols <- colorRampPalette(c("lightgoldenrod1", "gold","darkorange2", "firebrick1", "firebrick", "firebrick4",  "darkred"))(n = length(edited.breaks)-2)			#you can put in as many colors as you want here
my.col<-c("bisque", value.cols)

## breaks/colors, second choice:
my.breaks = c(seq(0,0.30,length=1),  								#just remember: your breaks must always = 1 more than the number of colors you are defining!
               seq(0.30,0.40,length=5),
               seq(0.40,0.50,length=5),
               seq(0.50,0.60,length=5),
               seq(0.60,0.70,length=5), 
               seq(0.70,0.80,length=5),
               seq(0.80,0.90,length=5),
               seq(0.90,1.00,length=5))
value.cols<-colorRampPalette(brewer.pal(11,"YlOrRd"))(n = length(edited.breaks)-2)
my.col<-c("lightyellow", value.cols)

## for legend:
my.leg.breaks = c(seq(0,0.30,length=1),  								#just remember: your breaks must always = 1 more than the number of colors you are defining!
               seq(0.30,0.40,length=1),
               seq(0.40,0.50,length=1),
               seq(0.50,0.60,length=1),
               seq(0.60,0.70,length=1), 
               seq(0.70,0.80,length=1),
               seq(0.80,0.90,length=1),
               seq(0.90,1.00,length=1))
value.leg.cols<-colorRampPalette(brewer.pal(11,"YlOrRd"))(n = length(my.leg.breaks)-2)
my.leg.col<-c("lightyellow", value.leg.cols)

# you can make a heatmap of each of these:
# Figure 2A:
heatmap.2(sums.theta.inter,
                    notecol="black",     
                    density.info="none",
                    key.xlab="",
                    key.title="",  
                    trace="none",         
                    margins =c(5,10),    
                    col=my.col,        
                    breaks=edited.breaks,    
                    key=T,  
                  	dendrogram="none",    
                    Colv=F,
                    Rowv=F,
                    srtCol=45,
                    cexRow= 1,
                    cexCol = 1,
                    keysize=1,
                    lwid = c(2,5),
                    lhei = c(2,5),
                    labCol=gsub("_", " ", colnames(sums.theta.inter)),
                    labRow=gsub("_", " ", rownames(sums.theta.inter)),
					#RowSideColors = as.character(rushsum2$color),
					#ColSideColors = as.character(phyla.col)
					na.color="grey80"
					)
legend("top",legend=c("0", "0.3", "0.5", "0.6", "0.7", "0.8", "1.0"), col=my.leg.col, cex=0.8, pch=19)

# Figure 2B:
heatmap.2(sums.theta.intra,
                    notecol="black",     
                    density.info="none",
                    key.xlab="",
                    key.title="",  
                    trace="none",         
                    margins =c(5,10),    
                    col=my.col,        
                    breaks=edited.breaks,    
                    key=T,  
                  	dendrogram="none",    
                    Colv=F,
                    Rowv=F,
                    srtCol=45,
                    cexRow= 1,
                    cexCol = 1,
                    keysize=1,
                    lwid = c(2,5),
                    lhei = c(2,5),
                    labCol=gsub("_", " ", colnames(sums.theta.intra)),
                    labRow=gsub("_", " ", rownames(sums.theta.intra)),
					#RowSideColors = as.character(rushsum2$color),
					#ColSideColors = as.character(phyla.col)
					na.color="grey80"
					)
#legend("top",legend=c("", "", ""), col=c("", "", ""), cex=0.8, pch=19)

### notes:
# if you look at the grids, you see that there are blank spots
# this is because not everyone provided samples for every region:
	# only one individual gave Distal_jejunum samples (so there is no 'inter-comparison' here)
	# For the individual who gave a distal jejunal sample, they did not provide stomach, proximal jejunum, or stool at that point (so no 'intra-comparisons' here)
	# similarly, there are no stool samples available for intra-comparison against proximal jejunum or distal jejunum samples
	

### What about similarity at the hour? 
# exclude stool from this--let's just look at ALL the samplings
d.nostool<-d[(!d$s1_site %in% c("Stool") | d$s2_site %in% c("Stool")), ]

d.intra.nostool<-d.nostool[d.nostool$type %in% c("intra"), ]

# let's also do by matching hour...
d.match.samples<-d.intra.nostool[d.intra.nostool$s1_Hour==d.intra.nostool$s2_Hour, ]
sums <- ddply(d.match.samples, c("comp2", "comp1"), summarise,
               N    = length(thetayc),
               mean = mean(thetayc),
               sd   = sd(thetayc),
               se   = sd / sqrt(N) )
sums<-sums[order(sums$comp1), ]
sums.matrix<-acast(sums, comp1~comp2, value.var="mean")
sums.samplematch<-sums.matrix

heatmap.2(sums.samplematch,
                    notecol="black",     
                    density.info="none",
                    key.xlab="",
                    key.title="",  
                    trace="none",         
                    margins =c(5,10),    
                    col=my.col,        
                    breaks=edited.breaks,    
                    key=T,  
                  	dendrogram="none",    
                    Colv=F,
                    Rowv=F,
                    srtCol=45,
                    cexRow= 1,
                    cexCol = 1,
                    keysize=1,
                    lwid = c(2,5),
                    lhei = c(2,5),
                    labCol=gsub("_", " ", colnames(sums.samplematch)),
                    labRow=gsub("_", " ", rownames(sums.samplematch)),
					#RowSideColors = as.character(rushsum2$color),
					#ColSideColors = as.character(phyla.col)
					na.color="grey80"
					)



```

####### Fig. 2C, D, intra vs. inter-individual distances
	- mesal_alldist.txt

```
# ### let's look at intra versus inter-individual differences, by body site:
data<-read.table(file="mesal_alldist.txt", header=TRUE)
data$s1_site<-factor(data$s1_site, levels = c("Stomach", "Duodenum", "Proximal_Jejunum", "Mid_Jejunum", "Distal_Jejunum", "Stool"))
data$s2_site<-factor(data$s2_site, levels = c("Stomach", "Duodenum", "Proximal_Jejunum", "Mid_Jejunum", "Distal_Jejunum", "Stool"))

# define intra or inter-individual:
# note: for this first pass, using subject M046, who did this 3x, as separate observances
data$comp[data$s1_pID==data$s2_pID]<-"intra"
data$comp[!data$s1_pID==data$s2_pID]<-"inter"

# define organ-to-organ:
# first, let's do nonspecific site and combine with intra/inter:
df<-data[data$s1_location==data$s2_location,]
df$COMP<-paste(df$comp, df$s1_location, sep="_")
	# not comparing stool:
df<-droplevels(df[!df$s1_location %in% c("Stool"),])

### graph:
# if you want to combine all these into the same graph:
#plotspots<-rbind(c(1, 1, 1, 2, 2, 2, 3, 3), c(1, 1, 1, 2, 2, 2, 3, 3), c(1, 1, 1, 2, 2, 2, 3, 3))
#layout(plotspots)

# A: by general site:
df$COMP<-as.factor(df$COMP)
df$COMP <- factor(df$COMP, levels = c("intra_Stomach", "intra_Duodenum", "intra_Jejunum", "inter_Stomach", "inter_Duodenum","inter_Jejunum"))

comp.col <- function(n) {
colorvec <- vector(mode="character", length=length(n))
for (i in 1:length(n)) {
colorvec[i] = "light grey"
if ( n[i] == "intra" ) {
colorvec[i] = "dodgerblue"
}
if ( n[i] == "inter" ) {
colorvec[i] = "red"
}
}
c(colorvec)
}
d<-df

### Figure 2C, D:
#plot(d$COMP, d$thetayc, type = "n", las=2, cex.axis=0.8, ylab=expression(paste(italic(theta [YC]))), xaxt="n")
#plot(d$COMP, d$braycurtis, type = "n", las=2, cex.axis=0.8, ylab="Dissimilarity (Bray-curtis distance)", xaxt="n")
# label example: xlab= expression(bla~bla~italic(bli~bli)~bla~bold(blom)~italic(bla)))
plot(d$COMP, d$thetayc, type = "n", las=2, cex.axis=0.8, ylab=expression("Dissimilarity,"~italic(theta)["YC"] ), xaxt="n", boxwex=0.5)
points(jitter(as.numeric(d$COMP), amount=0.2), d$thetayc, col = adjustcolor(comp.col(d$comp), alpha=0.5), pch = 16, cex=0.5)
names<-gsub(".*_", "", as.character(levels(d$COMP)))
text(x =  seq(1,length(unique(d$COMP)),by=1), y = par("usr")[3]-max(d$thetayc, na.rm=TRUE)/11, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=0.7)
text(x =  c(2,5), y = c(-.5, -.5), adj = 1, labels = c("intra", "inter"), xpd = TRUE, cex=0.7)

## stats:
# split by location:
d1<-d[d$COMP %in% c("intra_Stomach", "inter_Stomach"), ]
t.test(thetayc~COMP, data=d1)
wilcox.test(thetayc~COMP, data=d1)
#W = 29714, p-value = 1.599e-15
d2<-d[d$COMP %in% c("intra_Duodenum", "inter_Duodenum"), ]
t.test(thetayc~COMP, data=d2)
wilcox.test(thetayc~COMP, data=d2)
#W = 148670, p-value = 2.509e-05
d3<-d[d$COMP %in% c("intra_Jejunum", "inter_Jejunum"), ]
t.test(thetayc~COMP, data=d3)
wilcox.test(thetayc~COMP, data=d3)
#W = 52201, p-value = 9.941e-09

# inter or intra-similarity:
d4<-d[d$comp=="intra", ]
kruskal.test(thetayc~COMP, data=d4)
#Kruskal-Wallis chi-squared = 8.4329, df = 2, p-value = 0.01475
d5<-d[d$comp=="inter", ]
kruskal.test(thetayc~COMP, data=d5)
#Kruskal-Wallis chi-squared = 647.46, df = 2, p-value < 2.2e-16

### Dunn's posthoc test (since there was significance with Kruskal-Wallis)
library(FSA)

PT = dunnTest(thetayc ~ COMP,
              data=d,
              method="bh")   
# Dunn (1964) Kruskal-Wallis multiple comparison
#  p-values adjusted with the Benjamini-Hochberg method.
  
  
```



####### Fig. S2, PCOA: files used (not used in manuscript):
	- mesal_summary.txt

```
sums<-read.table(file="mesal_summary.txt", header=TRUE, sep="\t")
sums<-sums[sums$seq_status %in% c("OK"), ]
sums$site[sums$site %in% c("Antrum")]<-"Stomach"
rownames(sums)<-sums$seqID
sums$site<-factor(sums$site, levels = c("Stomach", "Duodenum", "Proximal_Jejunum", "Mid_Jejunum", "Distal_Jejunum", "Stool"))

# method 1:
site.col <- function(n) {
colorvec <- vector(mode="character", length=length(n))
for (i in 1:length(n)) {
colorvec[i] = "light grey"
if ( n[i] == "Distal_Jejunum" ) {
colorvec[i] = "firebrick4"
}
if ( n[i] == "Duodenum" ) {
colorvec[i] = "darkblue"
}
if( n[i] == "Mid_Jejunum") {
colorvec[i] = "orange"
}
if( n[i] == "Proximal_Jejunum") {
colorvec[i] = "yellow"
}
if( n[i] == "Stomach") {
colorvec[i] = "dodgerblue1"
}
if( n[i] == "Stool") {
colorvec[i] = "sienna4"
}
}
c(colorvec)
}
plot2<-plot(data$pcoa03_axis1, data$pcoa03_axis2, ylab="Axis 2 (6.4%)", xlab="Axis1 (13.4%)", data = data, pch=21, col="black", bg=site.col(data$site), cex=1, cex.lab=1, cex.axis=1, ylim=c(-0.6, 0.6), xlim=c(-0.6, 0.6))

# OR, method 2:
# ensure that your data is in the correct order levels:
levels(data$site)
# define a color block
site.col2<-c("dodgerblue1", "darkblue", "yellow", "orange", "firebrick3", "sienna4")
plot(data$pcoa03_axis1, data$pcoa03_axis2, ylab="Axis 2 (6.4%)", xlab="Axis1 (13.4%)", pch=21, col="black", bg=site.col2[data$site], cex=1, cex.lab=1, cex.axis=1)
legend("topright",legend=gsub("_", " ", levels(data$site)), col="black", pt.bg=site.col2, cex=0.8, pch=21)

# nmds, if want that instead:
plot(data$nmds03_axis1, data$nmds03_axis2, ylab="NMDS 2", xlab="NMDS 1", pch=21, col="black", bg=site.col2[data$site], cex=1, cex.lab=1, cex.axis=1)
legend("bottomright",legend=gsub("_", " ", levels(data$site)), col="black", pt.bg=site.col2, cex=0.7, pch=21)






```