#### Mesalamine Study Manuscript:  Figure 1, community description and diversity
##### Anna M. Seekatz
##### 08.10.18

###### Figures created:
	- Fig. 1A:  barplots of mean taxonomic abundance by organ
	- Fig. 1B:  boxplots of median diversity by organ


####### Files used:
	- mesal_summary.txt
	- mesal_genfrac2p.meta.txt
	- mesal_genfrac2p.txt


```

library(plyr)

par(mfrow-c(1,2))

### Figure 1A: barplots of community membership (genus-level)

# if plotting 2 on same graph:
par(mfrow=c(2,1), xpd=T, mar=c(5,4,2,5))

### let's graph general bargraphs first:
data<-read.table(file="mesal_genfrac2p.meta.txt", header=TRUE)
genera2<-read.table(file="mesal_genfrac2p.txt", header=TRUE, sep="\t", comment.char = "")
color<-c(as.character(genera2$color), "grey47")

# get rid of sequences that failed sequencing:
data<-data[data$seq_status %in% c("OK"), ]

# let's get the mean per group--let's use site
rownames(data)<-data$seqID
datam<-droplevels(data[!is.na(data$Prevotella), c(11, 14:ncol(data))])
# note: we are combining "antrum" with "stomach" --> stomach
datam$site[datam$site=="Antrum"]<-"Stomach"
datam<-droplevels(datam)

mean.data<-ddply(datam, c("site"), colwise(mean, is.numeric))
mean.data$site <- factor(mean.data$site, levels = c("Stomach", "Duodenum", "Proximal_Jejunum", "Mid_Jejunum", "Distal_Jejunum", "Stool"))
mean.data<-mean.data[order(mean.data$site), ]
barc<-as.data.frame(mean.data[,2:ncol(mean.data)])
rownames(barc)<-mean.data[,1]

# graph:
# graphing both sets:
bar_g<-as.matrix(t(barc))
par(mar=c(5,3,2,5))
par(xpd=T)
plot1<-barplot(bar_g, main="mean phylotype proportions by body site", ylab="", ylim=c(0,100), col=color, xlim=c(0,7), axisnames=FALSE, cex.axis=0.8, width=c(0.75,0.75), space=c(0.5, 0.5))
#legend(8,100,legend=rownames(bar_g),col=color,fill=color,cex=0.5)			# will probably have to generate this separately
text(x = plot1, y = par("usr")[3]-5, srt = 45, adj = 1, labels = unique(rownames(barc)), xpd = TRUE, cex=0.8)
mtext(side=2, text="relative abundance genera (%)", line=2, cex=0.8)


### Figure 1A: barplots of community membership (genus-level)

df<-read.table(file="mesal_summary.txt", header=T)
df<-df[df$seq_status %in% c("OK"), ]
rownames(df)<-df$seqID
df$site[df$site=="Antrum"]<-"Stomach"
df<-droplevels(df)
df$site<-factor(df$site, levels = c("Stomach", "Duodenum", "Proximal_Jejunum", "Mid_Jejunum", "Distal_Jejunum", "Stool"))

# plotting by defined factors:
data<-df
div.boxplot<- function (n, m) {
				plot(data[,n] ~ as.factor(data[,m]), data = data, ylab=names(data[n]), boxwex=0.7,
					xlab="", xaxt="n", outline=FALSE, mpg=c(2,1,0), ylim=c(0,max(data[,n]*1.1, na.rm=TRUE)), cex=1, cex.lab=1, cex.axis=1)
				points(data[,n] ~ jitter(as.numeric(as.factor(data[,m]), factor=0)), data = data, col=adjustcolor("grey45", alpha=0.5), pch=21, bg="grey80", cex=0.7)
				print(as.character(unique(data[,m])))
				names<-as.character(levels(data[,m]))
				text(x =  seq(1,length(unique(data[,m])),by=1), y = par("usr")[3]-max(data[n], na.rm=TRUE)/12, srt = 45, adj = 1, labels = paste0(names, "\n n = ", as.numeric(summary(data$site))), xpd = TRUE, cex=0.8)
				stat<-kruskal.test(data[,n]~data[,m], data)
				print(stat$p.value)
				text(length(unique(data[,m]))/2, max(data[,n], na.rm=TRUE), labels=paste0("Kruskal-Wallis test \np < ", signif(stat$p.value, digits=3)), cex=0.5)
				mtext(names(data[m]), 1, line=2)
				}
par(mar=c(5,3,2,5))
lapply(c("invsimpson_03"), FUN=div.boxplot, m=c("site"))


### graph both:
par(mfrow=c(2,1), xpd=T, mar=c(5,4,2,5))
	# 1A: phylotype by body site
bar_g<-as.matrix(t(barc))
par(mar=c(5,3,2,5))
par(xpd=T)
plot1<-barplot(bar_g, main="", ylab="", ylim=c(0,100), col=color, xlim=c(0,7), axisnames=FALSE, cex.axis=0.8, width=c(0.75,0.75), space=c(0.5, 0.5))
legend(8,100,legend=rownames(bar_g),col=color,fill=color,cex=0.5)
text(x = plot1, y = par("usr")[3]-5, srt = 45, adj = 1, labels = unique(rownames(barc)), xpd = TRUE, cex=0.8)
mtext(side=2, text="relative abundance genera (%)", line=2, cex=0.8)
	# 1B: diversity by body site:
lapply(c("invsimpson_03"), FUN=div.boxplot, m=c("site"))

```
