#### Mesalamine Study Manuscript:  Figure 4, streamplots
##### Anna M. Seekatz
##### 08.22.18

###### Figures created:
	- Fig. 4: streamplots of individuals (duodenum)
	- Fig S1: streamplots of individuals (stomach) 
	- Figure S2: streamplot of individual admitted 3X (duodenum and jejunum)


####### Files used:
	- mesal_genfrac2p.meta.txt
	- mesal_genfrac2p.txt


```
# subset to only the samples from p3, duodenum and jejunum:
# read in genus-level phylotype data
data<-read.table(file="mesal_genfrac2p.meta.txt", header=TRUE)
data<-data[data$seq_status %in% c("OK"), ]
genera2<-read.table(file="mesal_genfrac2p.txt", header=TRUE, sep="\t", comment.char = "")
color<-c(as.character(genera2$color), "grey47")

# patient with return visits:
levels(data$pID)
# "M047"   "M048"   "M053"   "M061"   "M062"   "M063"   "M064"
# "M046-A" "M046-B" "M046-C"

# code for streamplots:
plot.stacked <- function(x,y, ylab="Relative abundance (%)", xlab="", ncol=1, xlim=range(x, na.rm=T), ylim=c(0, 100), border = NULL, col=color){
    plot(x,y[,1], ylab=ylab, xlab=xlab, ylim=ylim, xaxs="i", yaxs="i", xlim=xlim, t="n", new=T)
    bottom=0*y[,1]
    for(i in 1:length(y[1,])){
        top=rowSums(as.matrix(y[,1:i]))
        polygon(c(x, rev(x)), c(top, rev(bottom)), border=border, col=col[i])
        bottom=top
    }
    abline(h=seq(0,200000, 10000), lty=3, col="grey")
    #legend("topleft", rev(colnames(y)), ncol=ncol, inset = 0, fill=rev(col), bty="0", bg="white", cex=0.8, col=col)
    par(new=F)
}
multi.stream.plot<- function (n) {
	for (i in length(dlist)) {
		x <- lapply(dlist, "[[", "Hour")
		y <- lapply(dlist, subset, select=c(14:ncol(dlist[[i]])))	
		y <- as.matrix(y)
		m <- nrow(y)
		n <- ncol(y)
		ph<- lapply(dlist, subset, select=c("Hour", "pH"))
		mes<- lapply(dlist, subset, select=c("Hour", "Mesalamineconcentration"))
		plot.stacked(x[[i]],y[[i]])
		title(names(dlist[i]))
			# if you want pH:
		par(new = T)
		with(ph[[i]], plot(Hour, pH, xlab=NA, ylab=NA, col="white", type="l", lwd=3, axes=F, , ylim=c(1,8.2)))	
		par(new = T)
		with(ph[[i]], plot(Hour, pH, pch=21, axes=F, xlab=NA, ylab=NA, cex=1.2, col="black", bg="white", ylim=c(1,8.2)))
		axis(side = 4)
		mtext(side = 4, line = 2, 'pH')
			# if you want mesalamine conc. in red:
		par(new = T)
		with(mes[[i]], plot(Hour, Mesalamineconcentration, xlab=NA, ylab=NA, col="red", type="l", lwd=3, axes=F, ylim=c(0,975000)))	
		par(new = T)
		with(mes[[i]], plot(Hour, Mesalamineconcentration, pch=21, axes=F, xlab=NA, ylab=NA, cex=1.2, col="black", bg="red", ylim=c(0,975000)))
		axis(side = 4, col="red", col.ticks="red", lty=1, tcl=-1, col.axis="red", line=0)
		mtext(side = 4, line = 3, 'Mesalamine', col="red")
	}
}

### before plotting, let's look at distribution of bodysites per patient, again:
sums<-data
sums$pID_site<-paste(sums$location, sums$pID, sep="_")
summary(as.factor(sums$pID_site))

# for Duodenum, graph: M047, M048, M061, M062, M063, M064
# for Stomach, graph: M048, M061, M062, M063, M064
# will also include graphs for M046, A/B/C: only stomach for C, though (include in stomach graph)

### Fig. 3: plots for 6 patients, duodenum:
par(mfrow=c(2,3))
#levels(data$pID)
# "M047"   "M048"   "M053"   "M061"   "M062"   "M063"   "M064"
# "M046-A" "M046-B" "M046-C"
# note: M053 has no Duodenal samples
d <-droplevels(data[data$location %in% c("Duodenum") & data$pID %in% c("M047"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)
d <-droplevels(data[data$location %in% c("Duodenum") & data$pID %in% c("M048"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)
d <-droplevels(data[data$location %in% c("Duodenum") & data$pID %in% c("M061"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)
d <-droplevels(data[data$location %in% c("Duodenum") & data$pID %in% c("M062"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)
d <-droplevels(data[data$location %in% c("Duodenum") & data$pID %in% c("M063"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)
d <-droplevels(data[data$location %in% c("Duodenum") & data$pID %in% c("M064"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)


### Fig. S1: plots for 6 patients, stomach:
# note: no stomach sample for M047
par(mfrow=c(2,3))
d <-droplevels(data[data$location %in% c("Stomach") & data$pID %in% c("M048"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)
d <-droplevels(data[data$location %in% c("Stomach") & data$pID %in% c("M061"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)
d <-droplevels(data[data$location %in% c("Stomach") & data$pID %in% c("M062"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)
d <-droplevels(data[data$location %in% c("Stomach") & data$pID %in% c("M063"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)
d <-droplevels(data[data$location %in% c("Stomach") & data$pID %in% c("M064"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)
d <-droplevels(data[data$location %in% c("Stomach") & data$pID %in% c("M046-C"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)


### Fig. S2: plots for patient 046, duodenum and jejunum:
# note: no stomach sample for M047
par(mfrow=c(2,3))
d <-droplevels(data[data$location %in% c("Duodenum") & data$pID %in% c("M046-A"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)
d <-droplevels(data[data$location %in% c("Duodenum") & data$pID %in% c("M046-B"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)
d <-droplevels(data[data$location %in% c("Duodenum") & data$pID %in% c("M046-C"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)
d <-droplevels(data[data$location %in% c("Jejunum") & data$pID %in% c("M046-A"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)

	# note: for the B admission in the jejunum, the mid- and distal jejunum were sampled at the same time
	# need to simplify it: let's choose mid-jejunum, since that was the site sampled multiple times
d <-droplevels(data[data$site %in% c("Mid_Jejunum") & data$pID %in% c("M046-B"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)

d <-droplevels(data[data$location %in% c("Jejunum") & data$pID %in% c("M046-C"), ])
dlist <- split(d, f=d$pID)
lapply(c("pID"), FUN=multi.stream.plot)



```