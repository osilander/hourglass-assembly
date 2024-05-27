### Figures for genome paper

library(RColorBrewer)
library(viridis)
plot.cols <- brewer.pal(9, "Set1")
fig2 <- 0 ## coverage, gc, size
fig3 <- 1 ## genome busco comparisons
fig4 <- 0 ## exon depletion
fig5 <- 0 ## phasing stats
figS1 <- 0 ## RAM time
par(mar = c(5.1, 4.1, 4.1, 2.1))
contigs <- read.table(file="./data/fig2-contig-data-q20.txt", header=T)

#################
## Fig2
## coverage, gc, size
#################
if(fig2) {
	##############
	### high quality mappers >q20
	contigs <- read.table(file="./data/fig2-contig-data-q20.txt", header=T)
	plot.col <- viridis(100)
	min.num <- min(contigs$uniq.kmer)
	max.num <- max(contigs$uniq.kmer)
	norm.num <- (contigs$uniq.kmer - min.num)/(max.num-min.num)
	col.num <- 100-floor(norm.num*length(plot.col)-1) + 1
	
	min.size <- min(log10(contigs$length))
	max.size <- max(log10(contigs$length))
	norm.size <- (log10(contigs$length) - min.size)/(max.size-min.size)

	
	pdf(file="./figures/Figure2.pdf", height=3., width=10)
	par(las=1)
	par(mar=c(5,5.2,1,1))
	par(mfrow=c(1,3))
	plot(contigs$gc, contigs$median.cov, pch=21, cex=(norm.size+0.001)*1.8, log="", xlim=c(32,55),ylim=c(0.,64.), xlab="% GC", ylab="Median depth calculated using\nonly high quality mappers", lwd=0.5, col=plot.col[col.num])
	abline(h=c(52,26), lty=2, col=plot.col[2])
	text(53,c(23.5,49.5), labels=paste("Depth=",c(26,52),sep=""), cex=0.8)

	# rectangles
	rect(rep(32.7,100), seq(4,34,length.out=101)[1:100], rep(33.2,100), seq(4,34,length.out=101)[2:101], col=rev(plot.col),border=NA)
	text(35,14,labels="% unique\n21-mers\nin contig", cex=0.7)
	text(rep(33.2,4), c(4,14,24,34), labels=c(10,40,70,100), cex=0.8, pos=2)
	
	mtext("A",side=2, line=3, adj=1.5, padj=-8., cex=1.4)
	
	#############################
	### all mappers
	contigs <- read.table(file="./data/fig2-contig-data.txt", header=T)
	plot(contigs$gc, contigs$median.cov, pch=21, cex=(norm.size+0.001)*1.8, log="y", xlim=c(32,55),ylim=c(10.,500.), xlab="% GC", ylab="Median depth calculated\nusing all mappers", lwd=0.5, col=plot.col[col.num])
	abline(h=c(52,26), lty=2, col=plot.col[2])
	text(53,c(22,44), labels=paste("Depth=",c(26,52),sep=""), cex=0.8)

	# rectangles
	rect(rep(32.7,100), 10^seq(1.9,2.7,length.out=101)[1:100], rep(33.2,100), 10^seq(1.9,2.7,length.out=101)[2:101], col=rev(plot.col),border=NA)
	text(35,10^2.3,labels="% unique\n21-mers\nin contig", cex=0.7)
	text(rep(33.2,4), 10^seq(1.9,2.7,length=4), labels=c(10,40,70,100), cex=0.8, pos=2)
	
	mtext("B",side=2, line=3, adj=1.5, padj=-8., cex=1.4)
	
	#############################
	### N regions
	contigs <- read.table(file="./data/fig2-contig-data-q20.txt", header=T)
	plot(contigs$n.ratio, contigs$median.cov, pch=21, cex=(norm.size+0.001)*1.8, xlim=c(0.3,0.8),ylim=c(0,64), xlab="Fraction of contig masked by RepeatMasker", ylab="Median depth calculated using\nonly high quality mappers", lwd=0.5, col=plot.col[col.num])
	text(53,c(22,44), labels=paste("Depth=",c(26,52),sep=""), cex=0.8)
	
	# rectangles
	rect(rep(0.32,100), seq(4,34,length.out=101)[1:100], rep(0.33,100), seq(4,34,length.out=101)[2:101], col=rev(plot.col),border=NA)
	text(0.37,14,labels="% unique\n21-mers\nin contig", cex=0.7)
	text(rep(0.33,4), c(4,14,24,34), labels=c(10,40,70,100), cex=0.8, pos=2)
	
	#rect(rep(0.31,100), seq(4,34,length.out=101)[1:100], rep(0.321,100), seq(4,34,length.out=101)[2:101], col=rev(plot.col),border=NA)
	#text(0.361,14,labels="% unique\n21-mers\nin contig", cex=0.7)
	#text(c(0.312,0.312), c(4,34), labels=c(10,100), cex=0.8, pos=4)
	
	mtext("C",side=2, line=3, adj=1.5, padj=-8., cex=1.4)
	
	dev.off()
	
	cat("Figure2\n")
}

#################
## Fig3
## busco comparisons
#################
if(fig3) {
	library(UpSetR)
	patterns <- c("miss")
	patterns.names <- c("Missing")
	#patterns <- c("miss","frag","dupl","single")
	#patterns.names <- c("Missing","Fragmented","Duplicated","Single")
	plot.list <- list()
	m.list <- list()
	for (i in 1:length(patterns)) {
		full.paths <- dir("./data/busco-data", pattern=patterns[i], full.names=T)
		files <- dir("./data/busco-data", pattern=patterns[i])
		# order better (in order of complete singles)
		correct.order <- c("aduncus","melas","albir","truncatus","obliq","orcinus","delphis", "coeru","cruciger","sinus")
		# Find the index of partial matches in the file_names vector
		correct.order.ind <- sapply(correct.order, function(x) grep(x, files))
		# Flatten the list of indices		
		correct.order.ind <- unlist(correct.order.ind)
		# Order the file_names vector based on the order of partial matches
		files <- files[correct.order.ind]
		full.paths <- full.paths[correct.order.ind]
		
		s.names <- gsub(paste("-", patterns[i],".tsv", sep=""), "", files)
		s.names  <- gsub("(^[[:alpha:]])", "\\U\\1", s.names, perl=TRUE)
		s.names <- gsub("-", " ", s.names)
		busco.list <- list()
		for (j in 1:length(files)) {
			busco.list[s.names[j]] <- read.table(full.paths[j])
		}
		plot.list[[i]] <- upset(fromList(busco.list), order.by = "freq", sets=s.names, mainbar.y.label=paste("Intersection of ", patterns.names[i]," BUSCOs",sep=""), nintersects=15, sets.x.label=paste(patterns.names[i], " BUSCOs",sep=""), keep.order=T,text.scale=1.2, set_size.numbers_size=T)
	}
	if(0) {
		### plots must be done out of the main file unknown reasons. Copy to the R command line:
		setwd("~/Documents/Manuscripts/Current/dolphin")
		# missing
		i <- 1
		pdf(file=paste("./figures/Figure3-",patterns[i],".pdf", sep=""), height=5, width=7)
		plot.list[[i]]
		dev.off()
		# single
		i <- 4
		pdf(file=paste("./figures/Figure3-",patterns[i],".pdf", sep=""), height=5, width=7)
		plot.list[[i]]
		dev.off()
		cat("Figure3\n")
	}
}


#################
## Fig4
## exon depletion
#################
if(fig4) {
	total.indels <- read.table(file="./data/indel-sizes-total.txt")
	intron.indels <- read.table(file="./data/indel-sizes-introns.txt")
	exon.indels <- read.table(file="./data/indel-sizes-exons.txt")
	
	total.indels.abs <- matrix(c(seq(1:50),total.indels[50:1,2]+total.indels[51:100,2]),nrow=50,ncol=2)
	total.indels.abs[50,2] <- 17 # missing value in introns and exons but fill here
	exon.indels.abs <- matrix(c(seq(1:50), exon.indels[50:1,2]+ exon.indels[51:100,2]),nrow=50,ncol=2)
	exon.indels.abs[50,2] <- 1 # missing value in introns and exons but fill here
	intron.indels.abs <- matrix(c(seq(1:50), intron.indels[50:1,2]+ intron.indels[51:100,2]),nrow=50,ncol=2)
	intron.indels.abs[50,2] <- 3 # missing value in introns and exons but fill here

	intron.enrich <- (intron.indels.abs[,2]/sum(intron.indels.abs[,2]))/(total.indels.abs[,2]/sum(total.indels.abs[,2]))
	exon.enrich <- (exon.indels.abs[,2]/sum(exon.indels.abs[,2]))/(total.indels.abs[,2]/sum(total.indels.abs[,2]))

	pdf(file="./figures/Figure4.pdf", height=8, width=12)
	par(las=1)
	par(mfrow=c(2,1))
	par(mar=c(2.5,7,4,1))
	
	### first the introns
	plot.col <- viridis(100)
	max.indels <- ceiling(log10(max(intron.indels.abs[,2])))
	min.indels <- floor(log10(min(intron.indels.abs[,2])))
	norm.indel <- (log10(intron.indels.abs[,2]) - min.indels)/(max.indels-min.indels)
	col.indel <- 100-floor(norm.indel*length(plot.col)-1) + 1
	col.indel[which(col.indel>100)] <- 100
	
	codons <- seq(3,48,by=3)
	barplot(intron.enrich, width=1, space=0, ylab="Intronic fold-enrichment\ncompared to identically\nsized indels genomewide",col=plot.col[col.indel], cex.lab=1.2)
	axis(1,at=c(codons)-0.5, labels=codons)
	
	legend.x <- 12
	x.width <- 12.5
	tot.labels <- 6
	x.locs <- seq(0,x.width*2,length=tot.labels) - x.width + legend.x
	legend.y <- 2.5

	
	# rectangles
	rect(seq(legend.x-x.width,legend.x+x.width,length.out=101)[1:100], rep(legend.y-0.45,100), seq(legend.x-x.width,legend.x+ x.width,length.out=101)[2:101], rep(legend.y-0.2,100), col=rev(plot.col),border=NA)
	text(legend.x,legend.y-0.6,labels="Total indels of size N", cex=1.2)
	### get a little fancy with the labels in log space
	ten.exp <- seq(min.indels, max.indels, length=tot.labels)
	ten.exp <- paste("^", ten.exp,sep="")
	# Create expressions with superscripted numbers
	legend.labels <- parse(text=paste(rep(10, tot.labels), ten.exp, sep=""))
		
	#legend.labels <- seq(min.indels, max.indels, length=tot.labels)
	text(x.locs+1, rep(legend.y-0.08,tot.labels), labels=legend.labels, cex=1.1, pos=2, offset=0.5)
	
	mtext("A",side=2, line=2, adj=2., padj=-7., cex=1.9)

	### then exons
	max.indels <- ceiling(log10(max(exon.indels.abs[,2])))
	nz.indels <- exon.indels.abs[which(exon.indels.abs[,2]>0),]
	min.indels <- floor(log10(min(nz.indels[,2])))
	norm.indel <- (log10(exon.indels.abs[,2]) - min.indels)/(max.indels-min.indels)
	col.indel <- 100-floor(norm.indel*length(plot.col)-1) + 1
	col.indel[which(col.indel>100)] <- 100

	lines(c(-2,50),c(1,1),lty=2)
	exon.enrich[22:50] <- NA
	par(mar=c(5,7,1.5,1))
	barplot(exon.enrich, width=1, space=0, beside=T, xlab="Length of indel", ylab="Exonic fold-enrichment\ncompared to identically\nsized indels genomewide",col=plot.col[col.indel], ylim=c(0,4.3), cex.lab=1.2)
	axis(1,at=c(codons)-0.5, labels=codons)
	text(c(codons-0.5), exon.enrich[codons]+0.2, labels=exon.indels.abs[codons,2], cex=1.)
	lines(c(-2,21),c(1,1),lty=2)
	
	# rectangles
	legend.x <- 35
	x.width <- 12.5
	tot.labels <- 4
	x.locs <- seq(0,x.width*2,length=tot.labels) - x.width + legend.x
	legend.y <- 3.5
	
	rect(seq(legend.x-x.width,legend.x+x.width,length.out=101)[1:100], rep(legend.y-0.6,100), seq(legend.x-x.width,legend.x+ x.width,length.out=101)[2:101], rep(legend.y-0.2,100), col=rev(plot.col),border=NA)
	
	text(legend.x,legend.y-0.85,labels="Total indels of size N", cex=1.2)
	### get a little fancy with the labels in log space
	ten.exp <- seq(min.indels, max.indels, length=tot.labels)
	ten.exp <- paste("^", ten.exp,sep="")
	legend.labels <- parse(text=paste(rep(10, tot.labels), ten.exp, sep=""))
	text(x.locs+1, rep(legend.y+0.0,tot.labels), labels=legend.labels, cex=1.1, pos=2, offset=0.5)
	
	mtext("B",side=2, line=2, adj=2., padj=-6.7, cex=1.9)
	
	dev.off()
	cat("Figure4\n")
}

#################
## Fig5
## phasing
#################
if(fig4) {
	library(viridis)
	d <- read.table(file="./phased-stats-no-all.tsv", header=T)
	
	pdf(file="./figures/Figure5.pdf", height=6, width=7)
	par(las=1)
	par(mar=c(5,5,1,1))
	par(mfrow=c(2,2))
	par(mgp=c(2.6,1,0))
	num.cats <- 8
	plot.col <- viridis(num.cats, alpha=0.4)
	
	xy.contigs <- contigs[which(contigs$median.cov<=30),]
	xy.contigs <- xy.contigs$contig[which(xy.contigs$uniq.kmer>=0.90)]
	xy.poly <- d[(d$chromosome %in% xy.contigs),]
	
	auto.contigs <- contigs[which(contigs$median.cov>=45),]
	auto.contigs <- auto.contigs$contig[which(contigs$uniq.kmer>=0.90)]
	auto.poly <- d[(d$chromosome %in% auto.contigs),]
	
	rep.contigs <- contigs[which(contigs$uniq.kmer<0.90),]
	rep.poly <- d[which(d$chromosome %in% rep.contigs$contig),]
	
	blocks <- auto.poly$blocks
	blocks[which(blocks>num.cats)] <- num.cats
	
	plot(auto.poly$contig_size, auto.poly$heterozygous_variants, log="xy", cex=0.7, pch=19, col=plot.col[blocks], xlab="Contig size", ylab="Total het. sites",xaxt="n",yaxt="n", ylim=c(3,100e3))
	
	# adjust colour scheme for xy contigs
	blocks <- xy.poly$blocks
	blocks[which(blocks>num.cats)] <- num.cats
	points(xy.poly$contig_size, xy.poly$heterozygous_variants, cex=0.7, pch=19, col=plot.col[blocks])

	# adjust colour scheme for rep contigs
	blocks <- rep.poly$blocks
	blocks[which(blocks>num.cats)] <- num.cats
	points(rep.poly$contig_size, rep.poly$heterozygous_variants, cex=0.7, pch=19, col=plot.col[blocks])

	axis(1, at=10^seq(0,7), labels=c("1","10","100","1Kbp","10Kbp","100Kbp","1Mbp", "10Mbp"))
	axis(2, at=10^seq(0,7), labels=c("1","10","100","1K","10K","100K","1M", "10M"))
	mtext("A",side=2, line=2, adj=2.5, padj=-6., cex=1.25)
	
	### legend
	### redo colours without transparency
	plot.col <- viridis(num.cats)
	rect(rep(9e6,num.cats), 10^seq(0.5,2.3,length.out=(num.cats+1))[1:num.cats], rep(14e6,num.cats), 10^seq(0.5,2.3,length.out=(num.cats+1))[2:(num.cats+1)], col=plot.col,border=NA)
	text(10e6,20,labels="Total\nphased\nblocks", cex=0.8, pos=2)	
	text(c(25e6,31e6), c(10^0.55, 10^2.15), labels=c("1",expression("">=8)),cex=0.8, pos=2)
	
	# switch back
	plot.col <- viridis(num.cats, alpha=0.4)
	# redo the plot but without the repetitive contigs
	
	blocks <- auto.poly$blocks
	blocks[which(blocks>num.cats)] <- num.cats
	
	plot(auto.poly$contig_size, auto.poly$heterozygous_variants, log="xy", cex=0.7, pch=19, col=plot.col[blocks], xlab="Contig size", ylab="Total het. sites",xaxt="n",yaxt="n",ylim=c(3,100e3))
	
	points(xy.poly$contig_size, xy.poly$heterozygous_variants, cex=0.7, pch=21, col="grey")
	
	axis(1, at=10^seq(0,7), labels=c("1","10","100","1Kbp","10Kbp","100Kbp","1Mbp", "10Mbp"))
	axis(2, at=10^seq(0,7), labels=c("1","10","100","1K","10K","100K","1M", "10M"))
	mtext("B",side=2, line=2, adj=2.5, padj=-6., cex=1.25)
	
	### Legend
	### redo colours without transparency
	plot.col <- viridis(num.cats)
	rect(rep(9e6,num.cats), 10^seq(0.5,2.3,length.out=(num.cats+1))[1:num.cats], rep(14e6,num.cats), 10^seq(0.5,2.3,length.out=(num.cats+1))[2:(num.cats+1)], col=plot.col,border=NA)
	text(10e6,20,labels="Total\nphased\nblocks", cex=0.8, pos=2)	
	text(c(25e6,31e6), c(10^0.55, 10^2.15), labels=c("1",expression("">=8)),cex=0.8, pos=2)
		
	hist(auto.poly$heterozygous_variants/auto.poly$contig_size*1e3, breaks=seq(0,20,by=0.1), xlim=c(0,5), xlab="SNPs per Kbp", main="", col=plot.col[4])
	mtext("C",side=2, line=2, adj=2, padj=-6., cex=1.25)
	
	hist(xy.poly$heterozygous_variants/xy.poly$contig_size*1e3, breaks=seq(0,20,by=0.1), xlim=c(0,5), xlab="SNPs per Kbp", main="", col=plot.col[4])
	mtext("D",side=2, line=2, adj=2, padj=-6., cex=1.25)

	dev.off()
	cat("Figure5\n")
}

#################
## FigS1
## RAM
#################
if(figS1) {
	mem <- read.table(file="./data/raven.track.txt", header=T, stringsAsFactors=F)
	mem.sub <- read.table(file="./data/raven.filtlong.track.txt", header=T, stringsAsFactors=F)
	mem$time <- as.POSIXct(mem$time, format = "%H:%M:%S")
	mem.sub$time <- as.POSIXct(mem.sub$time, format = "%H:%M:%S")
	# fix 24 clock
	mem[275:length(mem$time),1] <- mem[275:length(mem$time),1]+(24*60*60)
	mem.sub[253:length(mem.sub$time),1] <- mem.sub[253:length(mem.sub$time),1]+(24*60*60)
	
	mem$time <- mem$time-mem$time[1]
	mem.sub$time <- mem.sub$time-mem.sub$time[1]
	mem$time <- mem$time/(60*60)
	mem.sub$time <- mem.sub$time/(60*60)
	pdf(file="./figures/FigureS1.pdf", height=4, width=5)

	par(las=1)
	par(mar=c(5,6,1,2))
	plot(mem$time, mem$MEM, ty="l", xlim=c(0,39), ylim=c(0,65),xlab="Assembly runtime (hours) on 16 threads", ylab="Percent of RAM used\n(4 x 64Gb total)", lwd=2)
	lines(mem.sub$time, mem.sub$MEM, col="light blue", lwd=2)
	mem.last <- dim(mem)[1]
	mem.max <- which(mem$MEM==max(mem$MEM))
	mem.sub.max <- which(mem.sub$MEM==max(mem.sub$MEM))
	mem.sub.last <- dim(mem.sub)[1]
	text(c(mem$time[mem.last], mem.sub$time[mem.sub.last]), c(mem$MEM[mem.last], mem.sub$MEM[mem.sub.last])-0.5, pos=4, labels=c("142 Gb DNA\nsequence", "90 Gb DNA\nsequence"), cex=0.7, offset=0.1)
	
	mem.max.label <- round(mem$MEM[mem.max]*256/100,1)
	mem.sub.max.label <- round(mem.sub$MEM[mem.sub.max]*256/100,1)
	text(c(mem$time[mem.max], mem.sub$time[mem.sub.max]), c(mem$MEM[mem.max], mem.sub$MEM[mem.sub.max])+2.5, pos=4, labels=c(paste("Max. usage\n", mem.max.label, " Gb RAM", sep=""), paste("Max. usage\n", mem.sub.max.label, " Gb RAM", sep="")), cex=0.7, offset=0.1)
	
	dev.off()
	cat("FigureS1\n")
}	

