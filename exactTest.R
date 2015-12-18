#!/usr/bin/Rscript

#load Library
library("edgeR")

# load arguments : target file with name of count file & and working directory
args<-commandArgs(TRUE)
targetFile<-args[1]
workDir<-args[2]
threshold<-args[3] # minimal number of read count in all condition to keep data in the analysis (below this value, read count is considered as noise)
pv<-args[4] #p-value, usually equal to 0.05

# Set working directory
setwd(workDir)
getwd()

# load data (count table)
targets<-read.delim(file=targetFile, header=T,comment.char="") #get target file with count table name inside
d2<-readDGE(targets, columns=c(1,2), header=F) # Create the object use for DE with edgeR, merging the count tables
colnames(d2$counts)<-targets$description # set the proper column name ( experiment condition) of the count tables
d.filtrage<-d2[rowSums(d2$counts) >= threshold,] #filter data to keep only those which are not considered as noise (usually, data with read count in all conditions > ~15)


#recalculate libraries' size after filtering
d.filtrage$samples$lib.size<-colSums(d.filtrage$counts)

#create output directory
outDir<-paste(levels(d2$samples$group)[1],"_vs_",levels(d2$samples$group)[2],sep="")
dir.create(outDir)

# set output directory
setwd(outDir)

# write excel table with filtered count
write.table(d.filtrage$counts, paste(outDir,"comptage_filtres.xls",sep=""), sep='\t', col.names=NA, row.names=T, quote=F)


#####################
# Start of analysis #
#####################

##### Edger ExactTest #####

## NORMALIZATION ##
# Test 3 normalisation type = TMM, RLE & upperquartile

# TMM
d.TMM<-d.filtrage
d.TMM<-calcNormFactors(d.TMM)

# RLE
d.RLE<-d.filtrage
d.RLE<-calcNormFactors(d.RLE)

# UpperQuartile
d.upperquartile<-d.filtrage
d.upperquartile<-calcNormFactors(d.upperquartile)

## ESTIMATE DISPERSION ##
# estimation of common dispersion because we only have 3 replicates
# tagwise dispersion fo 4 or more replicates
d.TMM<-estimateCommonDisp(d.TMM)
d.RLE<-estimateCommonDisp(d.RLE)
d.upperquartile<-estimateCommonDisp(d.upperquartile)

# write excel table with filtered & normalized matrix count
write.table(d.TMM$pseudo.counts, "TMM_counts.xls", sep='\t', col.names=NA, row.names=T, quote=F)
write.table(d.RLE$pseudo.counts, "RLE_counts.xls", sep='\t', col.names=NA, row.names=T, quote=F)
write.table(d.upperquartile$pseudo.counts, "UpperQuartile_counts.xls", sep='\t', col.names=NA, row.names=T, quote=F)


## STATISTICAL TEST ##

# Experimental conditions
# To adapt according to the experimental plan
condition1<-levels(d2$samples$group)[1]
condition2<-levels(d2$samples$group)[2]
condition1<-as.character(condition1)
condition2<-as.character(condition2)

# TEST
de.TMM<-exactTest(d.TMM, pair=c(condition2,condition1), dispersion="common")
de.RLE<-exactTest(d.RLE, pair=c(condition2,condition1), dispersion="common")
de.upperquartile<-exactTest(d.upperquartile, pair=c(condition2,condition1), dispersion="common")


# genes with p-value < pv (p-value threshold set by users in launch parameters)
sum.TMM<-summary(decideTestsDGE(de.TMM, p.value=pv))
ntop.TMM<-sum.TMM[1]+sum.TMM[3]
top_com.TMM<-topTags(de.TMM, n=ntop.TMM)

sum.RLE<-summary(decideTestsDGE(de.RLE, p.value=pv))
ntop.RLE<-sum.RLE[1]+sum.RLE[3]
top_com.RLE<-topTags(de.RLE, n=ntop.RLE)

sum.upperquartile<-summary(decideTestsDGE(de.upperquartile, p.value=pv))
ntop.upperquartile<-sum.upperquartile[1]+sum.upperquartile[3]
top_com.upperquartile<-topTags(de.upperquartile, n=ntop.upperquartile)

## MA-PLOT ##
# TMM
png( file = paste("edgeR-TMM_ExactTest_",condition1,"_vs_",condition2,"_p-value_",pv,".png",sep='\t'), height=600, width=600)
plotSmear(d.TMM, pair=c(condition2,condition1), de.tags=rownames(top_com.TMM$table), main = paste("FC plot using common dispersion ",condition1," vs ",condition2,sep=""))
abline(h=c(-1,1), col="dodgerblue")
dev.off()

# RLE
png( file = paste("edgeR-RLE_ExactTest_",condition1,"_vs_",condition2,"_p-value_",pv,".png",sep='\t'), height=600, width=600)
plotSmear(d.RLE, pair=c(condition2,condition1), de.tags=rownames(top_com.RLE$table), main = paste("FC plot using common dispersion ",condition1," vs ",condition2,sep=""))
abline(h=c(-1,1), col="dodgerblue")
dev.off()

# upperquartile
png( file = paste("edgeR-UpperQuartile_ExactTest_",condition1,"_vs_",condition2,"_p-value_",pv,".png",sep='\t'), height=600, width=600)
plotSmear(d.upperquartile, pair=c(condition2,condition1), de.tags=rownames(top_com.upperquartile$table), main = paste("FC plot using common dispersion ",condition1," vs ",condition2,sep=""))
abline(h=c(-1,1), col="dodgerblue")
dev.off()


## SCATTER PLOT ##

# TMM
jpeg( filename = paste("edgeR-TMM_",condition1,"_vs_",condition2,"_scatterplot.jpeg",sep=""), width=1300, height=900, quality=100, bg="white", res=NA)
pairs( log2(d.TMM$pseudo.counts[,1:4]), xlim=c(0,20), ylim=c(0,20), pch=46)
dev.off()

# RLE
jpeg( filename = paste("edgeR-RLE_",condition1,"_vs_",condition2,"_scatterplot.jpeg",sep=""), width=1300, height=900, quality=100, bg="white", res=NA)
pairs( log2(d.RLE$pseudo.counts[,1:4]), xlim=c(0,20), ylim=c(0,20), pch=46)
dev.off()

# UpperQuartile
jpeg( filename = paste("edgeR-UpperQuartile_",condition1,"_vs_",condition2,"_scatterplot.jpeg",sep=""), width=1300, height=900, quality=100, bg="white", res=NA)
pairs( log2(d.upperquartile$pseudo.counts[,1:4]), xlim=c(0,20), ylim=c(0,20), pch=46)
dev.off()


## TABLE OF RESULTS ##
# TMM
write.table(file = paste("edgeR-TMM_ExactTest_",condition1,"_VS_",condition2,"_p-value",pv,".xls",sep=""), x=top_com.TMM$table, sep='\t', quote=F, row.names=T, col.names=NA)

# RLE
write.table(file = paste("edgeR-RLE_ExactTest_",condition1,"_VS_",condition2,"_p-value",pv,".xls",sep=""), x=top_com.RLE$table, sep='\t', quote=F, row.names=T, col.names=NA)

#UpperQuartile
write.table(file = paste("edgeR-UpperQuartile_ExactTest_",condition1,"_VS_",condition2,"_p-value",pv,".xls",sep=""), x=top_com.upperquartile$table, sep='\t', quote=F, row.names=T, col.names=NA)



















