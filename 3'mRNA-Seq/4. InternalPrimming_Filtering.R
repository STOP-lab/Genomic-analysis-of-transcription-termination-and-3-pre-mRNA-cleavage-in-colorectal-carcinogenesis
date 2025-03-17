library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd <- ("3Prime-Seq/Aligned/Stranded_BAMs/")

### loading internal priming mask with regions with >= 6 consectutive As/Ts, >6 A/Ts in 10 nt window with
### exceptions of 3' gene ends (GENCODE) and experimentally detected APAs (Dirty et al. 2012)

load("/home/micgdu/Analysis/Kinga/3endQuantSec/IntPrFilter/A_mask4.RData")
load("/home/micgdu/Analysis/Kinga/3endQuantSec/IntPrFilter/T_mask4.RData")

# listing stranded Fwd/Rev bam files with alignments
Fwd<-list.files(path="3Prime-Seq/Aligned/Stranded_BAMs/", pattern="Fwd", full.names = FALSE)  
Rev<-list.files(path="3Prime-Seq/Aligned/Stranded_BAMs/", pattern="Rev", full.names = FALSE)

# terminator function:
# - extends 3'ends that represent PAs by 10 ntds
# - detects the reads that overlap (min 5 nt) with potential internal priming regions from the mask
# - retrieves non-overlapping reads => genuine PAs
# - uses functions "PAs_catcher_Fwd" & "PAs_catcher_Rev"
# - for responsive starnded bam files they retrieve: 1.) PAs 2.) PAs coverage 3.) normalised bedgraph (rpm)
# - finally all three types of files are generted for 1.) PAs 2.) filtered out signal (control of the process; may be used as quality control?)

terminator<- function(i){

		F<-readGAlignments(Fwd[i]) 
		R<-readGAlignments(Rev[i]) 

		F<-keepStandardChromosomes(F, species="Homo sapiens", pruning.mode="coarse")
		F<-dropSeqlevels(F,"chrM",pruning.mode="coarse")
		Fg<-GRanges(F)
		Fg1<-Fg
		start(Fg)<-end(Fg)
		end(Fg)<-end(Fg)+10
		Fg<-trim(Fg)

		R<-keepStandardChromosomes(R, species="Homo sapiens", pruning.mode="coarse")
		R<-dropSeqlevels(R,"chrM",pruning.mode="coarse")
		Rg<-GRanges(R)
		Rg1<-Rg
		end(Rg)<-start(Rg)
		start(Rg)<-start(Rg)-10
		Rg<-trim(Rg)

		### overlap

		olapsF<-findOverlaps(Fg, A_mask2,minoverlap=5, ignore.strand=TRUE)
		Fhits<-unique(queryHits(olapsF))
		olapsR<-findOverlaps(Rg, T_mask2,minoverlap=5, ignore.strand=TRUE)
		Rhits<-unique(queryHits(olapsR))

		### towards bw

		Findex<-1:length(Fg)
		FcleanIndex<-setdiff(Findex, Fhits)

		Rindex<-1:length(Rg)
		RcleanIndex<-setdiff(Rindex, Rhits)

		Fclean<-Fg1[FcleanIndex] # only PAs
		Rclean<-Rg1[RcleanIndex]

		F_int<-Fg1[Fhits] # only internal priming
		R_int<-Rg1[Rhits]

		Fcs<-length(FcleanIndex)/1e6
		Rcs<-length(FcleanIndex)/1e6

		Fics<-length(F_int)/1e6
		Rics<-length(R_int)/1e6

		PAs_catcher_Fwd(Fclean,Fwd[i],Fcs)
		PAs_catcher_Rev(Rclean,Rev[i],Rcs)

		PAs_catcher_Fwd(F_int,paste0(Fwd[i],"_IntPr",""),Fics)
		PAs_catcher_Rev(R_int,paste0(Rev[i],"_IntPr",""),Rics)
}

PAs_catcher_Fwd<-function(gr, name,f){
		Path<- getwd()
		grr<-gr
		start(grr)<-end(gr)-1
		save(grr, file=paste(Path,"/PAS/",name, "_raw_PAS_Fwd.RData", sep=""))
		covG<-GRanges(coverage(grr),seqinfo=seqinfo(Hsapiens))
		save(covG, file=paste(Path,"/PAS/", name, "_raw_PAS_coverage_Fwd.RData", sep=""))
		end(covG)<-end(covG)+1
		covG<-trim(covG)
		a<-as.data.frame(covG)
		a<-a[,c(1:3,6)]
		a[,2]<-format(a[,2], scientific=FALSE)
		a[,3]<-format(a[,3], scientific=FALSE)
		a[,4]<-format(a[,4]/f, scientific=FALSE)
		write.table(a ,file=paste(Path,"/PAS/", name, "_rpm_PAS_coverage_Fwd.bedGraph", sep=""),quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)
}

PAs_catcher_Rev<-function(gr, name,f){
		Path<- getwd()
		grr<-gr
		end(grr)<-start(gr)+1
		save(grr, file=paste(Path,"/PAS/",name, "_raw_PAS_Rev.RData", sep=""))
		covG<-GRanges(coverage(grr),seqinfo=seqinfo(Hsapiens))
		save(covG, file=paste(Path,"/PAS/", name, "_raw_PAS_coverage_Rev.RData", sep=""))
		end(covG)<-end(covG)+1
		covG<-trim(covG)
		a<-as.data.frame(covG)
		a<-a[,c(1:3,6)]
		a[,2]<-format(a[,2], scientific=FALSE)
		a[,3]<-format(a[,3], scientific=FALSE) ## 1 addef for proper display in UCSC
		a[,4]<-format(-a[,4]/f, scientific=FALSE)
		write.table(a ,file=paste(Path,"/PAS/", name, "_rpm_PAS_coverage_Rev.bedGraph", sep=""),quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)
}
for (i in 1:length(Fwd)) terminator(i)
