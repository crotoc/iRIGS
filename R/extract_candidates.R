#' Extract candidates from snp regions.
#'
#' @param opt A list of parameters in Gibbs
#' @return A data.table containing the genes in snp regions
#' @import
#' data.table
#' magrittr
#' @importFrom GenomicRanges  mcols GRanges findOverlaps values resize width
#' @export
#' @examples
#' opt$gene<-extract_candidates(opt)

extract_candidates<-function(opt){
    gwas <- fread(opt$cand_loci,stringsAsFactors=FALSE)

    cat("Reading completed\n")
    if(is.element("identifier_pvalue",names(gwas))){
        gwas[identifier_pvalue==0,identifier_pvalue==1e-300]
    }

    ws <- opt$window_size

    if(!is.element("identifier_pvalue",colnames(gwas))){
        gwas$identifier_pvalue<-1
    }
    
    if(substr(gwas$chr[1],1,3)!="chr"){
        gwas$chr<-paste("chr",gwas$chr,sep="")
    }
    
    gwas[,gwas_start:=ifelse((pos_hg19-ws)>0,pos_hg19-ws,0)]
    gwas[,gwas_end:=pos_hg19+ws]

    cat("Convert snp to GenomicRange objects\n")
    Rle <- get("Rle", envir = asNamespace("GenomicRanges"), inherits = TRUE)
    IRanges <- get("IRanges", envir = asNamespace("GenomicRanges"), inherits = TRUE)
    gwas.gr <- GenomicRanges::GRanges(seqnames=Rle(gwas[["chr"]]),ranges = IRanges(start=gwas[["gwas_start"]],end=gwas[["gwas_end"]]))
    GenomicRanges::mcols(gwas.gr) <- gwas
    #print(gwas.gr)
    
    cat("Read gene info\n")
    gene <- fread(opt$gene_symbol,stringsAsFactors=FALSE)
    gene<-gene[!is.na(official_name),]
    gene[strand=="+",tss:=start_hg19]
    gene[strand=="-",tss:=end_hg19]
    gene$Name <- gsub("\\.\\d+",'',gene$Name,perl=T)
    gene <- gene[official_name %in% opt$nodes,]
    
    ## gene.gr <- GenomicRanges::GRanges(seqnames = Rle(gene[["chrom"]]),ranges =IRanges(gene[["start_hg19"]],end=gene[["end_hg19"]]),strand = Rle(gene[["strand"]]),ensembl=gene[["Name"]],tss=gene[["tss"]])
    gene.gr <- GenomicRanges::GRanges(seqnames = Rle(gene[["chrom"]]),ranges =IRanges(gene[["tss"]],end=gene[["tss"]]),strand = Rle(gene[["strand"]]),ensembl=gene[["Name"]],tss=gene[["tss"]])
    
    GenomicRanges::mcols(gene.gr) <- gene

    queryHits <- get("queryHits", envir = asNamespace("GenomicRanges"), inherits = TRUE)
    subjectHits <- get("subjectHits", envir = asNamespace("GenomicRanges"), inherits = TRUE)
    overlap.index <- GenomicRanges::findOverlaps(gwas.gr,gene.gr)
    output <- cbind(GenomicRanges::values(gwas.gr)[queryHits(overlap.index),],GenomicRanges::values(gene.gr[subjectHits(overlap.index),]))
    output$distance=abs(output$pos_hg19-output$tss)
    cols <- c("identifier","Name","official_name","chr","pos_hg19","identifier_pvalue","distance")
    #print(output)
    #get("as.data.frame",envir=asNamespace("GenomicRanges"),inherits=TRUE)
    output <- GenomicRanges::as.data.frame(output) %>% data.table %>% .[,..cols]
    cat("In",ws,"region of",output$identifier %>% unique %>% length,"loci ","there are",nrow(output),"genes\n",sep=" ")


####################################################################
#######		extract genes in loci regions with no genes in ws region
#####################################################################
    i=1
    while(length(gwas.gr[!gwas.gr$identifier %in% output$identifier,])!=0){
        i <- i+1
        gwas.gr.nogene <- gwas.gr[!gwas.gr$identifier %in% output$identifier,]
        gwas.gr.nogene <- GenomicRanges::resize(gwas.gr.nogene, width = 2*ws*i + width(gwas.gr.nogene), fix = "center")

        overlap.nogene.index <- GenomicRanges::findOverlaps(gwas.gr.nogene,gene.gr)
        output.nogene <- cbind(GenomicRanges::values(gwas.gr.nogene)[queryHits(overlap.nogene.index),],GenomicRanges::values(gene.gr[subjectHits(overlap.nogene.index),]))
        output.nogene$distance=abs(output.nogene$pos_hg19-output.nogene$tss)
        output.nogene <- GenomicRanges::as.data.frame(output.nogene) %>% data.table %>% .[,..cols]
        cat("In",2*ws*i,"region of",output.nogene$identifier %>% unique %>% length,"loci ","there are",nrow(output.nogene),"genes\n",sep=" ")

        if(nrow(output.nogene)!=0){
            output <- rbind(output,output.nogene)
        }
    }

    cat("In",ws,"and",ws*i,"region of",output$identifier %>% unique %>% length,"loci ","there are",nrow(output),"genes\n",sep=" ")
    
    #print(output)    
    ## output<-list()
    ## for(i in 1:nrow(gwas))
    ## {
    ##     start<-max(gwas[i,]$pos_hg19-ws,0)
    ##     end<-gwas[i,]$pos_hg19+ws
    ##     temp<-gene[gene$chrom==gwas[i,]$chr,]
    ##     index1<-temp$start_hg19<end & temp$start_hg19>start
    ##     index2<-temp$end_hg19<end & temp$end_hg19>start

    ##     if(sum(index1|index2)>0)
    ##     {
    ##         temp<-temp[index1|index2,]  
    ##         temp$identifier<-gwas[i,]$identifier
    ##         output<-rbind(output,temp)
    ##     }
    ## } 

    ## output<-output[,c("official_name","identifier")] 
    ## output$identifier_pvalue<-gwas[match(output$identifier,gwas$identifier),]$identifier_pvalue

    colnames(output)<-c("SNP","ensembl","official_name","chr","pos_hg19","SNP_pvalue","distance")  

    ## avoid same gene in same snp region but with different ensembl ID
    output <- output[!(duplicated(output[,c("SNP","official_name")])),]
    output[,index:=1:.N,by="SNP"]

    
    if(opt$max_gene != 0 && is.numeric(opt$max_gene)){
        output <- output[index<=opt$max_gene,]        
    }
    
    if(opt$evaluate_region){
        output
    } else{
        output[,-6]
    }
}


