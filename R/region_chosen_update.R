#' Whether a region is chosen as a conditional-on region in the next iteration. Combined the prior calculated based on a beta distribution of region pvalue and the network contribution.
#'
#' @param
#' pickup A number indicate the region index.
#' region_chosen A vector of TRUE/FALSE. The chosen region status in the last iteration
#' remaining A vector of genes. The candidates gene from last iteration.
#' opt
#' @return Logical. TRUE when the region is chosen.
#'
#' @export
#' @examples
#'  region_chosen<-unlist(foreach(pickup=1:length(opt$region)) %dopar% region_chosen_update(pickup,region_chosen,remaining,opt))
#' 
 region_chosen_update<-function(pickup,region_chosen,remaining,opt)
 {
     exclude_samegene<-T
     index<-region_chosen
     index[pickup]<-F
  
     pro_p.nodes <- names(opt$pro_p)
     cols <- pro_p.nodes[pro_p.nodes %in% unique(remaining[index])]
     rows <- pro_p.nodes[pro_p.nodes %in% unlist(opt$region[pickup])]

     ## conditional genes, exclude the same genes
     cols <- cols[!cols %in% rows]
     pickup_p<-opt$pro_p[opt$nodes %in% unlist(opt$region[pickup]),..cols]

     bg <-opt$pro_p[!opt$nodes %in% unlist(opt$region),..cols]
     bg<-bg[,p:=rowSums(.SD)][,"p"]
     bg<-bg[p>0,]

     if( !is.null(dim(pickup_p)) && ncol(pickup_p)==0 )  stop("Error: no conditional genes!\n")

     p <-data.table(official_name = rows,p = apply(pickup_p,1,sum))
     p <- merge(p,opt$extra_weight,by="official_name")[,p_joint:=p*extra_weight]

     #cat("BF MODE:",opt$bf_mode,"\n")
     if(opt$bf_mode=="quan"){
         #cat("BF MODE: AVG\n")
         p_M1 <- sum(p$p_joint)/length(p$p_joint)
         p_M0 <- (sample(bg$p,1)*sample(opt$extra_weight$extra_weight,1))
         BF <- p_M1/p_M0
     } else if(opt$bf_mode =="rui"){
         #cat("BF MODE: MAX\n")
         p_M1 <- max(p$p_joint) 
         p_M0 <- max(sample(bg$p,length(p$p_joint)) * sample(opt$extra_weight$extra_weight,length(p$p_joint)))
         BF <- p_M1/p_M0
     }

     opt$p_region <- opt$gene[,c("SNP","SNP_pvalue")][!duplicated(SNP),][,logp:=-log10(SNP_pvalue)][SNP == names(opt$region[pickup]) ,]
                                        #b<-min(-log10(snp_p$SNP_pvalue))
     b<-(-log10(0.05))
     a<-opt$p_region$logp
     prior<-beta(a+1,b)/beta(a,b+1) 

     BF<-BF*prior
     theta<- BF/(1+BF)
     
     if(rbinom(1,size=1,prob=theta)) T else F
 }
 

