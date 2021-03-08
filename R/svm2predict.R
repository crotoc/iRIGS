#' Predict a candidate gene for each region conditioning on the candidates selected in rest regions.
#'
#' @param
#' pickup Numeric. Region index.
#' region_chosen A vector of TRUE/FALSE indicates which regions are selected as conditioning regions
#' remaining A vector of genes chosen in the last iteration ((n-1)th round).
#' opt A list of parameter in Gibbs
#' @return Characters of a single gene 
#' @import
#' data.table
#' @export
#' @examples
#' predict_gene(pickup,region_chosen,remaining,opt,circle,count)
#' 
svm2predict <-function(pickup,region_chosen,remaining,opt,circle,count){
    
    exclude_samegene<-T
    index<-region_chosen
    index[pickup]<-F

    rows <- opt$nodes[opt$nodes %in% unlist(opt$region[pickup])]
    
    if(opt$svm & opt$step == "identification"){
        positive<- data.table(official_name=remaining[!remaining %in% rows],label="risk")
        negative <- data.table(official_name=sample(count[p<median(p),][!official_name %in% rows,]$official_name,size=length(remaining[index])),label="bg")
        tb <- rbind(positive,negative)
        tb$label <- as.factor(tb$label)
        training <- merge(opt$extra_processed,tb,by="official_name")

        pickup_p<-opt$pro_p[opt$nodes %in% unlist(opt$region[pickup]),cols,with=FALSE]
        ##remove consistency columns
        training[,lapply(training,function(x){ifelse(length(unique(x))==1,FALSE,TRUE)}) %>% unlist,with=F]

        
        p1 <- opt$svm_weight
        if(circle%%opt$svm_interval==1 & pickup==1){
            cat("\n use SVM weight","\n")
            print(p1[,.SD[order(-risk)[1]],by="SNP"][order(SNP),])
        }
    } else {
        p1 <- opt$extra_weight
        if(circle%%opt$svm_interval==1 & pickup==1){
            cat("\n use Extra weight","\n")
            print(p1)
        }
    }
    
    p.post <- merge(p,p1,by="official_name")[,p_joint:=p*extra_weight]

    ## avoid probabilit=0
    if(sum(p.post$p_joint)==0){
        pn <- 1/nrow(p.post)
        p.post[,p_joint:=pn]
    }

    ## if(circle%%opt$svm_interval==1 & pickup==1){
    ##     cat("\n use post weight:\n")
    ##     print(p.post)
    ## }
    
    sample(as.character(p.post$official_name),1,replace=T,prob=p.post$p_joint)
 }
