#' Predict a candidate gene for each region conditioning on the candidates selected in rest regions.
#'
#' @param
#' pickup Numeric. Region index.
#' region_chosen A vector of TRUE/FALSE indicates which regions are selected as conditioning regions
#' remaining A vector of genes chosen in the last iteration ((n-1)th round).
#' opt A list of parameter in Gibbs
#' @return Characters of a single gene 
#' @import
#' e1071
#' @export
#' @examples
#' predict_gene(pickup,region_chosen,remaining,opt,circle)
#' 
predict_gene<-function(pickup,region_chosen,remaining,opt_tmp,circle){
    
    exclude_samegene<-T
    index<-region_chosen
    index[pickup]<-F

    cols <- opt_tmp$nodes[opt_tmp$nodes %in% unique(remaining[index])]
    rows <- opt_tmp$nodes[opt_tmp$nodes %in% unlist(opt_tmp$region[pickup])]

    ### conditional genes, exclude the same genes
    cols <- cols[!cols %in% rows]
    
    pickup_p<-opt_tmp$pro_p[opt_tmp$nodes %in% unlist(opt_tmp$region[pickup]),cols,with=FALSE]
    
    if( !is.null(dim(pickup_p)) && ncol(pickup_p)==0 )  stop("Error: no conditional genes!\n")
    
    p <- data.table(official_name = rows,p = apply(pickup_p,1,sum))

    if(opt_tmp$method != "NULL" && opt_tmp$step == "identification"){
        p1 <- opt_tmp$method_weight
        if(circle%%opt_tmp$method_interval==1 && pickup==1 && opt_tmp$verbos){
            cat("\n use method weight","\n")
            print(p1[,.SD[order(-risk)[1]],by="SNP"][order(SNP),])
        }
    } else {
        p1 <- opt_tmp$extra_weight
        if(circle%%opt_tmp$method_interval==1 && pickup==1 && opt_tmp$verbos){
            cat("\n use Extra weight","\n")
            print(p1)
        }
    }
    
    p.post <- merge(p,p1[,c("official_name","extra_weight")],by="official_name")[,p_joint:=p*extra_weight]


    ## ask quan how to combine p value
    ## avoid probabilit=0
    if(sum(p.post$p_joint)==0){
        pn <- 1/nrow(p.post)
        p.post[,p_joint:=pn]
    }

    ## if(circle%%opt_tmp$method_interval==1 && pickup==1){
    ##     cat("\n use post weight:\n")
    ##     print(p.post)
    ## }
    
    sample(as.character(p.post$official_name),1,replace=T,prob=p.post$p_joint)
 }
