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
#' data.table
#' @export
#' @examples
#' predict_gene(pickup,region_chosen,remaining,opt,count))
#' 
svm2weight<-function(pickup,region_chosen,remaining,opt,count){
    
    exclude_samegene<-T
    index<-region_chosen
    index[pickup]<-F

    rows <- opt$nodes[opt$nodes %in% unlist(opt$region[pickup])]
    extra_processed.region <- opt$extra_processed[official_name %in% rows,]

    positive<- data.table(official_name=remaining[!remaining %in% rows],label="risk")
    negative <- data.table(official_name=sample(count[p<median(p),][!official_name %in% rows,]$official_name,size=length(remaining[index])),label="bg")
    tb <- rbind(positive,negative)
    tb$label <- as.factor(tb$label)
    
    training <- merge(opt$extra_processed,tb,by="official_name")
    ##remove consistency columns
    training[,lapply(training,function(x){ifelse(length(unique(x))==1,FALSE,TRUE)}) %>% unlist,with=F]
    model <- svm(label~.,training[,-1],probability=TRUE)
    pred <- predict(model ,extra_processed.region,probability=T)
    prob <- data.table(official_name=extra_processed.region$official_name,attr(pred, "probabilities"))
    prob[,SNP:=names(opt$region[pickup])]
    prob.label <- merge(prob,tb[,c("official_name","label")],by="official_name",all.x=T)
    #names(prob.label)[names(prob.label)=="risk"]="extra_weight"
    prob.label[,extra_weight:=risk^opt$svm_weight_power]
    ## if(sum(!is.na(p.post$label)) !=0){stop("Wrong positive and negative sets...")}
    return(prob.label)
 }
