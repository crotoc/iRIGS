
#' Mahalanobis transformation of evidence.
#'
#' @param dt a data.table containig the evidence data, without row names.
#' @return a data.table containing the decorred features
#' @import data.table
#' @export
#' @examples
#' weight <- mahalanobis_transformation(evi)
#' 
mahalanobis_transformation <- function(dt){
    extra <- dt[,-1]
    sigma<-cov(extra)
    s<-svd(sigma)
    s1<- s$u %*% diag((s$d)^0.5) %*% t(s$v)
    s2<- s$u %*% diag((s$d)^(-0.5)) %*% t(s$v)

    mu<-apply(extra,2,mean)
    extra_tran<-t(s2%*%(apply(extra,1,function(x) x-mu)))

    extra_p<-apply(extra_tran,2,pnorm,lower.tail=F)
    for(i in 1:ncol(extra_p))
        extra_p[,i][is.na(extra_p[,i])]<-median(extra_p[,i][!is.na(extra_p[,i])])

    ## extra_weight1<-apply(extra_p,1,function(x)  (-2*sum(log(x))))
    ## extra_weight2<-apply(extra_p,1,function(x) { pchisq(-2*sum(log(x)),df=2*length(x))})
    weight<-(-log(apply(extra_p,1,function(x) { pchisq(-2*sum(log(x)),df=2*length(x),lower.tail=F)})))
    weight<-data.table(dt[,1],extra_weight=weight)
    weight$extra_weight<-as.numeric(weight$extra_weight)
    weight[is.infinite(extra_weight),extra_weight:=max(weight[!is.infinite(extra_weight),2])]
    
    if(sum(weight$extra_weight==0)>0)
        weight[extra_weight==0,extra_weight:=min(weight[extra_weight>0,2])]

    weight
}


#v <- var(extra)
#e <- eigen(v,symmetric=TRUE)
#ginv <- e$vectors %*% diag(e$values^(-0.5)) %*% t(e$vectors)
#W <- sweep(extra,2,colMeans(extra)) %>% as.matrix

