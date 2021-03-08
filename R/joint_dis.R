#' Return joint distribution of two genes in sampling.
#'
#' @param opt A list of parameters in Gibbs script
#' @return A data.table of coocurrence frequence
#' 
#' @import data.table reshape2
#' @export
#' @examples
#' 
joint_dis <- function(opt){
    x <- melt(opt$pg$chosen[,round:=paste("round",.I,sep="")],id.vars="round")
    x[,name:=paste(variable,value,sep="_")]
    x <- dcast(x,name~round)
    y <- opt$gene[,name:=paste(SNP,official_name,sep="_")][,"name"]
    x <- merge(y,x,by="name",all.x=TRUE)

    for(i in names(x)[-1]){
        x[!is.na(get(i)),(i):="1"]
        x[is.na(get(i)),(i):="0"]
        x[[i]]=as.numeric(x[[i]])
    }

    co <- as.matrix(x[,-1]) %*% t(as.matrix(x[,-1])) %>% data.table
    names(co) <- x[[1]]
    co
}
