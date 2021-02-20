#' Process generic evidence and generate generic weight.
#'
#' @param opt A list of parameters in Gibbs script
#' @return A list containing generic_proccesed and generic_weight used in next step
#' @import data.table
#' @export
#' @examples
#' opt$generic <- generic_evi(opt)
 
generic_evi <- function(opt){
###########################################################################
# if generic_evi_path are null using equal weight
###########################################################################
    
    if((opt$generic_evi_path=="NULL" | length(dir(opt$generic_evi_path))==0) & opt$generic_evi_file == "NULL")
    {
        cat("No generic evi path and file used\n")
        weight<-unique(opt$gene[,"official_name"])[,extra_weight:=1]
        return(list(generic_weight=weight))
    } else if(opt$generic_evi_file!="NULL"){
        cat("READ generic evi file: ",opt$generic_evi_file,"\n")
        evi.origin <- fread(opt$generic_evi_file)
    } else{
        cat("READ generic evi path: ",opt$generic_evi_path,"\n")
        evi_file=dir(opt$generic_evi_path)[grep(opt$generic_evi_prefix,dir(opt$generic_evi_path),perl=T)]
        if(length(evi_file)==0){
            cat("No generic evi file readin !!! Use generic evi weight!!!\n")
            return(list(generic_weight=weight))
        }else{
            cat("Extra generic file readin: ",evi_file,"\n")
        }
        evi_file<-paste(opt$generic_evi_path,evi_file,sep="/")
        evi.list <- lapply(evi_file, FUN=function(x){fread(x,stringsAsFactors = FALSE)[official_name %in% gene$official_name,]})
        evi.origin <- Reduce(function(x,y){merge(x,y,by="official_name")},evi.list)
    }

    gene <- opt$gene[order(official_name,-distance),c("official_name","distance")]
    names(gene)[2] <- "distance_less"
    gene <- gene[!duplicated(official_name),]
    evi <- merge(gene,evi.origin,by="official_name",all.x = TRUE)


    if(!opt$distance){
        evi <- evi[,-2]
    }
    evi <- evi_direction(evi) 

    cat("Generic evi file filling NA: ",opt$fillNA,"\n")
    FUN <- match.fun(opt$fillNA)
    for(i in names(evi)[-1]){
        col.fill <- evi[!is.na(get(i)),get(i),] %>% FUN
        evi[[i]][is.na(evi[[i]])] <- col.fill
    }


    weight <- mahalanobis_transformation(evi)
    
    return(list(generic_processed=evi,generic_weight=weight))
}


