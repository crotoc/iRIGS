#' Processing extra evidence.
#'
#' @param opt A list of parameters 
#' @return A list containing extra_processed and extra_weight used in next step
#' @import data.table
#' @export
#' @examples
#'opt$extra <- extra_evi(opt)

extra_evi<-function(opt){

###########################################################################
# if extra_evi_path and generic_evi_path are null using equal weight
###########################################################################
    if((opt$extra_evi_path=="NULL" | length(dir(opt$extra_evi_path))==0) & opt$extra_evi_file == "NULL" ){
        return(list(evi_processed=opt$generic$generic_processed,evi_weight=opt$generic$generic_weight))
    }


###########################################################################
# if extra_evi_file exist, read it and if not, read the extra_evi_dir
###########################################################################
     training_gene <- opt$training_gene$count[!is.na(grp),c("official_name","grp")]
     training_gene <- training_gene[,.SD[order(official_name)]]
     training_gene <- training_gene[!duplicated(official_name),]


    if(opt$extra_evi_file!="NULL"){
        cat("Reading extra evidance file...",opt$extra_evi_file,"\n")
        evi <- fread(opt$extra_evi_file,stringsAsFactors = FALSE)
    } else{
        cat("Reading extra evidance path...",opt$extra_evi_path,"\n")
        evi_path<-opt$extra_evi_path
        evi_file=dir(opt$extra_evi_path)[grep(opt$extra_evi_prefix,dir(opt$extra_evi_path),perl=T)]
        if(length(evi_file)==0){
            cat("No extra evi file readin !!! Use generic evi weight!!!\n")
            return(list(evi_processed=opt$generic$generic_processed,evi_weight=opt$generic$generic_weight))
        }else{
            cat("Extra evi file readin: ",evi_file,"\n")
        }
        evi.list.name <- evi_file
        evi_file<-paste(evi_path,evi_file,sep="/")
        evi.list <- lapply(evi_file, FUN=function(x){y=fread(x,stringsAsFactors = FALSE);names(y)[grep("gene_symbol",names(y),perl=T)] <- "official_name";y})
        names(evi.list) <- evi.list.name

        
        ## cat("Extra evi file filling NA: ",opt$fillNA,"\n")
        ## FUN <- match.fun(opt$fillNA)
        ## evi.list <- lapply(evi.list,function(evi){
        ##     for(i in names(evi)[-1]){
        ##         col.fill <- evi[!is.na(get(i)),get(i),] %>% FUN
        ##         evi[[i]][is.na(evi[[i]])] <- col.fill
        ##     }
        ##     evi
        ## })
        
        if(opt$method=="rf"){
            ## evi is used for random forest
            evi.list <- lapply(evi.list,rf_evi_processing)
            return(list(evi_list=evi.list))
        }
        ## evi is used for logistic regresion
        evi <- Reduce(function(x,y){merge(x,y,by="official_name")},evi.list)
    }


    cat("Extra evi file filling NA: ",opt$fillNA,"\n")
    FUN <- match.fun(opt$fillNA)
    for(i in names(evi)[-1]){
        col.fill <- evi[!is.na(get(i)),get(i),] %>% FUN
        evi[[i]][is.na(evi[[i]])] <- col.fill
    }


    
###########################################################################
# Select the training_gene in evi
###########################################################################
     
    evi.training <- merge(training_gene,evi,by="official_name",all=FALSE)
    evi.training[grp=="positive",grp:="1"]
    evi.training[grp=="negative",grp:="0"]
    evi.training[,2] <- as.numeric(evi.training[[2]])
    names(evi.training)<-gsub("greater|less","",names(evi.training))
    
###########################################################################
# Calculate the significance and direction of features using logist regression
###########################################################################

    if(opt$method != "NULL"){
    
        cat("Calculate the significance and direction of features...\n")
        
        ft <- lapply(3:ncol(evi.training),function(x){summary(glm(as.formula(paste("grp~ ","`",names(evi.training)[x],"`",sep="")),evi.training,family="binomial"))})
        ft.pvalue <- data.table(cols=names(evi.training)[-1:-2],pvalue=lapply(ft,function(x){if(dim(x$coefficients)[1]>1){x$coefficients[2,4]}else{1}}) %>% unlist, coef = lapply(ft,function(x){if(dim(x$coefficients)[1]>1){x$coefficients[2,3]}else{0}}) %>% unlist)
        ft.pvalue[pvalue<opt$p_thres,direction:=ifelse(coef>0,paste(cols,"_lmgreater",sep=""),paste(cols,"_lmless",sep=""))]
        ft.pvalue[,index:=.I+1]
        evi.select <- evi[,c(1,ft.pvalue[!is.na(direction),]$index),with=F]
        names(evi.select)[-1] <- ft.pvalue[!is.na(direction),]$direction
    }else{
        evi.select <- evi
        ft.pvalue <- 0
    }
    
    cat("Select", ncol(evi.select)-1, "features...\n")
    gene <- opt$gene[order(official_name,-distance),c("official_name","distance")]
    names(gene)[2] <- "distance_less"
    gene <- gene[!duplicated(official_name),]
    evi.select <- merge(gene,evi.select,by="official_name",all.x = TRUE)
    if(!opt$distance){
        cat("ignoring distance...\n")
        evi.select <- evi.select[,-2]
    }

    if(ncol(evi.select)!=1){
        cat("Add sign according to the direction...\n")    
        evi.select <- evi_direction(evi.select)

        cat("Extra evi.select file filling NA: ",opt$fillNA,"\n")
        for(i in names(evi.select)[-1]){
            col.fill <- evi.select[!is.na(get(i)),get(i),] %>% FUN
            evi.select[[i]][is.na(evi.select[[i]])] <- col.fill
        }

        #if only one feature is left, don't need mahananobis
        if(ncol(evi.select)==2){
            weight <- data.table(official_name=evi.select$official_name,extra_weight = pnorm(scale(evi.select[[2]])[,1],lower.tail=FALSE))
            
        
        } else {
            cat("Mahalanobis transformation...\n")
            weight <- mahalanobis_transformation(evi.select)
        }
        return(list(evi_processed=evi.select,evi_weight=weight,evi_pvalue=ft.pvalue))
    } else{
        cat("No extra evidence passed the threshold\n")
        cat("Only the generic evidence used\n")
        return(list(evi_processed=opt$generic$generic_processed,evi_weight=opt$generic$generic_weight,evi_pvalue=ft.pvalue))
    }

    
    
}


rf_evi_processing <- function(evi){
    cat("evi feature dim: ",dim(evi),"\n")
    ## evi is used for random forest
    evi <- evi[!is.na(apply(evi[,-1],1,max)),]
    evi <- evi[!duplicated(official_name),]
    cat("evi feature dim after removing NAs: ",dim(evi),"\n")
    evi <- merge(opt$gene[,"official_name"],evi,by="official_name")
    cat("select local genes: ",dim(evi),"\n")
    evi <- evi[!duplicated(official_name),]
    cat("remove duplicated genes: ",dim(evi),"\n")
    return(evi)
}    
