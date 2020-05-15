library('data.table')

library('ggplot2')
library('cowplot')
library(gridExtra)
library(grid)

library('dtangle')
library('dtangle.data')
source('CIBERSORT_mod.R')
#devtools::install_github("Sun-lab/ICeDT")
library('ICeDT')
source('ICeDT.R')
library('DeconRNASeq')
library(MOMF)
library('EPIC')
#library('CellMix')
library('deconvSeq')
library('edgeR')
source('deconvSeq.R')
library('Rsolnp')

MTHS = c("dtangle","dtangle2","cibersort","icedt","deconrnaseq","momf","deconvseq","epic","lin","log")

deconv_run <- function(mix, ref, mrks,
                       methods=MTHS){

    print(methods)

    cat("running dtangle\n")
    dt.out = dtangle(Y=log2(1+mix),references=log2(1+ref),
                     markers = mrks)

    NA_dim <- c(nrow(mix),nrow(ref))
    NA_mtx = array(NA,NA_dim)


    dt2.phat = NA_mtx
    dt2.out = NULL    
    if("dtangle2" %in% methods){
        cat("running dtangle2\n")
        dt2.out = dtangle2(Y=mix,references=ref,
                           markers = mrks,
                           inv_scale=identity,
                           loss_smry="var",
                           dtangle_init=FALSE,
                           verbose=TRUE,
                           optim_opts=list(constr="box")
                           )
        dt2.phat = dt2.out$estimates
    }
    
    M = Reduce('c',lapply(dt.out$markers,names))

    Ym = t(mix[,M])
    Zm = t(ref[,M])

    cs.phat = NA_mtx
    if("cibersort" %in% methods){
        cat("running cibersort\n")
        cs.out = CIBERSORT(Zm,Ym)
        cs.phat = cs.out[,1:ncol(Zm)]
        rownames(cs.phat) = rownames(mix)
    }

    ice.phat = NA_mtx
    if("icedt" %in% methods){
        cat("running icedt\n")
        ice.phat = tryCatch({
            ice.out = ICeDT(Y = Ym, Z = Zm)#, tumorPurity = NULL, refVar = NULL)
            ice.phat = t(ice.out$rho)[,-1]
            ice.phat
        },error=function(e){
            print(e)
            NA_mtx
        })
    }

    drs.phat = NA_mtx
    if("deconrnaseq" %in% methods){
        cat("running deconRNAseq\n")
        drs.out = DeconRNASeq(datasets=data.frame(Ym), signatures=data.frame(Zm))
        drs.phat = drs.out$out.all
        rownames(drs.phat) = rownames(mix)
    }

    momf.phat = NA_mtx
    if("momf" %in% methods){
        cat("running momf\n")
        GList <- list(X1 = t(Zm), X2 = t(Ym))
        momf.out = momf.fit(DataX = GList,DataPriorU=Zm)
        momf.phat = momf.out$cell.prop
    }

    ds.phat = NA_mtx
    if("deconvseq" %in% methods){
        cat("running deconvSeq\n")
        ds.phat = tryCatch({
            dgeY = getdge(Ym,design=NULL)
            resultx1 = getx1.rnaseq(log(Zm),dgeY)
            resultx1$x1
        },error=function(e){
            print(e)
            array(NA,dim(cs.phat))
        })
    }

    
    #cat("running deconf")
    #deconf = t(coef(ged(object = Ym, 
    #                    x = MarkerList(markers), method = "deconf", verbose = verb)))

    #cat("running ssF")
    #ssFrobenius = t(coef(ged(object = Ym, 
    #                         x = MarkerList(markers), method = "ssFrobenius", verbose = verb)))

    #cat("running ssK")
    #ssKL = t(coef(ged(object = Ym, 
    #                         x = MarkerList(markers), method = "ssKL", verbose = verb)))

    #cat("running DSA")
    #DSA = t(coef(ged(object = Ym, 
    #                 x = MarkerList(markers), method = "DSA", verbose = verb)))

    #cat("running qprog")
    #qprog = t(coef(ged(object = Ym, 
    #                   x = Zm, method = "qprog", verbose = verb)))

    #cat("running lsfit")
    #lsfit = t(coef(ged(object = Ym, 
    #                   x = Zm, method = "lsfit", verbose = verb)))

    EPIC.phat = NA_mtx
    if("epic" %in% methods){
        cat("running epic\n")
        EPIC.phat = EPIC(bulk = Ym, reference = list(refProfiles = Zm, 
                                                     sigGenes = rownames(Zm)))$cellFractions[, 1:ncol(Zm)]
    }

    lin.phat = NA_mtx
    if("lin" %in% methods){
        cat("running lin reg\n")
        lin.phat = coef(lm(Ym~0+Zm))
        lin.phat = apply(lin.phat,c(1,2),function(x)max(0,x))
        lin.phat = t(lin.phat)/colSums(lin.phat)
        colnames(lin.phat) <- sapply(colnames(lin.phat),function(x)gsub("Zm","",x))
        lin.phat <- as.matrix(lin.phat)
        rownames(lin.phat) = rownames(mix)
    }

    log.phat = NA_mtx
    if("log" %in% methods){
        cat("running log reg\n")
        log.phat = coef(lm(log(1+Ym)~0+log(1+Zm)))
        log.phat = apply(log.phat,c(1,2),function(x)max(0,x))
        log.phat = t(log.phat)/colSums(log.phat)
        colnames(log.phat) = colnames(lin.phat)
        rownames(log.phat) = rownames(mix)
    }

    phl = list(
                dtangle=dt.out$estimates,
                    dtangle2=dt2.phat,
                cibersort=cs.phat,
                icedt=ice.phat,
                    linear=lin.phat,
                log=log.phat,
                drs=drs.phat,
        momf=momf.phat,
        ds=ds.phat,
                epic=EPIC.phat)

return(list(phats=phl,
                out=dt2.out)
           )
}

run_method = function(dset,N_markers=NULL,pct=NULL){

    ps = dset$annotation$pure_samples
    data = dset$data$log

    ref = 2^t(sapply(ps,function(x)colMeans(data[x,,drop=FALSE])))
    mix = 2^data[-unlist(ps),]

    mdf = data.table(top=apply(ref,2,which.max),
                     val=apply(ref,2,function(x)abs(diff(sort(x,decreasing=TRUE)[1:2]))),
                     g=1:ncol(ref),
                     gn=colnames(ref))
    smdf = mdf[order(top,-abs(val)),]
    mrks = lapply(unique(smdf$top),function(x){
        SD=smdf[top==x,]
        tp = SD$g
        names(tp)=SD$gn
        return(tp)
    })
    names(mrks) = rownames(ref)

    truth = as.matrix(dset$annotation$mixture[-unlist(ps),])
    rownames(truth) <- rownames(mix)

    if(is.null(N_markers))
        N_markers = max(floor(pct * ncol(ref) / nrow(ref)),1)

    print(N_markers)
    print(dset$name)
    
    mrks = lapply(mrks,head,n=N_markers)
    outd = deconv_run(mix,ref,mrks)

    err_df = function(out,truth){
        phats = out$phats
        phats = lapply(phats,as.matrix)
        #phats = phats[!sapply(phats,function(x)all(is.na(x)))]
        #phats = phats[which(sapply(sapply(phats,class),"[[",1)=="matrix")]
        edf = melt(lapply(phats,function(x)x-truth))
        colnames(edf) <- c("sample","ct","diff","method")
        tdf = melt(lapply(phats,function(x)truth))
        colnames(tdf) <- c("sample","ct","truth","method")
        adf = melt(lapply(phats,function(x)x))
        colnames(adf) <- c("sample","ct","estimate","method")
        #etdf = merge(tdf,edf,by=c('sample','ct','method'))
        #errs = data.table(merge(etdf,adf,by=c('sample','ct','method')))
        errs = cbind(edf,tdf$truth,adf$estimate)
        colnames(errs)[c(5,6)] = c('truth','estimate')
        errs=data.table(errs)
        stopifnot(all(errs$diff==-(errs$truth-errs$estimate),na.rm=TRUE))
        return(errs)
    }
    E = err_df(outd,truth)
    E$dset = dset$name
    E$N_markers = N_markers
    E[,.(median(abs(diff))),by=method]

    print("done.")
    return(E)
}

dsets = get_dtangle_data()
N_seq = c(1,2,5,10,20,50,100,200,500,1000,2000,5000)
out = lapply(N_seq,function(n)lapply(dsets,run_method,N_markers=n))
saveRDS(file="out_final.rds",object=out)
