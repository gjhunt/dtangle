library('ggplot2')
library('cowplot')
library('data.table')

pdir = 'plots/'
dir.create(pdir,showWarnings=FALSE)

out = readRDS('out_final.rds')

cord <- c(14,12,7,11,1,2,4,3,15,6,8,13,10,9,5)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[cord]
}

library('dplyr')
E = lapply(1:length(out[[1]]),function(i){
    tmp = bind_rows(lapply(out,"[[",i))
    tmp$method[tmp$method=="cibersort"] <- "CIBERSORT"
    tmp$method[tmp$method=="dtangle2"] <- "Hybrid Scale"
    tmp$method[tmp$method=="linear"] <- "Regression"
    tmp$method[tmp$method=="log"] <- "LogRegression"
    tmp$method[tmp$method=="ds"] <- "deconvSeq"
    tmp$method[tmp$method=="drs"] <- "deconRNASeq"
    tmp$method[tmp$method=="icedt"] <- "ICeDT"
    tmp$method[tmp$method=="momf"] <- "MOMF"
    tmp$method[tmp$method=="epic"] <- "EPIC"
    tmp$method <- factor(tmp$method,levels=c('Hybrid Scale','Regression','ICeDT','EPIC','CIBERSORT','dtangle','MOMF','deconvSeq','deconRNASeq','LogRegression'))
    return(data.table(tmp))
})
E=E[sapply(E,nrow)>0]

newbox <- function(values) {
    #values <- na.omit(values)
    data.frame(
        ymin = min(values),
        lower = quantile(values,.25),
        middle = median(values),
        upper = quantile(values,.75),
        ymax = max(values)
    )
}

make_box_plots = function(E){
    E_cpy = E
    E_cpy$N_markers <- paste0("M = ",E_cpy$N_markers)
    E_cpy$N_markers <- ordered(E_cpy$N_markers,
                               levels=unique(E_cpy$N_markers))
    box=ggplot(data=E_cpy,mapping=aes(x=method,y=abs(diff),color=method))+stat_summary(fun.data=newbox,geom="boxplot")+theme_bw()+facet_grid(dset~N_markers)+labs(x="Method",y="|truth - estimate|")+theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="bottom")+scale_color_manual(values=gg_color_hue(15))+ggtitle(E$dset)
    box2 =ggplot(data=E,mapping=aes(x=as.factor(N_markers),y=abs(diff),color=method))+stat_summary(fun.data=newbox,geom="boxplot")+theme_bw()+facet_grid(dset~method)+labs(x="Num. Markers",y="|truth - estimate|")+theme(axis.text.x = element_text(angle = 60, hjust = 1))+guides(color=FALSE)+scale_color_manual(values=gg_color_hue(15))+ggtitle(E$dset)
    return(list(box=box,box2=box2))
}

box_plts = lapply(E,function(x)tryCatch({make_box_plots(x)},error=function(e)NULL))

pdf(paste0(pdir,"box_plots.pdf"),width=10,height=6)
for(p in box_plts){
    print(plot_grid(p$box,p$box2,nrow=2,labels='AUTO'))
}
dev.off()

make_scatter_plots=function(E){
    E_cpy = E
    E_cpy$N_markers <- paste0("N = ",E_cpy$N_markers)
    E_cpy$N_markers <- ordered(E_cpy$N_markers,
                               levels=unique(E_cpy$N_markers))
    scatter=ggplot(data=E_cpy,mapping=aes(x=truth,y=estimate,color=method,shape=ct,group=method))+geom_point()+geom_smooth(method='lm',se=FALSE)+ylim(0,1)+xlim(0,1)+geom_abline(intercept=0,slope=1)+scale_shape_manual(values=c(0:25,44-47))+theme_bw()+facet_wrap(N_markers~.)+ggtitle(unique(E$dset))
    return(scatter)
}

scatter_plts = lapply(E,make_scatter_plots)

pdf(paste0(pdir,"scatter_plots.pdf"),width=20,height=12)
for(p in scatter_plts){
    print(p)
}
dev.off()

# var by n
var_n = function(i){
    #tdf = E[[i]][,.(diff=median(abs(diff))),by=.(dset,method,N_markers)]
    tdf = E[[i]][,.(diff=median(abs(diff),na.rm=TRUE)),by=.(dset,method)]
    #tdf[,.(var=var(medErr)),by=.(method,dset)]
    tdf
}
var_dt = Reduce('rbind',lapply(1:length(E),var_n))
#var_dt=var_dt[,.(method=method,diff=diff/median(diff)),by=dset]
vdt = var_dt[,.(md=median(abs(diff))),by=method]
var_dt$method = ordered(var_dt$method,levels=vdt$method[order(vdt$md)])
lbl = var_dt[,.(med=round(median(abs(diff)),3),q3=round(quantile(abs(diff),.75,na.rm=TRUE),3),q2=round(quantile(abs(diff),.25,na.rm=TRUE),3),var=round(var(abs(diff)),3),ypos=.0035+median(abs(diff),na.rm=TRUE)),by=method]
lbl$label = paste0(lbl$med," (",lbl$var,")")
meta = ggplot(data=var_dt,mapping=aes(x=method,y=abs(diff),group=interaction(method),color=method))+stat_summary(fun.data=newbox,geom="boxplot",position=position_dodge())+geom_jitter(mapping=aes(shape=dset))+scale_shape_manual(values=0:15)
meta = meta + geom_text(data=lbl,mapping=aes(label = label, y =ypos), 
               position = position_dodge(width = .75), 
               show.legend = FALSE )+labs(y="|truth - estimate|") + 
    geom_text(data=lbl,mapping=aes(label = q3, y = .0035+q3,x=as.integer(method)+.17), 
               position = position_dodge(width = .75), 
              show.legend = FALSE )+labs(y="|truth - estimate|")+
    geom_text(data=lbl,mapping=aes(label = q2, y = -.0035+q2,x=as.integer(method)+.17), 
               position = position_dodge(width = .75), 
              show.legend = FALSE )+labs(y="|truth - estimate|")+
    scale_color_manual(values=gg_color_hue(15))+coord_cartesian(ylim=c(0,.35))+theme_bw()+labs(color="Method",shape="Dataset")
ggsave(file=paste0(pdir,"meta.pdf"),plot=meta,width=15,height=10)

var_n = function(i){
    #tdf = E[[i]][,.(diff=median(abs(diff))),by=.(dset,method,N_markers)]
    tdf = E[[i]]#[,.(diff=median(abs(diff))),by=.(dset,method)]
    #tdf[,.(var=var(medErr)),by=.(method,dset)]
    tdf
}
var_dtn = Reduce('rbind',lapply(1:length(E),var_n))
var_dtn$method = ordered(var_dtn$method,levels=vdt$method[order(vdt$md)])

meta = ggplot(data=var_dtn,mapping=aes(x=as.factor(N_markers),y=abs(diff),color=method))+stat_summary(fun.data=newbox,geom="boxplot",position=position_dodge())+scale_color_manual(values=gg_color_hue(15))+facet_grid(.~method)+theme(axis.text.x = element_text(angle = 90,line=0))+coord_cartesian(ylim=c(0,.35))
ggsave(file=paste0(pdir,"meta_byn.pdf"),plot=meta,width=15,height=10)

meta = ggplot(data=var_dtn,mapping=aes(x=as.factor(N_markers),y=abs(diff),color=method))+stat_summary(fun.data=newbox,geom="boxplot",position=position_dodge())+scale_color_manual(values=gg_color_hue(15))+coord_cartesian(ylim=c(0,.35))
ggsave(file=paste0(pdir,"meta_byn2.pdf"),plot=meta,width=15,height=10)

## Variance of estimates
var_est = var_dtn[,.(var=quantile(abs(diff),.5,na.rm=TRUE)),by=.(method,dset,N_markers)]
var_est = var_est[,.(var=var(var,na.rm=TRUE)),by=.(method,dset)]
plt = ggplot(data=var_est,mapping=aes(x=method,y=var,color=method))+stat_summary(fun.data=newbox,geom="boxplot")+geom_jitter(mapping=aes(shape=dset))+scale_shape_manual(values=0:15)+scale_color_manual(values=gg_color_hue(15))+scale_y_sqrt()+expand_limits(y=0)+labs(color="Method",shape="Dataset",y="Variance of Median Error across number of markers",x="Method")+theme_bw()
ggsave(plt,file=paste0(pdir,"var.pdf"),width=15,height=10)



# Best by dataset
best_n = function(i){
    #tdf = E[[i]][,.(medErr=median(abs(diff))),by=.(dset,method,N_markers)]
    tdf = E[[i]][,.(medErr=median(abs(diff))),by=.(dset,method,N_markers)]
    tdf[,.(best_n=N_markers[which.min(medErr)],
           worst_n=N_markers[which.max(medErr)],
           best=round(min(medErr),10),worst=round(max(medErr),10)),by=.(method,dset)]
}

best = Reduce('rbind',lapply(1:length(E),best_n))
best$method <- factor(best$method)
#levels(best$method) <- c("Our Approach","Regression","LogRegression","cibersort",'icedt')

best_plt = ggplot(data=best,mapping=aes(x=dset,y=best,fill=method))+geom_bar(stat='identity',position=position_dodge())+ggtitle("Best Error")+theme_bw()+labs(x="method",y="MAD, median(|truth - estimate|)")+theme(axis.text.x = element_text(angle = 60, hjust = 1))+scale_color_manual(values=gg_color_hue(15))
worst_plt = ggplot(data=best,mapping=aes(x=dset,y=worst,fill=method))+geom_bar(stat='identity',position=position_dodge())+ggtitle("Worst Error")+theme_bw()+labs(x="method",y="MAD, median(|truth - estimate|)")+theme(axis.text.x = element_text(angle = 60, hjust = 1))+scale_color_manual(values=gg_color_hue(15))

#plot_grid(best_plt,worst_plt)

#best$bestr = round(best$best,2)
#best$worstr = round(best$worst,2)
table(best[,.(method[best==min(best)]),by=dset]$V1)
table(best[,.(method[worst==min(worst)]),by=dset]$V1)
#table(best[,.(method[bestr==min(bestr)]),by=dset]$V1)
#table(best[,.(method[worstr==min(worstr)]),by=dset]$V1)

best$dset[best$dset=="Newman PBMC"] <- "Newman\n PBMC"
best$dset[best$dset=="Newman FL"] <- "Newman\n FL"
mbest = melt(best,id.vars=c('method','dset','best_n','worst_n'))
bw = ggplot(data=best,mapping=aes(x=method,xend=method,y=best,yend=worst,color=method))+geom_segment()+facet_grid(~dset)+theme_bw()+geom_point(data=mbest,mapping=aes(x=method,y=value,shape=variable,color=method),inherit.aes=FALSE)+theme(axis.text.x = element_text(angle = 60, hjust = 1))+labs(x="method",y="MAD, median(|truth - estimate|)")+scale_shape_manual(name="Shape",values=c(15,17),labels=c("Best","Worst"),guide="legend")+ggtitle("Best and worst cases across number of markers")+theme(legend.position="bottom",strip.text.x = element_text(angle = 0),axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_color_manual(values=gg_color_hue(15))
ggsave(paste0(pdir,"best_and_worst_by_dset.pdf"),bw,width=8,height=5)


## Meta across datasets
all_E = bind_rows(E)
all_E_smry = all_E[,.(medErr = median(abs(diff))),by=.(dset,N_markers,method)]

mbox1 = ggplot(data=all_E_smry,mapping=aes(x=as.factor(N_markers),y=medErr,color=method))+stat_summary(fun.data=newbox,geom="boxplot",position=position_dodge())+theme(axis.text.x = element_text(angle = 60, hjust = 1))+labs(x="Num. Markers",y="Median MAD, median of median(|truth - estimate|)")+scale_shape_manual(name="Shape",values=c(0,2),labels=c("Best","Worst"),guide="legend")+ggtitle("Median MAD Across all datasets")+theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust = 1))
mbox2 = ggplot(data=all_E_smry,mapping=aes(x=as.factor(N_markers),y=medErr,color=method))+stat_summary(fun.data=newbox,geom="boxplot",position=position_dodge())+theme(axis.text.x = element_text(angle = 60, hjust = 1))+labs(x="Num. Markers",y="Median MAD, median of median(|truth - estimate|)")+scale_shape_manual(name="Shape",values=c(0,2),labels=c("Best","Worst"),guide="legend")+ggtitle("Median MAD Across all datasets")+theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust = 1))+facet_grid(~method)
ggsave(paste0(pdir,"meta_boxplot.pdf"),plot_grid(mbox1,mbox2,nrow=2),width=20,height=12)


tdf = all_E_smry[,.(medErr=median(medErr)),by=.(method,N_markers)]
tdf=tdf[,.(best_n=N_markers[which.min(medErr)],
       worst_n=N_markers[which.max(medErr)],
       best=min(medErr),worst=max(medErr)),by=.(method)]
mtdf = melt(tdf,id.vars=c('method','best_n','worst_n'))


bw = ggplot(data=tdf,mapping=aes(x=method,xend=method,y=best,yend=worst,color=method))+geom_segment()+theme_bw()+geom_point(data=mtdf,mapping=aes(x=method,y=value,shape=variable),inherit.aes=FALSE)+scale_shape_manual(values=c(0,2))+theme(axis.text.x = element_text(angle = 60, hjust = 1))+labs(x="method",y="Median MAD, median of median(|truth - estimate|)")+scale_shape_manual(name="Shape",values=c(0,2),labels=c("Best","Worst"),guide="legend")+ggtitle("Meta across all datasets: Best and worst cases across number of markers")
ggsave(paste0(pdir,"meta_best_and_worst.pdf"),bw,width=20,height=12)

## best gap
best_gap = best[,.(best=round((best-min(best,na.rm=TRUE))/median(best,na.rm=TRUE),3),worst=round((worst-min(worst,na.rm=TRUE))/median(worst,na.rm=TRUE),3),method),by=dset]
best_gap$dset <- gsub("\n","",best_gap$dset)

bc = ggplot(data=best_gap,mapping=aes(x=method,y=best,fill=method))+theme_bw()+geom_bar(stat="identity")+facet_grid(~dset)+scale_fill_manual(values=gg_color_hue(15))+scale_y_sqrt(breaks=seq(0,1,.01))+theme(axis.text.x = element_text(angle = 60, hjust = 1))+ylab("Best-Case performance gap\n compared with best method")
wc = ggplot(data=best_gap,mapping=aes(x=method,y=worst,fill=method))+theme_bw()+geom_bar(stat="identity")+facet_grid(~dset)+scale_fill_manual(values=gg_color_hue(15))+scale_y_sqrt(breaks=c(seq(0,.1,.01),seq(.2,1,.1)))+theme(axis.text.x = element_text(angle = 60, hjust = 1))+ylab("Worst-Case performance gap \n compared with best method")

plt = plot_grid(plotlist=list(bc,wc),nrow=2)
ggsave(paste0(pdir,"pgap.pdf"),plt,width=10,height=10)

bc = ggplot(data=best_gap,mapping=aes(x=dset,y=best,fill=method))+theme_bw()+geom_bar(stat="identity")+facet_grid(~method)+scale_fill_manual(values=gg_color_hue(15))+scale_y_sqrt(breaks=seq(0,1,.01))+theme(axis.text.x = element_text(angle = 60, hjust = 1))+ylab("Best-Case performance gap\n compared with best method")
wc = ggplot(data=best_gap,mapping=aes(x=dset,y=worst,fill=method))+theme_bw()+geom_bar(stat="identity")+facet_grid(~method)+scale_fill_manual(values=gg_color_hue(15))+scale_y_sqrt(breaks=c(seq(0,.1,.01),seq(.2,1,.1)))+theme(axis.text.x = element_text(angle = 60, hjust = 1))+ylab("Worst-Case performance gap \n compared with best method")
plt = plot_grid(plotlist=list(bc,wc),nrow=2)
ggsave(paste0(pdir,"pgap2.pdf"),plt,width=10,height=10)
