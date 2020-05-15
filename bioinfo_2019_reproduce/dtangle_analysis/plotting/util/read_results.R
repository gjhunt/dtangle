pbq_res <- readRDS("../Deconvolution/Analysis/pval_breadth_quantile/pval_breadth_quantile.rds")
pbq_E <- all_errors(pbq_res)

rbq_res <- readRDS("../Deconvolution/Analysis/ratio_breadth_quantile/ratio_breadth_quantile.rds")
rbq_E <- all_errors(rbq_res)

prbq_E <- rbind(pbq_E, rbq_E)

mq_res <- readRDS("../Deconvolution/Analysis/m_quantile/m_quantile.rds")
mq_E <- all_errors(mq_res)

our_res <- readRDS("../Deconvolution/Analysis/ours_quantile/ours_quantile.rds")
our_E <- all_errors(our_res)

s_res <- readRDS("../Deconvolution/Analysis/slope/slope.rds")
s_E <- all_errors(s_res)

# paper data
paper_res <- readRDS("../Deconvolution/Analysis/paper/paper.rds")
paper_E <- all_errors(paper_res)


# get rid of this
nd_res <- readRDS("../Deconvolution/Analysis/newdata/newdata.rds")
nd_E <- all_errors(nd_res)
