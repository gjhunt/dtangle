library("ggplot2")
library("dtangle.data")
library("data.table")

source("simulators.R")
source("simulators2.R")

run_sim = function(n_samples = 10, mu = 1, sigma = 1, type = "gaussian", n_markers = 100, 
    seed = NULL, gen_sim = gen_sim) {
    
    sig = paste("n_samples =", n_samples, "mu =", mu, "sigma =", sigma, "type =", 
        type, "n_markers =", n_markers, "seed = ", seed, sep = " ")
    
    sim = gen_sim(n_samples = n_samples, mu = mu, sigma = sigma, seed = seed, type = type)
    
    wnzero = which(apply(sim$dset$data$log, 2, function(x) sum(x > 0)) >= 1)
    nnzero = length(wnzero)
    
    mrks = find_markers(Y = sim$dset$data$log, pure_samples = sim$pure_samples, marker_method = "regression")
    stopifnot(all(lengths(mrks$L) >= n_markers))
    mrks = lapply(mrks$L, "[", 1:n_markers)
    stopifnot(all(lengths(mrks) <= n_markers))
    
    M = unlist(mrks)
    pure = unlist(sim$pure_samples)
    Y = sim$dset$data$count[-pure, ]
    Z = sim$dset$data$count[pure, ]
    
    snrm = sim$snr[M]
    snrm = snrm[is.finite(snrm)]
    sig = paste0(sig, "\n", "snr_markers = ", round(mean(snrm, na.rm = TRUE), 5), 
        "\n", "snr_all = ", round(mean(sim$snr, na.rm = TRUE), 5), "\n", "nnzero = ", 
        nnzero)
    
    cat("dtangle\n")
    dt_out = dtangle(Y = log2(1 + Y), references = log2(1 + Z), pure_samples = sim$pure_samples, 
        markers = mrks)
    cat("dtangle2boxvar\n")
    dt2boxvar = dtangle2(Y = log2(1 + Y), references = log2(1 + Z), pure_samples = sim$pure_samples, 
        markers = mrks, verbose = FALSE, dtangle_init = FALSE, optim_opts = list(constr = "box"), 
        loss_smry = "var")
    # cat('dtangle2eql2\n') dt2eql2 =
    # dtangle2(Y=log2(1+Y),references=log2(1+Z),pure_samples=sim$pure_samples,markers=mrks,verbose=FALSE,dtangle_init=FALSE,
    # optim_opts=list(constr='eq'),loss_smry='L2')
    
    dtE = data.frame(method = "dtangle", est = as.vector(dt_out$estimates), truth = as.vector(sim$P[-pure, 
        ]))
    dt2boxvarE = data.frame(method = "dtangle2", est = as.vector(dt2boxvar$estimates), 
        truth = as.vector(sim$P[-pure, ]))
    # dt2E =
    # data.frame(method='dtangle2',est=as.vector(dt2_out$estimates),truth=as.vector(sim$P[-pure,]))
    # dt2eql2E =
    # data.frame(method='dtangle2eq',est=as.vector(dt2eql2$diag$p_hat),truth=as.vector(sim$P[-pure,]))
    
    
    Ym = t(Y[, M])
    Zm = t(Z[, M])
    Ym[Ym <= 0] <- 0
    Zm[Zm <= 0] <- 0
    Zm = t(do.call("rbind", lapply(sim$pure_samples, function(x) rowMeans(Zm[, x, 
        drop = FALSE]))))
    
    lin_mod = lm(Ym ~ 0 + Zm)
    p_hat_lin = coef(lin_mod)
    p_hat_lin[p_hat_lin < 0] <- 0
    p_hat_lin = t(t(t(p_hat_lin)/colSums(p_hat_lin)))
    linE = data.frame(method = "Regression", est = as.vector(p_hat_lin), truth = as.vector(sim$P[-pure, 
        ]))
    
    loglin_mod = lm(log2(1 + Ym) ~ 0 + log2(1 + Zm))
    p_hat_loglin = coef(loglin_mod)
    p_hat_loglin[p_hat_loglin < 0] <- 0
    p_hat_loglin = t(t(t(p_hat_loglin)/colSums(p_hat_loglin)))
    loglinE = data.frame(method = "LogRegression", est = as.vector(p_hat_loglin), 
        truth = as.vector(sim$P[-pure, ]))
    
    source("CIBERSORT_mod.R")
    p_hat_cs = tryCatch({
        cs_mod = CIBERSORT(Zm, Ym)
        cs_mod[, 1:3]
    }, error = function(e) {
        print(e)
        array(NA, dim(p_hat_lin))
    })
    
    csE = data.frame(method = "cibersort", est = as.vector(p_hat_cs), truth = as.vector(sim$P[-pure, 
        ]))
    
    # icedt
    library("ICeDT")
    source("ICeDT.R")
    print("icedt")
    ice.phat = tryCatch({
        ice.out = ICeDT(Y = Ym, Z = Zm, tumorPurity = NULL, refVar = NULL)
        ice.phat = t(ice.out$rho)[, -1]
        ice.phat/rowSums(ice.phat)
    }, error = function(e) {
        print(e)
        array(NA, dim(p_hat_lin))
    })
    iceE = data.frame(method = "icedt", est = as.vector(ice.phat), truth = as.vector(sim$P[-pure, 
        ]))
    
    library("deconvSeq")
    source("deconvSeq.R")
    library("Rsolnp")
    ds.phat = tryCatch({
        dgeY = getdge(Ym, design = NULL)
        resultx1 = getx1.rnaseq(log(1 + Zm), dgeY)
        resultx1$x1
    }, error = function(e) {
        print(e)
        array(NA, dim(p_hat_lin))
    })
    dsE = data.frame(method = "ds", est = unlist(ds.phat), truth = as.vector(sim$P[-pure, 
        ]))
    
    constE = data.frame(method = "Constant", est = 1/ncol(sim$P), truth = as.vector(sim$P[-pure, 
        ]))
    
    E = rbind(dtE, dt2boxvarE, linE, loglinE, csE, iceE, dsE)
    
    E = data.table(E)
    E$err = abs(E$truth - E$est)
    # smry=E[,.(cor=cor(est,truth),err=median(abs(est-truth))),by=.(method)] smry
    # mod_summary = E[,.(b0 = lm(truth~1+est,data=.SD)$coef[1],b1 =
    # lm(truth~1+est,data=.SD)$coef[2]),by=.(method)] mod_summary$b1e =
    # mod_summary$b1-1 rmet = mod_summary[,.(metric = abs(b1e)+b0),by=.(method)]
    
    return(list(err = E))  #,smry=smry,M=M,snr=mean(sim$snr)))
}

