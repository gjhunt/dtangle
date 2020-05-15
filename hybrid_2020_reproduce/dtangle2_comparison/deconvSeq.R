getx1.rnaseq =function (resultb0, dge_s, MAXITER = 1000, 
          x0 = NULL) 
{
  #if (NB0 == "all") {
  #  NB0 = length(resultb0$pfstat)
  #}else if (NB0 == "top_bonferroni") {
  #  NB0 = length(which(sort(p.adjust(resultb0$pfstat, "bonferroni")) < 
  #                       0.05))
  #}else if (NB0 == "top_fdr") {
  #  NB0 = length(which(sort(p.adjust(resultb0$pfstat, "fdr")) < 
  #                       0.05))
  #}
  b0r.1 = resultb0#$b0
  #pfstat.1 = resultb0$pfstat
  b0r = b0r.1[which(rownames(b0r.1) %in% dge_s$genes[, 1]), 
              ]
  #pfstat = pfstat.1[which(rownames(b0r.1) %in% dge_s$genes[, 
                                                        #   1])]
  #b0r = b0r[order(pfstat)[1:NB0], ]
  ntypes = ncol(b0r)
  if (is.null(x0)) 
    x0 = rep(1/ntypes, ntypes)
  offset_s = edgeR::getOffset(dge_s)
  offset_s = offset_s - mean(offset_s)
  dispersion_s.index = c()
  for (i in 1:nrow(b0r)) {
    dispersion_s.index = c(dispersion_s.index, which(dge_s$genes[, 
                                                                 1] %in% rownames(b0r)[i]))
  }
  dispersion_s = dge_s$tagwise.dispersion[dispersion_s.index]
  nsamples_s = nrow(dge_s$samples)
  x1 = matrix(0, nrow = nsamples_s, ncol = ntypes)
  one_vec = rep(1, ntypes)
  gx <- function(x) {
    one_vec %*% x - 1
  }
  converged = c()
  for (i in 1:nsamples_s) {
    print(i)
    y1 = round(dge_s$counts[rownames(b0r), i],0)
    offset_s.i = offset_s[i]
    obj <- function(x) {
      meanv = exp(b0r %*% x + offset_s.i)
      out = -sum(dnbinom(y1, size = dispersion_s, mu = meanv, log = TRUE))
      if(!is.finite(out)){
        print("not finite objective")
        return(0)
      }
      return(out)
    }
    opt_result = tryCatch(solnp(pars = x0, fun = obj, eqfun = gx, 
                                eqB = 0, LB = rep(0, ntypes), UB = rep(1, ntypes), 
                                control = list(outer.iter = MAXITER, inner.iter = MAXITER, 
                                               trace = 0)), error = function(e) {
                                                 print(i)
                                                 print(e)
                                                 NA
                                               })
    
    #print(opt_result)
    if (!is.na(unlist(opt_result)[[1]])) {
      x1[i, ] = opt_result$pars
      converged = c(converged, opt_result$convergence)
    }
    else {
      x1[i, ] = rep(NA, ntypes)
      converged = c(converged, NA)
    }
  }
  colnames(x1) = colnames(b0r)
  rownames(x1) = rownames(dge_s$samples)
  x1 = as.data.frame(x1)
  return(list(x1 = x1, converged = converged))
}