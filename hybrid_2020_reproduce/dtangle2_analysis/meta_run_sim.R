source('run_sim_fn.R')
library('cowplot')

type = "gaussian"
gen_sim=gen_sim
grd = expand.grid(
    n_markers = c(1,10,50,100,1000),
    mu=1,
    sigma=c(1/5,1/2,1,2,5)
)
print(dim(grd))
fname = paste0('./',type,'_final_wide.rds')
source('run_sim.R')

type = "nbinom"
gen_sim=gen_sim
grd = expand.grid(
    n_markers = c(1,10,50,100,1000),
    mu=1,
    sigma=c(1/5,1/2,1,2,5)
)
print(dim(grd))
fname = paste0('./',type,'_final_wide.rds')
source('run_sim.R')

type = "gaussian"
gen_sim=gen_sim2
grd = expand.grid(
    n_markers = c(1,10,50,100,1000),
    mu=1,
    sigma=c(1/5,1/2,1,2,5)
)
print(dim(grd))
fname = paste0('./',type,'_final_wide_small.rds')
source('run_sim.R')
