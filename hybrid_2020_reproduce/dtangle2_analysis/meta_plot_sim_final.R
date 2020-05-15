## paper
kp_mths = c("Hybrid", 'Log Reg.','Regression','cibersort','deconvSeq')
type = "gaussian"
fl = readRDS(paste0("./", type, "_final_wide.rds"))
pdir = "gaussian_final/"
dir.create(pdir)
source("plot_sim.R")

type = "nbinom"
fl = readRDS(paste0("./", type, "_final_wide.rds"))
pdir = "nbinom_final/"
dir.create(pdir)
source("plot_sim.R")

kp_mths = c("Hybrid", 'dtangle','ICeDT')
type = "gaussian"
fl = readRDS(paste0("./", type, "_final_wide_small.rds"))
pdir = "gaussian_final_small/"
dir.create(pdir)
source("plot_sim.R")

## wide
kp_mths = c("Hybrid", 'dtangle','ICeDT','Log Reg.','Regression','cibersort','deconvSeq')
type = "gaussian"
fl = readRDS(paste0("./", type, "_final_wide.rds"))
pdir = "gaussian_wide_final/"
dir.create(pdir)
source("plot_sim.R")

type = "nbinom"
fl = readRDS(paste0("./", type, "_final_wide.rds"))
pdir = "nbinom_wide_final/"
dir.create(pdir)
source("plot_sim.R")

type = "gaussian"
fl = readRDS(paste0("./", type, "_final_wide_small.rds"))
pdir = "gaussian_wide_final_small/"
dir.create(pdir)
source("plot_sim.R")

