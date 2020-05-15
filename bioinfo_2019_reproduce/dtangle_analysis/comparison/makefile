SHELL := /bin/bash

analysis:
	Rscript -e "formatR::tidy_dir(recursive=TRUE,arrow=TRUE)"
	R CMD BATCH --no-save ./newman.R && tail newman.Rout -n 5	
	R CMD BATCH --no-save ./main.R && tail main.Rout -n 5
	R CMD BATCH --no-save ./ts.R && tail ts.Rout -n 5	
	R CMD BATCH --no-save ./qs.R && tail qs.Rout -n 5
	R CMD BATCH --no-save ./log_paper.R && tail log_paper.Rout -n 5
	R CMD BATCH --no-save ./ss.R && tail ss.Rout -n 5
	R CMD BATCH --no-save ./sim.R && tail sim.Rout -n 5
	R CMD BATCH --no-save ./sim_mix.R && tail sim_mix.Rout -n 5

send:
	Rscript -e "formatR::tidy_dir(recursive=TRUE,arrow=TRUE)"
	sleep 15
	~/Dropbox/util/submit main.R wt=1:00:00:00 proc=12 mem=2500mb > main.Rout 2>&1 &
	sleep 15	
	~/Dropbox/util/submit qs.R wt=1:00:00:00 proc=6 mem=10000mb > qs.Rout 2>&1 &
	sleep 15	
	~/Dropbox/util/submit shen_orr.R wt=1:00:00:00 proc=12 mem=2500mb > shen_orr.Rout 2>&1 &
	sleep 15	
	~/Dropbox/util/submit sim.R wt=1:00:00:00 proc=6 mem=10000mb > sim.Rout 2>&1 &
	sleep 15	
	~/Dropbox/util/submit sim_mix.R wt=1:00:00:00 proc=6 mem=10000mb > sim_mix.Rout 2>&1 &
	sleep 15
	~/Dropbox/util/submit ss.R wt=1:00:00:00 proc=6 mem=10000mb > ss.Rout 2>&1 &
	sleep 15	
	~/Dropbox/util/submit ts.R wt=1:00:00:00 proc=6 mem=10000mb > ts.Rout 2>&1 &
	sleep 15
	~/Dropbox/util/submit suppl.R wt=1:00:00:00 proc=12 mem=2500mb > suppl.Rout 2>&1 &
	sleep 15
	~/Dropbox/util/submit newman.R wt=1:00:00:00 proc=12 mem=2500mb > newman.Rout 2>&1 &
	sleep 15
	~/Dropbox/util/submit newman1.R wt=1:00:00:00 proc=12 mem=2500mb > newman1.Rout 2>&1 &
	sleep 15
	watch -n 1 tail -n 3 *.Rout
clean:
	rm *.Rout
