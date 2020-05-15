run:
	rm *.Rout || true
	Rscript -e "formatR::tidy_dir(recursive=TRUE,arrow=TRUE)"
	R CMD BATCH --no-save ./paper.R && tail paper.Rout -n 5
	R CMD BATCH --no-save ./log_paper.R && tail log_paper.Rout -n 5
	R CMD BATCH --no-save ./suppl_s.R && tail suppl_s.Rout -n 5
	R CMD BATCH --no-save ./sim.R && tail sim.Rout -n 5
	R CMD BATCH --no-save ./time_analysis.R && tail time_analysis.Rout -n 5
	R CMD BATCH --no-save ./newman.R && tail newman.Rout -n 5	
clean:
	rm *.Rout

