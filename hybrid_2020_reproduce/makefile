.PHONY: sim_plots data_plots sim_analysis data_analysis

sim_plots:
	cd ./dtangle2_analysis/ && make	plot
	mkdir -p /output/simulations
	rsync -av --include="*.pdf" --exclude="*.rds*" --exclude="*.R*" --exclude="makefile" ./dtangle2_analysis/ /output/simulations

sim_analysis:
	cd ./dtangle2_analysis/ && make analysis
	make sim_plots

data_plots:
	cd ./dtangle2_comparison/ && make plot
	mkdir -p /output/data
	rsync -av --include="*.pdf" --exclude="*.rds*" --exclude="*.R*" --exclude="makefile" ./dtangle2_comparison/ /output/data

data_analysis:
	cd ./dtangle2_comparison/ && make analysis
	make data_plots





