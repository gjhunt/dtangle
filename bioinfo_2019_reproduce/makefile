.PHONY: plots analysis

plots:
	find ./dtangle_analysis/ -name "*.pdf" -exec rm {} \;
	find ./dtangle_analysis/ -name "*.eps" -exec rm {} \;
	cd dtangle_analysis/plotting && make
	rsync -av --exclude='*.R*' /dtangle_analysis/plotting /output

analysis:
	find ./dtangle_analysis/ -name "*.pdf" -exec rm {} \;
	find ./dtangle_analysis/ -name "*.eps" -exec rm {} \;
	cd dtangle_analysis/comparison/Analysis && rm -r *
	cd dtangle_analysis/comparison && make analysis
	cd dtangle_analysis/plotting && make
	rsync -av --exclude='*.R*' /dtangle_analysis/plotting /output
