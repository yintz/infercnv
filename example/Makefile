

all:
	./run.R

clean:
	rm -rf ./output_dir



debug:
	../scripts/plot_hspike.by_num_cells.R --infercnv_obj output_dir/run.final.infercnv_obj	
	../scripts/plot_hspike.R --infercnv_obj output_dir/run.final.infercnv_obj
	../scripts/run_HMM_on_hspike.R --infercnv_obj output_dir/run.final.infercnv_obj
	../scripts/plot_hspike_vs_sample_chrs.R --infercnv_obj output_dir/preliminary.infercnv_obj




i3:
	cat run.R  | sed s/HMM=TRUE/HMM=TRUE,HMM_type=\'i3\'/ > run.i3.R
	Rscript ./run.i3.R

Bayes:
	cat run.R  | sed s/HMM=TRUE/HMM=TRUE,BayesMaxPNormal=0.35/ > run.Bayes.R
	Rscript ./run.Bayes.R

