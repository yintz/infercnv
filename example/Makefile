

all:
	./run.R

clean:
	rm -rf ./output_dir



debug:
	../scripts/plot_hspike.by_num_cells.R --infercnv_obj output_dir/run.final.infercnv_obj	
	../scripts/plot_hspike.R --infercnv_obj output_dir/run.final.infercnv_obj
	../scripts/run_HMM_on_hspike.R --infercnv_obj output_dir/run.final.infercnv_obj
	../scripts/plot_hspike_vs_sample_chrs.R --infercnv_obj output_dir/preliminary.infercnv_obj
