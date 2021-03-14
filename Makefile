create_env:
	conda env create -f ./environment.yml

remove_env:
	conda remove --name pipeline-umi-amplicon --all

install:
	conda install -y git-lfs
	git clone -b dev https://git.oxfordnanolabs.local/research/medaka.git && cd medaka && sed -i '/mini_align/d' setup.py && pip install . && cd .. && rm -rf medaka
	pip install -e lib/
	git clone https://github.com/rrwick/Filtlong.git && cd Filtlong && make -j && rm -rf Filtlong

test:
	snakemake -pr --cores 1 all --configfile config.yml
	umi_stats -m 20 -M 60 -n EGFR_917 "example_egfr_single_read_run":example_egfr_single_read_run/ > example_egfr_single_read_run/summary.tsv
	cat example_egfr_single_read_run/summary.tsv

clean:
	rm -rf example_egfr_single_read_run/
