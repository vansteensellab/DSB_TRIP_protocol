# DSB_TRIP_protocol
git repository to accompany Schep et al. 2022

# dependencies
Before running the pipeline, some dependencies need to be installed.
Most of them can be installed using conda, others need to be downloaded from git using a bash script.   
In addition SRA toolkit needs to be installed manually:   
https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit


```bash

bash download_dependencies.sh

conda env create -f config/DSB_TRIP_env.yml

```

# example run
Before running your own analysis, you can try an example dataset.

```bash
conda activate dsb_trip

bash example.sh -t 10
```

This runs the pipeline on an example dataset using 10 cores.

After downloading the right data from SRA, the following snakemake command is called:

```bash
conda activate dsb_trip

snakemake -j 10 -s src/dsb_trip.snake \
          --configfile config/LBR_example_config.yaml \
          --config indir=example/raw outdir=example/results 
```

By replacing the --configfile and the --config arguments you can run the pipeline on your own experiment.
