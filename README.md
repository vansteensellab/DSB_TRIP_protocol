# DSB_TRIP_protocol
git repository to accompany Schep et al. 2022

# install dependencies
Before running the pipeline, some dependencies need to be installed.
Most of them can be installed using conda, others need to be downloaded from git using a bash script.

```bash

bash download_dependencies.sh

conda env create -f config/DSB_TRIP_env.yml

```

# example run
Before running your own analysis, you can try an example dataset.

```bash
bash example.sh -t 10
```

This runs the pipeline on an example dataset using 10 cores.
After downloading the right data from SRA, the following snakemake command is called:

```bash

snakemake -j 10 -s src/crispr_trip.snake \
          --configfile config/LBR_example_config.yaml --use-conda \
          --config indir=example/raw outdir=example/results 
```

By replacing the --configfile and the --config arguments you can run the pipeline on your own experiment.
