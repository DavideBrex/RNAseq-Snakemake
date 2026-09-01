conda activate snakemake

snakemake -j 30 -p -n --use-singularity  --singularity-args "-B /shares,/home -e"  
