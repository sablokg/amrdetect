# amrdetect 

- complete AMR antimicrobial detection pipeline

```
# conda create -n amr-pipeline -c bioconda -c conda-forge \
#   fastp spades quast ncbi-amrfinderplus abricate -y
```

```
# Usage: ./amr_pipeline.sh -1 R1.fastq.gz -2 R2.fastq.gz -o outdir -s sample_name [-t threads]
```

Gaurav Sablok \
gsablok@proton.me