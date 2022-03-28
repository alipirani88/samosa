# SAMOSA - Strain Assignment from MetagenOme SNP Analysis 

SAMOSA is a Snakemake computational workflow for performing strain level detection from Metagenomic samples. The pipeline is being extensively developed to detect mixed colonization (any species) from metagenome samples.


![Workflow](https://github.com/alipirani88/samosa/blob/main/dag.svg)

### Usage

```snakemake -j 999 --cluster-config config/cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes}  -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem}" --use-conda```
