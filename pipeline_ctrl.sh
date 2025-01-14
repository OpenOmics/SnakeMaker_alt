#!/bin/bash
set -e

module load python/3.7
module load snakemake/5.24.1

R=$2

if [ $1 == "npr" ]
then
    snakemake -npr --snakefile $R/Snakefile -j 1 --configfile $R/config.yaml
fi

if [ $1 == "process" ]
then
### WORKING
snakemake --latency-wait 120 --configfile $R/config.yaml -s $R/Snakefile -d $R --printshellcmds --use-conda --cluster-config $R/cluster.json --keep-going --cluster "sbatch --gres {cluster.gres} -n {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname}" -j 500 --rerun-incomplete --stats $R/Reports/snakemake.stats | tee -a $R/Reports/snakemake.log
fi

cd $R
mkdir -p slurmfiles

mv $R/slurm-*.out $R/slurmfiles/
if [ -f $R/HPC_usage_table.txt ]; then
	modtime1=`stat -c %y $R/HPC_usage_table.txt|awk -F "." '{print $1}'|sed 's/ /_/g' -|sed 's/:/_/g'|sed 's/-/_/g' -`
	mv $R/HPC_usage_table.txt $R/HPC_usage_table.txt.${modtime1}
fi
perl ~/Scripts/summarize_usage.pl
python ~/Scripts/filter_usage_summary.py > $R/HPCusagetable.txt.tmp
mv $R/HPCusagetable.txt.tmp $R/HPC_usage_table.txt

