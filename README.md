# SnakeMaker
A repository for the Maker Snakemake pipeline, designed to be run on the completed assemblies, after assembly & scaffolding.


### Before running the pipeline
Edit the config.yaml as follows:

 * Update 'input_dir' to be the path where all assemblies you want to annotate are stored. Make sure that all assemblies end in '.fasta'.
 * Update 'results_dir' to the path where all results will be stored, this works best if it is this repository directory.
 * Change 'augustus_name' to the most appropriate, taken from the list from running ```funannotate species```. Set as human by default.
 * Change 'transcript_file' and 'protein_file' to the path of known transcripts/proteins for this species, or a closely related species. protein_file list can also include further away protein databases, and includes uniprot by default. These lists can have more than one file, comma separated.
 * If repeat sequences for this genome are known, change 'repeat_file' to the path with the fasta of these repeat sequences. If repeat sequences are not known, uncomment the repeatmodeler rule/output & change the input repeats for repeatmasker/fun_mask to be the repeatmodeler output.

### Test running the pipeline

The pipeline can be test run with the following command in an interactive session:

```sh pipeline_ctrl.sh npr $PWD```

Assuming that you are in the directory from this repository.

### Running the pipeline

The pipeline can be fully run with the following command:

```sbatch --time=10-00:00:00 pipeline_ctrl.sh process $PWD```

Assuming that you are in the directory from this repository.
