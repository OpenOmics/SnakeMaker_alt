# FunannotateMake
A repository for the Funannotate Snakemake pipeline


### Before running the pipeline
Edit the config.yaml as follows:

 * Update 'results_dir' to the path where all results will be stored.
 * Change 'lineage_name' to the Busco odb10 database which will be used in this analysis, chosen from https://busco-archive.ezlab.org/data/lineages/.
 * Change 'augustus_name' to the most appropriate, taken from the list from running ```funannotate species```.
 * Change 'transcript_file' and 'protein_file' to the path of known transcripts/proteins for this species, or a closely related species.
 * Change 'species' to be the species name of the assembly being annotated, or a unique identifier.
 * If repeat sequences for this genome are known, change 'repeat_file' to the path with the fasta of these repeat sequences. If repeat sequences are not known, uncomment the repeatmodeler rule/output & change the input repeats for repeatmasker/fun_mask to be the repeatmodeler output.
 * Change 'funannotate_dir' to be the path where busco odb10 databases will be stored.
