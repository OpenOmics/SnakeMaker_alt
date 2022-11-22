
###########################################################################
# Genome annotation pipeline using Funannotate in Snakemake
# Snakemake/5.13.0
###########################################################################
from os.path import join
from snakemake.io import expand, glob_wildcards

result_dir = config["result_dir"]
Lineage_name = config["lineage_name"]
busco_db = Lineage_name.split("_odb")[0]
augustus_name = config["augustus_name"]

repeat_file = config["repeat_file"]
protein_file = config["protein_file"]
transcript_file = config["transcript_file"]
funannotate_dir = config["funannotate_dir"]
species = config["species"]
species_id = species.replace(" ","_")


SAMPLE = ["pilon124_round3_chromosomesnumbered"]
ID = glob_wildcards(join(fastq_dir, "{ids}.fastq.gz"))

print(SAMPLE)

print(ID)

rule All:
    input:
        # Repeats
        #expand(join(result_dir,"funannotate/{samples}/{samples}-families.fa"),samples=SAMPLE),
        expand(join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta.masked"),samples=SAMPLE),
        expand(join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta.out.gff"),samples=SAMPLE),

        # Funannotate - preprocessing
        join(funannotate_dir,busco_db + ".tar.gz")
        expand(join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta"),samples=SAMPLE),
        expand(join(result_dir,"funannotate/{samples}/{samples}.CSM.fasta"),samples=SAMPLE),

        # Funannotate - predict
        expand(join(result_dir,"funannotate/{samples}/fun/predict_results/"+species_id+".gff3"),samples=SAMPLE),
       	expand(join(result_dir,"funannotate/{samples}/fun/predict_results/"+species_id+".proteins.fa"),samples=SAMPLE),
       	expand(join(result_dir,"funannotate/{samples}/fun/predict_results/"+species_id+".mrna-transcripts.fa"),samples=SAMPLE),

        # Funannotate - annotate
        expand(join(result_dir,"funannotate/{samples}/fun/annotate_results/"+species_id+".gff3"),samples=SAMPLE),
        #expand(join(result_dir,"funannotate/{samples}/fun/annotate_misc/iprscan.xml"),samples=SAMPLE),
        #expand(join(result_dir,"funannotate/{samples}/fun/train_results/"+species_id+".gff3"),samples=SAMPLE),

rule fun_setup:
    input:
        fa=transcript_file,
    output:
        dbout=join(funannotate_dir,busco_db + ".tar.gz")
    params:
        fun_dir=funannotate_dir,
        busco_db=busco_db,
    shell:
        """
        module load funannotate/1.8.9
        funannotate setup -b {params.busco_db} -d {params.fun_dir}
        """
        
#rule RepeatModeler:
#  input:
#    fa=join(result_dir,"all-assemblies/{samples}.fasta"),
#  output:
#    fa=join(result_dir,"all-assemblies/{samples}-families.fa"),
#  params:
#    rname="RepeatModeler",
#    dir=join(result_dir,"all-assemblies"),
#    id="{samples}.{assemblers}"
#  threads:
#    48
#  shell:
#    """
#    cd {params.dir}
#    module load repeatmodeler
#    BuildDatabase -name {params.id} {input.fa}
#    RepeatModeler -database {params.id} -pa {threads} -LTRStruct >& {params.id}.out
#    """

rule RepeatMasker:
  input:
    fa=join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta"),
  output:
    fa=join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta.masked"),
    gff=join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta.out.gff"),
  params:
    rname="RepeatMasker",
    dir=join(result_dir,"funannotate/{samples}"),
    rep=repeat_file,
    threads="48",
  shell:
    """
    cd {params.dir}
    module load repeatmasker
    RepeatMasker -u -s -poly -engine rmblast -pa {params.threads} -gff -no_is -gccalc -norna -lib {params.rep} {input.fa}
    """

rule fun_clean:
    input:
        fa=join(result_dir,"hifiasm/{samples}.ragtag.fasta"),
    output:
        fa=temp(join(result_dir,"funannotate/{samples}/{samples}.cleaned.fasta")),
    params:
        rname="fun_clean",
        dir=join(result_dir,"funannotate/{samples}"),
    run:
        from Bio import SeqIO
        out = open(output.fa,"w")
        for r in SeqIO.parse(open(input.fa,"r"),"fasta"):
            if len(str(r.seq)) >= 1000:
                out.write(">" + r.description + "\n" + str(r.seq))

rule fun_sort:
    input:
        fa=join(result_dir,"funannotate/{samples}/{samples}.cleaned.fasta"),
    output:
        fa=join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta"),
    params:
        rname="fun_sort",
        dir=join(result_dir,"funannotate/{samples}"),
    shell:
        """
        module load funannotate
        mkdir -p {params.dir}
        funannotate sort -i {input.fa} -b scaffold -o {output.fa}
        """

rule fun_mask:
    input:
        fa=join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta"),
    output:
        fa=join(result_dir,"funannotate/{samples}/{samples}.CSM.fasta"),
    params:
        rname="fun_mask",
        rep=repeat_file,
        dir=join(result_dir,"funannotate/{samples}"),
        threads="32",
    shell:
        """
        module load funannotate
        mkdir -p {params.dir}
        funannotate mask -i {input.fa} --cpus {params.threads} -o {output.fa} -l {params.rep}
        """

rule fun_predict:
    input:
        fa=join(result_dir,"funannotate/{samples}/{samples}.CSM.fasta"),
    output:
        gff=join(result_dir,"funannotate/{samples}/fun/predict_results/"+species_id+".gff3"),
        prot=join(result_dir,"funannotate/{samples}/fun/predict_results/"+species_id+".proteins.fa"),
        tran=join(result_dir,"funannotate/{samples}/fun/predict_results/"+species_id+".mrna-transcripts.fa"),
    params:
        rname="fun_predict",
        dir=join(result_dir,"funannotate/{samples}/fun"),
        protein=protein_file,
        transcript=transcript_file,
        species=species,
        busco=augustus_name,
        fun_dir=funannotate_dir,
        busco_db=busco_db,
        threads="48",
    shell:
        """
        module load funannotate/1.8.9
        mkdir -p {params.dir}
        funannotate predict -i {input.fa} -o {params.dir} -d {params.fun_dir} --busco_db {params.busco_db} \
        --species "{params.species}" --busco_seed_species {params.busco} --cpus {params.threads} --organism other \
        --protein_evidence {params.protein} --transcript_evidence {params.transcript}
        """

#rule fun_iprs:
#    input:
#        gff=join(result_dir,"funannotate/{samples}/fun/predict_results/"+species_id+".gff3"),
#        prot=join(result_dir,"funannotate/{samples}/fun/predict_results/"+species_id+".proteins.fa"),
#        tran=join(result_dir,"funannotate/{samples}/fun/predict_results/"+species_id+".mrna-transcripts.fa"),
#    output:
#        xml=join(result_dir,"funannotate/{samples}/fun/annotate_misc/iprscan.xml"),
#    params:
#        rname="fun_iprs",
#        dir=join(result_dir,"funannotate/{samples}/fun"),
#        iprpath="/usr/local/apps/interproscan/5.57-90.0/interproscan_app/interproscan.sh",
#        threads="32"
#    shell:
#        """
#        module load funannotate/1.8.9 interproscan
#        funannotate iprscan -c {params.threads} -i {params.dir} -m local --iprscan_path {params.iprpath}
#        """

rule fun_annotate:
    input:
        fasta=join(result_dir,"funannotate/{samples}/{samples}.CSM.fasta"),
        gff=join(result_dir,"funannotate/{samples}/fun/predict_results/"+species_id+".gff3"),
        prot=join(result_dir,"funannotate/{samples}/fun/predict_results/"+species_id+".proteins.fa"),
        tran=join(result_dir,"funannotate/{samples}/fun/predict_results/"+species_id+".mrna-transcripts.fa"),
    output:
        gff=join(result_dir,"funannotate/{samples}/fun/annotate_results/"+species_id+".gff3"),
    params:
        fun_dir=funannotate_dir,
        rname="fun_annotate",
        dir=join(result_dir,"funannotate/{samples}/fun"),
        species=species,
        busco_db=busco_db,
        threads="32",
    shell:
        """
        module load funannotate/1.8.9
        mkdir -p {params.dir}
        funannotate annotate -i {params.dir} --cpus {params.threads} --busco_db {params.busco_db} --fasta {input.fasta} --gff {input.gff} -s "{params.species}" -d {params.fun_dir}
        """

#rule fun_train:
#    input:
#        fa=join(result_dir,"funannotate/{samples}/{samples}.CSM.fasta"),
#        R1=expand(join(fastq_dir,"{ids}_1.fastq.gz"),ids=ID),
#        R2=expand(join(fastq_dir,"{ids}_2.fastq.gz"),ids=ID),
#    output:
#        fa=join(result_dir,"funannotate/{samples}/fun/train_results/"+species_id+".gff3"),
#    params:
#        rname="fun_train",
#        dir=join(result_dir,"funannotate/{samples}/fun"),
#        species=species,
#        threads="48",
#    shell:
#        """
#        module load funannotate trinity minimap2 samtools
#        mkdir -p {params.dir}
#        funannotate train -i {input.fa} -o {params.dir} \
#        --left {input.R1} \
#        --right {input.R2} \
#        --stranded RF --cpus {params.threads} --species "{params.species}" --TRINITYHOME .
#        """

