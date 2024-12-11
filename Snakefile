
###########################################################################
# Genome annotation pipeline using Funannotate in Snakemake
# Snakemake/5.13.0
###########################################################################
from os.path import join
from snakemake.io import expand, glob_wildcards
from Bio import SeqIO

result_dir = config["result_dir"]
input_dir = config["input_dir"]

repeat_file = config["repeat_file"]
protein_file = config["protein_file"]
transcript_file = config["transcript_file"]
alt_transcript_file = config["alt_transcript_file"]
funannotate_dir = config["funannotate_dir"]
augustus="caenorhabditis",

#SAMPLE = ["pilon124_round3_chromosomesnumbered"]
SAMPLE = list(glob_wildcards(join(input_dir, "{ids}.fasta")))[0]
#ID = glob_wildcards(join(fastq_dir, "{ids}.fastq.gz"))

print(SAMPLE)

rule All:
    input:
        # Repeats
        #expand(join(input_dir,"{samples}-families.fa"),samples=SAMPLE),
        #expand(join(result_dir,"{samples}/{samples}-families.fa"),samples=SAMPLE),
        #expand(join(result_dir,"{samples}/{samples}.cleaned.sorted.fasta.masked"),samples=SAMPLE),
        #expand(join(result_dir,"{samples}/{samples}.cleaned.sorted.fasta.out.gff"),samples=SAMPLE),

        # Funannotate - preprocessing
        #expand(join(result_dir,"{samples}/{samples}.cleaned.sorted.fasta"),samples=SAMPLE),
        #expand(join(result_dir,"{samples}/{samples}.CSM.fasta"),samples=SAMPLE),

        # Maker ctrl files
        expand(join(result_dir,"{samples}/maker_opts_rnd1.ctl"),samples=SAMPLE),
        expand(join(result_dir,"{samples}/maker_opts_rnd2.ctl"),samples=SAMPLE),
        expand(join(result_dir,"{samples}/maker_opts_rnd3.ctl"),samples=SAMPLE),
        expand(join(result_dir,"{samples}/rnd1.maker.output/rnd1_master_datastore_index.log"),samples=SAMPLE),
        expand(join(result_dir,"{samples}/rnd2.maker.output/rnd2_master_datastore_index.log"),samples=SAMPLE),
        expand(join(result_dir,"{samples}/rnd3.maker.output/rnd3_master_datastore_index.log"),samples=SAMPLE),

        # Maker GFFs
        expand(join(result_dir,"{samples}/rnd1.maker.output/rnd1.all.gff"),samples=SAMPLE),
        expand(join(result_dir,"{samples}/rnd1.maker.output/snap/rnd1.snap.hmm"),samples=SAMPLE),
        expand(join(result_dir,"{samples}/rnd2.maker.output/rnd2.all.gff"),samples=SAMPLE),
        expand(join(result_dir,"{samples}/rnd3.maker.output/rnd3.all.gff"),samples=SAMPLE),

        
#rule RepeatModeler:
#  input:
#    fa=join(input_dir,"{samples}.fasta"),
#  output:
#    fa=join(input_dir,"{samples}-families.fa"),
#    rep=join(result_dir,"{samples}/{samples}-families.fa"),
#  params:
#    rname="RepeatModeler",
#    dir=input_dir,
#    id="{samples}",
#    opdir=join(result_dir,"{samples}"),
#  threads:
#    48
#  shell:
#    """
#    mkdir -p {params.opdir}
#    cd {params.dir}
#    module load repeatmodeler
#    BuildDatabase -name {params.id} {input.fa}
#    RepeatModeler -database {params.id} -pa {threads} -LTRStruct >& {params.id}.out
#    cp {output.fa} {output.rep}
#    """

#rule RepeatMasker:
#  input:
#    fa=join(result_dir,"{samples}/{samples}.cleaned.sorted.fasta"),
#    rep=join(result_dir,"{samples}/{samples}-families.fa"),
#  output:
#    fa=join(result_dir,"{samples}/{samples}.cleaned.sorted.fasta.masked"),
#    gff=join(result_dir,"{samples}/{samples}.cleaned.sorted.fasta.out.gff"),
#  params:
#    rname="RepeatMasker",
#    dir=join(result_dir,"{samples}"),
#    #rep=repeat_file,
#    threads="48",
#  shell:
#    """
#    cd {params.dir}
#    module load repeatmasker
#    RepeatMasker -u -s -poly -engine rmblast -pa {params.threads} -gff -no_is -gccalc -norna -lib {input.rep} {input.fa}
#    """

rule maker_opts1:
    input:
        fa=join(result_dir,"{samples}/{samples}.fasta"),
    output:
        ctl=join(result_dir,"{samples}/maker_opts_rnd1.ctl"),
    params:
        rname="maker_opts1",
        protein=protein_file,
        transcript=transcript_file,
        alt_transcript=alt_transcript_file,
        augustus=augustus,
        scripts_path=join(result_dir,"param_files"),
        outdir=join(result_dir,"{samples}"),
    shell:
        """
        mkdir -p {params.outdir}
        python3 {params.scripts_path}/generate_opts1.py {output.ctl} {input.fa} {params.protein} {params.transcript} {params.augustus} {params.alt_transcript}
        """



rule maker_rnd1:
    input:
        fa=join(result_dir,"{samples}/{samples}.fasta"),
        ctl=join(result_dir,"{samples}/maker_opts_rnd1.ctl"),
    output:
        log=join(result_dir,"{samples}/rnd1.maker.output/rnd1_master_datastore_index.log"),
    params:
        rname="maker_rnd1",
        outid="rnd1",
        outdir=join(result_dir,"{samples}"),
        bopts=join(result_dir,"param_files/maker_bopts.log"),
        exe=join(result_dir,"param_files/maker_exe.log"),
    shell:
        """
        module load maker
        cd {params.outdir}
        mpiexec -np 80 maker -RM_off -base {params.outid} {input.ctl} {params.bopts} {params.exe}
        """

rule make_gff1:
    input:
        log=join(result_dir,"{samples}/rnd1.maker.output/rnd1_master_datastore_index.log"),
    output:
        gff=join(result_dir,"{samples}/rnd1.maker.output/rnd1.all.gff"),
        snap=join(result_dir,"{samples}/rnd1.maker.output/snap/rnd1.snap.hmm"),
    params:
        rname="make_gff1",
        outdir=join(result_dir,"{samples}"),
    shell:
        """
        module load maker snap
        cd {params.outdir}/rnd1.maker.output/
        gff3_merge -d {input.log}
        mkdir -p {params.outdir}/rnd1.maker.output/snap
        cd {params.outdir}/rnd1.maker.output/snap
        maker2zff -x 0.5 -l 50 -c 0 -e 0 -o 0 -d {input.log}
        fathom genome.ann genome.dna -gene-stats > gene-stats.log
        fathom genome.ann genome.dna -validate > validate.log
        fathom genome.ann genome.dna -categorize 1000 > categorize.log
        fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log
        mkdir -p params
        cd params
        forge ../export.ann ../export.dna > ../forge.log
        cd ..
        hmm-assembler.pl genome params > {output.snap}
        """

rule maker_opts2:
    input:
        fa=join(result_dir,"{samples}/{samples}.fasta"),
        gff=join(result_dir,"{samples}/rnd1.maker.output/rnd1.all.gff"),
        snap=join(result_dir,"{samples}/rnd1.maker.output/snap/rnd1.snap.hmm"),
    output:
        ctl=join(result_dir,"{samples}/maker_opts_rnd2.ctl"),
    params:
        rname="maker_opts2",
        protein=protein_file,
        transcript=transcript_file,
        scripts_path=join(result_dir,"param_files"),
        outdir=join(result_dir,"{samples}"),
    shell:
        """
        mkdir -p {params.outdir}
        python3 {params.scripts_path}/generate_opts2.py {output.ctl} {input.fa} {input.gff} {input.snap} {params.protein} {params.transcript}
        """

rule maker_rnd2:
    input:
        fa=join(result_dir,"{samples}/{samples}.fasta"),
        ctl=join(result_dir,"{samples}/maker_opts_rnd2.ctl"),
    output:
        log=join(result_dir,"{samples}/rnd2.maker.output/rnd2_master_datastore_index.log"),
    params:
        rname="maker_rnd2",
        outid="rnd2",
        outdir=join(result_dir,"{samples}"),
        bopts=join(result_dir,"param_files/maker_bopts.log"),
        exe=join(result_dir,"param_files/maker_exe.log"),
    shell:
        """
        module load maker
        cd {params.outdir}
        valgrind mpiexec -c 1 -np 40 maker -RM_off -base {params.outid} {input.ctl} {params.bopts} {params.exe}
        """

rule make_gff2:
    input:
        log=join(result_dir,"{samples}/rnd2.maker.output/rnd2_master_datastore_index.log"),
    output:
        gff=join(result_dir,"{samples}/rnd2.maker.output/rnd2.all.gff"),
    params:
        rname="make_gff2",
        outdir=join(result_dir,"{samples}"),
    shell:
        """
        module load maker snap
        cd {params.outdir}/rnd2.maker.output/
        gff3_merge -d {input.log}
        """

rule maker_rnd3:
    input:
        fa=join(result_dir,"{samples}/{samples}.fasta"),
        ctl=join(result_dir,"{samples}/maker_opts_rnd3.ctl"),
    output:
        log=join(result_dir,"{samples}/rnd3.maker.output/rnd3_master_datastore_index.log"),
    params:
        rname="maker_rnd3",
        outid="rnd2",
        outdir=join(result_dir,"{samples}"),
        bopts=join(result_dir,"param_files/maker_bopts.log"),
        exe=join(result_dir,"param_files/maker_exe.log"),
    shell:
        """
        module load maker
        cd {params.outdir}
        valgrind mpiexec -c 1 -np 40 maker -RM_off -base {params.outid} {input.ctl} {params.bopts} {params.exe}
        """

rule make_gff3:
    input:
        log=join(result_dir,"{samples}/rnd3.maker.output/rnd3_master_datastore_index.log"),
    output:
        gff=join(result_dir,"{samples}/rnd3.maker.output/rnd3.all.gff"),
    params:
        rname="make_gff3",
        outdir=join(result_dir,"{samples}"),
    shell:
        """
        module load maker snap
        cd {params.outdir}/rnd3.maker.output/
        gff3_merge -d {input.log}
        """
