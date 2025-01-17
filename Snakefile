###########################################################################
# Genome annotation pipeline using Maker in Snakemake
# Snakemake/5.13.0
###########################################################################
from os.path import join
from snakemake.io import expand, glob_wildcards
from Bio import SeqIO

# Load configuration paths and files
result_dir = config["result_dir"]  # Directory to store results
input_dir = config["input_dir"]   # Directory containing input files
repeat_file = config["repeat_file"]  # File containing repeat sequences
protein_file = config["protein_file"]  # Protein reference file(s)
transcript_file = config["transcript_file"]  # Transcript reference file(s)
augustus = config["augustus"]  # Augustus configuration file


SAMPLE = list(glob_wildcards(join(input_dir, "{ids}.fasta")))[0]
if not SAMPLE:
    raise ValueError("No FASTA files found in the input directory!")

PROT = list(protein_file.split(","))  # Split multiple protein files
PROT_NAME = [os.path.basename(x) for x in PROT]
PROT_NAME = [x.rstrip(".fasta").rstrip(".faa") for x in PROT_NAME]

print(SAMPLE)
print(PROT_NAME)

rule All:
    input:
        # Repeats
        expand(join(input_dir,"{samples}-families.fa"),samples=SAMPLE),
        expand(join(result_dir,"{samples}-families.fa"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.fasta"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.fasta.masked"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.fasta.out.gff"),samples=SAMPLE),

        # Maker ctrl files
        expand(join(result_dir,"{samples}.maker_opts_rnd1.ctl"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.maker_opts_rnd2.ctl"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.maker_opts_rnd3.ctl"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.rnd1.maker.output/{samples}.rnd1_master_datastore_index.log"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.rnd2.maker.output/{samples}.rnd2_master_datastore_index.log"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.rnd3.maker.output/{samples}.rnd3_master_datastore_index.log"),samples=SAMPLE),

        # Maker GFFs
        expand(join(result_dir,"{samples}.rnd1.maker.output/{samples}.rnd1.all.gff"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.rnd1.maker.output/snap/rnd1.snap.hmm"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.rnd2.maker.output/{samples}.rnd2.all.gff"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.rnd3.maker.output/{samples}.rnd3.all.gff"),samples=SAMPLE),

        # rename gffs
        expand(join(result_dir,"{samples}.rnd2.maker.output/{samples}.gff3"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.rnd2.maker.output/{samples}.gtf"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.rnd2.maker.output/{samples}.clean.gtf"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.rnd2.maker.output/{samples}.fna"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.rnd2.maker.output/{samples}.faa"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.rnd3.maker.output/{samples}.gff3"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.rnd3.maker.output/{samples}.gtf"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.rnd3.maker.output/{samples}.clean.gtf"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.rnd3.maker.output/{samples}.fna"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.rnd3.maker.output/{samples}.faa"),samples=SAMPLE),

        # functional annotation
        expand(join(result_dir, "{samples}.rnd2.maker.output/{samples}.{prots}.gff3"),samples=SAMPLE,prots=PROT_NAME),
        expand(join(result_dir,"{samples}.rnd2.maker.output/{samples}.{prots}.blast"),samples=SAMPLE,prots=PROT_NAME),
        expand(join(result_dir, "{samples}.rnd3.maker.output/{samples}.{prots}.gff3"),samples=SAMPLE,prots=PROT_NAME),
        expand(join(result_dir,"{samples}.rnd3.maker.output/{samples}.{prots}.blast"),samples=SAMPLE,prots=PROT_NAME),

rule RepeatModeler:
  input:
    fa=join(input_dir,"{samples}.fasta"),
  output:
    fa=join(input_dir,"{samples}-families.fa"),
    fa2=join(result_dir,"{samples}.fasta"),
    rep=join(result_dir,"{samples}-families.fa"),
  params:
    rname="RepeatModeler",
    dir=input_dir,
    id="{samples}",
    opdir=result_dir,
  threads:
    48
  shell:
    """
    mkdir -p {params.opdir}
    cd {params.dir}
    module load repeatmodeler
    BuildDatabase -name {params.id} {input.fa}
    RepeatModeler -database {params.id} -pa {threads} -LTRStruct >& {params.id}.out
    cp {output.fa} {output.rep}
    cp {input.fa} {output.fa2}
    """

rule RepeatMasker:
  input:
    fa=join(result_dir,"{samples}.fasta"),
    rep=join(result_dir,"{samples}-families.fa"),
  output:
    fa=join(result_dir,"{samples}.fasta.masked"),
    gff=join(result_dir,"{samples}.fasta.out.gff"),
  params:
    rname="RepeatMasker",
    dir=result_dir,
    #rep=repeat_file,
    threads="48",
  shell:
    """
    cd {params.dir}
    module load repeatmasker
    RepeatMasker -u -s -poly -engine rmblast -pa {params.threads} -gff -no_is -gccalc -norna -lib {input.rep} {input.fa}
    """

rule maker_opts1:
    input:
        fa=join(result_dir,"{samples}.fasta.masked"),
    output:
        ctl=join(result_dir,"{samples}.maker_opts_rnd1.ctl"),
    params:
        rname="maker_opts1",
        protein=protein_file,
        transcript=transcript_file,
        augustus=augustus,
        scripts_path=join(result_dir,"param_files"),
        outdir=result_dir,
    shell:
        """
        mkdir -p {params.outdir}
        python3 {params.scripts_path}/generate_opts1.py {output.ctl} {input.fa} {params.protein} {params.transcript} {params.augustus}
        """



rule maker_rnd1:
    input:
        fa=join(result_dir,"{samples}.fasta.masked"),
        ctl=join(result_dir,"{samples}.maker_opts_rnd1.ctl"),
    output:
        log=join(result_dir,"{samples}.rnd1.maker.output/{samples}.rnd1_master_datastore_index.log"),
    params:
        rname="maker_rnd1",
        outid="{samples}.rnd1",
        outdir=result_dir,
        bopts=join(result_dir,"param_files/maker_bopts.log"),
        exe=join(result_dir,"param_files/maker_exe.log"),
    shell:
        """
        module load maker
        cd {params.outdir}
        valgrind mpiexec -c 1 -np 40 maker -RM_off -base {params.outid} {input.ctl} {params.bopts} {params.exe}
        """

rule make_gff1:
    input:
        log=join(result_dir,"{samples}.rnd1.maker.output/{samples}.rnd1_master_datastore_index.log"),
    output:
        gff=join(result_dir,"{samples}.rnd1.maker.output/{samples}.rnd1.all.gff"),
        snap=join(result_dir,"{samples}.rnd1.maker.output/snap/rnd1.snap.hmm"),
    params:
        rname="make_gff1",
        outdir=result_dir,
    shell:
        """
        module load maker snap
        cd {params.outdir}/{wildcards.samples}.rnd1.maker.output/
        gff3_merge -d {input.log}
        mkdir -p {params.outdir}/{wildcards.samples}.rnd1.maker.output/snap
        cd {params.outdir}/{wildcards.samples}rnd1.maker.output/snap
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
        fa=join(result_dir,"{samples}.fasta.masked"),
        gff=join(result_dir,"{samples}.rnd1.maker.output/{samples}.rnd1.all.gff"),
        snap=join(result_dir,"{samples}.rnd1.maker.output/snap/rnd1.snap.hmm"),
    output:
        ctl2=join(result_dir,"{samples}.maker_opts_rnd2.ctl"),
        ctl3=join(result_dir,"{samples}.maker_opts_rnd3.ctl"),
    params:
        rname="maker_opts2",
        protein=protein_file,
        transcript=transcript_file,
        scripts_path=join(result_dir,"param_files"),
        outdir=result_dir,
    shell:
        """
        mkdir -p {params.outdir}
        python3 {params.scripts_path}/generate_opts2.py {output.ctl2} {input.fa} {input.gff} {input.snap} {params.protein} {params.transcript}
        python3 {params.scripts_path}/generate_opts3.py {output.ctl3} {input.fa} {input.gff} {input.snap} {params.protein} {params.transcript}
        """

rule maker_rnd2:
    input:
        fa=join(result_dir,"{samples}.fasta.masked"),
        ctl=join(result_dir,"{samples}.maker_opts_rnd2.ctl"),
    output:
        log=join(result_dir,"{samples}.rnd2.maker.output/{samples}.rnd2_master_datastore_index.log"),
    params:
        rname="maker_rnd2",
        outid="{samples}.rnd2",
        outdir=result_dir,
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
        log=join(result_dir,"{samples}.rnd2.maker.output/{samples}.rnd2_master_datastore_index.log"),
        fa=join(result_dir,"{samples}.fasta"),
    output:
        gff=join(result_dir,"{samples}.rnd2.maker.output/{samples}.rnd2.all.gff"),
        faa=join(result_dir,"{samples}.rnd2.maker.output/{samples}.rnd2.all.faa"),
        fna=join(result_dir,"{samples}.rnd2.maker.output/{samples}.rnd2.all.fna"),
    params:
        rname="make_gff2",
        outdir=result_dir
    shell:
        """
        module load maker snap
        cd {params.outdir}/{wildcards.samples}.rnd2.maker.output/
        gff3_merge -d {input.log}
        /data/OpenOmics/references/brakerMake/gffread/gffread -g {input.fa} -y {output.faa} -x {output.fna} {output.gff}
        """

rule rnd2_gff_rename:
    input:
        gff=join(result_dir,"{samples}.rnd2.maker.output/{samples}.rnd2.all.gff"),
        aa=join(result_dir,"{samples}.rnd2.maker.output/{samples}.rnd2.all.faa"),
        cds=join(result_dir,"{samples}.rnd2.maker.output/{samples}.rnd2.all.fna"),
    output:
        map=temp(join(result_dir, "{samples}.rnd2.maker.output/{samples}.map")),
        gff=join(result_dir, "{samples}.rnd2.maker.output/{samples}.gff3"),
        aa=join(result_dir, "{samples}.rnd2.maker.output/{samples}.faa"),
        cds=join(result_dir, "{samples}.rnd2.maker.output/{samples}.fna"),
    params:
        species_id="{samples}",
        rname="rnd2_gff_rename",
    shell:
        """
        module load maker
        maker_map_ids --prefix {params.species_id} --justify 5  {input.gff} > {output.map}
        cp {input.gff} {output.gff}
        map_gff_ids {output.map} {output.gff}
        cp {input.aa} {output.aa}
        cp {input.cds} {output.cds}
        map_fasta_ids {output.map} {output.aa}
        map_fasta_ids {output.map} {output.cds}
        """

rule rnd2_gff2gtf:
    input:
        gff=join(result_dir, "{samples}.rnd2.maker.output/{samples}.gff3"),
    output:
        gtf=join(result_dir, "{samples}.rnd2.maker.output/{samples}.gtf"),
        clean=join(result_dir, "{samples}.rnd2.maker.output/{samples}.clean.gtf"),
    params:
        rname="rnd2_gff2gtf",
    shell:
        """
        module load agat/1.2.0 python
        agat_convert_sp_gff2gtf.pl --gff {input.gff} -o {output.gtf}
        python /data/OpenOmics/references/brakerMake/clean_gtf.py {output.gtf} > {output.clean}
        """

rule rnd2_gff_annot:
    input:
        gff=join(result_dir, "{samples}.rnd2.maker.output/{samples}.gff3"),
        prot=join(result_dir, "{samples}.rnd2.maker.output/{samples}.faa"),
        cds=join(result_dir, "{samples}.rnd2.maker.output/{samples}.fna"),
    output:
        gff=join(result_dir, "{samples}.rnd2.maker.output/{samples}.{prots}.gff3"),
        blast=join(result_dir,"{samples}.rnd2.maker.output/{samples}.{prots}.blast"),
    params:
        rname="rnd2_gff_annot",
        threads=8,
        uniprot=protein_file,
    shell:
        """
        module load blast maker
        blastp -query {input.prot} -db {params.uniprot} -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out {output.blast} -num_threads {params.threads}
        maker_functional_gff {params.uniprot} {output.blast} {input.gff} > {output.gff}
        """


rule maker_rnd3:
    input:
        fa=join(result_dir,"{samples}.fasta.masked"),
        ctl=join(result_dir,"{samples}.maker_opts_rnd3.ctl"),
    output:
        log=join(result_dir,"{samples}.rnd3.maker.output/{samples}.rnd3_master_datastore_index.log"),
    params:
        rname="maker_rnd3",
        outid="{samples}.rnd3",
        outdir=result_dir,
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
        log=join(result_dir,"{samples}.rnd3.maker.output/{samples}.rnd3_master_datastore_index.log"),
        fa=join(result_dir,"{samples}.fasta"),
    output:
        gff=join(result_dir,"{samples}.rnd3.maker.output/{samples}.rnd3.all.gff"),
        faa=join(result_dir,"{samples}.rnd3.maker.output/{samples}.rnd3.all.faa"),
        fna=join(result_dir,"{samples}.rnd3.maker.output/{samples}.rnd3.all.fna"),
    params:
        rname="make_gff3",
        outdir=result_dir,
    shell:
        """
        module load maker snap
        cd {params.outdir}/{wildcards.samples}.rnd3.maker.output/
        gff3_merge -d {input.log}
        /data/OpenOmics/references/brakerMake/gffread/gffread -g {input.fa} -y {output.faa} -x {output.fna} {output.gff}
        """

rule rnd3_gff_rename:
    input:
        gff=join(result_dir,"{samples}.rnd3.maker.output/{samples}.rnd3.all.gff"),
        aa=join(result_dir,"{samples}.rnd3.maker.output/{samples}.rnd3.all.faa"),
        cds=join(result_dir,"{samples}.rnd3.maker.output/{samples}.rnd3.all.fna"),
    output:
        map=temp(join(result_dir, "{samples}.rnd3.maker.output/{samples}.map")),
        gff=join(result_dir, "{samples}.rnd3.maker.output/{samples}.gff3"),
        aa=join(result_dir, "{samples}.rnd3.maker.output/{samples}.faa"),
        cds=join(result_dir, "{samples}.rnd3.maker.output/{samples}.fna"),
    params:
        species_id="{samples}",
        rname="rnd3_gff_rename",
    shell:
        """
        module load maker
        maker_map_ids --prefix {params.species_id} --justify 5  {input.gff} > {output.map}
        cp {input.gff} {output.gff}
        map_gff_ids {output.map} {output.gff}
        cp {input.aa} {output.aa}
        cp {input.cds} {output.cds}
        map_fasta_ids {output.map} {output.aa}
        map_fasta_ids {output.map} {output.cds}
        """

rule rnd3_gff2gtf:
    input:
        gff=join(result_dir, "{samples}.rnd3.maker.output/{samples}.gff3"),
    output:
        gtf=join(result_dir, "{samples}.rnd3.maker.output/{samples}.gtf"),
        clean=join(result_dir, "{samples}.rnd3.maker.output/{samples}.clean.gtf"),
    params:
        rname="rnd3_gff2gtf",
    shell:
        """
        module load agat/1.2.0 python
        agat_convert_sp_gff2gtf.pl --gff {input.gff} -o {output.gtf}
        python /data/OpenOmics/references/brakerMake/clean_gtf.py {output.gtf} > {output.clean}
        """

rule rnd3_gff_annot:
    input:
        gff=join(result_dir, "{samples}.rnd3.maker.output/{samples}.gff3"),
        prot=join(result_dir, "{samples}.rnd3.maker.output/{samples}.faa"),
        cds=join(result_dir, "{samples}.rnd3.maker.output/{samples}.fna"),
    output:
        gff=join(result_dir, "{samples}.rnd3.maker.output/{samples}.{prots}.gff3"),
        blast=join(result_dir,"{samples}.rnd3.maker.output/{samples}.{prots}.blast"),
    params:
        rname="rnd3_gff_annot",
        threads=8,
        uniprot=protein_file,
    shell:
        """
        module load blast maker
        blastp -query {input.prot} -db {params.uniprot} -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out {output.blast} -num_threads {params.threads}
        maker_functional_gff {params.uniprot} {output.blast} {input.gff} > {output.gff}
        """
