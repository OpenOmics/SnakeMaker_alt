
###########################################################################
# Genome annotation pipeline using Funannotate in Snakemake
# Snakemake/5.13.0
###########################################################################
from os.path import join
from snakemake.io import expand, glob_wildcards
from Bio import SeqIO

result_dir = config["result_dir"]
Lineage_name = config["lineage_name"]
busco_db = Lineage_name.split("_odb")[0]
augustus_name = config["augustus_name"]
input_dir = config["input_dir"]

#repeat_file = config["repeat_file"]
protein_file = config["protein_file"]
transcript_file = config["transcript_file"]
funannotate_dir = config["funannotate_dir"]
species = config["species"]
species_id = species.replace(" ","_")
augustus="caenorhabditis",

#SAMPLE = ["pilon124_round3_chromosomesnumbered"]
SAMPLE = list(glob_wildcards(join(input_dir, "{ids}.fasta")))[0]
#ID = glob_wildcards(join(fastq_dir, "{ids}.fastq.gz"))

print(SAMPLE)

rule All:
    input:
        # Repeats
        expand(join(input_dir,"{samples}-families.fa"),samples=SAMPLE),
        expand(join(result_dir,"funannotate/{samples}/{samples}-families.fa"),samples=SAMPLE),
        expand(join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta.masked"),samples=SAMPLE),
        expand(join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta.out.gff"),samples=SAMPLE),

        # Funannotate - preprocessing
        join(funannotate_dir,"databases",busco_db + ".tar.gz"),
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

        # Maker ctrl files
        expand(join(result_dir,"maker/{samples}/maker_opts_rnd1.ctl"),samples=SAMPLE),
        expand(join(result_dir,"maker/{samples}/maker_opts_rnd2.ctl"),samples=SAMPLE),
        expand(join(result_dir,"maker/{samples}/rnd1.maker.output/rnd1_master_datastore_index.log"),samples=SAMPLE),
        expand(join(result_dir,"maker/{samples}/rnd2.maker.output/rnd2_master_datastore_index.log"),samples=SAMPLE),

        # Maker GFFs
        expand(join(result_dir,"maker/{samples}/rnd1.maker.output/rnd1.all.gff"),samples=SAMPLE),
        expand(join(result_dir,"maker/{samples}/rnd1.maker.output/snap/rnd1.snap.hmm"),samples=SAMPLE),
        expand(join(result_dir,"maker/{samples}/rnd2.maker.output/rnd2.all.gff"),samples=SAMPLE),

        # Merged GFFs
        expand(join(result_dir,"{samples}.combined.gtf"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.transdecoder.gtf"),samples=SAMPLE),
        expand(join(result_dir,"{samples}.final.gtf"),samples=SAMPLE),

rule fun_setup:
    input:
        fa=transcript_file,
    output:
        dbout=join(funannotate_dir,"databases",busco_db + ".tar.gz"),
    params:
        rname="fun_setup",
        fun_dir=join(funannotate_dir,"databases"),
        busco_db=busco_db,
    shell:
        """
        module load funannotate/1.8.9
        funannotate setup -b {params.busco_db} -d {params.fun_dir}
        """
        
rule RepeatModeler:
  input:
    fa=join(input_dir,"{samples}.fasta"),
  output:
    fa=join(input_dir,"{samples}-families.fa"),
    rep=join(result_dir,"funannotate/{samples}/{samples}-families.fa"),
  params:
    rname="RepeatModeler",
    dir=input_dir,
    id="{samples}",
    opdir=join(result_dir,"funannotate/{samples}"),
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
    """

rule RepeatMasker:
  input:
    fa=join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta"),
    rep=join(result_dir,"funannotate/{samples}/{samples}-families.fa"),
  output:
    fa=join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta.masked"),
    gff=join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta.out.gff"),
  params:
    rname="RepeatMasker",
    dir=join(result_dir,"funannotate/{samples}"),
    #rep=repeat_file,
    threads="48",
  shell:
    """
    cd {params.dir}
    module load repeatmasker
    RepeatMasker -u -s -poly -engine rmblast -pa {params.threads} -gff -no_is -gccalc -norna -lib {input.rep} {input.fa}
    """

rule fun_clean:
    input:
        fa=join(input_dir,"{samples}.fasta"),
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
                out.write(">" + r.description + "\n" + str(r.seq) + "\n")

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
        rep=join(result_dir,"funannotate/{samples}/{samples}-families.fa"),
    output:
        fa=join(result_dir,"funannotate/{samples}/{samples}.CSM.fasta"),
    params:
        rname="fun_mask",
        #rep=repeat_file,
        dir=join(result_dir,"funannotate/{samples}"),
        threads="32",
    shell:
        """
        module load funannotate
        mkdir -p {params.dir}
        funannotate mask -i {input.fa} --cpus {params.threads} -o {output.fa} -l {input.rep}
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
        fun_dir=join(funannotate_dir,"databases"),
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
        fun_dir=join(funannotate_dir,"databases"),
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

rule maker_opts1:
    input:
        fa=join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta"),
    output:
        ctl=join(result_dir,"maker/{samples}/maker_opts_rnd1.ctl"),
    params:
        rname="maker_opts1",
        protein=protein_file,
        transcript=transcript_file,
        augustus=augustus,
        scripts_path=join(result_dir,"param_files"),
        outdir=join(result_dir,"maker/{samples}"),
    shell:
        """
        mkdir -p {params.outdir}
        python3 {params.scripts_path}/generate_opts1.py {output.ctl} {input.fa} {params.protein} {params.transcript} {params.augustus}
        """



rule maker_rnd1:
    input:
        fa=join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta"),
        ctl=join(result_dir,"maker/{samples}/maker_opts_rnd1.ctl"),
    output:
        log=join(result_dir,"maker/{samples}/rnd1.maker.output/rnd1_master_datastore_index.log"),
    params:
        rname="maker_rnd1",
        outid="rnd1",
        outdir=join(result_dir,"maker/{samples}"),
        bopts=join(result_dir,"param_files/maker_bopts.log"),
        exe=join(result_dir,"param_files/maker_exe.log"),
    shell:
        """
        module load maker
        cd {params.outdir}
        mpiexec -np 40 maker -RM_off -base {params.outid} {input.ctl} {params.bopts} {params.exe}
        """

rule make_gff1:
    input:
        log=join(result_dir,"maker/{samples}/rnd1.maker.output/rnd1_master_datastore_index.log"),
    output:
        gff=join(result_dir,"maker/{samples}/rnd1.maker.output/rnd1.all.gff"),
        snap=join(result_dir,"maker/{samples}/rnd1.maker.output/snap/rnd1.snap.hmm"),
    params:
        rname="make_gff1",
        outdir=join(result_dir,"maker/{samples}"),
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
        fa=join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta"),
        gff=join(result_dir,"maker/{samples}/rnd1.maker.output/rnd1.all.gff"),
        snap=join(result_dir,"maker/{samples}/rnd1.maker.output/snap/rnd1.snap.hmm"),
    output:
        ctl=join(result_dir,"maker/{samples}/maker_opts_rnd2.ctl"),
    params:
        rname="maker_opts2",
        protein=protein_file,
        transcript=transcript_file,
        scripts_path=join(result_dir,"param_files"),
        outdir=join(result_dir,"maker/{samples}"),
    shell:
        """
        mkdir -p {params.outdir}
        python3 {params.scripts_path}/generate_opts2.py {output.ctl} {input.fa} {input.gff} {input.snap} {params.protein} {params.transcript}
        """

rule maker_rnd2:
    input:
        fa=join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta"),
        ctl=join(result_dir,"maker/{samples}/maker_opts_rnd2.ctl"),
    output:
        log=join(result_dir,"maker/{samples}/rnd2.maker.output/rnd2_master_datastore_index.log"),
    params:
        rname="maker_rnd2",
        outid="rnd2",
        outdir=join(result_dir,"maker/{samples}"),
        bopts=join(result_dir,"param_files/maker_bopts.log"),
        exe=join(result_dir,"param_files/maker_exe.log"),
    shell:
        """
        module load maker
        cd {params.outdir}
        mpiexec -np 40 maker -RM_off -base {params.outid} {input.ctl} {params.bopts} {params.exe}
        """

rule make_gff2:
    input:
        log=join(result_dir,"maker/{samples}/rnd2.maker.output/rnd2_master_datastore_index.log"),
    output:
        gff=join(result_dir,"maker/{samples}/rnd2.maker.output/rnd2.all.gff"),
    params:
        rname="make_gff2",
        outdir=join(result_dir,"maker/{samples}"),
    shell:
        """
        module load maker snap
        cd {params.outdir}/rnd2.maker.output/
        gff3_merge -d {input.log}
        """

rule gffmerge:
    input:
        maker=join(result_dir,"maker/{samples}/rnd2.maker.output/rnd2.all.gff"),
        fun=join(result_dir,"funannotate/{samples}/fun/annotate_results/"+species_id+".gff3"),
    output:
        merge=join(result_dir,"{samples}.combined.gtf"),
    params:
        rname="gffmerge",
        id="{samples}",
    shell:
        """
        module load gffcompare bedtools
        gffcompare -o {params.id} {input.maker} {input.fun}
        """

rule gff_annot:
    input:
        fa=join(result_dir,"funannotate/{samples}/{samples}.cleaned.sorted.fasta"),
        merge=join(result_dir,"{samples}.combined.gtf"),
    output:
        cdna=temp(join(result_dir,"{samples}.cdna.fa")),
        merge=temp(join(result_dir,"{samples}.combined.gff3")),
        gff1=join(result_dir,"{samples}.cdna.fa.transdecoder_dir/longest_orfs.cds.best_candidates.gff3"),
        gff2=temp(join(result_dir,"{samples}.transdecoder.gff3")),
    params:
        rname="gff_annot",
        id="{samples}",
    shell:
        """
        module load TransDecoder
        gtf_genome_to_cdna_fasta.pl {input.merge} {input.fa} > {output.cdna}
        gtf_to_alignment_gff3.pl {input.merge} > {output.merge}
        TransDecoder.LongOrfs -t {output.cdna}
        TransDecoder.Predict -t {output.cdna}
        cdna_alignment_orf_to_genome_orf.pl \
        {output.gff1} \
        {output.merge} \
        {output.cdna} > {output.gff2}
        """

rule gff_format:
    input:
        gff2=join(result_dir,"{samples}.transdecoder.gff3"),
    output:
        gff3=join(result_dir,"{samples}.transdecoder.gtf"),
        gff4=join(result_dir,"{samples}.final.gtf"),
    params:
        rname="gff_format",
        param_dir=join(result_dir,"param_files"),
    shell:
        """
        module load singularity python
        singularity exec -B $PWD {params.param_dir}/agat_0.8.0--pl5262hdfd78af_0.sif agat_convert_sp_gff2gtf.pl --gff {input.gff2} -o {output.gff3}
        python3 {params.param_dir}/clean_gtf.py {output.gff3} > {output.gff4}
        """
