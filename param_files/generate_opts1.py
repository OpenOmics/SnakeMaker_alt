import sys

out1 = open(sys.argv[1],"w")
genome = sys.argv[2]
prot_fa = sys.argv[3]
tran_fa = sys.argv[4]
augustus = sys.argv[5]

out1.write("""#-----Genome (these are always required)\n""")
out1.write("""genome=""" + genome + """ #genome sequence (fasta file or fasta embeded in GFF3 file)\n""")
out1.write("""organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic\n""")
out1.write("""\n""")
out1.write("""#-----Re-annotation Using MAKER Derived GFF3\n""")
out1.write("""maker_gff= #MAKER derived GFF3 file\n""")
out1.write("""est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no\n""")
out1.write("""altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no\n""")
out1.write("""protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no\n""")
out1.write("""rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no\n""")
out1.write("""model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no\n""")
out1.write("""pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no\n""")
out1.write("""other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no\n""")
out1.write("""\n""")
out1.write("""#-----EST Evidence (for best results provide a file for at least one)\n""")
out1.write("""est= #set of ESTs or assembled mRNA-seq in fasta format\n""")
out1.write("""altest=""" + tran_fa + """ #EST/cDNA sequence file in fasta format from an alternate organism\n""")
out1.write("""est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file\n""")
out1.write("""altest_gff= #aligned ESTs from a closly relate species in GFF3 format\n""")
out1.write("""\n""")
out1.write("""#-----Protein Homology Evidence (for best results provide a file for at least one)\n""")
out1.write("""protein=""" + prot_fa + """ #protein sequence file in fasta format (i.e. from mutiple organisms)\n""")
out1.write("""protein_gff=  #aligned protein homology evidence from an external GFF3 file\n""")
out1.write("""\n""")
out1.write("""#-----Repeat Masking (leave values blank to skip repeat masking)\n""")
out1.write("""model_org= #select a model organism for RepBase masking in RepeatMasker\n""")
out1.write("""rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker\n""")
out1.write("""repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner\n""")
out1.write("""rm_gff= #pre-identified repeat elements from an external GFF3 file\n""")
out1.write("""prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no\n""")
out1.write("""softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)\n""")
out1.write("""\n""")
out1.write("""#-----Gene Prediction\n""")
out1.write("""snaphmm= #SNAP HMM file\n""")
out1.write("""gmhmm= #GeneMark HMM file\n""")
out1.write("""augustus_species=""" +augustus +""" #Augustus gene prediction species model\n""")
out1.write("""fgenesh_par_file= #FGENESH parameter file\n""")
out1.write("""pred_gff= #ab-initio predictions from an external GFF3 file\n""")
out1.write("""model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)\n""")
out1.write("""run_evm=0 #run EvidenceModeler, 1 = yes, 0 = no\n""")
out1.write("""est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no\n""")
out1.write("""protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no\n""")
out1.write("""trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no\n""")
out1.write("""snoscan_rrna= #rRNA file to have Snoscan find snoRNAs\n""")
out1.write("""snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs\n""")
out1.write("""unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no\n""")
out1.write("""allow_overlap=0 #allowed gene overlap fraction (value from 0 to 1, blank for default)\n""")
out1.write("""\n""")
out1.write("""#-----Other Annotation Feature Types (features MAKER doesn't recognize)\n""")
out1.write("""other_gff= #extra features to pass-through to final MAKER generated GFF3 file\n""")
out1.write("""\n""")
out1.write("""#-----External Application Behavior Options\n""")
out1.write("""alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases\n""")
out1.write("""cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)\n""")
out1.write("""\n""")
out1.write("""#-----MAKER Behavior Options\n""")
out1.write("""max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)\n""")
out1.write("""min_contig=5000 #skip genome contigs below this length (under 10kb are often useless)\n""")
out1.write("""\n""")
out1.write("""pred_flank=200 #flank for extending evidence clusters sent to gene predictors\n""")
out1.write("""pred_stats=0 #report AED and QI statistics for all predictions as well as models\n""")
out1.write("""AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)\n""")
out1.write("""min_protein=0 #require at least this many amino acids in predicted proteins\n""")
out1.write("""alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no\n""")
out1.write("""always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no\n""")
out1.write("""map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no\n""")
out1.write("""keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)\n""")
out1.write("""\n""")
out1.write("""split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)\n""")
out1.write("""min_intron=20 #minimum intron length (used for alignment polishing)\n""")
out1.write("""single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no\n""")
out1.write("""single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'\n""")
out1.write("""correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes\n""")
out1.write("""\n""")
out1.write("""tries=2 #number of times to try a contig if there is a failure for some reason\n""")
out1.write("""clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no\n""")
out1.write("""clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no\n""")
out1.write("""TMP= #specify a directory other than the system default temporary directory for temporary files\n""")
