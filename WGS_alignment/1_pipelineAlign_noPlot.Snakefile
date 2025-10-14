report: "report/workflow.rst"

my_basedir = workflow.current_basedir

import glob
import re
import os
import shutil

if not os.path.isdir(config['WORKPATH']): 
	os.mkdir(config['WORKPATH']) 

#Set working directory

# If in the config PATH_BAM is empty (start from FASTQ)
if not config.get("PATH_BAM"):
	config['workdir'] = config['WORKPATH'] + "BAM" + "/"
	config['PATH_BAM'] = config['workdir']
else:
	config['workdir'] = config['PATH_BAM']

if not os.path.isdir(config['workdir']): 
	os.mkdir(config['workdir']) 

# If in the config SENTIEON_MODULE is empty, then add Sentieon to the PATH. Otherwise, load the specified module
if not config.get("SENTIEON_MODULE"):
    os.environ["PATH"] = f"{config['SENTIEON_PATH']}/bin:" + os.environ["PATH"]
    print(f"[INFO] Using Sentieon from PATH: {config['SENTIEON_PATH']}")
else:
    shell.prefix(f"module load {config['SENTIEON_MODULE']}; ")
    print(f"[INFO] Using Sentieon module: {config['SENTIEON_MODULE']}")


os.environ["SENTIEON_LICENSE"]= config['SENTIEON_LICENSE']
os.environ["BCFTOOLS_PLUGINS"]= config['BCF_DIR']
   
config['KNOWN_SITES_STR'] = " "
if config['KNOWN_SITES'] is not None:
	for key in config['KNOWN_SITES']:
		config['KNOWN_SITES_STR'] = config['KNOWN_SITES_STR'] +"-k " + config['KNOWN_SITES'][key] +" "

config['KNOWN_SITES_STR']

config['my_basedir'] = my_basedir


os.chdir(config['workdir'])

def clearSamp(sample):
	if not config.get("PATH_BAM"):
		sample = sample.replace(config['PATH_FASTQ'], '')
		if config['PAIRED']==True:
			sample = sample.replace("_"+config['FASTQ_SUFFIXES']["1_SUFFIX"], '')
		else:
			sample = sample.replace(config['FASTQ_SUFFIX'], '')
	else:
		sample = sample.replace(config['PATH_BAM'], '')
		sample = sample.replace(config['BAM_SUFFIX'], '')

	return sample


def get_sampleNames(wildcards):

	import glob

	if not config.get("PATH_BAM"):
		if config['PAIRED']==True:
			sampleNames = glob.glob(config['PATH_FASTQ']+"*"+config['FASTQ_SUFFIXES']["1_SUFFIX"])
		else:
			sampleNames = glob.glob(config['PATH_FASTQ']+"*"+config['FASTQ_SUFFIX'])
	else:
		sampleNames = glob.glob(config['PATH_BAM']+"*"+config['BAM_SUFFIX'])

	sampleNames = list(map(lambda sample: clearSamp(sample), sampleNames))

	print(sampleNames)

	return sampleNames


def get_fastqFile(wildcards):
	ls = []
	if config['PAIRED']==True:
		ls.append(config['PATH_FASTQ'] + wildcards.sample + "_"+ config['FASTQ_SUFFIXES']["1_SUFFIX"]) 
		ls.append(config['PATH_FASTQ'] + wildcards.sample + "_"+ config['FASTQ_SUFFIXES']["2_SUFFIX"])
	else:
		ls.append(config['PATH_FASTQ'] + wildcards.sample + config['FASTQ_SUFFIX']) 

	return ls

def get_fileOutputName(namefile, startPATH = ""):
	ls = []
	for sample in get_sampleNames(""):
		ls.append(startPATH + sample + namefile) 

	return ls

import re

def get_fileOutputName2(namefile, startPATH = ""):
	ls = []
	suffix_clear = re.sub(".bam", "", config['BAM_SUFFIX'])
	for sample in get_sampleNames(""):
	  ls.append(startPATH + sample + suffix_clear + "/" + sample + namefile)
		#ls.append(startPATH + sample + "_recal" + "/" + sample + namefile) 

	return ls


def get_baseRecalibration(wildcards):
	ls = get_fileOutputName("_recal.bam") 
	return ls

def get_variantCallerDNAscope(wildcards):
	ls = get_fileOutputName("-filtDNAscope.vcf.gz")
	return ls

def get_variantCallerTNscope(wildcards):
	ls = get_fileOutputName("_output_tnscope.filter.vcf.gz")
	return ls

def get_fastqc(wildcards):
	ls = get_fileOutputName("_fastqc.zip", "qc/fastqc/")
	return ls

def get_samtools_flagstat(wildcards):
	ls = get_fileOutputName(".bam.flagstat", "qc/samtools/")
	return ls

def get_qc_corr(wildcards):
  suffix_clear = re.sub(".bam", "", config['BAM_SUFFIX'])
	#ls = get_fileOutputName2("_recal.W_gc_outliers_removed_smoothed.heatmap.png", config['WORKPATH']+"GC_correction_output/")
  ls = get_fileOutputName2(suffix_clear+".W_gc_outliers_removed_smoothed.heatmap.png", config['WORKPATH']+"GC_correction_output/")
  return ls

def get_Fragmentomics_file(wildcards):
	ls = get_fileOutputName("_motif_bin_"+ str(config['BIN_SIZE']) + "_DF.RData", config['WORKPATH']+"Fragmentomics_features/")
	#ls = [k for k in ls if "LB-" in k]
	#ls = [k for k in ls if "IWI" in k]
	return ls

print(get_Fragmentomics_file(""))


##### ALIGNMENT #######

if not config.get("PATH_BAM"):
	rule all:
		input:
			get_baseRecalibration
			#get_Fragmentomics_file

	# 1 Mapping reads with BWA-MEM, sorting for normal sample
	rule mapping:
		input:
			get_fastqFile
		output:
			"{sample}_sorted.bam"
		params:
			fastapath = config['FASTA'],
			platform = config['PLATFORM'],
			nthreads = config['NUM_THREADS']
		threads: config['NUM_THREADS']
		resources:
			mem_mb=102400 
		shell:
			"""
			(sentieon bwa mem -M -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:{params.platform}' -t {params.nthreads} -K 10000000 {params.fastapath} {input} || echo -n error ) | sentieon util sort -o {wildcards.sample}_sorted.bam -t {params.nthreads} --sam2bam -i -
			"""

	rule metrics:
		input:
			"{sample}_sorted.bam"
		output:
			report("{sample}_gc_metrics.txt", category="metrics")
		threads: 1
		resources:
			mem_mb=102400 
		params:
			fastapath = config['FASTA'],
			nthreads = config['NUM_THREADS']
		shell:
			"""
			sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}_sorted.bam --algo MeanQualityByCycle {wildcards.sample}_mq_metrics.txt --algo QualDistribution {wildcards.sample}_qd_metrics.txt --algo GCBias --summary {wildcards.sample}_gc_summary.txt {wildcards.sample}_gc_metrics.txt --algo AlignmentStat --adapter_seq '' {wildcards.sample}_aln_metrics.txt --algo InsertSizeMetricAlgo {wildcards.sample}_is_metrics.txt
			#sentieon plot GCBias -o {wildcards.sample}_gc-report.pdf {wildcards.sample}_gc_metrics.txt
			"""

	rule removeDuplicate:
		input:
			"{sample}_sorted.bam"
		output:
			"{sample}_deduped.bam"
		params:
			nthreads = config['NUM_THREADS']
		threads: config['NUM_THREADS']
		resources:
			mem_mb=102400 
		shell:
			"""
			sentieon driver -t {params.nthreads} -i {wildcards.sample}_sorted.bam --algo LocusCollector --fun score_info {wildcards.sample}_score.txt
			sentieon driver -t {params.nthreads} -i {wildcards.sample}_sorted.bam --algo Dedup --rmdup --score_info {wildcards.sample}_score.txt --metrics {wildcards.sample}_dedup_metrics.txt {wildcards.sample}_deduped.bam
			"""

	rule indelRealigner:
		input:
			"{sample}_deduped.bam"
		output:
			"{sample}_realigned.bam"
		threads: config['NUM_THREADS']
		resources:
			mem_mb=102400 
		params:
			fastapath = config['FASTA'],
			knowsitestr = config['KNOWN_SITES_STR'],
			nthreads = config['NUM_THREADS']
		shell:
			"""
			sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}_deduped.bam --algo Realigner{params.knowsitestr}{wildcards.sample}_realigned.bam
			"""

	rule baseRecalibration:
		input:
			"{sample}_realigned.bam"
		output:
			"{sample}_recal.bam"
		threads: config['NUM_THREADS']
		resources:
			mem_mb=102400 
		params:
			fastapath = config['FASTA'],
			knowsitestr = config['KNOWN_SITES_STR'],
			dbsnppath = config['DBSNP'],
			nthreads = config['NUM_THREADS']
		shell:
			"""
			sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}_realigned.bam --algo QualCal -k {params.dbsnppath}{params.knowsitestr}{wildcards.sample}_recal_data.table
			sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}_realigned.bam -q {wildcards.sample}_recal_data.table --algo QualCal -k {params.dbsnppath}{params.knowsitestr}{wildcards.sample}_recal_data.table.post
			sentieon driver -t {params.nthreads} --algo QualCal --before {wildcards.sample}_recal_data.table --after {wildcards.sample}_recal_data.table.post {wildcards.sample}_recal.csv
			#sentieon driver -t {params.nthreads} --algo QualCal --plot --before {wildcards.sample}_recal_data.table --after {wildcards.sample}_recal_data.table.post {wildcards.sample}_recal.csv
			#sentieon plot QualCal -o {wildcards.sample}_recal_plots.pdf {wildcards.sample}_recal.csv
			sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}_realigned.bam -q {wildcards.sample}_recal_data.table --algo ReadWriter {wildcards.sample}_recal.bam
			"""
else:
	print("Start From BAM")
	rule all:
		input:
			get_Fragmentomics_file

SUFFIX_CLEAR = re.sub(".bam", "", config['BAM_SUFFIX'])

##### GC CORRECTION #######

rule install_GCparagon:
	output:
		#str(config['my_basedir']) + "/GCparagon/GCparagon_installed.txt"
		str(config['my_basedir']) + "/GCparagon/setup.py"
	conda:
		"conda_env/GCparagon_py3.10_env.yml"  
	params:
		my_basedir = config['my_basedir']
	shell:
		"""
		cd {params.my_basedir}
		git clone https://github.com/BGSpiegl/GCparagon
		cd GCparagon
		#pip install .
		#echo "GCparagon_installed" > GCparagon_installed.txt
		cd {params.my_basedir}
		"""

rule GC_corr:
	input:
		str(config['my_basedir']) + "/GCparagon/setup.py",
		str(config['PATH_BAM']) + "{sample}"+ str(config['BAM_SUFFIX'])
	output:
		config['WORKPATH'] + "GC_correction_output/{sample}"+SUFFIX_CLEAR+"/{sample}"+SUFFIX_CLEAR +".W_gc_outliers_removed_smoothed.heatmap.png"
	threads: config['NUM_THREADS']
	resources:
		mem_mb=102400 
	params:
		PATH_BAM = config['PATH_BAM'],
		OUTPUT_DIR = config['WORKPATH'],
		SUFFIX_BAM = config['BAM_SUFFIX'],
		nthreads = config['NUM_THREADS'],
		#my_basedir = "/home3/adefalco/Fate-AI/"
		my_basedir = config['my_basedir']
	conda:
		"conda_env/GCparagon_py3.10_env.yml"  
	log:
		"{sample}.log"
	shell:
		"python {params.my_basedir}/GCparagon/src/GCparagon/correct_GC_bias.py --bam {params.PATH_BAM}/{wildcards.sample}{params.SUFFIX_BAM} -rtb {params.my_basedir}/acc_files/genome2bit/hg38.2bit -rgcd {params.my_basedir}/GCparagon/accessory_files/hg38_reference_GC_content_distribution.tsv -c {params.my_basedir}/GCparagon/accessory_files/hg38_minimalExclusionListOverlap_1Mbp_intervals_33pcOverlapLimited.FGCD.bed -o {params.OUTPUT_DIR}GC_correction_output/ --threads {params.nthreads}"

##### VARIANT CALLING DNA SCOPE #######

rule HC_VariantCaller:
	input:
		str(config['PATH_BAM']) + "{sample}" + str(config['BAM_SUFFIX'])
	output:
		"{sample}-output-hc.vcf.gz"
	threads: config['NUM_THREADS']
	resources:
		mem_mb=102400,
		time="12:00:00" 
	params:
		fastapath = config['FASTA'],
		dbsnppath = config['DBSNP'],
		nthreads = config['NUM_THREADS']
	log:
		"logs/Variant_calling_DNAscope/{sample}_HC.log"
	shell:
		"sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}_recal.bam --algo Haplotyper -d {params.dbsnppath} --emit_conf=30 --call_conf=30 {wildcards.sample}-output-hc.vcf.gz"


rule DNAscope_variantCaller:
	input:
		"{sample}-output-hc.vcf.gz"
	output:
		"{sample}-filtDNAscope.vcf.gz"
	threads: config['NUM_THREADS']
	resources:
		mem_mb=102400,
		time="12:00:00" 
	params:
		fastapath = config['FASTA'],
		dbsnppath = config['DBSNP'],
		bcfpath = config['BCF_DIR'],
		MLmodel = str(config['my_basedir'])+"/acc_files/"+str(config['ML_MODEL_N']),
		nthreads = config['NUM_THREADS']
	log:
		"logs/Variant_calling_DNAscope/{sample}_DNAscope.log"
	shell:
		"""
		sentieon driver -t {params.nthreads} -r {params.fastapath} -i {wildcards.sample}_recal.bam --algo DNAscope -d {params.dbsnppath} --model {params.MLmodel} {wildcards.sample}-tmpDNAscope.vcf.gz
		sentieon driver -t {params.nthreads} -r {params.fastapath} --algo DNAModelApply --model {params.MLmodel} -v {wildcards.sample}-tmpDNAscope.vcf.gz {wildcards.sample}-DNAscope.vcf.gz
		{params.bcfpath} filter -s ML_FAIL -i INFO/ML_PROB > 0.81 {wildcards.sample}-DNAscope.vcf.gz -O z -m x -o {wildcards.sample}-filtDNAscope.vcf.gz
		"""

##### VARIANT CALLING TN SCOPE #######

rule Variant_calling:
	input:
		str(config['PATH_BAM']) + "{sample}" + str(config['BAM_SUFFIX'])
	output:
		"{sample}_output_tnscope.pre_filter.vcf.gz"
	threads: config['NUM_THREADS']
	resources:
		mem_mb=102400,
		time="12:00:00" 
	params:
		fastapath = config['FASTA'],
		dbsnppath = config['DBSNP'],
		nthreads = config['NUM_THREADS'],
		inteval = str(config['my_basedir'])+"/acc_files/GenomeSizeHg38.txt"
	log: out = "logs/Variant_calling_TN/{sample}.log",
		 err = "logs/Variant_calling_TN/{sample}.err"
	shell:
		"""
		sentieon driver -r {params.fastapath} -t {params.nthreads} -i {wildcards.sample}_recal.bam \
   		--interval {params.inteval} \
		--algo TNscope \
    	--tumor_sample {wildcards.sample} \
    	--dbsnp {params.dbsnppath} \
    	--disable_detector sv \
    	--min_tumor_allele_frac 0.1 \
    	--min_tumor_lod 2.0 \
    	--min_init_tumor_lod 0.5 \
    	--assemble_mode 4 \
    	--pcr_indel_model NONE \
    	--min_base_qual 20 \
    	--resample_depth 100 \
    	{wildcards.sample}_output_tnscope.pre_filter.vcf.gz
		"""

rule Variant_filtration:
	input:
		"{sample}_output_tnscope.pre_filter.vcf.gz"
	output:
		"{sample}_output_tnscope.filter.vcf.gz"
	threads: config['NUM_THREADS']
	resources:
		mem_mb=102400,
		time="12:00:00" 
	params:
		tnscope_filter = str(config['my_basedir'])+"/acc_files/"+str(config['TNSCOPE_FILTER'])
	log:
		"logs/Variant_calling_TN/{sample}_filt.log"
	shell:
		"""
		sentieon pyexec {params.tnscope_filter} \
    	-v {wildcards.sample}_output_tnscope.pre_filter.vcf.gz \
    	--tumor_sample {wildcards.sample} \
    	-x ctdna --min_tumor_af 0.1 --min_depth 3 \
    	{wildcards.sample}_output_tnscope.filter.vcf.gz
		"""

rule Do_Variant_Calling:
	input:
		get_variantCallerDNAscope,
		get_variantCallerTNscope

##### QC PIPELINE #######

rule fastqc:
	input:
		get_fastqFile
	output:
		html="qc/fastqc/{sample}.html",
		zip="qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
	params: "--quiet"
	log:
		"logs/fastqc/{sample}.log"
	threads: 1
	resources:
		mem_mb=102400,
		time="12:00:00" 
	wrapper:
		"v1.23.5/bio/fastqc"

rule samtools_flagstat:
	input:
		config['PATH_BAM'] + "{sample}"+config['BAM_SUFFIX']
	output:
		"qc/samtools/{sample}.bam.flagstat",
	log:
		"{sample}.log",
	threads: 1
	resources:
		mem_mb=102400,
		time="12:00:00" 
	params:
		extra="",  # optional params string
	wrapper:
		"v1.23.5/bio/samtools/flagstat"

rule multiqc_all:
	input:
		get_fastqc,
		get_samtools_flagstat
	threads: 1
	resources:
		mem_mb=102400,
		time="12:00:00" 
	output:
		"qc/multiqc_all.html"
	params:
		extra="",  # Optional: extra parameters for multiqc.
		use_input_files_only=True, # Optional, use only a.txt and don't search folder samtools_stats for files
	log:
		"logs/multiqc.log"
	wrapper:
		"v1.23.5/bio/multiqc"


##### FRAGMENTOMICS #######

rule all_Fragmentomics:
	input:
		get_Fragmentomics_file

rule get_bin_bed_file:
	input:
		str(config['my_basedir'])+"/acc_files/GenomeSizeHg38.txt"
	output:
		str(config['my_basedir'])+"/acc_files/genome_{bin_size}.bed"
	threads: 1
	resources:
		mem_mb=102400,
		time="12:00:00" 
	params:
		bed_tools = config['BED_TOOLS_DIR'],
		BIN_SIZE = config['BIN_SIZE'],
		BIN_SIZE_1 = config['BIN_SIZE']+1,
		my_basedir = config['my_basedir']
	shell:
		"""
		cd {params.my_basedir}
		{params.bed_tools} makewindows -g acc_files/GenomeSizeHg38.txt -w {params.BIN_SIZE} -s {params.BIN_SIZE_1} > acc_files/genome_{params.BIN_SIZE}.bed
		"""


rule get_Fragmentomics_fragm_motif:
	input:
		str(config['my_basedir'])+"/acc_files/genome_{bin_size}.bed",
		"{sample}"+str(config['BAM_SUFFIX']),
		config['WORKPATH'] + "GC_correction_output/{sample}"+SUFFIX_CLEAR+"/{sample}"+SUFFIX_CLEAR +".W_gc_outliers_removed_smoothed.heatmap.png"
	output:
		config['WORKPATH'] + "Fragmentomics_output/{sample}_{bin_size}_res_frag_motif.RData"
	threads: config['NUM_THREADS']
	resources:
		mem_mb=102400,
		time="12:00:00" 
	params:
		my_basedir = config['my_basedir'],
		suffix_bam = config['BAM_SUFFIX'],
		path_bam = config['PATH_BAM'],
		output_dir = config['WORKPATH'],
		NUM_THREADS= config['NUM_THREADS'],
		MAPQ = config['MAPQ'],
		BIN_SIZE = config['BIN_SIZE'],
		fastapath = config['FASTA'],
		samtoolspath = config['SAMTOOLS_DIR'],
		MAX_FRAG_LENGHT = config['MAX_FRAG_LENGHT'],
		SUFFIX_BAM = config['BAM_SUFFIX']
	log:
		"logs/fastqc/{sample}_{bin_size}.log"
	conda:
		"conda_env/envR.yml"
	shell:
		"Rscript {params.my_basedir}/Script/Fragmentomics_Extraction.R {wildcards.sample} {params.path_bam}/{wildcards.sample}{params.suffix_bam} {params.output_dir}/Fragmentomics_output/ {params.BIN_SIZE} {params.NUM_THREADS} {params.fastapath} {params.samtoolspath} {params.MAPQ} {params.MAX_FRAG_LENGHT} {params.SUFFIX_BAM}"

print(config['WORKPATH'])

rule get_Fragmentomics_features:
	input:
		config['WORKPATH'] + "Fragmentomics_output/{sample}_{bin_size}_res_frag_motif.RData"
	output:
		config['WORKPATH'] + "Fragmentomics_features/{sample}_fragm_bin_{bin_size}_DF.RData",
		config['WORKPATH'] + "Fragmentomics_features/{sample}_motif_bin_{bin_size}_DF.RData"
	threads: config['NUM_THREADS']
	resources:
		mem_mb=102400,
		time="12:00:00" 
	params:
		my_basedir = config['my_basedir'],
		output_dir = config['WORKPATH'],
		NUM_THREADS= config['NUM_THREADS'],
		BIN_SIZE = config['BIN_SIZE']
	log:
		"logs/fastqc/{sample}_{bin_size}.log"
	conda:
		"conda_env/envR.yml"  
	shell:
		"Rscript {params.my_basedir}/Script/Fragmentomics_2D_Profile.R {wildcards.sample} {params.output_dir}/Fragmentomics_output/ {params.output_dir}/Fragmentomics_features/ {params.BIN_SIZE} {params.NUM_THREADS} {params.my_basedir}"


##### IchorCNA #######
"""
config["samples"] = get_sampleNames("")

rule all_ichorCNA:
	input:
		expand("output_ichorCNA/ichorCNA/{tumor}/{tumor}.cna.seg", tumor=config["samples"])
	output:
		report3=report("output_ichorCNA/results.tsv", category="TF"),
	conda:
		"envs/envIchorCNA.yaml"		
	log:
		"logs/ichorCNA/plot.log"	
	shell:
		"Rscript scripts/writeTable.R ./"

rule install_readCounter:
	output:
		readCounter=config["readCounterScript"]
	log:
		"logs/install.log"
	shell:
		"./install_hmm_copy.sh"


rule read_counter:
	input:
		readCounter=config["readCounterScript"]
	output:
		"readDepth/{samples}.bin{binSize}.wig"		
	params:
		readCounter=config["readCounterScript"],
		binSize=config["binSize"],
		qual="20",
		chrs=config["chrs"],
		suffix = config["BAM_SUFFIX"],
		pathBAM = config["PATH_BAM"],
	resources:
		mem=4
	log:
		"logs/readDepth/{samples}.bin{binSize}.log"
	shell:
		"{params.readCounter} {params.pathBAM}{wildcards.samples}{params.suffix} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}"

rule ichorCNA:
	input:
		tum="readDepth/{samples}.bin" + str(config["binSize"]) + ".wig",
	output:
		cna="output_ichorCNA/ichorCNA/{samples}/{samples}.cna.seg",
		report=report("output_ichorCNA/ichorCNA/{samples}/{samples}/{samples}_genomeWide.pdf", category="genomeWide")
	params:
		outDir="output_ichorCNA/ichorCNA/{samples}/",
		rscript=config["ichorCNA_rscript"],
		id="{samples}",
		ploidy=config["ichorCNA_ploidy"],
		normal=config["ichorCNA_normal"],
		gcwig=config["ichorCNA_gcWig"],
		mapwig=config["ichorCNA_mapWig"],
		normalpanel=config["ichorCNA_normalPanel"],
		estimateNormal=config["ichorCNA_estimateNormal"],
		estimatePloidy=config["ichorCNA_estimatePloidy"],
		estimateClonality=config["ichorCNA_estimateClonality"],
		scStates=config["ichorCNA_scStates"],
		maxCN=config["ichorCNA_maxCN"],
		includeHOMD=config["ichorCNA_includeHOMD"],
		chrs=config["ichorCNA_chrs"],
		chrTrain=config["ichorCNA_chrTrain"],
		genomeBuild=config["ichorCNA_genomeBuild"],
		genomeStyle=config["ichorCNA_genomeStyle"],
		centromere=config["ichorCNA_centromere"],
		fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
		minMapScore=config["ichorCNA_minMapScore"],
		maxFracGenomeSubclone=config["ichorCNA_maxFracGenomeSubclone"],
		maxFracCNASubclone=config["ichorCNA_maxFracCNASubclone"],
		exons=config["ichorCNA_exons"],
		txnE=config["ichorCNA_txnE"],
		txnStrength=config["ichorCNA_txnStrength"],
		plotFileType=config["ichorCNA_plotFileType"],
		plotYlim=config["ichorCNA_plotYlim"],
		altFracThreshold=config["ichorCNA_altFracThreshold"], 
		libdir=config["ichorCNA_libdir"]
	conda:
		"envs/envIchorCNA.yaml"	
	log:
		"logs/ichorCNA/{samples}.log"	
	shell:
		"Rscript {params.rscript} --id {params.id} --libdir {params.libdir} --WIG {input.tum} --gcWig {params.gcwig} --mapWig {params.mapwig} --normalPanel {params.normalpanel} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --centromere {params.centromere} --exons.bed {params.exons} --txnE {params.txnE} --txnStrength {params.txnStrength} --minMapScore {params.minMapScore} --fracReadsInChrYForMale {params.fracReadsChrYMale} --altFracThreshold {params.altFracThreshold} --maxFracGenomeSubclone {params.maxFracGenomeSubclone} --maxFracCNASubclone {params.maxFracCNASubclone} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"
"""