configfile: "/home2/adefalco/LiquidBiopsy_Data/cfMedip/pipeline/config_Mel.yaml"

script_folder = "/home/tufano/projects/cfMeDIP/pipelines/MedipPipeline_Antonio"

my_basedir = workflow.current_basedir

# --- functions to create lists of snakemake target files ---

import glob
import re
import os

def get_fastqFile(wildcards):
	ls = []
	if config['PAIRED']==True:
		ls.append(config['PATH_FASTQ'] + wildcards.sample + "_"+ config['FASTQ_SUFFIXES']["1_SUFFIX"]) 
		ls.append(config['PATH_FASTQ'] + wildcards.sample + "_"+ config['FASTQ_SUFFIXES']["2_SUFFIX"])
	else:
		ls.append(config['PATH_FASTQ'] + wildcards.sample + config['FASTQ_SUFFIX']) 

	return ls

def get_sortedBAM(wildcards):
	ls = []
	for exp in config["EXP"]:
		for samptype in config["SAMPTYPE"][exp]:
			# print(samptype)
			for tool in ["fastqc","trimgalore", "bowtie2", "sortedbam", "sortedbam_dup", "multiqc"]:
				path = "%s/%s/%s/%s/" % (config["OUTDIR"],tool,exp,samptype)
				if not os.path.exists(path):
					os.makedirs(path)
			path = "%s/%s/%s/" % (config["DATDIR"],exp,samptype)
			print(path)
			#os.chdir(path)
			#samples = set([re.sub("_R[1-2]_001.fastq.gz", "", i) for i in glob.glob('*.fastq.gz')]) 
			print(glob.glob('%s*.fastq.gz' % path))
			samples = set([re.sub("_R[1-2]_001.fastq.gz", "", i) for i in glob.glob('%s*.fastq.gz' % path)]) 
			print(samples)
			samples = set([re.sub(path, "", i) for i in samples]) 

			# print("OLD")
			print(samples)

			# print("glob.glob")
			# print(glob.glob('%s*.fastq.gz' % path))
			
			samples_old = samples.copy()
			#r for sample in samples_old:
			#r	if sample in config["EXCLUDE_SAMPLES"]:
			#r		samples.remove(sample)

			# print("NEW")
			print(samples)

			#samples = set([re.sub("_[1-2].fq.gz", "", i) for i in glob.glob('*.fq.gz')]) 
			for sample in samples:
				ls.append("%s/sortedbam_dup/%s/%s/%s.sorted.bam" % (config["OUTDIR"],exp,samptype,sample))
	return ls

print(get_sortedBAM(""))

def get_multiqc(wildcards):
	ls = []
	for exp in config["EXP"]:
		for samptype in config["SAMPTYPE"][exp]:
			ls.append("%s/multiqc/%s/%s/fastqc_samtools_flagstat.html" % (config["OUTDIR"],exp,samptype))
	return ls

def get_MEDIPobj(wildcards):
	ls = []
	for exp in config["EXP"]:
		for samptype in config["SAMPTYPE"][exp]:
			ls.append("%s/MEDIPS_%s/%s/medip.%s.rds" % (config["OUTDIR"],config["WS"],exp,samptype))
	return ls

def getNumSampleExp(exper):
	sizeExp = 0
	for samptype in config["SAMPTYPE"][exper]:
		path = "%s/%s/%s/" % (config["DATDIR"],exper,samptype)
		#os.chdir(path)

		samples = set([re.sub("_R[1-2]_001.fastq.gz", "", i) for i in glob.glob('%s*.fastq.gz' % path)]) 
		# samples = set([re.sub("_R[1-2]_001.fastq.gz", "", i) for i in glob.glob('*.fastq.gz')]) 
		
		samples = set([re.sub(path, "", i) for i in samples]) 

		# print("n samples before exclusion")
		# print(len(samples))
		samples_old = samples.copy()
		#r for sample in samples_old:
		#r	if sample in config["EXCLUDE_SAMPLES"]:
		#r 		samples.remove(sample)

		#samples = set([re.sub("_[1-2].fq.gz", "", i) for i in glob.glob('*.fq.gz')]) 
		print(samples)
		sizeExp = sizeExp + len(samples)
	return sizeExp

def get_class(wildcards):
	ls = []
	for exp in config["EXP"]:
		if config["SAMPTYPE"][exp][0] == "control":
#		if config["SAMPTYPE"][exp][0] == "CC":
			lab_control = config["SAMPTYPE"][exp][0]
			lab_tumor = config["SAMPTYPE"][exp][1]
		else:
			lab_control = config["SAMPTYPE"][exp][1]
			lab_tumor = config["SAMPTYPE"][exp][0]	
		if config["MODELSELECTION"]=="LeaveOneOut":
			for i in range(1,getNumSampleExp(exp)+1):
				ls.append("%s/MEDIPS_%s/%s/results/sampleprob_table_%s_%s_%s_%s.txt" % (config["OUTDIR"],config["WS"],exp,lab_tumor,lab_control,config["METHOD"], i))
		elif config["MODELSELECTION"]=="ShuffleSplit":
			for i in range(1,config["N_SPLITS"]+1):
				ls.append("%s/MEDIPS_%s/%s/results/sampleprob_table_%s_%s_%s_%s.txt" % (config["OUTDIR"],config["WS"],exp,lab_tumor,lab_control,config["METHOD"], i))	
	return ls

def get_ROC(wildcards):
	ls = []
	for exp in config["EXP"]:
		ls.append("%s/MEDIPS_%s/%s/results/ROC.png" % (config["OUTDIR"],config["WS"],exp,))	
	return ls

# --- snakemake rules ---

#shell

# rule to run the whole workflow
rule workflow:
	input:
		get_sortedBAM
		#get_multiqc
		#get_ROC
		

rule fastqc_initial:
    input:
        #get_fastqFile
        r1 = lambda wildcards: config['PATH_FASTQ'] + wildcards.sample + "_" + config['FASTQ_SUFFIXES']["1_SUFFIX"],
        r2 = lambda wildcards: config['PATH_FASTQ'] + wildcards.sample + "_" + config['FASTQ_SUFFIXES']["2_SUFFIX"]
    output:
        r1_html = "{outdir}/fastqc/{group}/{sampletype}/{sample}_R1_001_fastqc.html",
        r2_html = "{outdir}/fastqc/{group}/{sampletype}/{sample}_R2_001_fastqc.html"
    threads: 1
    resources:
        mem_mb=10240*2,
        total_cpus=1,
    params:
        outdir = "{outdir}/fastqc/{group}/{sampletype}"
    shell:
        """
        mkdir -p {params.outdir}
        if [ ! -s {output.r1_html} ]; then
            fastqc -o {params.outdir} {input.r1}
        fi
        if [ ! -s {output.r2_html} ]; then
            fastqc -o {params.outdir} {input.r2}
        fi
        """

rule trim_galore:
    input:
        #get_fastqFile
        r1 = lambda wildcards: config['PATH_FASTQ'] + wildcards.sample + "_" + config['FASTQ_SUFFIXES']["1_SUFFIX"],
        r2 = lambda wildcards: config['PATH_FASTQ'] + wildcards.sample + "_" + config['FASTQ_SUFFIXES']["2_SUFFIX"]
    threads: 1
    resources:
        mem_mb=10240*2,
        total_cpus=1
    output:
        trimmed_r1 = "{outdir}/trimgalore/{group}/{sampletype}/{sample}_R1_001_val_1.fq.gz",
        trimmed_r2 = "{outdir}/trimgalore/{group}/{sampletype}/{sample}_R2_001_val_2.fq.gz",
        fastqc_r1 = "{outdir}/fastqc/{group}/{sampletype}/{sample}_R1_001_val_1_fastqc.html",
        fastqc_r2 = "{outdir}/fastqc/{group}/{sampletype}/{sample}_R2_001_val_2_fastqc.html"
    params:
        out_trim = "{outdir}/trimgalore/{group}/{sampletype}",
        out_fastqc = "{outdir}/fastqc/{group}/{sampletype}",
        path_cutadapt = config['CUTADAPTPATH'],
        path_grimgalore = config['TRIMGALOREPATH']
    shell:
        """
        mkdir -p {params.out_trim} {params.out_fastqc}
        if [ ! -s {output.trimmed_r2} ]; then
          {params.path_grimgalore} --fastqc --fastqc_args "-o {params.out_fastqc}" \
            -o {params.out_trim} --paired {input.r1} {input.r2} --path_to_cutadapt {params.path_cutadapt}
        fi
        """

rule bowtie2_map:
    input:
        r1 = "{outdir}/trimgalore/{group}/{sampletype}/{sample}_R1_001_val_1.fq.gz",
        r2 = "{outdir}/trimgalore/{group}/{sampletype}/{sample}_R2_001_val_2.fq.gz"
    output:
        bam = "{outdir}/bowtie2/{group}/{sampletype}/{sample}.bam"
    threads: 40
    resources:
        mem_mb=10240*5,
        total_cpus=40
    params:
        index = config["INDEX"],
        genome = config["GENOME_index"],   # set this in config or params from shell env
        sam = "{outdir}/bowtie2/{group}/{sampletype}/{sample}.sam",
        #logfold = "{outdir}/bowtie2_log/{sample}.log",
        bowtie2_path = config["BOWTIE2PATH"],
        #bowtie2_path = "bowtie2",
        n_threads = config["NUM_THREADS"],
        samtools_path = config["SAMTOOLSPATH"]
    log:
        "{outdir}/bowtie2/{group}/{sampletype}/{sample}.log"
    shell:
        """
        echo "bowtie2"
        export BOWTIE2_INDEXES={params.index}
        mkdir -p $(dirname {output.bam})
        {params.bowtie2_path} -x {params.genome} -p {params.n_threads} -S {params.sam} -1 {input.r1} -2 {input.r2}
        {params.samtools_path} view -S -b {params.sam} > {output.bam}
        rm {params.sam}
        """

rule sort_index_bam:
    input:
        bam = "{outdir}/bowtie2/{group}/{sampletype}/{sample}.bam"
    output:
        sorted_bam = "{outdir}/sortedbam/{group}/{sampletype}/{sample}.sorted.bam",
        bai = "{outdir}/sortedbam/{group}/{sampletype}/{sample}.sorted.bam.bai"
    threads: 1
    resources:
        mem_mb=10240*2,
        total_cpus=1
    params:
        samtools_path = config["SAMTOOLSPATH"]
    shell:
        """
        mkdir -p $(dirname {output.sorted_bam})
        if [ ! -s {output.sorted_bam} ] || ls {output.sorted_bam}.tmp.* 1> /dev/null 2>&1; then
          rm -f {output.sorted_bam}*
          {params.samtools_path} sort {input.bam} -o {output.sorted_bam}
        fi
        if [ ! -s {output.bai} ]; then
          {params.samtools_path} index {output.sorted_bam}
        fi
        """

rule mark_duplicates:
    input:
        sorted_bam = "{outdir}/sortedbam/{group}/{sampletype}/{sample}.sorted.bam"
    output:
        dedup_bam = "{outdir}/sortedbam_dup/{group}/{sampletype}/{sample}.sorted.bam",
        dedup_bai = "{outdir}/sortedbam_dup/{group}/{sampletype}/{sample}.sorted.bam.bai"
    threads: 1
    resources:
        mem_mb=10240*2,
        total_cpus=1
    params:
        samtools_path = config["SAMTOOLSPATH"]
    shell:
        """
        mkdir -p $(dirname {output.dedup_bam})
        rm -f {output.dedup_bam}.*  # cleanup old files

        # sort by name
        {params.samtools_path} sort -n -o {output.dedup_bam}.namesort.bam {input.sorted_bam}

        # fixmate
        {params.samtools_path} fixmate -m {output.dedup_bam}.namesort.bam {output.dedup_bam}.fixmate.bam
        rm {output.dedup_bam}.namesort.bam

        # sort by position
        {params.samtools_path} sort -o {output.dedup_bam}.positionsort.bam {output.dedup_bam}.fixmate.bam
        rm {output.dedup_bam}.fixmate.bam

        # mark duplicates and remove
        {params.samtools_path} markdup -r -s {output.dedup_bam}.positionsort.bam {output.dedup_bam}
        rm {output.dedup_bam}.positionsort.bam

        # index dedup
        {params.samtools_path} index {output.dedup_bam}
        """