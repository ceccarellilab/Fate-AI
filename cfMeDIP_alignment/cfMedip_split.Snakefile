my_basedir = workflow.current_basedir

# --- functions to create lists of snakemake target files ---

import glob
import re
import os

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
			print(glob.glob('%s*.fastq.gz' % path))
			samples = set([re.sub("_R[1-2]_001.fastq.gz", "", i) for i in glob.glob('%s*.fastq.gz' % path)]) 
			samples = set([re.sub(path, "", i) for i in samples]) 
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

def getNumSampleExp(exper):
	sizeExp = 0
	for samptype in config["SAMPTYPE"][exper]:
		path = "%s/%s/%s/" % (config["DATDIR"],exper,samptype)

		samples = set([re.sub("_R[1-2]_001.fastq.gz", "", i) for i in glob.glob('%s*.fastq.gz' % path)]) 
		
		samples = set([re.sub(path, "", i) for i in samples]) 
		
		print(samples)
		sizeExp = sizeExp + len(samples)
	return sizeExp


# --- snakemake rules ---


# rule to run the whole workflow
rule workflow:
	input:
		get_sortedBAM
		#get_multiqc
		

rule fastqc_initial:
    input:
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
        bowtie2_path = config["BOWTIE2PATH"],
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