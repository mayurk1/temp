// Samples
ch_fastq_star = Channel.fromFilePairs(params.samples_path).into{ch_fastq; ch_star}

// References 
star_dir = Channel.value(file(params.hg38STAR_path))
ref_fa = Channel.value(file(params.ref_fa_path))
ref_flat = Channel.value(file(params.ref_flat_path))
ref_ribo_intervals = Channel.value(file(params.ref_ribo_intervals_path))

// Script variables
EXPECTED_PAIR_ORIENTATIONS = params.EXPECTED_PAIR_ORIENTATIONS
STRAND_SPECIFICITY = params.STRAND_SPECIFICITY

process fastqc{
	tag "FASTQC ${sample}"
	cpus 2
	memory '2 GB'
	time 1.hour
	module 'java-jdk/1.10.0_1'
	module 'fastqc/0.11.7'

	input:
		tuple val(sample), file(fq_files) from ch_fastq

	output:
		path("${sample}.R*.html") into fastqc_html
		path("${sample}.R*.zip") into fastqc_zip

	script:
		"""
		fastqc ${fq_files}
		"""
}

process STAR{
	tag "STAR ${sample}"
	cpus 12
	memory '48 GB' 
	time 1.hour
	module 'gcc/6.2.0'
	module 'STAR/2.6.1d'

	input:
		tuple val(sample), file(fq_files) from ch_star
		file genome from star_dir

	output:
		path("${sample}.*.bam") into aligned_bam
		path("${sample}.Log.final.out") into star_log
		path("${sample}.ReadsPerGene.out.tab") into star_tab

	script:
		"""
		STAR \
		--runThreadN 12 \
		--genomeDir ${genome} \
		--readFilesIn ${fq_files} \
		--outSAMtype BAM Unsorted \
		--outSAMunmapped Within \
		--outFileNamePrefix ${sample}. \
		--outSAMattrRGline ID:${sample}\tPU:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:${sample} \
		--outSAMattributes NH HI AS nM NM \
		--quantMode GeneCounts
		"""
}

process sortAndIndex{
	tag "Sort and Index ${bam.getName()}"
	cpus 6
	memory '6 GB'
	time 1.hour
	module 'gcc/6.2.0'
	module 'samtools/1.10'
	publishDir "results/"

	input:
		path bam from aligned_bam

	output:
		file("*${bam.getSimpleName()}.sorted.bam") into sorted_bam
		file("*${bam.getSimpleName()}.sorted.bam.bai") into sorted_bai

	script:
		"""
		samtools sort \
		-O bam \
		-o ${bam.getSimpleName()}.sorted.bam \
		-@ 6 \
		${bam}

		samtools index ${bam.getSimpleName()}.sorted.bam
		"""
}

process validate{
	tag "Validate ${sortedBam.getName()}"
	cpus 2
	memory '4 GB'
	time 24.min
	module 'gcc/6.2.0'
	module 'java-jdk/1.8.0_92'
	module 'picard/2.18.29'
	module 'R/3.6.1'

	input:
		file sortedBam from sorted_bam
		file bai from sorted_bai
		file ref_fa from ref_fa
		file ref_flat from ref_flat
		file ref_ribo_intervals from ref_ribo_intervals
		
	output:
		path("*${sortedBam.getName()}.validation_metrics") into val_metrics
		path("*${sortedBam.getName()}.gc_bias_metrics")	into bias_metrics
		path("*${sortedBam.getName()}.gc_bias_chart.pdf") into bias_chart
		path("*${sortedBam.getName()}.gc_bias_summary") into bias_sum
		path("*${sortedBam.getName()}.alignment_summary") into aln_sum
		path("*${sortedBam.getName()}.rna_seq_metrics")into seq_metrics

	script:
		"""
		# Validate BAM File
		java -Xmx4g -jar \${PICARD} \
			ValidateSamFile \
			INPUT=${sortedBam} \
			OUTPUT=${sortedBam.getName()}.validation_metrics \
			REFERENCE_SEQUENCE=${ref_fa} \
			MODE=VERBOSE \
			VALIDATE_INDEX=true \
			INDEX_VALIDATION_STRINGENCY=EXHAUSTIVE
			
		# GC Bias Metrics
		java -Xmx4g -jar \${PICARD} \
			CollectGcBiasMetrics \
			INPUT=${sortedBam} \
			OUTPUT=${sortedBam.getName()}.gc_bias_metrics \
			REFERENCE_SEQUENCE=${ref_fa} \
			CHART_OUTPUT=${sortedBam.getName()}.gc_bias_chart.pdf \
			SUMMARY_OUTPUT=${sortedBam.getName()}.gc_bias_summary \
			IS_BISULFITE_SEQUENCED=false \
			METRIC_ACCUMULATION_LEVEL=ALL_READS \
			METRIC_ACCUMULATION_LEVEL=READ_GROUP \
			ALSO_IGNORE_DUPLICATES=true

		# Alignment Summary		
		java -Xmx4g -jar \${PICARD} \
			CollectAlignmentSummaryMetrics \
			INPUT=${sortedBam} \
			OUTPUT=${sortedBam.getName()}.alignment_summary \
			REFERENCE_SEQUENCE=${ref_fa} \
			MAX_INSERT_SIZE=100000 \
			EXPECTED_PAIR_ORIENTATIONS=${EXPECTED_PAIR_ORIENTATIONS} \
			IS_BISULFITE_SEQUENCED=false \
			METRIC_ACCUMULATION_LEVEL=ALL_READS \
			METRIC_ACCUMULATION_LEVEL=READ_GROUP

		# RNA Seq Metrics
		java -jar -Xmx4g \${PICARD} \
			CollectRnaSeqMetrics \
			INPUT=${sortedBam} \
			OUTPUT=${sortedBam.getName()}.rna_seq_metrics \
			REF_FLAT=${ref_flat} \
			STRAND_SPECIFICITY=${STRAND_SPECIFICITY} \
			RIBOSOMAL_INTERVALS=${ref_ribo_intervals} \
			METRIC_ACCUMULATION_LEVEL=ALL_READS \
			METRIC_ACCUMULATION_LEVEL=READ_GROUP \
			ASSUME_SORTED=true
		"""
}

process multiQC{
	tag "MultiQC"
	cpus 1
	memory '1 GB'
	time 24.min
	module 'gcc/6.2.0'
	module 'multiqc/1.8.0'
	publishDir "results/"

	input:
		file fastqc_html from fastqc_html.collect()
		file fastqc_zip from fastqc_zip.collect()

		file star_log from star_log.collect()
		file star_tab from star_tab.collect()

		file val_metrics from val_metrics.collect()
		file bias_metrics from bias_metrics.collect()
		file bias_chart from bias_chart.collect()
		file bias_sum from bias_sum.collect()
		file aln_sum from aln_sum.collect()
		file seq_metrics from seq_metrics.collect()

	output:
		path("multiqc_report.html")

	script:
		"""
		multiqc . -f
		"""
}
