workflow BamToFastq {
	File input_bams_file
	Array[File] input_bams = read_lines(input_bams_file)

	#File input_bam # gs://neurogap/high_coverage/WGS/AAJO2J0L-02/v1/AAJO2J0L-02.1.0.bam

	String? gatk_docker_override
	String gatk_docker = select_first([gatk_docker_override, "broadinstitute/gatk:4.1.0.0"])
	String? gatk_path_override
	String gatk_path = select_first([gatk_path_override, "/gatk/gatk"])

	String? preemptible_tries_override
    Int preemptible_tries = select_first([preemptible_tries_override, "3"])

    scatter (input_bam in input_bams) {
        String output_base = basename(input_bam, ".bam")
        call BamToFastqTask {
            input:
                input_bam = input_bam,
                output_base = output_base,
                gatk_path = gatk_path,
                docker = gatk_docker,
                preemptible_attempts = preemptible_tries
        }
    }

}


task BamToFastqTask {
	# Command parameters
	File input_bam
	String output_base # AAJO2J0L-02.1.0
	String sample = sub(output_base, "\\.\\d+\\.0$", "") # AAJO2J0L-02.1.0
	String cov = sub(output_base, "^.*\\.(\\d+\\.0)", "$1") # 0
	String gatk_path

	# Runtime parameters
	String docker
	Int? machine_mem_gb
	Int? disk_space_gb
	Boolean use_ssd = false
	Int? preemptible_attempts

	#Int disk_size = ceil(size(input_bam, "GB") * 3) + 10
	Int disk_size = 200

	command <<<
        set -e

#        echo ${sample}
#        echo ${cov}

#		gatk SamToFastq \
#		-I=${input_bam} \
#		--INCLUDE_NON_PRIMARY_ALIGNMENTS=true \
#		--INCLUDE_NON_PF_READS=true \
#		--FASTQ=${sample}_${cov}_R1.fastq.gz \
#		--SECOND_END_FASTQ=${sample}_${cov}_R2.fastq.gz

        samtools sort \
        -n ${input_bam} \
        ${input_bam}.sorted

        ls -l ${input_bam}.sorted.bam

        bedtools bamtofastq \
        -i ${input_bam}.sorted.bam \
        -fq  ${sample}_${cov}_R1.fastq \
        -fq2 ${sample}_${cov}_R2.fastq

        gzip ${sample}_${cov}_R1.fastq
        gzip ${sample}_${cov}_R2.fastq
	>>>
	runtime {
		docker: docker
		memory: select_first([machine_mem_gb, 12]) + " GB" # was 8G
		#cpu: 2
		disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
		preemptible: select_first([preemptible_attempts, 5])
 }
	output {
		File fastq1 = "${sample}_${cov}_R1.fastq.gz"
		File fastq2 = "${sample}_${cov}_R2.fastq.gz"
	}
}