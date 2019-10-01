workflow LowCovImputation {
	File input_vcf
	String refined_base
	String imputed_base
	Boolean run_refine_genotypes
	Boolean run_impute_genotypes

	File intervals_file
	Array[Array[String]] intervals = read_tsv(intervals_file)

	File intervals_chrom
	Array[String] intervals_chr = read_lines(intervals_chrom)

	File ref_file
	Array[Array[String]] ref_files = read_tsv(ref_file)

    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File eval_interval_list
	File dbsnp_vcf
	File dbsnp_vcf_index

	String? gatk_docker_override
	String gatk_docker = select_first([gatk_docker_override, "broadinstitute/gatk:4.1.0.0"])
	String? gatk_path_override
	String gatk_path = select_first([gatk_path_override, "/gatk/gatk"])
	String? beagle_docker_override
	String beagle_docker = select_first([beagle_docker_override, "armartin/beagle_refine_impute:0.6"])
	String? bcftools_path_override
	String bcftools_path = select_first([bcftools_path_override, "/bcftools/bcftools"])
	String? bgzip_path_override
	String bgzip_path = select_first([bgzip_path_override, "/htslib/bgzip"])
	String? tabix_path_override
	String tabix_path = select_first([tabix_path_override, "/htslib/tabix"])
	String? beagle4_path_override
	String beagle4_path = select_first([tabix_path_override, "/beagle.27Jan18.7e1.jar"])
	String? beagle5_path_override
	#String beagle5_path = select_first([tabix_path_override, "/beagle.12Jul19.0df.jar"]) # beagle 5.0
	String beagle5_path = select_first([tabix_path_override, "/beagle.21Sep19.ec3.jar"]) # beagle 5.1

	Int? small_disk_override
    Int small_disk = select_first([small_disk_override, "100"])
    Int? medium_disk_override
    Int medium_disk = select_first([medium_disk_override, "200"])
    Int? large_disk_override
    Int large_disk = select_first([large_disk_override, "300"])
    Int? huge_disk_override
    Int huge_disk = select_first([huge_disk_override, "400"])

    String? preemptible_tries_override
    Int preemptible_tries = select_first([preemptible_tries_override, "3"])

	if (run_refine_genotypes) {
        scatter (idx in range(length(intervals))) {
            String refine_phase_interval = intervals[idx][0] + ':' + intervals[idx][1] + '-' + intervals[idx][2]
            String refine_interval_trim = intervals[idx][0] + ':' + intervals[idx][3] + '-' + intervals[idx][4]
            Int refine_chrom_index = sub(intervals[idx][0], "chr", "")

            String refine_reference = ref_files[(refine_chrom_index-1)][0]
            String refine_recomb_map = ref_files[(refine_chrom_index-1)][2]

            String input_vcf_idx = input_vcf + '.tbi'

            call RefineGenotypeTask {
                input:
                    input_vcf = input_vcf,
                    input_vcf_idx = input_vcf_idx,
                    phased_ref_vcf = refine_reference,
                    interval = refine_phase_interval,
                    interval_trim = refine_interval_trim,
                    recomb_map = refine_recomb_map,
                    refined_base = basename(refined_base),
                    ref_fasta_index = ref_fasta_index,
                    docker = beagle_docker,
                    beagle_path = beagle4_path,
                    bcftools_path = bcftools_path,
                    bgzip_path = bgzip_path,
                    tabix_path = tabix_path
            }
        }

        # Gather the VCF shards, grabbing the middle of each shard (minus the ends)
        call GatherVcfs as GatherRefinedVcf {
            input:
                input_vcfs_fofn = write_lines(RefineGenotypeTask.refined_vcf),
                output_vcf_name = basename(refined_base) + ".vcf.gz",
                disk_size = huge_disk,
                docker = gatk_docker,
                gatk_path = gatk_path,
                preemptible_tries = preemptible_tries
        }

        # Collect metrics on gathered VCF
        call CollectVariantCallingMetrics as CollectMetricsRefinedVcf {
            input:
                input_vcf = GatherRefinedVcf.output_vcf,
                input_vcf_index = GatherRefinedVcf.output_vcf_index,
                metrics_filename_prefix = basename(refined_base),
                dbsnp_vcf = dbsnp_vcf,
                dbsnp_vcf_index = dbsnp_vcf_index,
                interval_list = eval_interval_list,
                ref_dict = ref_dict,
                disk_size = large_disk,
                docker = gatk_docker,
                gatk_path = gatk_path,
                preemptible_tries = preemptible_tries
        }
    }
		
	if (run_impute_genotypes) {
        scatter (idx in range(length(intervals_chr))) {
            String impute_chr_interval = intervals_chr[idx]
            # String impute_phase_interval = intervals[idx][0] + ':' + intervals[idx][1] + '-' + intervals[idx][2]
            # String impute_interval_trim = intervals[idx][0] + ':' + intervals[idx][3] + '-' + intervals[idx][4]
            Int impute_chrom_index = sub(intervals_chr[idx], "chr", "")

            String impute_reference = ref_files[(impute_chrom_index-1)][1]
            String impute_recomb_map = ref_files[(impute_chrom_index-1)][2]

            call ImputeTask {
                input:
                    input_vcf = if (run_refine_genotypes) then GatherRefinedVcf.output_vcf else input_vcf,
                    phased_ref_vcf = impute_reference,
                    #interval = impute_phase_interval,
                    #interval_trim = impute_interval_trim,
                    interval = impute_chr_interval,
                    recomb_map = impute_recomb_map,
                    imputed_base = basename(imputed_base),
                    ref_fasta_index = ref_fasta_index,
                    docker = beagle_docker,
                    beagle_path = beagle5_path,
                    bcftools_path = bcftools_path,
                    bgzip_path = bgzip_path,
                    tabix_path = tabix_path
            }
        }

        # Gather the VCF shards
        call GatherVcfs as GatherImputedVcf {
            input:
                input_vcfs_fofn = write_lines(ImputeTask.imputed_vcf),
                output_vcf_name = basename(imputed_base) + ".vcf.gz",
                disk_size = huge_disk,
                docker = gatk_docker,
                gatk_path = gatk_path,
                preemptible_tries = preemptible_tries
        }

        # Collect metrics on gathered VCF
        call CollectVariantCallingMetrics as CollectMetricsImputedVcf {
            input:
                input_vcf = GatherImputedVcf.output_vcf,
                input_vcf_index = GatherImputedVcf.output_vcf_index,
                metrics_filename_prefix = basename(imputed_base),
                dbsnp_vcf = dbsnp_vcf,
                dbsnp_vcf_index = dbsnp_vcf_index,
                interval_list = eval_interval_list,
                ref_dict = ref_dict,
                disk_size = large_disk,
                docker = gatk_docker,
                gatk_path = gatk_path,
                preemptible_tries = preemptible_tries
        }

        # Outputs that will be retained when execution is complete
        output {
            File? refined_vcf = GatherRefinedVcf.output_vcf
            File? refined_vcf_index = GatherRefinedVcf.output_vcf_index
            File? refined_detail_metrics_file = CollectMetricsRefinedVcf.detail_metrics_file
            File? refined_summary_metrics_file = CollectMetricsRefinedVcf.summary_metrics_file
            File? imputed_vcf = GatherImputedVcf.output_vcf
            File? imputed_vcf_index = GatherImputedVcf.output_vcf_index
            File? imputed_detail_metrics_file = CollectMetricsImputedVcf.detail_metrics_file
            File? imputed_summary_metrics_file = CollectMetricsImputedVcf.summary_metrics_file
        }
    }
}


task RefineGenotypeTask {
	# Command parameters
	File input_vcf
	File input_vcf_idx
	File phased_ref_vcf
	File recomb_map
	File ref_fasta_index
	String interval
	String interval_trim
	String refined_base
	String beagle_path
	String bcftools_path
	String bgzip_path
	String tabix_path

	# Runtime parameters
	String docker
	Int? machine_mem_gb
	Int? disk_space_gb
	Boolean use_ssd = false
	Int? preemptible_attempts

	Int disk_size = ceil(size(input_vcf, "GB") * 2) + 10
# was 8G
	command <<<
        set -e

		java -Xmx16g -Xms16g \
		-jar ${beagle_path} \
		gl=${input_vcf} \
		ref=${phased_ref_vcf} \
		map=${recomb_map} \
		chrom=${interval} \
		modelscale=2 \
		niterations=0 \
		impute=false \
		lowmem=true \
		out=${refined_base}.${interval}

        ${bcftools_path} reheader \
        -f ${ref_fasta_index} \
        ${refined_base}.${interval}.vcf.gz \
        -o reheader.vcf.gz

		${bcftools_path} index \
		reheader.vcf.gz

		${bcftools_path} convert \
		-r ${interval_trim} \
		-o ${refined_base}.${interval_trim}.vcf reheader.vcf.gz

		${bgzip_path} ${refined_base}.${interval_trim}.vcf

        ${tabix_path} ${refined_base}.${interval_trim}.vcf.gz
	>>>
	runtime {
		docker: docker
		memory: select_first([machine_mem_gb, 20]) + " GB" # was 11G
		#cpu: 2
		disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
		preemptible: select_first([preemptible_attempts, 5])
 }
	output {
		File refined_vcf = "${refined_base}.${interval_trim}.vcf.gz"
		File refined_vcf_index = "${refined_base}.${interval_trim}.vcf.gz.tbi"
		File refined_vcf_log = "${refined_base}.${interval}.log"
	}
}

task ImputeTask {
	# Command parameters
	File input_vcf
	File phased_ref_vcf
	File recomb_map
	File ref_fasta_index
	String interval
	#String interval_trim
	String imputed_base
	String beagle_path
	String bcftools_path
	String bgzip_path
	String tabix_path

	# Runtime parameters
	String docker
	Int? machine_mem_gb
	Int? disk_space_gb
	Boolean use_ssd = false
	Int? preemptible_attempts

	Int disk_size = ceil(size(input_vcf, "GB") * 2) + 10

	command <<<
	    set -e

		java -Xmx8g -Xms8g \
		-jar ${beagle_path} \
		gt=${input_vcf} \
		ref=${phased_ref_vcf} \
		map=${recomb_map} \
		chrom=${interval} \
		out=${imputed_base}.${interval}

		${bcftools_path} reheader \
        -f ${ref_fasta_index} \
        ${imputed_base}.${interval}.vcf.gz \
        -o reheader.vcf.gz

		${bcftools_path} index \
		reheader.vcf.gz

		${bcftools_path} convert \
		-r ${interval} \
		-o ${imputed_base}.${interval}.vcf reheader.vcf.gz

		${bgzip_path} -f ${imputed_base}.${interval}.vcf

        ${tabix_path} ${imputed_base}.${interval}.vcf.gz
	>>>
	runtime {
		docker: docker
		memory: select_first([machine_mem_gb, 11]) + " GB"
		disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
		preemptible: select_first([preemptible_attempts, 3])
	}
	output {
		#File imputed_vcf = "${imputed_base}.${interval_trim}.vcf.gz"
		File imputed_vcf = "${imputed_base}.${interval}.vcf.gz"
		#File imputed_vcf_index = "${imputed_base}.${interval_trim}.vcf.gz.tbi"
		File imputed_vcf_index = "${imputed_base}.${interval}.vcf.gz.tbi"
		File imputed_vcf_log = "${imputed_base}.${interval}.log"
	}
}

task GatherVcfs {
	File input_vcfs_fofn
	String output_vcf_name
	String base_vcf_name = basename(output_vcf_name)
	String gatk_path

	String docker
	Int disk_size
	Int preemptible_tries

	command <<<
		set -e

		# Now using NIO to localize the vcfs but the input file must have a ".list" extension
		mv ${input_vcfs_fofn} inputs.list

		# --ignore-safety-checks makes a big performance difference so we include it in our invocation.
		# This argument disables expensive checks that the file headers contain the same set of
		# genotyped samples and that files are in order by position of first record.
		${gatk_path} --java-options "-Xmx6g -Xms6g" \
		GatherVcfsCloud \
		--ignore-safety-checks \
		--gather-type BLOCK \
		--input inputs.list \
		--output ${base_vcf_name}

		${gatk_path} --java-options "-Xmx6g -Xms6g" \
		IndexFeatureFile \
		--feature-file ${base_vcf_name}
	>>>
	runtime {
		docker: docker
		memory: "7 GB"
		cpu: "1"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible_tries
	}
	output {
		File output_vcf = "${output_vcf_name}"
		File output_vcf_index = "${output_vcf_name}.tbi"
	}
}

# this needs the index!
task CollectVariantCallingMetrics {
	File input_vcf
	File input_vcf_index

	String metrics_filename_prefix
	File dbsnp_vcf
	File dbsnp_vcf_index
	File interval_list
	File ref_dict

	String gatk_path
	String docker
	Int disk_size
	Int preemptible_tries

	command {
		${gatk_path} --java-options "-Xmx6g -Xms6g" \
			CollectVariantCallingMetrics \
			--INPUT ${input_vcf} \
			--DBSNP ${dbsnp_vcf} \
			--SEQUENCE_DICTIONARY ${ref_dict} \
			--OUTPUT ${metrics_filename_prefix} \
			--THREAD_COUNT 8 \
			--TARGET_INTERVALS ${interval_list}
	}
	output {
		File detail_metrics_file = "${metrics_filename_prefix}.variant_calling_detail_metrics"
		File summary_metrics_file = "${metrics_filename_prefix}.variant_calling_summary_metrics"
	}
	runtime {
		docker: docker
		memory: "7 GB"
		cpu: 2
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible_tries
	}
}