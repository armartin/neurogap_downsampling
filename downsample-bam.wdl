## Copyright Broad Institute, 2018
##
## This WDL pipeline downsamples a bam file to a given coverage
##
## Requirements/expectations :
## - A file list with two columns: a BAM and its coverage
## - BAM files should be Analysis-ready for a single sample (as identified in RG:SM)
## - Coverage as indicated is Mean Coverage (Raw) Full
##
## Output :
## - A downsampled BAM file
##
## Cromwell version support
## - Successfully tested on v36
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script.

# WORKFLOW DEFINITION

workflow DownsampleBam {
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  Float target_depth
  File bam_files_list

  Array[Array[String]] bam_files = read_tsv(bam_files_list)

  String? gatk_docker_override
  String gatk_docker = select_first([gatk_docker_override, "broadinstitute/gatk:latest"])
  String? gitc_docker_override
  String gitc_docker = select_first([gitc_docker_override, "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"])
  String? samtools_path_override
  String samtools_path = select_first([samtools_path_override, "samtools"])

  Int? preemptible_attempts

  scatter (idx in range(length(bam_files))) {

    String input_bam = bam_files[idx][0]
    Float cov = bam_files[idx][1]

    #is the input a cram file?
    Boolean is_cram = sub(basename(input_bam), ".*\\.", "") == "cram"
    String sample_basename = if is_cram then basename(input_bam, ".cram") else basename(input_bam, ".bam")

    if ( is_cram ) {
      call CramToBamTask {
            input:
              input_cram = input_bam,
              sample_name = sample_basename,
              ref_dict = ref_dict,
              ref_fasta = ref_fasta,
              ref_fasta_index = ref_fasta_index,
              docker = gitc_docker,
              samtools_path = samtools_path
      }
    }

    call DownsampleBamTask {
      input:
        input_bam = select_first([CramToBamTask.output_bam, input_bam]),
        ref_fasta = ref_fasta,
        probability = target_depth / cov,
        output_name = basename(input_bam, ".bam")  + "." + target_depth + ".bam",
        output_index = basename(input_bam, ".bam")  + "." + target_depth + ".bai",
        docker = gatk_docker,
        preemptible_attempts = preemptible_attempts
    } 
  }
 }


task CramToBamTask {
  # Command parameters
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File input_cram
  String sample_name

  # Runtime parameters
  String docker
  Int? machine_mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts
  String samtools_path

  Float output_bam_size = size(input_cram, "GB") / 0.40
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Int disk_size = ceil(size(input_cram, "GB") + output_bam_size + ref_size) + 25

  command {
    set -e
    set -o pipefail

    ${samtools_path} view -h -T ${ref_fasta} ${input_cram} |
    ${samtools_path} view -b -o ${sample_name}.bam -
    ${samtools_path} index -b ${sample_name}.bam
    mv ${sample_name}.bam.bai ${sample_name}.bai
  }
  runtime {
    docker: docker
    memory: select_first([machine_mem_gb, 15]) + " GB"
    disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
 }
  output {
    File output_bam = "${sample_name}.bam"
    File output_bai = "${sample_name}.bai"
  }
}

# Downsamples bam/sam file
task DownsampleBamTask {
    # Command parameters
    File input_bam
    File ref_fasta
    Float probability
    String output_name
    String output_index

    # Runtime parameters
    String docker
    Int? machine_mem_gb
    Int machine_mem = select_first([machine_mem_gb, 18]) # 7
    Int? disk_space_gb
    Int? preemptible_attempts

    Int command_mem_gb = machine_mem - 1
    Int disk_size = ceil( size(input_bam, "GB") * 2) + 100 # 20

    command {
       set -e

       gatk DownsampleSam \
       --REFERENCE_SEQUENCE=${ref_fasta} \
       --INPUT=${input_bam} \
       --OUTPUT=${output_name} \
       --PROBABILITY=${probability} \
       --CREATE_INDEX=true
    }

    runtime {
        docker: docker
        memory: machine_mem + " GB"
        disks: "local-disk " + select_first([disk_space_gb, disk_size]) + " HDD"
        preemptible: select_first([preemptible_attempts, 3])
    }

    output {
        File output_bam = "${output_name}"
        File output_bai = "${output_index}"
    }
}

