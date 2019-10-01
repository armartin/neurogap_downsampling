# NeuroGAP downsampling project

This repository contains the code generated and run to systematically compare low-coverage versus GWAS accuracy and sensitivity. Starting data consisted of high coverage whole genomes (target: 30X) from the 5 study sites included in this project. Reads or variants were then downsampled to correspond with a depth of coverage or sites from a GWAS array. Variant call set quality was assessed and performance was compared with respect to the high coverage genomes as a "truth" set.

---

## Processing steps
**For sequencing data:**
- Downsample cram files to specific depths (see workflow [downsample-bam.wdl](https://github.com/armartin/neurogap_downsampling/blob/master/downsample-bam.wdl) with corresponding example inputs [downsample-bam.json](https://github.com/armartin/neurogap_downsampling/blob/master/beagle-refine-impute.wdl))
- Run [HaplotypeCaller](https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/haplotypecaller-gvcf-gatk4-nio.wdl) per-sample with appropriate arguments [(see example inputs)](https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/haplotypecaller-gvcf-gatk4.hg38.wgs.inputs.json)
- Create a joint call set using [best practices](https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/joint-discovery-gatk4.wdl) with appropriate arguments [(see example inputs)](https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/joint-discovery-gatk4.hg38.wgs.inputs.json)

**For array data:**
- Subset full call set to sites on specific arrays using [Hail](https://hail.is/), as described in [extract_array_sites.py](https://github.com/armartin/neurogap_downsampling/blob/master/extract_array_sites.py)

**For both:**
- Run BEAGLE using this workflow ([beagle-refine-impute.wdl](https://github.com/armartin/neurogap_downsampling/blob/master/beagle-refine-impute.wdl)) and with corresponding example inputs ([beagle-refine-impute.json](https://github.com/armartin/neurogap_downsampling/blob/master/beagle-refine-impute.json)).
Differences between sequencing vs array analyses: 
- Use a 2-step approach for the low-coverage data consisting of genotyping refinement using genotype likelihoods rather than hard genotype calls (fields called GL rather than GT in VCF parlance). After genotypes are refined, use the hard calls to run imputation. BEAGLE v4.1 was used for refinement, and BEAGLE v5.1 was used for imputation. Different reference panels were used for these steps, as refinement is much more computationally intensive (refinement reference panel subsetting script here: [refine_reference.py](https://github.com/armartin/neurogap_downsampling/blob/master/refine_reference.py)).
- No point to refinement for array analyses, just run imputation using BEAGLE v5.1

## Analytical steps
- In addition to standard variant call metrics from GATK (gathered routinely as part of the [beagle-refine-impute.wdl](https://github.com/armartin/neurogap_downsampling/blob/master/beagle-refine-impute.wdl) and standard GATK workflows), gather concordance and sensitivity metrics (relative to full dataset) as in [variant_info.py](https://github.com/armartin/neurogap_downsampling/blob/master/variant_info.py)
