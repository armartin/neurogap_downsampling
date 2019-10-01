import hail as hl
import argparse


def main(args):
    tgp_vcf = hl.import_vcf(args.input_vcf, force_bgz=True, header_file='gs://neurogap/reference_data/header', reference_genome='GRCh38', min_partitions=200)
    sample_info = hl.import_table(args.sample_info).key_by('sample')
    tgp_vcf = tgp_vcf.annotate_cols(**sample_info[tgp_vcf.s])
    afr_vcf = tgp_vcf.filter_cols((tgp_vcf.super_pop == 'AFR') & ~((tgp_vcf.pop == 'ASW') | (tgp_vcf.pop == 'ACB'))) #  NA18498 is missing
    sample_info.filter(hl.is_missing(tgp_vcf.cols()[sample_info['sample']])).show()
    hl.export_vcf(afr_vcf, args.output_filename)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_vcf', default='gs://neurogap/reference_data/ALL.chr1.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz')
    parser.add_argument('--sample_info', default='gs://neurogap/reference_data/integrated_call_samples_v3.20130502.ALL.panel')
    parser.add_argument('--output_filename', default='gs://neurogap/reference_data/AFR.chr1.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.bgz')
    args = parser.parse_args()
    main(args)
