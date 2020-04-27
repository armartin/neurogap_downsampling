import hail as hl
import argparse


def main(args):
    covs = ['0.5', '1', '2', '4', '6']

    # liftover gencove data
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)

    for cov in ['1']:
        all_vcf = hl.hadoop_open('gs://neurogap/high_coverage/gencove/fastq_keep.txt')
        header = all_vcf.readline().split()
        header_index = dict(zip(header, list(range(len(header)))))
        comb = None

        if args.merge_vcf:
            for line in all_vcf:
                # read and merge
                line = line.split()
                if line[header_index['depth']] == cov and line[header_index['gencove_qc']] == 'PASS':
                    vcf = hl.import_vcf(line[header_index['path']], force_bgz=True, reference_genome='GRCh37', min_partitions=100)
                    if comb is None:
                        comb = vcf
                    else:
                        comb = comb.union_cols(vcf)

            # write out mt
            #comb = comb.naive_coalesce(5000)
            comb.write('gs://neurogap/high_coverage/gencove/merge_' + cov + '_hg19.mt', overwrite=args.overwrite)

        if args.liftover:
            comb = hl.read_matrix_table('gs://neurogap/high_coverage/gencove/merge_' + cov + '_hg19.mt')

            comb = comb.annotate_rows(new_locus=hl.liftover(comb.locus, 'GRCh38'))
            comb = comb.filter_rows(hl.is_defined(comb.new_locus))
            comb = comb.key_rows_by(locus=comb.new_locus, alleles=comb.alleles)

            # write out mt
            comb.write('gs://neurogap/high_coverage/gencove/merge_' + cov + '_grch38.mt', overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--merge_vcf', action='store_true')
    parser.add_argument('--liftover', action='store_true')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    main(args)
