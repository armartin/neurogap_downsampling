import hail as hl
from pprint import pprint
import argparse
import numpy as np


def load_files(file_prefix, overwrite, gencove, mt):
    """
    loads VCFs, run sample QC and variant QC, writes matrix table
    :param file_prefix:
    :param overwrite:
    :return:
    """
    if gencove:
        ngap_downsample = hl.read_matrix_table(file_prefix + '_grch38.mt')
    else:
        ngap_downsample = hl.import_vcf(file_prefix + '.vcf.gz', force_bgz=True, reference_genome='GRCh38',
                                        min_partitions=200)
        ngap_downsample = hl.split_multi_hts(ngap_downsample)
    ngap_downsample = ngap_downsample.filter_cols((ngap_downsample .s != 'NGE0018') & (ngap_downsample .s != 'NGE0130'))
    ngap_sample_qc = hl.sample_qc(ngap_downsample)
    ngap_sample_variant_qc = hl.variant_qc(ngap_sample_qc)
    ngap_sample_variant_qc.write(file_prefix + '.mt', overwrite=overwrite)


def geno_stats(ht_dict, geno_prefix):
    """
    computes numbers of hom ref, het, and hom var variants present in data
    :param ht_dict:
    :param geno_prefix:
    :return:
    """
    ht_samples = list(ht_dict.values())
    ht_samples = [
        ht.annotate(sample_qc=ht.sample_qc.drop('gq_stats', 'dp_stats')) if 'gq_stats' in list(ht.sample_qc) else ht for
        ht in ht_samples]
    ht_joined = ht_samples[0].union(*ht_samples[1:], unify=True)
    geno_stats = ht_joined.group_by(ht_joined.cov).aggregate(
        n_hom_ref_stats=hl.int(hl.agg.mean(ht_joined.sample_qc.n_hom_ref)),
        n_het_stats=hl.int(hl.agg.mean(ht_joined.sample_qc.n_het)),
        n_hom_var_stats=hl.int(hl.agg.mean(ht_joined.sample_qc.n_hom_var)))
    
    geno_stats.show(40)
    geno_stats.export(geno_prefix + 'geno_stats.tsv')


def compare_full(full_vcf, downsample_dict, downsample_prefix):
    """
    Counts numbers of variants in downsampled set present in full dataset
    :param full_vcf:
    :param downsample_dict:
    :param downsample_prefix:
    :return:
    """
    full_downsample = full_vcf.annotate(covs=hl.struct(
        **{cov: hl.is_defined(ht[full_vcf.key]) for cov, ht in downsample_dict.items()}
    ))
    
    count_comparisons = full_downsample.group_by(
        freq=hl.case().when(full_downsample.variant_qc.AC[1] == 1, 'singleton')
            .when(full_downsample.variant_qc.AC[1] <= 5, 'intermediate')
            .default('common')
    ).aggregate(
        **{'cov_' + x: hl.agg.fraction(full_downsample.covs[x]) for x in downsample_dict.keys()}
        )
    count_comparisons.export(downsample_prefix + 'present_in_full.tsv')


def concordance_tables(full_vcf, downsample_dict, output, overwrite):
    """
    runs concordance between full vcf and downsampled vcf
    :param full_vcf:
    :param downsample_dict:
    :param output:
    :param overwrite:
    :return:
    """
    global_conc, cols_conc, rows_conc = hl.concordance(full_vcf, downsample_dict)
    pprint(global_conc)
    cols_conc.write(output + 'samples.ht', overwrite=overwrite)
    rows_conc.write(output + 'variants.ht', overwrite=overwrite)


def concordance_frequency(full_vcf, concordance_table, output):
    full_variant_qc = full_vcf.rows()
    concordance_qc = full_variant_qc.annotate(concordance=concordance_table[full_variant_qc.key])
    freqs = list(np.linspace(0.5, 0, num=91)) ## note, this will need to be updated
    concordance_stats = concordance_qc.group_by(
        freq=hl.array(freqs).find(lambda x: concordance_qc.variant_qc.AF[1] >= x),
        snp=hl.is_snp(concordance_qc.alleles[0], concordance_qc.alleles[1])
    ).aggregate(
        n_variants=hl.agg.count(),
        unique_variants=hl.agg.array_agg(
            lambda row: hl.agg.array_agg(
                lambda element: hl.agg.count_where(element > 0),
                row),
            concordance_qc.concordance.concordance),
        geno_concordance=hl.agg.array_agg(
            lambda row: hl.agg.array_agg(
                lambda element: hl.agg.sum(element),
                row),
            concordance_qc.concordance.concordance)
    )

    concordance_stats = concordance_stats.annotate(total_concordant=concordance_stats.geno_concordance[3][3] +
                                                                    concordance_stats.geno_concordance[4][4],
                                                   total_discordant=concordance_stats.geno_concordance[2][3] +
                                                                    concordance_stats.geno_concordance[2][4] +
                                                                    concordance_stats.geno_concordance[3][2] +
                                                                    concordance_stats.geno_concordance[3][4] +
                                                                    concordance_stats.geno_concordance[4][2] +
                                                                    concordance_stats.geno_concordance[4][3]
                                                   )
    concordance_stats = concordance_stats.annotate(non_ref_concordance=concordance_stats.total_concordant/
                                                                       (concordance_stats.total_concordant +
                                                                        concordance_stats.total_discordant))
    concordance_stats.export(output + 'variants.tsv')


def main(args):
    if args.chip or args.gencove:
        file_prefix = args.downsample_prefix + '{cov}'
    else:
        file_prefix = args.downsample_prefix + '{cov}X'
    downsample_prefix = args.downsample_prefix
    if args.refined:
        file_prefix = file_prefix + '.refined'
        downsample_prefix = downsample_prefix + '.refined'
    elif args.imputed:
        file_prefix = file_prefix + '.imputed'
        downsample_prefix = downsample_prefix + '.imputed'

    if args.load_files:
        for cov in args.covs:
            load_files(file_prefix=file_prefix.format(cov=cov), overwrite=args.overwrite, gencove=args.gencove, mt=args.mt)
        # full depth mt
        load_files(file_prefix=args.allreads_prefix, overwrite=args.overwrite, gencove=False, mt=args.mt)

    if args.geno_stats:
        ht_dict = {cov: hl.read_matrix_table(file_prefix.format(cov=cov) + '.mt').cols().annotate(cov=cov) for cov in args.covs}
        ht_dict['all'] = hl.read_matrix_table(args.allreads_prefix + '.mt').cols().annotate(cov='all')
        geno_stats(ht_dict, downsample_prefix)

    if args.compare_full:
        downsample_dict = {cov: hl.read_matrix_table(file_prefix.format(cov=cov) + '.mt').rows() for cov in args.covs}
        full_vcf = hl.read_matrix_table(args.allreads_prefix + '.mt').rows()

        compare_full(full_vcf, downsample_dict, downsample_prefix)

    if args.concordance_tables:
        full_vcf = hl.read_matrix_table(args.allreads_prefix + '.mt')
        for cov in args.covs:
            downsample_vcf = hl.read_matrix_table(file_prefix.format(cov=cov) + '.mt')
            concordance_tables(full_vcf, downsample_vcf, file_prefix.format(cov=cov) + '_concordance_', args.overwrite)

    if args.concordance_frequency:
        full_vcf = hl.read_matrix_table(args.allreads_prefix + '.mt')
        for cov in args.covs:
            downsample_concordance = hl.read_table(file_prefix.format(cov=cov) + '_concordance_variants.ht')
            concordance_frequency(full_vcf, downsample_concordance, file_prefix.format(cov=cov) + '_concordance_')

    if args.pop_concordance:
        full_vcf = hl.read_matrix_table(args.allreads_prefix + '.mt')
        pop_info = hl.import_table(args.pop_info, key='high_cov_id')
        full_vcf = full_vcf.annotate_cols(**pop_info[full_vcf.s])
        sites = pop_info.aggregate(hl.agg.collect_as_set(pop_info.site))

        for site in sites:
            pop_full_vcf = full_vcf.filter_cols(full_vcf.site == site)
            for cov in args.covs:
                downsample_vcf = hl.read_matrix_table(file_prefix.format(cov=cov) + '.mt')
                downsample_vcf = downsample_vcf.annotate_cols(**pop_info[downsample_vcf.s])
                pop_downsample_vcf = downsample_vcf.filter_cols(downsample_vcf.site == site)
                concordance_tables(pop_full_vcf, pop_downsample_vcf, file_prefix.format(cov=cov) + '_' + site + '_concordance_',
                                   args.overwrite)
                pop_downsample_concordance = hl.read_table(file_prefix.format(cov=cov) + '_' + site + '_concordance_variants.ht')
                concordance_frequency(pop_full_vcf, pop_downsample_concordance, file_prefix.format(cov=cov) + '_' + site + '_concordance_')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--covs', type=lambda x: x.split(','), default=['0.5', '1.0', '2.0', '4.0', '6.0', '10.0', '20.0'], help="A comma separated list IDs")
    parser.add_argument('--downsample_prefix', default='gs://neurogap/high_coverage/neurogap_')
    parser.add_argument('--gencove', action='store_true')
    parser.add_argument('--mt', action='store_true')
    parser.add_argument('--refined', action='store_true')
    parser.add_argument('--imputed', action='store_true')
    parser.add_argument('--chip', action='store_true')
    parser.add_argument('--allreads_prefix', default='gs://neurogap/high_coverage/NeuroGap_30x_Pilot_Callset')
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--load_files', action='store_true')
    parser.add_argument('--geno_stats', action='store_true')
    parser.add_argument('--compare_full', action='store_true')
    parser.add_argument('--concordance_tables', action='store_true')
    parser.add_argument('--concordance_frequency', action='store_true')
    parser.add_argument('--pop_info', default='gs://neurogap/high_coverage/sample_id_map.txt')
    parser.add_argument('--pop_concordance', action='store_true')
    args = parser.parse_args()
    main(args)

