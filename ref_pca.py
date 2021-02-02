import hail as hl
import argparse
from gnomad.sample_qc.ancestry import *
import re


def intersect_datasets(ref_panel: str, ref_panel_qc: str, my_data: str, checkpoint: str, overwrite: bool = False):
    """
    Intersect NeuroGAP WGS with reference panel of HGDP + 1kG
    :param ref_panel: Path to matrix table of reference panel
    :param ref_panel_qc: Path to hail table of QC info to be annotated onto mt
    :param my_data: Path to matrix table of data to project
    :param checkpoint: Path to temporary intermediate mt
    :param overwrite: if True, overwrites existing data
    :return:
    """

    mt = hl.read_matrix_table(ref_panel)
    mt_qc = hl.read_table(ref_panel_qc)
    mt = mt.annotate_cols(**mt_qc[mt.s])
    all_sample_filters = set(mt_qc['sample_filters'])
    # bad_sample_filters were removing whole populations (mostly AFR and OCE) that passed all other QC, bad filters
    bad_sample_filters = {re.sub('fail_', '', x) for x in all_sample_filters if x.startswith('fail_')}
    # this filters to only variants that passed all gnomad QC or only failed filters in bad_sample_filters
    mt_filt = mt.filter_cols(( mt['sample_filters']['qc_metrics_filters'].difference(bad_sample_filters).length() == 0))

    # Example:
    # sample1: {} (keep)
    # sample2: {bad_sample_filter1, bad_sample_filter2} (keep)
    # sample3: {good_sample_filter1, bad_sample_filter1} (remove)

    # starting point: 4151
    # gnomAD QC: 3840
    # rescuing fails that shouldn't have by pop: 4017

    my_data = hl.read_matrix_table(my_data)

    # intersect hgdp + 1kG + my_data variants
    mt_intersect_snps = mt_filt.filter_rows(hl.is_defined(my_data.rows()[mt_filt.row_key]))
    mt_intersect_snps2 = my_data.filter_rows(hl.is_defined(mt_filt.rows()[my_data.row_key]))
    mt_intersect_snps = mt_intersect_snps.select_entries(*list(mt_intersect_snps2.entry))
    mt_intersect = mt_intersect_snps.select_cols().union_cols(mt_intersect_snps2.select_cols())

    mt_intersect.checkpoint(checkpoint, overwrite = overwrite, _read_if_exists = not overwrite)


def ld_prune_filter(mt: hl.MatrixTable, mt_ld: str, overwrite: bool = False):
    """
    Runs variant QC and filters out rare variants, those with missingness, and LD prunes to independent variants
    :param mt: Matrix table to run variant QC on and filter variants from
    :param mt_ld: Path to write intermediate filtered mt
    :param overwrite: if True, overwrites existing data
    :return:
    """
    mt.describe()
    mt = hl.variant_qc(mt)
    # mt_filt = mt.filter_rows((mt.variant_qc.AF[0] > 0.01) & (mt.variant_qc.AF[0] < 0.99))
    mt_filt = mt.filter_rows((mt.variant_qc.AF[0] > 0.05) & (mt.variant_qc.AF[0] < 0.95) &
                             (mt.variant_qc.call_rate > 0.999))

    # pruned = hl.ld_prune(mt_filt.GT, r2=0.2, bp_window_size=500000)
    pruned = hl.ld_prune(mt_filt.GT, r2=0.1, bp_window_size=500000)
    mt_filt = mt_filt.filter_rows(hl.is_defined(pruned[mt_filt.row_key]))
    mt_filt.write(mt_ld, overwrite)


def run_pc_relate(mt: hl.MatrixTable, pca_prefix: str, overwrite: bool = False):
    """
    Runs PC-relate to identify relatives in a matrix table
    :param mt: Matrix table to run PC-relate on
    :param pca_prefix: Prefix to path to output relatedness information
    :param overwrite: if True, overwrites existing data
    :return:
    """
    relatedness_ht = hl.pc_relate(mt.GT, min_individual_maf=0.05, min_kinship=0.05, statistics='kin',
                                  k=20).key_by()
    relatedness_ht.write(pca_prefix + 'relatedness.ht', args.overwrite)
    relatedness_ht = hl.read_table(pca_prefix + 'relatedness.ht')

    # identify individuals in pairs to remove
    related_samples_to_remove = hl.maximal_independent_set(relatedness_ht.i, relatedness_ht.j, False)
    mt_unrel = mt.filter_cols(hl.is_defined(related_samples_to_remove[mt.col_key]), keep=False)
    mt_rel = mt.filter_cols(hl.is_defined(related_samples_to_remove[mt.col_key]), keep=True)

    mt_unrel.write(pca_prefix + 'unrel.mt', args.overwrite)
    mt_rel.write(pca_prefix + 'rel.mt', args.overwrite)


def run_pca(mt: hl.MatrixTable, out_prefix: str, overwrite: bool = False):
    """
    Runs PCA on a dataset
    :param mt: dataset to run PCA on
    :param out_prefix: directory and filename prefix for where to put PCA output
    :return:
    """

    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt.GT, k=20, compute_loadings=True)
    pca_mt = mt.annotate_rows(pca_af=hl.agg.mean(mt.GT.n_alt_alleles()) / 2)
    pca_loadings = pca_loadings.annotate(pca_af=pca_mt.rows()[pca_loadings.key].pca_af)

    pca_scores.write(out_prefix + 'scores.ht', overwrite)
    pca_scores = hl.read_table(out_prefix + 'scores.ht')
    pca_scores = pca_scores.transmute(**{f'PC{i}': pca_scores.scores[i - 1] for i in range(1, 21)})
    pca_scores.export(out_prefix + 'scores.txt.bgz')  # individual-level PCs

    pca_loadings.write(out_prefix + 'loadings.ht', overwrite)  # PCA loadings


def project_individuals(pca_loadings, project_mt, out_prefix: str, overwrite: bool = False):
    """
    Project samples into predefined PCA space
    :param pca_loadings: existing PCA space
    :param project_mt: matrix table of data to project
    :param project_prefix: directory and filename prefix for where to put PCA projection output
    :return:
    """
    ht_projections = pc_project(project_mt, pca_loadings)
    ht_projections = ht_projections.transmute(**{f'PC{i}': ht_projections.scores[i - 1] for i in range(1, 21)})
    ht_projections.export(out_prefix + 'projected_scores.txt.bgz')
    return ht_projections


def main(args):
    if args.intersect_datasets:
        intersect_datasets(args.ref_panel, args.ref_panel_qc, args.my_data, args.checkpoint, args.overwrite)

    if args.run_ld_prune_filter:
        mt = hl.read_matrix_table(args.checkpoint)
        ld_prune_filter(mt, args.mt_ld, args.overwrite)

    if args.run_pc_relate:
        mt = hl.read_matrix_table(args.mt_ld)
        num_snps = mt.count_rows()
        run_pc_relate(mt, args.pca_prefix, args.overwrite)

    if args.run_pca:
        unrel = hl.read_matrix_table(args.pca_prefix + 'unrel.mt')
        run_pca(unrel, args.pca_prefix, args.overwrite)

    if args.project_relateds:
        loadings = hl.read_table(args.pca_prefix + 'loadings.ht')
        rel = hl.read_matrix_table(args.pca_prefix + 'rel.mt')
        project_individuals(loadings, rel, args.pca_prefix, args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_panel', default='gs://hgdp_tgp/output/tgp_hgdp.mt')
    parser.add_argument('--ref_panel_qc', default='gs://hgdp_tgp/output/gnomad_v3.1_sample_qc_metadata_hgdp_tgp_subset.ht')
    parser.add_argument('--my_data', default='gs://neurogap/high_coverage/NeuroGap_30x_Pilot_Callset.mt')

    parser.add_argument('--intersect_datasets', action='store_true')
    parser.add_argument('--checkpoint', default='gs://neurogap/high_coverage/temp.mt')

    parser.add_argument('--run_ld_prune_filter', action='store_true')
    parser.add_argument('--mt_ld', default='gs://neurogap/high_coverage/hgdp_tgp_neurogap_prune_filt.mt')

    parser.add_argument('--run_pc_relate', action='store_true')
    parser.add_argument('--run_pca', action='store_true')
    parser.add_argument('--project_relateds', action='store_true')
    parser.add_argument('--pca_prefix', default='gs://neurogap/high_coverage/hgdp_tgp_neurogap_pca_')

    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    main(args)
