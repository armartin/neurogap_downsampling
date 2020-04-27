import hail as hl
import argparse


gencove_covs = ['0.5', '1', '2', '4', '6']
beagle_covs = ['0.5', '1.0', '2.0', '4.0', '6.0']

for cov in range(len(gencove_covs)):
    gencove = hl.read_matrix_table('gs://neurogap/high_coverage/gencove/merge_' + gencove_covs[cov] + '_grch38.mt') # (85154681, 73)
    beagle = hl.read_matrix_table('gs://neurogap/high_coverage/neurogap_' + beagle_covs[cov] + 'X.imputed.mt') # (73044119, 91)

    # filter gencove to sites in beagle
    gencove_in_beagle = gencove.filter_rows(hl.is_defined(beagle.rows()[gencove.row_key])) # 68264356, 73
    beagle_in_gencove = beagle.filter_rows(hl.is_defined(gencove.rows()[beagle.row_key])) # (68264084, 91)

    ##
    gencove_in_beagle.write('gs://neurogap/high_coverage/gencove/merge_intersect_' + gencove_covs[cov] + '_grch38.mt', overwrite=True)
    beagle_in_gencove.write('gs://neurogap/high_coverage/neurogap_intersect_' + beagle_covs[cov] + 'X.imputed.mt', overwrite=True)
