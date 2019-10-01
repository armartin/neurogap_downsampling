import hail as hl
from pprint import pprint
import argparse


def main(args):
    full_vcf = hl.read_matrix_table(args.allreads_prefix + '.mt')
    
    # liftover chains
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)
    
    chips = hl.hadoop_open(args.chip_loci)
    chip_dict = {}
    for chip in chips:
        chip = chip.strip().split()
        chip_pos = hl.import_table(chip[1], filter='\[Controls\]', skip_blank_lines=True)
        chip_pos = chip_pos.filter(hl.array(list(map(str, range(1, 23))) + ['X', 'Y']).contains(chip_pos.chr))
        chip_pos = chip_pos.key_by(locus=hl.locus(chip_pos.chr, hl.int(chip_pos.pos)))

        #  liftover chip position info
        chip_pos = chip_pos.annotate(new_locus=hl.liftover(chip_pos.locus, 'GRCh38'))
        chip_pos = chip_pos.filter(hl.is_defined(chip_pos.new_locus))
        chip_pos = chip_pos.key_by(locus=chip_pos.new_locus)

        # filter full vcf to sites in genotype data
        geno_vcf = full_vcf.filter_rows(hl.is_defined(chip_pos[full_vcf.locus]))
        hl.export_vcf(geno_vcf, 'gs://neurogap/high_coverage/NeuroGap_30x_' + chip[0] + '.vcf.bgz')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--allreads_prefix', default='gs://neurogap/high_coverage/NeuroGap_30x_Pilot_Callset')
    parser.add_argument('--chip_loci', default='gs://neurogap/reference_data/chips/chip_loci.txt')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    main(args)
