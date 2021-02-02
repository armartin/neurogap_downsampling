import gzip
import argparse


def open_vcf(args, current_chrom):
    chr_vcf = gzip.open(args.vcf.split('.vcf.gz')[0] + '_' + str(current_chrom) + '.vcf.gz', 'wt')
    return chr_vcf


def main(args):
    vcf = gzip.open(args.vcf)
    header = []
    last_chrom = ''
    for line in vcf:
        current_line = line.strip().decode('utf-8')
        current_chrom = current_line.split()[0]
        if current_line.startswith('#'):
            header.append(current_line)
            continue
        if current_chrom != last_chrom:
            if last_chrom.startswith('chr'):
                chr_vcf.close()
            chr_vcf = open_vcf(args, current_chrom)
            for header_line in header:
                chr_vcf.write((header_line + '\n'))
        chr_vcf.write((current_line + '\n'))
        last_chrom = current_chrom
    chr_vcf.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf')
    args = parser.parse_args()
    main(args)
