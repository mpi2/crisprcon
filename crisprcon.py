#!/usr/bin/env python

from crisprcon.lib import crispFa
from crisprcon.lib import crispRNA
import argparse
import sys
import os


def main(args):
    if args.subparser_name == 'crispFa':
        chrom = args.coords.split(':')[0]
        ts = args.coords.split(':')[1].split('-')[0]
        te = args.coords.split(':')[1].split('-')[1]
        crispFa(input_fa=args.input_fa,
                dir=args.out_dir,
                temp_dir=args.temp_dir,
                output_vcf=args.out_dir + args.output_vcf,
                ref_genome=args.ref_genome,
                prefix=args.prefix,
                gff=args.gff,
                chrom=chrom,
                ts=ts,
                te=te,
                grna=args.grna)

    elif args.subparser_name == 'crispRNA':
        chrom = args.grna[0][1].split(':')[0]
        ts = int(args.grna[0][1].split(':')[1].split('-')[0]) - 1000
        te = int(args.grna[-1][1].split(':')[1].split('-')[1]) + 1000

        crispRNA(dir=args.out_dir,
                 temp_dir=args.temp_dir,
                 output_vcf=args.out_dir + args.output_vcf,
                 ref_genome=args.ref_genome,
                 prefix=args.prefix,
                 gff=args.gff,
                 grna=args.grna,
                 donor=args.donor)
        print ('Predicting consequence')
        crispFa(input_fa=args.temp_dir + args.prefix + '.crispr.fasta',
                dir=args.out_dir,
                temp_dir=args.temp_dir,
                output_vcf=args.out_dir + args.output_vcf,
                ref_genome=args.ref_genome,
                prefix=args.prefix,
                gff=args.gff,
                chrom=chrom,
                ts=ts,
                te=te,
                grna=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='crispCon', description='Variant analysis for CRISPR edited alleles.')
    subparsers = parser.add_subparsers(dest='subparser_name')

    # subparser for crisprCon
    parser_crispfa = subparsers.add_parser('crispFa', help='Variant caller for CRISPR edited alleles')
    parser_crispfa.add_argument('-i', '--input', dest='input_fa', help='Input FASTA or BAM file (format autodetected). Default: stdin', default=sys.stdin)
    parser_crispfa.add_argument('-o', '--output', dest='output_vcf', help='Output VCF file. Default: stdout', default=sys.stdout)
    parser_crispfa.add_argument('-f', '--fasta-ref', dest='ref_genome', help='Reference file in fasta format (must be .fai indexed)')
    parser_crispfa.add_argument('-r', '--region', dest='coords', help='Coordinates of CRISPR/Cas target region (ie: chr:start-end)')
    parser_crispfa.add_argument('-g', '--gff-annot', dest='gff', help='gff3 annotation file')

    parser_crispfa.add_argument('-d', '--dir', dest='out_dir', help='Optional: Directory for the output files. If not present, it will be created. Default: current directory', default=os.getcwd())
    parser_crispfa.add_argument('-t', '--temp', dest='temp_dir', help='Optional: Directory for temporary and intermediate files. Default: creates "temp" directory at the output directory')
    parser_crispfa.add_argument('-p', '--prefix', dest='prefix', help='Optional: Prefix for temporary files. Default: input_name')
    parser_crispfa.add_argument('-c', '--crispr-grna', dest='grna', nargs=2, type=int, action='append', help='Optional: gRNA pair coordinates used for CRISPR. Advised for validation of indels and potential false positives. Use different flags for more than one pair. Eg: -c 72206978 72207000 -c 72206982 72207004')

    # subparser for crispRNA
    parser_crisprna = subparsers.add_parser('crispRNA', help='gRNA consequence predictor. It takes the coordinates of gRNA and it outputs a VCF with the predicted consequences')
    parser_crisprna.add_argument('-o', '--output', dest='output_vcf', help='Output VCF file. Default: stdout', default=sys.stdout)
    parser_crisprna.add_argument('-f', '--fasta-ref', dest='ref_genome', help='Reference file in fasta format (must be .fai indexed)')
    parser_crisprna.add_argument('-g', '--gff-annot', dest='gff', help='gff3 annotation file')
    parser_crisprna.add_argument('-c', '--crispr-grna', dest='grna', nargs=2, action='append', help='gRNA pair coordinates and sequence used for CRISPR (sequence chr:start-end). Use different flags for more than one pair. Eg: -c GATTCTGGCATCATCTATGTGGG 1:72206978-72207000 -c GCTTCTGCCAATATCTATTTGGG 1:72206982-72207004')
    parser_crisprna.add_argument('-s', '--seq-donor', dest='donor', action='append', help='Sequence of the oligo donor')
    parser_crisprna.add_argument('-d', '--dir', dest='out_dir', help='Optional: Directory for the output files. If not present, it will be created. Default: current directory', default=os.getcwd())
    parser_crisprna.add_argument('-t', '--temp', dest='temp_dir', help='Optional: Directory for temporary and intermediate files. Default: creates "temp" directory at the output directory')
    parser_crisprna.add_argument('-p', '--prefix', dest='prefix', help='Optional: Prefix for temporary files. Default: input_name')

    args = parser.parse_args()

    if args.subparser_name is None:
        parser.print_help()
        exit(1)

    if args.subparser_name == 'crispFa':
        if type(args.input_fa) is not str:
            if sys.stdin.isatty():
                print ('\n**ERROR: input fasta or bam file required**\n\n')
                parser_crispfa.print_help()
                exit(1)
        elif not os.path.exists(args.input_fa):
            print ('\n**ERROR: input fasta or bam file {} not found**\n\n'.format(args.input_fa))
            parser_crispfa.print_help()
            exit(1)

    if args.ref_genome is None:
        print ('\n**ERROR: reference fasta file required**\n\n')
        parser_crispfa.print_help() if args.subparser_name == 'crispFa' else parser_crisprna.print_help()
        exit(1)
    elif not os.path.exists(args.ref_genome):
        print ('\n**ERROR: reference fasta file {} not found**\n\n'.format(args.ref_genome))
        parser_crispfa.print_help() if args.subparser_name == 'crispFa' else parser_crisprna.print_help()
        exit(1)

    args.out_dir = args.out_dir + '/' if args.out_dir[len(args.out_dir)-1] is not '/' else args.out_dir

    if args.temp_dir is None:
        if not os.path.exists(args.out_dir + 'temp'):
            os.makedirs(args.out_dir + 'temp')
            args.temp_dir = args.out_dir + 'temp'
        elif os.path.exists(args.out_dir + 'temp'):
            args.temp_dir = args.out_dir + 'temp'

    if args.subparser_name == 'crispFa':
        if args.prefix is None:
            args.prefix = args.input_fa.rsplit('.', 1)[0]

    if args.grna is None:
        if args.subparser_name == 'crispFa':
            args.grna = False
        if args.subparser_name == 'crispRNA':
            print ('\n**ERROR: no gRNA coordinates given**\n\n')
            parser_crisprna.print_help()
            exit(1)

    args.temp_dir = args.temp_dir + '/' if args.temp_dir[len(args.temp_dir)-1] is not '/' else args.temp_dir

    main(args)
