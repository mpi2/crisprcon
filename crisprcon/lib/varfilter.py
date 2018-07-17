
# -*- coding: utf-8 -*-


import csv
import os
import sys
from .infervariants import get_ref
from .vcfheader import vcf_header


try:
    csv.field_size_limit(sys.maxsize)
except:
    csv.field_size_limit(100000000)


def varFilter(dictBam, inVCF, outVCF, grna, ref_genome):
    lf = [c for c in dictBam if c['rpos'] == 0][0]
    rf = [c for c in dictBam if c['rpos'] == len(dictBam)-1][0]
    homArmDown = lf if lf['gs'] < rf['gs'] else rf
    homArmUps = rf if rf['ge'] > lf['ge'] else lf

    if homArmDown['strand'] != homArmUps['strand']:
        exit('ERROR: Exit variant calling: homologous regions have different orientation.')

    if homArmDown['chrom'] != homArmUps['chrom']:
        exit('ERROR: Exit variant calling: homologous regions are in different chromosome.')

    go = homArmDown['strand']
    gc = homArmDown['chrom']

    hs = [c['gs'] for c in dictBam if c['rpos'] == 0][0] if go == '+' else [c['ge'] for c in dictBam if c['rpos'] == 0][0]
    he = [c['ge'] for c in dictBam if c['rpos'] == len(dictBam)-1][0] if go == '+' else [c['gs'] for c in dictBam if c['rpos'] == len(dictBam)-1][0]

    with open(inVCF) as reader:
        # Parse file
        dictreader = _parse_vcf(reader)

        # Write out file
        _filter_vcf(dictreader, outVCF, dictBam, grna, hs, he, ref_genome)

    os.rename(outVCF + 'temp', outVCF)


def _parse_vcf(readable, fieldNames=None):

    prev, curr = 0, 0
    while True:
        line = readable.readline()
        if not line.startswith('#'):
            # lets start from prev # line, without the hash sign
            readable.seek(prev + 1)
            break
        else:
            prev = curr
            curr = readable.tell()

    # Determine dialect
    curr = readable.tell()
    # dialect = csv.Sniffer().sniff(readable.read(3000))
    dialect = 'excel-tab'
    readable.seek(curr)

    # Read file
    dictreader = csv.DictReader(readable, dialect=dialect, fieldnames=fieldNames)
    return dictreader


def _filter_vcf(dictreader, outVCF, dictBam, grna, hs, he, ref_genome):

    _vcf_fields = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')

    with open(outVCF + 'temp', mode='w') as writer:
        header = vcf_header(ref_genome)
        writer.write(header.h)
        # writer.write('#{}\n'.format('\t'.join(_vcf_fields)))

        output_vcf = []

        for line in dictreader:
            if _isSNP(line) or _isIndel(line):  # Only filter SNP and small indels
                # Filter by position: too close to homology arms or (optional) too far from a gRNA coordinate
                if abs(int(line['POS']) - hs) < 20 or abs(int(line['POS']) - he) < 20:
                    print ('Variant at ', line['POS'], 'filtered out. Too close to homology arms')
                    continue
                # # Filter by repetitive region
                # elif _isRepSnp(line, ref_genome) and not grna:
                #     print ('Variant at ', line['POS'], 'filtered out. SNP adjacent to repetitive region')
                #     continue
                # # Filter by repetitive region with optional gRNA
                # elif _isRepSnp(line, ref_genome) and grna:
                #     if _isgRna(line, grna):
                #         output_vcf.append([line['CHROM'], line['POS'], line['ID'], line['REF'],
                #                           line['ALT'], line['QUAL'], line['FILTER'], line['INFO']])
                #     else:
                #         print ('Variant at ', line['POS'], 'filtered out. SNP adjacent to repetitive region')
                #         continue
                # Filter by optional gRNA prox
                else:
                    if grna:
                        if _isgRna(line, grna):
                            output_vcf.append([line['CHROM'], line['POS'], line['ID'], line['REF'],
                                               line['ALT'], line['QUAL'], line['FILTER'], line['INFO']])
                        else:
                            print ('Variant at ', line['POS'], 'filtered out. Too far from gRNA position')
                            continue
                    else:
                        output_vcf.append([line['CHROM'], line['POS'], line['ID'], line['REF'],
                                           line['ALT'], line['QUAL'], line['FILTER'], line['INFO']])

            # Don't filter structural variants
            else:
                output_vcf.append([line['CHROM'], line['POS'], line['ID'], line['REF'],
                                   line['ALT'], line['QUAL'], line['FILTER'], line['INFO']])

        # Sort all results
        output_vcf.sort()
        output = "\n".join(["\t".join(map(str, vcf_row)) for vcf_row in output_vcf])

        # Write record
        writer.write(output + '\n')


def _isgRna(variant, grna):

    dist = 20

    grna_flat = [i for sublist in grna for i in sublist]

    return any([abs(int(variant['POS']) - c) <= dist for c in grna_flat])


def _isRepSnp(variant, ref_genome):
    # Filter variants (so far only SNPs) that follow or preceed a repetition of two or more bases of the same nucleotide
    if len(variant['REF']) == 1 and len(variant['ALT']) == 1:
        pseq = get_ref(ref_genome, variant['CHROM'], int(variant['POS'])-2, int(variant['POS'])+2)
        if pseq[0] == pseq[1] and pseq[0] == variant['ALT']:
            return (True)
        elif pseq[3] == pseq[4] and pseq[3] == variant['ALT']:
            return (True)
    else:
        return (False)


def _isSNP(variant):
    return (len(variant['REF']) == 1 and len(variant['ALT']) == 1)


def _isIndel(variant):
    # return ('INDEL' in variant['INFO'].split(';'))
    return (abs(len(variant['REF']) - len(variant['ALT'])) < 30 and abs(len(variant['REF']) - len(variant['ALT'])) > 0)
