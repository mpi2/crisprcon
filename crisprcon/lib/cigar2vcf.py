

import cigarParser
from itertools import islice


vcf_header = '''##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##reference=file:///nfs/production/mousegenomes/projects/REL-17/reference/GRCm38_68.fa
##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
##contig=<ID=8>
##contig=<ID=9>
##contig=<ID=10>
##contig=<ID=11>
##contig=<ID=12>
##contig=<ID=13>
##contig=<ID=14>
##contig=<ID=15>
##contig=<ID=16>
##contig=<ID=17>
##contig=<ID=18>
##contig=<ID=19>
##contig=<ID=X>
##contig=<ID=Y>
##contig=<ID=MT>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP:TANDEM,Description="Tandem duplication">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=CNV,Description="Copy number variable region">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
'''


def infer_variants(in_bam, out_vcf, ref_genome):
    with open(out_vcf, mode='w') as writer:
        writer.write(vcf_header)

        dictBam = cigarParser.parse_bam(in_bam)
        pos = dictBam[0]['start']
        read_pos = 0
        del_count = 0
        large_del = False
        for read in dictBam:
            cigarParsed = cigarParser.parse_cigar(read['cigar'])
            seq = list(read['seq'])
            if cigarParsed[0][0] == 'H' or cigarParsed[0][0] == 'S':
                if not large_del:
                    read_pos = cigarParsed[0][1]
            if cigarParsed[0][0] == 'H' or cigarParsed[0][0] == 'S' and large_del:
                pos += 1
                svend = read['start']
                svtype = 'DEL'
                svlen = pos - svend
                ref = get_ref(ref_genome, int(read['chrom']), (pos - 1), int(svend))
                alt = ''.join(seq[(read_pos - 1)])
                variant = '{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};SVLEN={};END={}\n'.format(read['chrom'], (pos - 1), ref, alt, svtype, svlen, svend)
                writer.write(variant)
                pos = read['start']
                read_pos = cigarParsed[0][1]
                large_del = False
            for event in cigarParsed:
                if event[0] == 'M' or event[0] == '=':
                    pos += event[1]
                    read_pos += event[1]
                elif event[0] == 'D':
                    del_count += 1
                    svtype = 'DEL'
                    svlen = -event[1]
                    svend = (pos + event[1])
                    alt = ''.join(seq[(read_pos - 1)])
                    deleted = cigarParser.get_deletion(read['tags'], del_count)
                    ref = ''.join([alt + deleted if deleted != 'N' else deleted])
                    variant = '{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};SVLEN={};END={}\n'.format(read['chrom'], pos, ref, alt, svtype, svlen, svend)
                    writer.write(variant)
                    pos += event[1]
                elif event[0] == 'I':
                    svtype = 'INS'
                    svlen = event[1]
                    svend = pos
                    ref = ''.join(seq[read_pos-1])
                    alt = ''.join(seq[(read_pos-1):(read_pos + event[1])])
                    variant = '{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};SVLEN={};END={}\n'.format(read['chrom'], pos, ref, alt, svtype, svlen, svend)
                    writer.write(variant)
                    # pos += 1
                elif event[0] == 'H' or event[0] == 'S':
                    large_del = True
                    continue


def get_ref(fasta, chr, start, end):
    with open(fasta + '.fai') as index:
        for line in index:
            if line.startswith('{}\t'.format(chr)):
                bytes = int(line.split('\t')[2])
                lines = round(bytes / int(line.split('\t')[4]))
                break
    with open(fasta) as f:
        copy = False
        refSeq = []
        for line in islice(f, round(lines-0.5), None):
            if copy and line[0].isalpha():
                refSeq += list(line.strip('\n'))
            elif line.startswith('>{} '.format(chr)):
                copy = True
            elif line.startswith('>'):
                break
        return (''.join(refSeq[(start - 1):end]))


def get_ref3(fasta, chr, start, end):
    with open(fasta + '.fai') as index:
        for line in index:
            if line.startswith('{}\t'.format(chr)):
                bytes = int(line.split('\t')[2])
                lines = round(bytes / int(line.split('\t')[4]))
                break
    with open(fasta) as f:
        copy = False
        refSeq = []
        for i in range(round(lines-0.5)):
            next(f)
        for line in f:
            if copy and line[0].isalpha():
                refSeq += list(line.strip('\n'))
            elif line.startswith('>{} '.format(chr)):
                copy = True
            elif line.startswith('>'):
                break
        return (''.join(refSeq[(start - 1):end]))
