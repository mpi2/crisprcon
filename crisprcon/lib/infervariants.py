
import subprocess
from .cigarParser import parse_cigar, parseMD
from sys import exit
from .alignment import flatten
from .ssr import _isSSR, _getSsrType, _ssrReadPos, TARGET_SEQ, mergeINS
from .vcfheader import vcf_header
import math


def main(dictBam, out_vcf, ref_genome, in_fa, temp, prefix):

    minSV = 1

    with open(out_vcf, mode='w') as writer:
        out = []
        header = vcf_header(ref_genome)
        writer.write(header.h)

        if _isSSR(in_fa):
            with open(in_fa, 'r') as fasta_ssr:
                for line in fasta_ssr:
                    if line.startswith('>') or line.startswith('@'):
                        continue
                    else:
                        seq = line
                        fasta_ssr_line = line
                        break
            out.append(list(flatten(inferSSR(dictBam=dictBam, ref_genome=ref_genome, in_fa=in_fa, ssr_type=_getSsrType(fasta_ssr_line)[0], ssr_pos=_ssrReadPos(fasta=fasta_ssr_line, dictBam=dictBam)))))
            in_fa = temp + '{}.noSSR.fasta'.format(prefix)

        with open(in_fa, 'r') as fasta:
            while True:
                header = fasta.readline()
                if not header.rstrip():
                    break
                seq_line = fasta.readline()
                if header.startswith('>'):
                    in_fa = '{}{}'.format(header, seq_line).rstrip()
                elif header.startswith('@'):
                    qual_header = fasta.readline()
                    qual_line = fasta.readline()
                    in_fa = '{}{}{}{}'.format(header, seq_line, qual_header, qual_line).rstrip()
                else:
                    exit('ERROR: Header {} seems malformed. Not fasta or fastq'.format(header))
        try:
            seq
        except:
            seq = seq_line

        out.append(list(flatten(inferVariants(dictBam=dictBam, ref_genome=ref_genome, fasta=in_fa, minSV=minSV))))
        out.append(list(flatten(inferIndels(dictBam=dictBam, ref_genome=ref_genome))))
        out.append(list(flatten(inferSNP(dictBam=dictBam))))
        out = list(flatten(out))
        if len(out) == 0:
            exit('ERROR: Failing to detect variants. No variants found.')

        out = mergeINS(out, seq, dictBam)

        writer.write(''.join(flatten(out)))


def inferVariants(dictBam, ref_genome, minSV, fasta):

    output = []

    if len(dictBam) == 1:
        return (output)

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

    i = 0
    while i < len(dictBam)-1:

        read = [c for c in dictBam if c['rpos'] == i][0]
        if len(dictBam) == 2:
            nread = [c for c in dictBam if c['rpos'] == i+1][0]
            output.append(inferSimpleSV(read=read, nread=nread, ncread=nread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=False, dictBam=dictBam))
            i += 1
        elif len(dictBam) > 2:
            nread = [c for c in dictBam if c['rpos'] == i+1][0]
            ncread = nread
            try:
                # ncread = [c for c in dictBam if c['gs'] >= read['ge']] if go == '+' else [c for c in dictBam if c['ge'] <= read['gs']]
                # ncread = min(ncread, key=lambda k: (k['gs'])) if go == '+' else max(ncread, key=lambda k: (k['ge']))
                ncread = [c for c in dictBam if c['ge'] > read['ge']]
                ncread = min(ncread, key=lambda k: (k['gs']))

            except:
                i += 1
                continue

            if nread == ncread:
                output.append(inferSimpleSV(read=read, nread=nread, ncread=ncread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=False, dictBam=dictBam))
            elif nread['ge'] <= homArmDown['gs'] or nread['gs'] >= homArmUps['ge'] or nread['chrom'] != gc:
                output.append(inferDUP(read=read, nread=nread, ncread=ncread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=False))
                output.append(inferSimpleSV(read=read, nread=ncread, ncread=ncread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=True, dictBam=dictBam))
                i += 1
            # elif (go == '+' and nread['gs'] >= read['ge'] and nread['ge'] <= ncread['gs']) or (go == '-' and nread['ge'] <= read['gs'] and nread['gs'] >= ncread['ge']):
            elif (nread['gs'] >= read['ge'] and nread['ge'] <= ncread['gs']):
                output.append(inferSimpleSV(read=read, nread=nread, ncread=ncread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=False, dictBam=dictBam))
            elif (nread['gs'] >= read['gs'] and nread['ge'] <= read['ge']) or (nread['gs'] >= ncread['gs'] and nread['ge'] <= ncread['ge']):
                output.append(inferSimpleSV(read=read, nread=nread, ncread=ncread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=False, dictBam=dictBam))
                output.append(inferSimpleSV(read=read, nread=ncread, ncread=ncread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=True, dictBam=dictBam))
                i += 1
            # elif go == '+' and (nread['ge'] <= read['gs'] or nread['gs'] >= ncread['ge']):
            #     output.append(inferDUP(read=read, nread=nread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=False))
            # elif go == '-' and (nread['gs'] >= read['ge'] or nread['ge'] <= ncread['gs']):
            #     output.append(inferDUP(read=read, nread=nread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=False))
            elif (nread['ge'] <= read['gs'] or nread['gs'] >= ncread['ge']):
                output.append(inferDUP(read=read, nread=nread, ncread=ncread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=False))
            i += 1
            continue

    return (output)


def inferIndels(dictBam, ref_genome):
    indels = []
    for i in range(len(dictBam)):
        read = dictBam[i]
        cigarParsed = parse_cigar(read['cigar'])
        seq = list(read['seq'])
        read_pos = 0
        pos = int(dictBam[i]['gs'])
        del_count = 0
        for event in cigarParsed:
            if event[0] == 'M' or event[0] == '=':
                pos += event[1]
                read_pos += event[1]
            elif event[0] == 'D':
                del_count += 1
                svtype = 'DEL'
                svlen = -event[1]
                svend = (pos-1) + event[1]
                ref = get_ref(ref_genome, read['chrom'], (pos-1), (pos-1)+abs(svlen))
                alt = ref[0]
                indels.append('{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};SVLEN={};END={}\n'.format(read['chrom'], pos-1, ref, alt, svtype, svlen, svend))
                pos += event[1]
            elif event[0] == 'I':
                svtype = 'INS'
                svlen = event[1]
                svend = pos - 1
                ref = ''.join(seq[read_pos-1])
                alt = ''.join(seq[(read_pos-1):(read_pos + event[1])])
                ins_phred = ins_qual(read['qual'][(read_pos):(read_pos + event[1])]) if read['qual'] != '*' else '.'
                indels.append('{}\t{}\t.\t{}\t{}\t{}\t.\tSVTYPE={};SVLEN={};END={}\n'.format(read['chrom'], pos-1, ref, alt, ins_phred, svtype, svlen, svend))
                read_pos += event[1]
            elif event[0] == 'S' or event[0] == 'H':
                continue
    return (indels)


def inferSNP(dictBam):
    snps = []
    for read in dictBam:
        mm = []
        read_snps = []
        cigarParsed = parse_cigar(read['cigar'])
        mdParsed = parseMD(read['tags'].split(':')[-1])
        g = 0
        r = 0
        for m in mdParsed:
            if m == 'D':
                g += 1
                continue
            elif m == 'M':
                r += 1
                g += 1
                continue
            else:
                mm.append([r, g, m])
                r += 1
                g += 1
        if any('I' in c for c in cigarParsed):
            pos = 0
            for event in cigarParsed:
                if event[0] == 'S' or event[0] == 'H' or event[0] == 'D':
                    continue
                if event[0] == 'M':
                    pos += event[1]
                if event[0] == 'I':
                    for i in range(len(mm)):
                        if mm[i][0] > pos:
                            mm[i][0] += event[1]
                    pos += event[1]
        for m in mm:
            chrom = read['chrom']
            st = read['gs'] + m[1]
            ref = m[2]
            alt = read['seq'][m[0]]
            qual = get_phred(read['qual'])[m[0]] if read['qual'] != '*' else '.'
            if alt not in ['A', 'C', 'G', 'T']:
                continue
            read_snps.append('{}\t{}\t.\t{}\t{}\t{}\t.\t.\n'.format(chrom, st, ref, alt, qual))
        snps.append(read_snps)
    return (snps)


def inferSimpleSV(read, nread, ncread, go, gc, ref_genome, minSV, fasta, nested, dictBam):
    sv = []

    header = fasta.split('\n')[0]
    seq_line = fasta.split('\n')[1]
    qual_header = ''
    qual_line = ''
    if header.startswith('@'):
        qual_header = fasta.split('\n')[2]
        qual_line = fasta.split('\n')[3]

    if nread['gs'] - read['ge'] > minSV:
        svstart = read['ge'] - 1
        svend = nread['gs'] - 1
        svtype = 'DEL'
        svlen = svstart - svend
        ref = get_ref(ref_genome, read['chrom'], (svstart), int(svend))
        alt = ref[0]
        sv.append('{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};SVLEN={};END={}\n'.format(read['chrom'], svstart, ref, alt, svtype, svlen, svend))
    if (max(nread['rs'], read['rs']) - min(read['re'], nread['re'])) - [nread['gs'] - read['ge'] if nread['gs'] - read['ge'] > 0 else 0][0] > minSV and nested is False:
        svstart = read['ge'] - 1
        svend = svstart
        svtype = 'INS'
        ref = get_ref(ref_genome, read['chrom'], (svstart), int(svend))
        alt = ref + seq_line[min(read['re'], nread['re'])-1:max(nread['rs'], read['rs'])-1]
        qual = ins_qual(qual_line[min(read['re'], nread['re'])-1:max(nread['rs'], read['rs'])-1]) if read['qual'] != '*' else '.'
        sv.append('{}\t{}\t.\t{}\t{}\t{}\t.\tSVTYPE={};END={}\n'.format(read['chrom'], svstart, ref, alt, qual, svtype, svend))

    elif min(read['re'], nread['re']) < max(nread['rs'], read['rs']) and nested is False:
        svstart = read['ge'] - 1
        svend = svstart
        svtype = 'INS'
        ref = get_ref(ref_genome, read['chrom'], (svstart), int(svend))
        alt = ref + seq_line[min(read['re'], nread['re'])-1:max(nread['rs'], read['rs'])-1]
        qual = get_phred(qual_line[min(read['re'], nread['re'])-1:max(nread['rs'], read['rs'])-1]) if read['qual'] != '*' else '.'
        sv.append('{}\t{}\t.\t{}\t{}\t{}\t.\tSVTYPE={};END={}\n'.format(read['chrom'], svstart, ref, alt, qual, svtype, svend))

    if min(nread['re'], ncread['re']) < max(nread['rs'], ncread['rs']) and nested is False and nread != ncread and ncread['rpos'] == nread['rpos'] + 1:
        svstart = read['ge'] - 1
        svend = svstart
        svtype = 'INS'
        ref = get_ref(ref_genome, read['chrom'], (svstart), int(svend))
        alt = ref + seq_line[min(nread['re'], ncread['re'])-1:max(ncread['rs'], nread['rs'])-1]
        qual = get_phred(qual_line[min(read['re'], nread['re'])-1:max(nread['rs'], read['rs'])-1]) if read['qual'] != '*' else '.'
        sv.append('{}\t{}\t.\t{}\t{}\t{}\t.\tSVTYPE={};END={}\n'.format(read['chrom'], svstart, ref, alt, qual, svtype, svend))

    if min(read['ge'], nread['ge']) > max(nread['gs'], read['gs']) + minSV:
        if read['ge'] == nread['gs']:
            # TODO: Maybe we can allow for more than 1 tandem duplication
            cs = max(nread['gs'], read['gs'])
            ce = min(read['ge'], nread['ge'])
            dupSeq = get_ref(ref_genome, read['chrom'], cs, ce - 1)
            svstart = max(read['ge'], nread['ge']) - 1
            svend = svstart
            svlen = len(dupSeq)
            svtype = 'DUP:TANDEM' if nread['strand'] == go else 'INVDUP:TANDEM'
            ref = get_ref(ref_genome, read['chrom'], svstart, svend)
            alt = ref + dupSeq
            qual = nread['mq']
            coords = '{}:{}-{}'.format(nread['chrom'], cs, ce - 1)
            sv.append('{}\t{}\t.\t{}\t{}\t{}\t.\tSVTYPE={};SVLEN={};END={};DUPCOORDS={}\n'.format(read['chrom'], svstart, ref, alt, qual, svtype, svlen, svend, coords))
            # return (sv)

        else:
            cs = max(nread['gs'], read['gs'])
            ce = min(read['ge'], nread['ge'])
            svstart = max(read['ge'], nread['ge']) - 1
            svend = svstart
            dupSeq = get_ref(ref_genome, read['chrom'], cs, ce - 1)

            svlen = len(dupSeq)
            svtype = 'DUP' if nread['strand'] == go else 'INVDUP'
            ref = get_ref(ref_genome, read['chrom'], int(svstart), int(svend))
            alt = ref + dupSeq
            qual = nread['mq']
            coords = '{}:{}-{}'.format(nread['chrom'], cs, ce - 1)
            sv.append('{}\t{}\t.\t{}\t{}\t{}\t.\tSVTYPE={};SVLEN={};END={};DUPCOORDS={}\n'.format(read['chrom'], svstart, ref, alt, qual, svtype, svlen, svend, coords))
            # return(sv)

    if nread['strand'] != go:
        svstart = nread['gs']
        svend = nread['ge'] - 1
        svtype = 'INV'
        svlen = nread['ge'] - nread['gs']
        ref = get_ref(ref_genome, read['chrom'], int(svstart), int(svend))
        alt = ref[::-1]
        sv.append('{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};SVLEN={};END={}\n'.format(read['chrom'], svstart, ref, alt, svtype, svlen, svend))

    return (sv)


def inferDUP(read, nread, ncread, go, gc, ref_genome, minSV, fasta, nested):
    sv = []
    svstart = read['ge'] - 1
    svend = svstart
    if (max(nread['rs'], read['rs']) - min(read['re'], nread['re'])) <= minSV and (nread['gs'] >= read['gs'] and nread['ge'] <= read['ge']):
        svtype = 'DUP:TANDEM' if nread['strand'] == go else 'INVDUP:TANDEM'
    else:
        svtype = 'DUP' if nread['strand'] == go else 'INVDUP'
    ref = get_ref(ref_genome, read['chrom'], (svstart), int(svend))
    alt = ref + nread['seq']
    qual = nread['mq']
    coords = '{}:{}-{}'.format(nread['chrom'], nread['gs'], nread['ge'] - 1)
    sv.append('{}\t{}\t.\t{}\t{}\t{}\t.\tSVTYPE={};END={};DUPCOORDS={}\n'.format(read['chrom'], svstart, ref, alt, qual, svtype, svend, coords))

    if min(read['re'], nread['re']) < max(nread['rs'], read['rs']) and nested is False:
        svstart = read['ge'] - 1
        svend = svstart
        svtype = 'INS'
        ref = get_ref(ref_genome, read['chrom'], (svstart), int(svend))
        alt = ref + fasta[min(read['re'], nread['re'])-1:max(nread['rs'], read['rs'])-1]
        sv.append('{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};END={}\n'.format(read['chrom'], svstart, ref, alt, svtype, svend))

    if min(nread['re'], ncread['re']) < max(nread['rs'], ncread['rs']) and nested is False and nread != ncread:
        # svstart = read['ge'] - 1 if go == '+' else nread['ge'] - 1
        svstart = ncread['gs'] - 1
        svend = svstart
        svtype = 'INS'
        ref = get_ref(ref_genome, nread['chrom'], (svstart), int(svend))
        alt = ref + fasta[min(nread['re'], ncread['re'])-1:max(ncread['rs'], nread['rs'])-1]
        sv.append('{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};END={}\n'.format(nread['chrom'], svstart, ref, alt, svtype, svend))

    return(sv)


def get_ref(fasta, chrom, start, end):
    p1 = subprocess.Popen(['samtools', 'faidx', fasta, '{}:{}-{}'.format(chrom, start, end)], stdout=subprocess.PIPE)
    refSeq = p1.communicate()[0].decode('utf-8')
    return(''.join(refSeq.split('\n')[1:]))


def inferSSR(ssr_pos, dictBam, ssr_type, in_fa, ref_genome):
    ssr_ins = []
    for pos in ssr_pos:
        chrom = pos.split(':')[0]
        svstart = int(pos.split(':')[1].split('-')[0]) - 1
        svend = svstart
        svtype = 'INS:LOXP' if ssr_type == 'LOXP' or ssr_type == 'LOXPRC' else 'INS:FRT'
        ref = get_ref(ref_genome, chrom, (svstart), int(svend))
        alt = ref + TARGET_SEQ[ssr_type]

        ssr_ins.append('{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};END={}\n'.format(chrom, svstart, ref, alt, svtype, svend))
    return(ssr_ins)


# def get_phred(qual_str):
#     # There might still be rare situations where this will fail. There is no right way to cover all cases
#     vals = [ord(c) for c in qual_str]  # Get the values for each symbol
#     val_range = (min(vals), max(vals))
#     if val_range[0] < 64:  # If the min value is smaller than 64, it must be Phred+33
#         return([x-33 for x in vals])
#     elif val_range[0] >= 64:  # If the min value is higher than 64, it must be Phred+64
#         return([x-64 for x in vals])

def get_phred(qual_str):
    # There might still be rare situations where this will fail. There is no right way to cover all cases
    vals = [ord(c) for c in qual_str]  # Get the values for each symbol
    return([x-33 for x in vals])


def calc_phred(prob):
    return (round(-10*math.log10(prob), 2))


def ins_prob(phred_list):
    probs = [1 - (10**(-x/10)) for x in phred_list]
    fragment_prob = 1
    for x in probs:
        fragment_prob *= x
    return(1 - fragment_prob)


def ins_qual(qual_str):
    phred_list = get_phred(qual_str)
    ins_probs = ins_prob(phred_list)
    return(calc_phred(ins_probs))
