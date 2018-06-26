
import subprocess
from .cigarParser import parse_cigar
from sys import exit
from .alignment import flatten, reverseComplement
from .ssr import _isSSR, _getSsrType, _ssrReadPos, TARGET_SEQ, mergeINS
from .vcfheader import vcf_header


def main(dictBam, out_vcf, ref_genome, in_fa, temp, prefix):

    minSV = 15

    with open(out_vcf, mode='w') as writer:
        out = []
        writer.write(vcf_header)
        if _isSSR(in_fa):
            with open(in_fa, 'r') as fasta_ssr:
                for line in fasta_ssr:
                    if not line.startswith('>'):
                        seq = line
                        fasta_ssr_line = line
            out.append(list(flatten(inferSSR(dictBam=dictBam, ref_genome=ref_genome, in_fa=in_fa, ssr_type=_getSsrType(fasta_ssr_line)[0], ssr_pos=_ssrReadPos(fasta=fasta_ssr_line, dictBam=dictBam)))))
            in_fa = temp + '{}.noSSR.fasta'.format(prefix)

        with open(in_fa, 'r') as fasta:
            for line in fasta:
                if not line.startswith('>'):
                    try:
                        seq
                    except:
                        seq = line
                    fasta_line = line
        out.append(list(flatten(inferVariants(dictBam=dictBam, ref_genome=ref_genome, in_fa=in_fa, minSV=minSV, fasta=fasta_line))))
        out.append(list(flatten(inferIndels(dictBam=dictBam, ref_genome=ref_genome))))

        out = mergeINS(list(flatten(out)), seq, dictBam)

        writer.write(''.join(flatten(out)))


def inferVariants(dictBam, ref_genome, in_fa, minSV, fasta):

    output = []

    if len(dictBam) == 1:
        return (output)

    lf = [c for c in dictBam if c['rpos'] == 0][0]
    rf = [c for c in dictBam if c['rpos'] == len(dictBam)-1][0]
    homArmDown = lf if lf['gs'] < rf['gs'] else rf
    print ('Homology arm downstream')
    print (homArmDown)
    homArmUps = rf if rf['ge'] > lf['ge'] else lf
    print ('Homology arm upstream')
    print (homArmUps)
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
            output.append(inferSimpleSV(read=read, nread=nread, ncread=nread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=False))
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
            print ('The iteration number is')
            print (i)
            print ('The read is')
            print (read)
            print ('The nread is')
            print (nread)
            print ('The ncread is')
            print (ncread)

            if nread == ncread:
                print ('one')
                output.append(inferSimpleSV(read=read, nread=nread, ncread=ncread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=False))
            elif nread['ge'] <= homArmDown['gs'] or nread['gs'] >= homArmUps['ge'] or nread['chrom'] != gc:
                print ('two')
                output.append(inferDUP(read=read, nread=nread, ncread=ncread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=False))
                output.append(inferSimpleSV(read=read, nread=ncread, ncread=ncread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=True))
                i += 1
            # elif (go == '+' and nread['gs'] >= read['ge'] and nread['ge'] <= ncread['gs']) or (go == '-' and nread['ge'] <= read['gs'] and nread['gs'] >= ncread['ge']):
            elif (nread['gs'] >= read['ge'] and nread['ge'] <= ncread['gs']):
                print ('three')
                output.append(inferSimpleSV(read=read, nread=nread, ncread=ncread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=False))
            elif (nread['gs'] >= read['gs'] and nread['ge'] <= read['ge']) or (nread['gs'] >= ncread['gs'] and nread['ge'] <= ncread['ge']):
                print ('four')
                output.append(inferSimpleSV(read=read, nread=nread, ncread=ncread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=False))
                output.append(inferSimpleSV(read=read, nread=ncread, ncread=ncread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=True))
                i += 1
            # elif go == '+' and (nread['ge'] <= read['gs'] or nread['gs'] >= ncread['ge']):
            #     output.append(inferDUP(read=read, nread=nread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=False))
            # elif go == '-' and (nread['gs'] >= read['ge'] or nread['ge'] <= ncread['gs']):
            #     output.append(inferDUP(read=read, nread=nread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=False))
            elif (nread['ge'] <= read['gs'] or nread['gs'] >= ncread['ge']):
                print ('five')
                output.append(inferDUP(read=read, nread=nread, ncread=ncread, go=go, gc=gc, ref_genome=ref_genome, minSV=minSV, fasta=fasta, nested=False))
            i += 1
            print ('And the output so far is:')
            print (output)
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
                indels.append('{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};SVLEN={};END={}\n'.format(read['chrom'], pos-1, ref, alt, svtype, svlen, svend))
                read_pos += event[1]
            elif event[0] == 'S' or event[0] == 'H':
                continue
    return (indels)

# def inferSNP(dictBam):
#     snps = []
#     for i in range(len(dictBam)):


def inferSimpleSV(read, nread, ncread, go, gc, ref_genome, minSV, fasta, nested):
    sv = []
    if nread['gs'] - read['ge'] > minSV:
        svstart = read['ge'] - 1
        svend = nread['gs'] - 1
        svtype = 'DEL'
        svlen = svstart - svend
        ref = get_ref(ref_genome, read['chrom'], (svstart), int(svend))
        alt = ref[0]
        sv.append('{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};SVLEN={};END={}\n'.format(read['chrom'], svstart, ref, alt, svtype, svlen, svend))

    if (max(nread['rs'], read['rs']) - min(read['re'], nread['re'])) - [nread['gs'] - read['ge'] if nread['gs'] - read['ge'] > 0 else 0][0] > minSV:
        # svstart = read['ge'] - 1 if go == '+' else nread['ge'] - 1
        svstart = read['ge'] - 1
        svend = svstart
        svtype = 'INS'
        ref = get_ref(ref_genome, read['chrom'], (svstart), int(svend))
        alt = ref + fasta[min(read['re'], nread['re'])-1:max(nread['rs'], read['rs'])-1]
        sv.append('{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};END={}\n'.format(read['chrom'], svstart, ref, alt, svtype, svend))

    elif min(read['re'], nread['re']) < max(nread['rs'], read['rs']) and nested is False:
        # svstart = read['ge'] - 1 if go == '+' else nread['ge'] - 1
        svstart = read['ge'] - 1
        svend = svstart
        svtype = 'INS'
        ref = get_ref(ref_genome, read['chrom'], (svstart), int(svend))
        alt = ref + fasta[min(read['re'], nread['re'])-1:max(nread['rs'], read['rs'])-1]
        sv.append('{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};END={}\n'.format(read['chrom'], svstart, ref, alt, svtype, svend))

    if min(nread['re'], ncread['re']) < max(nread['rs'], ncread['rs']) and nested is False and nread != ncread:
        # svstart = read['ge'] - 1 if go == '+' else nread['ge'] - 1
        svstart = read['ge'] - 1
        svend = svstart
        svtype = 'INS'
        ref = get_ref(ref_genome, read['chrom'], (svstart), int(svend))
        alt = ref + fasta[min(nread['re'], ncread['re'])-1:max(ncread['rs'], nread['rs'])-1]
        sv.append('{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};END={}\n'.format(nread['chrom'], svstart, ref, alt, svtype, svend))


    # if (nread['gs'] - read['ge']) > minSV and (nread['rs'] - read['re']) < minSV:
    if min(read['ge'], nread['ge']) > max(nread['gs'], read['gs']) + minSV and (max(nread['rs'], read['rs']) - min(read['re'], nread['re'])) <= minSV:
        # svstart = nread['gs'] - 1 if go == '+' else read['gs'] - 1
        # svend = read['ge'] - 1 if go == '+' else nread['ge'] - 1
        cs = max(nread['gs'], read['gs'])
        ce = min(read['ge'], nread['ge'])
        # dupSeq = get_ref(ref_genome, int(read['chrom']), (int(nread['gs'])), int(read['ge']))
        dupSeq = get_ref(ref_genome, read['chrom'], cs - 1, ce - 1)
        svstart = max(read['ge'], nread['ge']) - 1
        svend = svstart
        # ndups = 2  # TODO: Allow for more than 1 tandem duplication
        # ndups = 0
        # nnread = read[i + 1 + ndups]
        # while read['ge'] == nnread['ge']:
        #     ndups += 1
        #     nnread = read[i + 1 + ndups]
        svlen = len(dupSeq)
        svtype = 'DUP:TANDEM' if nread['strand'] == go else 'INVDUP:TANDEM'
        ref = get_ref(ref_genome, read['chrom'], svstart, svend)
        alt = ref + dupSeq
        coords = '{}:{}-{}'.format(nread['chrom'], cs, ce - 1)
        sv.append('{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};SVLEN={};END={};DUPCOORDS={}\n'.format(read['chrom'], svstart, ref, alt, svtype, svlen, svend, coords))
        return (sv)

    elif min(read['ge'], nread['ge']) > max(nread['gs'], read['gs']) + minSV and (max(nread['rs'], read['rs']) - min(read['re'], nread['re'])) > minSV:
        cs = max(nread['gs'], read['gs'])
        ce = min(read['ge'], nread['ge'])
        svstart = max(read['ge'], nread['ge']) - 1
        svend = svstart
        dupSeq = get_ref(ref_genome, read['chrom'], cs - 1, ce - 1)

        svlen = len(dupSeq)
        svtype = 'DUP' if nread['strand'] == go else 'INVDUP'
        ref = get_ref(ref_genome, read['chrom'], int(svstart), int(svend))
        alt = ref + nread['seq']
        coords = '{}:{}-{}'.format(nread['chrom'], cs, ce - 1)
        sv.append('{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};SVLEN={};END={};DUPCOORDS={}\n'.format(read['chrom'], svstart, ref, alt, svtype, svlen, svend, coords))
        return(sv)

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
    coords = '{}:{}-{}'.format(nread['chrom'], nread['gs'], nread['ge'] - 1)
    sv.append('{}\t{}\t.\t{}\t{}\t.\t.\tSVTYPE={};END={};DUPCOORDS={}\n'.format(read['chrom'], svstart, ref, alt, svtype, svend, coords))

    if min(read['re'], nread['re']) < max(nread['rs'], read['rs']) and nested is False:
        # svstart = read['ge'] - 1 if go == '+' else nread['ge'] - 1
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
