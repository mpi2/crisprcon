
from .cigarParser import parse_cigar
# from .alignment import reverseComplement

TARGET_SEQ = {
            'LOXP': 'ATAACTTCGTATAGCATACATTATACGAAGTTAT',
            'LOXPRC': 'ATAACTTCGTATAATGTATGCTATACGAAGTTAT',
            'FRT': 'GAAGTTCCTATTCCGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTC'
}


def reverseComplement(pattern):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                  'Y': 'R', 'R': 'Y', 'W': 'W', 'S': 'S', 'K': 'M',
                  'M': 'K', 'D': 'H', 'V': 'B', 'H': 'D', 'B': 'V',
                  'X': 'X'}
    pattern = pattern.strip()
    return ''.join([complement[base] for base in pattern[::-1]])


def _ssr_analysis(fasta):
    # Identifies a LoxP or FRT (Site Specific Recombination) sequence in the input fasta (a line with the sequence) and:
        # 1) Locates the insertion position in the recombinate sequence
        # 2) Returns the fasta sequence without the inserted sequence

    recombinate = _getSsrType(fasta)

    if not recombinate:
        return (fasta)
    elif len(recombinate) > 1:
        print ('WARNING: {} found. No site specific recombination considered.'.format(' and '.join(recombinate)))
        return (fasta)
    else:
        recombinate = recombinate[0]

    if recombinate == 'LOXP' or recombinate == 'LOXPRC':
        return(_loxp_analysis(fasta, recombinate))
    if recombinate == 'FRT':
        return(_frt_analysis(fasta, recombinate))


def _loxp_analysis(fasta, recombinate):
    fasta_clean = fasta.replace(TARGET_SEQ[recombinate], '')
    return (fasta_clean)


def _frt_analysis(fasta, recombinate):
    fasta_clean = fasta.replace(TARGET_SEQ[recombinate], '')
    return (fasta_clean)


def _isSSR(in_fa):
    with open(in_fa, 'r') as fasta:
        for line in fasta:
            if not line.startswith('>'):
                if any(line.find(TARGET_SEQ[c]) >= 0 for c in TARGET_SEQ):
                    return True
        return False


def _ssrReadPos(fasta, dictBam):
    # Input: fasta line where there is a SSR Insertion
    #        dictBam with the reads
    # Output: list of read positions were the insertions occured
    go = [c['strand'] for c in dictBam if c['rpos'] == 0][0]

    try:
        ssr_seq = [TARGET_SEQ[c] for c in TARGET_SEQ if fasta.find(TARGET_SEQ[c]) >= 0][0]
    except:
        exit('ERROR: No recombinate sequence')
    if go == '-':
        fasta = reverseComplement(fasta)
        ssr_seq = reverseComplement(ssr_seq)
    all_pos = [i for i in range(len(fasta)) if fasta.startswith(ssr_seq, i)]
    ssr_pos = []
    for pos in all_pos:
        pseq = get_pseq(fasta, pos, dictBam, ssr_seq)
        print ('This is the previous sequence')
        print (pseq)
        ssr_read = get_ssrRead(pseq, dictBam, pos)
        gpos = get_ssrPos(ssr_read, pos, fasta, ssr_seq, pseq, dictBam)
        ssr_pos.append(gpos)
    return (ssr_pos)


def get_pseq(fasta, pos, dictBam, ssr_seq):
    pseq = fasta[(pos - 20):pos]
    if not any([c['seq'].find(pseq) >= 0 for c in dictBam]):
        pseq = fasta[(pos + len(ssr_seq)):(pos + len(ssr_seq) + 20)]
    if not any([c['seq'].find(pseq) >= 0 for c in dictBam]):
        pseq = False
    return (pseq)


def get_ssrRead(pseq, dictBam, pos):
    if pseq is False:
        read_len = 0
        ssr_read = [c for c in dictBam if c['rpos'] == (len(dictBam)-1)][0]
        for read in dictBam:
            read_len += len(read['seq'])
            if read_len >= pos:
                ssr_read = read
    else:
        ssr_read = [read for read in dictBam if read['seq'].find(pseq) >= 0][0]
    return (ssr_read)


def get_ssrPos(ssr_read, pos, fasta, ssr_seq, pseq, dictBam):
    go = [c['strand'] for c in dictBam if c['rpos'] == 0][0]
    cigarParsed = parse_cigar(ssr_read['cigar'])
    if pseq is False:
        # gpos = ssr_read['gs'] if ssr_read['strand'] == '+' else ssr_read['ge']
        gpos = ssr_read['gs']
        gpos = '{}:{}-{}'.format(ssr_read['chrom'], gpos, gpos + len(ssr_seq))
        return (gpos)
    # pseq_pos = fasta.find(pseq) if [c['strand'] for c in dictBam if c['rpos'] == 0][0] == '+' else fasta.find(reverseComplement(pseq))
    pseq_pos = fasta.find(pseq)
    if pos > pseq_pos:
        rpos = ssr_read['seq'].find(pseq) + 20
    else:
        rpos = ssr_read['seq'].find(pseq) - 20 - len(ssr_seq)
    if rpos < 0:
        # gpos = ssr_read['gs'] if ssr_read['strand'] == '+' else ssr_read['ge']
        gpos = ssr_read['gs']
        gpos = '{}:{}-{}'.format(ssr_read['chrom'], gpos, gpos + len(ssr_seq))
        return (gpos)
    if rpos > ssr_read['re']:
        # gpos = ssr_read['ge'] if ssr_read['strand'] == '+' else ssr_read['gs']
        gpos = ssr_read['ge']
        gpos = '{}:{}-{}'.format(ssr_read['chrom'], gpos, gpos + len(ssr_seq))
        return (gpos)
    gpos = 0
    gposM = 0
    for c in cigarParsed:
        if c[0] == 'M' and gposM + c[1] >= rpos:
            break
        if c[0] == 'S' or c[0] == 'H':
            continue
        elif c[0] == 'M' or c[0] == 'I':
            gposM += c[1]
            gpos += c[1]
        else:
            gpos += c[1]

    gpos = gpos + (rpos - gposM) + ssr_read['gs']
    gpos = '{}:{}-{}'.format(ssr_read['chrom'], gpos, gpos + len(ssr_seq))
    return (gpos)


def _getSsrType(fasta):
    return ([c for c in TARGET_SEQ if fasta.find(TARGET_SEQ[c]) > 0])


def mergeINS(vcf, read, dictBam):

    go = [c['strand'] for c in dictBam if c['rpos'] == 0][0]
    vcf_out = []
    ssr_ins = []
    vcf_list = [c.split('\t') for c in vcf]
    # Check if there is a SSR and INS
    if any(['INS:LOXP' in c[7] or 'INS:FRT' in c[6] for c in vcf_list]) and any(['SVTYPE=INS' in x for x in [j for j in [c[7].split(';') for c in vcf_list]]]):
        for var in vcf_list:
            if 'INS' not in var[7]:
                vcf_out.append(var)
            else:
                ssr_ins.append(var)
        ssr_ins = sorted(ssr_ins, key=lambda k: (k[0], k[1]))
        i = 0
        ssr_n = 0
        while i < len(ssr_ins):
            var = ssr_ins[i]
            if i != len(ssr_ins) - 1:
                nvar = ssr_ins[i + 1]
            else:
                vcf_out.append(var)
                break
            # Look if INS and LOXP share coordinates
            if [c for c in var[7].split(';') if 'SVTYPE' in c][0].split('=')[1] == 'INS:LOXP' or [c for c in var[7].split(';') if 'SVTYPE' in c][0].split('=')[1] == 'INS:FRT':
                if [c for c in nvar[7].split(';') if 'SVTYPE' in c][0].split('=')[1] == 'INS' and [var[0], var[1]] == [nvar[0], nvar[1]]:
                    # Check which one comes first in the read
                    print ('This is the read')
                    print (read)
                    print ('And these the sequences')
                    print (var[4][1:])
                    print (nvar[4][1:])
                    varSeqPos = [j for j in range(len(read)) if read.startswith(var[4][1:], j)][ssr_n] if go == '+' else [j for j in range(len(read)) if reverseComplement(read).startswith(reverseComplement(var[4][1:]), j)][ssr_n]
                    nvarSeqPos = read.find(nvar[4][1:]) if go == '+' else reverseComplement(read).find(nvar[4][1:])
                    # if nvarSeqPos < 0:
                    #     print ('It is smaller than 0')
                    #     nvarSeqPos = reverseComplement(read).find(reverseComplement(var[4][1:]))

                    # Merge them in the same order
                    print ('These are the positions of the sequences. Should be higher than 0')
                    print (varSeqPos)
                    print (nvarSeqPos)
                    if varSeqPos < nvarSeqPos:
                        seq = var[4] + nvar[4][1:]
                    else:
                        seq = nvar[4] + var[4][1:]
                    vcf_out.append([var[0], var[1], var[2], var[3], seq, var[5], var[6], var[7]])
                    i += 1
                else:
                    vcf_out.append(var)
                ssr_n += 1
            if [c for c in var[7].split(';') if 'SVTYPE' in c][0].split('=')[1] == 'INS':
                if ([c for c in nvar[7].split(';') if 'SVTYPE' in c][0].split('=')[1] == 'INS:LOXP' or [c for c in nvar[7].split(';') if 'SVTYPE' in c][0].split('=')[1] == 'INS:FRT') and [var[0], var[1]] == [nvar[0], nvar[1]]:
                    # Check which one comes first in the read
                    print ('This is the read')
                    print (read)
                    varSeqPos = read.find(var[4][1:]) if go == '+' else reverseComplement(read).find(var[4][1:])
                    # if varSeqPos < 0:
                    #     varSeqPos = read.find(reverseComplement(var[4][1:]))
                    nvarSeqPos = [j for j in range(len(read)) if read.startswith(var[4][1:], j)][ssr_n] if go == '+' else [j for j in range(len(read)) if reverseComplement(read).startswith(reverseComplement(var[4][1:]), j)][ssr_n]
                    print ('These are the positions of the sequences. Should be higher than 0')
                    print (varSeqPos)
                    print (nvarSeqPos)
                    # Merge them in the same order
                    if varSeqPos < nvarSeqPos:
                        seq = var[4] + nvar[4][1:]
                    else:
                        seq = nvar[4] + var[4][1:]
                    vcf_out.append([var[0], var[1], var[2], var[3], seq, var[5], var[6], var[7]])
                    i += 1
                else:
                    vcf_out.append(var)
            i += 1
        return(['\t'.join(c) for c in vcf_out])
    else:
        return (vcf)
