
from .cigarParser import parse_cigar
import subprocess
import os.path
import collections
import shutil
import fileinput
import os
import glob
from .ssr import _isSSR, _ssr_analysis
from difflib import SequenceMatcher

FNULL = open(os.devnull, 'w')


def main(input_fa, out_bam, ref, mref, chrom, ts, te, temp, prefix):

    # Check for site specific recombination and remove it if present
    with open(input_fa) as fasta:

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
        print ('Instead of aligning ')
        print (in_fa)

        if _isSSR(input_fa):
            header = in_fa.split('\n')[0] + '\n'
            with open(temp + '{}.noSSR.fasta'.format(prefix), 'w') as writer:
                seq_line_noSSR = _ssr_analysis(seq_line)

                if header.startswith('>'):
                    seq_line = seq_line_noSSR
                    in_fa = '{}{}'.format(header, seq_line).rstrip()
                    print ('We are aligning')
                    print(in_fa)
                if header.startswith('@'):
                    matches = difflib.SequenceMatcher(None, seq_line, seq_line_noSSR).get_matching_blocks()
                    ssr_start = matches[0].a + matches[0].size
                    ssr_end = matches[1].a

                    seq_line = seq_line_noSSR
                    qual_line = qual_line[0:start] + qual_line[end:]
                    in_fa = '{}{}{}{}'.format(header, seq_line, qual_header, qual_line).rstrip()
                writer.write(in_fa)
            input_fa = temp + '{}.noSSR.fasta'.format(prefix)

        # Index reference file
        bwaIndex(ref)
        refIndex(ref)

        # Perform the alignment using BWA in a recursive manner
        sam = "\n".join(flatten(bwaAlign(in_fa, ref, mref, temp, chrom, ts, te, prefix)))

        # Exit pipeline of no alignments are found
    if len(sam) == 0:
        exit('ERROR: Failing to align. No alignments found.')

    # Create sam file (for visualization or evidence. Not used in the pipeline)
    samHeader(ref=ref, sam=sam, out_bam=out_bam)

    # Create a dictionary of the alignments that will be used in the variant calling
    dictBam = sam2dict(sam=sam, in_fa=in_fa)
    dictBam = sorted(dictBam, key=lambda k: (k['rpos']))

    dictBam = dictBam_QC(dictBam, in_fa)  # Just check if the order of the reads is correct (non_overlapping starts and ends)

    # Invert the order of the reads if the homology arms aligned to the reverse strand

    if [c['strand'] for c in dictBam if c['rpos'] == 0][0] == '-':
        for i in range(len(dictBam)):
            dictBam[i]['rpos'] = len(dictBam) - 1 - i

    # Remove temporary files
    for files in glob.glob(mref + '*'):
        os.remove(files)

    for files in glob.glob(input_fa + '.*'):
        os.remove(files)
    return (dictBam)


def bwaAlign(in_fa, ref, mref, temp, chrom, ts, te, prefix):
    """\
    Use BWA to align the read and all it's soft clipped fragments reiteratively

    """
    sam = []
    interval = 10000

    header = in_fa.split('\n')[0]
    seq_line = in_fa.split('\n')[1]
    qual_header = ''
    qual_line = ''
    if header.startswith('@'):
        qual_header = in_fa.split('\n')[2]
        qual_line = in_fa.split('\n')[3]

    if len(seq_line) >= 90:
        read = largeAlign(in_fa, ref)
        if not read:
            return(sam)
        sam.append(read)
        cigar = parse_cigar(read.split('\t')[5])
        # If there is any clipping, recover the clipped fragment and call bwaAlign with it.
        if any(['S' in c for c in cigar]):
            readPos = 0
            it = 1
            for event in cigar:
                if event[0] == 'S':
                    if list(flagAsBin(int(read.split('\t')[1])))[-5] == '1':
                        if header.startswith('>'):
                            line_rev = reverseComplement(seq_line)
                            out_fa = header + '\n' + reverseComplement(line_rev[readPos:event[1]]) if it != len(cigar) else header + '\n' + reverseComplement(line_rev[readPos:len(line_rev)])
                        elif header.startswith('@'):
                            line_rev = reverseComplement(seq_line)
                            qual_rev = qual_line[::-1]
                            out_fa = header + '\n' + reverseComplement(line_rev[readPos:event[1]]) + '\n' + qual_header + '\n' + qual_rev[readPos:event[1]][::-1] if it != len(cigar) else header + '\n' + reverseComplement(line_rev[readPos:len(line_rev)]) + '\n' + qual_header + '\n' + qual_rev[readPos:len(line_rev)][::-1]
                    else:
                        if header.startswith('>'):
                            out_fa = header + '\n' + seq_line[readPos:event[1]] if it != len(cigar) else header + '\n' + seq_line[readPos:len(seq_line)]
                        elif header.startswith('@'):
                            out_fa = header + '\n' + seq_line[readPos:event[1]] + '\n' + qual_header + '\n' + qual_line[readPos:event[1]] if it != len(cigar) else header + '\n' + seq_line[readPos:len(seq_line)] + '\n' + qual_header + '\n' + qual_line[readPos:len(seq_line)]
                    sam.append(bwaAlign(in_fa=out_fa, ref=ref, mref=mref, temp=temp, chrom=chrom, ts=ts, te=te, prefix=prefix))
                    readPos += event[1]
                else:
                    readPos += event[1]
                it += 1
    if len(seq_line) >= 50 and len(seq_line) < 90:
    # if len(seq_line) > 50:
        if not os.path.exists(mref):
            trimRef(ref, chrom, int(ts) - interval, int(te) + interval, mref)
        read = largeAlign(in_fa, mref)
        if not read:
            return(sam)
        sam.append(read)
        cigar = parse_cigar(read.split('\t')[5])
        if any(['S' in c for c in cigar]):
            readPos = 0
            it = 1
            for event in cigar:
                if event[0] == 'S':
                    if list(flagAsBin(int(read.split('\t')[1])))[-5] == '1':
                        if header.startswith('>'):
                            line_rev = reverseComplement(seq_line)
                            out_fa = header + '\n' + reverseComplement(line_rev[readPos:event[1]]) if it != len(cigar) else header + '\n' + reverseComplement(line_rev[readPos:len(line_rev)])
                        elif header.startswith('@'):
                            line_rev = reverseComplement(seq_line)
                            qual_rev = qual_line[::-1]
                            out_fa = header + '\n' + reverseComplement(line_rev[readPos:event[1]]) + '\n' + qual_header + '\n' + qual_rev[readPos:event[1]][::-1] if it != len(cigar) else header + '\n' + reverseComplement(line_rev[readPos:len(line_rev)]) + '\n' + qual_header + '\n' + qual_rev[readPos:len(line_rev)][::-1]
                    else:
                        if header.startswith('>'):
                            out_fa = header + '\n' + seq_line[readPos:event[1]] if it != len(cigar) else header + '\n' + seq_line[readPos:len(seq_line)]
                        elif header.startswith('@'):
                            out_fa = header + '\n' + seq_line[readPos:event[1]] + '\n' + qual_header + '\n' + qual_line[readPos:event[1]] if it != len(cigar) else header + '\n' + seq_line[readPos:len(seq_line)] + '\n' + qual_header + '\n' + qual_line[readPos:len(seq_line)]
                    sam.append(bwaAlign(in_fa=out_fa, ref=ref, mref=mref, temp=temp, chrom=chrom, ts=ts, te=te, prefix=prefix))
                    readPos += event[1]
                else:
                    readPos += event[1]
                it += 1

    if len(seq_line) < 50 and len(seq_line) >= 20:
        if not os.path.exists(mref):
            trimRef(ref, chrom, int(ts) - interval, int(te) + interval, mref)
        read = smallAlign(in_fa, mref, temp, prefix)
        if not read:
            return(sam)
        sam.append(read)
        cigar = parse_cigar(read.split('\t')[5])
        if any(['S' in c for c in cigar]):
            readPos = 0
            it = 1
            for event in cigar:
                if event[0] == 'S':
                    if list(flagAsBin(int(read.split('\t')[1])))[-5] == '1':
                        if header.startswith('>'):
                            line_rev = reverseComplement(seq_line)
                            out_fa = header + '\n' + reverseComplement(line_rev[readPos:event[1]]) if it != len(cigar) else header + '\n' + reverseComplement(line_rev[readPos:len(line_rev)])
                        elif header.startswith('@'):
                            line_rev = reverseComplement(seq_line)
                            qual_rev = qual_line[::-1]
                            out_fa = header + '\n' + reverseComplement(line_rev[readPos:event[1]]) + '\n' + qual_header + '\n' + qual_rev[readPos:event[1]][::-1] if it != len(cigar) else header + '\n' + reverseComplement(line_rev[readPos:len(line_rev)]) + '\n' + qual_header + '\n' + qual_rev[readPos:len(line_rev)][::-1]
                    else:
                        if header.startswith('>'):
                            out_fa = header + '\n' + seq_line[readPos:event[1]] if it != len(cigar) else header + '\n' + seq_line[readPos:len(seq_line)]
                        elif header.startswith('@'):
                            out_fa = header + '\n' + seq_line[readPos:event[1]] + '\n' + qual_header + '\n' + qual_line[readPos:event[1]] if it != len(cigar) else header + '\n' + seq_line[readPos:len(seq_line)] + '\n' + qual_header + '\n' + qual_line[readPos:len(seq_line)]
                    sam.append(bwaAlign(in_fa=out_fa, ref=ref, mref=mref, temp=temp, chrom=chrom, ts=ts, te=te, prefix=prefix))
                    readPos += event[1]
                else:
                    readPos += event[1]
                it += 1

    return (sam)


def bwaIndex(ref):
    """\
    Check for bwa mem indexes for the reference genome and build them if they don't exist

    """
    index_list = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    if not all(os.path.exists(ref + idx) for idx in index_list):
        subprocess.call(['bwa', 'index', ref], stderr=FNULL, stdout=FNULL)


def refIndex(ref):
    """\
    Check for the fasta .fai index for the reference genome and build it if it doesn't exist

    """

    if not os.path.exists(ref + '.fai'):
        subprocess.call(['samtools', 'faidx', ref], stderr=FNULL, stdout=FNULL)


def largeAlign(in_fa, ref):
    '''\
    Call the BWA mem aligner for long fragments

    '''
    p1 = subprocess.run(['bwa', 'mem', '-c', '250', '-v', '1', '-k', '14', '-W', '20', '-r', '10',
                        '-A', '1', '-L', '0', '-B', '4', '-O', '6', '-E', '6', ref, '-'],
                        stdout=subprocess.PIPE, input=in_fa.encode(), stderr=FNULL)
    p2 = subprocess.run(['samtools', 'view', '-F', '2048', '-F', '4', '-bS', '-'], input=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.run(['samtools', 'sort', '-O', 'sam'], input=p2.stdout, stdout=subprocess.PIPE)  # TODO: do we need this?
    sam = p3.stdout.decode("utf-8")

    read = '\n'.join([c for c in sam.split('\n') if not c.startswith('@') and len(c) is not 0])
    if read:
        if len(read.split('\t')[2].split(':')) > 1:
            sread = read.split('\t')
            chrom = sread[2].split(':')[0]
            pos = int(sread[2].split(':')[1].split('-')[0]) - 1 + int(sread[3])
            sread[2] = str(chrom)
            sread[3] = str(pos)
            read = '\t'.join(sread)
    return (read)


def smallAlign(in_fa, ref, temp, prefix):
    '''\
    Call the BWA aln aligner for small fragments

    '''

    p1 = subprocess.run(['bwa', 'aln', '-l', '20', ref, '-', '-f', '{}{}.aln'.format(temp, prefix)],
                        input=in_fa.encode(), stdout=subprocess.PIPE, stderr=FNULL)
    p2 = subprocess.run(['bwa', 'samse', ref, '{}{}.aln'.format(temp, prefix), '-'],
                        input=in_fa.encode(), stdout=subprocess.PIPE, stderr=FNULL)
    # p2 = subprocess.run(['bwa samse reference/C57BL_6J.fa <(bwa aln -l 20 reference/C57BL_6J.fa -) -'],
    #                     shell=True, executable="/bin/bash", input=p1.stdout, stdout=subprocess.PIPE, stderr=FNULL)
    p3 = subprocess.run(['samtools', 'view', '-q', '1', '-F', '2048', '-F', '4', '-bS', '-'],
                        input=p2.stdout, stdout=subprocess.PIPE)
    p4 = subprocess.run(['samtools', 'sort', '-O', 'sam'],
                        input=p3.stdout, stdout=subprocess.PIPE)
    sam = p4.stdout.decode("utf-8")

    read = '\n'.join([c for c in sam.split('\n') if not c.startswith('@') and len(c) is not 0])
    if read:
        if len(read.split('\t')[2].split(':')) > 1:
            sread = read.split('\t')
            chrom = sread[2].split(':')[0]
            pos = int(sread[2].split(':')[1].split('-')[0]) - 1 + int(sread[3])
            sread[2] = str(chrom)
            sread[3] = str(pos)
            read = '\t'.join(sread)

    return (read)


def trimRef(ref, chrom, ts, te, mref):
    '''\
    Use samtools to extract the region between start and end

    '''

    region = '{}:{}-{}'.format(chrom, ts, te)
    f = open(mref, 'w')
    subprocess.call(['samtools', 'faidx', ref, region], stdout=f)
    f.close
    refIndex(mref)
    bwaIndex(mref)


def reverseComplement(pattern):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                  'Y': 'R', 'R': 'Y', 'W': 'W', 'S': 'S', 'K': 'M',
                  'M': 'K', 'D': 'H', 'V': 'B', 'H': 'D', 'B': 'V',
                  'X': 'X'}
    pattern = pattern.strip()
    return ''.join([complement[base] for base in pattern[::-1]])


def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


def samHeader(ref, sam, out_bam):
    '''\
    Add header to sam file and ouput a bam file

    '''

    p1 = subprocess.run(['cat'], stdout=subprocess.PIPE, input=sam.encode())
    p2 = subprocess.run(['samtools', 'view', '-ht', ref + '.fai', '-', '-bS'], stdout=subprocess.PIPE, input=p1.stdout)
    p3 = subprocess.run(['samtools', 'sort', '-O', 'bam', '-o', out_bam], input=p2.stdout, stdout=subprocess.PIPE)


def flagAsBin(n):
    """\
    convert a samtools flag to its binary rep (padded with 0's)

    """
    return str(bin(n))[2:].zfill(17)


def getOrder(in_fa, substr, strand, p_st):
    for line in in_fa.split('\n'):
        if line.startswith('>') or line.startswith('@'):
            continue
        else:
            line = list(line)
            for i in range(len(line)):
                if not line[i] in ['A', 'G', 'C', 'T']:
                    line[i] = 'N'
            line = ''.join(line)
            if strand == '+':
                return(line.find(substr, p_st) + 1)
            else:
                substr_rev = reverseComplement(substr)
                return(line.find(substr_rev, p_st) + 1)


def getParsedSeq(seq, cigar):
    if seq == '*':
        return('*')
    cigar = parse_cigar(cigar)
    start = int(cigar[0][1]) if cigar[0][0] == 'S' or cigar[0][0] == 'H' else 0
    end = -int(cigar[len(cigar) - 1][1]) if cigar[len(cigar) - 1][0] == 'S' or cigar[len(cigar) - 1][0] == 'H' else len(seq)

    pseq = seq[start:end]
    return (pseq)


def sam2dict(sam, in_fa):
    '''\
    It reads every line of the sam variable and creates a list of dictionaries with format:
    d = {'chrom': chromosome, 'strand': strand, 'gs': genome start, 'ge': genome end,
        'cigar': cigar string, 'seq': read sequence, 'tags': tags, 'rs': read start, 're': read end, 'rpos': read position}
    '''
    dictBam = []
    sam_list = [c.split('\t') for c in sam.split('\n')]
    n = 0
    for read in sam_list:
        cigar = read[5]
        cigarParsed = parse_cigar(cigar)
        flag = read[1]
        strand = '-' if list(flagAsBin(int(flag)))[-5] == '1' else '+'
        seq = getParsedSeq(read[9], cigar)
        qual = getParsedSeq(read[10], cigar)
        rs = getOrder(in_fa, seq, strand, 0)
        re = int(rs) + len(seq)
        chrom = read[2]
        gs = int(read[3])
        ge = int(read[3]) + len(seq)
        mq = read[4]
        for c in cigarParsed:
            if c[0] == 'I':
                ge -= c[1]
            if c[0] == 'D':
                ge += c[1]
        tags = read[12]

        d = {'chrom': chrom, 'strand': strand, 'gs': gs, 'ge': ge, 'cigar': cigar, 'seq': seq, 'qual': qual, 'mq': mq, 'tags': tags, 'rs': rs, 're': re}
        dictBam.append(d)
        n += 1
    dictBam = sorted(dictBam, key=lambda k: k['rs'])
    for i in range(len(dictBam)):
        dictBam[i]['rpos'] = i
    return (dictBam)


def dictBam_QC(dictBam, in_fa):
    qc = is_qc(dictBam)
    it = 0
    while not qc:
        it += 1
        dictBam = dictBam_reOrder(dictBam, in_fa)
        qc = is_qc(dictBam)
    return (dictBam)


def is_qc(dictBam):
    qc = True
    i = 0
    limit = len(dictBam)
    while qc and i < limit:
        if i == 0:
            pst = dictBam[i]['rs']
            pend = dictBam[i]['re']
            i += 1
            continue
        else:
            st = dictBam[i]['rs']
            end = dictBam[i]['re']
            if st < pend or [st, end] == [pst, pend]:
                qc = False
            i += 1
    return (qc)


def dictBam_reOrder(dictBam, in_fa):
    re = 0
    for i in range(len(dictBam)):
        seq = dictBam[i]['seq']
        strand = dictBam[i]['strand']
        st = getOrder(in_fa, seq, strand, re)
        dictBam[i]['rs'] = st
        dictBam[i]['re'] = int(st) + len(seq)
        re = int(st) + len(seq) - 1
    dictBam = sorted(dictBam, key=lambda k: (k['rs']))
    for i in range(len(dictBam)):
        dictBam[i]['rpos'] = i
    return (dictBam)
