
from .cigarParser import parse_cigar
import subprocess
import os.path
import collections
import shutil
import fileinput
import os
import glob
from .ssr import _isSSR, _ssr_analysis

FNULL = open(os.devnull, 'w')


def main(in_fa, out_bam, ref, mref, chrom, ts, te, temp, prefix):

    # Check for site specific recombination and remove it if present
    if _isSSR(in_fa):
        with open(in_fa) as fasta:
            with open(temp + '{}.noSSR.fasta'.format(prefix), 'w') as writer:
                for line in fasta:
                    if line.startswith('>'):
                        writer.write(line)
                    else:
                        writer.write(_ssr_analysis(line))
        in_fa = temp + '{}.noSSR.fasta'.format(prefix)

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

    for files in glob.glob(in_fa + '.*'):
        os.remove(files)
    print (dictBam)
    return (dictBam)


def bwaAlign(in_fa, ref, mref, temp, chrom, ts, te, prefix):
    """\
    Use BWA to align the read and all it's soft clipped fragments reiteratively

    """
    sam = []
    interval = 10000
    with open(in_fa) as fasta:
        for line in fasta:
            if line.startswith('>'):
                fasta_header = line
                continue
            if len(line) >= 90:
                read = largeAlign(in_fa, ref)
                if not read:
                    continue
                sam.append(read)
                cigar = parse_cigar(read.split('\t')[5])
                # If there is any clipping, recover the clipped fragment and call bwaAlign with it.
                if any(['S' in c for c in cigar]):
                    readPos = 0
                    it = 1
                    for event in cigar:
                        if event[0] == 'S':
                            output = open(temp + '{}.{}To{}.fasta'.format(prefix, readPos, event[1]), 'w')
                            if list(flagAsBin(int(read.split('\t')[1])))[-5] == '1':
                                line_rev = reverseComplement(line)
                                output.write(fasta_header + reverseComplement(line_rev[readPos:event[1]]) + '\n') if it != len(cigar) else output.write(fasta_header + reverseComplement(line_rev[readPos:len(line_rev)]))
                                output.close()
                            else:
                                output.write(fasta_header + line[readPos:event[1]] + '\n') if it != len(cigar) else output.write(fasta_header + line[readPos:len(line)])
                                output.close()
                            sam.append(bwaAlign(in_fa=temp + '{}.{}To{}.fasta'.format(prefix, readPos, event[1]), ref=ref, mref=mref, temp=temp, chrom=chrom, ts=ts, te=te, prefix=prefix))
                            readPos += event[1]
                        else:
                            readPos += event[1]
                        it += 1

            if len(line) > 50 and len(line) < 90:
            # if len(line) > 50:
                if not os.path.exists(mref):
                    trimRef(ref, chrom, int(ts) - interval, int(te) + interval, mref)
                read = largeAlign(in_fa, mref)
                if not read:
                    continue
                sam.append(read)
                cigar = parse_cigar(read.split('\t')[5])
                if any(['S' in c for c in cigar]):
                    readPos = 0
                    it = 1
                    for event in cigar:
                        if event[0] == 'S':
                            output = open(temp + '{}.{}To{}.fasta'.format(prefix, readPos, event[1]), 'w')
                            if list(flagAsBin(int(read.split('\t')[1])))[-5] == '1':
                                line_rev = reverseComplement(line)
                                output.write(fasta_header + reverseComplement(line_rev[readPos:event[1]]) + '\n') if it != len(cigar) else output.write(fasta_header + reverseComplement(line_rev[readPos:len(line_rev)]))
                                output.close()
                            else:
                                output.write(fasta_header + line[readPos:event[1]] + '\n') if it != len(cigar) else output.write(fasta_header + line[readPos:len(line)])
                                output.close()
                            sam.append(bwaAlign(in_fa=temp + '{}.{}To{}.fasta'.format(prefix, readPos, event[1]), ref=ref, mref=mref, temp=temp, chrom=chrom, ts=ts, te=te, prefix=prefix))
                            readPos += event[1]
                        else:
                            readPos += event[1]
                        it += 1

            if len(line) < 50 and len(line) >= 20:
                if not os.path.exists(mref):
                    trimRef(ref, chrom, int(ts) - interval, int(te) + interval, mref)
                read = smallAlign(in_fa, mref)
                if not read:
                    continue
                sam.append(read)
                cigar = parse_cigar(read.split('\t')[5])
                if any(['S' in c for c in cigar]):
                    readPos = 0
                    it = 1
                    for event in cigar:
                        if event[0] == 'S':
                            output = open(temp + '{}.{}To{}.fasta'.format(prefix, readPos, event[1]), 'w')
                            if list(flagAsBin(int(read.split('\t')[1])))[-5] == '1':
                                line_rev = reverseComplement(line)
                                output.write(fasta_header + reverseComplement(line_rev[readPos:event[1]]) + '\n') if it != len(cigar) else output.write(fasta_header + reverseComplement(line_rev[readPos:len(line_rev)]))
                                output.close()
                            else:
                                output.write(fasta_header + line[readPos:event[1]] + '\n') if it != len(cigar) else output.write(fasta_header + line[readPos:len(line)])
                                output.close()
                            sam.append(bwaAlign(in_fa=temp + '{}.{}To{}.fasta'.format(prefix, readPos, event[1]), ref=ref, mref=mref, temp=temp, chrom=chrom, ts=ts, te=te, prefix=prefix))
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

    p1 = subprocess.Popen(['bwa', 'mem', '-c', '250', '-v', '1', ref, in_fa], stdout=subprocess.PIPE, stderr=FNULL)
    p2 = subprocess.Popen(['samtools', 'view', '-F', '2048', '-F', '4', '-bS', '-'], stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    p3 = subprocess.Popen(['samtools', 'sort', '-O', 'sam'], stdin=p2.stdout, stdout=subprocess.PIPE)
    p2.stdout.close()
    sam = p3.communicate()[0].decode('utf-8')
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


def smallAlign(in_fa, ref):
    '''\
    Call the BWA aln aligner for small fragments

    '''

    p1 = subprocess.Popen(['bwa', 'aln', '-l', '20', ref, in_fa], stdout=subprocess.PIPE, stderr=FNULL)
    p2 = subprocess.Popen(['bwa', 'samse', ref, '-', in_fa], stdin=p1.stdout, stdout=subprocess.PIPE, stderr=FNULL)
    p1.stdout.close()
    p3 = subprocess.Popen(['samtools', 'view', '-q', '1', '-F', '2048', '-F', '4', '-bS', '-'], stdin=p2.stdout, stdout=subprocess.PIPE)
    p2.stdout.close()
    p4 = subprocess.Popen(['samtools', 'sort', '-O', 'sam'], stdin=p3.stdout, stdout=subprocess.PIPE)
    p3.stdout.close()
    sam = p4.communicate()[0].decode('utf-8')
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
    with open(in_fa) as fasta:
        for line in fasta:
            if not line.startswith('>'):
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
        rs = getOrder(in_fa, seq, strand, 0)
        re = int(rs) + len(seq)
        chrom = read[2]
        gs = int(read[3])
        ge = int(read[3]) + len(seq)
        for c in cigarParsed:
            if c[0] == 'I':
                ge -= c[1]
            if c[0] == 'D':
                ge += c[1]
        tags = read[12]

        d = {'chrom': chrom, 'strand': strand, 'gs': gs, 'ge': ge, 'cigar': cigar, 'seq': seq, 'tags': tags, 'rs': rs, 're': re}
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
        print ('This is the qc iteration')
        print (it)
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
