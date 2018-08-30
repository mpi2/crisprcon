
# import pysam
import re
import numpy as np
# Parse the cigar string from the .bam file


# def parse_bam(in_bam):
#     # It reads every line of the BAM file and creates a list of dictionaries with format 'start' = start of read 'cigar' = cigar string
#     bam_list = pysam.Samfile(in_bam, 'rb')
#     dictBam = []
#     for read in bam_list:
#         cigar = ''.join(read.cigarstring)
#         chrom = read.reference_name
#         flag = read.flag
#         start = read.pos
#         end = read.aend
#         seq = read.seq
#         tags = read.tags
#         strand = '-' if int(flag) == 16 else '+'
#         d = {'chrom': chrom, 'strand': strand, 'start': start, 'end': end, 'cigar': cigar, 'seq': seq, 'tags': tags}
#         dictBam.append(d)
#     return (dictBam)


def parse_cigar(cigarString):
    # It parses a cigar string and makes a list with format ['character', numeric]
    cigarParsed = []
    cigarList = list(cigarString)
    cigarList = [int(i) if i.isalpha() is False else i for i in cigarList]
    prev = 0
    for i in range(len(cigarList)):
        if isinstance(cigarList[i], (int, float, complex)):
            continue
        else:
            cigarParsed.append([cigarList[i], int(''.join([str(c) for c in cigarList[prev:i]]))])
            prev = i + 1
    return (cigarParsed)


def RepresentsInt(i):
    '''
    Helper function that returns True if `i` is an int.
    '''
    try:
        int(i)
        return True
    except ValueError:
        return False


def parseMD(mdstr):
    '''
    Returns a numpy array with positions for where a `fromMM` base is, in the
    reference genome. This can then be used to intersect with the query string
    to find all desired mutations.
    '''
    # Split MD String at every single [ACGT] or ^:
    mdSub = re.sub(r'([\\^]*[ACGT]+)[0]*', ' \\1 ', mdstr)
    mdSplit = re.split('[ ]+', mdSub)
    mutArr = []
    # Iterate over Array and replace all mutations from the MD string with the
    # letter of the corresponding reference.
    # eg: 2G1 will produce the numpy array: 'M M G M'
    # All ^[ACGT]* by "D" and the number with a corresponding stretch of "M"
    for i in range(len(mdSplit)):
        mdPos = mdSplit[i]
        if len(mdPos) > 0 and RepresentsInt(mdPos):
            mutArr.extend(list('M'*int(mdPos)))
        elif re.match('\\^', mdPos):
            mutArr.extend(list('D'*(len(mdPos) - 1)))
        elif len(mdPos) == 1:
            mutArr.extend(mdPos)
        else:
            # I'm not yet quite sure, if this won't just break at some point.
            # In my BAM/SAM files I have seen rare cases with two consecutive
            # mismatches in the MD tag causing this series of ifs to report incorrect
            # positions if I don't catch this.
            mutArr.extend(list(mdPos))
    return(mutArr)

# def parse_md(md):
#     md_parsed = []
#     mdList = list(md)
#     mdList = [int(i) if i.isalpha() is False else i for i in mdList]
#     prev = 0
#     for i in range(len(mdList)):
#         if md[i].isdigit():
#             continue
#         else:
#             md_parsed.append([mdList[i], int(''.join([str(c) for c in mdList[prev:i]]))])
#             prev = i + 1
#     return (md_parsed)




def get_deletion(tag, del_count):
    try:
        del_seq = ''.join([c[1].split('^')[del_count] for c in tag if c[0] == 'MD'])
        del_seq = re.split('[0-9]', del_seq)[0]
        return (del_seq)
    except:
        return ('N')
