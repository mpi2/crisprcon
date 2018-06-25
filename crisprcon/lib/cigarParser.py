
# import pysam
import re
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


def parse_md(md):
    md_parsed = []
    mdList = list(md)
    mdList = [int(i) if i.isalpha() is False else i for i in mdList]
    prev = 0
    for i in range(len(mdList)):
        if md[i].isdigit():
            continue
        else:
            md_parsed.append([mdList[i], int(''.join([str(c) for c in mdList[prev:i]]))])
            prev = i + 1
    return (md_parsed)


def get_deletion(tag, del_count):
    try:
        del_seq = ''.join([c[1].split('^')[del_count] for c in tag if c[0] == 'MD'])
        del_seq = re.split('[0-9]', del_seq)[0]
        return (del_seq)
    except:
        return ('N')
