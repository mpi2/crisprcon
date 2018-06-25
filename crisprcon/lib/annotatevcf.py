
import os
import subprocess
from .varfilter import _parse_vcf
from .infervariants import vcf_header


def annotatevcf(vcf, gff, ref, temp):
    vcf_intersect = temp + vcf.split('/')[-1] + '.intersect'
    vcf_temp = temp + vcf.split('/')[-1] + '.temp'
    vcf_ann = temp + vcf.split('/')[-1] + '.csq'
    csq(vcf, vcf_ann, gff, ref)
    add_exon_csq(vcf_ann, gff, vcf_intersect)
    with open(vcf_intersect) as reader:
        dictreader = _parse_vcf(reader)
        # Write out file
        vcf_format(dictreader, vcf_temp)

    subprocess.call(['bcftools', 'view',
                     vcf_temp,
                     '-Oz',
                     '-o',
                     vcf_temp])
    subprocess.call(['tabix',
                     '-p',
                     'vcf',
                     vcf_temp])
    subprocess.call(['bcftools', 'concat',
                     vcf_temp,
                     vcf_ann,
                     '-D',
                     '-a',
                     '-Oz',
                     '-o', vcf])


def csq(vcf_in, vcf_out, gff, ref):
    '''\
    Bcftools csq to annotate the variants calls from the merged VCF


    Options:
        vcf                         VCF generated by Split Read Alignment
        MM_vcf                      VCF generated by Alignment Read Mismatching
        out_vcf                     VCF with merged calls
    '''

    subprocess.call(['bcftools', 'csq',
                     '-g', gff,
                     '-f', ref,
                     '-Oz',
                     '-o', vcf_out,
                     vcf_in])

    subprocess.call(['tabix',
                     '-p', 'vcf',
                     vcf_out])


def add_exon_csq(vcf, gff, output):
    if os.path.isfile(output):
        os.remove(output)
    log = open(output, 'w')
    log.write('\t'.join(['#{}'.format(n) for n in range(1, 18)]) + '\n')
    log.flush()

    p1 = subprocess.Popen(['bedtools', 'intersect',
                           '-a', vcf,
                           '-b', gff,
                           '-loj'],
                          stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['awk', '$11 == "exon" || $11 == "five_prime_UTR" || $11 == "three_prime_UTR"'],
                     stdin=p1.stdout,
                     stdout=log,
                     stderr=log)
    ret_code = p2.wait()
    log.flush()
    p1.stdout.close()


def vcf_format(dictreader, vcf_out):

    _vcf_fields = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')

    with open(vcf_out, 'w') as writer:

        writer.write(vcf_header)
        # writer.write('#{}\n'.format('\t'.join(_vcf_fields)))

        output_vcf = []
        for line in dictreader:
            # Take vcf fields from the dictreader
            CHROM, POS, ID, REF, ALT, QUAL, FILTER = line['1'], line['#2'], line['#3'], line['#4'], line['#5'], line['#6'], line['#7']

            # Check if current line is already in the output_vcf
            if len(output_vcf) > 0:
                if [output_vcf[-1][0], output_vcf[-1][1], output_vcf[-1][3], output_vcf[-1][4]] == [CHROM, POS, REF, ALT]:
                    info = output_vcf[-1][7].split(';')
                else:
                    info = line['#8'].split(';')
            else:
                info = line['#8'].split(';')

            bcsq = info.pop().split('=')[1].split(',')

            if any('@' in c for c in bcsq):
                if len(output_vcf) == 0:
                    output_vcf.append([CHROM, POS, ID, REF, ALT, QUAL, FILTER, line['#8']])
                    continue
                else:
                    if [output_vcf[-1][0], output_vcf[-1][1], output_vcf[-1][3], output_vcf[-1][4]] != [CHROM, POS, REF, ALT]:
                        output_vcf.append([CHROM, POS, ID, REF, ALT, QUAL, FILTER, line['#8']])
                        continue
                    elif [output_vcf[-1][0], output_vcf[-1][1], output_vcf[-1][3], output_vcf[-1][4]] == [CHROM, POS, REF, ALT]:
                        continue

            bcsq = [c.split('|') for c in bcsq]

            biotype = line['#11']
            bt_info = line['#17'].split(';')
            if biotype == 'exon' and any('rank' in c for c in bt_info):
                biotype = '{}_{}'.format(biotype, [c for c in bt_info if 'rank' in c][0].split('=')[1])
            st = int(line['#2'])
            try:
                end = int([c for c in info if 'END' in c][0].split('=')[1])
            except:
                end = st
            btst = int(line['#12'])
            btend = int(line['#13'])
            bt_transcript = [c for c in bt_info if 'Parent=transcript' in c][0].split(':')[1]

            if st <= btst and end >= btend:
                p_mut = str(1.0)
            else:
                p_mut = str((min(end, btend) - max(st, btst))/(btend - btst))

            # Introduce the biotype information in the right transcript inside the BCSQ field
            for n, j in enumerate(bcsq):
                if bt_transcript in j and j[-2:] != [biotype, p_mut]:
                    if 'exon' in j[-2] or 'UTR' in j[-2]:
                        bcsq[n][-2] = '{}&{}'.format(j[-2], biotype)
                        bcsq[n][-1] = '{}&{}'.format(j[-1], p_mut)
                    else:
                        bcsq[n] = j + [biotype, p_mut]

            # Reinsert the BCSQ field into the INFO field
            bcsq = 'BCSQ=' + ','.join(['|'.join(c) for c in bcsq])
            info.append(bcsq)

            # Append to output_vcf if line is not present or subsitute the line if already present
            if len(output_vcf) == 0:
                output_vcf.append([CHROM, POS, ID, REF, ALT, QUAL, FILTER, ';'.join(info)])
            else:
                if output_vcf[-1][0:7] == [CHROM, POS, ID, REF, ALT, QUAL, FILTER]:
                    output_vcf[-1] = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, ';'.join(info)]
                else:
                    output_vcf.append([CHROM, POS, ID, REF, ALT, QUAL, FILTER, ';'.join(info)])

        # Sort all results
        output_vcf.sort()
        output = "\n".join(["\t".join(map(str, vcf_row)) for vcf_row in output_vcf])
        # Write record
        writer.write(output + '\n')