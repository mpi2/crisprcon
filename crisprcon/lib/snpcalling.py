
import os
import subprocess

FNULL = open(os.devnull, 'w')


def snpCalling(in_bam, out_vcf, ref):
    '''\
    Samtools mpileup and Bcftools call to produce a sort and normalized VCF

    usage:
        snpCalling -r <reference> -i <input> [-o <output>]

    Options:
        -r, --reference <reference>  Reference file in fasta format
        -i, --input <input>          The input file in BAM format
        -o, --output <output>        The output file in VCF format
                                     [default: stdout]
    '''

    p1 = subprocess.Popen(['samtools', 'mpileup',
                           '-g',
                           '-E',
                           '-Q', '0',
                           '-p',
                           '-f', ref,
                           in_bam],
                          stdout=subprocess.PIPE,
                          stderr=FNULL)
    p2 = subprocess.Popen(['bcftools', 'call',
                           '-mv',
                           '-V', 'indels',
                           '-Oz'],
                          stdin=p1.stdout,
                          stdout=subprocess.PIPE,
                          stderr=FNULL)
    p1.stdout.close()
    p3 = subprocess.Popen(['bcftools', 'norm',
                           '-c', 's',
                           '-f', ref,
                           '-Oz'],
                          stdin=p2.stdout,
                          stdout=subprocess.PIPE,
                          stderr=FNULL)
    p2.stdout.close()
    p4 = subprocess.Popen(['bcftools', 'view',
                           '-G',
                           '-Oz',
                           '-o', out_vcf],
                          stdin=p3.stdout,
                          stdout=subprocess.PIPE,
                          stderr=FNULL)
    p3.stdout.close()
    p4.communicate()[0]

    subprocess.call(['tabix',
                     '-p', 'vcf',
                     out_vcf],
                    stderr=FNULL)
