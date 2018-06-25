import subprocess


def svNorm(in_vcf, out_vcf, ref):
    '''\
    Bcftools sort and Bcftools sort to produce a sort and normalized VCF

    '''

    p1 = subprocess.Popen(['bcftools', 'view',
                           '-Oz',
                           in_vcf],
                          stdout=subprocess.PIPE)

    p2 = subprocess.Popen(['bcftools', 'sort', '-Oz', '-o', out_vcf], stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    p2.communicate()[0]

    subprocess.call(['tabix',
                     '-p', 'vcf',
                     out_vcf])


# def svNorm(in_vcf, out_vcf, ref):
#     '''\
#     Bcftools sort and Bcftools sort to produce a sort and normalized VCF
#
#     '''
#
#     p1 = subprocess.Popen(['bcftools', 'view',
#                            '-Oz',
#                            in_vcf],
#                           stdout=subprocess.PIPE)
#
#     p2 = subprocess.Popen(['bcftools', 'sort', '-Oz'], stdin=p1.stdout, stdout=subprocess.PIPE)
#     p1.stdout.close()
#     p3 = subprocess.Popen(['bcftools', 'norm',
#                            '-c', 's',
#                            '-f', ref,
#                            '-Oz',
#                            '-o', out_vcf],
#                           stdin=p2.stdout,
#                           stdout=subprocess.PIPE)
#     p2.stdout.close()
#     p3.communicate()[0]
#
#     subprocess.call(['tabix',
#                      '-p', 'vcf',
#                      out_vcf])
