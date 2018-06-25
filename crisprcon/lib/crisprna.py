

from .resolvegrna import resolvegRNA


def crispRNA(output_vcf, ref_genome, gff, dir, temp_dir, prefix, grna, donor):

    '''
    gRNA consequence prediction. It takes the coordinates and sequence of the CRISPR gRNA and, if present,
    the sequence of the donor, and predicts the consequence that the cut will have in the genome.
    '''

    mutFa = resolvegRNA(grna=grna,
                        ref=ref_genome,
                        donors=donor
                        )
    mfa = open(temp_dir + prefix + '.crispr.fasta', 'w')
    mfa.write('>{}_CRISPR\n{}\n'.format(prefix, mutFa))
    mfa.close
