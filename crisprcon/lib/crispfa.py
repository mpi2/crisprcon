from .alignment import main as alignment
from .snpcalling import snpCalling as snpCalling
from .infervariants import main as inferVariants
from .svnorm import svNorm as svNorm
from .mergevar import mergeVar as mergeVar
from .varfilter import varFilter as varFilter
from .annotatevcf import annotatevcf as annotateVCF


def crispFa(input_fa, output_vcf, ref_genome, gff, dir, temp_dir, prefix, chrom, ts, te, grna):

    """
    Variant caller for CRISPR edited alleles. The input is in FASTA format.
    """

    output_bam = temp_dir + prefix + '.bam'
    vcf_temp = temp_dir + prefix + 'vcf.gz.temp'

    mref = list(ref_genome.rpartition('/'))
    mref.insert(-1, prefix + '.masked.')
    mref = ''.join(mref)

    # Performing alignment
    dictBam = alignment(input_fa=input_fa,
                        out_bam=output_bam,
                        ref=ref_genome,
                        mref=mref,
                        chrom=chrom,
                        ts=ts,
                        te=te,
                        temp=temp_dir,
                        prefix=prefix)

    # Infering structural variants, SNPs, and indels
    inferVariants(dictBam=dictBam,
                  out_vcf=vcf_temp,
                  ref_genome=ref_genome,
                  in_fa=input_fa,
                  temp=temp_dir,
                  prefix=prefix)

    # Variant normalization
    svNorm(in_vcf=vcf_temp,
           out_vcf=vcf_temp,
           ref=ref_genome)

    # Filtering SNPs and small indels
    varFilter(dictBam=dictBam,
              inVCF=vcf_temp,
              outVCF=output_vcf,
              grna=grna,
              ref_genome=ref_genome)

    # Annotating variant consequence
    annotateVCF(vcf=output_vcf,
                gff=gff,
                ref=ref_genome,
                temp=temp_dir)
