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
    mm_vcf = temp_dir + prefix + '.MM.vcf.gz'
    sr_vcf = temp_dir + prefix + '.SR.vcf'

    mref = list(ref_genome.rpartition('/'))
    mref.insert(-1, prefix + '.masked.')
    mref = ''.join(mref)

    # Performing alignment
    dictBam = alignment(in_fa=input_fa,
                        out_bam=output_bam,
                        ref=ref_genome,
                        mref=mref,
                        chrom=chrom,
                        ts=ts,
                        te=te,
                        temp=temp_dir,
                        prefix=prefix)
    # SNP calling
    snpCalling(in_bam=output_bam,
               out_vcf=mm_vcf,
               ref=ref_genome)

    # Infering structural variants and indels
    inferVariants(dictBam=dictBam,
                  out_vcf=sr_vcf,
                  ref_genome=ref_genome,
                  in_fa=input_fa,
                  temp=temp_dir,
                  prefix=prefix)

    # Variant normalization
    svNorm(in_vcf=sr_vcf,
           out_vcf=sr_vcf + '.gz',
           ref=ref_genome)

    # Merging variants
    mergeVar(SR_vcf=sr_vcf + '.gz',
             MM_vcf=mm_vcf,
             out_vcf=output_vcf)

    # Filtering SNPs and small indels
    varFilter(dictBam=dictBam,
              inVCF=output_vcf,
              outVCF=output_vcf,
              grna=grna,
              ref_genome=ref_genome)

    # Annotating variant consequence
    annotateVCF(vcf=output_vcf,
                gff=gff,
                ref=ref_genome,
                temp=temp_dir)
