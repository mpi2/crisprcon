#!/usr/bin/env python

from .infervariants import get_ref


def resolvegRNA(grna, ref, donors):
    '''
    Attached is the mutagenesis designs.
    1. Find the coordinates of the homology arms of the oligos.
    2. Construct a fasta file of the predicted mutant allele by removing the region between gRNA or by inserting the oligos sequence. For single gRNA designs create a 1 bp indel.

    Input: Reference genome
           gRNA Coordinates
           Oligos sequences?
           Output: Fasta file with mutation:
    1. Remove region between gRNAs
    2. Insert oligo sequence if available and produce and extra fasta

    '''

    chrom = grna[0][1].split(':')[0]
    cuts = inferCutCoord(grna)
    start = min(cuts) - 5000
    end = max(cuts) + 5000
    cuts = list(map(lambda x: x-start, cuts))
    fa = get_ref(ref, chrom, start, end)
    cutFa = trimFasta(fa, cuts)
    if donors:
        mutFa = inferDonors(cutFa, donors)
    else:
        mutFa = inferDel(cutFa)

    # TODO: Get smaller read length??
    # s = int((len(mutFa)/2) - 200)
    # e = int((len(mutFa)/2) + 200)

    # return(mutFa[s:e])
    return(mutFa)


def trimFasta(fasta, cuts):
    cutFa = []
    lf = list(fasta)
    for i in range(len(cuts)):
        dcut = cuts[i]
        ucut = cuts[i] + 1
        if i == 0 and i == len(cuts) - 1:
            dcut = cuts[i]
            ucut = cuts[i] + 1
            fasta_ls = ''.join(lf[:dcut])
            fasta_rs = ''.join(lf[ucut:])
            cutFa.append(fasta_ls)
            cutFa.append(fasta_rs)
        elif i == 0:
            ucut = cuts[i]
            fasta_ls = ''.join(lf[:ucut])
            cutFa.append(fasta_ls)
        elif i == len(cuts) - 1:
            dcut = cuts[i - 1] + 1
            ucut = cuts[i]
            fasta_ls = ''.join(lf[dcut:ucut])
            fasta_rs = ''.join(lf[ucut+1:])
            cutFa.append(fasta_ls)
            cutFa.append(fasta_rs)
        else:
            dcut = cuts[i - 1] + 1
            ucut = cuts[i]
            fasta_ls = ''.join(lf[dcut:ucut])
            cutFa.append(fasta_ls)
    return (cutFa)


def inferCutCoord(gRNA):
    cuts = []
    for rna in gRNA:
        seq = rna[0]
        coords = rna[1]
        ds = seq[-2:]
        us = seq[:2]
        sc = int(coords.split(':')[1].split('-')[0])
        ec = int(coords.split(':')[1].split('-')[1])
        if ds == 'GG':
            cut = ec - 4
        if us == 'CC':
            cut = sc + 4
        cuts.append(cut)
    return (cuts)


def inferDel(cutFa):
    return (cutFa[0] + cutFa[-1])


def inferDonors(cutFa, donors):
    if len(donors) > (len(cutFa)-1):
        exit('There are more donors than cuts in the DNA. Exiting.')

    if (len(cutFa)-1) == len(donors):
        j = 0
        mutFa = []
        prevBrkp = 0
        for donor in donors:
            fastaPaired = cutFa[j:j+2]
            for i in range(len(donor)):
                hau = donor[i:i+10]
                c = fastaPaired[0].find(hau)
                if c > 0:
                    ucf = c
                    ucd = donor.find(hau)
                    break
            mutFa.append(fastaPaired[0][prevBrkp:ucf])
            for i in range(len(donor)):
                had = donor[-(i+10):(len(donor)-i)]
                c = fastaPaired[1].find(had)
                if c > 0:
                    dcf = c + 10
                    dcd = donor.find(had) + 10
                    break
            mutFa.append(donor[ucd:dcd])
            prevBrkp = dcf
            j += 1
            if j == len(cutFa) - 1:
                mutFa.append(fastaPaired[1][dcf:])
        return (''.join(mutFa))

    if (len(cutFa)-1) == 2 and len(donors) == 1:
        del cutFa[1]
        return (inferDonors(cutFa, donors))

    if (len(cutFa)-1) > len(donors):
        mCutFa = []
        mCutFa.append(cutFa[0])
        mCutFa.append(''.join(cutFa[1:-1]))
        mCutFa.append(cutFa[-1])
        return (inferDonors(mCutFa, donors))
