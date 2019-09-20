from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from Bio import SeqIO
import primer3
from Bio.Alphabet import generic_dna,generic_rna,generic_protein
# from Bio.Blast import NCBIWWW
import csv
# from pydna.dseq import Dseq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import random


def RandomDNA(length):
    str_dna= ''.join(random.choice("A"*25+"C"*25+"G"*25+"T"*25) for _ in range(length))
    seq = Seq(str_dna, generic_dna)
    return seq
    
def RandomDNA_without_site(length, site_to_ruleout):
    handler = 0
    while handler == 0:
        candi=RandomDNA(length)
        if (site_to_ruleout not in candi) and (site_to_ruleout not in candi.reverse_complement()):
            handle = 1
            return candi
def RandomDNA_without_sites(length, sites_list_to_ruleout=[]):
    handler = 0
    while handler == 0:
        candi = RandomDNA(length)
        if any(site in candi for site in sites_list_to_ruleout) or any(site in candi.reverse_complement() for site in sites_list_to_ruleout):
            handler = 0
        else:
            handler = 1
            return candi


def rc(seq):
    return( str(Seq(seq).reverse_complement()))

# Return GC percentage*100
def gc_counter(mly_primer):
    gc=100*((mly_primer.count('G') + mly_primer.count('C'))/len(mly_primer))
    return gc
def gc(s):
    return((s.count('G')+s.count('C'))/len(s))


#### no more than 4 runs
def runs_counter(mly_primer):
    run_out =     mly_primer.count('GGGG') +     mly_primer.count('CCCC') +     mly_primer.count('AAAA') +     mly_primer.count('TTTT')
    if run_out > 0:
        return (False)
    else:
        return(True)


### no more than 5 di-repeats
def repeat_counter(mly_primer):
    repeats =     mly_primer.count('ATATATATAT') + mly_primer.count('ACACACACAC') + mly_primer.count('AGAGAGAGAG') +     mly_primer.count('TATATATATA') + mly_primer.count('TCTCTCTCTC') + mly_primer.count('TGTGTGTGTG') +     mly_primer.count('CACACACACA') + mly_primer.count('CGCGCGCGCG') + mly_primer.count('GAGAGAGAGA') +     mly_primer.count('GTGTGTGTGT') + mly_primer.count('GCGCGCGCGC')
    
    if repeats > 0:
        return (False)
    else:
        return(True)


### end stability
def end_3(mly_primer):
    xx= mly_primer[-5:].count('C') + mly_primer[-5:].count('G')
    e = False
    if mly_primer[-1]=='G' or mly_primer[-1]=='C':
        e = True

    if xx ==3 and e:
        return (True)
    
    


### generate unique list of primers 
### Added temp tuple 
### Add end3_off option
def primer_generator(length,digestion_site, tests,end_CG=True):
    mly_primer_20 = list()
    i = 0
    rc_digestion_site=str(Seq(digestion_site).reverse_complement())
    bp=length-5-len(digestion_site)
    while i <= tests:
        i=i+1
        mly_primer = str(RandomDNA_without_site(bp,digestion_site)) + digestion_site + str(RandomDNA_without_site(5,digestion_site))
        s=primer3.calcHairpin(mly_primer)
        if end_CG:
            if 53<primer3.calcTm(mly_primer)<55 and mly_primer.count(digestion_site) + mly_primer.count(rc_digestion_site) == 1\
            and not s.structure_found and 50 <= gc_counter(mly_primer) <= 60 and end_3(mly_primer) \
            and runs_counter(mly_primer) and repeat_counter(mly_primer):
                mly_primer_20.append(mly_primer)
        else:
            if 53<primer3.calcTm(mly_primer)<55 and mly_primer.count(digestion_site) + mly_primer.count(rc_digestion_site) == 1\
            and not s.structure_found and 50 <= gc_counter(mly_primer) <= 60 \
            and runs_counter(mly_primer) and repeat_counter(mly_primer):
                mly_primer_20.append(mly_primer)           
    return (list(tuple(mly_primer_20)))  


### with function pop(0), I manage to generate unique pairs from one list passing HETERODIMERTM
def primer_pair_generator(mly_primer_20):
    mly_primer_pair=list()
    mly_primer_20_a=list(mly_primer_20)
    for i in mly_primer_20:
        mly_primer_20_a.pop(0)
        for j in mly_primer_20_a:
            if primer3.calcHeterodimerTm(i,j)< 0:
                mly_primer_pair.append ((i,j))
    return mly_primer_pair


### with function pop(0), I manage to generate unique pairs from one list
def uniqu_pair(mly_primer_20):
    mly_primer_pair=list()
    mly_primer_20_a=list(mly_primer_20)
    for i in mly_primer_20:
        mly_primer_20_a.pop(0)
        for j in mly_primer_20_a:
            mly_primer_pair.append ((i,j))
    return (mly_primer_pair)


### Find top primer pairs based on HeterodimerTm value(smallest the best), 
### Each primer can only be used in one pair of primers
### similarity between pairs: use score of alignment to tell apart

def optimal_combination(primers_list):
    primer_candidates=primers_list[:]
    optimal_primers=list()
    while len(primer_candidates)>1:
        top=top_com(primer_candidates)
        primer_candidates.remove(top[1])
        primer_candidates.remove(top[2])
        
        alignments = pairwise2.align.globalms(top[1], top[2],1, -1, -0.5, -0.1)
        alignments_end3 = pairwise2.align.globalms(top[1][-5:], top[2][-5:],1, -1, -0.5, -0.1)

        if alignments_end3[0][2] < 3 and alignments[0][2] <10 :
            top.append(alignments[0])
            top.append(alignments_end3[0])
            optimal_primers.append(top)
#         optimal_primers.append(top)

    return(optimal_primers)            

## yield will create a generator instead of a list
def substring_indexes(substring, string):
    """ 
    Generate indices of where substring begins in string

    >>> list(find_substring('me', "The cat says meow, meow"))
    [13, 19]
    """
    last_found = -1  # Begin at -1 so the next position to search from is 0
    while True:
        # Find next index of substring, by starting after its last known position
        last_found = string.find(substring, last_found + 1)
        if last_found == -1:  
            break  # All occurrences have been found
        yield last_found
# ### Spacer  and MIPs



Mly1_F='TATGAGTGTGGAGTCGTTGC'
Mly1_R='GCTTCCTGATGAGTCCGATG'

constant_adaptor_r = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'  # 5->3
constant_adaptor_f = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'  # 5->3




om6_spacer = 'NNN' + str(Seq(constant_adaptor_f[-27:], generic_dna).reverse_complement()) + constant_adaptor_r[-27:] +'NNN'




umi_12_spacer = 'NNNNNN' + str(Seq(constant_adaptor_f[-26:], generic_dna).reverse_complement()) + constant_adaptor_r[-26:] +'NNNNNN'



##### OM6 Structure
# fw_watson='CCTGCTGAATTAGCCACCAAGTA'
# rev_watson='AACTTTTCAGAGGGAGCTTGCAA'
# ordered='TATGAGTGTGGAGTCGTTGCTACTTGGTGGCTAATTCAGCAGGNNNAGATCGGAAGAGCACACGTCTGAACTCTTTCCCTACACGACGCTCTTCCGATCTNNNTTGCAAGCTCCCTCTGAAAAGTTCATCGGACTCATCAGGAAGC'
# seq_ref ='CCTGCTGAATTAGCCACCAAGTACGCAAACTTTTCAGAGGGAGCTTGCAA'
# to_valid=Mly1_F + str(Seq(fw_watson).reverse_complement()) +'NNN' + str(Seq(constant_adaptor_f[-27:]).reverse_complement()) + constant_adaptor_r[-27:] +'NNN'+ str(Seq(rev_watson).reverse_complement()) + str(Seq(Mly1_R).reverse_complement())  
# to_valid==ordered





("ACGCGTTTGAGTCAGTGCTCACCATCATGTCACAAAGCACA")[::-1]





####control oligos
fw_watson='AGAATAGGCATCTGAGGACAGCC'
rev_watson='TCACCATCATGTCACAAAGCACA'
digestion_site = 'GAGTC'
rc_digestion_site=str(Seq(digestion_site).reverse_complement())





def digestion_sites_counter (seq, digestion_site):
    rc_digestion_site=str(Seq(digestion_site).reverse_complement())
    print (seq.count(digestion_site) + seq.count(rc_digestion_site))




def DigestionAnalysis(OM6_control_oligo,MlyI):
    rb=RestrictionBatch([MlyI])
    ana=Analysis(rb,Seq(OM6_control_oligo))
    # ana.print_as('number')
    # ana.print_as('alpha')
    ana.print_as('map')
    ana.print_that()
    for seq in MlyI.catalyse(Seq(OM6_control_oligo)):
        dsDNA(str(seq))




####control oligos

fwd_watson='AGAATAGGCATCTGAGGACAGCC'
rev_crick='TCACCATCATGTCACAAAGCACA'
Mly1_F='TATGAGTGTGGAGTCGTTGC'
Mly1_R='GCTTCCTGATGAGTCCGATG'
constant_adaptor_f
constant_adaptor_r
UMI_length=6
def control_om_parts(fwd_watson,rev_crick,Mly1_F,Mly1_R,constant_adaptor_f,constant_adaptor_r,UMI_length):

    OM6_control_oligo=Mly1_F + str(Seq(fwd_watson).reverse_complement()) +'N'* (int(UMI_length/2)) + str(Seq(constant_adaptor_f[-27:]).reverse_complement()) + constant_adaptor_r[-27:] +'N'* (int(UMI_length/2))+ str(Seq(rev_crick).reverse_complement()) + str(Seq(Mly1_R).reverse_complement())  
    OM6_control_probe=str(Seq(fwd_watson).reverse_complement()) +'N'* (int(UMI_length/2)) + str(Seq(constant_adaptor_f[-27:]).reverse_complement()) + constant_adaptor_r[-27:] +'N'* (int(UMI_length/2))+ str(Seq(rev_crick).reverse_complement())
    OM6_control_spacer='N'* (int(UMI_length/2)) + str(Seq(constant_adaptor_f[-27:]).reverse_complement()) + constant_adaptor_r[-27:] +'N'* (int(UMI_length/2))
    Mly1_F_fwd = Mly1_F + str(Seq(fwd_watson).reverse_complement())
    Mly1_R_rev = str(Seq(rev_crick).reverse_complement()) + str(Seq(Mly1_R).reverse_complement())
    oligo=Mly1_F + str(Seq(fwd_watson).reverse_complement()) +'N'* (int(UMI_length/2)) + str(Seq(constant_adaptor_f[-27:]).reverse_complement()) + constant_adaptor_r[-27:] +'N'* (int(UMI_length/2))+ str(Seq(rev_crick).reverse_complement()) + str(Seq(Mly1_R).reverse_complement())
#     print (Mly1_F_fwd)
#     print ('|'* len (Mly1_F_fwd))
#     print (oligo)
#     print (' '* (len(oligo)-len(Mly1_R_rev))+'|'* len (Mly1_R_rev))
#     print (' '* (len(oligo)-len(Mly1_R_rev))+Mly1_R_rev)
    return (OM6_control_oligo,OM6_control_probe,OM6_control_spacer,Mly1_F_fwd, rc(Mly1_R_rev))




#### om_parts 

def om_parts(fw_watson,rev_watson,Mly1_F,Mly1_R,constant_adaptor_f,constant_adaptor_r,UMI_length):

    OM6_control_oligo=Mly1_F + str(Seq(fw_watson).reverse_complement()) +'N'* (int(UMI_length/2)) + str(Seq(constant_adaptor_f).reverse_complement()) + constant_adaptor_r+'N'* (int(UMI_length/2))+ str(Seq(rev_watson).reverse_complement()) + str(Seq(Mly1_R).reverse_complement())  
    OM6_control_probe=str(Seq(fw_watson).reverse_complement()) +'N'* (int(UMI_length/2)) + str(Seq(constant_adaptor_f).reverse_complement()) + constant_adaptor_r +'N'* (int(UMI_length/2))+ str(Seq(rev_watson).reverse_complement())
    OM6_control_spacer='N'* (int(UMI_length/2)) + str(Seq(constant_adaptor_f).reverse_complement()) + constant_adaptor_r +'N'* (int(UMI_length/2))
    Mly1_F_fwd = Mly1_F + str(Seq(fw_watson).reverse_complement())
    Mly1_R_rev = str(Seq(rev_watson).reverse_complement()) + str(Seq(Mly1_R).reverse_complement())
    oligo=Mly1_F + str(Seq(fw_watson).reverse_complement()) +'N'* (int(UMI_length/2)) + str(Seq(constant_adaptor_f).reverse_complement()) + constant_adaptor_r +'N'* (int(UMI_length/2))+ str(Seq(rev_watson).reverse_complement()) + str(Seq(Mly1_R).reverse_complement())
#     print (Mly1_F_fwd)
#     print ('|'* len (Mly1_F_fwd))
#     print (oligo)
#     print (' '* (len(oligo)-len(Mly1_R_rev))+'|'* len (Mly1_R_rev))
#     print (' '* (len(oligo)-len(Mly1_R_rev))+Mly1_R_rev)
    return (OM6_control_oligo,OM6_control_probe,OM6_control_spacer,Mly1_F_fwd, rc(Mly1_R_rev))



ill_sequence_spacer ='NNNNNN' + str(Seq(constant_adaptor_f, generic_dna).reverse_complement()) + constant_adaptor_r +'NNNNNN'


# illumina_p5 = "AATGATACGGCGACCACCGAGATCTACAC[UMI][i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
# illumina_p7 = "CAAGCAGAAGACGGCATACGAGAT[UMI][i7]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
illumina_p5 = "AATGATACGGCGACCACCGAGATCTACAC"+'N'*8 +'AAAAAAAA'+"ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
illumina_p7 = "CAAGCAGAAGACGGCATACGAGAT"+'N'*8 +'TTTTTTTT'+"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
whole_ill_spacer =str(Seq(illumina_p5, generic_dna).reverse_complement()) + illumina_p7




illumina_p5 = "AATGATACGGCGACCACCGAGATCTACAC"
illumina_p7 = "CAAGCAGAAGACGGCATACGAGAT"


def dsDNA(oligo):
    print ("dsDNA style:")
    while len(oligo) > 50:
        piece=oligo[:50]
        print ('5-'+piece+'-3')
        print (' '*2 +'|'*len(piece)+' '*2 )
        print ('3-'+rc(piece)[::-1]+'-5')
        print ()
        oligo=oligo[50:]
    print ('5-'+oligo+'-3')
    print (' '*2 +'|'*len(oligo)+' '*2 )
    print ('3-'+rc(oligo)[::-1]+'-5')
    print ()


# def dsDNA(oligo):
#     output=list()
#     output.append ("dsDNA style:")
#     while len(oligo) > 50:
#         piece=oligo[:50]
#         ws='5-'+piece+'-3'
#         output.append (ws)
#         bonds=' '*2 +'|'*len(piece)+' '*2
#         output.append (bonds)
#         cs='3-'+rc(piece)[::-1]+'-5'
#         output.append (cs)
#         output.append ('\n')
#         oligo=oligo[50:]
#     ws='5-'+oligo+'-3'
#     output.append (ws)
#     bonds=' '*2 +'|'*len(oligo)+' '*2
#     output.append (bonds)
#     cs='3-'+rc(oligo)[::-1]+'-5'
#     output.append (cs)
#     output.append ('\n')
#     return output  
# s=dsDNA(whole_ill_spacer)


### find a error when break into two parts, the direction of secondpart should only be reversed, instead of rc
###Define every concept in a clear seperate variable
def break2shorts(input_long_oligo, overlap=30):
#     print ('5-'+input_long_oligo+'-3')
#     print ('|'*(len(input_long_oligo)+4))
#     print ('3-'+rc(input_long_oligo)[::-1]+'-5')

    print ("original long:",len(input_long_oligo))
    dsDNA(input_long_oligo)
    
    waston = input_long_oligo
    crick = rc(input_long_oligo)  ### should not use 3-->5 define a DNA seq
    
    l= int((len(input_long_oligo)+overlap)/2)
    
    print ("\n two shorts:")
    print ('5-'+waston[:l]+'-3')
    print ('|'*(len(input_long_oligo)+4))
    print (' '*(l-overlap)+ '3-'+ crick[::-1][l-overlap:]+'-5')

    part_A= waston[:l]
    Part_B= crick[::-1][l-overlap:][::-1]    
    return (part_A,len(part_A),Part_B,len(Part_B))

# 