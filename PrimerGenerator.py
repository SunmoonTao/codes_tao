
# coding: utf-8

# ###

# In[ ]:


from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from Bio import SeqIO
import primer3
from Bio.Alphabet import generic_dna,generic_rna,generic_protein
from Bio.Blast import NCBIWWW
import csv
from pydna.dseq import Dseq


# In[ ]:


import random
from Bio.Seq import Seq

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


# In[ ]:


# Return GC percentage*100
def gc_counter(mly_primer):
    gc=100*((mly_primer.count('G') + mly_primer.count('C'))/len(mly_primer))
    return gc


# In[ ]:


#### no more than 4 runs
def runs_counter(mly_primer):
    run_out =     mly_primer.count('GGGG') +     mly_primer.count('CCCC') +     mly_primer.count('AAAA') +     mly_primer.count('TTTT')
    if run_out > 0:
        return (False)
    else:
        return(True)


# In[ ]:


### no more than 5 di-repeats
def repeat_counter(mly_primer):
    repeats =     mly_primer.count('ATATATATAT') + mly_primer.count('ACACACACAC') + mly_primer.count('AGAGAGAGAG') +     mly_primer.count('TATATATATA') + mly_primer.count('TCTCTCTCTC') + mly_primer.count('TGTGTGTGTG') +     mly_primer.count('CACACACACA') + mly_primer.count('CGCGCGCGCG') + mly_primer.count('GAGAGAGAGA') +     mly_primer.count('GTGTGTGTGT') + mly_primer.count('GCGCGCGCGC')
    
    if repeats > 0:
        return (False)
    else:
        return(True)


# In[ ]:


### end stability
def end_3(mly_primer):
    xx= mly_primer[-5:].count('C') + mly_primer[-5:].count('G')
    e = False
    if mly_primer[-1]=='G' or mly_primer[-1]=='C':
        e = True

    if xx ==3 and e:
        return (True)
    
    
    


# In[ ]:


### generate unique list of primers 
### Added temp tuple 
def primer_generator(length,digestion_site, tests):
    mly_primer_20 = list()
    i = 0
    rc_digestion_site=str(Seq(digestion_site).reverse_complement())
    bp=length-5-len(digestion_site)
    while i <= tests:
        i=i+1
        mly_primer = str(RandomDNA_without_site(bp,digestion_site)) + digestion_site + str(RandomDNA_without_site(5,digestion_site))
        s=primer3.calcHairpin(mly_primer)
        if 53<primer3.calcTm(mly_primer)<55 and mly_primer.count(digestion_site) + mly_primer.count(rc_digestion_site) == 1        and not s.structure_found and 50 <= gc_counter(mly_primer) <= 60 and end_3(mly_primer)         and runs_counter(mly_primer) and repeat_counter(mly_primer):
            mly_primer_20.append(mly_primer)
    return (list(tuple(mly_primer_20)))     


# In[ ]:


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


# In[ ]:


### with function pop(0), I manage to generate unique pairs from one list
def uniqu_pair(mly_primer_20):
    mly_primer_pair=list()
    mly_primer_20_a=list(mly_primer_20)
    for i in mly_primer_20:
        mly_primer_20_a.pop(0)
        for j in mly_primer_20_a:
            mly_primer_pair.append ((i,j))
    return (mly_primer_pair)


# In[ ]:


# mly_primer_18 = list()
# i = 0
# while i < 1000000:
#     i=i+1
#     mly_primer = str(RandomDNA_without_site(8,'GAGTC')) + 'GAGTC' + str(RandomDNA_without_site(5,'GAGTC'))
#     s=primer3.calcHairpin(mly_primer)

#     if 53<primer3.calcTm(mly_primer)<55 and mly_primer.count('GAGTC') + mly_primer.count('GACTC') == 1\
#     and not s.structure_found and 50 <= gc_counter(mly_primer) <= 60 and end_3(mly_primer) \
#     and runs_counter(mly_primer) and repeat_counter(mly_primer):
#         mly_primer_18.append(mly_primer)
# len(mly_primer_18)
        
# for i in mly_primer_18:
#         print (i)


# In[ ]:



### Find top primer pair based on HeterodimerTm value(smallest the best)
def top_com(primers_list):
    temp_list=list()
    primers_list_a=primers_list[:]
    for i in primers_list:
        primers_list_a.pop(0)
        for j in primers_list_a:
            temp_list.append((primer3.calcHeterodimerTm(i,j),i,j))
    temp_list.sort(key=lambda x:x[0])
    top = temp_list[0]
    return (top)


# In[ ]:


### Find top primer pairs based on HeterodimerTm value(smallest the best), 
### Each primer can only be used in one pair of primers

def optimal_combination(primers_list):
    primer_candidates=primers_list[:]
    optimal_primers=list()
    while len(primer_candidates)>1:
        top=top_com(primer_candidates)
        primer_candidates.remove(top[1])
        primer_candidates.remove(top[2])
        optimal_primers.append(top)

    return(optimal_primers)            


# In[ ]:


def rc(str_seq):
    return str(Seq(str_seq).reverse_complement())


# ### Spacer  and MIPs

# In[ ]:


# OM6_12K = Mly1_F + rc(te.left.referencevalue.sequence) +'NNN' + rc(constant_adaptor_f[-27:]) + constant_adaptor_r[-27:] +'NNN'+ rc(te.right.referencevalue.sequence) + rc(Mly1_R)  


# In[ ]:


Mly1_F='TATGAGTGTGGAGTCGTTGC'
Mly1_R='GCTTCCTGATGAGTCCGATG'

constant_adaptor_r = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'  # 5->3
constant_adaptor_f = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'  # 5->3


# In[ ]:


om6_spacer = 'NNN' + str(Seq(constant_adaptor_f[-27:], generic_dna).reverse_complement()) + constant_adaptor_r[-27:] +'NNN'


# In[ ]:


om6_spacer


# In[ ]:


umi_12_spacer = 'NNNNNN' + str(Seq(constant_adaptor_f[-26:], generic_dna).reverse_complement()) + constant_adaptor_r[-26:] +'NNNNNN'


# In[ ]:


umi_12_spacer


# In[ ]:


##### OM6 Structure
# fw_watson='CCTGCTGAATTAGCCACCAAGTA'
# rev_watson='AACTTTTCAGAGGGAGCTTGCAA'
# ordered='TATGAGTGTGGAGTCGTTGCTACTTGGTGGCTAATTCAGCAGGNNNAGATCGGAAGAGCACACGTCTGAACTCTTTCCCTACACGACGCTCTTCCGATCTNNNTTGCAAGCTCCCTCTGAAAAGTTCATCGGACTCATCAGGAAGC'
# seq_ref ='CCTGCTGAATTAGCCACCAAGTACGCAAACTTTTCAGAGGGAGCTTGCAA'
# to_valid=Mly1_F + str(Seq(fw_watson).reverse_complement()) +'NNN' + str(Seq(constant_adaptor_f[-27:]).reverse_complement()) + constant_adaptor_r[-27:] +'NNN'+ str(Seq(rev_watson).reverse_complement()) + str(Seq(Mly1_R).reverse_complement())  
# to_valid==ordered


# In[ ]:


("ACGCGTTTGAGTCAGTGCTCACCATCATGTCACAAAGCACA")[::-1]


# In[ ]:


####control oligos
fw_watson='AGAATAGGCATCTGAGGACAGCC'
rev_watson='TCACCATCATGTCACAAAGCACA'
digestion_site = 'GAGTC'
rc_digestion_site=str(Seq(digestion_site).reverse_complement())

OM6_control_oligo=Mly1_F + str(Seq(fw_watson).reverse_complement()) +'NNN' + str(Seq(constant_adaptor_f[-27:]).reverse_complement()) + constant_adaptor_r[-27:] +'NNN'+ str(Seq(rev_watson).reverse_complement()) + str(Seq(Mly1_R).reverse_complement())  
OM6_control_probe=str(Seq(fw_watson).reverse_complement()) +'NNN' + str(Seq(constant_adaptor_f[-27:]).reverse_complement()) + constant_adaptor_r[-27:] +'NNN'+ str(Seq(rev_watson).reverse_complement())
for i in oc:
    Mly1_F=i[1]
    Mly1_R=i[2]
    Mly1_F_fwd = Mly1_F + str(Seq(fw_watson).reverse_complement())
    Mly1_R_rev = str(Seq(rev_watson).reverse_complement()) + str(Seq(Mly1_R).reverse_complement())

    oligo=Mly1_F + str(Seq(fw_watson).reverse_complement()) +'NNN' + str(Seq(constant_adaptor_f[-27:]).reverse_complement()) + constant_adaptor_r[-27:] +'NNN'+ str(Seq(rev_watson).reverse_complement()) + str(Seq(Mly1_R).reverse_complement())
#     print (Mly1_F_fwd)
#     print ('|'* len (Mly1_F_fwd))
#     print (oligo)
#     print (' '* (len(oligo)-len(Mly1_R_rev))+'|'* len (Mly1_R_rev))
#     print (' '* (len(oligo)-len(Mly1_R_rev))+Mly1_R_rev)
    print (Mly1_F_fwd, rc(Mly1_R_rev),Mly1_F_fwd.count(digestion_site) + Mly1_F_fwd.count(rc_digestion_site),Mly1_R_rev.count(digestion_site) + Mly1_R_rev.count(rc_digestion_site))


# In[ ]:


def digestion_sites_counter (seq, digestion_site):
    rc_digestion_site=str(Seq(digestion_site).reverse_complement())
    print (seq.count(digestion_site) + seq.count(rc_digestion_site))


# In[ ]:


from Bio.Restriction import *
MlyI.search(Seq(OM6_control_oligo))


# In[ ]:


# digestion_sites_counter(OM6_control_oligo,"GAGTC")
digestion_sites_counter(OM6_control_oligo,MlyI.site)


# In[ ]:


MlyI.catalyse(Seq(OM6_control_oligo))


# In[ ]:


from Bio.Restriction import *

def DigestionAnalysis(OM6_control_oligo,MlyI):
    rb=RestrictionBatch([MlyI])
    ana=Analysis(rb,Seq(OM6_control_oligo))
    # ana.print_as('number')
    # ana.print_as('alpha')
    ana.print_as('map')
    ana.print_that()
    for seq in MlyI.catalyse(Seq(OM6_control_oligo)):
        dsDNA(str(seq))


# In[ ]:


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


# In[ ]:


control_om_parts(fw_watson,rev_watson,Mly1_F,Mly1_R,constant_adaptor_f,constant_adaptor_r,UMI_length)


# In[ ]:


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


# In[ ]:


ill_sequence_spacer ='NNNNNN' + str(Seq(constant_adaptor_f, generic_dna).reverse_complement()) + constant_adaptor_r +'NNNNNN'


# In[ ]:


illu=(om_parts(fw_watson,rev_watson,Mly1_F,Mly1_R,constant_adaptor_f,constant_adaptor_r,UMI_length=12))


# In[ ]:


ill_sequence_spacer ==illu[2]


# In[ ]:


umi12=om_parts(fw_watson,rev_watson,Mly1_F,Mly1_R,constant_adaptor_f[-26:],constant_adaptor_r[-26:],12)


# In[ ]:


umi12_control_probe=umi12[1]


# In[ ]:


# illumina_p5 = "AATGATACGGCGACCACCGAGATCTACAC[UMI][i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
# illumina_p7 = "CAAGCAGAAGACGGCATACGAGAT[UMI][i7]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
illumina_p5 = "AATGATACGGCGACCACCGAGATCTACAC"+'N'*8 +'AAAAAAAA'+"ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
illumina_p7 = "CAAGCAGAAGACGGCATACGAGAT"+'N'*8 +'TTTTTTTT'+"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
whole_ill_spacer =str(Seq(illumina_p5, generic_dna).reverse_complement()) + illumina_p7


# In[ ]:


illumina_p5 = "AATGATACGGCGACCACCGAGATCTACAC"
illumina_p7 = "CAAGCAGAAGACGGCATACGAGAT"


# In[ ]:


whole_ill_spacer


# In[ ]:


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
        
s=dsDNA(whole_ill_spacer)


# In[ ]:


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


# In[ ]:


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

# test_1313231="TATGAGTGTGGAGTCGTTGCGGCTGTCCTCAGATGCCTATTCTNNNAGATCGGAAGAGCACACGTCTGAACTCTTTCCCTACACGACGCTCTTCCGATCTNNNTGTGCTTTGTGACATGATGGTGACATCGGACTCATCAGGAAGC"
# break2shorts(test_1313231)


# In[ ]:


break2shorts(OM6_control_probe,50)


# In[ ]:


break2shorts(umi12_control_probe)


# In[ ]:


###validate with snapgene
OM6_control_oligo=="TATGAGTGTGGAGTCGTTGCGGCTGTCCTCAGATGCCTATTCTNNNAGATCGGAAGAGCACACGTCTGAACTCTTTCCCTACACGACGCTCTTCCGATCTNNNTGTGCTTTGTGACATGATGGTGACATCGGACTCATCAGGAAGC"


# ## Blast online

# In[ ]:


#cd OMs


# In[ ]:


from Bio.Blast import NCBIXML


# In[ ]:


seq=Seq(fw_watson,generic_dna)
fasta_string = open("seq.fasta").read()


# In[ ]:


with open('seq.fasta','w') as f:
    f.write('>\r')
    f.write(str(seq))
f.close()   
    


# In[ ]:


with open('seq.fasta','r') as f:
    s=f.read()
    print (s)
f.close()


# In[ ]:


seq


# ## Available database
# Name	Title	Type
# nt	Nucleotide collection	DNA
# nr	Non-redundant	Protein
# refseq_rna	NCBI Transcript Reference Sequences	DNA
# refseq_protein	NCBI Protein Reference Sequences	Protein
# swissprot	Non-redundant UniProtKB/SwissProt sequences	Protein
# pdbaa	PDB protein database	Protein
# pdbnt	PDB nucleotide database	DNA

# In[ ]:


result_handle = NCBIWWW.qblast("blastn", "nt", str(seq))
# result_handle = NCBIWWW.qblast("blastn", "nt", "test_5.fasta", megablast=True)

br = NCBIXML.read(result_handle)


# In[ ]:


E_VALUE_THRESH = 0.001
# for m in br:
for alignment in br.alignments:
    for hsp in alignment.hsps:
        print("sequence:", alignment.title)
        print("length:", alignment.length)
        print("e value:", hsp.expect)


# In[ ]:


br.alignments


# ## Blast locally

# In[ ]:


from Bio.Blast.Applications import NcbiblastnCommandline


# In[ ]:


cline = NcbiblastnCommandline(query="/Users/ltao/gitsss/biopython/Doc/examples/m_cold.fasta", db="/Users/ltao/ref/blastdb/nt",
                              evalue=0.001, out="m_cold.xml", outfmt=5)
cline
print(cline)


# In[ ]:


stdout, stderr = blastn_cline()


# ### Parsing NCBIXML

# In[ ]:


from Bio.Blast import NCBIXML

result_handle = open("m_cold.xml")
blast_record = NCBIXML.read(result_handle)


# In[ ]:


blastn_cline = NcbiblastnCommandline(query='seq.fasta', db="/Users/ltao/ref/blastdb/nt", evalue=0.1,
                                     outfmt=5, out="seq_blstn.xml", )
blastn_cline
print(blastn_cline)
stdout, stderr = blastn_cline()


# In[ ]:


from Bio.Blast import NCBIXML

result_handle = open("seq_blstn.xml")
blast_record = NCBIXML.read(result_handle)


# In[ ]:


blast_record.alignments

