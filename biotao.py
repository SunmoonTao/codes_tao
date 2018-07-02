
# coding: utf-8

# In[2]:


from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from Bio import SeqIO
import primer3
from Bio.Alphabet import generic_dna,generic_rna,generic_protein
from Bio.Restriction import *


# In[3]:


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


# In[4]:




# In[5]:




# In[6]:


# Return GC percentage*100
def gc_counter(mly_primer):
    gc=100*((mly_primer.count('G') + mly_primer.count('C'))/len(mly_primer))
    return gc


# In[7]:


#### no more than 4 runs
def runs_counter(mly_primer):
    run_out =     mly_primer.count('GGGG') +     mly_primer.count('CCCC') +     mly_primer.count('AAAA') +     mly_primer.count('TTTT')
    if run_out > 0:
        return (False)
    else:
        return(True)


# In[8]:


### no more than 5 di-repeats
def repeat_counter(mly_primer):
    repeats =     mly_primer.count('ATATATATAT') + mly_primer.count('ACACACACAC') + mly_primer.count('AGAGAGAGAG') +     mly_primer.count('TATATATATA') + mly_primer.count('TCTCTCTCTC') + mly_primer.count('TGTGTGTGTG') +     mly_primer.count('CACACACACA') + mly_primer.count('CGCGCGCGCG') + mly_primer.count('GAGAGAGAGA') +     mly_primer.count('GTGTGTGTGT') + mly_primer.count('GCGCGCGCGC')
    
    if repeats > 0:
        return (False)
    else:
        return(True)


# In[9]:


### end stability
def end_3(mly_primer):
    xx= mly_primer[-5:].count('C') + mly_primer[-5:].count('G')
    e = False
    if mly_primer[-1]=='G' or mly_primer[-1]=='C':
        e = True

    if xx ==3 and e:
        return (True)
    
    
    


# In[10]:


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


# In[11]:


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


# In[12]:


### with function pop(0), I manage to generate unique pairs from one list
def uniqu_pair(mly_primer_20):
    mly_primer_pair=list()
    mly_primer_20_a=list(mly_primer_20)
    for i in mly_primer_20:
        mly_primer_20_a.pop(0)
        for j in mly_primer_20_a:
            mly_primer_pair.append ((i,j))
    return (mly_primer_pair)


# In[13]:


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


# In[14]:



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


# In[15]:


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


# In[37]:




# In[42]:





# In[18]:


def rc(str_seq):
    return str(Seq(str_seq).reverse_complement())


# ### Spacer  and MIPs

# In[19]:


# OM6_12K = Mly1_F + rc(te.left.referencevalue.sequence) +'NNN' + rc(constant_adaptor_f[-27:]) + constant_adaptor_r[-27:] +'NNN'+ rc(te.right.referencevalue.sequence) + rc(Mly1_R)  


# In[20]:


#str(my_dna)[::-1]


# In[21]:




# In[22]:


Mly1_F='TATGAGTGTGGAGTCGTTGC'
Mly1_R='GCTTCCTGATGAGTCCGATG'

constant_adaptor_r = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'  # 5->3
constant_adaptor_f = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'  # 5->3
###use 18 bp adaptors; use 26bp illumina_spacer; UMI = 12bp on both side.


# In[23]:


om6_spacer = 'NNN' + str(Seq(constant_adaptor_f[-27:], generic_dna).reverse_complement()) + constant_adaptor_r[-27:] +'NNN'


# In[49]:





# In[24]:


umi_12_spacer = 'NNNNNN' + str(Seq(constant_adaptor_f[-26:], generic_dna).reverse_complement()) + constant_adaptor_r[-26:] +'NNNNNN'


# In[50]:




# In[25]:


##### OM6 Structure
# fw_ref='CCTGCTGAATTAGCCACCAAGTA'
# rev_ref='AACTTTTCAGAGGGAGCTTGCAA'
# ordered='TATGAGTGTGGAGTCGTTGCTACTTGGTGGCTAATTCAGCAGGNNNAGATCGGAAGAGCACACGTCTGAACTCTTTCCCTACACGACGCTCTTCCGATCTNNNTTGCAAGCTCCCTCTGAAAAGTTCATCGGACTCATCAGGAAGC'
# seq_ref ='CCTGCTGAATTAGCCACCAAGTACGCAAACTTTTCAGAGGGAGCTTGCAA'
# to_valid=Mly1_F + str(Seq(fw_ref).reverse_complement()) +'NNN' + str(Seq(constant_adaptor_f[-27:]).reverse_complement()) + constant_adaptor_r[-27:] +'NNN'+ str(Seq(rev_ref).reverse_complement()) + str(Seq(Mly1_R).reverse_complement())  
# to_valid==ordered


# In[115]:


####control oligos
fw_ref='AGAATAGGCATCTGAGGACAGCC'
rev_ref='TCACCATCATGTCACAAAGCACA'
digestion_site = 'GAGTC'
rc_digestion_site=str(Seq(digestion_site).reverse_complement())

OM6_control_oligo=Mly1_F + str(Seq(fw_ref).reverse_complement()) +'NNN' + str(Seq(constant_adaptor_f[-27:]).reverse_complement()) + constant_adaptor_r[-27:] +'NNN'+ str(Seq(rev_ref).reverse_complement()) + str(Seq(Mly1_R).reverse_complement())  
OM6_control_probe=str(Seq(fw_ref).reverse_complement()) +'NNN' + str(Seq(constant_adaptor_f[-27:]).reverse_complement()) + constant_adaptor_r[-27:] +'NNN'+ str(Seq(rev_ref).reverse_complement())

# In[127]:


def digestion_sites_counter (seq, digestion_site):
    rc_digestion_site=str(Seq(digestion_site).reverse_complement())
    print (seq.count(digestion_site) + seq.count(rc_digestion_site))


# In[132]:


# digestion_sites_counter(OM6_control_oligo,"GAGTC")


# In[136]:




# In[137]:




# In[234]:



def DigestionAnalysis(OM6_control_oligo,MlyI):
    rb=RestrictionBatch([MlyI])
    ana=Analysis(rb,Seq(OM6_control_oligo))
    # ana.print_as('number')
    # ana.print_as('alpha')
    ana.print_as('map')
    ana.print_that()
    for seq in MlyI.catalyse(Seq(OM6_control_oligo)):
        dsDNA(str(seq))


# In[235]:




# In[107]:


####control oligos

fw_ref='AGAATAGGCATCTGAGGACAGCC'
rev_ref='TCACCATCATGTCACAAAGCACA'
Mly1_F
Mly1_R
constant_adaptor_f
constant_adaptor_r
UMI_length=6
def control_om_parts(fw_ref,rev_ref,Mly1_F,Mly1_R,constant_adaptor_f,constant_adaptor_r,UMI_length):

    OM6_control_oligo=Mly1_F + str(Seq(fw_ref).reverse_complement()) +'N'* (int(UMI_length/2)) + str(Seq(constant_adaptor_f[-27:]).reverse_complement()) + constant_adaptor_r[-27:] +'N'* (int(UMI_length/2))+ str(Seq(rev_ref).reverse_complement()) + str(Seq(Mly1_R).reverse_complement())  
    OM6_control_probe=str(Seq(fw_ref).reverse_complement()) +'N'* (int(UMI_length/2)) + str(Seq(constant_adaptor_f[-27:]).reverse_complement()) + constant_adaptor_r[-27:] +'N'* (int(UMI_length/2))+ str(Seq(rev_ref).reverse_complement())
    OM6_control_spacer='N'* (int(UMI_length/2)) + str(Seq(constant_adaptor_f[-27:]).reverse_complement()) + constant_adaptor_r[-27:] +'N'* (int(UMI_length/2))
    Mly1_F_fwd = Mly1_F + str(Seq(fw_ref).reverse_complement())
    Mly1_R_rev = str(Seq(rev_ref).reverse_complement()) + str(Seq(Mly1_R).reverse_complement())
    oligo=Mly1_F + str(Seq(fw_ref).reverse_complement()) +'N'* (int(UMI_length/2)) + str(Seq(constant_adaptor_f[-27:]).reverse_complement()) + constant_adaptor_r[-27:] +'N'* (int(UMI_length/2))+ str(Seq(rev_ref).reverse_complement()) + str(Seq(Mly1_R).reverse_complement())
#     print (Mly1_F_fwd)
#     print ('|'* len (Mly1_F_fwd))
#     print (oligo)
#     print (' '* (len(oligo)-len(Mly1_R_rev))+'|'* len (Mly1_R_rev))
#     print (' '* (len(oligo)-len(Mly1_R_rev))+Mly1_R_rev)
    return (OM6_control_oligo,OM6_control_probe,OM6_control_spacer,Mly1_F_fwd, rc(Mly1_R_rev))


# In[108]:


control_om_parts(fw_ref,rev_ref,Mly1_F,Mly1_R,constant_adaptor_f,constant_adaptor_r,UMI_length)


# In[46]:


ill_sequence_spacer ='NNNNNN' + str(Seq(constant_adaptor_f, generic_dna).reverse_complement()) + constant_adaptor_r +'NNNNNN'


# In[28]:


# illumina_p5 = "AATGATACGGCGACCACCGAGATCTACAC[UMI][i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
# illumina_p7 = "CAAGCAGAAGACGGCATACGAGAT[UMI][i7]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
illumina_p5 = "AATGATACGGCGACCACCGAGATCTACAC"+'N'*8 +'AAAAAAAA'+"ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
illumina_p7 = "CAAGCAGAAGACGGCATACGAGAT"+'N'*8 +'TTTTTTTT'+"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
whole_ill_spacer =str(Seq(illumina_p5, generic_dna).reverse_complement()) + illumina_p7


# In[39]:


illumina_p5 = "AATGATACGGCGACCACCGAGATCTACAC"
illumina_p7 = "CAAGCAGAAGACGGCATACGAGAT"


# In[51]:





# In[247]:


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
        


# In[248]:





# In[242]:


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


# In[245]:





# In[90]:


def break2shorts(input_long_oligo, overlap=30):
#     print ('5-'+input_long_oligo+'-3')
#     print ('|'*(len(input_long_oligo)+4))
#     print ('3-'+rc(input_long_oligo)[::-1]+'-5')

    print ("original long:",len(input_long_oligo))
    dsDNA(input_long_oligo)
    
    fwd = input_long_oligo
    rev = rc(input_long_oligo)[::-1]
    
    l= int((len(input_long_oligo)+overlap)/2)
    
    print ("\n two shorts:")
    print ('5-'+fwd[:l]+'-3')
    print ('|'*(len(input_long_oligo)+4))
    print (' '*(l-overlap)+ '3-'+ rev[l-overlap:]+'-5')
    
    return (fwd[:l],len(fwd[:l]),rc(rev[l-overlap:]),len(rc(rev[l-overlap:])))


# In[214]:




# In[31]:




# In[32]:




