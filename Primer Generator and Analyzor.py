
# coding: utf-8

# In[3]:


from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from Bio import SeqIO
import primer3
from Bio.Alphabet import generic_dna,generic_rna,generic_protein


# In[4]:


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


# In[5]:


primer3.calcHomodimer("TAGAGATAGTGTGGTGGCTCA")


# In[6]:


primer3.calcHeterodimer("TAGAGATAGTGTGGTGGCTCA", str(Seq("TAGAGATAGTGTGGTGGCTCA").reverse_complement()))


# In[7]:


# Return GC percentage*100
def gc_counter(mly_primer):
    gc=100*((mly_primer.count('G') + mly_primer.count('C'))/len(mly_primer))
    return gc


# In[8]:


#### no more than 4 runs
def runs_counter(mly_primer):
    run_out =     mly_primer.count('GGGG') +     mly_primer.count('CCCC') +     mly_primer.count('AAAA') +     mly_primer.count('TTTT')
    if run_out > 0:
        return (False)
    else:
        return(True)


# In[9]:


### no more than 5 di-repeats
def repeat_counter(mly_primer):
    repeats =     mly_primer.count('ATATATATAT') + mly_primer.count('ACACACACAC') + mly_primer.count('AGAGAGAGAG') +     mly_primer.count('TATATATATA') + mly_primer.count('TCTCTCTCTC') + mly_primer.count('TGTGTGTGTG') +     mly_primer.count('CACACACACA') + mly_primer.count('CGCGCGCGCG') + mly_primer.count('GAGAGAGAGA') +     mly_primer.count('GTGTGTGTGT') + mly_primer.count('GCGCGCGCGC')
    
    if repeats > 0:
        return (False)
    else:
        return(True)


# In[10]:


### end stability
def end_3(mly_primer):
    xx= mly_primer[-5:].count('C') + mly_primer[-5:].count('G')
    e = False
    if mly_primer[-1]=='G' or mly_primer[-1]=='C':
        e = True

    if xx ==3 and e:
        return (True)
    
    
    


# In[90]:


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


# In[89]:


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


# In[14]:


### with function pop(0), I manage to generate unique pairs from one list
def uniqu_pair(mly_primer_20):
    mly_primer_pair=list()
    mly_primer_20_a=list(mly_primer_20)
    for i in mly_primer_20:
        mly_primer_20_a.pop(0)
        for j in mly_primer_20_a:
            mly_primer_pair.append ((i,j))
    return (mly_primer_pair)


# In[26]:


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


# In[77]:



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


# In[82]:


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


# In[94]:


sss=primer_generator(22, 'GAGTC', 5000)
print (len(sss))


# In[95]:


oc=optimal_combination(sss)
len(oc)


# In[96]:


oc


# In[15]:


str(my_dna)[::-1]


# In[ ]:


str_rna=''


# In[ ]:


str_protein=''


# In[5]:


str(my_dna.complement())


# In[4]:


str(my_dna.reverse_complement())


# In[146]:


my_dna.find("ACT")


# In[149]:


my_dna.count("AT")


# In[ ]:


my_dna.transcribe()


# In[ ]:


my_rna = my_dna.transcribe() 
my_rna

my_rna.back_transcribe()


# In[ ]:


messenger_rna = my_rna
messenger_rna.translate()


# In[ ]:


coding_dna=my_dna
coding_dna.translate(table=2, to_stop=True)


# In[12]:


import sys
sys.path.append("/home/liming/s/LimingTao/Bioinformatics/code/clineage/")
import clineage.wsgi
from utils.SequenceManipulations import randseq


# In[14]:


randseq(15)


# In[66]:


Mly1_F='TATGAGTGTGGAGTCGTTGC'
Mly1_R='GCTTCCTGATGAGTCCGATG'
Mly_F_a="TGCGCATGGAGTCTCACC"
Mly_R_a="TGCTTGCGGAGTCAATGC"
Mly_F_g="TGCATGCAGAGTCAACGC"
Mly_R_g="TCATTCGCGAGTCGCAAC"
constant_adaptor_r = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'  # 5->3
constant_adaptor_f = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'  # 5->3
###use 18 bp adaptors; use 26bp illumina_spacer; UMI = 12bp on both side.
sss='NNNNNN' + str(Seq(constant_adaptor_f, generic_dna).reverse_complement()) + constant_adaptor_r +'NNNNNN'


# In[69]:


om6_spacer = 'NNN' + str(Seq(constant_adaptor_f[-27:], generic_dna).reverse_complement()) + constant_adaptor_r[-27:] +'NNN'


# In[72]:


umi_12_spacer = 'NNNNNN' + str(Seq(constant_adaptor_f[-26:], generic_dna).reverse_complement()) + constant_adaptor_r[-26:] +'NNNNNN'


# In[70]:


ill_sequence_spacer ='NNNNNN' + str(Seq(constant_adaptor_f, generic_dna).reverse_complement()) + constant_adaptor_r +'NNNNNN'


# In[82]:


# illumina_p5 = "AATGATACGGCGACCACCGAGATCTACAC[UMI][i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
# illumina_p7 = "CAAGCAGAAGACGGCATACGAGAT[UMI][i7]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
illumina_p5 = "AATGATACGGCGACCACCGAGATCTACAC"+'N'*8 +'AAAAAAAA'+"ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
illumina_p7 = "CAAGCAGAAGACGGCATACGAGAT"+'N'*8 +'TTTTTTTT'+"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
whole_ill_spacer =str(Seq(illumina_p5, generic_dna).reverse_complement()) + illumina_p7


# In[79]:


len(whole_ill_spacer)


# In[83]:


whole_ill_spacer

