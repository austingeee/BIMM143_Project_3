#!/usr/bin/env python
# coding: utf-8

# # BIMM143 Project 3

# # Intoduction
# scientific question: How similar are the genes for Casein Kinase 1 Delta, a protein that helps control circadian rhythm, between humans and different species?
# 
# Casein Kinase 1 Delta (CK1D) is a protein kinase that phosphorylates a number of different proteins in the metabolic pathways for circadian rhytm. CK1D is known to specifically target PER1/2 genes, which controls circadian cycle length. Additionally, CK1D can also itself be phosphorylated, which directly effects its function, and thus the function of proteins downstream in the pathway (https://www.pnas.org/content/115/23/5986).
# 
# It seems that most complex eukaryotes contain genes for Casein Kinase 1 Delta. This project aims to determine wether or not amino acid sequecnes for CK1D is conserved across different species, and if so to what degree. 
# 
# Hypothesis: If an organism is closer in evolutionary relation to Humans (Homo Sapiens), then the alignment score between those two sequences will be higher. 
# 
# Protein sequneces for CK1D were collected for multible organisms: Humans, Fruit Bats, Zebrefish, House Mouse, Chimpanzees, and Clawed Frogs. A pairwise alignment method was used to compare each of the organisms aino acid sequences, using Humans as a reference for all comparisons. The sequences used for the alignment were found on the National Center for Biotehnology Information (NCBI, https://www.ncbi.nlm.nih.gov/). The sequences were saved indivisually as .txt files to be read in the actual code. 
# . 

# # Part 1: Load the Packages.
# 
# The two packages that are needed for this code are:
# 
# BioPython: a bioinformatics package that contains different functions for analyzing biological data. For this aplication, the pairwise2 function was imported to be used to align the different protein sequences. 
# 
# matplotlib: a package that contains functions for creating a multitude of different graphs and tables for displaying data. For this application, the pyplot subpackage was used to access the function to create a table, which is how the data will be presented. 
# 

# # Part 2: load in the data and perform analyses
# 
# text files are simple files containing any desired text. For this project the text files used contained the amino acid sequences for CK1D for various different species. In the code below, a function is defined called "align()" which takes two sequecnes and aligns them. The .txt files are read as strings, and the strings are then passed through the pairwise2 alingment function from BioPython. Four scoring values are used in the alignment: +1 for matches in the sequence, -1 for a mismatch, -1 for a gap, and -2 for a gap 2 amino acids or longer. The 'align()' function is then passed through each organism, using Humans as the reference

# In[9]:


import matplotlib.pyplot as plt
from Bio import pairwise2

def align(ref_seq, seq1):
    "Align and score two amino acid sequences using two .txt files containing the sequences"
    
    # open sequence files - two text files containing desired amino acid sequences to compare
    
    ref_seq = open(ref_seq, 'r') 
    seq1 = open(seq1, 'r')
    
    # set alignment function: pairwise alignment of two sequecnes
    
    alignment = pairwise2.align.globalms(ref_seq.read(), seq1.read(), 1, -1, -2, -1, score_only=True)
    result = alignment
    # arguments 1, -1, -2, -1 are scoring parameters for a match, mismatch, gap extention, and gap opneing, respectively
    return(result)

# define all desired comparisons and store as seperate variables

human_mice = align('Homo_Sapiens_Seq.txt', 'House_Mice_seq.txt') 
human_frog = align('Homo_Sapiens_Seq.txt', 'Clawed_Frog_Seq.txt')
human_fish = align('Homo_Sapiens_Seq.txt', 'Zebra_Fish_Seq.txt')
human_bat = align('Homo_Sapiens_Seq.txt', 'Fruit_Bat_Seq.txt')
human_chimp = align('Homo_Sapiens_Seq.txt', 'Chimp_Seq.txt')

# creation of table for display of scores

plt.table(cellText=[['Zebrafish', str(human_fish)],
                    ['Clawed Frog', str(human_frog)],
                    ['House Mouse', (human_mice)],
                   ['Fruit Bat', str(human_bat)],
                   ['Chimp', str(human_chimp)]],
          cellLoc='center', 
          colLabels=['Comparison w/ Human Reference', 'Alignment Score'], loc = 'top', 
          colColours = ['lightblue', 'lightblue'])

# lines of code to remove plot that is automatically created using the table function (not needed)

x = plt.gca() 
x.get_xaxis().set_visible(False)
x.get_yaxis().set_visible(False)
plt.box(on=None)

# Pritning Table with description

print('   Alignment Scores for Casein Kinase 1 Delta')
print('      Between Humans and Different Species')
plt.show()


# #  Part 3: displaying scores in a table
# 
# The result of this code gives a table showing the 5 different species comparisons, along with their respective alignment scores. 

# #  Part 4: Analyzing results 
# 
# The data in the table defintely shows a trend in alignment scores across teh different species. For example, starting with the zebrea fish, the socre for this alignment is 304. The Chimp alignment however gave back a score of 404, which is defintely higher than the zebrafish score. Since Chimpanzees are much closely realted to humans than the zebra fish, this data supports the hypothesis that evolutionariy closer organisms to humans will have higher alignment scores compared to more distant relativs. This hypothesis is further supported by the trend seen in the data. As species get closer to humans on the evolutionary tree, the scores increse, indicating a more conerved sequence. 
# Another interesting thing to note is that although different species were tested, the amino acid sequence scores are still pretty high compared to the closest match (the chimp). This shows that the amino acid sequecnes for CK1D seems relatively conserved across all the species tested. 

# 
