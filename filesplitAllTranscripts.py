from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

#Put each sequence into a dictionary,with the key being the header in the FASTA
#http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec63

utr_dict = SeqIO.to_dict(SeqIO.parse("All660UTRsWithTranscriptIDs.fasta", "fasta"))

'''
This block of code filters for 3' UTR sequences at least 15 nt and unavailable sequences
del_items = [] #contain keys of sequences that will be deleted
'''
def filter_unavailable(dict_seq):
    del_items = []
    for gene in dict_seq:
        seq = dict_seq[gene].seq
        if  seq == "Sequenceunavailable":  #can add or len(seq) < 15 to filter for nucleotides 
            del_items.append(gene)
    filtered_dict = dict_seq
    for key in del_items:
        del filtered_dict[key]
    return filtered_dict

utr_dict = filter_unavailable(utr_dict)

'''
#can make histogram function later - given a dictionary where values are strings, plot a histogram of all the string lengths.

#plot histogram - summary statistics for unfiltered UTR regions
#make a dictionary the same as utr_dict called utr_dict_hist except the values of the keys are the lengths of the 3' UTR regions
utr_dict_hist = {}
for gene, utr in utr_dict.items():
    num_seq = len(str(utr_dict[gene].seq))
    utr_dict_hist[gene] = num_seq

#get list of values from the histogram dictionary
values = list(utr_dict_hist.values())

#compute min and max of 3' UTR lengths 
min_utr_len = min(values)
max_utr_len = max(values)
#algorithm for binning: https://stackoverflow.com/questions/6855710/how-to-have-logarithmic-bins-in-a-python-histogram ; I chose arbitrary # of bins at 30
#numpy linspace package returns evenly spaced numbers over a specified interval 
plt.hist(values, bins = 10 ** np.linspace(np.log10(min_utr_len), np.log10(max_utr_len), 30), histtype = "bar")
plt.xscale('log')
plt.xlabel("Sequence Length")
plt.ylabel("Frequency")
plt.show()

'''


#check if there are duplicates in 270ntsequences.txt by first adding all the 270 nt sequences to a dictionary, before writing to 270ntsequences.txt

#make a function that converts a dictionary with sequences in Biopython's SeqRecord format to string format. Input: dictionary w values in SeqRecord for,at.
#output: dictionary with String values.
def seqdict_to_strdict(dict1):
    strdict = {}
    for gene in dict1:
        sequence = dict1[gene].seq
        string_seq = str(sequence)
        strdict[gene] = string_seq
    return strdict

utr_dict = seqdict_to_strdict(utr_dict)

#return dictionary with keys being in the format:
#Gene name | transcript id | 270 nt sequence number for this id | sequence
#and being the sequence. Will use this to check for, and remove, duplicates in the next step; and after, write the sequences to a file.  
def dict_with_270nt_seq(dict1):
    dict_270nt_seq = {}
    for gene in dict1:
            gene_header_list = gene.split("|")
            gene_name = gene_header_list[0]
            gene_transcript_id = gene_header_list[1]
            string_seq = dict1[gene]
            #now handle this just like a string, and write to a file just like a string.
            num_lines = 1
            count = 0
            len_seq = len(string_seq)
            while count < len_seq:
                key = gene_name + " | " + gene_transcript_id + " | " + str(num_lines)
                value = string_seq[count:count+270]
                dict_270nt_seq[key] = value
                count = count + 220
                #50-overlap is created by just increasing by 220 nt each iteration over the string, but writing 270 nt
                num_lines = num_lines + 1
    return dict_270nt_seq

dict_270_nt = dict_with_270nt_seq(utr_dict)


#given a dictionary that can have keys with duplicate values, return a dictionary such that there is only one key per value, and other keys with duplicate values
#are removed. Write to a file called filename the keys that are removed from the dictionary. 
def remove_duplicates(dict1, filename):
    dict_270_nt_no_dupl = {}
    count = 0
    with open(filename, "w") as filevar:
        filevar.write("Deleted Sequences\nGene Name | Transcript ID | Transcript #")                    
        non_duplicate_values_list = []
        keys_not_deleted = []
        for gene,utr in dict1.items():
            if utr not in non_duplicate_values_list:
                non_duplicate_values_list.append(utr)
                keys_not_deleted.append(gene)
        for gene in dict1:
            if gene not in keys_not_deleted:
                filevar.write("\n" + gene)
                count = count + 1
            if gene in keys_not_deleted:
                dict_270_nt_no_dupl[gene] = dict1[gene]
    print(count)
    print("len no dupl" + str(len(dict_270_nt_no_dupl)))
    

    
#remove keys from original dictionary
no_dupl_dict = remove_duplicates(dict_270_nt, "deleted_transcripts.txt")




#format:
#Gene name | transcript id | 270 nt sequence number for this id | sequence
#There are gene paralogs with same gene id + gene name but different 3' UTR regions
def write_270_nt(dict1):
    with open("270ntsequences.txt", "w") as filevar:
        filevar.write("Gene name | transcript id | 270 nt sequence number for this id | sequence\n") 
        for gene in dict1:
            filevar.write(gene + dict1[gene])
            #now handle this just like a string, and write to a file just like a string.


write_270_nt(no_dupl_dict)


#check if there are duplicates - print length of set vs length of
#No duplicates


#can plot number of fragments do that
#can check for genes with multiple transcript - look at if they have overlapping sequences. IF two isoforms have the same
#manually check some of these big genes to see if you have same UTR sequence
#just look at like 1 or 2 genes

#algorithm to write sequences to a file
#make a dictionary where keys are FASTA headers and values are SeqRecord objects that you could get the sequence values
#Ex: ["AAACCACAC...", "AATTGTGTTGC...", ...]
#make a for loop that just iterates over the list and writes the name of the gene
    #+ the position in the list + sequence at that position + new line character
    #ex: write("actin1 - AAA...\n")
#1. get all the keys of the dictionary


#test - for COMT [0:270] should be TGTGGACACACTGGACATGGTCTTCCTCGACCACTGGAAGGACCGGTACCTGCCGGACACGCTTCTCTTGGAGGAATGTGGCCTGCTGCGGAAGGGGACAGTGCTACTGGCTGACAACGTGATCTGCCCAGGTGCGCCAGACTTCCTAGCACACGTGCGCGGGAGCAGCTGCTTTGAGTGCACACACTACCAATCGTTCCTGGAATACAGGGAGGTGGTGGACGGCCTGGAGAAGGCCATCTACAAGGGCCCAGGCAGCGAAGCAGGGCC
    #       COMT [220:490] should be  GACGGCCTGGAGAAGGCCATCTACAAGGGCCCAGGCAGCGAAGCAGGGCCCTGACTGCCCCCCCGGCCCCCCTCTCGGGCTCTCTCACCCAGCCTGGTACTGAAGGTGCCAGACGTGCTCCTGCTGACCTTCTGCGGCTCCGGGCTGTGTCCTAAATGCAAAGCACACCTCGGCCGAGGCCTGCGCCCTGACATGCTAACCTCTCTGAACTGCAACACTGGATTGTTCTTTTTTAAGACTCAATCATGACTTCTTTACTAACACTGGCTA
    #       COMT [440:710] should be  GATTGTTCTTTTTTAAGACTCAATCATGACTTCTTTACTAACACTGGCTAGCTATATTATCTTATATACTAATATCATGTTTTAAAAATATAAAATAGAAATTAAGAATCTAAA
    