from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

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


#make a function that converts a dictionary with sequences in Biopython's SeqRecord format to string format. Input: dictionary w values in SeqRecord for,at.
#output: dictionary with String values.
def seqdict_to_strdict(dict1):
    strdict = {}
    for gene in dict1:
        sequence = dict1[gene].seq
        string_seq = str(sequence)
        strdict[gene] = string_seq
    return strdict

#Find gene IDs that aren't
#Input: CSV with gene names of all genes you want data for, and dictionary with header in FASTA format where first thing before | is the name of the gene
#plan: convert all of these to a list
#will manually get data for the last 19 genes that ENSEMBL couldn't get MANE 3' information for.
#output: write to file called add_manually with a list of genes I need to get 3' UTR info for.
def find_genes(dict1):
    with open("CSV of gene names of all 660 haploinsufficient genes.csv", "r") as filevar:
        first_line_with_names = filevar.readline()
        gene_list = first_line_with_names.split(",")
    gene_list_dict = []
    for header in dict1:
        gene_name = header.split("|")[0]
        gene_list_dict.append(gene_name)
    genes_manually_add = []
    for gene in gene_list:
        if gene not in gene_list_dict:
            genes_manually_add.append(gene)
    with open("Gene Names to add.txt", "w") as filevar2:
        for gene in genes_manually_add:
            filevar2.write(gene)
            filevar2.write("\n")


#plot_hist - plots summary statistics of sequence lengths of all transcripts (assuming dictionary values are sequences in string format). Second parameter is title of histogram
def plot_hist(dict1):
    dict_hist = {}
    for transcript, utr in dict1.items():
        num_seq = len(utr)
        dict_hist[transcript] = num_seq
    values = list(dict_hist.values())
    min_utr_len = min(values)
    max_utr_len = max(values)
    #algorithm for binning: https://stackoverflow.com/questions/6855710/how-to-have-logarithmic-bins-in-a-python-histogram ; I chose arbitrary # of bins at 30
    #numpy linspace package returns evenly spaced numbers over a specified interval
    plt.hist(values, bins = 10 ** np.linspace(np.log10(min_utr_len), np.log10(max_utr_len), 30), histtype = "bar")
    plt.xscale('log')
    plt.xlabel("Sequence Length")
    plt.ylabel("Frequency")
    plt.show()


#return dictionary with keys being in the format:
#Gene name | transcript id | 270 nt fragment number for this transcript id
#and the value is the sequence. Will use this to check for, and remove, duplicates in the next step; and after, write the sequences to a file in FASTA format.
#Addition: add additional thing you return, a dictionary with the number of 270nt sequences needed for 3' region for each gene. Useful for summary statistics.
def dict_with_270nt_seq(dict1):
    dict_270nt_seq = {}
    dict_number_seq = {}
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
                key = gene_name + "|" + gene_transcript_id + "|" + str(num_lines)
                value = string_seq[count:count+270]
                dict_270nt_seq[key] = value
                count = count + 220
                if count > len_seq:#summary statistics
                    dict_number_seq[gene] = num_lines#summary statistics
                #50-overlap is created by just increasing by 220 nt each iteration over the string, but writing 270 nt
                num_lines = num_lines + 1
    return dict_270nt_seq, dict_number_seq

#write sequences before removal to file - debugging function
def write_seqs_before_dup_removal(dict1):
    with open("seqs_before_duplication_removal.fasta", "w") as filevar:
        for key in dict1:
            filevar.write(">" + key + "\n" + dict1[key] + "\n")


#write dict_number_seqs to a new file called seqs_per_transcript.txt. Can do analysis in Python, but easier in Excel. In Excel, just import this text file as CSV to open it.
def write_seqs_per_transcript(dict1):
    with open("seqs_per_transcript.txt", "w") as filevar:
        filevar.write("Gene name | transcript id | # of sequences\n")
        for gene in dict1:
            filevar.write(gene + "|" + str(dict1[gene]) + "\n")




#given a dictionary that can have keys with duplicate values, return a dictionary such that there is only one key per value, and other keys with duplicate values
#are removed. Write to a file called filename the keys that are removed from the dictionary.
#return another dictionary with numbers of each transcript id and fragment, to compare later for quality control testing.
#second dictionary returns gene_name|transcript_id as key and value as #times it appears, will compare to dictionary before adding.
def remove_duplicates(dict1, filename):
    dict_270_nt_no_dupl = {}
    dict_qc_test = {}
    count_deleted_duplicates = 0
    count_original = 0
    with open(filename, "w") as filevar:
        filevar.write("Deleted Sequences\nGene Name | Transcript ID | Transcript #\n")
        non_duplicate_values_list = []
        keys_not_deleted = []
        for header,utr in dict1.items():
            if utr not in non_duplicate_values_list:
                non_duplicate_values_list.append(utr)
                keys_not_deleted.append(header)
                #qc part
                split_list = header.split("|")
                gene_name = split_list[0]#weird numbers here to make same format to compare to other dictionary in qc test
                transcript_id = split_list[1]
                qc_key = gene_name +"|" + transcript_id
                if qc_key not in dict_qc_test:
                    dict_qc_test[qc_key] = 0
                dict_qc_test[qc_key] = dict_qc_test[qc_key] + 1

        for header in dict1:
            if header not in keys_not_deleted:
                filevar.write("\n" + header)
                count_deleted_duplicates = count_deleted_duplicates + 1
            if header in keys_not_deleted:
                dict_270_nt_no_dupl[header] = dict1[header]
                count_original = count_original + 1
        print("Number of duplicate sequences: " + str(count_deleted_duplicates))
        print("Number of non-duplicates: " + str(count_original))
    return dict_270_nt_no_dupl, dict_qc_test


#alternative way to remove duplicates: remove_duplicates2
def remove_duplicates2(dict1, filename):
    dict_270_nt_no_dupl = {}
    with open(filename, "w") as filevar:
        filevar.write("Deleted Sequences\nGene Name | Transcript ID | Fragment # | Replaced by Fragment \n") # add replaced by
        for key in dict1:
            if dict1[key] in dict_270_nt_no_dupl:
                filevar.write(key + " | " + dict_270_nt_no_dupl[dict1[key]] + "\n" + dict1[key] + "\n")
            dict_270_nt_no_dupl[dict1[key]] = key
    print("Number of non-duplicates remove_duplicates2: " + str(len(dict_270_nt_no_dupl)))
    return dict_270_nt_no_dupl



'''
Make a file that shows the number of fragments per transcript before duplicate
removal and after (inputs: original dictionary, dictionary after removal);
Output: writing to file in CSV format. Both dictionaries have the same keys if they both have the item.
Format:
Gene Id | Transcript ID | Number of fragments before | Number of fragments after
'''
def make_fragment_deleted_table(original_dict_nums, dict_after_removal_nums):
    #print("dict_after_removal_nums keys", list(dict_after_removal_nums.keys()))
    with open("deleted_transcripts_data.tsv", "w") as filevar:
        for key in original_dict_nums:
            filevar.write(key  + "|" + str(original_dict_nums[key]) + "|")
            if key in dict_after_removal_nums:
                filevar.write(str(dict_after_removal_nums[key]))
            else:
                filevar.write("0")
            filevar.write("\n")

#filter 1: If 1 fragment is subset of another, then remove the smaller one.
def filter_subsequence(dict1):
    values_list = list(dict1.values())
    new_dict = {}
    for key in dict1:
        add_or_not = True
        for value in values_list:
            if dict1[key] != value and dict1[key] in value:
                add_or_not = False
        if add_or_not == True:
            new_dict[key] = dict1[key]
    print("Length of new dictionary: " + str(len(new_dict)))
    return new_dict

#count number of fragments smaller than a particular size
def count_15(dict1):
    count200 = 0
    count10 = 0
    for x in dict1:
        if len(dict1[x]) <= 200:
            count200 = count200 +1
        if len(dict1[x]) <= 10:
            count10 = count10 +1
    print("Number < 200:  ", count200, "\n")
    print("Number < 10:  ", count10, "\n")

#Write Fragments to file in Fasta format:
#header: >Gene name | transcript id | 270 nt fragment number for this id
#sequence on next ine
def write_270_nt(dict1):
    with open("270ntsequences.fasta", "w") as filevar:
        #filevar.write("Gene name | transcript id | 270 nt sequence number for this id | sequence\n")
        for gene in dict1:
            filevar.write(">" + gene + "\n" + dict1[gene] + "\n")
            #now handle this just like a string, and write to a file just like a string.

#for 3' UTR sequences initially smaller than 270 nt, add GFP to end of that sequence until
#it reaches 270 nt. Return a new dictionary with all sequences initially larger than 270 nt, and those
#that are smaller than 270 nt corrected to be 270 nt by appending GFP to the end.
def add_GFP_to_270(dict1):
    gfp = "ATGGTGAGCAAGGGCGCCGAGCTGTTCACCGGCATCGTGCCCATCCTGATCGAGCTGAATGGCGATGTGAATGGCCACAAGTTCAGCGTGAGCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCTGTGCCCTGGCCCACCCTGGTGACCACCCTGAGCTACGGCGTGCAGTGCTTCTCACGCTACCCCGATCACATGAAGCAGCACGACTTCTTCAAGAGCGCCATGCCT"
    new_dict = {}
    for seq in dict1:
        if len(dict1[seq]) > 270:
            new_dict[seq] = dict1[seq]
        else:
            len_utr = len(dict1[seq])
            new_seq = dict1[seq] + gfp[:270-len_utr]
            new_dict[seq] = new_seq
    return new_dict

def run_program():
    utr_dict = SeqIO.to_dict(SeqIO.parse("UTRs of haploinsufficient genes MANE curated 641 out of 660.fasta", "fasta"))
    find_genes(utr_dict)
    utr_dict = filter_unavailable(utr_dict)
    utr_dict = seqdict_to_strdict(utr_dict)
    #plot_hist(utr_dict)
    utr_dict = filter_subsequence(utr_dict) #new 11/9
    utr_dict = add_GFP_to_270(utr_dict)
    dict_270_nt, dict_number_seqs = dict_with_270nt_seq(utr_dict)
    write_seqs_before_dup_removal(dict_270_nt)
    write_seqs_per_transcript(dict_number_seqs)
    print("Number of total 270 nt transcripts: " + str(len(dict_270_nt))) #32561 transcripts/stuff currently in the dictionary
    no_dupl_dict, qc_dict = remove_duplicates(dict_270_nt, "deleted_transcripts.txt")
    no_dupl_dict2 = remove_duplicates2(dict_270_nt, "deleted_transcripts2.txt")
    #assertion test: fragments in no_dupl_dict = fragments in no_dupl_dict2 (they're just flipped so can't do a direct comparison). Keys are different though.
    assert set(no_dupl_dict.values()) == set(no_dupl_dict2.keys())
    make_fragment_deleted_table(dict_number_seqs, qc_dict)
    subsequence_dict = filter_subsequence(no_dupl_dict)
    count_15(subsequence_dict)
    write_270_nt(subsequence_dict)

run_program()
