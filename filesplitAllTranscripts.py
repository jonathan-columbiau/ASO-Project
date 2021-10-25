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

#print(len(utr_dict))
#print(type(list(utr_dict.keys())[0]))
#format of key: NFKB2|ENST00000189444; type: str
#plot_hist - plots summary statistics of sequence lengths of all transcripts (assuming dictionary values are sequences in string format). Second parameter is title of histogram
'''
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

plot_hist(utr_dict)
'''
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
                key = gene_name + " | " + gene_transcript_id + " | " + str(num_lines)
                value = string_seq[count:count+270]
                dict_270nt_seq[key] = value
                count = count + 220
                if count > len_seq:#summary statistics
                    dict_number_seq[gene] = num_lines#summary statistics
                #50-overlap is created by just increasing by 220 nt each iteration over the string, but writing 270 nt
                num_lines = num_lines + 1
    return dict_270nt_seq, dict_number_seq



dict_270_nt, dict_number_seqs = dict_with_270nt_seq(utr_dict)
#write dict_number_seqs to a new file called seqs_per_transcript.txt. Can do analysis in Python, but easier in Excel. In Excel, just import this text file as CSV to open it.
def write_seqs_per_transcript(dict1):
    with open("seqs_per_transcript.txt", "w") as filevar:
        filevar.write("Gene name | transcript id | # of sequences\n")
        for gene in dict1:
            filevar.write(gene + "|" + str(dict1[gene]) + "\n")
write_seqs_per_transcript(dict_number_seqs)


print("Number of total 270 nt transcripts: " + str(len(dict_270_nt))) #32561 transcripts/stuff currently in the dictionary




#given a dictionary that can have keys with duplicate values, return a dictionary such that there is only one key per value, and other keys with duplicate values
#are removed. Write to a file called filename the keys that are removed from the dictionary.
def remove_duplicates(dict1, filename):
    dict_270_nt_no_dupl = {}
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
        for header in dict1:
            if header not in keys_not_deleted:
                filevar.write("\n" + header)
                count_deleted_duplicates = count_deleted_duplicates + 1
            if header in keys_not_deleted:
                dict_270_nt_no_dupl[header] = dict1[header]
                count_original = count_original + 1
        print("Number of duplicate sequences: " + str(count_deleted_duplicates))
        print("Number of non-duplicates: " + str(count_original))
    return dict_270_nt_no_dupl




#output sequences should be FASTA file. Then can take fragments for a single gene, copy and paste 10 entries, do Blast. Try before and after removing duplicates.

#create a dictionary where sequence is key, use id as the value. Then can loop through dictionary where utr is new key and header is new value, and output that as the final.
#if key already exists, just
#a lot of transcripts have identical protein coding sequence so same stop codon

#alternative way to remove duplicates after meeting with Xuebing: remove_duplicates2
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


#remove keys from original dictionary
no_dupl_dict = remove_duplicates(dict_270_nt, "deleted_transcripts.txt")
no_dupl_dict2 = remove_duplicates2(dict_270_nt, "deleted_transcripts2.txt")

count = 0
for key in no_dupl_dict2:
    if len(key) <= 15:
        count = count + 1
print("Number of fragments <= 15 nt " + count)

'''
#assertion test: fragments in no_dupl_dict = fragments in no_dupl_dict2 (they're just flipped so can't do a direct comparison). Keys are different though.
assert set(no_dupl_dict.values()) == set(no_dupl_dict2.keys())
#assertion test works - these have the same values
'''


'''
Make a file that shows the number of fragments per transcript before duplicate
removal and after (inputs: original dictionary, dictionary after removal);
output: writing to file in CSV format
'''

def make_fragment_deleted_table(original_dict, removal_dict):
    with open("deleted_transcripts_data.csv", "w") as filevar:
        split_key = origina



#Fasta format:
#header: >Gene name | transcript id | 270 nt fragment number for this id
#sequence on next ine
#There are gene paralogs with same gene id + gene name but different 3' UTR regions
def write_270_nt(dict1):
    with open("270ntsequences.fasta", "w") as filevar:
        #filevar.write("Gene name | transcript id | 270 nt sequence number for this id | sequence\n")
        for gene in dict1:
            filevar.write(">" + gene + "\n" + dict1[gene] + "\n")
            #now handle this just like a string, and write to a file just like a string.
write_270_nt(no_dupl_dict)
