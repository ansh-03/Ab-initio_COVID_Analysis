import csv
from sys import stdin, stdout
from Bio import Entrez
from Bio import SeqIO
from time import time
from selenium import webdriver

accession_no = list()

# exporting accession numbers of the genomes
def getting_list():
    with open('R 26_3.csv') as f:
        reader = csv.reader(f)
        a = 0
        for row in reader:
            if a == 0: a+=1
            else: accession_no.append(row[0])

# retrieving data from GenBank
def data_gathering(accession_no_list, starting_point):
    Entrez.email = 'anshulg0303@gmail.com'
    count = starting_point

    for i_d in accession_no_list:
        with Entrez.efetch(db = "nucleotide", id = i_d, rettype = "gb", retmode = "csv") as net_handle:
            seq_r = SeqIO.read(net_handle, 'gb')
            features = seq_r.features
            row = [i_d+'\n', str(int(features[0].location.end))+'\n', str(seq_r.seq)+'\n']
            number = 0

            for a in range(len(features)):
                if features[a].type == 'CDS':
                    start, stop = int(features[a].location.start) + 1, int(features[a].location.end)
                    location = '%i - %i\n'%(start, stop)
                    sequence = features[a].qualifiers['translation'][0]
                    row += [location, sequence+'\n']

                elif features[a].type == 'gene':
                    number += 1
            
            row.insert(3, str(number)+'\n')

            with open('genebank\%i'%(count), 'w') as f:
                f.writelines(row)
        
        print('%i'%(count) + ' Done')
        count += 1

# i = time()
getting_list()
data_gathering(accession_no[66573:], 66573)
# f = open('runtime.txt', 'a')
# f.write('Gathering GeneBank Data: '+'%i'%(time()-i)+''+'s')
# f.close()