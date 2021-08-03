from Bio import pairwise2
import pandas as pd
import xlsxwriter
import os
from time import time

class alignment_score:
    #constructor
    def __init__(self, seq1, seq2):
        self.X = seq1
        self.Y = seq2
        self.max_len = max(len(seq1), len(seq2))
        self.global_alg()
        # self.local_alg()
    
    def global_alg(self):
        alignments = pairwise2.align.globalms(self.X, self.Y, 1, -0.5, -0.25, -0.1)
        self.score_global = self.max_score(alignments)

    def local_alg(self):
        alignments = pairwise2.align.localms(self.X, self.Y, 1, -0.5, -0.25, -0.1)
        self.score_local = self.max_score(alignments)
    
    def max_score(self, alignment_list):
        _score = []
        for a in alignment_list: _score.append(a.score/self.max_len)
        return (max(_score)*100, min(_score)*100, (sum(_score)/len(_score))*100)

class align_softberry_to_genebank:
    def __init__(self, gb_f, sb_f):
        self.sb_loc, self.sb_aa_seq, self.sb_aa_len, self.sb_n = [], [], [], 0
        self.gb_loc, self.gb_aa_seq, self.gb_aa_len, self.gb_n = [], [], [], 0
        self.align_completely, self.align_comp_no, self.not_align_completely = [], 0, 0

        self.sb_gene_list(sb_f)
        self.gb_gene_list(gb_f)
        self.align_sb_gb()
    
    def sb_gene_list(self, sb_f_h):
        seq = []
        for i, line in enumerate(sb_f_h):
            if i < 3: continue
            elif line.startswith('>Potential'):
                if len(seq) != 0: self.sb_aa_seq.append(''.join(seq))
                l = line.split()
                l_ = l[2].split(':')[1].split('.')
                self.sb_loc.append((int(l_[0]), int(l_[2].split(',')[0])))
                self.sb_aa_len.append(int(l[-2]))
                seq = []       
            else:
                seq.append(line.split('\n')[0])
        self.sb_aa_seq.append(''.join(seq))

    def gb_gene_list(self, gb_f_h):
        for i, line in enumerate(gb_f_h):
            if i < 3: continue
            elif i == 3: self.gb_n = int(line.split('\n')[0])
            elif i % 2 != 0:
                seq = line.split('\n')[0]
                self.gb_aa_seq.append(seq)
                self.gb_aa_len.append(len(seq))
            else:
                l = line.split()
                self.gb_loc.append((int(l[0]), int(l[2])))

    def align_sb_gb(self):
        for sb_ind, l in enumerate(self.sb_loc):
            if l not in self.gb_loc: 
                self.not_align_completely += 1
                continue
            gb_ind = self.gb_loc.index(l)
            # print(l, self.sb_aa_seq[sb_ind], self.gb_aa_seq[gb_ind], sep='\n')
            a_s = alignment_score(self.sb_aa_seq[sb_ind], self.gb_aa_seq[gb_ind])
            self.align_completely.append(l)
            self.align_completely.append(a_s.score_global)
            self.align_comp_no += 1

files = sorted(list(map(int, (os.listdir('vgas/vgas_s/')))))
data = []
for i, a in enumerate(files):
    f1, f2 = open('genebank_1/'+str(a), 'r'), open('vgas/vgas_s/'+str(a), 'r')
    align = align_softberry_to_genebank(f1, f2)
    # print(time()-start)
    row = align.align_completely
    row.insert(0, align.not_align_completely)
    row.insert(1, align.align_comp_no)
    row.insert(0, str(a)+'.txt')
    data.append(row)
    print(i+1, 'Done')

df = pd.DataFrame(data)

writer = pd.ExcelWriter('Vgas.xlsx', engine='xlsxwriter')

# Convert the dataframe to an XlsxWriter Excel object.
df.to_excel(writer, sheet_name='Sheet1', index=False)

# Close the Pandas Excel writer and output the Excel file.
writer.save()