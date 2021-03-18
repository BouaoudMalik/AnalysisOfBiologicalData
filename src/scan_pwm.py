import sys
from Bio import motifs, Entrez
from Bio.Seq import Seq

from utils import find_cds, mrna_to_gen, upstream_gene_seq, download_promotors
from pwm import pfm2pssm, scan_seq, scan_all_seqs, scan_all_sequences, personalized_best_window, \
    best_Of_All_window

Entrez.email = "malik.bouaoud.etu@univ-lille.fr"
jasparMatrix = "../data/" + str(sys.argv[1])  # MA0056.1
genbankArnm = str(sys.argv[2])  # "NM_007389"
tailleSequence = int(sys.argv[3])  # 1000
score = int(sys.argv[4])  # -2
ps = 0.1

gene = mrna_to_gen(genbankArnm)
seq = upstream_gene_seq(gene, tailleSequence)

with open(jasparMatrix) as handle:
    m = motifs.read(handle, "jaspar")
pssm = pfm2pssm(m, ps)

scanSeqList=scan_seq(pssm, seq, score)

l = ["NM_007389", "NM_079420", "NM_001267550", "NM_002470", "NM_003279", "NM_005159"
    , "NM_003281", "NM_002469", "NM_004997", "NM_004320", "NM_001100", "NM_006757"]
download_promotors(l,1000,"data")
res = scan_all_seqs(pssm, l, tailleSequence, score)
TFname=m.name

for x in res:
    for y in x[0]:
        if  y[1] >= score and y[0]<tailleSequence:
             print(str(TFname)+" "+str(y[0]) +" "+str(y[1]))

