import argparse
import sys

from Bio import motifs

from pwm import pfm2pssm, scan_all_seqs, scan_all_sequences, best_Of_All_window, personalized_best_window, \
    personalized_scan_sequences, best_of_One_window
from utils import download_promotors
# Parseur
parser = argparse.ArgumentParser()
parser.add_argument('l', nargs='+', help='<Required> Set flag')
parser.add_argument("-m", "--pfm", help="Position frequency matrix should be entred", action='store')
parser.add_argument("-t", "--threshold", type=float, help="Float threshold", action="store")
parser.add_argument("-l", "--promotor_lenght", type=int, help="Correct int lenght should be entred", action="store")
parser.add_argument("-w", "--window_size", type=int, help="Correct int window size should be entred", action="store")
parser.add_argument("-s", "--window_threshold", type=float, help="Correct int window threshold should be entred",
                    action="store")
parser.add_argument("-p", "--pseudocount", type=float, help="Correct float pseudo count should be entred",
                    action="store")

args = parser.parse_args()
# les if pour les valeurs par defaut et renommage
if args.promotor_lenght:
    longueur = args.promotor_lenght
else:
    longueur = 1000
if args.window_size:
    windowsSize = args.window_size
else:
    windowsSize = 40
if args.pseudocount:
    ps = args.pseudocount
else:
    ps = 0.1

jasparMatrix = "../data/" + args.pfm
score = args.threshold
windowsScore = args.window_threshold
listSeq = args.l
#listSeq=["NM_007389", "NM_079420", "NM_001267550", "NM_002470", "NM_003279", "NM_005159"
 #  , "NM_003281", "NM_002469", "NM_004997", "NM_004320", "NM_001100", "NM_006757"]

# telecharger les sequence en amont du gene
download_promotors(listSeq,longueur,"data")


with open(jasparMatrix) as handle:
    m = motifs.read(handle, "jaspar")

pssm = pfm2pssm(m, ps)
TFname = m.name
res = scan_all_seqs(pssm, listSeq, longueur, score)
count = 1
startCoordinate, endCoordinate = (0, windowsSize)
# affichage des fenetre ayant le meilleur score
for x in res:
    if (count < len(res) - 1):
        bestWindowX = best_of_One_window(res, startCoordinate, endCoordinate, longueur, windowsScore)
        print(str(count) + " " + str(TFname) + " [" + str(bestWindowX[0]) + ":" + str(bestWindowX[1]) + "] " + str(
            bestWindowX[2]))
        startCoordinate = endCoordinate
        endCoordinate += windowsSize
        count += 1

print("---------------------------------- \n List d'occurence de chaque fênetre"
      "\n-----------------------------------")

# liste des occurrences du motif de chaque fenetre
startCoordinate, endCoordinate = (0, windowsSize)
count = 0
for i in res:
    count += 1
    for y in i[0]:
     x = i[1]
     if startCoordinate <= y[0] <= endCoordinate and y[1] >= windowsScore and x<len(res):
            print(str(count) + " " + listSeq[x] + " " + str(TFname) + " " + str(y[0]) + " " + str(y[1]))

print("Meilleure fenêtre avec notre solution personalisée :\n")
print(personalized_best_window(res, windowsSize, longueur, windowsScore))
print("Meilleur fenêtre avec la solution d'occurence")
print(best_Of_All_window(res,windowsSize,longueur,windowsScore))
