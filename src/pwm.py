from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq
from numpy import inf

def pfm2pssm(pfm, pseudoPoids):
    pwm = pfm.pwm
    return pwm.log_odds()


def scan_seq(pssm, seq, seuil):
    l = []
    for x in pssm.search(seq, seuil):
        if x[1] > seuil:
            l.append((x[0], x[1]))
    return l


# la structure sera une liste contenant un tuple : le premier element est la liste retourner
# par scan_seq et le deuxième est le numero de la sequence suivant l'ordre de la list donné en entrée
# Exemple list[tuple(List,Int = numéro de séquence)]

def scan_all_seqs(pssm, seqList, longueur, seuil):
    seqNumber = 1
    listOfAllSeq = []
    for x in seqList:
        path = "../data/" + x + "_" + str(longueur) + ".fa"
        with open(path, 'r') as file:
            data = file.read().replace('\n', '')
        sequence = Seq(data)

        # on récuère tout grace à - inf et non pas seuil
        l = scan_seq(pssm, sequence, seuil)
        listOfAllSeq.append((l, seqNumber))
        seqNumber = seqNumber + 1
    return listOfAllSeq


def scan_all_sequences(listOfScanAllSeqs, startCoordinate, endCoordinate, scoreWindow):
    score = 0
    for i in listOfScanAllSeqs:
        for y in i[0]:
            if startCoordinate <= y[0] <= endCoordinate and y[1] >= scoreWindow:
                score = score + 1
    return score


# return a tuple of three element start of the widows its end and its score
def best_Of_All_window(listOfScanAllSeqs, sizeWindow, longueur, scoreWindow):
    savedStart = -1
    savedEnd = -1
    startCoordinate = 0
    endCoordinate = sizeWindow
    l = []
    while startCoordinate < longueur and endCoordinate <= longueur:
        l.append(best_of_One_window(listOfScanAllSeqs, startCoordinate, endCoordinate, longueur, scoreWindow))
        startCoordinate=endCoordinate
        endCoordinate+=sizeWindow
    return getMax(l)

# retourn le promoteur de la fentre avec le meilleur score
def best_of_One_window(listOfScanAllSeqs, startCoordinate, endCoordinate, longueur, windowScore):
    max = scan_all_sequences(listOfScanAllSeqs, startCoordinate, endCoordinate, windowScore)
    for i in listOfScanAllSeqs:
        if startCoordinate < longueur and endCoordinate <= longueur:
            min = scan_all_sequences(listOfScanAllSeqs, startCoordinate, endCoordinate, windowScore)
            if min > max:
                max = min

    return startCoordinate, endCoordinate, max


# calcul la moyenne de la fenetre
def personalized_scan_sequences(listOfScanAllSeqs, startCoordinate, endCoordinate, windowScore):
    sum = 0
    count = 1
    for i in listOfScanAllSeqs:
        if startCoordinate <= i[0] <= endCoordinate and i[1] >= windowScore:
            sum = sum + i[1]
            count = count + 1
    return sum / count


# En resumé on calcul la moyenne de chaque sequence fenetrée, puis on la stoque dans une liste,
# On suite on calcul la moyenne generale de toutes les sequence fenetré
# le score serra alors le max des moyenne de chaques fenetre

def personalized_best_window(listOfScanAllSeqs, sizeWindow, longueur, windowScore):
    startCoordinate = 0
    endCoordinate = sizeWindow
    l = []
    for i in listOfScanAllSeqs:
        moyOne = 0
        if startCoordinate < longueur and endCoordinate <= longueur:
            moyOne = personalized_scan_sequences(i[0], startCoordinate, endCoordinate, windowScore)
            l.append((startCoordinate, endCoordinate, moyOne))
            startCoordinate = endCoordinate
            endCoordinate = endCoordinate + sizeWindow

    return getMax(l)


# Moyenne de toutes les listes
# recherche de la fenetre qui a une moyenne sup à moyAll et qui soit max
def getMax(listOfMoy):
    savedStart = -1
    savedEnd = -1
    maxMoy = listOfMoy[0][2]
    for x in listOfMoy:
        if x[2] > maxMoy:
            maxMoy = x[2]
            savedStart = x[0]
            savedEnd = x[1]
    return savedStart, savedEnd, maxMoy
