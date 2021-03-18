import os

from Bio import SeqIO
from Bio import Entrez

Entrez.email = "malik.bouaoud.etu@univ-lille.fr"

def find_cds(seqRec):
    l = []
    for x in seqRec.features:
        if x.type == 'CDS':
            x1 = int(x.location.start)
            x2 = int(x.location.end)
            l.append((x1, x2))

    return l


def mrna_to_gen(idMRNA):
    try:
        handle2 = Entrez.elink(dbfrom="nucleotide", id=idMRNA, db="gene")
        record2 = Entrez.read(handle2)
        handle2.close()
        linked = [link["Id"] for link in record2[0]["LinkSetDb"][0]["Link"]]
        return linked[0]
    except ValueError:
        print("ERREUR DE ELINK")


def upstream_gene_seq(idGen,longueur):
    handle=Entrez.esummary(db="gene", id=idGen)
    record = Entrez.read(handle)
    handle.close()
    start = record["DocumentSummarySet"]["DocumentSummary"][0]['GenomicInfo'][0]['ChrStart']
    stop = record["DocumentSummarySet"]["DocumentSummary"][0]['GenomicInfo'][0]['ChrStop']
    NumeroAccession = record["DocumentSummarySet"]["DocumentSummary"][0]['GenomicInfo'][0]['ChrAccVer']
    if start>stop:
        valueSTRAND=2 #brin nge
    else :
        valueSTRAND=1 #brin pos
    s=SeqIO.read(
             Entrez.efetch(db="nucleotide", id=NumeroAccession, seq_start=start,retmax=longueur
                      , seq_stop=stop, rettype="gb", strand=valueSTRAND,retmode="text")
                          , "gb")
    return s.seq[:longueur]


def download_promotors(listMRna,sizeSeq,repository):
    path="../"+str(repository)+'/'
    for x in listMRna:
        idgen=mrna_to_gen(x)
        seq=upstream_gene_seq(idgen,sizeSeq)
        filepath = os.path.join(path, x+str('_')+str(sizeSeq)+str(".fa"))
        if not os.path.exists(path):
            os.makedirs(path)
        f=open(filepath,"a")
        f.write(str(seq))


