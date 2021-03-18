import unittest

from src.pwm import scan_all_seqs, personalized_best_window, best_Of_All_window
from src.utils import find_cds, mrna_to_gen, upstream_gene_seq
from Bio import SeqIO
from Bio import Entrez

Entrez.email = "malik.bouaoud.etu@univ-lille.fr"


class MyTestCase(unittest.TestCase):
    def testIfGenAndFastaSequencesEquals(self):
        fasta = SeqIO.read("data/sequence.fasta", "fasta")
        gb = SeqIO.read("data/sequence.gb", "genbank")
        seqFasta = str(fasta.seq)
        seqGb = str(gb.seq)
        self.assertTrue(fasta.seq == gb.seq)

    def testFindCds(self):
        gb = SeqIO.read("data/sequence.gb", "genbank")
        l = find_cds(gb)
        self.assertTrue([(51, 1425)] == l)

    def testIfContenentImportedFromSeqIoReadSameAsLocal(self):
        gb = SeqIO.read("data/sequence.gb", "genbank")
        gbwithHandle = SeqIO.read(Entrez.efetch(db="nucleotide", id="NM_007389", rettype="gb", retmode="text")
                                  , "gb")
        self.assertTrue(find_cds(gbwithHandle) == find_cds(gb))

    def testIfAspecifiedIDIsInElinkReturnedList(self):
        handle2 = Entrez.elink(dbfrom="nucleotide", id="NM_007389", db="gene")
        record2 = Entrez.read(handle2)
        handle2.close()
        linked = [link["Id"] for link in record2[0]["LinkSetDb"][0]["Link"]]

        self.assertTrue("11435" in linked)

    def testMrnaTOGEn(self):
        self.assertTrue(int(mrna_to_gen("NM_007389")) == 11435)

    def testForEssumary(self):
        handle3 = Entrez.esummary(db="gene", id="11435")
        record3 = Entrez.read(handle3)
        handle3.close()
        print(record3)
        genInfo = NumeroAccession = record3["DocumentSummarySet"]["DocumentSummary"][0]['GenomicInfo'][0]
        NumeroAccession = record3["DocumentSummarySet"]["DocumentSummary"][0]['GenomicInfo'][0]['ChrAccVer']
        pos1 = record3["DocumentSummarySet"]["DocumentSummary"][0]['GenomicInfo'][0]['ChrStart']
        pos2 = record3["DocumentSummarySet"]["DocumentSummary"][0]['GenomicInfo'][0]['ChrStop']
        self.assertTrue(str('NC_000068.8') == str(NumeroAccession))
        self.assertTrue((73410681, 73393624) == (int(pos1), int(pos2)))


if __name__ == '__main__':
    unittest.main()
