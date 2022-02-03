import math


class Bioinformatics(object):


    def __init__(self):
        pass

    def findAminoAcidFromCodon(self,codons):
        codon_full = {
            'UUU': 'Phe',
            'UUC': 'Phe',
            'uua': 'Leu',
            'UUG': 'Leu',
            'CUU': 'Leu',
            'CUC': 'Leu',
            'CUA': 'Leu',
            'CUG': 'Leu',
            'AUU': 'Ile',
            'AUC': 'Ile',
            'AUA': 'Ile',
            'AUG': 'Met',
            'GUU': 'Val',
            'GUC': 'Val',
            'GUA': 'Val',
            'GUG': 'Val',
            'UCU': 'Ser',
            'UCC': 'Ser',
            'UCA': 'Ser',
            'UCG': 'Ser',
            'CCU': 'Pro',
            'CCC': 'Pro',
            'CCA': 'Pro',
            'CCG': 'Pro',
            'ACU': 'Thr',
            'ACC': 'Thr',
            'ACA': 'Thr',
            'ACG': 'Thr',
            'GCU': 'Ala',
            'GCC': 'Ala',
            'GCA': 'Ala',
            'GCG': 'Ala',
            'UAU': 'Tyr',
            'UAC': 'Tyr',
            'UAA': 'STOP',
            'UAG': 'STOP',
            'CAU': 'His',
            'CAC': 'His',
            'CAA': 'Gln',
            'CAG': 'Gln',
            'AAU': 'Asn',
            'AAC': 'Asn',
            'AAA': 'Lys',
            'AAG': 'Lys',
            'GAU': 'Asp',
            'GAC': 'Asp',
            'GAA': 'Glu',
            'GAG': 'Glu',
            'UGU': 'Cys',
            'UGC': 'Cys',
            'UGA': 'STOP',
            'UGG': 'Trp',
            'CGU': 'Arg',
            'CGC': 'Arg',
            'CGA': 'Arg',
            'CGG': 'Arg',
            'AGU': 'Ser',
            'AGC': 'Ser',
            'AGA': 'Arg',
            'AGG': 'Arg',
            'GGU': 'Gly',
            'GGC': 'Gly',
            'GGA': 'Gly',
            'GGG': 'Gly'
        }
        codon_three_to_one = {
            "Lys" : "K",
            "Asn" : "N",
            "Thr" : "T",
            "Ser" : "S",
            "Ile" : "I",
            "Gln" : "Q",
            "Met" : "M",
            "His" : "H",
            "Pro" : "P",
            "Arg" : "R",
            "Leu" : "L",
            "Glu" : "E",
            "Asp" : "D",
            "Ala" : "A",
            "Gly" : "G",
            "Val" : "V",
            "Stp" : "O",
            "Tyr" : "Y",
            "Cys" : "C",
            "Trp" : "W",
            "Phe" : "F"
        }
        aminoAcid = ""
        for i in range(0,len(codons),3):
            actualCodon = codons[i:i+3]
            aminoAcid += codon_three_to_one [codon_full[codons[i:i+3]]]
        return aminoAcid

    def entropyCalc2(self,profile):
        entropy_matrix = 0.0
        for i in profile:
            print (i)
            entropy_matrix += self.entropyCalcColumn(i)
            print (entropy_matrix)
        return (entropy_matrix)


    def entropyCalc(self,profile):
        entropy_matrix = 0.0
        for i in profile:
            for element in i:
                if element != 0:
                    entropy_matrix += element*math.log(element, 2)
            print (entropy_matrix)
        return (-1*entropy_matrix)

    def entropyCalcColumn(self, column):
        entropyColumn = 0
        for element in column:
            if element != 0:
                entropyColumn += element*math.log(element, 2)
        return -1*entropyColumn

    def entropyCalcColumn_old(self, column):
        entropyColumn = 0
        for i in range(0, len(column)):
            element = column[i]
            if element != 0:
                entropyColumn += element*math.log(element, 2)
        return -1*entropyColumn

    def count(self,motifs):
        k = len(motifs[0])
        count = {}
        for symbol in "ACGT":
            count[symbol]=[]
            for j in range(k):
                count[symbol].append(0)
        t = len(motifs)
        for i in range(t):
            for j in range(k):
                symbol = motifs[i][j]
                count[symbol][j] +=1
        return count

    def profile(self,motifs):
        countMatrix = self.count(motifs)
        profileMatrix = {}
        k = len(motifs[0]) # number of columns
        t = len(motifs) # number of rows
        for symbol in "ACGT":
            profileMatrix[symbol]=[]
            for j in range(k):
                profileMatrix[symbol].append(0)
        for symbol in "ACGT":
            for j in range(k):
                profileMatrix[symbol][j]=countMatrix[symbol][j]/t
        return profileMatrix

    def consensus(self,motifs):
        consensusString = ""
        k = len(motifs[0])  # number of columns
        countMatrix = self.count(motifs)
        for j in range(k):
            m = 0
            frequentSymbol = ""
            for symbol in "ACGT":
                if countMatrix[symbol][j]> m:
                    m = countMatrix[symbol][j]
                    frequentSymbol = symbol
            consensusString += frequentSymbol
        return consensusString

    def score(self,motifs):
        consensusString = self.consensus(motifs)
        countMatrix = self.count(motifs)
        score = 0
        k = len(motifs[0])
        for symbol in "ACGT":
            for j in range (k):
                if symbol != consensusString[j]:
                    score += countMatrix[symbol][j]
        return score

    def kmerProb(self,kmer,profileMatrix):
        k = len(kmer)
        probability = 1
        for j in range(k):
            symbol = kmer[j]
            probability *= profileMatrix[symbol][j]
        return probability

    def mostProbableKmer(self,dna,k,profileMatrix):
        max = 0
        maxKmer =dna[0:k] # return first kmer if probability for all is zero
        totalLenght=len(dna)
        for i in range(len(dna)-k+1):
            kmer = dna[i:i+k]
            prob = self.kmerProb(kmer,profileMatrix)
            if prob > max:
                max = prob
                maxKmer = kmer
        return maxKmer

    def greedyMotifSearch(self,dna,k,t):
        bestMotifs = []
        for i in range (0,t):
            bestMotifs.append(dna[i][0:k]) # dna es una lista de strings (dna1, dna2, etc) y t es el largo de esa lista, con esto agrego el primer kmer de cada dna1,dna2,... de la lista
        n = len(dna[0])
        for i in range (n-k+1):
            motifs = [] # la lsita de strgins motifs
            motifs.append(dna[0][i:i+k]) # agrego primero el promer kmer despues corro una letra y agrego otro
            for j in range (1,t): #starting at 1 because i already have dna[0] with the sliding window
                p =self.profile(motifs[0:j]) # create a profile matrix first of motif 0 then motifs 0 and 1 them motif 0,1 and 2
                motifs.append(self.mostProbableKmer(dna[j],k,p))
            if self.score(motifs) < self.score(bestMotifs):
                bestMotifs = motifs
        return bestMotifs




def main():

    bio = Bioinformatics()
    profile_matrix = [
        [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
        [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
        [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
    ]
    column = [0.2, 0.6, 0.0, 0.2]
    #print (bio.entropyCalcColumn(column))
    #print (bio.entropyCalc(profile_matrix))
    # motifs = ["AACGTA","CCCGTT","CACCTT","GGATTA","TTCCGG"]
    # print (bio.count(motifs))
    # print (bio.profile(motifs))
    # print(bio.consensus(motifs))
    # print(bio.score(motifs))
    # kmer="TCGTGGATTTCC"
    # kmer2 = "ACGGGGATTACC"
    # profileMatrix = {"A":[0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
    #     "C":[0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
    #     "G":[0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
    #     "T":[0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
    # }
    # print(bio.kmerProb(kmer2,profileMatrix))
    # profileMatrix =  {
    #       'A': [0.2, 0.2, 0.3, 0.2, 0.3],
    #       'C': [0.4, 0.3, 0.1, 0.5, 0.1],
    #       'G': [0.3, 0.3, 0.5, 0.2, 0.4],
    #       'T': [0.1, 0.2, 0.1, 0.1, 0.2]
    #     }
    # A: 0.4,0.3,0.0,0.1,0.0,0.9
    #
    # C: 0.2,0.3,0.0,0.4,0.0,0.1
    #
    # G: 0.1,0.3,1.0,0.1,0.5,0.0
    #
    # T: 0.3,0.1,0.0,0.4,0.5,0.0
    # dna = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
    # k=5
    # print(bio.mostProbableKmer(dna,k,profileMatrix))
    # dna = [
    #     "GGCGTTCAGGCA",
    #     "AAGAATCAGTCA",
    #     "CAAGGAGTTCGC",
    #     "CACGTCAATCAC",
    #     "CAATAATATTCG"
    # ]
    # k=3
    # t=5
    # dna =  [
    # "GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC",
    # "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG",
    # "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC",
    # "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC",
    # "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG",
    # "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA",
    # "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA",
    # "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG",
    # "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG",
    # "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"
    # ]
    # k=15
    # t=10
    # bestMotifs=bio.greedyMotifSearch(dna,k,t)
    # score = bio.score(bestMotifs)
    # print(bestMotifs)
    # print(score)
    profileMatrix =  {
          'A': [0.4,0.3,0.0,0.1,0.0,0.9],
          'C': [0.2,0.3,0.0,0.4,0.0,0.1],
          'G': [0.1,0.3,1.0,0.1,0.5,0.0],
          'T': [0.3,0.1,0.0,0.4,0.5,0.0]
        }
    print(bio.)

    # print(bio.kmerProb("AAGTTC",profileMatrix))
    # print(bio.findAminoAcidFromCodon("CCAAGUACAGAGAUUAAC"))
    # print (bio.findAminoAcidFromCodon("CCCAGGACUGAGAUCAAU"))
    # print (bio.findAminoAcidFromCodon("CCGAGGACCGAAAUCAAC"))
    # print (bio.findAminoAcidFromCodon("CCCAGUACCGAAAUUAAC"))

if __name__ == "__main__":
    main()