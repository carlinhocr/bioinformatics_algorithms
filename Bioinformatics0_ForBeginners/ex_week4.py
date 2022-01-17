import math


class Bioinformatics(object):


    def __init__(self):
        pass

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

    def countWithPseudocounts(self,motifs):
        k = len(motifs[0])
        count = {}
        for symbol in "ACGT":
            count[symbol]=[]
            for j in range(k):
                count[symbol].append(1) # to create the lapalce pseudocoutns intead of initializing with zeroes I do it with ones
        t = len(motifs)
        for i in range(t):
            for j in range(k):
                symbol = motifs[i][j]
                count[symbol][j] += 1
        return count

    def profileWithPseudocounts(self,motifs):
        countMatrix = self.countWithPseudocounts(motifs)
        profileMatrix = {}
        k = len(motifs[0]) # number of columns
        t = len(motifs) # number of rows
        totatDivCount = 0 # laPlace
        for symbol in "ACGT":
            totatDivCount += countMatrix[symbol][0]
            profileMatrix[symbol]=[]
            for j in range(k):
                profileMatrix[symbol].append(0)
        for symbol in "ACGT":
            for j in range(k):
                profileMatrix[symbol][j]=countMatrix[symbol][j]/totatDivCount
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

    def greedyMotifSearchWithPseudocounts(self,dna,k,t):
        bestMotifs = []
        for i in range (0,t):
            bestMotifs.append(dna[i][0:k]) # dna es una lista de strings (dna1, dna2, etc) y t es el largo de esa lista, con esto agrego el primer kmer de cada dna1,dna2,... de la lista
        n = len(dna[0])
        for i in range (n-k+1):
            motifs = [] # la lsita de strgins motifs
            motifs.append(dna[0][i:i+k]) # agrego primero el promer kmer despues corro una letra y agrego otro
            for j in range (1,t): #starting at 1 because i already have dna[0] with the sliding window
                p =self.profileWithPseudocounts(motifs[0:j]) # create a profile matrix first of motif 0 then motifs 0 and 1 them motif 0,1 and 2
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
    dna = [
        "GGCGTTCAGGCA",
        "AAGAATCAGTCA",
        "CAAGGAGTTCGC",
        "CACGTCAATCAC",
        "CAATAATATTCG"
    ]
    k=3
    t=5

    dna =  [
    "GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC",
    "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG",
    "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC",
    "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC",
    "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG",
    "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA",
    "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA",
    "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG",
    "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG",
    "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"
    ]
    k=15
    t=10
    print(bio.greedyMotifSearchWithPseudocounts(dna,k,t))
    # bestMotifs=bio.greedyMotifSearch(dna,k,t)
    # score = bio.score(bestMotifs)
    # print(bestMotifs)
    # print(score)
    # profileMatrix =  {
    #       'A': [0.4,0.3,0.0,0.1,0.0,0.9],
    #       'C': [0.2,0.3,0.0,0.4,0.0,0.1],
    #       'G': [0.1,0.3,1.0,0.1,0.5,0.0],
    #       'T': [0.3,0.1,0.0,0.4,0.5,0.0]
    #     }
    # print(bio.kmerProb("AAGTTC",profileMatrix))
    # print(bio.findAminoAcidFromCodon("CCAAGUACAGAGAUUAAC"))
    # print (bio.findAminoAcidFromCodon("CCCAGGACUGAGAUCAAU"))
    # print (bio.findAminoAcidFromCodon("CCGAGGACCGAAAUCAAC"))
    # print (bio.findAminoAcidFromCodon("CCCAGUACCGAAAUUAAC"))
    # motifs = (
    #     'AACGTA',
    #     'CCCGTT',
    #     'CACCTT',
    #     'GGATTA',
    #     'TTCCGG')
    # print(bio.countWithPseudocounts(motifs))
    # print(bio.profileWithPseudocounts(motifs))

if __name__ == "__main__":
    main()