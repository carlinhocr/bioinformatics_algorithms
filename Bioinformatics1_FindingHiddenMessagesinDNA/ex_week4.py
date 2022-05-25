import math
import random

# comentario de test loco
class Bioinformatics(object):


    def __init__(self):
        pass

    def randomMotifs(self,dna,k,t):
        randomizedMotifs=[]
        for row in dna:
            startNum = random.randint(0,len(row)-1)
            if startNum + k < len(row):
                randomizedMotifs.append(row[startNum:startNum+k])
            else:
                randomizedMotifs.append(row[startNum-k:startNum])
        return randomizedMotifs

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

    def motif(self, profileMatrix, dna,k=4):
        mostProbableKmerDna = []
        for row in dna:
            mostProbableKmerDna.append(self.mostProbableKmer(row,k,profileMatrix))
        return mostProbableKmerDna

    def consensus(self,motifs):
        consensusString = ""
        k = len(motifs[0])  # number of columns
        countMatrix = self.countWithPseudocounts(motifs)
        for j in range(k):
            m = 0
            frequentSymbol = ""
            for symbol in "ACGT":
                if countMatrix[symbol][j]> m:
                    m = countMatrix[symbol][j]
                    frequentSymbol = symbol
            consensusString += frequentSymbol
        return consensusString

    def scorewithPseudocounts(self,motifs):
        consensusString = self.consensus(motifs)
        countMatrix = self.countWithPseudocounts(motifs)
        score = 0
        k = len(motifs[0])
        for symbol in "ACGT":
            for j in range (k):
                if symbol != consensusString[j]:
                    score += countMatrix[symbol][j]
        return score

    def randomizedMotifSearch(self, dna, k, t):
        m = self.randomMotifs(dna, k, t)
        bestMotifs = m
        while True:
            profile = self.profileWithPseudocounts(m)
            m = self.motif(profile, dna, k)
            if self.scorewithPseudocounts(m) < self.scorewithPseudocounts(bestMotifs):
                bestMotifs = m
            else:
                return bestMotifs

    def randomizedMotifSearchNtimes(self,dna,k,t,n):
        m= self.randomMotifs(dna,k,t)
        bestMotifs = m
        for i in range(0,n):
            m = self.randomizedMotifSearch(dna,k,t)
            if self.scorewithPseudocounts(m) < self.scorewithPseudocounts(bestMotifs):
                bestMotifs = m[:]
        return bestMotifs


# ==============================

    # def score(self,motifs):
    #     consensusString = self.consensus(motifs)
    #     countMatrix = self.count(motifs)
    #     score = 0
    #     k = len(motifs[0])
    #     for symbol in "ACGT":
    #         for j in range (k):
    #             if symbol != consensusString[j]:
    #                 score += countMatrix[symbol][j]
    #     return score
    #
    # def count(self, motifs):
    #     k = len(motifs[0])
    #     count = {}
    #     for symbol in "ACGT":
    #         count[symbol] = []
    #         for j in range(k):
    #             count[symbol].append(0)
    #     t = len(motifs)
    #     for i in range(t):
    #         for j in range(k):
    #             symbol = motifs[i][j]
    #             count[symbol][j] += 1
    #     return count
    #
    #
    #
    #
    # def greedyMotifSearch(self,dna,k,t):
    #     bestMotifs = []
    #     for i in range (0,t):
    #         bestMotifs.append(dna[i][0:k]) # dna es una lista de strings (dna1, dna2, etc) y t es el largo de esa lista, con esto agrego el primer kmer de cada dna1,dna2,... de la lista
    #     n = len(dna[0])
    #     for i in range (n-k+1):
    #         motifs = [] # la lsita de strgins motifs
    #         motifs.append(dna[0][i:i+k]) # agrego primero el promer kmer despues corro una letra y agrego otro
    #         for j in range (1,t): #starting at 1 because i already have dna[0] with the sliding window
    #             p =self.profile(motifs[0:j]) # create a profile matrix first of motif 0 then motifs 0 and 1 them motif 0,1 and 2
    #             motifs.append(self.mostProbableKmer(dna[j],k,p))
    #         if self.score(motifs) < self.score(bestMotifs):
    #             bestMotifs = motifs
    #     return bestMotifs
    #
    # def greedyMotifSearchWithPseudocounts(self,dna,k,t):
    #     bestMotifs = []
    #     for i in range (0,t):
    #         bestMotifs.append(dna[i][0:k]) # dna es una lista de strings (dna1, dna2, etc) y t es el largo de esa lista, con esto agrego el primer kmer de cada dna1,dna2,... de la lista
    #     n = len(dna[0])
    #     for i in range (n-k+1):
    #         motifs = [] # la lsita de strgins motifs
    #         motifs.append(dna[0][i:i+k]) # agrego primero el promer kmer despues corro una letra y agrego otro
    #         for j in range (1,t): #starting at 1 because i already have dna[0] with the sliding window
    #             p =self.profileWithPseudocounts(motifs[0:j]) # create a profile matrix first of motif 0 then motifs 0 and 1 them motif 0,1 and 2
    #             motifs.append(self.mostProbableKmer(dna[j],k,p))
    #         if self.score(motifs) < self.score(bestMotifs):
    #             bestMotifs = motifs
    #     return bestMotifs
    #
    # def randomizedMotifSearch_exam(self):
    #     dna = ['AAGCCAAA','AATCCTGG','GCTACTTG','ATGTTTTG']
    #     k=3
    #     m =['CCA','CCT','CTT','TTG']
    #     profile = self.profileWithPseudocounts(m)
    #     m = self.motif(profile, dna, k)
    #     print (m)
    #
    # def normalize(self,probabilities):
    #     sum = 0
    #     for element in probabilities.values():
    #         sum += element
    #     for element2 in probabilities.keys():
    #         probabilities[element2]=probabilities[element2]/sum
    #     return probabilities
    #
    # def weightedDie(self,probabilities):
    #     kmer = ''
    #     number = random.uniform(0,1)
    #     acum = 0
    #     for element in probabilities.keys():
    #         if number >= acum and number < acum + probabilities[element]:
    #             kmer = element
    #         acum += probabilities[element]
    #     return kmer
    #
    # def profileGeneratedString(self, Text, profile, k):
    #     n = len(Text)
    #     probabilities = {}
    #     probabilities_normalized = {}
    #     probabilities_weighted = {}
    #     for i in range(0,n-k+1):
    #         kmer = Text[i:i+k]
    #         probabilities[kmer] = self.kmerProb(kmer,profile)
    #     probabilities_normalized = self.normalize(probabilities)
    #     string_probabilities_weighted = self.weightedDie(probabilities_normalized)
    #     return string_probabilities_weighted
    #
    # def cloning(self, li1):
    #     li_copy = li1[:]
    #     return li_copy
    #
    # def gibbsSampler(self,dna, k, t, N):
    #     bestMotifs = []
    #     motifs = self.randomMotifs(dna,k,t)
    #     bestMotifs = motifs
    #     for j in range(1, N):
    #         i = random.randint(0, t-1)
    #         motifs_without_i= self.cloning(motifs)
    #         motifs_selected = motifs[i]
    #         motifs_without_i.pop(i)
    #         profile = self.profileWithPseudocounts(motifs_without_i)
    #         dna_i = dna[i]
    #         kmer_profile_weighted = self.profileGeneratedString(dna_i,profile,k)
    #         motifs.pop(i) # it is more clear to use an auxliary list motifs_without_i for debugging
    #         motifs.insert(i, kmer_profile_weighted)
    #         if self.score(motifs) < self.score(bestMotifs):
    #             bestMotifs = motifs
    #     return bestMotifs
    #
    # def gibbsSamplerRepeatXtimes(self,dna,k,t,n,x):
    #     motifs = self.randomMotifs(dna,k,t)
    #     bestMotifs = motifs
    #     for i in range(0,x):
    #         m = self.gibbsSampler(dna,k,t,n)
    #         if self.score(m) < self.score(bestMotifs):
    #             bestMotifs = m
    #     return bestMotifs
    #

def main():

    bio = Bioinformatics()
    # #bio.randomizedMotifSearch_exam()
    # probabilities = {'A': 0.22, 'C': 0.54, 'G': 0.58, 'T': 0.36}
    # print (bio.normalize(probabilities))
    # probabilites = {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1}
    # print(bio.normalize(probabilites))
    # probabilites = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
    # print(bio.weightedDie(probabilites))
    # k = 2
    # Text = 'AAACCCAAACCC'
    # profile = {'A': [0.5, 0.1], 'C': [0.3, 0.2], 'G': [0.2, 0.4], 'T': [0.0, 0.3]}
    # print(bio.profileGeneratedString(Text,profile,k))
    # profile_matrix = [
    #     [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
    #     [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
    #     [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
    #     [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
    # ]
    # column = [0.2, 0.6, 0.0, 0.2]
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
    #
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
    # x=20
    # n=100
    # print(bio.gibbsSamplerRepeatXtimes(dna,k,t,n,x))
    # print(bio.greedyMotifSearchWithPseudocounts(dna,k,t))
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
    # profileMAtrix= {'A': [0.8, 0.0, 0.0, 0.2], 'C': [0.0, 0.6, 0.2, 0.0], 'G': [0.2, 0.2, 0.8, 0.0],
    #                     'T': [0.0, 0.2, 0.0, 0.8]}
    # dna = ['TTACCTTAAC', 'GATGTCTGTC', 'ACGGCGTTAG', 'CCCTAACGAG', 'CGTCAGAGGT']
    #print(bio.motif(profileMAtrix,dna))
    # k = 3
    # t = 5
    # print(bio.randomMotifs(dna,k,t))
    # dna = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
    #        'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
    # k = 8
    # t = 5
    # n = 100
    # print(bio.gibbsSampler(dna,k,t,n))
    # # print(bio.randomizedMotifSearchNtimes(dna,k,t,100))
    # profileMatrix =  {
    #       'A': [0.4,0.3,0.0,0.1,0.0,0.9],
    #       'C': [0.2,0.3,0.0,0.4,0.0,0.1],
    #       'G': [0.1,0.3,1.0,0.1,0.5,0.0],
    #       'T': [0.3,0.1,0.0,0.4,0.5,0.0]
    #     }
    # ---------------------------------------------------------------------------------------------------------------
    # Exercise 1 Randomized motif search n time
    # Debug Set
    # dna = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
    #        'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
    # k = 8
    # t = 5
    # n = 1000
    # print(bio.randomizedMotifSearchNtimes(dna,k,t,n))
    #File real exercise set
    filename = 'dataset_161_5.txt'
    with open(filename, "r") as dataset:
        data = []
        for line in dataset:
            data.append(line.strip())
        k_string,t_string = data[0].split()
        k = int(k_string)
        t = int(t_string)
        dna_sequences = data[1].strip().split() # several motifs
        n = 1000
        # print (dna_sequences,k,t,n)
        list = bio.randomizedMotifSearchNtimes(dna_sequences,k,t,n)
        print(*list) # separate the elements with one space only (using *
    # ---------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    main()