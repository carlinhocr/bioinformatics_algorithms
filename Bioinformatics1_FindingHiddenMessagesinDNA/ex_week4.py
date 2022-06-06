import math
import random
from scipy.stats import binom

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

    def kmerProb(self, kmer, profileMatrix):
        k = len(kmer)
        probability = 1
        for j in range(k):
            symbol = kmer[j]
            probability *= profileMatrix[symbol][j]
        return probability

    def mostProbableKmer(self,dna,k,profileMatrix):
        maxProb = 0
        maxKmer = dna[0:k]  # return first kmer if probability for all is zero
        for i in range(len(dna)-k+1):  # total Length=len(dna)
            kmer = dna[i:i+k]
            prob = self.kmerProb(kmer, profileMatrix)
            if prob > maxProb:
                maxProb = prob
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

    def randomizedMotifSearchNTimes(self, dna, k, t, n):
        m_initial = self.randomMotifs(dna, k, t)  # only to have an initial motif set
        bestMotifs = m_initial
        for i in range(0, n):
            m = self.randomizedMotifSearch(dna, k, t)   # a full run of randomized motif search that finishes with a score that is not less than the last one
            if self.scorewithPseudocounts(m) < self.scorewithPseudocounts(bestMotifs):
                bestMotifs = m[:]
        return bestMotifs

    def probabilityAllKmersFoundNStrings(self, numberOfKmers, kmersLenght, nucleotidesLengh):
        # probability of not capturing the  kmer in the first nucleotide string:
        probNotFoundFirstNucleotide = (nucleotidesLengh - kmersLenght) / (nucleotidesLengh - kmersLenght + 1)
        probAllKmersNotFound = pow(probNotFoundFirstNucleotide, numberOfKmers)
        probAllKmersFound = 1 - probAllKmersNotFound
        return probAllKmersFound

    def probabilityOfAtLeast2Found(self, numberOfKmers, kmersLenght, nucleotidesLengh):
        # This is a binomial probability question. Think about the fact that you're flipping a weighted coin 10 times and you
        # want to calculate the probability of flipping heads at least two times. This is the same as the probability of NOT flipping heads either zero times or one time, which can be calculated as 1 - P(H=0) - P(H=1).
        # There are 600 - 15 + 1 positions in a single string, of which a single one is correct. That suggests our coin is
        # weighted such that the probability of heads is 1/586.
        # The probability of zero heads, P(H=0), is what was calculated previously. That means the 10 coin flips all
        # came up tails, at a probability of 585/586 each time.
        # Now we have to calculate the probability of heads coming up only once. That is the probability of 9 tails and 1 head.
        # BUT, that head flip can come in any order, of which there are 10 possible positions (10 choose 1 is simply 10).
        # So after you multiply the independent probabilities (9 of the tails and one of the heads), you have to multiply
        # that by the 'binomial coefficient', which in this case is 10 for the ten possible ways those flips can be ordered.
        # This will give you P(H=1).
        # Now we have P(H=0) and P(H=1), but we want the probability that this did NOT happen,
        # so the full answer is 1 - P(H=0) - P(H=1).
        # You
        # can use scipy.stats module
        # with cumulative distribution function in binom
        # p(success) = 1 / 586
        # number
        # of
        # trials = 10
        # binom.cdf(1, 10, 1 / 586) is the
        # probability
        # of
        # 1 or 0
        # success
        # event.
        # from scipy.stats import binom
        probOneFound = 1 / (nucleotidesLengh - kmersLenght + 1)
        ProbabilityOneOrZeroEvents = 1 - binom.cdf(1, numberOfKmers, probOneFound)
        return ProbabilityOneOrZeroEvents



# ==============================

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

    def count(self, motifs):
        k = len(motifs[0])
        count = {}
        for symbol in "ACGT":
            count[symbol] = []
            for j in range(k):
                count[symbol].append(0)
        t = len(motifs)
        for i in range(t):
            for j in range(k):
                symbol = motifs[i][j]
                count[symbol][j] += 1
        return count
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
    def normalize(self,probabilities):
        sum = 0
        for element in probabilities.values():
            sum += element
        for element2 in probabilities.keys():
            probabilities[element2]=probabilities[element2]/sum
        return probabilities


    def weightedDie (self, probabilities):
        kmer = ''
        number = random.uniform(0,1)
        acum = 0
        for element in probabilities.keys():
            if number >= acum and number < acum + probabilities[element]:
                kmer = element
            acum += probabilities[element]
        return kmer

    def profileGeneratedString(self, Text, profile, k):
        n = len(Text)
        probabilities = {}
        probabilities_normalized = {}
        probabilities_weighted = {}
        for i in range(0,n-k+1):
            kmer = Text[i:i+k]
            probabilities[kmer] = self.kmerProb(kmer,profile)
        probabilities_normalized = self.normalize(probabilities)
        string_probabilities_weighted = self.weightedDie(probabilities_normalized)
        return string_probabilities_weighted

    def cloning(self, li1):
        li_copy = li1[:]
        return li_copy

    def gibbsSampler(self,dna, k, t, N):
        bestMotifs = []
        motifs = self.randomMotifs(dna,k,t)
        bestMotifs = motifs
        for j in range(1, N):
            i = random.randint(0, t-1)
            motifs_without_i= self.cloning(motifs)
            motifs_selected = motifs[i]
            motifs_without_i.pop(i)
            profile = self.profileWithPseudocounts(motifs_without_i)
            dna_i = dna[i]
            kmer_profile_weighted = self.profileGeneratedString(dna_i,profile,k)
            motifs.pop(i) # it is more clear to use an auxliary list motifs_without_i for debugging
            motifs.insert(i, kmer_profile_weighted)
            if self.score(motifs) < self.score(bestMotifs):
                bestMotifs = motifs
        return bestMotifs

    def gibbsSamplerRepeatXtimes(self,dna,k,t,n,x):
        motifs = self.randomMotifs(dna,k,t)
        bestMotifs = motifs
        for i in range(0,x):
            m = self.gibbsSampler(dna,k,t,n)
            if self.score(m) < self.score(bestMotifs):
                bestMotifs = m
        return bestMotifs


def main():

    bio = Bioinformatics()
    # filename = 'dataset_161_5.txt'
    # with open(filename, "r") as dataset:
    #     data = []
    #     for line in dataset:
    #         data.append(line.strip())
    #     k_string,t_string = data[0].split()
    #     k = int(k_string)
    #     t = int(t_string)
    #     dna_sequences = data[1].strip().split() # several motifs
    #     n = 1000
    #     # print (dna_sequences,k,t,n)
    #     list = bio.randomizedMotifSearchNTimes(dna_sequences,k,t,n)
    #     print(*list) # separate the elements with one space only (using *
    # ---------------------------------------------------------------------------------------------------------------
    # print(bio.probabilityAllKmersFoundNStrings(numberOfKmers=10,kmersLenght=15,nucleotidesLengh=600))
    # print(bio.probabilityOfAtLeast2Found(numberOfKmers=10, kmersLenght=15, nucleotidesLengh=600))
    # k = 8
    # t = 5
    # n = 100
    # dna_sequences = ['CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA',
    # 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
    # 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
    # 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
    # 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
    # list = bio.gibbsSampler(dna_sequences,k,t,n)
    # print(*list)
    filename = 'dataset_163_4.txt'
    with open(filename, "r") as dataset:
        data = []
        for line in dataset:
            data.append(line.strip())
        k_string,t_string,n_string = data[0].split()
        k = int(k_string)
        t = int(t_string)
        n = int(n_string)
        dna_sequences = data[1].strip().split() # several motifs
        print(dna_sequences,k,t,n)
        random_starts = 20
        lista = bio.gibbsSamplerRepeatXtimes(dna_sequences,k,t,n,random_starts)
        print(*lista) # separate the elements with one space only using *
    # ---------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()