import math


class Bioinformatics(object):


    def __init__(self):
        pass

    def cloningList(li1):
        li_copy = li1[:]
        return li_copy

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

    def consensusFromProfile(self,probMatrix):
        numColumns = len(probMatrix["A"])
        consensusMax={}
        for i in range(numColumns):
            localMax = 0
            for symbol in "ACGT":
                if probMatrix[symbol][i] > localMax:
                    localMax = probMatrix[symbol][i]
            consensusMax[i] = []
            for symbol in "ACGT":
                if probMatrix[symbol][i] == localMax:
                    consensusMax[i] += symbol
        return consensusMax


    def stringListInStringList(self,listElements,listToSearch):
        found = []
        notFound = []
        for element in listElements:
            if element in listToSearch:
                found.append(element)
            else:
                notFound.append(element)
        return found, notFound

    def buildString(self,list):
        n = len(list)
        letter = []
        first_run = True
        count_letter = 0
        for x in range(len(list),0,-1):
            letter.append([]) # initialized letter
            if first_run == True:
                for element in list[x-1]:
                    letter[0].append(element)
                first_run = False
            else:
                for element in list[x-1]:
                    for elementString in letter[count_letter]:
                        letter[count_letter+1].append(element+elementString)
                count_letter += 1
        return (letter[n-1])

    def move_cursor_ale(self,cur_item, cur_indexes, max_indexes, len_lista):
        found = False
        while not found:
            if cur_indexes[cur_item] < max_indexes[cur_item]:
                cur_indexes[cur_item] += 1
                for i in range(cur_item + 1, len_lista):
                    cur_indexes[i] = 0
                cur_item = len_lista - 1
                found = True
            else:
                cur_item -= 1
        return cur_item

    def get_join_ale(self,lista):
        len_lista = len(lista)
        cur_indexes = [0] * len_lista
        max_indexes = [len(subl) - 1 for subl in lista]
        cur_item = len_lista - 1

        final = []

        while cur_indexes != max_indexes:
            final.append("".join([lista[i][cur_indexes[i]] for i in range(0, len_lista)]))
            cur_item = self.move_cursor_ale(cur_item, cur_indexes, max_indexes, len_lista)

        final.append("".join([lista[i][cur_indexes[i]] for i in range(0, len_lista)]))
        return final

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
#------------------------------------------------------


    def profileAsList(self,motifs):
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
        profileMatrixList = []
        for numbers in profileMatrix.values():
            profileMatrixList.append(numbers)
        return profileMatrixList

    def probabiltyKmers(self, k, strechesArrays, bases):
        likelihood = 1 / (4 ** k)  # number of nucleotides x kmer * each nucleitide acgt if they are all equally probable
        sequences = bases - k + 1  # number of nucleotides of lenght k in bases number of nucloetides
        oneStrechProbability = sequences * likelihood
        allStrechesProbability = oneStrechProbability * strechesArrays
        print(allStrechesProbability)
        return allStrechesProbability


    def inmediateNeighbors(self,pattern):
        neighborhood = []
        neighborhood.append(pattern)
        nucleotides = ["A","C","G","T"]
        for i in range(len(pattern)):
            symbol= pattern[i]
            for nucleo in nucleotides:
                if nucleo != symbol:
                    if i == 0:
                        neighbor = nucleo + pattern[i+1:]
                    elif i == len(pattern)-1:
                        neighbor = pattern[0:i]+nucleo
                    else:
                        neighbor = pattern[0:i]+nucleo+pattern[i+1:]
                    neighborhood.append(neighbor)
        return neighborhood

    def iterativeNeighbors(self,pattern,d):
        neighborhood = [pattern]
        result = []
        if d == 0:
            pass
        else:
            for j in range(1,d+1):
                for neighbor in self.iterativeNeighbors(pattern,j-1):
                    neighborhood = neighborhood + self.inmediateNeighbors(neighbor)
        for element in neighborhood:
            if element not in result:
                result.append(element)
        result.sort()
        return result

    def hammingDistance(self,p,q):
        if len(p) >= len(q):
            mayor = p
            menor = q
        else:
            mayor = q
            menor = p
        count = 0
        for i in range (0,len(menor)):
            if p[i] != q[i]:
                count += 1
        count += len(mayor)-len(menor)
        return count

    def motifEnumeration(self,dna,k,d):
        patterns = []
        encontrado = {}
        patternNeighboors = {}
        for dna_string in dna:
            n = len(dna_string)
            for i in range(n-k+1):
                pattern = dna_string[i:i+k]
                patternNeighboors[pattern] = self.iterativeNeighbors(pattern,d)
                encontrado[pattern]=[]
                for patternNeig in patternNeighboors[pattern]:
                    count = 0
                    for dna_string_2 in dna:
                        for j in range(n-k+1):
                            distance = self.hammingDistance(patternNeig,dna_string_2[j:j+k])
                            if distance<=d:
                                count += 1
                                break # so I only count once per dna
                    if count >= len(dna):
                        patterns.append(patternNeig)
        return list(dict.fromkeys(patterns))

    def maximumScore(self,t,k):
        nucloetides ="ACGT"
        maximumRights = math.ceil(t/len(nucloetides)) # rounds to the next integer
        print (maximumRights)
        maximunMismatchesperColumn = t - maximumRights
        totalMistmatechs = maximunMismatchesperColumn * k
        print(totalMistmatechs)

    def distancePatternAndStrings(self,pattern,dna):
        k = len(pattern)
        distance = 0
        for motif in dna:
            hammingDistance = len(pattern)+len(motif) #just like making it a very big number
            for i in range(len(motif)-k+1):
                pattern2 = motif[i:i+k]
                distancePatterns = self.hammingDistance(pattern,pattern2)
                if hammingDistance > distancePatterns:
                    hammingDistance = distancePatterns
            distance += hammingDistance
        return distance

    def allStrings(self,k):
        allStringsLenghtK = []
        pattern =""
        for i in range(k):
            pattern += "A"
        allStringsLenghtK = self.iterativeNeighbors(pattern,k)
        return allStringsLenghtK

    def medianString(self,dna,k):
        median = ""
        distance = k + len(dna[0]) # like settign distance to infinity or a maximum value you can use also sys.maxsize
        patternNeighbours = self.allStrings(k)
        for pattern in patternNeighbours:
            distancePatternsDNA = self.distancePatternAndStrings(pattern, dna)
            if distance > distancePatternsDNA:
                distance = distancePatternsDNA
                median = pattern
        return median

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
    #print(bio.allStrings(4))
    # k = 3
    # t = 5
    # dna =["GGCGTTCAGGCA","AAGAATCAGTCA","CAAGGAGTTCGC","CACGTCAATCAC","CAATAATATTCG"]
    # filename = 'dataset_160_9.txt'
    # with open(filename, "r") as dataset:
    #     data = []
    #     for line in dataset:
    #         data.append(line.strip())
    #     k_string,t_string = data[0].split()
    #     k = int(k_string)
    #     t = int(t_string)
    #     dna_sequences = data[1].strip().split() # several motifs
    # list = bio.greedyMotifSearchWithPseudocounts(dna_sequences,k,t)
    # for element in list:
    #     print(element)
    # filename = 'dataset_159_3.txt'

    # with open(filename, "r") as dataset:
    #     data = []
    #     for line in dataset:
    #         data.append(line.strip())
    #     dna = data[0]
    #     k = int(data[1])
    #     raw_profile = data[2:]
    #     bases = ['A', 'C', 'G', 'T']
    #     prof = [list(map(float, raw_profile[i].split())) for i in range(len(raw_profile))]
    #     prof_dict = dict(zip(bases, prof))
    # print(bio.mostProbableKmer(dna,k,prof_dict))
    # profile = {
    #     'A': [0.25,0.289,0.25,0.224,0.276,0.224,0.303,0.316,0.289,0.197,0.276,0.237,0.224],
    #     'C': [0.342,0.303,0.382,0.316,0.342,0.316,0.289,0.263,0.237,0.382,0.355,0.211,0.263],
    #     'G': [0.197,0.197,0.197,0.158,0.158,0.237,0.132,0.289,0.263,0.224,0.197,0.276,0.263],
    #     'T': [0.211,0.211,0.171,0.303,0.224,0.224,0.276,0.132,0.211,0.197,0.171,0.276,0.25]
    # }
    #
    # dna = "AATAGACCGCGATCGAAGCCTATTATCTTCAGGGGAGCAATAAGGCTCCATGCCCTCCTTTATAATATGTACTCCCTATGAAAGCCTTGGACCAACCACACGCCCCGTGCGCCCCGTTTACAGACTAGAAGGAGGTGACCTTTTCGGCCGATGTCTCAGGAGGCGGACACGCAATACCGGGGCTGAGCTTAAATGCGCAGATTAGTGTGGTTCGGATGGTGTTGAGCTGCTGGGTTTGGAGGGACTATGCCACCTGGGAGTGACTGAAGACACACCGACGTCCTCCGAGTGGGAACTACGCGTCAGATACCCCACATGTGATGGGACAGTTATTCGCCAAACACGCCGGCCGCGCCAGGAAACCCACCTCACTCGAAAGTAGGGTACTGCTCAAGCCCATGGCTATCTTTCATTCAGCTTCTACTTGTAATCCATTCGCCATATTTGGAATAACTCGCACGCGCAGGGGTGACGGCTGCGTTAAAGACGGTAGACATTCCATAGAATAAACTCAATCGAGTGCTATGCGGTAGGACGGGCCGTTTAGTGGACTCCCATATAATATGATATTTTAGTCGCTTCATGTTCCAACGAGCGTCAGCGCTGATAACTTCGCAGATACAACAGACAAGTCTGGGAAAAACATAGATGAAAGATATAGGTGTAACCGAATTCCCATGGCACGGCTAGTTTGGGTTTCCCATTGTCCGTTAAGCAATAGTCTAAGAGATTATCAATTTCCGGTGAGCCTGAGTTACCTTCCCATGTGCACCGAAAGGCTTAAACGCCCACCTTAGAAGGCGGGTCATATATGTGGAGGCCAATATGAATGCTAGGGCGGATCGGCCAATCATGCTACAGAGCTTAGTCCACATATGAGCTTAAGCGGTGGCTTATGCCCGCATCTGCACTCTAAGGGGCCGTCGCTACCAGACGGGTGGTTATTAGTGCTGTCGAGGTTTGCGGGTACAATACCAGAACCAGTTCCGGCTGTCGCCAA"
    # k = 13
    # print(bio.mostProbableKmer(dna,k,profile))
    # profile = {
    #
    #     'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
    #     'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
    #     'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
    #     'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
    # }
    # pattern = "TCGGGGATTTCC"
    # print (bio.kmerProb(pattern,profile))
    # with open("dataset_158_9.txt") as inp:
    #     input_items = inp.read().strip().splitlines()
    #     k_string = input_items[0].strip()
    #     k = int(k_string)
    #     dna_sequences = input_items[1].strip().split() # several motifs
    # print(bio.medianString(dna_sequences,k))
    # dna = ['TTACCTTAAC', 'GATATCTGTC', 'ACGGCGTTCG', 'CCCTAAAGAG', 'CGTCAGAGGT']
    # pattern = "AAA"
    # with open("dataset_5164_1.txt") as inp:
    #     input_items = inp.read().strip().splitlines()
    #     pattern = input_items[0].strip()
    #     dna_sequences = input_items[1].strip().split() # several motifs
    # print(bio.distancePatternAndStrings(pattern,dna_sequences))
    #bio.probabiltyKmers(9,500,1000)
    # dna = [
    #     "ATTTGGC",
    #     "TGCCTTA",
    #     "CGGTATC",
    #     "GAAAATT",
    # ]
    # k=3
    # d=1
    # dna = [
    #     "TGATGTTCACCACTAGCTCAAAGCT",
    #     "TGAAAGCCGATCACTTGATAAGCGA",
    #     "TAAGCCGTTGTGATTCGTAATTATC",
    #     "CGGTCTGATAACAGTCGTCTAATAA",
    #     "TGATTGAGATCATAGAGTATACCCA",
    #     "TGTGCCCTGATGATATACCGGTAAA",
    # ]
    # k=5
    # d=1
    # print(bio.motifEnumeration(dna,k,d))
    #bio.maximumScore(10,15)

    # motifs = [
    #     "TCGGGGGTTTTT",
    #     "CCGGTGACTTAC",
    #     "ACGGGGATTTTC",
    #     "TTGGGGACTTTT",
    #     "AAGGGGACTTCC",
    #     "TTGGGGACTTCC",
    #     "TCGGGGATTCAT",
    #     "TCGGGGATTCCT",
    #     "TAGGGGAACTAC",
    #     "TCGGGTATAACC"
    # ]
    # profileMAtrix = bio.profileAsList(motifs)
    # print (profileMAtrix)
    # entropy = bio.entropyCalc(profileMAtrix)
    # print(entropy)



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
    #print(bio.kmerProb("CAGTGA",profileMatrix))
    consensusList = bio.consensusFromProfile(profileMatrix)
    print(consensusList)
    print("Carlos",bio.buildString(consensusList))
    print ("Ale   ", bio.get_join_ale(list(consensusList.values())))
    totalList = bio.buildString(consensusList)
    tofind =['AAGTGA', 'AAGTTA', 'ACGCGA','AATTGA', 'ATTTTA', 'TCGCGA']
    found, notfound = bio.stringListInStringList(tofind,totalList)
    print("Found: ",found)
    print("Not Found", notfound)

    # print(bio.findAminoAcidFromCodon("CCAAGUACAGAGAUUAAC"))
    # print (bio.findAminoAcidFromCodon("CCCAGGACUGAGAUCAAU"))
    # print (bio.findAminoAcidFromCodon("CCGAGGACCGAAAUCAAC"))
    # print (bio.findAminoAcidFromCodon("CCCAGUACCGAAAUUAAC"))
    # dna = ["CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC",
    #        "GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC",
    #        "GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"]
    # print(bio.medianString(dna,7))
    # a =[0.5,0,0,0.5]
    # b = [0.25,0.25,0.25,0.25]
    # c=[0,0,0,1]
    # d=[0.25,0,0.5,0.25]
    # print(bio.entropyCalcColumn(a), "A")
    # print(bio.entropyCalcColumn(b), "B")
    # print(bio.entropyCalcColumn(c), "C")
    # print(bio.entropyCalcColumn(d), "D")
    # motifs = ["AAGAGA","AAGCGA","AAGCGA","ACGCGA","CCGCGA","CCGGTA","GGGTTA","TGGTTA","TGGTTA","TTGTTA"]
    #
    # Consider
    # the
    # following
    # motif
    # matrix:
    #
    # CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC
    #
    # GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC
    #
    # GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG
    #
    # Which
    # of
    # the
    # following
    # 7 - mers is a
    # median
    # string
    # for this motif matrix? (Select all that apply.)
    #
    # 2
    # points
    #
    # TCTGAAG
    #
    # CGTGTAA
    #
    # GTAGGAA
    #
    # GATGAGT
    #
    # ATAACGG
    #
    # TAGTTTC
    #print(bio.consensus(motifs))
if __name__ == "__main__":
    main()