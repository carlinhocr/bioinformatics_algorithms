import math


class Bioinformatics(object):
    def __init__(self):
        pass


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
            #print (element)
            if element != 0:
                entropyColumn += element*math.log(element, 2)
            #print(entropyColumn)
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
    motifs = ["AACGTA","CCCGTT","CACCTT","GGATTA","TTCCGG"]
    print (bio.count(motifs))
    print (bio.profile(motifs))
if __name__ == "__main__":
    main()