class Bioinformatics(object):
    def __init__(self):
        pass

    def skewArray_v2(self,Genome):
        skew = {}
        skew[0] = 0
        for i in range(0,len(Genome)):
            if Genome[i] == "C":
                skew[i+1]=skew[i]-1
            elif Genome[i] == "G":
                skew[i+1]=skew[i]+1
            elif Genome[i] == "A":
                skew[i+1]=skew[i]
            elif Genome[i] == "T":
                skew[i+1]=skew[i]
        return list(skew.values()) #returns a list of the values of the dict

    def skewArray(self,Genome):
        skew = []
        skew.append(0)
        for i in range(0,len(Genome)):
            if Genome[i] == "C":
                skew.append(skew[i]-1)
            elif Genome[i] == "G":
                skew.append(skew[i]+1)
            elif Genome[i] == "A":
                skew.append(skew[i])
            elif Genome[i] == "T":
                skew.append(skew[i])
        return skew

    def MinimumSkew(self,Genome):
        positions =[]
        skew = self.skewArray_v2(Genome)
        min_skew = min(skew)
        print (min_skew)
        for i in range (0,len(skew)):
            if skew[i] == min_skew:
                positions.append(i)
        return positions

#-------------------------------------------------------------------------------


def main():

    bio = Bioinformatics()
    # genome = "GAGCCACCGCGATA"
    # print(bio.skewArray_v2(genome).values())
    with open('dataset_7_10.txt','r') as file7_10:
        text = file7_10.read()
    print(bio.MinimumSkew(text))



if __name__ == "__main__":
    main()