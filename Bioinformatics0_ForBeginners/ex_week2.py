class Bioinformatics(object):
    def __init__(self):
        pass

    def PatternCount(self, text, pattern):
        count = 0
        for i in range(0,(len(text)-len(pattern))+1):
            subtext = text[i:len(pattern)+i]
            if subtext == pattern:
                count +=1
        return count

    def SymbolArray(self, Genome, symbol):
        array = {}
        n = len(Genome)
        ExtendedGenome = Genome + Genome[0:n//2]
        for i in range(n):
            array[i] = self.PatternCount(ExtendedGenome[i:i+(n//2)],symbol)
        return array

    def FasterSymbolArray(self, Genome, symbol):
        array = {}
        n = len(Genome)
        ExtendedGenome = Genome + Genome[0:n//2]

        # look at the first half of Genome to compute first array value
        array[0] = self.PatternCount(Genome[0:n//2],symbol)

        for i in range(1, n):
            # start by setting the current array value equal to the previous array value
            array[i] = array[i-1]
            # the current array value can differ from the previous array value by at most 1
            if ExtendedGenome[i-1] == symbol:
                array[i] = array[i]-1
            if ExtendedGenome[i+(n//2)-1] == symbol:
                array[i] = array[i]+1
        return array

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
        return skew

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
        skew = self.skewArray(Genome)
        min_skew = min(skew)
        print (min_skew)
        for i in range (0,len(skew)):
            if skew[i] == min_skew:
                positions.append(i)
        return positions

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

    def approximatePatternMatching(self, genome, pattern, d):
        positions = []  # output variable
        for i in range(0, (len(genome) - len(pattern)) + 1):
            subtext = genome[i:len(pattern) + i]
            if self.hammingDistance(subtext,pattern) <=d:
                positions.append(i)
        return positions

    def approximatePatternCount(self, text, pattern,d):
        count = 0
        for i in range(0,(len(text)-len(pattern))+1):
            subtext = text[i:len(pattern)+i]
            if self.hammingDistance(subtext,pattern) <=d:
                count +=1
        return count

def main():

    bio = Bioinformatics()
    text = "AAATTAGAGTACCAAATTAGAAATTAGAACAAATTAGTCAGGAGTCGTTTGTTCGAAAATTAGCTAAATTAGATTTCAAATTAGAAATTAGGAAATTAGTCAAATTAGGGAAATTAGCCAAAAATTAGAAAATTAGAAATTAGAGAAATTAGAAAAATTAGAGATAATTAGAAATTAGCGTGGCAAATTAGGAAATTAGAAATTAGATTCTAAATTAGAAATTAGGAAATTAGAAATTAGATGAAATTAGTAAATTAGCAAATTAGGAAATTAGAAATTAGAAATTAGAAATTAGCGTATAAATTAGGAAATTAGTTATATGAGAAATTAGCCACAGTAAAATTAGCGAAATTAGGAAATTAGACAAGCTTGACGCGCAAATTAGTTAAATTAGAATAAAATTAGGTGCTCTGAAGTGCACTAAATTAGTGAAAAATTAGGAGAAATTAGGAAAATTAGTAAATTAGCAAATTAGGGCAAATTAGAAATTAGATCTGCTAAATTAGTAACCAAATTAGAAATTAGCAACAAATTAGAAATTAGTAAATTAGTATCGTTAAATTAGAAATTAGTTAAATTAGCAAATTAGGACTAAATTAGATCGTAAAATTAGCAAATTAGCAAATTAGGCGATATAAATTAGAAATTAGATAAATTAGATAAAATTAGTTAAATTAGTAAATTAGATCAAATTAGGACCAAATTAGTGAAATTAGGGAAAATTAGCCATTAAATTAGTTGAAATTAGGTAAATTAGTAAAAATTAGCAAATTAGGGTAGATTGGAAATTAGGTAAATTAGAGACCGCCAAAATTAGGGTAAATTAGGTAGAAATTAGTAAAATTAGGCATGAAATTAGAGATAAATTAGTAAAAATTAGCCGAAATTAGAAAATTAGCAAATTAGAAATTAGGTGCCAAAATTAGATGGATAAATTAGAAATTAGTGTAAATTAGGAAATTAGAAATTAGAAATTAGTCTGAAAATTAGAGAAATTAGTCGATTAAAATTAGAAATTAGTTGTGGTCAAATTAGAAATTAG"
    pattern = "AAATTAGAA"
    genome = "AAAAGGGG"
    symbol = "A"
    # print(bio.SymbolArray(genome,symbol))
    # print(bio.FasterSymbolArray(genome,symbol))
    # Genome = "CATGGGCATCGGCCATACGCC"
    # print(bio.skewArray(Genome))
    # Genome = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
    # Genome = "GATACACTTCCCGAGTAGGTACTG"
    # print(bio.MinimumSkew(Genome))
    string1 = "CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA"
    string2 = "CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG"
    print(bio.hammingDistance(string1,string2))
    # pattern ="ATTCTGGA"
    # genome="CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
    # d=3
    # print(bio.approximatePatternMatching(genome,pattern,d))
    # pattern ="GAGG"
    # genome="TTTAGAGCCTTCAGAGG"
    # d=2
    # print(bio.approximatePatternCount(genome,pattern,d))
    a=list(range(5))
    b=a
    a[2]=12
    print (b)

if __name__ == "__main__":
    main()