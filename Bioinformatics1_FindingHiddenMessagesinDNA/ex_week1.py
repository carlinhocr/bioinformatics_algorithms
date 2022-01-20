class Bioinformatics(object):
    def __init__(self):
        pass

    def patternCount(self, text, pattern):
        count = 0
        for i in range(0,(len(text)-len(pattern))+1):
            subtext = text[i:len(pattern)+i]
            if subtext == pattern:
                count +=1
        return count

    def frequentWords(self,Text, k):
        words = []
        freq = self.frequencyMap(Text, k)
        m = max(freq.values())
        for key in freq:
            if freq[key] == m:
                words.append(key)
        return words

    def frequencyMap(self,Text, k):
        freq = {}
        n = len(Text)
        for i in range(n - k + 1):
            Pattern = Text[i:i + k]
            freq[Pattern] = 0
        for i in range(n - k + 1):
            Pattern = Text[i:i + k]
            freq[Pattern] += 1
        return freq

    def reverse_v2(self,Pattern):
        return Pattern[::-1]

    def reverse(self,Pattern):
        reverso = ""
        for i in range(len(Pattern) - 1, -1, -1):  # voy hasta el menos 1 para    que range incliuya al cero
            reverso += Pattern[i]
        return reverso

    def complement(self,Pattern):
        complemento = ""
        for i in range(0, len(Pattern), 1):
            if Pattern[i] == "A":
                lc = "T"
            elif Pattern[i] == "T":
                lc = "A"
            elif Pattern[i] == "G":
                lc = "C"
            elif Pattern[i] == "C":
                lc = "G"
            complemento += lc
        return complemento

    def reverseComplement(self,Pattern):
        Pattern = self.reverse_v2(Pattern)
        Pattern = self.complement(Pattern)
        return Pattern

    def patternMatching(self,Pattern, Genome):
        positions = []  # output variable
        for i in range(0, (len(Genome) - len(Pattern)) + 1):
            subtext = Genome[i:len(Pattern) + i]
            if subtext == Pattern:
                positions.append(i)
        return positions

#    def findClumps(self,text,k,l,t):

#-------------------------------------------------------------------------------

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

    def motifEnumeration(self, dna, pattern, d):
        positions = []  # output variable
        for i in range(0, (len(dna) - len(pattern)) + 1):
            subtext = dna[i:len(pattern) + i]
            if self.hammingDistance(subtext,pattern) <=d:
                positions.append(i)
        return positions

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
    # text = 'TCTATGGATCTGTCTACATACCTCACTGATATCATACCTCATACCTTCTATGGATCTGTCTATCTATGGATCTGTCTATCTATGGATTGGACGGTGCACTGATATTCTATGGATCTGTCTATGGACGGTGCATACCTTGGACGGTGCACTGATATCATACCTCTGTCTATCTATGGATCTGTCTATGGACGGTGTGGACGGTGTCTATGGATTGGACGGTGTCTATGGATCACTGATATCTGTCTATCTATGGATTCTATGGATCATACCTTGGACGGTGCTGTCTACTGTCTACACTGATATCATACCTCACTGATATTGGACGGTGCATACCTCTGTCTACACTGATATTGGACGGTGTGGACGGTGCACTGATATCTGTCTACACTGATATCACTGATATTCTATGGATCTGTCTACTGTCTATCTATGGATCATACCTCATACCTTCTATGGATCTGTCTACACTGATATCTGTCTACACTGATATCACTGATATTGGACGGTGCATACCTTCTATGGATCATACCTCACTGATATCATACCTTCTATGGATTGGACGGTGCACTGATATTGGACGGTGCATACCTCTGTCTATCTATGGATTCTATGGATCATACCTTCTATGGATCTGTCTACTGTCTATGGACGGTGCACTGATATTCTATGGATCTGTCTATGGACGGTGCACTGATATTCTATGGATCACTGATATTCTATGGATCACTGATATTCTATGGATCACTGATATTGGACGGTGTGGACGGTGCACTGATATCACTGATATCACTGATATCATACCTTCTATGGATCTGTCTATCTATGGATTGGACGGTGCATACCTTGGACGGTGCATACCTTCTATGGAT'
    # k = 11
    # print(bio.frequentWords(text,11))
    # text = 'CCATATGGAACGCCTGTGGGTCTGTTCCCTTACGTCTTTAGGGTGTCGGATATTTTTACCGGTGCACTTAGCTAAATAGTACCTTCGGAACTTTGCCTCCCGCAACCTCGTCCGGGCGGAAAGCACCATATCTGACCGCCAATTAACTGCGGATCCTAGCGTATCCGTAATATAGATTAGTACTTTACGGACACGCTACCACGACCGAGTGTGAGTCAGGCAAGCTCGATTTACATGCGCGACCATTCGATGTCTAGGGTTTGGTGGGTCGCTCACCACCTGACAGGTCGCATTACTATACACCTGAAAGCGGACACGCCATCGGGAGAATCTGTTTAACCGTACACGGCGGGTATACTCCCAATCACACCCCTTTTACTTGATTATGGGATTGCGCTCTAAGCCTCATCTCGTAATGCATGACCGGATGATTATTACGACCGGGGGATAGGCGGGCTATGAAATGCTTTACATCGTGGAGCAAATCGTTCCTATGTACCGGTTCGTCGGCCGTCTGGAGGTTGAGAGGCTTAATGCATACGGCTTCAATAACCGTCCCAGCGTCCCGCTTGTGATCTGAAAACTGCATCGCACTAGTACTACGGTAGGTCCTCGTTTGTATATGACCCACTCTGACAGCGTAGGCTCTAAGATCTTCGTCGACTCACTCTATGGAGTTGACGGCACTAAGGGAGCCGCGGTGTATATGACACGCAGTACTCCGTCCGTGAAACATTCCGGCCTGTCTGTCGAAAGAGCATGATAGCGAAGTATGTCTTCTTCGTTCGTCGTACACTACAATTGTTGTTTTCTATGGCTGGACGGGGACAGTGTCGTCCCCGGTCAAACTTGGGACTGTCCTCCGGCCGTCTAAGCTGCCCCTCCTTAACAATCTCTAGATAAAGGTAACTCCGAGGCATGTTCCGGCAAGTCCTACGCGATCCGCTGAAACTGTTTGGACAGTATGAGTCCCGACGAATTGGAGGATGTATATCACAAGTCGGAAAGTGTAAGATACGAAAGACACGAACTAATCCTATGGCGATATGCTAACTGGGCAGGATTGGCAGAATGGACGTACTGTAAATCCGCACCAATGTGCTCCAAACCGAAACATTGATCACCACTACTAATATCGAAGGTCCTTCTTCTCGCGAAGTGGTATAACCATAGATCTAGAGTTTTGGGACAATCGGTGATGGCTAGCACAATAGTTTCCTATCCACTACAGTTGGTTCTTCTCATTCCCCCCTCGTCGTAAGCAACACTAGGTAGCCGTCATAATGCCTCGGAGTTAGTAGTGCCATAAACATACCAGAAGTGCGACTGATCTGGTGGAGTCGCATGCAAAGCAGTCGTCAAACCGTCTACCGGCACATCAGCCGGCTCGTAATGTGGGTAATTAGAATACGAATAAAATATAGCTTAACGGACTGAGCAGGTCGGCAAGGAAAGGCTGAGGTCCTGAACCGTTTACAACGGATATGCGTAAGTGCTAAACTGACAAAATTAAAACGCGATCGGTCCTGTAGGGGAGCACGGACAAAAATAATCTTGCTGTTGTTTCCTCTATAGCTTGAGAGATCTGCTCATTACCCTGCTAACAGTGGGCATGAAATTAATGGCTCGTATACTTCAAGATGTGTACCTTACTTTCGTACCCGAGGGCCCGCTCTATGCTTCAACCAGCCTTTGCTTCTGTTGAGCCACATAAAGTGATCAAACTGCTCTAGCGGTTCACCCCTAAGACCATTCTAGCTTTAAGATATTAGGCTAGGTCCACATATGTCGTTTCGCAGGCGGGTGACGTCACTGTCTATTGTCCATTTATTGTAGCTACTTGCGGTTGTTTCCAACTTGCACGTTCGGCCGTCTCGTGTTGCTAAGACACGTAGCCGCTTCTTAGAAAGATTATGAGGACGAAGTATCAAGCGTCTAACACCCTATCAGGTATATTCGAGTTTATCAGTTACAGAAACATATTCACTATGTTCAACCTACAAACGCTAGCAAACAGTAGGATAGCACTTGCGAAGACGAAAGGTTTATGAGCCAAATAATGCCCATGCAACGGCTATCCTTCGCAATCTCAGGGGTGCGGCAGATTTTTCTAGTGTTGCGTCACTTCGGTATGTACTGCGTCATGATTAGTGCCTTATTAAGACGTCTCTCGTTGGGCGTACCACGGTATTCTTGTAAAGCCCACCGTACGCCTTCCACAGTCGAGCCCCCTTTAAGCTCCTGAGGCCGCACTGACTTCTTCCACAAGAGCAACTATCACCCGTTCCCATGGATACACCTTCTCTTTTAACGAGCGCAGATTTAATTTGGGGGATCGGCCGACACTTTAGCGTCATAAGTACTTTGATTTTTGTCTATACCAATAGATATATATAATAGCCCCAAGTGGCCGGCCACACCACACTCTCCCCCTTGTTAACCCGCGCCTTATTGTACCTATTCCAGATTCATACTCAGATGTTAAGTGCCTTTCTAGGATACACGAGGCACTTCTAGTTTCACGAACGTGAGGTGGCGTCGATAATCAGCACCCGCAAGCCCTCTAGTCCTATAAGGTGGATGTGCTGGGGGTATCGCCAACCATACCGTTGCGCTTTGAAAGCACGCATGTGCGGAACGGGCCCCGGTTGAGGGTTTAGTATGGTCTGCCAAAGCTTGTACGTACCGTGGGTCAACGAATCCCCCGCTGGATAAAGAGTCAGGAGTGAGGTGACTTATCAGTTTAGGCCTACCCACAGCGACCTCATTTAATCATGTCCGAGCCAGGACACCCCCCCCGCTTAAAGCGTAGGCTGCTAGCTACCACTTAAGGACTTCAATCGATAAACCTATATCACAAGTGTGAGCAACAGTCTAAATAGGACAACACATGTTAGAACTTTCGCAGGGTATCCATGCACTCCCAGTGGTACGGGTTGAGCTCGCGCGCGCAAACGGTATTGTACCACCTAAATATGCCTTAGAACATCCGAGGGAGGTAGCAATCGTGGGCAGAGGGTCCTCGGGATTACTTGTACTAGTTACGGGGGGTGGGTGACTTGCCAGCGTTAGTGCTTATTGATTTAAGCTAGAACCTGGGGGGTGGAAAGTTGACTCAAGGGCACCTAATTCAGCTAGTGTACCTTTAGTTGGCCAGCCATCTCATAGGTGCTGATAATGACGGATTTTTCGATTGGCCGCGGCATACTCTTTCTCCTCGCTCGGTAGGGCGTTCCTTGCGGTGTCGGACTCGACCTCGCACACCTCCGGTATCGTAGAAATACTTCGGCTCGCTAGATCTGGCAAAAACCTCGTATTCCACTGAGATGAGCGCTCTACGCAGAACGCAGCACGCGTAAACTGGGGACGGATAGTCTCGGTGTGGCCACGCGTTTCCAACAGATACAGCGACGAATAGCGTCTGGCTCTCGGCAAATTCTACGGTACTATGAACCACTATGACAACGCTAAAGACGTCAATGTTGTCCAAATCGTCCGCACCGATACATGGTCGAACAGTAGCTAACGAGCGTCTGCCGTGGGTTGTTTGAGATGCTCACTACCCGTCTGAAATACCGAAGGCCCGCAGATACGACCTCCAAGTTCGGTCTTAGACTGTTATTGAAAAAAATGACGATGTGTAGTTGGCCTTCAGTAGGCGGCGTATCCCCTTACGCGCGTTAAGTTCCCGGTAATAATTATATGCCGGCTCACGATGCTGACTGACTAGGCCCTAGGTATTGTTACGATTTCAGTTGTCTATGGCCCGGGAAAACAAATCTCACAGCACGAGGTGTACGCTTTGAGCTCCCCAATCGCAAGAGCCCCCCTGGAGCATTGGCCTTACTGATTAATCATGCGGTCTCAACTATTAAATCAGGGTATGTCGTGAATACAGCCGGCCGGGTGTGCATGTCAATATATGGTAGGGCGATACGCCCACAACCTGGCCAAAGGTAGGTGCCGACTAGAGGTGACCCAATTATCGCAGAAGTTAGTTTTAGCGGGAGGACCAAACCGGGGCCTTACCGCGCGATCTAACAAGTTACGTAAATCAACACCCCTGATACTTTGATAACTGGCGTATTCATGGATGGTGTCCGAATTCCAGTGGAGATTTAAGTACGATTGACGATCTGTGAATGTCTACCCGCTTTAAAACGGACTACTGTGTCAGCATCAAACTCTGTACAAATACAACTTAAGCCCAACTGATATGGCGACAGTATTTAGTCAAATGAATTCACCTGGCTCTCTCCGATTCACGTTTCTTCCTGACGCGAAGAGTCCTCCTGGCGAACCTGAGTTAAAGGCTTAACCAGCTGTATAGTGAGGGTTAAGTATCAACTGTGCTCTAGGCCAAGCTCCACGTGGCGACCCGGTCCCAACCTCCCCGATATGCCTCCCTATAAACTATGTCCCGCATCCAGACAGGGGGCAGTAGACATCCAACGGGACGCTCGAAGGCGCCGGTGGTACAAATTATGATGCGGAATGGTCAGGAGGTCTTAGGGTGTGCCCAGAACACCGCCGGGACAAATCAGAAGATAAGAAAATATACACATACATATGCACTCGTCGTGAAATCACTTTCGTCTAAGGCTCAGCGTCCTTCACATTACTACCGGCGCTGCCGCTCGAATTGCTCGCGCTAGTGTGCGGTACCCTGGTCGTTTGCTAAGGGCGACATACAATAAGTGCCACGAACCTTTTACGCGCGCCTCTCCGCGAAGACTAATCTACTCCACGTGACCAATGGAGTCGAAGTGACCTAAATTGCAAGATTAAATAAGAGCCAACCCTAACTTAACCTACTAATCTATTCAAAAGCAAACGAGTTAATTTCAGCGAACCTATGATCCGGCCAGAGCTAGGAAGGGGGGCCCGTAGTCTTATGTGCATCAAACGCATGGGTTACAATAAAAACCCCGACCCCGGCCCTGTTTTCTACGGTTAAGTATTATACCTAACATGTTACGCTCCGTGCAAAAAACCCCTGAACCTCTAAGTGTTCACTACCTGCAGGTGCGGCGACTTCATGTCGTCTTCGCACGTGGCTCTACCGGACAGAGTGGCCACTTAAGGATATGCGATGACATACGAAGGTTGGTAGCTCCTAAACCGGTTGTGCGCTTGGTGAGCTTATATATTGCCGGTATGCGAGTCAGGATATGAATCATGCAAAACAATGTCGAGGGCCCCACATAGTCTGTTCCCCCTAGCAACACGAGCCTGACAAGCGGTGTCTTATCGGCCTCAATGTCGGGTGCAGACCAGGTGGCCGGAGATAAACTTGGGGGTAACGTAAAAAAGAGAAGAACACGACAGGACGATTTTATTTTGGTTCGCCATGGTAGATTCTAGGTTGGAGGATCAATTTCTTCAATCGTTACTCAATACCGGCATAGGGCCGCTGAAGACGAGTGCACTCCGGAGCCTACCCGGTCACTGGTCTCCCAGCCTGCAAGCTTCTAATAAGACCAGTAGACTTCGGTGGATCGCCCAGGCGAAATCAATTGATGGATCGCCCACACGTCGCAACAGCATCTCACGCTTGAAAGTCACCTTCTGGATGCCGCAACCTGGATCGTCCCAATACCAAATGCCCTGTATCTTGTGGATCCTTCAAGACCGGGCGTTCATCGGTTTAAATAGAGGGTGACTCAAGCTTATCCCGAGACCCATAACCGACACCCACGTAACAGAGTTTTCGTAGAGTGATCTTGCCTAAAACTGTCTGCGAGCCCAGCACAATTAAAGGTGTGGCCCCACGAGAAAGAATGCGCGCAAGTGTTCACATGCTGCACGAGTGGCGTCTCCGCGTCAAGCATTGACCTGCGGCTTGCAACAACAGAGAAACCAGCGAGCACGCGAGAACACCGAGACGCTGACCCCACCAGGGCAGTGCTCTGAAATATTTACTCTTTGCTTAGCGCCGCGTCTTCACCCGTAAAGGGAATTTCTTTAGTGCGCAACGGCCATTGGACCTCGTGGTCGCGCATATTATAGTGTTCGGAAATTAATCGGCACATCTTTATGGATCCCCCTCGCCTGATAAATGCCTCGACTGTGCCTAAGGGCCGATCGGGCTTCCCTAGCTTTTATTGCGTAGGGGACTATGCAACCAGATGATGTACATAATGTGGGCATGCCCCGTTCGGCCCATGACCCTCGTCCCGCCTGTTGTCCGGATAGAGAGGGATGCAGGAGGAGGCAAGCCCGAACACGTCGCCCACTTTGAGATGATCACCAGAACCCCCTATGAGCTCTGTACTAAGACTGGCTTCTAAAGCCTCGCCTTATCGTCCCAAATAGCGTTACTGTAAAGCGCCACGATAACATCAATTCAGACAATGTTGCGCACCGTACTGGGCACTCACAGGAGCGTGATTGGCAGGGAGTCAGATATCAAAATCTACCTTCGAGGGAGCTTCACTAAGAAATACTAACTCTGTATATGGCTAACCAGTTGTACTCCGGGTACATGGCTGTTCACGGATCTGCGGTAAGTCTGCGCGATCGTCTGCTATCCCTGAGACGGCTTGGGCGATAGGTGGTCGTAATCCGCTTGATCAATTTACCTCGCTGTCTGCACTTGTGTCGCTACCACATGACTGCTTCCTAACAGTCGCCACTCTTTGATGCCCGGGGGGAGTCTAAGGGTGCCAGGTTTGGGGCTCACATCATTGGCACAGATAACAGCTGGAGTTAGTAGACGATTCTGTACGTACATTCTAAGGCGCTCTGAAGAATCGCCATATAAATTCCATTGACGTCGTCCGATGCTGCTACCGAATCTGATAGGCTACTACCACCCAAACCAGCCTAGGTATCCCAGTTAGCACATTTAACGTACTGGTAGACATGGCAACGCTGCAGCAGTTACAATACGTCTATTAGGTGCAATCGGTCGCGCCATAGTGCAGGAAGTTTTCGGCCGAGCTCATTTTGGGCAACCTGCGACGACTCGGGAACCGTAAGTATGCAAGACTCTTGTATTGGTTTATCAATGTCACGTACTGCTCAGCTTTGATTTTAACACTCGTATTATATGGCTTTGCACTAGAGTTAAACTATATTGGTCCTAATAGCGGCTTACCCTGCGTAGCTCTCACTCGAACGCAGCTCGTCACTATGGACCGGTCATCGTCATGCCATGGCCCTCCCAAATTTCGCCTTCTAGGGCCCCTCGAGAGTGCTCTCACGTGAGTCTACAGGCGACAAATTAGCCAAGGCGTGCAGCGTCGGGAGTACAATTGCATTCACTTAGTCCCTACCATGTTGAGGAGCCCCGTCATGATCCTTAGATGAGAACTAGCGCGTGCGGCGCGTACTAGGACACGTTTGTTAAATGGGTGACCACATTAACGGGGAATATGAAATTCTACAACATTTAGAAAGAATCAGCAATTGATGTGTTGTAGAATCATTGCCCGCTGTCACCTGATATCTAGATTTTTTGACGGCAAACCACCACGTGGTCACGTGGGCGGTCAACTGTCGTACGGGTTTGTCGTACATTGGTACAATGGACACTGGTTTACGACTACTGCCCAGCCGCGGTTGCATCTAATAGTTGACGCTACTAATCGCAGGGCTGTGCTTTCATGCTTGTTCCCCCGGTATCCTCAAAATGCCGAGCCCTACCAAGCTGCCGATCCAAACCTCCAAACTGCGCTTGATTACCTCCATGACATGTATTGCAGAATCCCCAATTTACCAGAATTTTAGCAACTAAATACCGGTGGAACACAGCCTTAAGGTATTATGACGCGCCCAAGTAGCGGTGTCATTTAGAAACATAGTCTATACTAGGGGGTACCCCATCCCATGCAGTGAAACTATGCCGTTCGATTTGAAACTGGGTCTGATTTCGTTAATACCGCTGCGCGCTAAATGGAATAACACAAGTGTCCTCTACGTCCAATACTCATCCGGGGCAAGAACATCTCGGTGCATCGGAGGTCCCAATGCCATCTCTGGAAGTCAACTCATGGTGGATGACCTCGTAATGCTGGGGGAGGCAACGCGGGCGCTGATGG'
    # print(bio.reverseComplement(text))
    # text = 'CCCCACTCTCTGACCCTGTGCCACTCTAAGCCACTCTTGCCACTCTTCCACTCTCCACTCTCCACTCTCCACTCTGCCTCCACTCTCTATCGCGCCACTCTCGCGCCACTCTAGCCACTCTCTCCCCACTCTGAGGTACCACTCTTCCACTCTACAACCACTCTCCACTCTCCACTCTGCCACTCTTTCCCCACTCTCCACTCTGAACCACTCTCGACCACTCTCCACTCTTCCCACTCTCCACTCTATTCCCACTCTCCACTCTATGCAGCCACTCTGCGTGGCCCACTCTCCACTCTGGTCCACTCTACCACTCTCCACTCTCCCACTCTCTACCACTCTTTTCCACTCTGCCACTCTTCCACTCTGGTTAACTCCTATTCGCCACTCTGACCACTCTAGCGACCACTCCACTCTTGAATCACCACTCTACCACTCTCCACTCTACCACTCTACCCACTCTCCACTCTCTACCACTCTCCACTCTCCGGCCACTCTTCCACTCTGGCCACTCTCCCACTCTAAAGCCCACTCTCCACTCTCCACTCTGCACCACTCTTCCCACTCTATAAGGCCACTCTCCACTCTCCACTCTCTCCCACTCTTACCACTCTGCCACTCTGTTGCGAAAAACAACCACTCTGGTAACTCCCCACTCTTGCCCACTCTCCACTCTGATTCCACTCTAGGTTCGGACCACTCTCTGATGCCACTCTTCCACTCTCCACTCTGACACCACTCTCCACTCTCCACTCTTGCCACTCTTACACCACTCTCCCACTCTGTCCCCCACTCTCCCACTCTCCACTCTCCACTCTGCCACTCTGTACCACTCTAGCACAGCCACTCTTCCTGACACCCACTCTCCCACTCTACAGCCACTCTCCCACTCTCCCACTCTGGCCACTCTAACCCACTCTCAGGCCGCCACTCTGAACGCCAGCCCTCCACTCTCCACTCTTCCACTCTTCCACTCTATTTCCACTCTATCTCCACTCTCCACTCTTTCCCCACTCTCCACTCTCCACTCTGTTTCCCACTCTACCCCCCACTCTCGCCACTCTTCATAGCTCTACCACTCTATCCCCACTCTTTTGACCAGCTAAGCCACTCTTCCACTCTTCCACTCTCCACTCTGTGATCACCACTCTATCAGAGGTCCCACTCTAACCCACTCTTTAGCACCACTCTCCCACTCTCCACTCTTCACCACTCTCCACTCTAACACCACTCTATACAAAATCCCACTCTCCCACTCTCCACTCTCCACTCTCGTCCACTCTGCGCCCACTCTCCACTCTCCACTCTCCCACTCTGTCCCACTCTGTGATGTGCCACTCTCTCCAACCACTCTCCACTCTCCCGCCACTCTCCACTCTAAAACGCTCCACTCTTCTGGCCACTCTTAATCCACTCTACCACTCTTCCACTCTGGCCACTCTTGCTCCACTCTACATCGGCCACTCTTCATGACGCCACTCTAGCCACTCTCCCTTGCCACTCTACCACTCTCCACTCTCCACTCTCAGCCACTAGGCCACTCTCCACTCTTCCACTCTACCACTCTGTTTCCACTCTCCACTCTCCACTCTCCACTCTGCCACTCTCCACTCTCCCACTCTCTCCACTCTCCACCACTCTCCACTCTTGCCACTCTCCACTCTGCCACTCTGTATGCCACTCTCCACTCTACGCCACTCTCCACTCTTCCACTCTCCACTCTCCACTCTGCCACTCTCCACTCTCCACTCTGTAATCCACTCTACCACTCTCTGTTCAAATACCACTCTGTTGTCCACTCTCCCACTCTACCACTCTACCACTCTGACCACTCTCCACTCTCCACTCTGCCACTCTCCACTCTACCACTCTTCCACTCTCCCACTCTGGCACCATCCACTCTCTACCCCACTCTCAGACCCACTCTCCACTCTCCACTCTGCCCACTCTTCGAGCCCACTCTCACCACTCTGCCACTCTTTCTACCCACTCTTAGGGCCACTCTCCCACTCTGAGTTCCCACTCTCCACTCTAGGCCTAGCCACTCTGGACCACTCTCCACTCTCCACTCTCCACTCTCGGCCACTCTGTGATTCCACTCTCCACTCTGTCCCCACTCTTACGCCACTCTCCACTCTCGCCCCACTCTCCACTCTAATGGCCACTCTCCACTCTCGCCACTCTACGGCTGCCACTCTTTCCACTCTCTCCACTCTGCCACTCTTTCCACTCTCCACTCTCGCCCACTCTCTTCGCGCCACTCTTCCACTCTCGCCACTCTCCACTCTCCACTCTTCTACGCCACTCTAATCTATGCACTCCCACTCTCCACTCTAGTCCCACTCTCCACTCTGTGAATGAATATGCCACTCTATCCACTCTCACCACTCTCCACTCTCCACTCTATCCCACTCTCCACTCTCCACTCTCCTGGTGCCACTCTGCCCACTCTCGCCCACTCTCCCGTGACCACTCTCCACTCTCCACTCTCCACTCTCCTGCGTCCAAACCACTCTACTCAGGCCACTCTCAGCACCACTCTGTCAACGCCCACTCTGTGACCCCACTCTAGCTCGGCGTCCCACTCTCCACTCTTGCCCACTCTACCGACCACTCTCTCCACTCTCAAGCGCCACTCTTTCTCCCCACTCTACCACTCTTCCACTCTTTAGCATGGGGTACCACTCTGTCCACTCTCCACTCTGGCCACTCTTGCCACTCTCCCCACTCTCCGCCACTCTTCCCACTCTGGCCACTCTGGAACCACTCTCCACTCTTGATTCCACTCTCCACTCTGATTCGTCCACTCTCCACTCTATCCACTCTCCCACTCTCCTCCACTCTATCCACTCTGTACCCACTCTTCGCCACTCTCCACTCTGCCACTCTTCCACTCTACCACTCTCCACTCTACTTCCACTCTCCACTCTCCCCACTCTCCACTCTCAGAGCCCACTCTCCACTCTACCCACTCTACCACTCTCCACTCTGTCCACTCTCCACTCTCCACTCTAAACCACTCTCCACTCTTCGTTCCACTCTGTCCACTCTCCACTCTCCACAACCACTCTCTCTCCCACTCTCCACTCTGTCCCACTCTTGTCTTCCACTCTCCCCTGGTGCCCACTCTCCACTCTCCACTCTTGCAAACCCACTCTCCCGCCACTCTCGACTCCACTCTTACCACTCTGTCCACTCTAGTTCACCACTCTCCACTCTTCCACTCTGGGACCACTCTCCCACTCTCCACTCTGCCACTCTCCACTCTCCACTCTCCACTCTGCTTAGGCGCCACTCTCACACCACTCTTTGACCACTCTGCCACTCTCATTCCCGATGTACCACTCTGGGGTCACACGGCCACTCTCCACCACTCTCTCAATCCACTCTATCCACTCTACCCACTCTTTGAGGATTATTGCCACTCTGCCACTCTCTCCACTCTGGGCCACTCTGCCACTCTGCCACTCTCATGCCACTCTTTCCACTCTCCACTCTGGTCCACTCTACCACTCTAGGGCACAGCCACTCTCCCCACTCTTCCACTCTCCCACTCTGCTGCCCACTCTCAAACCACTCTCCACTCTCGCCACTCTGGCCCACTCTAGACCACTCTCTTACCACTCTATACCCACTCTGCCACTCTGCCACTCTCCGATACCACTCTCCACTCTATCCACTCTTCCACTCTCCACTCTAGATTTTCCACTCTCCCACTCTGACCACTCTTTTCGAAGCACCACTCTCCACTCTCCACTCTGCCACTCTGACCCACTCTACCACTCTCCACTCTGCCACTCTCCACTCTTCCACTCTCCACTCTGCTGCCACTCTCCACTCTCCACTCTCACCACTCTTCCACTCTCCCACTCTTTCCACTCTCGCCCACTCTGCACGCCAGACGTGACAACCACTCTTGGCCACTCTAATCCCACTCTACCACTCTAGCGATTCCCACTCTATCACCACTCTCCACTCTGATTTCCACTCTTGCATCCCACTCTCCACTCTAGAACTCCCTGCTCATACCACTCTCCACTCTTCCCACTCTCCACTCTGCTCGTGCCCACTCTATTCAAGCGCACAGAAGCCCACTCTGTCCACTCTACCACTCTACCACTCTTTTCGGTCCACTCTACCACTCTGACCACTCTTCCACTCTTGTCACCACTCTCCACTCTCGCCACTCTCCACTCTAGAGGTTCCACTCTCCACTCTAAGCCACTCTATCCACTCTTCCACTCTATCCACTCTCCACTCTCCACTCTCCCACTCTTCCACTCTGCCACTCTGAATGAGGTGCCACTCTGGATCCACTCTCCACTCTCCCCACTCTAGGAGGGCCCACTCTCTCGCTAGTCCACTCTCCCACTCTTCCACTCTTGAGTAGATCATCAACCACTCTTGAAGCCCCACTCTTACGAACCACCACTCTAGACCACTCTCGCCACTCTCGTAACCACTCTCGCCTCCACTCTTCCACTCTGGCCCACTCTCCCACTCTTTCTCCACTCTACAGAATCCACTCTTACCACTCTCCCACTCTGCATTCGTTAATTTCCACTCTATACCACTCTTCCACTCTACCACTCTTTCTCTGTAACTGATCCACTCTTTAGAGGCCACTCTGGTACCACTCTGTCAGTTCCACTCTACCACTCTCCCCACTCTCGCCACTCTGACCACTCTCGCCCACTCTGGACCCACTCTCCACTCTCCACTCTCCACTCTAAGCTCCACCACTCTGCCACTCTTATGACCACTCTTCCACTCTGGCCACTCTCCACTCTCCACTCTTCCACTCTCCACTCTCCATCCGCACCACTCTTAGCCACTCTCCACTCTACCCACTCTTGTGATCCACTCTAGCCACTCTTGCCACTCTTTCCACTCTATTCCACTCTCCACTCTGCCACTCTCTTTGGCATATCACGACCACTCTCCCACTCTAACCACTCTACCACTCTTCATACCTTCCACTCTAACCACTCTGCTCCACTCTGGTATAACCACTCTTACCACTCTCATCCACTCTTTACCACTCTTCCCACCACTCTACCACTCTCCACTCTGCCACTCTTTGACCACTCTCCACTCTAATACCACTCTTCCACTCTCATTCATCCCCACTCTTCCACTCTCCACTCTTACCACTCTGGGGTGGGCGCCCACTCTCCCACTCTCCACTCTCCACTCTCCACTCTCCACTCTCCACTCTTCTGCCACTCTCCACTCTCCACTCTCAAACGGGCCCCACTCTCATCCCCACTCTACCACTCTCACGACATAGCAACCGTCCACTCTCCACTCTCCACTCTTCGTCCACTCTCCCACTCTACGGCACACCACTCTATTCGCGCCACTCTACACCACTCTTTGGGCGGTCCACTCTCTAAGGCCACTCTCCCCACTCTATGAGCCACTCTTCATATCCCCCACTCTCCCACTCTCCACTCTCCACTCTGGCCACTCTGTCCACCACTCTCACCACTCTCCACTCTCTAACTCCACTCTCCACTCTCGTACCACTCTAGGATCCCCCACTCTCCCACTCTCCACTCTCCTTTCCACTCTATCCACTCTCCGGGCGTACCACTCTTTTAGCCACTCTGTAGTACTCCCACTCTTCCACTCTGCCACTCTTGCCACTCTTAGGTGCCCCCACTCTAAAAGCTAACCACTCTCCACTCTCCACTCTGTGCACCACCCACTCTGCCACTCTATACCACTCTTCCCACTCTGGTCCACTCTGTCCACTCTTGGCCACTCTAACCACTCTAAACCACTCTGGTTCCCACTCTGATAGACCCCCACTCTCCACTCTGTCGCCACTCTCGCCACTCTAGCCACTCTCCACTCTCCACTCTTGCCACTCTGCCACTCTTCCACTCTACCCACTCTGATCCCACTCTACCACTCTATCCACTCTTCCACTCTCTCCACTCTTACCCCACTCTTTTAACCACTCTCACCACTCTCCCACTCTGCCACTCTCACCTTGACCACCACTCTTATGATGAATAACCCACTCTCCACTCTGGCCCCACTCTATCCACTCTCACCACTCTGCCCCACTCTCCACTCTGGAGATATCCACTCTTCCACTCTAAATCCACTCTTATCGTTGGTAACCACTCTATACCACTCTCCACTCTCATGCCACCACTCTCCACTCTCCACTCTAAGCGAATTGCCACTCTATCCACTCTCCACTCTCCACTCTCTGTGCCACTCTGGCTCGACCACTCTGAACCACTCTTCCACTCTTACCACTCTCCCACTCTCCTTACGGCCACTCTCCACTCTTGTTCCACTCTCCACTCTCCACTCTAACCACTCTCCACTCTACCCACTCTATGCCACTCTTCCACTCTCCACTCTCAGTGCCACTCTCCACTCTATCCACTCTCCACTCTCCACTCTTTCCCCACTCTCCACTCTCCACTCTCCGAACCACTCTGCCACTCTGGGGGGCCACTCTCCACTCTTTCCGGAACCACTCTGGCTACCACTCTGCGCCACTCTCCACTCTCCACTCTCGCCACTCTGCCACTCTCCACTCTTCCACTCTGCCACTCTCCACTCTTCCACTCTGTGGGTTAACCACTCTCCACTCTATAGACCCACTCTACCACTCTACACCACTCTAGAGGCCACTCTTACCACTCTCGCCACTCTTGATAAACCACTCTCCCACTCTGCCCACTCTTCATCGGTCCACTCTCCACTCTTCAAGCTGAGCTCTTAGGATTCCACTCTTCCACTCTTCTCGCCACTCTGGGCCACTCTAGGTCCACTCTTCCACTCTCAGCCCACTCTCGATCCCACTCTTCCCACTCTCTGGCCAAAAGGTTTACTGACCACTCTGCCACTCTTTCCCACTCTAGTGTCCACTCTGGACCACTCTGCCACTCTACTCCACTCTCCACTCTCCACTCTCCACTCTCCACTCTGCTCCCACTCTTACCGACTCCACTCTGCCCACTCTTTCCTCCACTCTGTGTGCAGGCCACTCTTCCACTCTGGCCCACTCTTGCCACTCTGCCCACTCTATTCCACTCTGCCACTCTGTCGGACCAGGCCCACTCTCCACTCTGGTCCACTCTTGAGGCAGCTATCCACTCTGCCACTCTCCGAGGCATCCACTCTTCCACTCTGTTCCACTCTCCACTCTCCACTCTGCCACTCTGAACAACCGGCCACTCTCCACTCTCCACTCTCCACTCTATGTTGTGGTGCTACCACTCTCACCACTCTTCCACTCTCCACTCTTAATGCCCACTCTCTGACCACTCTGACCACTCTAGGTTACCACTCTGCCACTCTAACCACTCTCCACTCTTGCCACTCTCCACTCTCCCACTCTCCACTCTCCACTCTACCACTCTCCACTCTTCCACTCTACCACTCTCCACTCTCCACTCTGTCCACTCTTTACGTCCACTCTCTCCACTCTCCACTCTTGTCCACTCTTCGCCCACTCTTTATTCCCACTCTAGCCACTCTCCGCCCCACTCTACCACTCTCCACTCTCCACTCTCCACTCTCCACTCTCCACTCTGCCACTCTCCACTCTAGTACGCCACTCTGCCACTCTACCAGCGCCACTCTGACCCACTCTCCACTCTCCACTCTCGTTACCACTCTGCCACTCTAGCCACTCTACTACCACTCTCCACTCTACCACTCTCACCACTCTCCACTCTCCACTCTTCCACTCTTTGAGGCCCACTCTTTGCCACTCTGCCACTCTATCCACTCTTCCACTCTCCACTCTCCCCACTCTTACCACTCTCCACTCTTGCCCACTCTCTTGCCACTCTGTCCACTCTTGAGATTACCCACTCTTACTCCCCACTCTACCACTCTCCACTCTGGTCCAGACTTGCCACTCTCCACTCTTCTCCCCCACTCTGGCCCACTCTCCACTCTCTCCACTCTCCACTCTATTGTCCACTCTCCACTCTCGCCACTCTAGACCACTCTCCCCACTCTCCACTCTCCACTCTTCCACTCTTCATTACCACTCTCCCACTCTGTCCCCCACTCTCCGCCACTCTCCCACTCTCGTAGCCCACTCTCCACTCTATTCACAGCCCACTCTGCTACCCATCCACTCTATGCCCACTCTCCCACTCTAACCCACTCTGAATAGACCACGCCACTCTCCACTCTCCACTCTACCACTCTCCCACTCTCCACTCTTATGACCACTCTCCACTCTAACCACTCTGTCGTGCCACTCTCCCTATATCCACTCTATAAGTTCGACCACTCTTCCACTCTCCACTCTTCCACTCTATTGTGGTGCCACTCTCGCCACTCTGACCACTCTGAGACCACTCTCTCCCACTCTCCTTCCACTCTCCGACCCACTCTTGGCCAACCACTCTCCGGTCCACTCTCCACTCTGGCTAGCCCACTCTAGGCCACTCTGAACCACTCTCTACTCCGTATATCCACTCTAGTCCACTCTGATCCACTCTGCCACTCTGAACCACTCTGGACCACTCTCCACTCTCGCCACTCTCCACTCTGGCCACTCTCCACTCTCTACCACTCTCCACTCTCCACTCTGAATTGATCTGTGGGCCACTCTCCACTCTCCGCCCACTCTGCCACTCTCCACTCTACTGGCCCACTCTACGCCTTCCACTCTGAATGCCACTCTCCCACTCTCCACTCTCGCCATATACCAAGCCACTCTCTTCCACTCTCCACTCTCCCACTCTTCGTTGTCCCACTCTGCATGTACCACTCTCCACTCTACCGTGGGTTCTGGACTACCCACTCTCGCCACTCTCCACTCTCCGTCGCCACTCTTGGATTCTGCTCCACTCTCCACTCTAGTCCCACTCTTCGACCACTCTCGCCACTCTCCACTCTCGAATCCACTCTCCACTCTAGGAGTCCACTCTCCCACTCTATCCACTCTGCCACTCTCCCACTCTTCATACCCACTCTCCCACTCTCCCACTCTCCACTCTAGACCACTCTCCACTCTCGCCACTCTCCACTCTCCCACTCTGTCCACTCTCCACTCTCCACTCTCCCACTCTTCCACTCTCCACTCTAGCCACTCTCGCCACTCTAAGCCCACTCTAGTATCCCACTCTAATTTGGCCACTCTCACCACTCTGCCACTCTCCACTCTGGTGAATGCCCACTCTTCCACTCTCCACTCTCCCCCACGTCCACTCTACCCACTCTCCACTCTCACCACTCTAGCTCCACTCTTCACCACTCTCCACTCTCCACTCTCCACTCTCCACTCTACCACTCTCTCCCACTCTCGACCCCACTCTCCACTCTCCACTCTCCCCCACTCTCACCACTCTGAACCACTCTCCACTCTAGTCCCACTCTCCACTCTCCCCACTCTTTCCACTCTCCTTTGATCCACTCTCTAACCACTCTGGCCCACTCTCCACTCTAAACCCACTCTTCCACTCTTTCCACTCTCTTTTTCCACTCTCCACTCTACGACCACTCTAGACCACTCTCCACTCTCCACTCTTGCCACTCTCTCTCCACTCTCCTTTACGGGTCGCCACTCTGCCACTCTCCACTCTCCACTCTTCCCACTCTATTTAACGTACCACTCT'
    # pattern= 'CCACTCTCC'
    # print(bio.patternMatching(pattern,text))
    with open('Vibrio_cholerae.txt','r') as fileVibrio:
        text = fileVibrio.read()
    pattern = 'CTTGATCAT'
    print(bio.patternMatching(pattern,text))


if __name__ == "__main__":
    main()