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

    def frequencyTable(self, text="hola", k=0):
        freqMap = {}
        n = len(text)
        for i in range(0, n - k):
            pattern = text[i:k+i]
            if not pattern in freqMap:
                freqMap[pattern]=1
            else:
                freqMap[pattern] += 1
        return freqMap

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

    def frequencyTable(self, text="hola", k=0):
        freqMap = {}
        n = len(text)
        for i in range(0, n - k):
            pattern = text[i:k+i]
            if not pattern in freqMap:
                freqMap[pattern]=1
            else:
                freqMap[pattern] += 1
        return freqMap

    def MaxMap(self,freqMap):
        maximum = 0
        for element in freqMap:
            if abs(freqMap[element]) > maximum:
                maximum = abs(freqMap[element])
        return maximum

    def BetterFrequentWords(self,text,k):
        frequentPatterns = []
        freqMap = self.frequencyTable(text,k)
        max = self.MaxMap(freqMap)
        for pattern in freqMap:
            # the patern is the key in the dictionary, the value it is how often it is found in text
            if freqMap[pattern] == max:
                frequentPatterns.append(pattern)
        return frequentPatterns



    def findClumps(self,text,k,l,t):
        patterns = []
        n = len(text)
        for i in range (0,n-l+1):
            window = text[i:i+l]
            freqMap = self.frequencyMap(window,k)
            for element in freqMap:
                if freqMap[element] >= t:
                    patterns.append(element)
            #remove duplicates by converting to a dictionary and then the keys to a list
            patterns = list(dict.fromkeys(patterns))
        return patterns




#-------------------------------------------------------------------------------


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
    # with open('Vibrio_cholerae.txt','r') as fileVibrio:
    #     text = fileVibrio.read()
    # pattern = 'CTTGATCAT'
    # print(bio.patternMatching(pattern,text))
    text = 'CCGGTATATACAGATGACGGCTACCGACTTCACAGTTGTTGCGCCAGAACAGTGGGACGTTATTACAATCGGGTGCGGCATGCCCCGAAGGACAATGTCTCTTTTGCGGTGACTCGTGTGGCTCGGCATTATCGTAGAAAGGGCTTTTCTTTTGCGGTGTTTGCGGTGGCGGTGGTGGTGATACGTAAATCCGGAATCCGGTCAGCCTTTTGCGGTGGGACGCTTTTGCGGTGGACCGTTTTAACCGTTGTGTAAATTTTTGACCCTTTTGCGGTGTGCGGTGATCAATAGCCTGTCTTTTGCGGTGCTTTTGCGGTGACATCAACGGACCAGCTATAGGACATGACCTTTTGCGGTGCTATTTCCACCGGACTAAGTTCATAGTCATGTCTTTTGCGGTGTTTGCGGTGACCGTTTTGATGCGCCCTTTTGCGGTGTTTGCGGTGCCGCTTTTGCGGTGTACCTTTTGCGGTGTTTGCGGTGCCCCACCCTGTTCCCTTTTGCGGTGTTATTCTTTTGCGGTGGATGCCACGCCGTTCCAAGGCCCTACCGCGTAATCTGTCTTTTGCGGTGGAATAAGGATCCAGTATGTCTGGTGTATGCACCTCGGCCTTTTGCGGTGTTTTGCGGTGGTGAAGACCTTTTGCTTTTGCGGTGATACCTTTTGCGGTGACCGGTAGCCAGGATGCGAAGCTTCCCAAGTCTATGGCCCAAGTCTAGCGTCCGATTTCTCTTTACCCACATGGCATTATAAGTCAATTATTGACATGGTCTAAAGGATTTGCGGGCTCCGGGGCAATAAATCGTCAATCTCGGTCACAGTCTGCAAAATATTAATAACTATTATTGCCTTTGGTGGATAACCGATCACGACTATGTTGCAGAAGGTTTACGCTCTTTTTAAGCGGCCGACAACATAGCTGATCCTAGGTAATGGTAGAAATCGGATCGGACGCCATTTTTTTGTGCCTTGAGGTGGAACAATAGGGCAAAGTGAGGGGGCCGAGTGGCCATCAGCGTCACCATTGTTTCTCGCGTGTGTATGCAGGTGACTATGGTGTTCGACAAACAGAGTACTTGAGGGATGAAGGAATCTCACAACGCCGAATTATGTCATGACGCTTTGTCACGCACACTGTGGTCCAACCAGCTTCCTTTAAGGTGGTGAACACCTAGAGGACCATCAGATTTCTTACAGGACCTCAAAGTCTGTTTCGCGTGGACTGCGGTTACTGTTCTGAGAAGTAGGGGGAGTGGAACATAAGTAGAAGTAGGGGGATCCCAAAAGTAGGGGGACCGTAGGCCGGGCATTGAAGTAGGGGGATAAGTAGGGGGAAGTAGGGGGAATGGTCTAAGTAGGGGGAATGCTAAGTAGGGGGACTAGAGTGTAAGTAGGGGGACGATCTACTGCACCCTGTCCGGGGAAGTAGGGGGAGGGATTAAACACATTTAGCCTTAAGCAGGGATGGGGCTAAAGTAGGGGGAACGATAACTATGTGATATCTGTAAGTAGGGGGAAGCGTTCCATGGAGAAGTAAGTAGGGGGAGGGGGAAGGGGGATAGGGGGAGGTATCGAAGGAACACGTCCCATGCAAATCATTACTGAACTTACACCCGAAGACCAACCTTCAGAAGTAGGGGGAATGAAGTAGGGGGAGCGGGTTTAAGGGCACCAAGAAGTAGGGGGAGACCACTTTGCTCATGTGCTCATGCGGGATGCTCATGCGGCTCATGCTCATGCTGAAGTAGGGGGAGCGGAAGTAGGGGGAGAGAAGTAGGGGGAGAGCGGCAAGTAGGGGGAAATAAGTAGGGGGAAGGGGGACTTGAAAATTGCTCATGCGGTGCTCATGCTCATGATAATGAACGGATAAATAATGAACGGCATGCGGCATGCGGAAATGCCATCTTTTGACGGAATATAATGAACGGGATAATGAACGGCGTCTTGATGCTCATGCGGCATGCGGATCTTTGCTCATGCGGTGGTCATAATAATGAATAATGAACGGGAACGGGGCTCTGCTCATGCGGATGCTGCATAATGAACGGTGATAATATAATGAACGGAAATAAATAATGAACGGGTCTCTTTCCATGCTCATAATGAACGGGAACGGCATGCGGGGAATAATGAACGGCATGCGGATGCAAGTGCTCATATAATGAACGGGTGCTCATGCGGTTAGATCTTGAAAATGCTACATGTGCTCATGCGGTTCTTGAAAATGCTCATGCGGATGAATAATGAACGGGCATTGCTCAATAATGAACATAATGAACGGTGAATTTCAACGTATTAAAGTGATAATGAACGGTGATAATGAACGGAACGGATAATGAACGGTTCATAATGAACGGGGGAAAATGTTAAATAATGAAATAATGAACGGAGACATGCTATATTACTCTGTTGATAATGAAATAATGAACGGATAATGAACGGGGACCCCAGCATATGCTATTGATCTGGTTAACTGCCCTTGTTGCAGAGCCCTCATGCTCTCACCAGTACGGGAATTTGCATACCGGCGAAGTGGCTTCATTTGTAAATGAAGTGATCAGCGTACTTGGAGGCGTCCCCTAGGGACCGATTTCGTGACCCTGTACACGATCGGAGTTCTGGTTCCCGAAGCGCATCGCTTTTTTTCAAGCCAGGTTCGAGGTGGATTTAAGCTTACAAGGGCTGGTAGTCGATGGATCAAAGTGACTCATTACACACTCAATCCACCTGACATAGGCATACTCGCATTGCGAACGCCAGCTGTTGAATGCTGTCTACTCAGACTGAGGACAGACGGTCTAAGGACTAAGGACAACAACAAGGAGCGGCGCCTAAGTTACGCGCTAAGGCTAAGGACAACGGACAACTGTATGCCGCTTCTAAGGACAACAGGACAACACGAACTAAGGACAACACAGCCTTGTTGGAGAGGGGGCGGACACTGTCCAGTGTTATGAGCTAAGGACAACGTATCTAAGGACAACGACAACATCCGTATCGGTTAGGCGGTGTGCCGCTTTTCGCCTTTTGGCATCTAAGGACAACCTAAGGACAACGACAACCCTACTAAGGACAACACAACCGTTTCCACGCAGGTTGTCTCCGGTAGAACTAAGGACAACTCTACTAAGGACAACGCTAAGGACAACCAAGCCGGATAGCACACGCAGCTAAGGACAACGCTAAGGACAACACAGAATCGGCTAAGGACAACAGGTGAAATGGCGCTAAGGACAACTCTCTACTAAGGACAACCAGCTACTAAGGACAACCATAGGTAGTAGGCAGCCTCCGGGACCCACTTTGTTTCTAAGGACAACCACCTAAGGACAACTAGGATAGTTCTAAGGACAACATCCGTGTATCTTCGTGCTAAGGACAACCGATTATAGCACGTCACCAAAGGCTAGAGTCACGCCTTTCGGAGATATCAGAGAGCTTACCCACCACTCCATTCACCTTACTTAGAAACGACATGTTTAGCTAAATATCGCTGATGGGGCTCCTTATCATGGGGTAACGGATTACCTGTGATAGACAAATAAGGCAGAGACACCCCCCAATAGTGAACATTCTTTAATATTCTGGGAAAGTGGATAAACCAGGCACCATTTGCTAGCATTAAGGTCGCCGCTGGGATGCTCCTCGCGCGCGTATCCTCGTTGGCTAAGCTCGACAAAGATGCAGGAGGAGCCGTCATAACAAAGCCGAAACGGTCCTAGATTTACTGACATCTGCTTGTTCGGGGAATGCGCAAAGGTCTGTAGCTTAGCAAACATATAAAGATCGCTAAGATCGCTCCACTGGTGGATGTGCATTAGGAATGTATAAAGATCGCTGTAAAGATCGCTGATAAAGATCGCTAGCAACAGCAAGCCTCGTAAAGATCGCTCGTAAAGATCGCTTAGCATCGTCCCTGTCGCAAATTATAAAGATCGCTCCTGGCCCGTATCAACCTCCGCACATCCTTGCTTCGCTAAAGATCGCTAAAAAGGATAAAGATCGCTTCGTAAAGATCGCTGGCCATCCTTATTACGACCCGCGCCCGTATGTAAAGATCGCTCACATAATAGTGGGTGCCGTGCAGCCACCCTTTACCACACAGGTACTAAAGATCGCTCTAAAGATCGCTCAATGCGACATAAATAAAGATCGCTAAAGATCGCTCTCCTAGGCTAAAGATTAAAGATCGCTGAATACGCTCGTCATAATGATCCCTAAATAAAGATCGCTAAGATCGCTTTTTTTTAAAGTAAAGATCGCTCTTTGCTTGCTATTGTCTTTCGTTTGCACACTCTTAAAGATCGCTAAGATCGCTAAGATCGCTGCTTAAAGATCGCTCCGTACGTAGGTCTGGTGACGATTGCATAAGAGAACGTCTCAGGTAAGGCCGTAAGTACAGGGGGCACTATTGGTCATCTCGTGACTCGCTTGGAGAGGCTGCTACATCCGCGCATTAGTAGATCTCGGTGGTCATTGGAGCGATACAAAAATTACAGCCGCTCTGAGCTACAAGAGATCCGATTTGTCTCGCTAAAGTCCAGCCCCTTCGAATCGTATCTGTAAGTTACAGGTCATGATGTTTAGCTTGATTCGGAGGAAGTAGATCGGCGTACCAGTTCCAGGGTCCTTTTCCGGAAAGCCCTCGTGTATCCCTATCCGAAATTACATGTGTCTAAATCGCTGGTATTCCTAGACCCGTTCCCACTCTCTCTCAAGACGTCTGCAGTCCTGAACTCGCAGCTATAAATTGCCATACACTTAGGTGAAGCCGCACCCTAACCAACGTTCTGTCTGAAGGTGACTCTTGCTCGTGCCAAAGACATCCGACGCTCCGAATGGGGTTACCCCTGCCTCATTCTCACCCTACACTCGAGGTTCATCCGGTGTCTCTATGCTCCTATTGGCAATCACGGCTTTCTGGCGGCCTTCGGGGATGTCAGGGAAATGGAGACAAACAGAGCAGGCGAGGGATGCCAAGACGTTTGGCTAAGAATGGCAACAATTAATTTGGAAGAACTGATACTTACTAAACGGCGATCGGGACATACTCGTTTGTCTCCCGAGGTTACTCGGATGGTCAACCTATTTGATTTAGTGTCGTGCTGTCGCTAATCCGCATCACCCCTTCTGCCGTGCCAGACCTTTTCCCCCGGCGTGATATATTCTCGTACAAGTTAAATACATAACTGAAGGAGTAACTCGGAAGGGGATAACCCCAGTCTCCTCGCGTACGACAATAGGCGCCGGTGCCAGCAGTGTAACAAAGGGCTGGCATATCAGGTAGGTCTATCCTTAAAGGGACCATGGAATATGCATAAATCCCAGGTTGTCAGGTGGACTCGCCTGAGTAATAAGTCGCCAAGATTCAATAGAAACAGCTGCATATTATTGCAGGAGCAGTCGACGGAGGAGCAGTCGCAGTCGTCCTTGACTAGGAGCAGTCGGGACTTAGGAGCAGTCGAGGAGCAGTCGCAGAGGAGCAGTCGACGGCTCCAGGAGCGACGTGGTAAGGAGCAGTCGAAGCGTTCACAGGGTGAAGAAGAAAGGAGCAGTCGATTAAACATTAGGAGCAGTCGCATTACCGCTTTGCGGAGGAGCAGTCGACTCACTTCGTTCTCGCTGACAACTGGAGGAGGAGCAGTCGGTAAGGAGCAGAGGAGCAGTCGCGGTAAGGAGCAGTCGGGAGCAGTAGGAGCAGTCGCATGCCCTCTACTCAAGGAGCAGTCGTTGGGCTCCAATTACAAGGAGCAGTCGAGACGTACTCCAGGAGCAAGGAAGGAGCAGTCGAGGAGCAGTCGGAGAATGAAAGGAGCAGTCGAGGAGCAGTCGATTAGGAGCAGTCGTCTCATTTTGACCAGGAGCAGTCGAGGAGCAGTCGGAATATTGTTCCACTACTTACAGAGACCAGGAGCAGTCGCAGAAACGATTACACCGGTACGGTAGAAGGAGCAGTCGCTACGTTTAAGGAGCAGTCGACCGGTTACACCGGTAACCGGTTACACCGGTACCTTACACCGGTACTTACACCGGTAGGCGCCAATTACACCGGTACCGGTAGACCACCAGTTACACCGGTTACACCGGTAGGTATGACCACCAGAATTACACCGGTAACCTTACACCGGTAGGTAGAACCACCAGAACATTACACCGGTAAGGGCAACGTTACACCGGTACGACTTTTTTACACCGGTACCTTCGTCCGATTTTCATCGTTGGCGAAATATGTTACACCGGTAGGCCTATTATAACGAAGACGTACGTACTCCTTCCATGAAGATCACATCGTACGTAACGTGGACTTACACCGGTACCTTACACCGGTACCGGTATATGATTACACCGGTATAAAGGAGACCACCAGAACCTTACACCGGTACTAATTACACCGGTAAATATTTACACCGGTAAATATGGCCTAATATGGCCTATGACCTAAAATATTTACACCGGTAGAAAAATATGGCCTAGGAGCAATTACACCGGTACACCGGTAAAAGGGGGCATACGTATCTCCTCACAATATGGCCTACCCTTTCCTGCTCGAATATGGCCTAATGGCCTAATAATATGGCCTAAATATGGCCTAATATGGCCTAGGCCAAATATGGCCTACTAGGAGCCTGGTGATTGGTCCTCGGAGGGAAAAATATGGCCTATGGCCTATGCAGAACACAGGTAGAACCCAGAAATATAATATGGCCTAGGCTTGTATTGATGACTGTAAGATCTTAGGCTAGATCCGCATAGCAAGCTTATAATATGGCCTAGCAACGACCTAAGAATATGGCCTAGATAATATGGCCTACTAGACATGAAATTCCTATCGGCACTAGATAAACTACTTTGGCAGCTGTGTATTGTTATGATGTGCGCAGTACAATGCGCAGTACATGCGCAGTACATAGTGAAGCAATTGTGCATCACTGGCTCATGTTTATTGGATAGTATCCTCCGGATGCGCAGTATGCGCAGTACAGCTGCGCAGTACAAATCGTCTGCGCAGTACACTGCGCAGTACAGTTCTACTGCGCACCCAAAGGCTGCGCAGTACAACATACCATGCGCAGTACAGTCGCCGGTCGTTGCGCAGTACAACGGTGATGAGCGAGAAAAATAATGCATGCGCAGTTGCGCAGTACACATGTAGGAATGCGCAGTACATCGCATAGCGATGTGCGCAGTACACTGCGCAGTACAATGCGCAGTACACCTGCGTGTATGTGCGCAGTACATGCGCAGTACACGCAGTGCGCAGTACACAGATTTGCGCAGTACAGAGTTCGTACATCACCCCACGAATGCTGCGCAGTACAATGCGCAGTACATATTCCTAACAGTGAGTCCTCGAGTAATGCGCAGTACAAGCTAGGCGTGGGGATGCGCAGTACAAACACCTTTGAAAGCATGCGCAGTACACCCGATATCACGGCAAACCTGCGCAGTACAGCGCAGTACAGGACCGCGCGCAGCCTATCTGTATTACGTCTCCATCTAAGACCTAACCGGGCATAGTCCGCAGTTTACGAGGGTATGGGTGGTTAAGATTCAACAGGATGCATCGACTTCCTCCTTAGATGCCTCCTGATCCCCTCCAGGGTCTTTATAAGCGCGGCCTTGCGCTGCCTTGGTCGTGTCCGCGGTTCCTGTTGAGCGGCCTAAGACCTGTTCCATGCCATGTTAAATCAGTCAGTGAACGCTCAAAGGTGTGTCCGGAAGCAAAACTTTTACCGCCCTCTGCAACAGAAATTCCCGATTCCGGGGATGACTGTATATGACGCGATACCGTGAGGCACCCCTGCCATTGGTGACATAGTAATCGGGGTGAGATGTGGAGATGCAAGCGATGTAACAACGGAATTGGTGCCCTACCACGAATACTGAATATTCACTGAATATTACTGAATATTCTTCGATTCCGGTGTCACTGAATATTCGACTGAATATTCTGAATATTCCATCCTAGAGACACCGATAGTATGTCCGGGCCAAAACTACTGAATATTCAATATTCAGTGGACCCGGTTAATAGTTATGCCGTATATATACTGAATATTCCTGAATATTCGGGCACTGTGTGTTTCAGCACTGAATACTGAATATTCCCGTATGATCGGGTACGAATACTGAATATTCCTGAACGGGAGCAACTGAATATTCCAATAGAACTAGACTGAATATTCGACTGAATATTCTGTACAGGCGCACTGAATATTCTAGAACTGGATTGATTACGATCATCCTGAAACGATTACGATCGAGATTACGATCAACGATCAGATCAGTTAAGGCACTGAATATTCACGGCACGATTACGATCATCGATTACGATCATCACTGAATAGATTACGATCACGCCATCCGGTCTAACGACATGTTACTCGGAACTGACGATTACGATCACCGAAATCAGGATTACGATCACACGGTCACTGAATAGATTACGATCAATTACGATGATTACGATCAAGATTACGATCGATTACGATCAACGATCATGTCGCGGTCGACCTCATAACTGTACCAGATTACGATCAGCTGATTACGATCAATCCGTCCAGATTACGATCATTACGATCATATAAGGCAATGAACACTATAGTGATTACGATCACGAACACATCTGATTACGATCATAAGAGGCCGACTTCTGATTACGATCACTGTTGAACGATTACGATGATTACGATCAGAAGGCAGTCCTATTGGATTACGATCACATACACGATTACGATCACAGTTGCCCCTGACCGAGGATTACGATCATAGCCAGGCGATTACGATCAATTACGATCATCTCCAACATGCTAGATTCAGCTTAGCAAG'
    k = 11
    l = 535
    t = 18
    # text = 'CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA'
    # k = 5
    # l = 50
    # t = 4
    print(bio.findClumps(text,k,l,t))

if __name__ == "__main__":
    main()