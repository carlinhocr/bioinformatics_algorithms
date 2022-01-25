class Bioinformatics(object):
    def __init__(self):
        pass
    def patterncount(self, text, pattern):
        count = 0
        for i in range(0,(len(text)-len(pattern))+1):
            print(i)
            print(len(pattern))
            subtext = text[i:len(pattern)+i]
            print(subtext)
            if subtext == pattern:
                count +=1
        print("It repeats",count)
        return count

    def FrequentWords(self,Text, k):
        words = []
        freq = self.FrequencyMap(Text, k)
        m = max(freq.values())
        for key in freq:
            if freq[key] == m:
                words.append(key)
        return words

    def FrequencyMap(self,Text, k):
        freq = {}
        n = len(Text)
        for i in range(n - k + 1):
            Pattern = Text[i:i + k]
            freq[Pattern] = 0
        for i in range(n - k + 1):
            Pattern = Text[i:i + k]
            freq[Pattern] += 1
        return freq

    def Reverse_v2(Pattern):
        return Pattern[::-1]

    def Reverse(Pattern):
        reverso = ""
        for i in range(len(Pattern)-1,-1,-1): #voy hasta el menos 1 para    que range incliuya al cero
            reverso += Pattern[i]
        return reverso

    def Complement(Pattern):
        complemento = ""
        for i in range(0,len(Pattern),1):
            if Pattern [i] == "A":
                lc = "T"
            elif Pattern [i] == "T":
                lc = "A"
            elif Pattern [i] == "G":
                lc = "C"
            elif Pattern [i] == "C":
                lc = "G"
            complemento += lc
        return complemento

    def ReverseComplement(Pattern):
        Pattern = Reverse(Pattern)
        Pattern = Complement(Pattern)
        return Pattern

    def PatternMatching(Pattern, Genome):
        positions = [] # output variable
        for i in range(0,(len(Genome)-len(Pattern))+1):
            subtext = Genome[i:len(Pattern)+i]
            if subtext == Pattern:
                positions.append(i)
        return positions

def main():

    bio = Bioinformatics()
    text = "AAATTAGAGTACCAAATTAGAAATTAGAACAAATTAGTCAGGAGTCGTTTGTTCGAAAATTAGCTAAATTAGATTTCAAATTAGAAATTAGGAAATTAGTCAAATTAGGGAAATTAGCCAAAAATTAGAAAATTAGAAATTAGAGAAATTAGAAAAATTAGAGATAATTAGAAATTAGCGTGGCAAATTAGGAAATTAGAAATTAGATTCTAAATTAGAAATTAGGAAATTAGAAATTAGATGAAATTAGTAAATTAGCAAATTAGGAAATTAGAAATTAGAAATTAGAAATTAGCGTATAAATTAGGAAATTAGTTATATGAGAAATTAGCCACAGTAAAATTAGCGAAATTAGGAAATTAGACAAGCTTGACGCGCAAATTAGTTAAATTAGAATAAAATTAGGTGCTCTGAAGTGCACTAAATTAGTGAAAAATTAGGAGAAATTAGGAAAATTAGTAAATTAGCAAATTAGGGCAAATTAGAAATTAGATCTGCTAAATTAGTAACCAAATTAGAAATTAGCAACAAATTAGAAATTAGTAAATTAGTATCGTTAAATTAGAAATTAGTTAAATTAGCAAATTAGGACTAAATTAGATCGTAAAATTAGCAAATTAGCAAATTAGGCGATATAAATTAGAAATTAGATAAATTAGATAAAATTAGTTAAATTAGTAAATTAGATCAAATTAGGACCAAATTAGTGAAATTAGGGAAAATTAGCCATTAAATTAGTTGAAATTAGGTAAATTAGTAAAAATTAGCAAATTAGGGTAGATTGGAAATTAGGTAAATTAGAGACCGCCAAAATTAGGGTAAATTAGGTAGAAATTAGTAAAATTAGGCATGAAATTAGAGATAAATTAGTAAAAATTAGCCGAAATTAGAAAATTAGCAAATTAGAAATTAGGTGCCAAAATTAGATGGATAAATTAGAAATTAGTGTAAATTAGGAAATTAGAAATTAGAAATTAGTCTGAAAATTAGAGAAATTAGTCGATTAAAATTAGAAATTAGTTGTGGTCAAATTAGAAATTAG"
    pattern = "AAATTAGAA"
    count = bio.patterncount(text,pattern)

if __name__ == "__main__":
    main()