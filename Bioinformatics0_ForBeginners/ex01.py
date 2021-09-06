class Bioinformatics(object):
    def __init__(self):
        pass
    def patterncount(self, text, pattern):
        count = 0
        for i in range(0,(len(text)-len(pattern))):
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

def main():

    bio = Bioinformatics()
    text = "AAATTAGAGTACCAAATTAGAAATTAGAACAAATTAGTCAGGAGTCGTTTGTTCGAAAATTAGCTAAATTAGATTTCAAATTAGAAATTAGGAAATTAGTCAAATTAGGGAAATTAGCCAAAAATTAGAAAATTAGAAATTAGAGAAATTAGAAAAATTAGAGATAATTAGAAATTAGCGTGGCAAATTAGGAAATTAGAAATTAGATTCTAAATTAGAAATTAGGAAATTAGAAATTAGATGAAATTAGTAAATTAGCAAATTAGGAAATTAGAAATTAGAAATTAGAAATTAGCGTATAAATTAGGAAATTAGTTATATGAGAAATTAGCCACAGTAAAATTAGCGAAATTAGGAAATTAGACAAGCTTGACGCGCAAATTAGTTAAATTAGAATAAAATTAGGTGCTCTGAAGTGCACTAAATTAGTGAAAAATTAGGAGAAATTAGGAAAATTAGTAAATTAGCAAATTAGGGCAAATTAGAAATTAGATCTGCTAAATTAGTAACCAAATTAGAAATTAGCAACAAATTAGAAATTAGTAAATTAGTATCGTTAAATTAGAAATTAGTTAAATTAGCAAATTAGGACTAAATTAGATCGTAAAATTAGCAAATTAGCAAATTAGGCGATATAAATTAGAAATTAGATAAATTAGATAAAATTAGTTAAATTAGTAAATTAGATCAAATTAGGACCAAATTAGTGAAATTAGGGAAAATTAGCCATTAAATTAGTTGAAATTAGGTAAATTAGTAAAAATTAGCAAATTAGGGTAGATTGGAAATTAGGTAAATTAGAGACCGCCAAAATTAGGGTAAATTAGGTAGAAATTAGTAAAATTAGGCATGAAATTAGAGATAAATTAGTAAAAATTAGCCGAAATTAGAAAATTAGCAAATTAGAAATTAGGTGCCAAAATTAGATGGATAAATTAGAAATTAGTGTAAATTAGGAAATTAGAAATTAGAAATTAGTCTGAAAATTAGAGAAATTAGTCGATTAAAATTAGAAATTAGTTGTGGTCAAATTAGAAATTAG"
    pattern = "AAATTAGAA"
    count = bio.patterncount(text,pattern)

if __name__ == "__main__":
    main()