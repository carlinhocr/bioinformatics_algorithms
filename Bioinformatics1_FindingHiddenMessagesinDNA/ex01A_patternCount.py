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

def main():

    bio = Bioinformatics()
    text = "AAATTAGAGTACCAAATTAGAAATTAGAACAAATTAGTCAGGAGTCGTTTGTTCGAAAATTAGCTAAATTAGATTTCAAATTAGAAATTAGGAAATTAGTCAAATTAGGGAAATTAGCCAAAAATTAGAAAATTAGAAATTAGAGAAATTAGAAAAATTAGAGATAATTAGAAATTAGCGTGGCAAATTAGGAAATTAGAAATTAGATTCTAAATTAGAAATTAGGAAATTAGAAATTAGATGAAATTAGTAAATTAGCAAATTAGGAAATTAGAAATTAGAAATTAGAAATTAGCGTATAAATTAGGAAATTAGTTATATGAGAAATTAGCCACAGTAAAATTAGCGAAATTAGGAAATTAGACAAGCTTGACGCGCAAATTAGTTAAATTAGAATAAAATTAGGTGCTCTGAAGTGCACTAAATTAGTGAAAAATTAGGAGAAATTAGGAAAATTAGTAAATTAGCAAATTAGGGCAAATTAGAAATTAGATCTGCTAAATTAGTAACCAAATTAGAAATTAGCAACAAATTAGAAATTAGTAAATTAGTATCGTTAAATTAGAAATTAGTTAAATTAGCAAATTAGGACTAAATTAGATCGTAAAATTAGCAAATTAGCAAATTAGGCGATATAAATTAGAAATTAGATAAATTAGATAAAATTAGTTAAATTAGTAAATTAGATCAAATTAGGACCAAATTAGTGAAATTAGGGAAAATTAGCCATTAAATTAGTTGAAATTAGGTAAATTAGTAAAAATTAGCAAATTAGGGTAGATTGGAAATTAGGTAAATTAGAGACCGCCAAAATTAGGGTAAATTAGGTAGAAATTAGTAAAATTAGGCATGAAATTAGAGATAAATTAGTAAAAATTAGCCGAAATTAGAAAATTAGCAAATTAGAAATTAGGTGCCAAAATTAGATGGATAAATTAGAAATTAGTGTAAATTAGGAAATTAGAAATTAGAAATTAGTCTGAAAATTAGAGAAATTAGTCGATTAAAATTAGAAATTAGTTGTGGTCAAATTAGAAATTAG"
    pattern = "AAATTAGAA"
    count = bio.patterncount(text,pattern)

if __name__ == "__main__":
    main()