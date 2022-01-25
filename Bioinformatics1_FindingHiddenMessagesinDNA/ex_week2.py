class Bioinformatics(object):
    def __init__(self):
        pass

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
        for i in range (0,len(skew)):
            if skew[i] == min_skew:
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

    def prefix(self,pattern):
        if len(pattern)<=1:
            return_value = ""
        else:
            return_value = pattern[:-1]
        return return_value

    def firstSymbol(self,pattern):
        return_value = pattern[0]
        return return_value

    def sufix(self,pattern):
        if len(pattern)<=1:
            return_value = pattern
        else:
            return_value = pattern[1:]
        return return_value

    def neighbors(self,pattern,d):
        neighborhood = []
        suffixNeighbors = []
        nucleotides = ["A","C","G","T"]
        if d == 0:
            return_value = [pattern]
        else:
            if len(pattern) == 1:
                return_value = nucleotides
            else:
                suffixNeighbors = self.neighbors(self.sufix(pattern), d)
                for item in suffixNeighbors:
                    if self.hammingDistance(self.sufix(pattern),item) < d:
                        for nucleo in nucleotides:
                            neighborhood.append(nucleo + item)
                    else:
                        neighborhood.append(self.firstSymbol(pattern)+item)
                return_value = neighborhood
        return_value.sort()
        return return_value

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

    def frequentWordsWithMismatches(self,text,k,d):
        patterns = []
        freqMap = {}
        n = len(text)
        for i in range(n-k+1):
            pattern = text[i:i+k]
            neighborhood = self.neighbors(pattern,d)
            for j in range (len(neighborhood)-1):
                neighbor = neighborhood[j]
                if neighbor not in freqMap.keys():
                    freqMap[neighbor] = 1
                else:
                    freqMap[neighbor] +=1
        m = max(freqMap.values())
        for item in freqMap.keys():
            if freqMap[item] == m:
                patterns.append(item)
        return patterns

    def frequentWordsWithMismatchesAndReverseComplements(self,text,k,d):
        patterns = []
        freqMap = {}
        n = len(text)
        for i in range(n-k+1):
            pattern = text[i:i+k]
            neighborhood = self.neighbors(pattern,d)
            for j in range (len(neighborhood)-1):
                neighbor = neighborhood[j]
                if neighbor not in freqMap.keys():
                    freqMap[neighbor] = 1
                else:
                    freqMap[neighbor] +=1
        count = 0
        for item in freqMap.keys(): # finding the max of item + reverser complement
            itemReverseComplement = self.reverseComplement(item)
            if itemReverseComplement in freqMap.keys():
                sum_item = freqMap[item] + freqMap[itemReverseComplement]
            else:
                sum_item = freqMap[item]
            if sum_item > count:
                count = sum_item
        for item in freqMap.keys(): # finding the max items
            itemReverseComplement = self.reverseComplement(item)
            if itemReverseComplement in freqMap.keys():
                sum_item = freqMap[item] + freqMap[itemReverseComplement]
            else:
                sum_item = freqMap[item]
            if sum_item == count:
                patterns.append(item)
        return(patterns)

    def findSalmonellaOri(self):
        with open('Salmonella_enterica.txt', 'r') as fileSalmonella:
            text = fileSalmonella.read().replace('\n', '')
        minimum_skew = self.MinimumSkew(text)[0]
        print ("minimum skew ",minimum_skew)
        oriWindow_start = minimum_skew-500
        oriWindow_end = minimum_skew+500
        oriwindowText = text[oriWindow_start:oriWindow_end]
        d = 1
        k = 9
        kmers = self.frequentWordsWithMismatchesAndReverseComplements(oriwindowText,k,d)
        return(kmers)



#-------------------------------------------------------------------------------


def main():

    bio = Bioinformatics()
    # genome = "GAGCCACCGCGATA"
    # print(bio.skewArray_v2(genome).values())
    # with open('dataset_7_10.txt','r') as file7_10:
    #     text = file7_10.read()
    # print(bio.MinimumSkew(text))
    # string1="GGGCCGTTGGT"
    # string2="GGACCGTTGAC"
    # string1="CCGAGGTAATGTGCTATGCCTGACGTTACACGGGTATAGGGGTTGTGTAATCCCAAGTCTATACCGACTGCGGACTCCGAGACCACACAGATGTAGCCAGGCCGCACGTAAGCATGGCAAAAGGGTAGACAACTCACCGCCTCAGTATAATTCTGTGAGCTCGGGAAGTAGGTGCCCGTTCAAAGAACTGGTGCTTGTCTTTCTACGGAGTCAATCATGATACAACGGAAGCCCTTTCGTTCACGACCGGACGGATGTGTTTTGTTTTAGGTGGTCGCCGTATCCCAGCACTCTGCGCGGCGCACAGGGTACCATTCAGGCAGCCGACAAGGCAACCTGGGATATAAAGGTACTGGTGCTCGGAATGTTTGGCATGCGACTTCGGAAATCCAAATCACATATTCGGTACAACGGCGAATCTGGTTGCGATCGCTATAGTTACGTATATTGCCACCAGCGAACTTATGCGGCGGCGGTCACGGTAGTTCTCGACAGGATTATGGTGAGGATCTGTCGCGCTATGAAATGCCAGGGGATCCCATAGGAGTCCCGTGTTCCTGCAGGCTTAGTCAGTACCCGATGAATCTATACGGGTCTGGATATACTACATCACCGCAACCATTTGTCTAGCAAGAACAGCAATGACGAACCAAAGCTGTTTAGGATAAGTCGCAAGCCAGAACAGCGTTCATCGTTCGGACGGCTGGACAAGGATAGAAAACCTTACCTAACGATATCCAGGAAGCGTCGCTCCGCTGATATGCCTAAGATTTTTGAAGTTTCTCAGACTGGTATCAGCTGATGTAAGACGACCCCAAGCGCCCCAACGGGAGTCTTAATAATTCGACAGCTTCCCCGTTCCGAATATTTTTCGTGGCCCCAGAAGAAAGTCGCGGCTTATATAGTCGACCTGTATGAGATTTGACTTTTCCAAAGTTCGACTGAGTCGAACGCGAGCCGAGTCGCACCCTGTTAACACAGCTGCGTCATAACACTTTGTAGTCCTGTACAAAACTCTCTTGACCCGCGTATATGTACATGTCGTGGCAACGCGAACGATGCTTCTCTGACACAGCTGGAGTATATAGCCCCTCTCGGCTCTAGTGTTCTATCGCTCTTA"
    # string2="CCATGTCGAGACCGCACCGGTGCGTCTTAGTGCCCATGTACTAATGCCTGTCGTTAGCTAAGACTGGCTCAATTCAGGTGGCCCTTCATTGCTATTCGACCCAGTGCCTGCGTACAGACAAGTCGATATTGTTAGCGCCAAAGTAACTCTGAGGGCACGTAGAGTCGTCACTTTCTGACAGCATGAATAACACTGTGCTTTCGGGGCCTTTAAAGCAGGACGAACAACGCGTCTAAGTGGTCTCCGAACCGCGGGCATACATAACGGTGCGCAGGGGGGCTCTTCCATGCCCGTAAATCGACACGGGCCACTAGGTTATGTGAAGCATCATTTAGCCGAGCGACAAGGTAGCTGCCGGGCTGTGGTGCCTTACAACTCGCCCGTGTTGCGGGTGGCATTTGTTGGACATGCCCCCTAACGAGATCGTAAACACACTACTCCTTGAATACGTCGATTCCGTTCCGCGCGGTGAAACTGCCCGCTCAGGAGGACCATCGATAACCGGTGCATCGCATAAGCTACCTCCGCGTGAAAACGTGTTCATGTACGCAACTATCGACTCATTAGCGAGTTGTCGAGCGCGCGCGTACGAAATCTTGAAGGTATTATCCCCGCGTCTGGGTGCACTCACGTTATAGTCTTTTTTACTGTGAGATCCCAAGGCGCCCACCAGTGATGTCAAGCAGGTGTCTTGTCTGAAACACTCCGTCGCTAGGTCGAGAATCATCAATCTCGTCTAGCGCGGCGCTCTGATGAACTTTTAGGGAGAGTACCAATCACGAGACGTGGTGTTACCCCGGCACAGAGCGTACATAGAGCCACCGATAGGAACATCTCTTTCCCAGGAGTGGTTTGGGGAATTGATCCGACATACGCAGTTATCATGGAGATTTATAACGTCCCAGCTGTGAACGCCTAAAATGGGGATCGTGTTCTCCGTTCCAGCTCTCACTTACTAGGAGACTCTCTTGTAGTCGATAGAAGACCCAAGCATGAACGTTCAGAGACCACAGAATCGTAGGCTTCTATAACGGTGTGCGCTGATTCGTCCTATTTTCGCTAATTCATTTCTTAGTAGAATCAAGGTTCACATTTCAGATCGCGGCCTACGCAGGCAGGC"
    # print(bio.hammingDistance(string1,string2))
    # pattern = "ATTCTGGA"
    # text = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
    # d = 3 # mismatches
    # pattern = "GCCGACAGATG"
    # text = "CTCTCTGGCGATTACTCGTTGCTCCTAGATTCGGACGCGTAGGGAGAGGTATGAGGGGGAGCAACACTGGTACGCAATCGTAAGGACGACGGAGGCGTTTTTAATCCATCCAGGCTGCGTATGCTGAAAGCATAGAGGACTAGAGATTTGTAACCAGTCGGACGTTCCTGCTGTATTAATCTGTGATACCCTCGCGACTTTCTATAAATTACGGGAACAACGCACTAGGGAGAGGCGACTGGCCGCTAAGTCTTGCGCCATATGAGCTAACCCACCCAATTGCCGATCCACATTGCCGGAAGCGTGATAAGGCAGCATCAATTAGCCTCACGCCGGTGTAGTGGCGGTGATCTAATGATGGAATCTACTGGTTATCGTTGGCTGTATAACCCCGAGGATTACGCGTCACGCGGTGACATGGATGCCGGGTGACTACAAATAATGGTTGCTGGGCGGGCTATGCATCTCATAGTCTAACCCGCCTACATAAAGGAGACGCTCAACGCGCCCCTACATAATAGCGGATGTTTTGGTAATTCGATTTAGTCATAGTTCTTGATCCCCCAGCGATTATGGATTCAGAGCCTCTCCGAACTACCTCAACGGGGGGATTGTCGCAAGGCTATCGTCCTCTTCTGCAAAAGCGAAGCGAGATCATTGATGACGCGCTTGGGCAGTAACCCCTCCTTAGCTCTCGGTAAATCTTCACAAACTAATGTCATCTGCTGGATTCTAACGCGAATTGTTAGTGCGATATTCGGTCTCCTCGACTGATGAATACCTGTCGCTCCGGACAGCGAATCAGGCATAAGATGCAGCTTTATATTTCTTCATTTCAGGGTAGCAAGCGGTATCTGAACGAATGGTTTCCACAGAAGACTTCTGTAGGGTACAGGGTTTGAGAGGGACGATCCGTCTACAAATAGAAACGCTAACCACATTATATGAACGAGACTACAGTGGTAGGCGTGGTTCCTTTTACCATACCTCGGGTATTACCACTACCTTTCCTGTTCCGAGGTATTATAATGCGCAAGTCGCGAGTGGAATTGGTTAATGTAAGCCAGTCCTGCAACTTGGGTCCTTTGCCAGGTGCCGCAACGTGACTAGTCGGTCAACCCCCGTATAGGTCGTTACCTTCAGGTCCGGAGCAATCAGTACGGGACTAAGCTATGCAGAGACCTGAACTGTATGGCACAAATGGCCACCGTGATGAGTCCCAAGGTGGACAAGCTCAGAAGTAGCCCATGTTTGATAAATTGAACGCTCTAACCTCGTTGAGATTTGAAGACGGTGGACAAGGAAACGCCTGGTAAAGGGTATCCGGTAAGTTATGTTTGAGCCCTAGGATGCCAAGCTCGCCAACTAAAGCCGCGTGAAGGAGGCCATGGGGGCGTCAATGACAAATAAAGTATGCCTAGACTTAAATCTAGTCGGACATGTTTTGGTGACAATGTCGTTCACCTCGACGCCCGATAAGGGGATGACAACGACGAATAGCGATCTAGTGTGTCCATCAGATAACACTGTACCCGGGACCCGTACGACTATGGATATCATCCGTTTGGTGCATTGTGCACCTTCCTACCTGGAATAACAATTCGCCTTAAAAACCCGCGTGTGCCCTGTTTGGCCATAGTGCGGATGAACCACGAAATAACTAGTACTGCGATGTTCAGCCAACTTCAATCAATGGTTAATCAGAGTGAATTAGGTCAAGTATGCCTTTATGAAACGAACACTGCTTCGACGCGATTTCGGTCACGGTCCGGCTGCTCACAGTCTGATTTTAAGCTCGCATCGTTTGCGTACATGAGAAGGTCTAAGTGTCCAAGGCCAGCGTAAGATCATGGTCGTGACAGACTCTCATGTTCGTAACAGCACAGATTGGATCAGATCAATTCATAGGTTGACTATGGCCGATAGCACAGAGGGGGAGGGGCAAATTGGAGCAGCGCATTGAGGAGGAGTCCGAGTTGAACAGATAAGGTGCAGATGCGTCGAAAATTGTACACAGCGGATCTCCTAATTCACTGTGAATTTCTCCAGGCGTGCGTGAATTATCTATTGAAATAATACTAGAATCGCGGACATTAGGCCCGTATAGTACTGCCACACTAGGTAAAACTACTAATGTCCGACAGCAGTGGCCTTAACGCCATACATACATGGATTCCCTCGCGTGATCTGACGTGGCTTAGTAAACTGTTGTTAATCGGTAGTAGGCCAGGGTTTCGCACGCCGTCACCAGAGGTACGAGCAAATGTTGAATATCTTTTAGTGTGTTGATAGACTGCAGCAGTGTGGAGTACTGAAAAGCTGCGGTAGATCCGTGTGGACCACGTGCCGTGGCTCATAGATAGGGTGCTTTCAGGCCTCGACAGAACGGATCAACATCTCAGCCGCCCGTATAGTGAGACGTCTTATAAAGACGCTCTTAGAACTCCATGGTTTCCTTGGTCTTTGCGGTTAGCACAATCTCGCCCTTAGAGGTAACTGTGCGCATAATCATAAACGTCCCCATGACAGCGACGAACACTCCATAACATTCCAATGTTGCTGTTGTAACAGGAGATATATTGTCTCGGAAGGGGACGGTACACACAATCACAGTCCGTTCGCAGTCCTTTTTTCCCAGCTTCCGTAAATTACGGGCTAACGTGGCACTTTGGTAGCCCGCACGAAAGAATGTGAGAGTCTGGGAAAAACTCACATCGCCTCCGTCAACTCAGTGTGCGTGGAGCTGAAGACCGAATGCTTCGCTCACTAGCGGGAGCATACCTGTGTCGTTATCCGGGTGTGAACACTGAAGGGAATCCATGGGCCGGTATATGTGTACAGGCGGCATTTCGGGGATAGTGTACAACGTTTTGCCTCTAAGTTGACGCCACTAGAACGCCCTTAGATAATCTTCGCACTAACTCCGGTTTACTGTGCTAAGCCACTATTTGCCGCAGGTTGTGGTGGAACGTGTGGCAGGATAGGGCCCTGTGCGAGAGTCCGCGATTGTGCCTGGAATCTCGATAATTGTCTATACCTGCCGACGTGTCATACCTATAATTCGATCGCCCGGATCCCATCCGGGATAGTCTTACTGCCGTTGGTTTCCTCCTCGGAGATCAAGTCTTGCCGTAAAAGGAAATTCGACATGCGTCTTAGTTAACAAGATCAAGCGAAGACAAGTCCCTAGCGTGCTCGTATCGGAACTTCCGGCGAACTGGTGGCCCAAAGGCCATACAACTTGGTCTCCAGAGTACGCTTTAAAATACAAGTGACAGCACGGGGAGAGGAAGCTTGTTTCCCGACTGATTGCGCACTGGATGCCCTAGGTTGGTTCAGGCGTAACTCAGGGTTCTCACTAGGGCTTTTCTTGGCGGACTATCTGGGGGGCGTTAATTCATGACTGTGATAATTTATGCCTGTAACGCAATATGCGATTCGACAGCGCAAATTCCCGCCGCAAAGATCCAACACGCTGTAGACATATATTGTTGATCTTACGGTCTTTGAGTGTTTCGTCGTCCTGAGGGCTAGGCGAGCGGCGCAGAATCTCAGCCAGACAACACCGTCTCGGATGTTCGTCCGTGTCTAAGTCGTAGTACCAGAATTCGGCCATCATCGCCCTTCTCTTGGACTTAGGAAATACACATGTACTCTAGTAATAGAGGTGATTTAAAATTAAACTTGCGTAGCAGTTGTGAGCCATGGCCTGGCGACGTTCTCATTATGACGGGCATGTAGGAGCGGTAGTGTTCCGACTAGAGCCATACGTCAGGAAGTGTGGCCGCTCCATTTTAGAGGTCCCTACGTAAGCAACCCGGATTCCACGATTAGTGTCTCCCCAGGTACCTGTGTTAGGGGCTTCGCTCTATAATGCCAATCTCGTAAGGGCAGGGAACAACTTAACTCTGTAATTAGATTCGCGTTCGAATTTCCCAGCTAAGAACTAACTCATCTCATGTACGTGGTATGAACAGGGGAGTAGTTTGAAAGGGTGTATCCGTTGGTAGATATCCTCGTGCAATGCCGTTTGGGGAGAACATGCGAACATCTCAAATAACCCTAAAATAGGGATTTGAGGTTGGGCGTGCCAACCGATATCTGGTACCGGCGTTCGCTAAGGGATGCTGCATCCCTAGGTTGCAGATTAAGAGATTTACCCCTAACTATTGTGCTCGTTGTGCACGCGGCGAGAGTAGCGGCCGTACGCCAAATATTGAGGCATAAAGAGAACTGAATTGACCAAAGAAACACGTTTAAAAGCCCTCAATAGTGAGAGCTTCCAGACGTCTATCGAAGAATATTCGGACAAATACGGTTTGGTTGGACACGACGTAGTCCAGAATTAAGTTCAGCTGAGACATGCACGGAATGACGGAACGGCTGTCGTGCTTGTCTACTCGCTACACGATCTTCCGGGCGGACAGCAAAGTTATGAAACACCCTTAATATCTCTTTCTTGGCCTGCACTACATTGCTGTCGGACGCTACCCTGCTAGACTAATCTTGCCAAGAACTTCGTCCCCAGACATATATTAGGTGTCCGAATGCTTGTATCGCACGTAATGGTTCGTGTCACTGAAACTGCGTGGCTAGAACTCCATATGCCGGGCAAGAGTACACAATGCTTGACTTTGTCCCACTCCGCCGGGCATAGGCTAACGGTTCTGTGCGTTCGAGCGATTATGAACATAGCCCAAGGCATATCAGTCCCATCGGACGGGGATAGCGTTGTACGCAATGGCGTATCCAGCTATAGGGTGTATGTTGCGAATCGGAGGCACCGCAGGCCAGGCGCGGGTGGTAGTTTACACGGTTTACCATTTGCCACCAACGCATAGTACCTTGGTGCGTGGTCTAGAGGTGCTCGGATACATCAAAGTCATGCCACTGTACAGGGTCGCGCGTTAATCCAAGCTTCGCTGTTCTGTGAGTTGCCCTTAATGCGATTGCAGCACTGGTCGGAGGAACGGGCAACTAACTGGGCCATTGTCCATAAGGACTCAGGAACGGGTTCACCAGTCGCTAGGGAACTGTTAGAGCGACACCGTGTACACGGAAAATCCGGCCCGAATAAGTGTGCTGAGACCACCCAGTGACGACAGGACGATAGTGGCAGGGAGTAGCGATATAGACGGGACCGTCGCAATGCGGTTGGTACCGCTATCAAGATGGCGGGTTCTATGACAGTAGTCGGAACTAACTTTGCTAGGTCCCTCTACCACCCGTGACCAAGCTCTGTGTAACCTCAGTGTTTTGAACAGATAACGGCCTGACACTAGCGTTACGGAATAGGTACCAACCGTCGCAATCAGTTGTAGCAAGGCGCGTGATCCAACACCCTCGCCTTCATGATAGACCGTTGTGTTTCCTAAACGACAAACCCCGTATAGCCGGAGTTAGTCAAGCCTATCGCGGATGCCTTCTCCGCTGATCCTCTATGAGCCGTACCGCACTGAGCGTAGACTATCTGAACGCAGCAACAAATTGGTCAAGGACCTCCTTTTACGCCTCCTCAGTGAGGAGTAGGAGCTCTTTAGCCATGGACATTCAAAGGCCCACGGCAAGATGGGAGGCCCATGGGTCCGCATAGCGGCACTCTAGAATATCCCCTCTGCCACAGGTGCCCCCTGCGGGTGACTCGATTTTAGTCAAGCATTAGCACATAGTTATCCGAACCTATGCCCACCCCCTAACCGGCCCCATCAATCCACTTTCAACGTATCCCAGAACGGGATGCTAGTCGTCGTCCAACCATTTAGCTCATAATGTGCGAATACGCGGAGTTGATCCGTGTGGATAACCACATATTTCCTAAAGAAGTCCAAGGTAACCCAAGCCTCCAACAGTTTGGTGGGGAGTTAACCGACGGATCAATGCGAAGCCCGTATGCGCCACATTCAACTTGTGCCATGTGCAAGTATACATTCTTGTATGTTGACGTCGACTGCAAATTGCACTCGAGGCTCCCAGAGCGTAGGTCGGGTTTTTTACTTAAGATGGCTAATAGGTGGTTCCGTAAATATGCTTGAGGCAGTGAAAGCAAAGATATGATCAGTGCGGCACCTGACTCCAGTGCTGTTGGTGGTTTGGCTCCTCTTGGTCCCACGCGCAAATTTAGTCAGATGGGGCCTACCCGATGTCTGAGAATAGTTAACACAGTTCATGACGCAATTCCTGTACCCTGCTAGAAATTATCTACCCTGTACATAGGTTTTTAGCATATATTTTACCAGAGTGTGTTTCAATTGTTCTAAGATTGGCGAACGGAGAATGTTGAGTAACCAAAGTTCTTCGTTGCACCAAAATTTCTACCGTCATTGGATTCCAACTGAAAGTTTCGTCTTGTCGAAACAGGTTGGCCCTTCCCTACTTCTAGCGACCTTTAATGAGTCGGGGTACGTAGTATGGATGTCATTGTACAGTGTAGCCCGTCGCAATCGCGAAGGATCCTGATGTTACCGCGACCCTCGGATCGACTCATTCCAGCAGTTAGATCCGGATAACCTCAGACACGTTGCCCTCTGTCGAAGGCTAATGACGTAATGTGTAAGGTGGGCTTGATGAACCGCCGGGAAATATTGTCCCAGTTTCGTGCCCAAGCTACCACAGTCGTCCACAATCGAAAACGTCAGTGATTGTCACATTTGTTCACGTGCTTACGGTTTAAGGTCTCAGAAGAGCACAGCTAACAATCAAAGTGACAGAATTTTATAGGGCGTTACCGATAGTGCATATTTTTGAGAGATGCATCCGTACTGGGACAGCTGAAACTAGTCGTTTGAGCACCATGTGGAATACAGTTCCCAGAGTCGGACGCCTCCTAAAAATAGGTTAACTCCGTGCAAACCGGTGAACAACTCCTTTGGGGGCGAGTGATCGTAAAGTGAGGCAGCAGCTTCTTATGTACTAAGTGGTTAGTAGCCTCGTTAACTTGGCTACCGTGGTCCGTGGCCCGGCTCGGGTTAGTCACCCAGGGTGCGCCTTGTCTGCAATCATCCTTTGTAAACTGCTTCCTACAGTGGCCGTCAGTGAGACTGCCAAGACATCATGATTGAAACCCTTTGGCCCTAGATAACGAGCATAGGCAAGAGGCGAAGACGGCTGCCAAGCAACGACTAACTCTTACAGCTTTTCAGGCTTGGAATAATGCATAGTAATCCCTCCACGGGTACAATTTTGTCTCTTAGTGTTGTGCGACCCGCTTCCGGCGCAGGCCTCGCTGTTAGCGGGAGACGGGACTCGCGCTGTGTTAAGCCCAGCTGATGGTTCGAGACAGCCCTATGGGCAATCAGTTATACCTCAATGGTAATTTCACCTCGACCCTTCGCCTAAGTCACTTACATCCTGATCAAATTACTGGACTTGTATCCTCGCTGTGAACCTCGATGCCAATTTGGCGTATTTAGTGGGGAATGCCTCCGAGTCCATAGTTAGAATTTTTTTCATGACAGACCACTGACCGAGGCCCGTTGTTGTCCTGTATTTGTGAGATTCTCGTCAGAGACGACCACAGTAAAGTCCGTATACGGACACGTGGTCTCCCGGATTCGCAACATCGGCTCGTTATGTGTCGACGACGATATGATCTGTCTGGCATGCGCTTGAACACGCCCCACCACTGTCCAATTGCCACACCCGCTCTCCTGACAAGGCTGTTCGTGCCAATACCCCGTCCTTGTCATCTTACACAAGGTAGTAGTATGATTCGACAGCTGCACCGAGACATAGCAATGGCTCAGGCCGGCGAACCACACCACGACATTGCTATTCTACTGATGGAAAGGCCGAACCCCAGAGGCGGCGATTGGCCTTCCATTAACCATCTGCCAAATTTCAGCCTCAGAGGATCTCTAACCGTAGCTACTCGTATTTATATATGCCTCACAATTACACGTATCGAGACTTGGAAAGCCAAAAGCAGGGTTGAAGTTTTCCTTGATTCTAAATAGGCAAACCATGTGACCTCTCTACATTTACCTCTGAAGCCAGAGCGTATGGACCGGGGTACTGTGAGCCGTAGGACAGGCCTACGCCACGGGAACTACGTTATAGTGGACAGGTACGACTTCGCGGGTTACACGTACTCAGCAGCGGCGTCGGTACTGATTCCTAAAGTGTAAAGTTGATTGTCAATTCTTATTCCGTATAGTACAGACTAGCTACGTCAACGACGATGTGTAGCACCTGAGGCGCACCGATTTCATAATAGCTCCCACGAACTACATCGTAATAGCTTAACAGTTTTGCGAAAAAACTTCTCGGCATGTCAAGTTATAGTTACTCGCACTCGAAGTTGAAACAGGGCGTCGTTCTCCGCCAAATACCTGTCCTGAAGGAAGGGGTCGGAAAAGGGGTGGCTCATGGCGATGAGGACACAGTGGGACTCTAAGCAGTTTGAGAAAGACCCTCGGGAGCACGCGACCTTTGTTATACTGTCATAAGATTATCGGTCGCTGGTGATTTAGACGCCTTACCTTCCAATATTAGTCTACGTCCAATGCAGAATAACGGTAACTCCACTAATAGCACCGCGGTTTTGGTAAGGGCCGTCTGTTGTTGGCTGAGTGGTAGGTTGTTTGACATAAGCTAAGGACCCGCAGTCCTACAGCGTCGGACTCACCATGCTCAGTCTAAATTTTGTAATCCACGCCCCCAGCAGGGTCGTGCGCGTACCTTGGGGACTTAGTCAGGCCTCAATAAGGTAGGTTCCAGCGCGAAGCGTAACTTGCCATGAAGCAGCCACGTACCTGAGGTCATTGGTGCGTTCCGAACGTCCTTCAAGTAGACGGTTGCTGGGTATAGTCGCGTGTTATAGGTTTAGGATCGCCTTTAAATGAGTTCAACTGATGAGAGCGGACCACCTGGGTGGCGACTGATGTATCCGGAGCCTTCAAGAATACGCGAAACGGGGCCCACGGGGACGGTTGACTACCCGATTGGAGAAAAAGGTTAGCCGAGTGGATCAAAAGGAAGAGGTAAATCCTTCGGACAAAGGCAATAGATGATTACTTAGTGGCTTGGTCACGATTAGTGTTCAACATATGCCCACGCGACCACGTTGCACGATAAGCGGAGTTCCAGTGGATACCCTCCCGGGGTAACCCGGCGACCATACTAACGAAACCTGAACCCACCCGGCTCGAATAAACAATCGACCATAATTTTCCAAGCAGTCGTTCAGCCCCTGTCTACATATATGTTTTTTCATACAGCCTTCCGCATGGCATTTGATGTGTCTGCCACCTGTGTAATAGATGACAGTCCACTAATTTTTCTAGTGAACGCTCATGAGTTTGTTTAATGAACAGAAAGGGAGAGATTTCTGCTAGCCTGAATACCCTTATCAGTAACCGACTTCTCACTTCGATGCTCCAGGGGAGGATTTTACGCCTCGATCTCGGATTCTATGGCGATGATACCGGTATATCCCGCATCAGGGGCCCGGCATTTATAGTCTGTCTCTAAACTTTCACCACTAAATTCCCTGTTATTAGGTACTCGGGCTTACTCAGGGGGGGTGGCGACGTGCGCTTAAGTCATTCGACCAGCTGTATAGATAATCGCGCGTTCAATATAAATAGTCTTCAGTGTCTATCTATTTCAAGAAAATAAGATGTGCTATATATATCTGCTCGACAAGGTTCATAGGTAACAGAGAAGAACAAGTTAGCCTGGCTAGCTTGTCGATATCGCGAGAAGAGTCGAGATGACATGCGGCGCCGTCCTAGTCTCAACCACCATTTGTAGGGTTAAACCAATCTGACAGGGGGCGTGTGGAAGGTTTTCCGAAGTGGACTGACATTGAGATATCCTAGTAAACTGATCGAGGCAACACGAAGTTTAAAATCGCCTGCTATGCTTTTTTGACGGATTACGGTACTATCATTTCTTAGTCCACCGTATGAGCACGAGGACTTCCTAACGAACATGCGTCGTGCCGAGTCCAGATAGTTAAAACATGCGATACAATGATACGAACACTGAGAGATCCTTTAGGCCCTCAGTTAACCGTACGTTCCAATCTTCACCCAGGTAACCCCGCGCCCTCGTAGCGATATTCATGTTCGATCCCGTGATTGTGGAGGCCGATTCCATCGACCCGTCGAGACTCTACATATAACCTACTAAACGCCTTGTAGGCTTGGCGGTCTTTGTGTGTCCCGAGCAGCAGTAGCTTCGAAAAGCTCGGTGCGCTCCACACGACAGTGACGGCTTTGTATTTCAAATCAGTACCCTAAAAGTGTACGGATCGACAATCTGCCCGGTCGCTCTTAGCCTTCAGGTTGTACTAGCGGGCAGTCAATACAAGGGAACGGTTGACCCGCGGGATTGTTGCGAGTTAGCAGCTAGTGCTGATATGCACCCTACGGTAGCAAGATATAGGCTACCACGTAATGGAGTGGTCAAAGCGCTGGGGGTTCAGATTGTGCTTGCGCCGGCTTACGATACCGCAGGATTGCTACATAACAAGAGCGTGCGATTCGGTGGCACCTGTGGCAGGGTGGGTAGACGTCATGCCCGTATACAAAACACAAATTCGATCTTCCCGGGCTAGCTGCGATTAGGGAAAAGTAGGCCAAATGCTCGTAGTTGCCACTTACATGGCCGGCTATGTTTCTGTTTCGAAGTTATCCGGTTGACCCTCCCCCCAGCCCATTGTAAAGAATATTTCCTTGGGGGCGCACTGGGTTACTTTCAGAGCCGGCCGGCTACGCTCAAAGGATCTGCGAACCATAACTGCGCGACATGGAAGTCTATCCATTGTGTTGCAGTAATTCAAGGATGGAGGCGTGAAGTACAAATTGCTAATACATCCTGCTTCGGTCGTCGCTAACTAACTCAGGACCCGCGGAGCTCGCGATTTTTTAGAGGATACGTGAACCCACTGGGGAGTACAACCCGAGACCGCTGCCCGACTGACTTGCTTTGTCGGGTAGGACACGGTCGGGCCAAGCGCGATGCTTTGATTAGGTGTTAAGTCAGCGAGTCGGGTGGTGTCCGAGCGCGACACGCTTCGGTGACAGCGTCAGTTCGTCTTCGTTAGCAGCATGAGTGCCGCTCCTGATACCTAGAGTTGTTCCGGGCCACATTCCACACTTTACGCGGCTCTCTGGGTCTTGTACACTTATTCCTGCCTCTCGGGCACCAAATAGGGCAACCATCGAACGAATCCTAGCTCTTGGACATCCCAATTGCCACCGCAAACTGTGCTAAACAGCCGCTCCGGTAGGGTGAACTGCACACGGCTAACCTCGACTCGAAGACAGATACGAATGTTGGGGAAAAGGGACCTCGATGGGGTGCGGCAGACACTCGGCTAAATGTTCTCCGCCCTGTATGATTCGGTATAGAGCATTTTCAGGCCATGCTTACAGTTGGTTGGGGGCTGTCTGTACTTAGAAGGGAGGTAGTTTCGACTTTGAGCTGTTCTTGTGGTCTGACTAAACATTCGGCGCAGCGAGGACGCAGTTTGTTCAGTGTCCATTCCAATACTAGAGCTGCCGTCCAGTGGAAAGAAGTCTTGTCACAGCTTTCTGGGCGCCAAGGGCATGTGATGATCCAACGAACCGGATACATTCAGGCGCCGTACGTGTTGCTCGGACGAGTAAGAGTAATGATTCACCCCAAAAAAACAGTCGGGGGGACATGCCGACGTGCACCTAATATCTACTTTGTACCCGTGCAGAGAGCGGATGGCCGGCACTTCTATTGCAGCGACCACAGATAGCCCGGTTACTGAGGGGGCCGTCAGTACAGATAACCCTAACGGGCTCCACAATTCATCGGCGCTTGGGAAATCCGAAACTGGGATGGGCTAGTTCCCCCGTAAAGTGCTTGGCCTCAGCTAAGAAATATAGTATACACTTCGCCAAATCGCCTTGCCCAGTCGGCTTATTGTAGTACTGGAGGCCTAGCGTATGAAGACGCGTAAGTCAGGTGGGATAGGTTATGAGTAACAGGGCTTCCCACCAAAGCACTGATGCTACCGATTTATGGATAACTAACTGACACCTCGAAATGTTAGGGGTTAATTCTCAACCTTCCGGGCCTAGTAATCTACTCAAAACGTCGCGGGTCATCTGGGCTTCAGCGACGCAGAGGCACCACGACATTATATCACGTTGTAGGTTCAGAGGCATGGGCGGTAGTTCGGCTTTGGACGTGGAGGCATCGCACACCAGTGCGGAGCTTTTAGTCTGAACAATTGGTGCAGTTATTTTGTCCACAACTAATACAGTCAGATGATCTCGGATCCTATTGTCTAAAACCTTAGAAGGTTCCTGACAAGGACGCAGGCAGCGTTTCTTTGAAATCGCCTCAAAAAGCTAGCTTTTTTCGCGAGTACCGGAAAAGTTATCTCACTGCGGGGGGAACCCGGATCAATGTTATCTTCATGATACTTGAGAGATGGGAGTATTGAGTCGCAGCGACTTCTATCGCAGATCTTACGGCAACGGAAAGTCTCTCTTTAGTGGACGTATGCCAGAACCTTCTCTGGGGCGGTAAACGCCACCGGTCGTTTGGACAGAGAGTGCCTTCACGCAAGGAGTGCATACAAGGGTCACCCTGGGGCTTGAGGAATTTTTCAGATAGCAGAAGTACGTGCGTCGAAAGCGTGGACGAATATCGACGACCCCAAATCATACAACTCTTATGGCCAAGGCCGCATAAATTTAACCAGTACGTCGCCCATCTTCAGCTGTCCTCATGTTGTTCTACCCTCCGGGAACTCCGGCGAGAACCTATCCTAGCTGGCAGGTCCTGTTATCAGGCCTAAGATGATGGTTCTGGTTTCGCCGATTTCACCCCCGCGAGGTGCGGTTTCAAAACGGTTAATGCTACTTACGAAGTCGGTAGCTATGGACCCACATGCGTGCGCCGGTGGGCTGGTAGTCTTCGGCGCGACGACAATGCCCTAAAGGTCAACTTACAGGTACGAGCATGGGTCAAGTATTCTACGATAGTTATACAGCCTAACCCGTGAATTCCTGGAGTATACATGTCAAAGTCGTAATCTTAAAGCACATATCTTTGCTTTTTTGTGATCCATAATACCGCGCGGGGACAAAGGAGACCAGAGCCTCTCACAACGTGGATTCGCACGGCCATATGTCAGACCGATGGTCAGGGATCTTTCAGAAATGGGGGTAATTGAGTAACGCTCCATAGGTCAGATTAACCCTGCCTGTCGACCACCTGCGACATGCATTTAATTTTTGTAACTCGGCGGCTTGCAGAAACGTCCTGACGTCCGATAACGGCCCTTTTGTGTACAAAAAAAAAATACAAGGAGTTTTCGGGGTCCCGGCGGACATGCCGCGGACCAGCTTATGTATCTGGCCCTTGAAGGTTCGCTAGAAATTAGCCCTGCTGCCTATCATCTTAATGTGCGTTGTACTACATACATGTTGGGAGCTAGGTAGAAATACGCTATCCGGCTACACGTAAACTGCTTGTTTTCCAGGGCAAAGACGGGCGTGACGCTAACCAGAGCTGTCGATTTAGCGGCTGCTAAGATGAGTCACGTCTCCCTGTTTACACATAAGTTATGGTTACGGGGCATCGAGAGGTTACCACGGGTAGATGCCTTAGTTGACGACAACATGAGGTGCTCGCGACGGGTGCTTTAATTCAAAATTGGAACATACTAAGAGACACCGCACGCTCGCAGTAGAGCACTGCCGCACCTCTGACGCTCGCCAAGGTAGGACATAGATCAATGTACCTGACGACCGCCTGGAGACATCTGAAAGAGGAATAGATGGTTAGAAACTCAACGCCCTATTAGAATGGCCAGGGGGGGGGAGAGCAGTTTTGAGAGCATTTATGCCACCTCTTTAAAAAAGGGTCAGTGCGAACTAGACGGGATGGTCCGATGATTGATGCGAAATCTTCCTAGCACGATGCGCCTAGGCTCCTGGTCCCGTTCGAAACGACAACAAAACCGCGATCAAGACCTGAACGGTAATTGCTATACGCGGACAAGCGCATGACCAACCCGCCCTCAATGCTATGGGGGCGCGCAGAGCGCACTTGTGACCTAAAAAATGTAGCCGCCTGTGGCTTTAAGCGGCGGGTTCTAAGCCGTCTCGACTGGCCAATTGGGAACTAGCGGGCGGAGGGGGTATAGTCCGCGGATATGGACGGACTTCCCCAGAGGCTAGGAGTCCTCCTGAACCGCAGACAGTTGTTAGTGTACTCTGCCAGTAGACGGCTTTCGCGAACGGCAATAGCACCCAAGTCCTGGGTGGATTGGACCGTACTCGAATAAGGTCCACACCGGTTTTTAATTGAACGGTTTGCATACCACGTTCTAGCTCAAGCCACTGATAAGTATCCAACCTGGACCGGCTGCCCACACCCCGCCGACATCGTTACCCAATAGCGGCTGGAGGGGGCCGTTGCATAGCCCACGTGCTATCATACCAGACGGACTATGTTTGAACAGCCGCTTGGACGTGGCGGCAAACCCAATTCGTGACCTGGGGCCAACAACATACGTCTTTTTTGACCTAGTCTCTGGGGTGACAATCGAGAGTGTTTAAACTCTCACTAGTGCTTAATCGTACGAAATATGTCATTCCGTCTAGGTGGTCCTGAGGCTCGAGTGTGCTGCGATAGATGCCGCAATCATTCCGATATGGTCCGACGAGAATGAGACTAGCGTTTTCGGGAAAGGTGGTATTTACCTATACTCGCGCAATCTGAGGAACTACCCGTAGCGTATCTGCTATTTAAAAACATGCAAGACGCCAACTCACTGTTGCTATACACTACCCCACCCGTCCGGCCCTAGCTGCCTTCAGTCGGTGCTGTGACACAAGTAACAGGGACGCCAGCGGGCCCGCCAGGAAAGTTGCAGAGTAACTGTTCCCGCCTTGAAATTGCCCGTGCGTTAACGCGTCCAGTAGCCGAGCCTTAGGAACGGTTGACAACGACCGGCCCCTTGGGAACATGTTGGGCAGCACGAGAACAAAAGATCTGCCTGAGACCCCGCAGGTACTCACCCTCTGCCCGGTGGATGGACACGGCCACCGAAGTATATGAAATCGAAGCAGATTACCTCATTATTGGTATCTTGCAACCTCTAAGGGGCTGGCTGTTGTGGAGTCTCTGTTCAGAGTCTTACCTCAAATTGTCAGATTTTGCTACCAGCTATGCCGTCCCGCGGCTTTATCTCGGGCGTTGTCTACTAGGAGGTCCGTTATCAAAGAAGCACAAATTGGAGAGAATTGGTTCCAAACTTTCGGGAGTGGACAAAATACGGTTCGTTGAAGACCGTTTAGAGTGACACTTATACATAGCCACACCTGCACTATTATTGCCTGAGAATCAGCGTACTAAAATATCCGCTTAACTAGCCGGTCTGGCTCCAATGAGGACTACTATACGGTGTCATGTGTAACGACGTCACGTAACAGCTCTCCCACTTGCGAGGTACAACATTGTTGTAAATTCGGTGAAGCTTTCTACCCGCATATAAACGAAACACACTGAGTGCGTGACATACGATGGTAACTACGTTTCGCCACGGGACGCGAGGACGTTTTGGAATCATTAGTGATTTGAGTCGACTCCCTTGCACATGGTGGCACAACGGAATAGCGTACTAGGGGGGTAGGGCACACTAACTCGATTGTAGGTGAGCGTCTAGCGTGCTCCGCACTACAAACAATAACCCGTAGGACCATTTTTCATGTCCTACGAAAATGGGCTTCGTGCGGCCCGACGGCATCCGACGTCACATATCTCCAGGAGCCCGGAGAAATGCGGGATCGAGCCAAAATTCGGAACTTGGCCCTCTCTCCTGCAAATAAGTAACCCAAAAAAAAGTAAAACCAGAACTAGGACTAACCGCCTTTTGCGTGAGCATTCACACTCCATCCGTTTGGGGAGCAGTGTTGGCCTCCATAACCGGTCATCTATTGCACAAATGGATGGTTACTAGAAGTTACCTGATCACGGAAATTAAGCCATAAGTCCTTAGAGTGCCCAGATTCATGTCCTTACCGTGGCAGGCTACTCGGCCGCAACTATTTGCATGTAATTTTCTAGTAGATCCTGCTTGATCTATGACGTGTATCCGGTCTTATGGGCAAGCCCGTTGCTTAGAAGCTAGTGCAAATCAGTCGTGGTACGCGTTTCGCTCCACAGACGTTTGCGAAGTAAGATCTTTCCCTCAGGGTTTCCCAATTATGTACTAGGAATAAGAGTAGGACTTGATATCAGAGACGTCGATGGTAGGGCGTGGGATGCGAGACGTCGTCTGCTGGAACGCGCCCTTATCGCCGCGTGACGCTCACGCGTTCGAACAGATGACAAACGTCCACTTGGTCGGCAGACCTAGTCCAAGGCGGAGCGCCCCCGCACCATCAACGTTTGCAGCGAAGTACCCCCTGTATCCTGCAGATTGAGGGCCCACAGACGCGATCCAAACTAGTGGCGGGGTCGTCGGTTGCAAGTTGCTGGATTTATCCCGAATCACGTGATACGTAAAGCTCTCACGTCAATGCCGAGAAGGATGAACAAGATGATCGCTAAAAACTGCGACTAATTTCGAGTAATTGTCCTCTCAGTAACGCTCTTTGGTTTTCCTCGAATCATTTGCAAACACGCCAGTCATAACTTCCGCTAAACTGGGCATAAGGACTGTGTGTGTGGTCCAGCATAAATCGGTGCCCAAGGAACGGCTAGACGATACGACCCGGTAAGAGTACGTTCCGTGGCTATTTCAGTGGACCTACGACAGGCGAATCTCGCGAACCCGGGGCTATGACGTCCGTTGAGCTGGAGTACCTGCTACGTCACGGTACCAGTAATCAGACTAACGCGGCGAAGTCGCGCGTGGGTGAAGTTGCCTAGCGGCCACTGACCGACGGGTGAGCGTGGCTGAGAAGGTAGATTTTGGATACCGTCGAACTGTTGCTAGCCTCTTAGCTGTTGGGCATACAGAAACTTGAACCCAGGCTTTACGAATACTGGAGTCTGTACAAGGAACTAATGCGCCACGACCAGTCGTCCTAGCAAGCCCGTTCGCCAGTGTGGGCCTATCTTACAGACGCATGCCTGCCTCGTCGCTCAGTCCTGGACCGCCTAACTGAGGATCCTAGACACACGTGGACCTCAGATCCGTAATCATCTAGGCTATCGTGGGCCCGAAGCGTCTTTTCAGGATTGTTCAAAGCTCTAGAAGCTTTAATAGAGCGGCTGACGAACTATTGGGGCGGCGTCTGGCCTTGTCCCTCGCCTCGTAACTCTGACTGAGCTGTCGTGTGTGTTTCTGTACTCTGTAATCCTGTAGTCCCAATGATACACCACGACGAGCCGGGCTAAACTCGATTGCAACCGCTCAGCGTGTAGTTTGAGGGTGCAAGCAGGGTAGTTCTCCGGGTAAACGCTTGGGGGTTTTAGCATATGATGTAGGTTCGCTTCCAGTGCAGTTTCATACAGGTTCCATTCGATACTCCCGTGGACAAATAATGGCGTACTAAATCCTTAGAACCGTCATGGCTCTGCTACTCTGAGACTAGTCTTCCATACACCTCCAGATATGGTTGATAAATAGGGTTCTGAAATCGTAATCCTCAACGTTGCGCTGTATGCCGTCAGGCGTGGTTTCCTCCTACATTCGGTTGATGGATTATGTTGATCTAACTCTTAAAATTGCAATACGCATTTTATTCGAGGCCATGTGCCGACAGATG"
    # d = 4 # mismatches
    # print(bio.approximatePatternMatching(text,pattern,d))
    # text = "ATGCTTTCTCGGTTTATTATCCCGCCGTCAGACGCTAGGGCTTGCCCGCGTTCAAGGGAGGTTATGTACAACGGCCCTCAGCATGAGTGTACGAGTGGCTTTTTTAGCCCGGGCCCGACGGATTGACATTTTGAAGAATGCAGAGGTATCAGTTATCTTGCTCTGTTCGACGTTTCTCTGCGTGCATGATCGTGAGTATCGACGAGATTTGAGGCCACCTTCCTACGAGTGGCGTTAGTCCCTGCCTTCCTGGTCGTTGTGTTCCGTCCCAGCCGGCGCTCCCATGAGGCTCAGGATAATTCGTAGCATATAATAC"
    # pattern = "TGAGGC"
    # d=2
    # print(bio.approximatePatternCount(text,pattern,d))
    #print (bio.inmediateNeighbors("AA"))
    # pattern = "GTGCATCC"
    # print("resultado final",bio.neighbors(pattern,2))
    # print("resultado final",bio.iterativeNeighbors(pattern,2))
    # text = "GTTAACGATGATGAACGATGATGGTTGTTAACGAACGGTTAACGATGCCCGTTATGATGCTGGCTGGATGCCCATGATGAACGCTGGCTGGAACGAACGGTTATGCCCATGAACGAACGCTGGCCCAACGAACGCTGGCTGGCCCCTGGAACGGTTAACGATGCTGGATGGTTGTTGTTCTGGCTGGATGCCCCTGGCCCCTGGAACGAACGATGATGCTGGAACGCCCATGGTTGTTGTTCCCATGATGCTGGAACGCCCCTGGATGAACGCTGGAACGATGGTTATGATGCTGGCTGGATGATGGTTCTGGCTGGAACGCCCAACGAACGCTGGAACGAACGCTGGCCCCCCGTTATG"
    # k = 5
    # d = 3
    # print(bio.frequentWordsWithMismatches(text,k,d))
    # text = "GAGAACGCAACGGACAACGCAGAACGGAGATGCAGCTGTGGCGCTGTGGAACGGAACGACGTGTGCAGACAGAACGTGGAGCCAACGCAGCTGCACAACGACGGACATGGCCATGCAACGCACAACGGAACGTGACGGAGACAGATGACGGAGACAACGGCGATGGAGATGCATGGCCAACGACGTGGCGCCAGCGAGCACGACGTGGAGAACGGATG"
    # k = 6
    # d = 3
    # text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
    # k = 4
    # d = 1
    # print(bio.frequentWordsWithMismatchesAndReverseComplements(text,k,d))
    print(bio.findSalmonellaOri())

if __name__ == "__main__":
    main()