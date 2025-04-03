import gffutils as gff
import os
import re
from pathlib import Path

#https://www.sciencedaily.com/releases/2020/09/200925113447.htm
# to do:
# clean-up
# better scoring
# start.bat
# use data
# AT content
# inverted repeats
# ambiguities
# multiple nonas
# use program to find cruci iterons and make database
# 1stem-loop inverted repeats start right after nona end and a little before nona start, always longer than 5?

#fix print messages, capitalization w list method

knownIterons = ['TTGTCCAC','AGTGGGA', 'GCCACCC', 'GGGGA', 'TCTGA', 
                'TTGAGAA', 'GGCCGGG', 'GTCCCG', 'TATCTCGCT', 'GTTACAT',
                'AGGCGCA', 'CTGGTCTAAA', 'CACACACACACA', 'AGCGTT', 'GGGTTAGG', 
                'TGGAA', 'AAAT', 'TGTCTC', 'TAAAC', 'AGTGTCC', 
                'GGCTC', 'GATCAG', 'TGTAAC', 'CGTATC', 'GTTTA', 
                'GAAA', 'TACCCATT', 'AAAC', 'AATT', 'TAAAAAA', 
                'AAAT', 'GAGAAT', 'TCCACGA', 'TAACC', 'ATTTTTT', 
                'GCCG', 'AGTGTC', 'GATAGAAG', 'CCTT', 'TCTCAC', 
                'AAAAT', 'GGAACA', 'TCGCTCCC', 'CTTTTATA', 'TTGCAC', 
                'ATGAGCC', 'GGGACAC', 'AGGGACA', 'GACCC', 'GTGTCCT',
                'GGGACACA', 'GTGTCA', 'GTCCCT', 'TCTGA']

knownIteronsNoAT = ['TTGTCCAC','AGTGGGA', 'GCCACCC', 'GGGGA', 'TCTGA', 
                    'TTGAGAA', 'GGCCGGG', 'GTCCCG', 'TATCTCGCT', 'GTTACAT or GTTAAAT',
                    'AGGCGCA', 'CTGGTCTAAA', 'CACACACACACA', 'AGGGTT or AGCGTT', 'GGGTTAGG', 
                    'TGGAA', 'TGTCTC', 'AGTGTCC', 'GGCTC', 'GATCAG', 'TGTAAC', 'CGTATC', 'GTTTA', 
                    'GAAA', 'TACCCATT', 'GAGAAT', 'TCCACGA', 'TAACC', 
                    'GCCG', 'AGTGTC', 'GATAGAAG', 'CCTT', 'TCTCAC', 'GGAACA', 'TCGCTCCC', 'TTGCAC', 
                    'ATGAGCC', 'GGGACAC', 'AGGGACA', 'GACCC', 'GTGTCCT',
                    'GGGACACA', 'GTGTCA', 'GTCCCT']

# Update the IUPAC dictionary to handle both cases
iupac_to_nuc = {
    'a': 'a', 't': 'tu', 'g': 'g', 'c': 'c',
    'r': 'ag', 'y': 'ctu', 's': 'gc', 'w': 'atu',
    'k': 'gtu', 'm': 'ac', 'b': 'cgtu', 'd': 'agtu',
    'h': 'actu', 'v': 'acg', 'n': 'acgtu', 'u': 'ut',
    'A': 'a', 'T': 'tu', 'G': 'g', 'C': 'c',
    'R': 'ag', 'Y': 'ctu', 'S': 'gc', 'W': 'atu',
    'K': 'gtu', 'M': 'ac', 'B': 'cgtu', 'D': 'agtu',
    'H': 'actu', 'V': 'acg', 'N': 'acgtu', 'U': 'ut'
}

# Update complement dictionary to handle both cases
complementDNA = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
    'R': 'Y', 'Y': 'R', 'M': 'K', 'K': 'M',
    'r': 'y', 'y': 'r', 'm': 'k', 'k': 'm',
    'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D',
    'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd',
    'N': 'N', 'n': 'n', 'W': 'W', 'w': 'w'
}

tagIteronsDict =    {'TGGTCA': 'RDHV', 'CGGCAG': 'PCV1/PCV2 Iteron',
                    'GGGGCACC': 'BFDV1/2 Iteron', 'GGTGTCTGGAGTC': 'ToLCV-NDE A1 Iteron',
                    'GGTGTC': 'ToLCV-NDE A1 5` Iteron', 'GGCGTCTGGCGTC': 'ToLCV-NDE A2 Iteron',
                    'GGAGCCAC': 'StCV Iteron', 'GGAACCAC': 'FiCV Iteron', 'GGAGCCAC': 'RACV/CaCV/CoCV Iteron',
                    'GGGGCCAT': 'GuCV Iteron', 'GTACTCC': 'DuCV Iteron', 'GTACTCC': 'GoCV Iteron',
                    'GGCGT': 'Tomato Leaf Curl New Delhi Cucumber Iteron',
                    'GGTGT': 'Tomato Leaf Curl New Delhi Iteron',
                    'GGAGT': 'Tomato Motile VIrus Iteron',
                    'GGTGTC': 'Soybean Crinkle Leaf Virus Iteron',
                    'TATTAC': 'RDHV Minimal Nick Site'}

tagIteronsDict = {}



def createDB(gffFile):
    db = gff.create_db(gffFile, ":memory:", merge_strategy = "create_unique", keep_order = True)
    listOfNonas = list(db.features_of_type("nonanucleotide"))
    if len(listOfNonas) < 1:
        listOfNonas = []
        for line in open(gffFile, 'r'):
            if "Nona" in line:
                listOfNonas.append(gff.feature.feature_from_line(line))
    if len(listOfNonas) < 1:
        #print('No nonanucleotides to search around!')
        return None
    nnlist = listOfNonas
    listOfLoops = list(db.features_of_type("stem_loop"))
    if len(listOfLoops) < 1:
        listOfLoops = []
        for line in open(gffFile, 'r'):
            if "tem-loop" in line:
                listOfLoops.append(gff.feature.feature_from_line(line))
    if len(listOfLoops) < 1:
        #print('No stem-loops to search around!')
        return None
    stem_looplist = listOfLoops
    regionlength = list(db.features_of_type("region"))[0].end
    return nnlist, regionlength, stem_looplist

def getOutputFile(inputFile, outputDir):
    """Get output file path, handling both string and Path objects correctly"""
    # Convert to Path objects if they aren't already
    input_path = Path(inputFile)
    output_dir = Path(outputDir)
    
    # Get the base name without .gff extension
    input_filename = input_path.stem
    
    # Create output path
    output_file = output_dir / f"{input_filename}.gff"
    return str(output_file)  # Return as string for compatibility


class Iteron(object):
    def __init__(self, sequence, positions = None, dist = None, aveDist = None, closeDist = None, distFromNona = None, score = 1000, stemLoop = None, knownIteron = None, perfect = None, tag = False):
        if positions is None:
            self.positions = []
        else:
            self.positions = positions
        if dist is None:
            self.dist = []
        else:
            self.dist = dist
        if distFromNona is None:
            self.distFromNona = []
        else:
            self.distFromNona = distFromNona
        if stemLoop is None:
            self.stemLoop = []
        else:
            self.stemLoop = stemLoop
        self.dist = dist
        self.aveDist = aveDist
        self.closeDist = closeDist
        self.sequence = sequence
        self.score = score
        self.stemLoop = stemLoop
        self.knownIteron = knownIteron
        self.perfect = perfect
        self.tag = tag

    def fixPositions(self, buffer, nn, regionlength):
        for x in range(len(self.positions)):
            if nn.start-buffer < 0:
                overlap = buffer-nn.start
            # elif nn.end+buffer > regionlength:
            #     overlap = regionlength-nn.start+buffer
            else:
                overlap = False
            if not overlap == False:
                if self.positions[x] < overlap:
                    self.positions[x] = nn.start - buffer + self.positions[x]
                else:
                    self.positions[x] = nn.start - buffer + self.positions[x]
            else:
                self.positions[x] = nn.start - buffer + self.positions[x]
            #if nn.start < self.positions[x] < nn.end:
                #print('tagged iteron?')
                #print(self.sequence)

    def getDistance(self):
        distances = []
        for y in range(len(self.positions)-1):
            distances.append(self.positions[y+1]-self.positions[y]-len(self.sequence))
        self.dist = distances
        
    def averageDistance(self):
        '''Adds the average distance between all repeats to the info list at index 2'''
        indexsum = 0
        for x in self.positions:
            indexsum += x-self.positions[0]
        if len(self.positions) > 1:
            indexave = indexsum/(len(self.positions)-1) - len(self.sequence)
            self.aveDist = round(indexave, 2)
        else:
            self.aveDist = 0
        
    def closestDistance(self):
        '''Adds the closest distance between any two repeats to the info list at index 3'''
        positions = sorted(self.dist)
        self.closeDist = positions[0]

    def distanceFromNona(self, nn):
        '''Uses the position of each repeat to calculate the distance from the nona'''
        for y in range(len(self.positions)):
            if self.positions[y] < nn.start:
                if self.positions[y] - nn.start + len(self.sequence) <= 0:
                    self.distFromNona.append(self.positions[y] - nn.start + len(self.sequence))
                else:
                    self.distFromNona.append(0)
            elif self.positions[y] > nn.end:
                self.distFromNona.append(self.positions[y] - nn.end - 1)
            else:
                self.distFromNona.append(0)

    def checkStemLoop(self, stem_loop, nn, regionlength):
        start = 0
        end = 0
        self.stemLoop = []
        if nn.start > stem_loop.start:
            for i in self.positions:
                if nn.start > i >= stem_loop.start:
                    self.stemLoop.append('Start')
                elif nn.end < i + len(self.sequence) - 1 <= stem_loop.end:
                    self.stemLoop.append('End')
                else:
                    self.stemLoop.append(None)
        elif nn.start < stem_loop.start:
            for i in self.positions:
                if nn.start + regionlength > i + regionlength >= stem_loop.start:
                    self.stemLoop.append('Start')
                elif nn.end + regionlength < i + len(self.sequence) + regionlength - 1 <= stem_loop.end:
                    self.stemLoop.append('End')
                else:
                    self.stemLoop.append(None)
        for x in self.stemLoop:
            if x == 'Start':
                start += 1
            elif x == 'End':
                end += 1
        if start < 1 or end < 1:
            self.stemLoop = None
        
    def calcScore(self, wiggle, buffer, goodLength, bestDist):
        """Evaluates the iteron candidates based on distance between the iterons,
        distance from the nonanucleotide, whether the iteron is all A's and T's,
        whether the entire candidate is one nucleotide, and candidate length. 
        Can compare to a database of known iterons"""
        score = 0
        distcheck = False
        seconddistcheck = False
        AT = ['A', 'T']
        ATcount = 0
        sameBaseCount = 0
        # dist between iterons
        for x in self.dist:
            if abs(len(self.sequence) - x) <= wiggle:
                distcheck = True
            if x <= bestDist:
                seconddistcheck = True
        if not distcheck:
            score += self.closeDist
        if not seconddistcheck:
            score += (self.closeDist-bestDist)*3
        
        # dist from nona check
        distances = sorted([abs(x) for x in self.distFromNona])
        score += distances[0]
            
        # AT content check
        for i in range(len(self.sequence)):
            if self.sequence[i] in AT:
                ATcount += 1
        if ATcount == len(self.sequence):
            score += 30
        elif ATcount == len(self.sequence) - 1:
            score += 15

        # same base check
        firstBase, secondBase = self.sequence[0], self.sequence[1]
        if firstBase == secondBase:
            # baseball?
            for i in range(len(self.sequence)):
                if self.sequence[i] == firstBase:
                    sameBaseCount += 1
        elif firstBase == self.sequence[2] or secondBase == self.sequence[2]:
            for i in range(len(self.sequence)):
                if self.sequence[i] == self.sequence[2]:
                    sameBaseCount += 1
        if sameBaseCount == len(self.sequence):
            score += 100
        elif sameBaseCount == len(self.sequence) - 1:
            score += 50

        # length check
        if len(self.sequence) > goodLength:
            score += (len(self.sequence) - goodLength)*1.5
        if len(self.sequence) <= 5:
            score += 50
        # 3' check (later)

        # known check
        if not self.knownIteron == None:
            score -= 10
        
        self.score = score

    def secondCheck(self, samplestring, buffer):
        """Uses current candidate sequences to find potential imperfect third+ repeats"""
        isPerfect = []
        listOfPositions = self.positions
        i=0
        while i < len(self.positions):
            isPerfect.append(True)
            i+=1
        if self.sequence in tagIteronsDict:
            for n in range(len(samplestring)-len(self.sequence)):
                thisRepeat = samplestring[n:n+len(self.sequence)]
                if getInvertedRepeat(thisRepeat) == self.sequence:
                    check = True
                    for x in listOfPositions:
                        if n > x + len(self.sequence) or n < x - len(self.sequence):
                            pass
                        else:
                            check = False
                    if check:
                        for y in range(len(listOfPositions)):
                            if n < listOfPositions[y]:
                                index = y
                                break
                            else:
                                index = -1
                        isPerfect.insert(index, True)
                        listOfPositions.insert(index, n)
        if len(self.sequence) > 5 and len(self.positions) > 2:
            listOfPositions = self.positions
            for n in range(len(samplestring)-len(self.sequence)):
                thisRepeat = samplestring[n:n+len(self.sequence)]
                count = 0
                for y in range(len(thisRepeat)):
                    if not thisRepeat[y] == self.sequence[y]:
                        count += 1
                if count == 1:
                    if n+len(self.sequence) <= buffer or n > len(samplestring)-buffer:
                        check = True
                        for x in listOfPositions:
                            if n > x + len(self.sequence) or n < x - len(self.sequence):
                                pass
                            else:
                                check = False
                        if check:
                            for y in range(len(listOfPositions)):
                                if n < listOfPositions[y]:
                                    index = y
                                    break
                                else:
                                    index = -1
                            isPerfect.insert(index, False)
                            listOfPositions.insert(index, n)
        self.perfect = isPerfect
        self.positions = listOfPositions

    # def checkKnown(self):
    #     """Checks a database of iterons (warning: slows program significantly)"""
    #     for i in range(len(self.sequence)):
    #         regex = self.sequence 
    #         for x in knownIterons:
    #             if not re.search(regex, x) == None:
    #                 self.knownIteron = re.search(regex, x)
    #             elif not re.search(regex, getInvertedRepeat(x)) == None:
    #                 self.knownIteron = re.search(regex, getInvertedRepeat(x))
    #             elif not re.search(regex, getReverseRepeat(x)) == None:
    #                 self.knownIteron = re.search(regex, getReverseRepeat(x))

    #         l = list(self.sequence)
    #         l.insert(i, '[ATCG]')
    #         regex = ''.join(l)
    #         for x in knownIterons:
    #             if not re.search(regex, x) == None:
    #                 self.knownIteron = re.search(regex, x)
    #             elif not re.search(regex, getInvertedRepeat(x)) == None:
    #                 self.knownIteron = re.search(regex, getInvertedRepeat(x))
    #             elif not re.search(regex, getReverseRepeat(x)) == None:
    #                 self.knownIteron = re.search(regex, getReverseRepeat(x))

    #         l = list(self.sequence)
    #         del l[i]
    #         regex = ''.join(l)
    #         for x in knownIterons:
    #             if not re.search(regex, x) == None:
    #                 self.knownIteron = re.search(regex, x)
    #             elif not re.search(regex, getInvertedRepeat(x)) == None:
    #                 self.knownIteron = re.search(regex, getInvertedRepeat(x))
    #             elif not re.search(regex, getReverseRepeat(x)) == None:
    #                 self.knownIteron = re.search(regex, getReverseRepeat(x))

    #         l = list(self.sequence)
    #         l[i] = '[ATCG]'
    #         regex = ''.join(l)
    #         for x in knownIterons:
    #             if not re.search(regex, x) == None:
    #                 self.knownIteron = re.search(regex, x)
    #             elif not re.search(regex, getInvertedRepeat(x)) == None:
    #                 self.knownIteron = re.search(regex, getInvertedRepeat(x))
    #             elif not re.search(regex, getReverseRepeat(x)) == None:
    #                 self.knownIteron = re.search(regex, getReverseRepeat(x))






def getInvertedRepeat(string):
    return ''.join([complementDNA[char] for char in string[::-1]])

def getReverseRepeat(string):
    return ''.join([char for char in string[::-1]])

def getSampleString(gffFile, nn, buffer, regionlength):
    """Get sequence from GFF file, preserving case"""
    samplestringprec = []
    start = False
    gffFile = open(gffFile, 'r')
    for line in gffFile:
        if ">" in line:
            start = True
        elif start:
            line = line.strip('\n')
            samplestringprec.append(line)
    samplestring = ''.join(samplestringprec)
    
    if nn.end + buffer > len(samplestring):
        samplestringprec.append(samplestring[:buffer - regionlength + nn.end])
        samplestring = ''.join(samplestringprec)
    elif nn.start - buffer < 0:
        splitstringlist = []
        splitstringlist.append(samplestring[nn.start-buffer-1:len(samplestring)])
        splitstringlist.append(samplestring[0:nn.end+buffer])
        return ''.join(splitstringlist)
    
    # Return without forcing case
    return samplestring[nn.start-buffer-1:nn.end+buffer]

def buildStringDict(samplestring, buffer, minLength, maxLength):
    '''Builds a dictionary of every substring inside the sequence string between the 
    minimum and maximum length as keys, storing number of times repeated and locations
    in an info list at indices 0 and 1'''
    samplestringdict = {}
    for i in range(len(samplestring)-minLength+1):
        # i is index, n is distance from index, a is built string
        a = samplestring[i]
        n = 1
        while 1 <= len(a) < minLength:
            a += samplestring[i+n]
            n += 1
        while minLength <= len(a) <= maxLength:
            check = False
            if a in samplestringdict:
                if i - len(a) >= samplestringdict[a][1][-1]:
                    samplestringdict[a][0] += 1
                    samplestringdict[a][1].append(i)
                else:
                    samplestringdict[a][1].append(i)
                check = True
            # else:
            #     if getInvertedRepeat(a) in samplestringdict:
            #         if i+1 - len(a) >= samplestringdict[getInvertedRepeat(a)][1][-1]:
            #             samplestringdict[getInvertedRepeat(a)][0] += 1
            #             samplestringdict[getInvertedRepeat(a)][1].append(i+1)
            #         else:
            #             samplestringdict[a][1].append(i+1)
            #         check = True
            #     if getReverseRepeat(a) in samplestringdict:
            #         if i+1 - len(a) >= samplestringdict[getReverseRepeat(a)][1][-1]:
            #             samplestringdict[getReverseRepeat(a)][0] += 1
            #             samplestringdict[getReverseRepeat(a)][1].append(i+1)
            #         else:
            #             samplestringdict[a][1].append(i+1)
            #         check = True
            if not check:
                if a in tagIteronsDict:
                    samplestringdict[a] = [1001, [i]]
                else:
                    samplestringdict[a] = [1, [i]]
            # don't overflow sequence
            if i+n >= len(samplestring):
                break
            a += samplestring[i+n]
            n += 1
    return samplestringdict

def removeJunk(samplestringdict, samplestring, buffer, maxDist):
    '''Removes all substrings that only repeat once, then removes all substrings 
    that are entirely contained inside another substring. Then, does preliminary distance
    and location analysis on candidates'''
    # remove repeats (lower count than length of positions)
    repeatdict = {x: samplestringdict[x] for x in samplestringdict if samplestringdict[x][0] > 1}
    repeatdictkeys = list(repeatdict.keys())
    for possub in repeatdict:
        for possuper in repeatdictkeys:
            # check if:
            #   1) the possible superstring is longer than the possible substring,
            #   2) the substring is contained in the superstring,
            #   3) both strings appear the same number of times, and
            #   4) the possible substring has not already been removed
            if     len(possuper) > len(possub) \
               and possub in possuper \
               and (repeatdict[possuper][0] == repeatdict[possub][0] or repeatdict[possuper][0]+1000 == repeatdict[possub][0]) \
               and possub in repeatdictkeys:
                repeatdictkeys.remove(possub)
    
    finalrepeatdict = {x: repeatdict[x] for x in repeatdictkeys}

    for x in finalrepeatdict:
        if x not in tagIteronsDict:
            positions = finalrepeatdict[x][1][:]
            y = 0
            while y < len(positions):

                if len(positions) <= 1:
                    break

                if y == 0:
                    thisDistance = positions[y+1]-positions[y]-len(x)
                    if thisDistance <= maxDist:
                        y += 1
                    else:
                        positions.pop(0)
                            
                elif y == len(positions)-1:
                    previousDistance = positions[y]-positions[y-1]-len(x)
                    if 0 <= previousDistance <= maxDist:
                        y += 1
                    else:
                        positions.pop(y)
                        y -= 1

                else:
                    thisDistance = positions[y+1]-positions[y]-len(x)
                    previousDistance = positions[y]-positions[y-1]-len(x)
                    if 0 <= thisDistance <= maxDist or 0 <= previousDistance <= maxDist:
                        y += 1
                    elif thisDistance < 0:
                        if previousDistance < 0:
                            positions.pop(y)
                        else:
                            y += 1
                    elif previousDistance < 0:
                        positions.pop(y)
                        y -= 1
                    else:
                        positions.pop(y)

            finalrepeatdict[x][1] = positions
        else:
            positions = finalrepeatdict[x][1][:]
            y = 0
            while y < len(positions):

                if len(positions) <= 1:
                    break

                if y == 0:
                    thisDistance = positions[y+1]-positions[y]-len(x)
                    y += 1
                            
                elif y == len(positions)-1:
                    previousDistance = positions[y]-positions[y-1]-len(x)
                    if 0 <= previousDistance:
                        y += 1
                    else:
                        positions.pop(y)
                        y -= 1

                else:
                    thisDistance = positions[y+1]-positions[y]-len(x)
                    previousDistance = positions[y]-positions[y-1]-len(x)
                    if 0 <= thisDistance or 0 <= previousDistance:
                        y += 1
                    elif thisDistance < 0:
                        if previousDistance < 0:
                            positions.pop(y)
                        else:
                            y += 1
                    elif previousDistance < 0:
                        positions.pop(y)
                        y -= 1
                    else:
                        positions.pop(y)

            finalrepeatdict[x][1] = positions
    finalrepeatdict = {x:finalrepeatdict[x] for x in finalrepeatdict if len(finalrepeatdict[x][1]) > 1 or finalrepeatdict[x][0] > 1000}

    for x in finalrepeatdict:
        if x not in tagIteronsDict:
            startcount = 0
            endcount = 0
            for y in finalrepeatdict[x][1]:
                if y+len(x) <= buffer:
                    startcount += 1
                elif y > len(samplestring) - buffer:
                    endcount += 1
            midcount = len(finalrepeatdict[x][1]) - startcount - endcount
            if startcount < 2 and endcount < 2:
                finalrepeatdict[x][0] = 0
            elif startcount < 2:
                finalrepeatdict[x][0] = endcount
                w = 0
                while w < startcount + midcount:
                    finalrepeatdict[x][1].pop(0)
                    w += 1
            elif endcount < 2:
                finalrepeatdict[x][0] = startcount
                w = 0
                while w < endcount + midcount:
                    finalrepeatdict[x][1].pop(-1)
                    w += 1
            else:
                finalrepeatdict[x][0] = startcount + endcount
                w = 0
                while w < midcount:
                    finalrepeatdict[x][1].pop(startcount)
                    w += 1
    repeatdict = {x: finalrepeatdict[x] for x in finalrepeatdict if finalrepeatdict[x][0] > 1}
    repeatdictkeys = list(repeatdict.keys())
    for possub in repeatdict:
        for possuper in repeatdictkeys:
            # check if:
            #   1) the possible superstring is longer than the possible substring,
            #   2) the substring is contained in the superstring,
            #   3) both strings appear the same number of times, and
            #   4) the possible substring has not already been removed
            if     len(possuper) > len(possub) \
               and possub in possuper \
               and (repeatdict[possuper][0] == repeatdict[possub][0] or repeatdict[possuper][0]+1000 == repeatdict[possub][0]) \
               and possub in repeatdictkeys:
                repeatdictkeys.remove(possub)
    finalrepeatdict = {x: repeatdict[x] for x in repeatdictkeys}
    for x in finalrepeatdict:
        if x in tagIteronsDict:
            finalrepeatdict[x].append(tagIteronsDict[x])
        else:
            finalrepeatdict[x].append(False)
    return finalrepeatdict

def fetchIteron(iteronslist, sequence):
    for x in iteronslist:
        if x.sequence == sequence:
            return x
    return None

def outputInfo(iteronslist, inputFile, outputFile, regionlength, doStemLoop, doKnownIterons, notFirstNona):
    '''Takes iteron candidates and outputs it into a copy of the original GFF
    file, for import into Geneious. Distinguishes between iterons 
    and stem-loop repeats'''
    newiteronfeatures = []
    newFile = []
    start = True
    if len(iteronslist) > 0:
        if notFirstNona:
            inputFile = open(outputFile, 'r')
            db = gff.create_db(outputFile, ":memory:", merge_strategy = "create_unique", keep_order = True)
            listOfIterons = list(db.features_of_type("iteron"))
            for iteron in iteronslist:
                for x in range(len(iteron.positions)):
                    for it in listOfIterons:
                        posList = [iteron.positions[x], iteron.positions[x]+regionlength, iteron.positions[x]-regionlength]
                        if it.start in posList and it.start+len(iteron.sequence)-1 == it.end:
                            iteron.positions[x] = 'REMOVE'
                            break
                iteron.positions = [x for x in iteron.positions if not x == 'REMOVE']
        else:
            inputFile = open(inputFile, 'r')
        a = inputFile.name
        a = a.split('/')[-1]
        a = a.strip('.gff')
        if notFirstNona:
            a = a.strip(' - copy')
        inputFileName = a
        identifier = 1
        sIdentifier = 1
        for x in iteronslist:
            if not x.stemLoop == None and doStemLoop == True:
                for y in range(len(x.positions)):
                    perfectTag = ""
                    if not x.perfect[y]:
                        perfectTag = "$"
                    if x.positions[y] <= 0:
                        start = x.positions[y] + regionlength
                        end = start + len(x.sequence)
                    else:
                        start = x.positions[y]
                        end = start + len(x.sequence)
                    if x.sequence in tagIteronsDict:
                        input("hey check this out")
                    elif not x.stemLoop[y] == None:
                        newiteronfeatures.append(str(inputFileName) + "\tCRUISE\tstem_loop_repeats" + "\t" + str(start) + "\t" + str(end-1) + "\t.\t.\t.\tName=" + perfectTag + "stem-loop repeats" + str(sIdentifier) + '\n')
                    else:
                        newiteronfeatures.append(str(inputFileName) + "\tCRUISE\tstem_loop_repeats" + "\t" + str(start) + "\t" + str(end-1) + "\t.\t.\t.\tName=" + perfectTag + "S iteron" + str(sIdentifier) + '\n')
                sIdentifier += 1
            # elif not x.knownIteron == None and x.stemLoop == None and doKnownIterons == True:
            #     for y in x.positions:
            #         if y <= 0:
            #             start = y + regionlength
            #             end = start + len(x.sequence)
            #         else:
            #             start = y
            #             end = y + len(x.sequence)
            #         newiteronfeatures.append(str(inputFileName) + "\tCRUISE\tknown_iteron" + str(identifier) + "\t" + str(start) + "\t" + str(end-1) + "\t.\t.\t.\tName=K iteron" + '\n')
            #     identifier += 1
            else:
                for i in range(len(x.positions)):
                    y = x.positions[i]
                    perfectTag = ""
                    if not x.perfect[i]:
                        perfectTag = "$"
                    if y <= 0:
                        start = y + regionlength
                        end = start + len(x.sequence)
                    else:
                        start = y
                        end = y + len(x.sequence)
                    if x.sequence in tagIteronsDict:
                        newiteronfeatures.append(str(inputFileName) + "\tCRUISE\t" + 'tagIteron' + "\t" + str(start) + "\t" + str(end-1) + "\t.\t.\t.\tName=" + perfectTag + tagIteronsDict[x.sequence] + str(identifier) + '\n')
                    else:
                        newiteronfeatures.append(str(inputFileName) + "\tCRUISE\titeron" + "\t" + str(start) + "\t" + str(end-1) + "\t.\t.\t.\tName=" + perfectTag + "iteron" + str(identifier) + '\n')
                identifier += 1
        newiteronfeatures = ''.join(newiteronfeatures)
        for line in inputFile:
            skip = False
            if not skip:
                newFile.append(line)
            if '##sequence' in line and start:
                newFile.append(newiteronfeatures)
                start = False
        newFileS = ''.join(newFile)
        inputFile.close()
        outputFile = open(outputFile, 'w')
        outputFile.write(newFileS)
        return True
    #input("no iterons")
    return False

def findPosIterons(gffFileIn, gffFileOut, minLength, maxLength, buffer, wiggle, rank, numberTopIterons, maxScore, goodLength, doStemLoop, doKnownIterons, maxDist, bestDist, scoreRange):
    '''Takes a gffFile and returns a copy with annotated iterons'''
    if not createDB(gffFileIn) == None:
        Total, RealTotal, foundIterons, newIteronCount = 0, 0, 0, 0
        add = 0
        nnlist, regionlength, stem_looplist = createDB(gffFileIn)
        notFirstNona = False
        for nn in nnlist:
            stem_loop = None
            for posStemLoop in stem_looplist:
                if posStemLoop.start < nn.start and posStemLoop.end > nn.end:
                    stem_loop = posStemLoop
                    break
            if stem_loop == None:
                for posStemLoop in stem_looplist:
                    if posStemLoop.start < regionlength+nn.start and posStemLoop.end > regionlength+nn.end:
                        stem_loop = posStemLoop
                        break
            if stem_loop == None:
                for posStemLoop in stem_looplist:
                    if posStemLoop.start+regionlength < nn.start and posStemLoop.end+regionlength > nn.end:
                        stem_loop = posStemLoop
                        break
            if not stem_loop == None:
                finaliteronslist = []
                sequence = getSampleString(gffFileIn, nn, buffer, regionlength)
                # print('Nona: ' + sequence[nn.start - 1:nn.end])
                sequencedict = removeJunk(buildStringDict(sequence, buffer, minLength, maxLength), sequence, buffer, maxDist)
                iteronslist  = [Iteron(sequence = x, positions = sequencedict[x][1], tag = sequencedict[x][2]) for x in sequencedict]

                for x in iteronslist:
                    if not x.tag == False:
                        x.secondCheck(sequence, buffer)
                        x.fixPositions(buffer, nn, regionlength)
                        x.score = 1000
                    else:
                        x.secondCheck(sequence, buffer)
                        x.fixPositions(buffer, nn, regionlength)
                        x.getDistance()
                        x.averageDistance()
                        x.closestDistance()
                        x.distanceFromNona(nn)
                        # x.checkKnown()
                        x.calcScore(wiggle, buffer, goodLength, bestDist)
                        # print(x.sequence, x.positions, x.dist, x.distFromNona, x.score)
                        x.checkStemLoop(stem_loop, nn, regionlength)
                        if not x.stemLoop == None:
                            x.score = 0
                            add += 1
                        if rank:
                            pass
                        else:
                            Total += 1
                            if x.score <= maxScore:
                                finaliteronslist.append(x)
                                RealTotal += 1
                    #print(x.sequence,x.score)
                if rank:
                    iteronslist.sort(key=lambda x:x.score)
                    finaliteronslist = iteronslist[:numberTopIterons+add]
                    if len(finaliteronslist) > 0:
                        finaliteronslist = [x for x in finaliteronslist if x.score < finaliteronslist[0].score + scoreRange]
                # else:
                #     # print(len(finaliteronslist))
                for x in iteronslist:
                    if x.tag:
                        finaliteronslist.append(x)
                #for x in finaliteronslist:
                    #print(x.sequence, x.positions, x.stemLoop)
                anyIterons = outputInfo(finaliteronslist, gffFileIn, gffFileOut, regionlength, doStemLoop, doKnownIterons, notFirstNona)
                if anyIterons:
                    foundIterons = True
                    notFirstNona = True
                    newIteronCount = len(finaliteronslist)
            #else:
                #print('No stem-loop?')
            #return finaliteronslist
        return (foundIterons, newIteronCount)