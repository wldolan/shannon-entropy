import sys, os, re
import matplotlib.pyplot as plt
import numpy as np
import math
from sets import Set

#AA frequencies from Proteins and Proteomics: A Laboratory Manual. RJ Simpson CSHL Press 2008
natFreq={'A':0.075,
         'R':0.052,
         'N':0.046,
         'D':0.052,
         'C':0.018,
         'E':0.041,
         'Q':0.063,
         'G':0.071,
         'H':0.022,
         'I':0.055,
         'L':0.091,
         'K':0.058,
         'M':0.028,
         'F':0.039,
         'P':0.051,
         'S':0.074,
         'T':0.060,
         'W':0.013,
         'Y':0.033,
         'V':0.065,
         'X':1.0,
         '-':1.0}

##readfile##
def readfile(filename):
    with open(filename) as file:
        seq=[(part[2].replace('\n',''))	
        for part in
            [entry.partition('\n')
        for entry in file.read().split('>')[1:]]]
    seqs=seq
    return (seqs)

	
#seqs= readfile(filename)
#print seqs

##get columns##
def getcolumns(seqs):
	n=len(seqs[0])
	#print n
	columns =[]
	alphabet =[]
	for j in range (n)  :
	   characters=[]
	   for i in range (len(seqs[0:])) :
		   characters.append(seqs[i][j])
	   columns.append(characters)
	   alphabet.append(list(set(columns[j])))
	return (columns, alphabet)


def freqs(alphabet,columns,natFreq):
    #print ('Alphabet of symbols in the string:')
    #print (alphabet)
    # calculate the frequency of each symbol in the column
    freqList = []
    adjFList = []
    for i in range(len(columns)):
        freq = []
        adjF = []
        for symbol in alphabet[i]:
            ctr = 0
            for sym in columns[i]:
                if sym == symbol:
                    ctr += 1
            charFreq=(float(ctr) / len(columns[i]))
            freq.append(charFreq)
            adjF.append(charFreq/natFreq[symbol])
        freqList.append(freq)
        adjFList.append(adjF)
        #print ('Frequencies of alphabet symbols:')
    #print freqList
    return (freqList,adjFList)


## Shannon entropy ##
def shannon(freqList, adjFList):
    ent = []
    for i in range(len(columns)):
        charEnt=[]
        for j in range(len(freqList[i])):
            #print freq
            charEnt.append((freqList[i][j] * (math.log(adjFList[i][j], 2))))
        ent.append(-sum(charEnt))
        
    #print ('Shannon entropy:')
    #print (ent)
    return ent 


#freqList=freqs(alphabet,columns)
#ent=shannon(freqList)

#plot the Shannon entropies
def plotent(ent):
	plt.plot(range(len(ent)), ent)
	plt.title('Position-specific Shannon Entropies')
	plt.xlabel('Position')
	plt.ylabel('Entropy')
	plt.show()

#### optional functions ####
	
#find motif coordinates
def findmotif(seqs, motif):
    p = re.compile(motif)
    matches=[]
    spans=[]
    for string in seqs:
        iterator = p.finditer(string)  
        matches.append(iterator)

        for match in iterator:
            #print match.span()
            spans.append(match.span())
            #print spans
    print (set(spans))

##plot a range##
def plotrange(begin, end):
	plt.plot(range(begin,end), ent[begin:end])
	plt.title('Position-specific Shannon Entropies')
	plt.xlabel('Position')
	plt.ylabel('Entropy')
	plt.show()

## aa frequencies for a range##
def freqrange(begin, end):
    for pos in range(begin,end):
        print "Position: ",pos, ent[pos]
        for char in range(len(freqList[pos])):
            print "\t{0}\t{1}%".format(alphabet[pos][char],100*(freqList[pos][char]))

## find positions below an entropy cutoff ##
def entcutoff(begin, end, cutoff):
    for pos in range(begin,end):
        if ent[pos] < cutoff:
            print ("Position: ",pos, ent[pos])
            for char in range(len(freqList[pos])):
                print ("\t{0}\t{1}%".format(alphabet[pos][char],100*(freqList[pos][char])))

## find positions above a frequency cutoff ##
def freqcutoff(begin, end, cutoff):
    for pos in range(begin, end):
        for char in range(len(freqList[pos])):
            if 100*freqList[pos][char] >= cutoff and alphabet[pos][char] is not "-": 
                print ("Position: ",pos, ent[pos])
                for char in range(len(freqList[pos])):
                    print ("\t{0}\t{1}%".format(alphabet[pos][char],100*(freqList[pos][char])))

def altFunc(command):
    if command == 'findmotif':
        motif = sys.argv[3]
        result = findmotif(seqs,motif)

    elif command == 'plotrange':
        begin = int(sys.argv[3])
        end = int(sys.argv[4])
        result = plotrange(begin, end)

    elif command == 'freqrange':
        begin = int(sys.argv[3])
        end = int(sys.argv[4])
        result = freqrange(begin, end)

    elif command == 'entcutoff':
        begin = int(sys.argv[3])
        end = int(sys.argv[4])
        cutoff = int(sys.argv[5])
        result = entcutoff(begin, end, cutoff)

    elif command == 'freqcutoff':
        begin = int(sys.argv[3])
        end = int(sys.argv[4])
        cutoff = int(sys.argv[5])
        result = freqcutoff(begin, end, cutoff)

    else: print ("commands: findmotif <motif>, plotrange <begin> <end>, freqrange <begin> <end>, entcutoff <begin> <end> <cutoff>, freqcutoff <begin> <end> <cutoff>")

	
#### MAIN ####
filename=sys.argv[1]
seqs=readfile(filename)
columns, alphabet= getcolumns(seqs)
freqList,adjFList=freqs(alphabet,columns,natFreq)
ent=shannon(freqList,adjFList)
if len(sys.argv)>2:
    command=sys.argv[2]
    altFunc(command)
else: plotent(ent)

