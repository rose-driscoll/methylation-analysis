import sys

#run as "python find_CpGs_B2.py gapped_refseq_B2.txt > CpGs_B2.txt"

inFile = open(sys.argv[1],'r')

gapped_ref = ""
for line in inFile:
	gapped_ref = line

ungapped=0
B2174=0
B2187=0
B2213=0
B2238=0
B2253=0
B2302=0
B2314=0
for character in range(len(gapped_ref)):
	if gapped_ref[character] == "-":
		pass
	else:
		ungapped = ungapped + 1
		if ungapped == 174:
			B2174 = character+1
		if ungapped == 187:
			B2187 = character+1
		if ungapped == 213:
			B2213 = character+1
		if ungapped == 238:
			B2238 = character+1
		if ungapped == 253:
			B2253 = character+1
		if ungapped == 302:
			B2302 = character+1
		if ungapped == 314:
			B2314 = character+1


print str(B2174)+","+str(B2187)+","+str(B2213)+","+str(B2238)+","+str(B2253)+","+str(B2302)+","+str(B2314)