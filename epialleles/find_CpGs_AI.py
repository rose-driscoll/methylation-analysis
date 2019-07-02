import sys

#run as "python find_CpGs_AI.py gapped_refseq_AI.txt > CpGs_AI.txt"

inFile = open(sys.argv[1],'r')

gapped_ref = ""
for line in inFile:
	gapped_ref = line

ungapped=0
AI130=0
AI133=0
AI135=0
AI143=0
AI150=0
AI163=0
AI181=0
AI187=0
AI228=0
AI252=0

for character in range(len(gapped_ref)):
	if gapped_ref[character] == "-":
		pass
	else:
		ungapped = ungapped + 1
		if ungapped == 130:
			AI130 = character+1
		if ungapped == 133:
			AI133 = character+1
		if ungapped == 135:
			AI135 = character+1
		if ungapped == 143:
			AI143 = character+1
		if ungapped == 150:
			AI150 = character+1
		if ungapped == 163:
			AI163 = character+1
		if ungapped == 181:
			AI181 = character+1
		if ungapped == 187:
			AI187 = character+1
		if ungapped == 228:
			AI228 = character+1
		if ungapped == 252:
			AI252 = character+1

print str(AI130)+","+str(AI133)+","+str(AI135)+","+str(AI143)+","+str(AI150)+","+str(AI163)+","+str(AI181)+","+str(AI187)+","+str(AI228)+","+str(AI252)
