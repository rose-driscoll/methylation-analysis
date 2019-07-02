import sys

#run as "python find_CpGs_AP.py gapped_refseq_AP.txt > CpGs_AP.txt"

inFile = open(sys.argv[1],'r')

gapped_ref = ""
for line in inFile:
	gapped_ref = line

ungapped=0
# AP117=0 This is the one that is inside the primer annealing site so I'm leaving it out
AP145=0
AP159=0
AP188=0
AP279=0
AP295=0
AP317=0
AP326=0

for character in range(len(gapped_ref)):
	if gapped_ref[character] == "-":
		pass
	else:
		ungapped = ungapped + 1
		if ungapped == 145:
			AP145 = character+1
		if ungapped == 159:
			AP159 = character+1
		if ungapped == 188:
			AP188 = character+1
		if ungapped == 279:
			AP279 = character+1
		if ungapped == 295:
			AP295 = character+1
		if ungapped == 317:
			AP317 = character+1
		if ungapped == 326:
			AP326 = character+1

print str(AP145)+","+str(AP159)+","+str(AP188)+","+str(AP279)+","+str(AP295)+","+str(AP317)+","+str(AP326)

