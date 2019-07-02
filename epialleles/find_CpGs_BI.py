import sys

#run as "python find_CpGs_BI.py gapped_refseq_BI.txt > CpGs_BI.txt"

inFile = open(sys.argv[1],'r')

gapped_ref = ""
for line in inFile:
	gapped_ref = line

ungapped=0
BI217=0
BI333=0
BI358=0
BI370=0

for character in range(len(gapped_ref)):
	if gapped_ref[character] == "-":
		pass
	else:
		ungapped = ungapped + 1
		if ungapped == 217:
			BI217 = character+1
		if ungapped == 333:
			BI333 = character+1
		if ungapped == 358:
			BI358 = character+1
		if ungapped == 370:
			BI370 = character+1

print str(BI217)+","+str(BI333)+","+str(BI358)+","+str(BI370)