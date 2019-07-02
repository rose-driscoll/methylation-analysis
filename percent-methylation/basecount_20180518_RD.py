import sys

#run as "python basecount_20180518_RD.py sample.mpileup > sample.counts"

inFile = open(sys.argv[1],'r')

print 'sample\tamp_bp\tamplicon\tbp\tA\tG\tC\tT\tdel\tins\tinserted\tambiguous'
outtable = []
for line in inFile:
        sample = sys.argv[1][0:3]
        # I need to remember to start all of my filenames with the sample number!!!
        # That way, the first three characters of the filename indicate sample
        # This will be useful if I decide to concatenate all of my CpG files before importing into R, which I think would be a good idea
        # Also then I don't have to add a sample column in R        
        data = line.strip().split('\t')
        amplicon = data[0][0:2]
        bp = data[1]
        if len(data) < 5:
        	bases = ''
        else:
        	bases = data[4].upper()
        ref = data[2].upper()
        
        types = {'A':0,'G':0,'C':0,'T':0,'-':0,'+':[],'X':[]}

        i = 0
        while i < len(bases):
                base = bases[i]
                if base == '^' or base == '$':
                        i += 1
                elif base == '-':
                        i += 1
                elif base == '*':
                        types['-'] += 1
                elif base == '+':
                        i += 1
                        addNum = int(bases[i])
                        addSeq = ''
                        for a in range(addNum):
                                i += 1
                                addSeq += bases[i]

                        types['+'].append(addSeq)
                elif base == '.' or base == ',':
                        types[ref] += 1
                else:
                        if types.has_key(base):
                                types[base] += 1
                        else:
                                types['X'].append(base)

                i += 1

        adds = '.'
        if len(types['+']) > 0:
                adds = ','.join(types['+'])

        amb = '.'
        if len(types['X']) > 0:
                amb = ','.join(types['X'])

        out = [sample,amplicon+bp,amplicon,bp,types['A'],types['G'],types['C'],types['T'],types['-'],len(types['+']),adds,amb]
        print '\t'.join([str(x) for x in out])