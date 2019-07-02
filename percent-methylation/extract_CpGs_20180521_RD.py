import sys

#run as "python extract_CpGs_20180521_RD.py basecount_file > CpGs_file"

basecountFile = open(sys.argv[1],'r')

#print 'sample\tamp_bp\tamplicon\tbp\tCpG\tposition\tA\tG\tC\tT\tdel\tins\tinserted\tambiguous\tprop_meth\tpercent_meth'

for line in basecountFile:
        data = line.strip().split('\t')
        if data[1] in ['AI130', 'AI133', 'AI135', 'AI143', 'AI150', 'AI163', 'AI181', 'AI187', 'AI228', 'AI252', 'AP117', 'AP145', 'AP159', 'AP188', 'AP279', 'AP295', 'AP317', 'AP326', 'B2174', 'B2187', 'B2213', 'B2238', 'B2253', 'B2302', 'B2314', 'BI217', 'BI333', 'BI358', 'BI370']:
                data[3] = int(data[3])
                data[4] = int(data[4])
                data[5] = int(data[5])
                data[6] = int(data[6])
                data[7] = int(data[7])
                data[8] = int(data[8])
                data[9] = int(data[9])
                prop_meth = float(data[6])/(float(data[6])+float(data[7]))
                percent_meth = prop_meth*100
                data = data + [prop_meth, percent_meth]
                if data[2] == "AI":
                	position = data[3] - 1472
                elif data[2] == "AP":
                	position = data[3] - 391
                elif data[2] == "B2":
                	position = data[3] - 246
                elif data[2] == "BI":
                	position = -data[3] - 1079
                CpG = str(position)+data[2][0]
                data = data[:4] + [CpG] + [position] + data[4:]
                print '\t'.join([str(x) for x in data])