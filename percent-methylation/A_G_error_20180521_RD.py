import sys

#run as "python A_G_error_20180521_RD.py basecount_file > A_G_error_file"

basecountFile = open(sys.argv[1],'r')

AIref = 'TTTTAACAGCTCATTTGCGAATAACTCCCTAAACATGTACTAGTGATAAAGATGCGGTTGACCACAACACAGAGCCCAGAGACGAGGTTACGCACACTTCTGCTTTTCAGGAGAAACAAAGGGGGCTGTCGCCGCGTAATCACGCAGAACGATTTAGGTTTACGCCATCACACAGGGCACCGTTATCGAGGCAGGGACACAGCACACAGTCAGTCATTTGAAGTGTGCGCATTCTTGGACTGGAGTGTGGGCGGTGAGCTTGCAGGGTTCTTTATGCATAATGAGCCAAGGAAACAAAGGCAGGCACAAGAGCAGACAGCATGGGTTATGTTGAGGATTTTTAGAATTTATAATGCTCAGAGGAGAACTGAAAGTGGGTAAAAGTCACTGATAATGCCAATGTTCACTGAGTTGTTCTGTTCATTTTGTTTGTTTGTTTAAATCAGAGACCAGCCTTAATACTGATGACCTCTTTTT'
APref = 'ATGTGCATTTATTGAAAAGCTGATTTTTAACAGACAGCTAGAGCAAATGACACTGATGTTCTCAAACTCAAAGCAAAGCAACTTCTTTTGAAGTTTAAAGTAGTTGAAAACCTTCTCGTCAACAAAATGCCAAATACAAGCTCACGCCAACATGCAACCGCTGAACTAGGTCCTGTAAACCCAAGGGCGCAAGCACAAACCCAATCCCTCAAGGTTACAGCAGTATTGTTTACCTTTTTCCTGCAGTGTTGTGGCTTTTGCATTACCCTGACCTGGCTCGTAACCAGCTCAGACCGCCTATAAAGCTGCCAGCATTCGGACCCTCCGGGTCTGTGCAGGCTGTTCTACATCACCCTTCTTATGGATCTGATCTCTGCTTGTGAACGGGCAATGAGTCCTGTACGCTTAGATGCCGAGGTGTCAGATCTATCCGTTACCTCAAATGCCATCCAATCGCATGGGATATCAATGGCAACCAGA'
B2ref = 'TTAATCATAAAGCTTGATGTGGGTGTAAATGTTTGGTTTAAATGTGGCTCGGTGCTATGTTACGTCCCAGACTGAAGGACAGACTCCACAGAAGATTAGAAAGAAGGTAACAGGAAGGTTACCCACAGGCTGATCATCTTCTGTTTTGTCTGCTGCTTATTAAAAGATCTGAACGCAGACTGCAGACGCACTGTGTGTGTGTGTGTGTGCTCCGTTGGGTGTACCTGGCTAGAGAAACGCAGAGGATGCTGCCGGTAGAGGAGCTCACTGCTGCTGGCCCCNTGGTGGCTGATACAGTGTCCGAAGTCACTGCCGTCCTGCTTCTTCTGCTGCTGTTGTTGCTTTTCACCACCTGGAGACAAAGAAAACAGTCACACATACCAGGTGAGGCACTAATAAATATTTCCCTGTAAGTAACTATATACTGTATACACACACATAACAAATATCATATATATA'
BIref = 'ACTTCACATTTAGCTTTGTTGCACCCAGAGAACCCCCACTGTCCAATTCATCAACCCCGGAGGCTTTTGAAAGGGAAATTAATTAGCTGGGTCAGACTGACCTGTGTCAAAGGAGGCAAGGACTTACTCAGGAGTTCAAAAAGGGTAAAGGAGAATTAAAAAATGTCATTACTATATTCTTAAAAACTCAGTCTTAATAATTAAATCTAAAAGTGCCGTCAAAGTGACTGATTGCATGTTAGTTAGTCCAGTGTAAGAAATCATAGTTGCAAACTTCAAAGGAAAAAACACTGTTCATGCAAAGGTTTTGTGAATGAGCTGGCTTTCAAGCTCGGTACCTGTGAATAGCTGGTATGCCGTGTTTAAAGACGCTGGCCACACATTTAGGGGGGGGATGTGTTGGTCTGGCTCTGAAATTAGCACTAGAAGCTACCCTGACCTCTGTTTCCTCCCCAAAAAGCAGCAAATAAATACAAATTGAGAAAAAAAATCTAATGCTGACTATTTACCACAGTTCAAA'

#print 'sample\tamp_bp\tamplicon\tbp\toriginal_base\tA\tG\tC\tT\tdel\tins\tinserted\tambiguous\tprop_error\tpercent_error'

for line in basecountFile:
	data = line.strip().split('\t')
	if data[1] not in ['amp_bp', 'AI130', 'AI133', 'AI135', 'AI143', 'AI150', 'AI163', 'AI181', 'AI187', 'AI228', 'AI252', 'AP117', 'AP145', 'AP159', 'AP188', 'AP279', 'AP295', 'AP317', 'AP326', 'B2174', 'B2187', 'B2213', 'B2238', 'B2253', 'B2302', 'B2314', 'BI217', 'BI333', 'BI358', 'BI370']:
		data[3] = int(data[3])
		data[4] = int(data[4])
		data[5] = int(data[5])
		data[6] = int(data[6])
		data[7] = int(data[7])
		data[8] = int(data[8])
		data[9] = int(data[9])
		if data[2] == 'AI':
			if AIref[data[3]-1] == 'A':
				if float(data[4])+float(data[5])+float(data[6])+float(data[7]) == 0:
					prop_error = 0
				else:
					prop_error = (float(data[5])+float(data[6])+float(data[7]))/(float(data[4])+float(data[5])+float(data[6])+float(data[7]))
				percent_error = prop_error*100
				data = data[0:4] + ['A'] + data[4:] + [prop_error, percent_error]
				print '\t'.join([str(x) for x in data])
			elif AIref[data[3]-1] == 'G':
				if float(data[4])+float(data[5])+float(data[6])+float(data[7]) == 0:
					prop_error = 0
				else:
					prop_error = (float(data[4])+float(data[6])+float(data[7]))/(float(data[4])+float(data[5])+float(data[6])+float(data[7]))
				percent_error = prop_error*100
				data = data[0:4] + ['G'] + data[4:] + [prop_error, percent_error]
				print '\t'.join([str(x) for x in data])
		elif data[2] == 'AP':
			if APref[data[3]-1] == 'A':
				if float(data[4])+float(data[5])+float(data[6])+float(data[7]) == 0:
					prop_error = 0
				else:
					prop_error = (float(data[5])+float(data[6])+float(data[7]))/(float(data[4])+float(data[5])+float(data[6])+float(data[7]))
				percent_error = prop_error*100
				data = data[0:4] + ['A'] + data[4:] + [prop_error, percent_error]
				print '\t'.join([str(x) for x in data])
			elif APref[data[3]-1] == 'G':
				if float(data[4])+float(data[5])+float(data[6])+float(data[7]) == 0:
					prop_error = 0
				else:
					prop_error = (float(data[4])+float(data[6])+float(data[7]))/(float(data[4])+float(data[5])+float(data[6])+float(data[7]))
				percent_error = prop_error*100
				data = data[0:4] + ['G'] + data[4:] + [prop_error, percent_error]
				print '\t'.join([str(x) for x in data])
		elif data[2] == 'B2':
			if B2ref[data[3]-1] == 'A':
				if float(data[4])+float(data[5])+float(data[6])+float(data[7]) == 0:
					prop_error = 0
				else:
					prop_error = (float(data[5])+float(data[6])+float(data[7]))/(float(data[4])+float(data[5])+float(data[6])+float(data[7]))
				percent_error = prop_error*100
				data = data[0:4] + ['A'] + data[4:] + [prop_error, percent_error]
				print '\t'.join([str(x) for x in data])
			elif B2ref[data[3]-1] == 'G':
				if float(data[4])+float(data[5])+float(data[6])+float(data[7]) == 0:
					prop_error = 0
				else:
					prop_error = (float(data[4])+float(data[6])+float(data[7]))/(float(data[4])+float(data[5])+float(data[6])+float(data[7]))
				percent_error = prop_error*100
				data = data[0:4] + ['G'] + data[4:] + [prop_error, percent_error]
				print '\t'.join([str(x) for x in data])
		elif data[2] == 'BI':
			if BIref[data[3]-1] == 'A':
				if float(data[4])+float(data[5])+float(data[6])+float(data[7]) == 0:
					prop_error = 0
				else:
					prop_error = (float(data[5])+float(data[6])+float(data[7]))/(float(data[4])+float(data[5])+float(data[6])+float(data[7]))
				percent_error = prop_error*100
				data = data[0:4] + ['A'] + data[4:] + [prop_error, percent_error]
				print '\t'.join([str(x) for x in data])
			elif BIref[data[3]-1] == 'G':
				if float(data[4])+float(data[5])+float(data[6])+float(data[7]) == 0:
					prop_error = 0
				else:
					prop_error = (float(data[4])+float(data[6])+float(data[7]))/(float(data[4])+float(data[5])+float(data[6])+float(data[7]))
				percent_error = prop_error*100
				data = data[0:4] + ['G'] + data[4:] + [prop_error, percent_error]
				print '\t'.join([str(x) for x in data])



