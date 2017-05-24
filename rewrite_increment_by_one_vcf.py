import pdb
from vcf_tools import bpf,htag,nhfvcf,header_info,bof_vcf
import argparse

# argument parser
parser = argparse.ArgumentParser(description="""""")
parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True)
parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True)
parser.add_argument('-d','--delimiter',type=str,default='\t')
args=parser.parse_args()


fin=args.input
fout=args.output
delim=args.delimiter

# figure out the header information for the file so can replicate a blank record
[headers,num_samples,header_line]=header_info(fin, delim)

template = ['name', 0, '.', '.', '<NON_REF>', '.', '.', 'MLEAC=.;MLEAF=.', 'GT:AD:DP:RGQ']
blank_record='.:0,0:0:0' # make sure that this is true for VCF standard (or at least for the data being used now)
[template.append(blank_record) for x in range(0,num_samples)] # append correct number of blank records

count=0
num=[]
errors={1: []}

for x in range(0, header_line+1): # skip the header lines, get to data
	fout.write(fin.next())

while 1:
	try:
		count+=1
		line = fin.next().strip().split(delim) # read in as csv
		if count==1:
			template[0]=line[0]
		num.append(int(line[bpf-1]))
		if len(num) > 1 and  num[-1] - num[-2] != 1:
			diff = num[-1] - num[-2] - 1
			print('Error at vcf bp: ' + str(num[-2]) + ' ' + str(num[-1]) + '\n')	
			for x in range(1, diff+1):
				template[1]=str(num[-2]+x) # generated record position
				fout.write(delim.join(template)+'\n')
		
			fout.write(delim.join(line)+'\n') # write the next record position
		
			if diff in errors.keys():
				errors[diff].append(diff)
				count+=diff
			else:
				errors[diff] = []
				errors[diff].append(diff)
				count+=diff
		else:
			fout.write(delim.join(line)+'\n') # write existing line to file

		if len(num) > 1:
			num.pop(-2)
	except (EOFError, StopIteration) as e:
		for x in errors.keys():
			print('Number of basepairs between read: ' + str(x) + ', number of times occured ' + str(len(errors[x])))
		break
