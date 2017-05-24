import argparse

def verify_vcf(f, delim='\t'):
	count=0
	num=[]
	errors={1: []}
	while 1:
		temp=f.next().strip()[0:6]
		if temp != '#CHROM':
			pass
		else:
			break
	while 1:
		try:
			count+=1
			num.append(int(f.next().strip().split(delim)[1]))
			if len(num) > 1 and  num[-1] - num[-2] != 1:
				diff = num[-1] - num[-2] - 1
				if diff in errors.keys():
					errors[diff].append(diff)
					count+=diff
				else:
					errors[diff] = []
					errors[diff].append(diff)
					count+=diff
			if len(num) > 1:
				num.pop(-2)
		except (EOFError, StopIteration) as e:
			message = []
			for x in errors.keys():
				if len(errors[x]) > 0:
					message.append('Number of basepairs between reads: ' + str(x) + ', number of times occured ' + str(len(errors[x])))
					message=['\n'.join(message)]
			message='\n'.join(message)
			break
	return message

def main():
	parser = argparse.ArgumentParser(description=""".""")
	parser.add_argument('-i','--input',help="""Input the multisample VCF file.""",required=True)
	parser.add_argument('-d','--delim',help="""Input the delimiter for the VCF file.""",required=False,default='\t')
	args=parser.parse_args()

	with open(args.input, 'r') as f:
		count=0
		num=[]
		errors={1: []}
		while 1:
			temp=f.next().strip()[0:6]
			if temp != '#CHROM':
				pass
			else:
				break	
		while 1:
			try:
				count+=1
				num.append(int(f.next().strip().split(args.delim)[1]))
				if len(num) > 1 and  num[-1] - num[-2] != 1:
					diff = num[-1] - num[-2] - 1
					#print('Error at vcf bp: ' + str(num[-2]) + ' ' + str(num[-1]) + '\n')
					#print('\tFile line: ' + str(count) + '\n')
					if diff in errors.keys():
						errors[diff].append(diff)
						count+=diff
					else:
						errors[diff] = []
						errors[diff].append(diff)
						count+=diff
				if len(num) > 1:
					num.pop(-2)
			except (EOFError, StopIteration) as e:
				for x in errors.keys():
					print('Number of basepairs between reads: ' + str(x) + ', number of times occured ' + str(len(errors[x])))
				break

if __name__ == "__main__":
	main()
