"""
Tools for reading in vcf files.
Frequently Needed Functions
- find the number of samples
	- given the file name and the assumption that the data is organized by 1:
- find the number of 
"""
from collections import deque
import pdb

global bpf,htag,nhfvcf
# by definition
bpf=2; # basepair field
htag='#CHROM' # header tag
nhfvcf=9 # number of header fields in vcf

def open_file(fin):
	"""Check if the file is open, if not open and return the handle."""
	if hasattr(fin, 'closed'):
		fin=open(fin.name, 'r')
	else:
		fin=open(fin, 'r')
	return fin

def header_info(fin, delim):
	"""Return the headers for the file, the number of samples in the file, and the line number. Line by line read of file to prevent running out of memory."""
	global htag,nhfvcf
	fin=open_file(fin)
	it=0
	for line in fin:
		if line[0:len(htag)] != htag:
			it+=1
		else:
			headers=line.strip().split(delim)
			return [headers, len(headers) - nhfvcf, it]

def bof_vcf(fin, delim, header_line):
	global bpf
	fin=open_file(fin)
	it=0
	for line in fin:
		if it <= header_line:
			it+=1
		else:
			break
	return int(line.strip().split(delim)[bpf-1])

def eof_vcf(fin, delim, header_line):
	"""Get the upper bound of the data from the file."""
	global bpf
	fin=open_file(fin)
	def file_len(fin):
		for i, l in enumerate(fin):
			pass
		return i + 1
	eof_line=deque(fin, 1)
	eof=int(eof_line[0].strip().split(delim)[bpf-1])
	eof_check=file_len(open_file(fin)) - header_line
	if eof_check <= eof:
		return eof
	else:
		return eof_check

class windowReader():
	"""Iterate through the previous reading frame of data before reading in new line to the frame (can easily deal with any kind of step size and window size this way). """
	def __init__(self, fin, pdata, pfd, flank, delim, dbp, pfdbp):
		self.fin = fin
		self.pdata = pdata
		self.dbp = deque(dbp, len(pdata)) # holds the data basepair info
		self.pfdbp = [deque(pfdbp[0], flank), deque(pfdbp[1], flank)] # holds the flank data basepair info
		self.pfd = [deque(pfd[0], flank), deque(pfd[1], flank)]
		self.flank = flank
		self.flag = 0
		self.delim = delim

	def increment(self):
		if self.flag == 0: # return the head flanker
			try:
				return self.pfd[0].popleft(),self.pfdbp[0].popleft()
			except IndexError: # until empty
				self.flag = 1
		if self.flag == 1: # return the window data
			try:
				return self.pdata.popleft(),self.dbp.popleft()
			except IndexError: # until empty
				self.flag = 2
		if self.flag == 2: # return the tail flanker
			try:
				return self.pfd[1].popleft(),self.pfdbp[1].popleft()
			except IndexError: # until empty
				self.flag = 3
		if self.flag == 3: # return new line from file
			try:
				return self.readline()
			except EOFError: # until empty? (should be taken care of elsewhere in region_search.py)
				return -1 # not sure what to do specfically

	def readline(self):
		new = self.fin.next().strip().split(self.delim)
#		if new[1] == '2483447':
#			pdb.set_trace()
#			pass
		newbp = new[3:5] # grab the REF / ALT
		new = [x[0] for x in new] # grab the SNPs
		new = new[nhfvcf:]
		return new,newbp 

def read_vcf_records(fin, delim, window, basepair, cp, header_line, data, dbp, flank, bof_line, flank_data=[], fdbp=[]):
	"""Read the vcf file line by line described by the basepair bounds with the python csv module. Return useful data from the file."""
	# works based on keeping the same file open and just grabbing + basepair[2] lines each time
	global htag,nhfvcf
	if len(data) == 0: # first iteration
		for i,l in enumerate(fin): # cut through all the header info on first run
			if i == header_line:
				break
		
		# if want to read a range that is not from the beginning of the file, then skip to that point here
		skip_range = range(bof_line, basepair[0]-flank)
		if len(skip_range) > 0: # if enough room to start before the first basepair
			for x in skip_range: # skip through the unwanted data
				fin.next()
		else: # should be taken care of elsewhere, must have preceding conserved region
			pass # do nothing because specified start coincides with the beginning of the file
			#error='Initial skip_range for the conserved region is of zero length.'
			#sys.exit(error)

		fdbp = [[], []]
		flank_data = [[],[]]

		# read in head flank region (before window region)
		for x in range(0, flank):
			flank_data[0].append(fin.next().strip().split(delim))

		# read in window region
		for x in range(0, window): # fill the entire window, not just increment
			data.append(fin.next().strip().split(delim))
		# grab the REF / ALT
		dbp = deque([x[3:5] for x in data])
		# remove the unnecessary data - just want 1, 0, or .
		data=[x[nhfvcf:] for x in data]
		data=[[el[0] for el in row] for row in data]
		data=deque(data, window)
		
		# read in the tail flank (after window region)
		for x in range(0, flank): # get the after flank
			flank_data[1].append(fin.next().strip().split(delim))
		# grab the REF / ALT here
		fdbp[0] = deque([x[3:5] for x in flank_data[0]])
		fdbp[1] = deque([x[3:5] for x in flank_data[1]])

		flank_data[0]=[x[nhfvcf:] for x in flank_data[0]]
		flank_data[0]=[[el[0] for el in row] for row in flank_data[0]]

		flank_data[1]=[x[nhfvcf:] for x in flank_data[1]]
		flank_data[1]=[[el[0] for el in row] for row in flank_data[1]]
		
		"""
		for el in flank_data[0]:
			print('Head: ' + '\t' + str(el) + '\n')
		print('\n')
		for el in data:
			print('New Window: ' + '\t' + str(el) + '\n')
		print('\n')
		for el in flank_data[1]:
			print('Head: ' + '\t' + str(el) + '\n')
		print('\n')
		"""
	else: # increment the window
		pfd = flank_data # store the previous flank data
		pfdbp = fdbp # store the previous flank_data basepairs
		wr = windowReader(fin, data, pfd, flank, delim, dbp, pfdbp)
		fdbp = [[], []]
		flank_data = [[], []]

		flank_data[0] = deque([], flank)
		flank_data[1] = deque([], flank)
		fdbp[0] = deque([], flank)
		fdbp[1] = deque([], flank)

		count = 0
		while count < basepair[2] + flank: # remove any gap between old window and new and read in the new flank
			next_snp,next_bp = wr.increment()
			flank_data[0].append(next_snp)
			fdbp[0].append(next_bp)
			#print('Head: ' + '\t' + str(flank_data[0][-1]) + '\n')
			count+=1
		#print('\n')
		count = 0
		window_data=deque([],window)
		window_data_bp=deque([],window)
		while count < window: # read in the new window
			next_snp,next_bp = wr.increment()
			window_data.append(next_snp)
			window_data_bp.append(next_bp)
			#print('New Window: ' + '\t' + str(window_data[-1]) + '\n')
			count+=1
		count = 0
		while count < flank: # read in the end of the new flank
			next_snp,next_bp = wr.increment()
			flank_data[1].append(next_snp)
			fdbp[1].append(next_bp)
			#print('Tail: ' + '\t' + str(flank_data[1][-1]) + '\n')
			count+=1
		data=window_data
		dbp=deque(window_data_bp, len(window_data_bp))
	return data,flank_data,dbp,fdbp
