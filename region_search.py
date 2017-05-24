"""PCR Primer Discrimination Window Analysis Tool"""
__author__="Brody DeSilva"
__email__="bdesilva@uab.edu"

import argparse
import sys
import pdb
import vcf_tools as vcf_tools
from collections import deque
from math import ceil
from verify_vcf_increment_by_one import verify_vcf

"""
Implement the deque data holder for the window view (pop/pull each old/new element)
Identify way to store the various unique patterns
"""
global fin,delim,fout,header_line,eof_line,bof_line,unknown,percent,output_limit

def check_values(**kwargs):
	global fin,delim,eof_line,bof_line
	# input a variable to test format uniformity, returns nothing, just exits w/error if no pass
	keys=sorted(kwargs.keys())
	for kw in kwargs:
		if kw == 'basepair':
			if kwargs[kw] == None:
				error='Basepair argument cannot be \'None\'. Ensure that basepair was input correctly.\n'
				sys.exit(error)
			if len(kwargs[kw]) > 3:
				error='Input at most 2 numbers for the basepair.\n'
				sys.exit(error)
			if kwargs[kw][0] < 0:
				error='Starting basepair must be at or beyond the first basepair of the sequence.\n'
				sys.exit(error)
			if kwargs[kw][1] > eof_line:
				error='Final basepair must be either at or before the last basepair of the sequence.\n'
				sys.exit(error)
		if kw == 'window':
			if kwargs[kw] > eof_line:
				error='Window cannot be larger than the file to be read in.\n'
				sys.exit(error)
			try:
				if kwargs['basepair'][1] < kwargs[kw]:
					error='Window cannot be larger than the range in the file to be read in.\n'
					sys.exit(error)
			except SystemExit:
				sys.exit(error)

		if kw == 'unknown':
			if type(kwargs[kw]) is not str:
				error='Ensure that the unknown value is of type string.'
				sys.exit(error)
		if kw == 'output_limit':
			if kwargs[kw] < 1:
				error='Ensure that the output limit is one number and is greater than or equal 1.'
				sys.exit(error)
		if kw == 'flank':
			if kwargs[kw] < 0:
				error='Ensure that the flank is one number and is greater than or equal to 0.'
				sys.exit(error)
#			if kwargs['basepair'][-1] < kwargs[kw]:
#				error='Ensure that the flank is less than the step size.'
#				sys.exit(error)
			if kwargs['basepair'][1] - kwargs['basepair'][0] < kwargs[kw]:
				error='Ensure that the flank is less than the size of the region to be searched.'
				sys.exit(error)
			if kwargs['basepair'][0] - bof_line - kwargs[kw] < 0:
				warning='Must have space for the conserved region before the start of the window. Changing the starting position to be at least ' + str(bof_line + kwargs[kw]) + '.'
				print(warning)
				return ['preflank', bof_line + kwargs[kw]]
		if kw == 'percent':
			if kwargs[kw] < 0:
				error='Ensure that the entered percent is greater than 0.'
				sys.exit(error)
			if kwargs[kw] > 100:
				error='Ensure that the entered percent is less than 100.'
				sys.exit(error)

def set_values(**kwargs):
	global fin,delim,eof_line,bof_line
	keys=sorted(kwargs.keys())

	for kw in kwargs:
		if kw == 'basepair':
			# input will be a list of two elements for each argument: example basepair=[[0, 0, 0], [220, 100, 300]]
			# insert check to this specific keyword here
			if kwargs[kw][1] != None:
				if kwargs[kw][1][0] == -1: # set to beginning of file
					kwargs[kw][1][0]=bof_line
				if kwargs[kw][1][1] == -1: # set to end of file
					kwargs[kw][1][1]=eof_line
				if kwargs[kw][1][2] == -1: # set to end of file, not useful, just for completeness
					kwargs[kw][1][2]=eof_line
				if len(kwargs[kw][1]) > 0:
					kwargs[kw][0][0]=kwargs[kw][1][0]
				if len(kwargs[kw][1]) > 1:
					kwargs[kw][0][1]=kwargs[kw][1][1]
				if len(kwargs[kw][1]) > 2:
					kwargs[kw][0][2]=kwargs[kw][1][2]
					if kwargs[kw][0][2] == -1:
						kwargs[kw][0][2]=eof_line
				return kwargs[kw][0] # usage is the returned value will always
			else:
				return kwargs[kw][0] # return the default value as the new value for default
		elif kw == 'window' or kw == 'delimiter' or kw == 'unknown' or kw == 'output_limit' or kw == 'flank' or kw == 'percent':
			if kwargs[kw] != None:
				return kwargs[kw]

class Album:
	# collects data about which snapshot has what features, etc.
	def __init__(self):
		self.data = []
		pass
	def add_snapshot(self, snapshot):
		# add snapshot to album and organize / update the album
			# updating will entail at least a record of all the ranges for each of the types of data
			self.data.append(snapshot)
	def organize(self):
		# organize the album, allow for various organization methods
		# organize each individual snapshot by largest data
		data=[[x.score for x in self.data],[x.unfilt_score for x in self.data]]
		index=sorted(range(len(data[0])), key=data[0].__getitem__)
		index.reverse()
		data=[[data[0][x] for x in index], [data[1][x] for x in index]]
		return [index, data]
	def limit_output(self):
		global output_limit
		remove_queue=[]
		index=self.organize()[0]
		count=0
		for i in index:
			count+=1
			if count > output_limit:
				try:
					remove_queue.append(i)
				except:
					pass
		remove_queue.sort(reverse=True) # sort the elements so that the index is preserved for earlier elements
		for i in remove_queue: # remove the elements in place (index changes, but already accounted for), circumvent copy each time element is deleted
			del self.data[i]

class Snapshot:
	def __init__(self, cp, size, data, dbp):
		# use the window as data ref to capture the snapshot
		self.bounds=[cp, cp+size]
		self.score,self.unfilt_score,self.data,self.filtered_total=self.read_data_and_calculate(zip(*list(data)))
		self.dbp = deque(dbp)
		self.bpm = [] # basepair matrix
	def read_data_and_calculate(self, data):
		"""Data is given by each list is an individual sample."""
		global unknown,percent
		# concatenate the data within each sample
		for x in range(0, len(data)):
			data[x] = ''.join(data[x])
		# convert to set
		# unique the set and count how many are unique
		new_data = list(data)
		for el in new_data:
			if ceil(float(el.count(unknown))/len(el)*100) > percent:
				new_data = [y for y in new_data if y != el]
		
		unique_data = list(set(data))
#		unique_data_copy = list(unique_data)
		unfilt_len = len(unique_data)
		# post-processing for '.'
#		for el in unique_data:
#			test_el = list(el[flank:-flank])
#			if ceil(float(test_el.count(unknown))/len(test_el)*100) > percent:
#				unique_data_copy = [y for y in unique_data_copy if y != el]
#		unique_data = unique_data_copy
#		if len(unique_data) > 10:
#			pdb.set_trace()
#		if len(set(new_data)) > 45:
#			pdb.set_trace()
		return [len(set(new_data)), unfilt_len, data, len(new_data)]

	def snpToBp(self):
		if len(self.data[0]) != len(self.dbp):
			sys.exit('Regions not the same length, error occurring somewhere.')
		
		data_t = zip(*self.data) # data transposed
		
		self.bpm = [[x[0].replace('<NON_REF>', '.') for x in self.dbp], [x[1].replace('<NON_REF>', '.') for x in self.dbp], []] # REF, ALT, sample-specific
		
		# fix the dbp if a 2 or 3 SNP response
		for x in range(0, len(self.bpm[1])):
			if self.bpm[1][x].count(','):
				self.bpm[1][x] = self.bpm[1][x].split(',')
		
		for sample in range(0, len(self.data)):
			self.bpm[-1].append([])
	
		# store in the x-basepair, y-sample format for ease of printing
		for bp in range(0, len(self.dbp)):
			for sample in range(0, len(self.data)):
				try:
					cur = int(data_t[bp][sample])
					if cur == 0: # option \t0:, no SNP
						self.bpm[2][sample].append(self.bpm[0][bp])
					else: # dealing with multiple other options
						self.bpm[2][sample].append(self.bpm[1][bp][cur-1]) # if option of \t1:, \t2:, or \t3:, SNP
				except ValueError:
					self.bpm[2][sample].append('.')

	def printData(self, fout, flank, sample_names=[], out=0):
		"""Print out the snap data, optionally specify what data to be printed."""

		def printObjVCF(snap, fout, flank, printobj, bpm=0):
			fout.write('Window Bounds (flank inclusive): ' + str(snap.bounds[0] - flank) + '\t' + str(snap.bounds[1] + flank) + '\n')
			for bp in range(0, len(snap.data[0])):
				fout.write(str(snap.bounds[0] + bp - flank))
				if bpm == 1:
					fout.write('\t' + str(snap.bpm[0][bp]) + '\t' + str(','.join(snap.bpm[1][bp])))
				for sample in range(0, len(snap.data)):
					fout.write('\t' + printobj[sample][bp])
				fout.write('\n')
				if bp == flank-1 or bp == len(snap.data[0])-1 - flank:
					fout.write('\n')
			fout.write('\n')

		def printObjRow(snap, fout, flank, printobj, sample_names, bpm=0):
			fout.write('Bounds (flank inclusive): ' + str(snap.bounds[0]-flank) + '\t' + str(snap.bounds[1]+flank) + '\n')
			max_sample_name = max([len(x) for x in sample_names]) + 1
			if max_sample_name < 4:
				max_sample_name = 4
			fout.write('REF' + ' '*(max_sample_name-3))
			for bp in range(0, len(snap.data[0])):
				fout.write(str(snap.bpm[0][bp]))
				if bp == flank-1 or bp == len(snap.data[0])-1 - flank:
					fout.write('\t')

			fout.write('\nAlt' + ' '*(max_sample_name-3))
			for bp in range(0, len(snap.data[1])):
				if len(snap.bpm[1][bp]) > 1:
					fout.write('[' + ','.join(snap.bpm[1][bp]) + ']')
				else:
					fout.write(str(snap.bpm[1][bp]))
				if bp == flank-1 or bp == len(snap.data[0])-1 - flank:
					fout.write('\t')

			fout.write('\n')
			for sample in range(0, len(snap.data)):
				fout.write('\n' + str(sample_names[sample]) + ' '*(max_sample_name-len(sample_names[sample])))
				for bp in range(0, len(snap.data[0])):
					fout.write(str(snap.bpm[2][sample][bp]))
					if bp == flank-1 or bp == len(snap.data[0])-1 - flank:
						fout.write('\t')
			fout.write('\n\n')

		if out == 0: # all data
			# snp
			printObjVCF(self, fout, flank, self.data)
			# bp
			printObjVCF(self, fout, flank, self.bpm[-1], bpm=1)
		elif out == 1: # snp data
			printObjVCF(self, fout, flank, self.data)
		elif out == 2: # bp data
			printObjVCF(self, fout, flank, self.bpm[-1], bpm=1)
		elif out == 3: # snp, sample-based as rows
			printObjRow(self, fout, flank, self.bpm[-1], sample_names, bpm=1)
		elif out == 4: # FASTA REF Output with flanking regions included
			fout.write('FASTA REF Output with Flanking Regions\n')
			fout.write(str(self.bounds[0]-flank) + '-' + str(self.bounds[1]+flank) + '\n')
			for bp in range(0, len(self.data[0])):
				fout.write(str(self.bpm[0][bp]))
			fout.write('\n\n')
		elif out == 5: # FASTA REF Output without flanking region
			fout.write('FASTA REF Output without Flanking Regions\n')
			fout.write(str(self.bounds[0]) + '-' + str(self.bounds[1]) + '\n')
			for bp in range(flank, len(self.data[0])-flank):
				fout.write(str(self.bpm[0][bp]))
			fout.write('\n\n')
	
	def addFlankToSnapshot(self, head, tail, head_bp, tail_bp):
		global unknown
		head = dataToString(head)
		tail = dataToString(tail)
		# since not getting results, need to more specifically limit the data
		# basically transpose the data, then look at each basepair and identify if there is more than one number present
		# ignore the unknown for this as some samples have very little data in general
		if not flank_data_is_conserved(head) or not flank_data_is_conserved(tail):
			self.score = 0
			self.unfilt_score = 0
			del self.data
			return
		# set the data
		for x in range(0, len(self.data)):
			self.data[x] = head[x] + self.data[x] + tail[x]
		# set the basepair data
		for x in range(0, len(head_bp)):
			self.dbp.appendleft(head_bp.pop())
		for x in range(0, len(tail_bp)):
			self.dbp.append(tail_bp.popleft())


def dataToString(data):
	for x in range(0, len(data)):
		data[x] = ''.join(data[x])
	return data

def flank_data_is_conserved(fd):
	global unknown
	line=[]
	for bp in fd[0]:
		line.append([])
	for sample in fd:
		for bp in range(0, len(sample)):
			line[bp].append(sample[bp])
	for bp in line:
		unique_data = []
		for sample in bp:
			if sample != unknown and unique_data.count(sample) == 0:
				unique_data.append(sample)
			if len(unique_data) > 1: # any size 2 combination of 0, 1, or 2
				return False
	return True


def main():
	global fin,delim,fout,eof_line,bof_line,unknown,percent,output_limit
	# Use argparse to grab the data from the command line
	parser = argparse.ArgumentParser(description=""".""")
	parser.add_argument('-i','--input',help="""Input the multisample VCF file.""",type=argparse.FileType('r'),required=True)
	parser.add_argument('-o','--output',help=""".""",type=argparse.FileType('w'),required=True)
	parser.add_argument('-v','--verbose',help=""".""",action='store_true')
	parser.add_argument('-bp','--basepair',help="""Choose the beginning and (optional) ending base pair for the window to begin searching (use -1 if you want the last basepair in the file). If no ending base pair is chosen, will default to the final base pair in the genome. If no increment is chosen, will default to default increment. With one option, will default to 1:end by default_increment.\nEx: '-bp 0 300' or '-bp 0' or '-bp 0 300 100'""",nargs='+',type=int)
	parser.add_argument('-w','--window',help="""Choose the length (in basepairs) of the window.\nEx: '-w 300'""",type=int)
	parser.add_argument('-d','--delimiter',help="""Specify a delimiter for your input and output files, otherwise use tab as delimiter.\nEx: '-d \t'""",type=str,default='\t')
	parser.add_argument('-u','--unknown',help="""Specify the unknown data character.\nEx: '.'""",type=str,default='.')
	parser.add_argument('-l','--output_limit',help="""Limit output controls the number of entries that are printed out to the output file.""",type=int,default=50)
	parser.add_argument('-f','--flank',help="""Input the length in basepairs of the flankers for post-identification of the conservative regions for the data.\nEx: '-f 100' lead of 100 bp and flank of 100 bp.""", type=int, default=100)
	parser.add_argument('-p','--percent',help="""Input the percent of the unknown value allowed in a single sample's current window. Anything below this percent will be filtered out and the sample's data will not contribute to the score.\nEx: 20""",type=int,default=0)
	parser.add_argument('-s','--verify',help="""Verify the VCF file before running, exit upon missing records.""",type=int,default=0)

	args=parser.parse_args()

	# set the input and output files (already checked by default in the parser)
	fin=args.input
	fout=args.output
	fin=vcf_tools.open_file(fin)

	if args.verify != 0:
		# Check the VCF File for missing records
		message=verify_vcf(fin)
		if len(message) > 0:
			sys.exit(message)

	# FIX THE DELIM FROM BEING '\\t'
	# set delim (must be before the basepair)
	if hasattr(args, 'delimiter') and args.delimiter != None:
		delim=set_values(delimiter=args.delimiter)
	else:
		delim='\t'
	# check delim
	check_values(delimiter=delim)
	# set header_line
	header=vcf_tools.header_info(fin, delim)
	header_line=header[2]
	# set eof_vcf
	eof_line=vcf_tools.eof_vcf(fin, delim, header_line)
	# set bof_vcf
	bof_line=vcf_tools.bof_vcf(fin, delim, header_line)
	# set bp
	if hasattr(args, 'basepair') and args.basepair != None:
		bp=set_values(basepair=[[-1, eof_line, 1], args.basepair])
	else:
		bp=[0, eof_line, 1]
	# check bp
	check_values(basepair=bp)
	# set window
	if hasattr(args, 'window') and args.window != None:
		window=set_values(window=args.window)
	else:
		window=500
	# check window with the basepair
	check_values(window=window,basepair=bp)
	
	if hasattr(args, 'unknown') and args.unknown != None:
		unknown=set_values(unknown=args.unknown)
	else:
		unknown='.'
	check_values(unknown=unknown)
	
	if hasattr(args, 'output_limit') and args.output_limit != None:
		output_limit=set_values(output_limit=args.output_limit)
	else:
		output_limit=50
	check_values(output_limit=output_limit)

	if hasattr(args, 'flank') and args.flank != None:
		flank=set_values(flank=args.flank)
	else:
		flank=100
	
	value=check_values(flank=flank, basepair=bp, window=window)
	if value is not None and value[0] == 'preflank':
		bp[0] = value[1]
	
	if hasattr(args, 'percent') and args.percent != None:
		percent=set_values(percent=args.percent)
	else:
		percent=0
	check_values(percent=percent)
	# not sure if necessary
	#if flank > bp[-1]:
	#	error='Flank must be smaller than the step size.'
	#	sys.exit(error)
	
	# initialize other variables
	descrip_tag='##' # unique tag that identifies description lines
	header_tag='#CHROM' # the tag that identifies the header line

	it=0 #iterator
	alt_it = 0
	# ensure fresh file - reading from the beginning
	fin=vcf_tools.open_file(fin)
	data = deque([], window)
	dbp = deque([], window) # ensure is set up correctly
	
	[data, flank_data, dbp, fdbp]=vcf_tools.read_vcf_records(fin, delim, window, bp, it, header_line, data, dbp, flank, bof_line) # the SNP data and basepair letter
	album=Album()
	while it + bp[0] + window + bp[2] + flank < bp[1]:
		# inspect the fdbp and dbp variables
		album.add_snapshot(Snapshot(it+bp[0], window, data, dbp)) # + bp[0] bc bounds must be run independent
		album.data[-1].addFlankToSnapshot(zip(*list(flank_data[0])), zip(*list(flank_data[1])), deque(fdbp[0], flank), deque(fdbp[1], flank)) # if conserved region, snapshots added to the end, else delete snap
		if album.data[-1].score == 0:
			album.data.pop()
		if alt_it > output_limit*2:
			alt_it = 0
			album.limit_output() # deletes the data for the snapshots that are not in the top output_limit
		else:
			alt_it+=1
		it+=bp[2]
		# update the current position and the data, then create new snapshot with the new "window" (data incremented)
		try:
			[data, flank_data, dbp, fdbp]=vcf_tools.read_vcf_records(fin, delim, window, bp, it, header_line, data, dbp, flank, bof_line, flank_data, fdbp)
			#pdb.set_trace()
		except StopIteration:
			# end of file
			pdb.set_trace()
			break
	# rearrange the snapshots in order of highest filtered data
	album.data = [album.data[x] for x in album.organize()[0]]
	# prepare the basepairs for printing out
	# if 1 '<NON-REF>' and 1 bp
		# change the non to '.'
		# the above is only if the assumption that '<NON-REF>' occurs only when there is no SNP at a certain bp position, to check, make print out ! if non-ref and a SNP for a sample
	 #if 2 '<NON-REF>' (i.e. vcf verification script filled in the basepair)
	 	# put two '.'s
	# each snapshot,
		# each sample
			# each bp
				# print out bp based on SNP or not

	# write to file
	fout.write('#' + fin.name + '\n')
	fout.write('\t'.join(header[0])+'\n')
	fout.write('#Number of samples:' + '\t' + str(header[1])+'\n')
	fout.write('#start|stop|increment\t' + str(bp) + '\n')
	fout.write('#window size:\t' + str(window) + '\n')
	fout.write('#flank size:\t' + str(flank) + '\n')
	fout.write('#unknown vcf data character:\t"' + str(unknown) + '"\n')
	fout.write('#filter percent:\t' + str(percent) + '\n')
	fout.write('#filtered score\tfiltered total samples\tunfiltered score\tbounds (flank inclusive)\n')
	try:
		i=0
		for el in album.data:
			fout.write(str(el.score) + '\t' + str(el.filtered_total) + '\t' + str(el.unfilt_score) + '\t' + str(el.bounds[0] - flank) + '\t' + str(el.bounds[1] + flank) + '\n')
			if i >= output_limit:
				break
			else:
				i+=1
		fout.write('\n')

		# write the stored data to file in the same form of the vcf file : y axis is basepair, x axis is sample
		# separate by each snap by the bounds of the snap
		i=0
		for snap in album.data:
			if i >= output_limit:
				break
			else:
				# write the basepair matrix to file
				snap.snpToBp()
				snap.printData(fout, flank)
				snap.printData(fout, flank, out=5)
				snap.printData(fout, flank, out=4)
				snap.printData(fout, flank, header[0][-len(snap.data):], 3)

				i+=1
	except TypeError: # deal with no output returning valid results
		pass
			
	fin.close()
	fout.close()

if __name__ == "__main__":
	main()
