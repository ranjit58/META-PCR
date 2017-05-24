import argparse
import pdb
import os
import sys
from collections import deque

global file_offset

def make_folder(path):
	try:
		os.makedirs(path)
		return 1
	except OSError as e:
		if e[0] == 17:
			return 1
		else:
			return 0
def write_to_file(fin, data_file, fout, delim='\t'):
	end_bounds=[]
	window_size=0

	with open(fin, 'r') as f:
		while 1: # get window size
			line=f.next().strip()
			try:
				if len(line) >= 7 and line[0:7] == '#window':
					window_size=int(line.split(delim)[-1])
					break
			except:
				sys.exit('Failure to get window.')
		while 1: # get flank size
			line=f.next().strip()
			try:
				if len(line) >= 6 and line[0:6] == '#flank':
					flank_size=int(line.split(delim)[-1])
					break
			except:
				sys.exit('Failure to get flank.')
		while 1: # skip headers
			line=f.next().strip()
			try:
				if len(line) >= 9 and line[0:9] == '#filtered':
					break
			except:
				sys.exit('Failure to find bound header tag.')
		while 1: # read until no more bounds
			line=f.next().strip()
			if len(line) == 0:
				break
			else:
				end_bounds.append(int(line.split(delim)[-1]))
	with open(fout, 'w+') as f: # write to file
		header=retrieve_header_from_vcf(data_file)
		for h in header: # write each line to output file
			f.write(h)
		indices=set()
		for bound in end_bounds: # write out the snapshot reference bounds to file
			lines=retrieve_lines_from_file(data_file, bound-window_size-2*flank_size-1, bound-1, indices)
			for line in lines:
				f.write(line)

def retrieve_lines_from_file(fin, start, stop, exclusions):
	global file_offset
	lines=[]
	with open(fin, 'r') as f:
		for i, line in enumerate(f):
			if i >= start + file_offset and i < stop + file_offset and i not in exclusions: # prevent repetition of lines
				lines.append(line)
				exclusions.add(i)
			elif i >= stop + file_offset:
				break
			else:
				pass
	return lines

def retrieve_header_from_vcf(fin, file_type='vcf'):
	global file_offset
	header=[]
	if file_type == 'vcf':
		with open(fin, 'r') as f:
			for i, line in enumerate(f):
				if i < file_offset:
					header.append(line)
				else:
					break
	return header

def main():
	global file_offset
	parser = argparse.ArgumentParser(description="""Enter the paths to the output from region_search.py to export them to an IGV formatted track file.""")
	parser.add_argument('-i','--input',help="""Input the path to the topmost level directory for which all files underneath it (in any subdirectory) are to be converted into IGV tracks.""",type=str,required=True)
	parser.add_argument('-e','--exclusions',help="""Input the path (from the topmost level directory) to any subdirectories for which IGV output should not be made. Multiple paths can be input at once.""",nargs='+',type=str)
	parser.add_argument('-o','--output',help="""Input the path to the output folder. If it does not exist, it will be created.""",type=str,required=True)
	parser.add_argument('-d','--data',help="""Input the data folder (with data files within) for which the output files will be matched - commonly a VCF file.""",type=str,required=True)
	parser.add_argument('-f','--offset',help="""Input the number of lines reserved for headers in the data file, commonly 125 for VCF files.""",type=int,required=False,default=125)
	parser.add_argument('-t','--type',help="""Input the data file type, (example: .vcf).""",type=str,required=False,default='.vcf')

	args=parser.parse_args()
	
	in_dir=args.input
	out_dir=args.output
	excl=args.exclusions
	data_dir=args.data
	file_offset=args.offset
	filetype=args.type
	
	if filetype[0] != '.':
		filetype='.'+filetype
	
	if data_dir[-1] != os.sep:
		data_dir=data_dir+os.sep

	filtered_excl = list(excl)
	flag = 0
	for i, e in enumerate(excl):
		if e.find(os.sep) == 0:
			e=e[1:]
			filtered_excl[i]=e
		if e.rfind(os.sep) == len(e)-1:
			e=e[:-1]
			filtered_excl[i]=e
	excl=list(filtered_excl)
	del filtered_excl

	files=[]
	dirs=[]
	for root, dirnames, filenames in os.walk(in_dir):
		for d in dirnames:
			dirs.append(root.split(in_dir)[-1] + os.sep + d)
		for filename in filenames:
				files.append(os.path.join(root, filename))
	
	# remove from the dirs any dir that has one of the excluded dirs in it (prevent it from being created
	filtered_dirs=[]
	flag=0
	for d in dirs:
		if d.find(os.sep) == 0:	# remove the '/' at the beginning
			d = d[1:]
		if d.find(os.sep) == len(d)-1: # remove a '/' if it occurs at the end
			d = d[:-1]
		for e in excl:
			if d == e:
				flag=1
		if flag==0:
			filtered_dirs.append(d)
		flag=0
	dirs=list(filtered_dirs)
	del filtered_dirs
	for d in dirs: # make output directories
		make_folder(os.path.join(out_dir, d))
	make_folder(out_dir) # in case no subdirs
	
	# remove from files any file that has one of the excluded dirs in it	
	filtered_files=[]
	flag=0
	for f in files:
		for e in excl:
			if f.find(e) != -1:
				flag=1
		if flag==0:
			filtered_files.append(f)
		flag=0
	files=list(filtered_files)
	del filtered_files

	for f in files:
		out_file=f.split(in_dir)[-1]
		basename=out_file.split(os.sep)[-1].split('.')[0]
		if out_file.find(os.sep) == 0:
			out_file=out_file[1:]
		out_file=os.path.join(out_dir, basename+filetype)
		data_file=data_dir+basename+os.sep+basename+filetype
		write_to_file(f, data_file, out_file)
	
if __name__ == "__main__":
	main()
