#!/usr/bin/env python

import re, sys, cPickle, os

def check_opts(opts):
	# check files sanity
	if opts.input is None and opts.combine is None:
		die_with_help()
	if opts.input is not None and opts.combine is not None:
		die_with_message("--input and --combine options cannot be used together")

	if opts.input is not None and not os.path.exists(os.path.realpath(opts.input)):
		die_with_message("Cannot find file '%s'" % opts.input)
	if opts.combine is not None and not os.path.exists(os.path.realpath(opts.combine)):
		die_with_message("Cannot find file '%s'" % opts.combine)

	if opts.fasta is None:
		die_with_message("Please provide the original fasta file (--fasta)")
	if not os.path.exists(os.path.realpath(opts.fasta)):
		die_with_message("Cannot find fasta file %s" % opts.fasta)

	# check format option sanity
	f = opts.format
	if f is None or f not in ['0', '1', '2', '3', '4']:
		die_with_message("Invalid -f value: %s" % f)
	if opts.bins is not None:
		try:
			int(opts.bins)
		except ValueError:
			die_with_message("Invalid -b option '%s'. Please provide an integer" % opt.bins)
	if opts.exclude_only is not None and opts.include_only is not None:
		die_with_message("--include_only and -exclude_only options are mutually exclusive")
	


def bin(rate_d, bin_no):
	rates = [rate_d[k]["rate"] for k in rate_d.keys()]

	upper = float(max(rates))
	lower = float(min(rates))

	step = (upper-lower)/int(bin_no)
	if step == 0:
		die_with_message("The data is too homogeneous to bin! Congratulations...??")
	bin_divs = []
	d = lower
	while d <= upper:
		bin_divs.append(d)
		d += step
	if bin_divs[-1] != upper:
		bin_divs.append(upper)

	bin_d = rate_d.copy()
	for k in rate_d.keys():
		bin_d[k]["bin"] = get_bin(bin_divs, bin_d[k]["rate"])

	return bin_d

def get_bin(parts, rate):
	for i in range(len(parts)-1):
		if rate >= parts[i] and rate <= parts[i+1]:
			return (len(parts)-1)-i

	
#def histogram(num_list, name_list):
def histogram(binned_data):
	histo = ''
	width = 60

	# format the data
	bin_counts = {}
	[total, max_bin] = [0, 0];

	#print binned_data.values()
	for site in binned_data.values():
	 	bin_no = site['bin']
	 	count = site['count']
	 	total += int(count)

	 	if bin_no > max_bin: max_bin = bin_no

		try:
			bin_counts[bin_no] += int(count)
		except KeyError:
			bin_counts[bin_no] = int(count)

	max_bin_len = len(str(max_bin))
	for i in range(1, max_bin+1):
		histo += "bin: %d%s|" % (i, ' ' * (max_bin_len - len(str(i))))
		try:
			this_count = bin_counts[i]
			num_cols = int((float(this_count)/total) * width)
		except KeyError: # bin empty
			[this_count, num_cols] = [0,0]

		pad_cols = width - num_cols
		histo += "%s%s| %d\n" % ('=' * num_cols, ' ' * pad_cols, this_count)

	return histo

def generate_nexus(binned_data, seq_data, excl_bins, mask, format_opt):
	nexus = ''

	species_order = sorted(seq_data.keys())

	nexus += nexus_header(binned_data, species_order, seq_data[species_order[0]])

	padded_speceies = pad_str(seq_data.keys())
	positions = split_fasta_into_positions(seq_data, species_order)

	if format_opt == 1 or format_opt == 3: # sites require sorting
		pos_order = sorted_position_order(binned_data)
	else:
		pos_order = range(len(positions))

	# create comments
	# print matrix
	# print charsets

	return nexus

def sorted_position_order(binned_data):
	sorted_pos = []

	rate_map = {}
	for x in binned_data:
		rate_map[x['rate']] = x['sites']

	sorted_rates = sorted(rate_map.keys())
	for r in sorted_rates:
		for p in sorted(rate_map[r]):
			sorted_pos.append(p)

	return sorted_pos

def split_fasta_into_positions(seq_data, species_order):
	pos_list = []

	for i in range(len(seq_data[species_order[0]])):
		pos = ''
		for species in species_order:
			pos += seq_data[species][i]
		pos_list.append(pos)

	return pos_list

def nexus_header(binned_data, species, rep_seq):
	nex_header = ''
	nex_header += "#NEXUS\n\n[This file contains data that has been analysed for site specific rates]\n"
	nex_header += "[using TIGER, developed by Carla Cummins in the laboratory of]\n"
	nex_header += "[Dr James McInerney, National University of Ireland, Maynooth]\n\n"
	
	nex_header += "[Histograms of number of sites in each bin:]\n"        
	nex_header += histogram(binned_data)

	nex_header += "\n\n\nBEGIN TAXA;\n"
	nex_header += "\tDimensions NTax = %d;\n" % len(species_order)
	nex_header += "\tTaxLabels %s;\nEND;\n" % " ".join(species_order)
	
	nex_header += "BEGIN CHARACTERS;\n"
	nex_header += "\tDimensions nchar = %d;\n" % len(rep_seq)
	nex_header += "\tFormat datatype = %s gap = - interleave;\nMatrix\n\n" % detect_datatype(rep_seq)

	return nex_header

def detect_datatype(seq):
	seq = seq.upper()
	oLen = float(len(seq))

	seq_C = re.sub('[ACGT]+', '', seq)
	nLen = float(len(seq_C))
	perc = (nLen/oLen)*100

	if perc < 20.0:
		return "DNA"
	else:
		seq_C = seq.replace("[01]+", "", seq)
		if len(seq_C) == 0:
			return "standard"
		else:
			return "protein"

def pad_str(taxon_names):
	max_len = 0;
	for name in taxon_names:
		if len(name) > max_len:
			max_len = len(name)

	padded = {}
	for name in taxon_names:
		pad_len = (max_len - len(name)) + 3
		pad_str = ' ' * pad_len
		padded[name] = "%s%s" % (name, pad_str)

	return padded

def generate_fasta(binned_data, seq_data, excl_bins, mask):
	bin_map = map_bins_to_positions(binned_data)
	fasta_str = ''

	for species in seq_data.keys():
		masked_seq = ''
		for pos,base in enumerate(seq_data[species]):
			if bin_map[pos] in excl_bins:
				masked_seq += 'X'
			else:
				masked_seq += base
		
		if not mask:
			final_seq = re.sub('X', '', masked_seq)
		else:
			final_seq = masked_seq

		fasta_str += ">%s\n" % species
		fasta_str += "%s\n" % final_seq

	return fasta_str

def map_bins_to_positions(binned_data):
	bin_map = {}

	for x in binned_data.values():
		cur_bin = x['bin']
		for pos in x['sites']:
			bin_map[pos] = cur_bin

	return bin_map

def combine_rates(combine_file):
	combine_fh = open(combine_file, 'r')

	all_rates = {}
	for rate_file in combine_fh:
		rate_file = rate_file.rstrip()
		with open(rate_file, 'rb') as fh:
			these_rates = cPickle.load(fh)
			all_rates.update(these_rates)

	return all_rates

def bins_to_exclude(opts):
	excl_bins = []
	if ( opts.exclude_only is not None ):
		excl_bins = [int(b) for b in opts.exclude_only.split(',')]
	elif ( opts.include_only is not None ):
		incl_bins = [int(b) for b in opts.include_only.split(',')]
		for b in range(1,int(opts.bins)):
			if b not in incl_bins:
				excl_bins.append(b)

	return excl_bins

def parse_fasta(fasta_file):
	seq_data = {}
	seq_fh = open(fasta_file, 'r')
	for line in seq_fh:
		line = line.rstrip()
		if line[0] == '>':
			cur_species = line[1:]
			seq_data[cur_species] = ''
		else:
			seq_data[cur_species] += line

	return seq_data

def run(opts):
	check_opts(opts)

	if ( opts.input is not None ):
		with open(opts.input, 'rb') as fh:
			rates = cPickle.load(fh)
	elif ( opts.combine is not None ):
		rates = combine_rates(opts.combine)

	binned_data = bin(rates, opts.bins)

	# create list of bins to mask/remove
	excl_bins = bins_to_exclude(opts)

	# parse original fasta for seq data
	seq_data = parse_fasta(opts.fasta)

	# write output in specified format
	if ( opts.format == '4' ): # output fasta
		output_str = generate_fasta(binned_data, seq_data, excl_bins, opts.mask)
		suffix = '.fa'
	else:
		output_str = generate_nexus(binned_data, seq_data, excl_bins, opts.mask, opts.format)
		suffix = '.nex'

	if ( opts.output is not None ):
		try:
			outfile = open("%s%s" % (opts.output, suffix), 'w')
			outfile.write(output_str)
		except IOError:
			die_with_message("Cannot open outfile %s%s" % (opts.output, suffix))
	else:
		print output_str


def die_with_help():
	print """
****************
TIGER  v2.0 Help:
****************

tiger output Options:
	
	-i|input            Specify input file. Must be in .gr format.

	-c|combine          Specify input file. This file should contain a list of .gr files to be combined.

	-fa|fasta           Provide original .fa file for sequence data.

	-o|output           Specify prefix name for output files (prints to stdout if not provided)

	-f|format           Changes formatting options.
						NEXUS, with comments:
						-f 0: output bin numbers, sites unsorted (default)
						-f 1: output bin number, sites sorted on rank
						-f 2: displays rank values rather than bin numbers
						-f 3: displays rank values and sites sorted on rank
						FastA:
						-f 4

	-inc|include_only   Give list of bins to include (only) in output
						-inc 3,4,5,6 (Note: No spaces, just commas)

	-exc|exclude_only   Give list of bins to exclude from output
						-exc 1,2,9,10

	-m|mask             Mask -inc/-exc sites. (Default is to remove them)

	-b|bins             Set the number of bins to be used.
						-b <int>: Sites will be placed into <int> number of bins. <int> is a whole number.

						Default is 10.
	Examples:
		1.  Write a FastA file, masking site that fall into Bin1, Bin2, Bin9 and Bin10 of 10 bins:
			tiger output -i sample.gr -fa my_data.fa -f 4 -exc 1,2,9,10 -b 10 --mask

		2. Write a NEXUS file combining test.0.gr, test.1.gr, test.2.gr with sites sorted on rank
			tiger output -c list_of_gr_files.txt -fa my_data.fa -f 3
   
	 """
	sys.exit(1)

def die_with_message(message):
	print message
	sys.exit(1)