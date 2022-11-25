#!/usr/bin/env python

import re, sys, cPickle, os, math

def check_opts(opts):
	if opts.input is None:
		die_with_help()
	if not file_exists(opts.input):
		die_with_message("Cannot find file '%s'" % opts.input)
	if opts.unknowns is not None:
		if len(opts.unknowns) < 1:
			die_with_message("Invalid -u option.")
	if opts.split is not None:
		try:
			int(opts.split)
		except ValueError:
			die_with_message("Please provide an integer value for --split option")

def file_exists(f):
	return os.path.exists(os.path.realpath(f))

def check_aln(seqs):
	base = len(seqs[0])
	for s in seqs:
		if len(s) != base:
			die_with_message("Sequences are not of even lengths. Please ensure the data is aligned.")
	return base

def parse_fasta(input):
	data = []
	fa_str = "\n".join(open(input).readlines())
	parts = fa_str.split(">")
	parts = parts[1:]

	for p in parts:
		spl = p.split("\n")
		name = spl.pop(0)
		seq = ""
		for s in spl:
			if re.search("\w", s):
				seq += s.rstrip()
		data.append([name, seq])

	return data

def patterns(data):
	seqs = [i[1] for i in data]

	l = len(seqs[0])
	sites = []
	for x in range(l):
		s = [y[x] for y in seqs]
		sites.append(s)

	return [site_pattern(s) for s in sites]

def site_pattern(site):
	pat = {}
	order = []
	for i, base in enumerate(site):
		try:
			pat[base].append(i)
		except KeyError:
			pat[base] = [i]
			order.append(base)

	return "|".join( [",".join([str(x) for x in pat[o]]) for o in order ])

def pattern_counts_sets(pats):
	uniq = {}
	for x, p in enumerate(pats):
		try:
			uniq[p]["count"] += 1
		except KeyError:
			uniq[p] = {"count": 1}

		try:
			uniq[p]["sites"].append(x)
		except KeyError:
			uniq[p]["sites"] = [x]


	return uniq

def run(opts):
	check_opts(opts)
	seq_data = parse_fasta(opts.input)
	seq_len = check_aln(seq_data)
	pats = patterns(seq_data)

	uniq_pats = pattern_counts_sets(pats)

	prefix = opts.output
	if opts.output is None:
		prefix = gen_prefix(opts.input)
	out = "%s.ref.ti" % prefix
	with open(out, 'wb') as fh:
		#cPickle.dump(seq_data, fh)
		cPickle.dump(uniq_pats, fh)
	
	if opts.split is not None:
		write_subsets(uniq_pats, int(opts.split), prefix)
		
def write_subsets(pats, num, prefix):
	step = math.ceil(float(len(pats))/num)
	c = 1
	subs = []

	pat_keys = pats.keys()
	while len(pat_keys) > 0:
		subs.append({})

		for x in range(int(step)):
			try:
				pat = pat_keys.pop()
				subs[-1][pat] = pats[pat]
				if len(pat_keys) < step:
					pat = pat_keys.pop()
					subs[-1][pat] = pats[pat]
			except IndexError:
				break

	for i in range(len(subs)):
		with open("%s.%d.ti" % (prefix, i), 'wb') as fh:
			cPickle.dump(subs[i], fh)

def gen_prefix(input):
	if "/" in input:
		spl = input.split('/')
		i = spl.pop()
	else:
		i = input

	if '.' in i:
		spl = i.split('.')
		o = spl.pop(0)
	else:
		o = i
	return o

def die_with_help():
    print("""
    ****************
    TIGER  v2.0 Help:
    ****************
    
    tiger index Options:
    
        -i|input    Specify input file. File must be in FastA format and must be aligned prior.
                    Datasets with uneven sequence lengths will return an error.
    
        -s|split    Split dataset across multiple files to run simultaneously. Takes int argument.
    
        -o|output   Specify the prefix name of output files.
    
        -u|unknowns Specify unknown characters in the alignment. Unknown characters are omitted from 
                    site patterns and so are not considered in the analysis.
                    -u ?,-,*: defines ?, - and * as unknown characters. (*Be sure to put only a comma
                    between characters, NO SPACE!!)
             
                    Default is ? only
    
        Examples:
            1. Generate a .tgr file for complete sequence named full_seq.tgr & set unknowns to ? and - :
                    tiger index -i my_file.aln -o full_seq -u ?,-
    
            2. Generate 10 subsets of the data with an output prefix of tiger_split and a reference:
                    tiger index -i my_file.aln -o tiger_split -s 10
                ** Results in files named tiger_split.1.tgr, tiger_split.2,tgr, and so on, along with
                   tiger_split.ref.tgr
       
         """)
    sys.exit(1)

def die_with_message(message):
	print(message)
	sys.exit(1)
	