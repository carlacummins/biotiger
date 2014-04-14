def check_opts(opts):
	if opts.input is None:
		die_with_help()
	if not os.path.exists(os.path.realpath(opts.input)):
		die_with_message("Cannot find file '%s'" % opts.input)
	f = opts.format
	if f is None or f not in ['0', '1', '2', '3']:
		die_with_message("Invalid -f value: %s" % f)
	if opts.bins is not None:
		try:
			int(opts.bins)
		except ValueError:
			die_with_message("Invalid -b option '%s'. Please provide an integer")
	if opts.run_ptp is not None:
		if opts.rands is not None:
			try:
				int(opts.rands)
			except ValueError:
				die_with_message("Invalid -z option '%s'. Please provide an integer")
		if opts.p_value is not None:
			try:
				float(opts.p_value)
			except ValueError:
				die_with_message("Invalid -p option '%s'. Please provide a floating point number")

def bin(bin_no, rate_d):
	print rate_d
	rates = [rate_d[k]["rate"] for k in rate_d.keys()]

	print rates

	upper = float(max(rates))
	lower = float(min(rates))

	step = (upper-lower)/bin_no
	print "upper: %s, lower: %s, step: %s" % (upper, lower, step)
	if step == 0:
		die_with_message("Something is up")
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
	print parts
	for i in range(len(parts)-1):
		if rate >= parts[i] and rate <= parts[i+1]:
			return i+1

	
def histogram(num_list, name_list):
    upper = float(max(num_list))
    pad = len(str(upper))
    parts = []
    for i in range(1,61):
        parts.append((upper/60)*i)

    for m, n in enumerate(num_list):
        pr = name_list[m] + "|"
        low = 0.0
        if n == 0:
            pr = pr + " "*61
        for p, hi  in enumerate(parts):
            if n > low and n <= hi:
                pr = pr + "="*(p+1) + (" "*(60 - p))
                break
            low = hi
        print "[" + pr + "|" + str(n) + " "*(pad-len(str(n))) + "]"



def print_nexus(taxon_names, formatted_pats, formatted_comm):
	print "#NEXUS\n\n[This file contains data that has been analysed for site specific rates]"
	print "[using TIGER, developed by Carla Cummins in the laboratory of]"
	print "[Dr James McInerney, National University of Ireland, Maynooth]\n\n"
	
	print "[Histograms of number of sites in each category:]"        
	
	histogram(hcounts, hnames)
	filled_names = pad_str(taxon_names)

	print "\n\n"
	
	print "\n\n\nBEGIN TAXA;"
	print "\tDimensions NTax = ", len(seqs), ";"
	print "\tTaxLabels ", " ".join(filled_names), ";\nEND;\n"
	
	print "BEGIN CHARACTERS;"
	print "\tDimensions nchar = ", len(seqs[0]), ";"
	print "\tFormat datatype = ", datatype, " gap = - interleave;\nMatrix\n"

def run(opts):
	check_opts(opts)
	bins = bin(opts.bins, rates)

def die_with_help():
    print """
****************
TIGER  v2.0 Help:
****************

tiger output Options:
    
    -i|input            Specify input file. Must be in .gr format.

    -c|combine          Specify input file. This file should contain a list of .gr files to be combined.

    -f|fasta            Provide original .fa file for sequence data.

    -o|output           Specify prefix name for output files.

    -f|format           Changes formatting options.
                        NEXUS, with comments:
                        -f 0: output bin numbers, sites unsorted (default)
                        -f 1: output bin number, sites sorted on rank
                        -f 2: displays rank values rather than bin numbers
                        -f 3: displays rank values and sites sorted on rank
                        FastA:
                        -f 4

    -inc|include_only   Give list of charsets to include
                        -inc Bin3,Bin4,Bin5,Bin6 (Note: No spaces, just commas)

    -exc|exclude_only   Give list of charsets to exclude
                        -exc Bin1,Bin2,Bin9,Bin10

    -m|mask             Mask -inc/-exc sites. (Default is to remove them)

    -b|bins             Set the number of bins to be used.
                        -b <int>: Sites will be placed into <int> number of bins. <int> is a whole number.

                        Default is 10.
    Examples:
        1.  Write a FastA file, masking site that fall into Bin1, Bin2, Bin9 and Bin10 of 10 bins:
            tiger output -i sample.gr -f my_data.fa -exc Bin1,Bin2,Bin9,Bin10 -b 10 --mask

        2. Write a NEXUS file combining test.0.gr, test.1.gr, test.2.gr with sites sorted on rank
            tiger output -c list_of_gr_files.txt -fa my_data.fa -f 3
   
     """
    sys.exit(1)

def die_with_message(message):
	print message
	sys.exit(1)