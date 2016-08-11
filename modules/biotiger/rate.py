import cPickle, sys, os, re

def check_opts(opts):
    if opts.input is None:
        die_with_help()
    if not os.path.exists(os.path.realpath(opts.input)):
        die_with_message("Cannot find input file '%s'" % opts.input)

    if opts.reference is None:
        opts.reference = opts.input

    if not os.path.exists(os.path.realpath(opts.reference)):
        die_with_message("Cannot find reference file '%s'" % opts.reference)

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


def rate_sites(pat_counts, ref_counts):
    rate_d = pat_counts.copy()
    for k in pat_counts.keys():
        if "|" not in k:
            rate_d[k]["rate"] = 1.0
        else:
            rate_d[k]["rate"] = site_rate(k, ref_counts)

    # rates = []
    # for p in pat_counts.keys():
    #   rates.append(rate_d[p])

    return rate_d


def site_rate(a, ref_counts):
    pat_rates = []
    dividand = 0

    patA = set_pattern(a)

    for p in ref_counts.keys():
        sub = 0
        if '|' in p:       # only score against non-constant sites
            patB = set_pattern(p)
            pr = 1.0 - score(patA, patB)

            reps = ref_counts[p]["count"]
            if a == p:
                reps -= 1
        
            for i in range(reps):
                #print "%s v %s = %f" % (patA, p, pr)
                pat_rates.append(pr)
                dividand += 1

    #print "%f/%d = %f" % (sum(pat_rates), dividand, (sum(pat_rates)/dividand))
    return 1.0 - (sum(pat_rates)/dividand)

def score(a, b):
    found = 0
    sets = len(b)

    for set_b in b:
        for set_a in a:
            if set_b.issubset(set_a):
                found += 1
                break

    return float(found)/sets

def set_pattern(p):
    return [set([int(y) for y in x.split(',')]) for x in p.split('|')]

def sort(rates, patterns):
    if len(rates) != len(patterns):
        print "Something's weird here. len(rates) != len(patterns)"
        sys.exit(1)

    rate_d = {}
    for i, r in enumerate(rates):
        rate_d[i] = r

    sort_order = sorted(rate_d, key=rate_d.get, reverse=True)

    return [ [ rates[o] for o in sort_order ], [ patterns[p] for p in sort_order ] ]

def rate_list(pats):
    rates = {}
    for k in pats.keys():
        r = pats[k]['rate']
        for s in pats[k]['sites']:
            rates[s] = r

    rate_list = []
    for i in range(len(pats)+1):
        rate_list.append(rates[i])

    return rate_list


def run(opts):
    check_opts(opts)

    with open(opts.input, 'rb') as in_h:
        pat_counts = cPickle.load(in_h)
    if opts.reference is None:
        ref_counts = pat_counts.copy()
    else:
        with open(opts.reference, 'rb') as ref_h:
            ref_counts = cPickle.load(ref_h)

    rates = rate_sites(pat_counts, ref_counts)
    print rates

    if opts.output is None:
        prefix = gen_prefix(opts.input)
    else:
        prefix = opts.output

    # write out .gr pickle
    with open("%s.gr" % prefix, 'wb') as fh:
        cPickle.dump(rates, fh)

    # write rate list if required
    if opts.rate_list:
        ratel = rate_list(rates)
        if opts.rate_list is None:
            rl = "%s.rates" % prefix
        else:
            rl = opts.rate_list
        rlh = open(rl, 'w')
        rlh.write( "\n".join([str(x) for x in ratel]) )
        rlh.close()

    # do PTP test if required
    # if opts.run_ptp:
    #   continue


def gen_prefix(input):
    if "/" in input:
        spl = input.split('/')
        i = spl.pop()
    else:
        i = input

    if '.' in i:
        spl = i.split('.')
        spl.pop()
        o = ".".join(spl)
    else:
        o = i
    return o

def die_with_help():
    print """
****************
TIGER  v2.0 Help:
****************

tiger rate Options:

    -i|input            Specify input file. File should be in .ti format.

    -r|reference        Specify reference sequence (.ti). -i file is used as default if none is provided.

    -o|output           Specify prefix name for output files.
    
    -rl|rate_list       A list of the rate at each site may be optionally written to a specified
                        file. 
                        -rl <file.txt> : writes list of the rates at each site to file.txt.

    -ptp                Specifies that a PTP test should be run. 
                        * Note * this option has a huge effect on running time!

    -z|randomisations   Number of randomisations to be used for the PTP test. 
                        -z <int>: each site will be randomised <int> times. <int> is a whole number.

                        Default is 100

    -p|p_value          Specify p-value which denotes significance in PTP test.
                        -p <float>: site will be denoted as significant if p-value is better than <float>.
                        <float> is a floating point number.

                        Default is 0.05

    -pl|pval_list       Write a list of p-values to a specified file.
                        -pl <file.txt>: writes list of p-values for each site to file.txt.

    Examples:
        1. Calculate rates for file test.ref.ti against itself, with a list of rates:
            tiger rate -i test.ref.ti -rl
        2. Calculate rates for file test.0.ti against test.ref.ti with a PTP test and a list of p values
            tiger rate -i test.0.ti -r test.ref.ti -ptp -pl
      
     """
    sys.exit(1)

def die_with_message(message):
    print message
    sys.exit(1)



