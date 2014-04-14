import sys, os
sys.path.append(os.path.realpath('./modules'))
sys.path.append(os.path.realpath('.'))

import biotiger
import pytest, filecmp, tempfile, shutil

# create temp dir and write config pointing to it
tmpdir = tempfile.mkdtemp(dir=os.getcwd())
confpath = "%s/test.conf" % tmpdir
conf = open(confpath, 'w')
conf.write("db_location\t%s\n" % tmpdir)
conf.close()

base_args = [ '-d', 'resistome_testing', '-c', confpath]

fasta = ""
patterns = []
unique_patterns = []

def test_parse():
    global fasta
    input_file = "tests/data/sample.fa"
    exp_fasta = [['tax1', 'AAAAAAAAAA'],
    	         ['tax2', 'AAATTAAAAA'],
    	         ['tax3', 'AAATTCAAAT']]
    fasta = biotiger.parse_fasta(input_file)
    
    assert fasta == exp_fasta

def test_pattern():
    site = ['A', 'A', 'A', 'G', 'G', 'T']
    exp = "0,1,2|3,4|5"
    got = biotiger.site_pattern(site)

    assert got == exp

def test_pattern_list():
    global patterns
    exp = [
        '0,1,2',
        '0,1,2',
        '0,1,2',
        '0|1,2',
        '0|1,2',      
        '0,1|2',
        '0,1,2',
        '0,1,2',
        '0,1,2',
        '0,1|2'
    ]
    patterns = biotiger.patterns(fasta)

    assert patterns == exp

def test_pattern_counts():
    global unique_patterns
    exp = {
        '0,1,2' : 6,
        '0|1,2' : 2,
        '0,1|2' : 2    
    }
    unique_patterns = biotiger.pattern_counts(patterns)

    assert unique_patterns == exp

#def test_rates():
#    exp = [
#        1.0,
#        1.0,
#        1.0,
#        0.75,
#        0.75,
#        0.75,
#        1.0,
#        1.0,
#        1.0,
#        0.75
#    ]
#    rates = biotiger.rate_list(patterns)
#
#    assert exp == rates

def test_set():
    exp = {
    '0,1,2' : [set([0,1,2])],
    '0|1,2' : [set([0]), set([1,2])]
    }
    pats = ['0,1,2', '0|1,2']
    got = biotiger.set_pattern(pats)

    print got

    assert got == exp

def test_sort():
    patterns = ['A', 'B', 'C', 'D']
    rates = [0.4, 0.1, 0.3, 0.2]

    exp_patterns = ['A', 'C', 'D', 'B']
    exp_rates = [0.4, 0.3, 0.2, 0.1]

    [got_rates, got_patterns] = biotiger.sort(rates, patterns)

    assert got_rates == exp_rates
    assert got_patterns == exp_patterns

    cleanup()
        
#def test_save():
#        args = ['-f', "%s/test.csv" % tmpdir]  + base_args
#        opts = parse_args(args)
#        manage_pickle.run('save', opts)
#        
#        assert filecmp.cmp("%s/test.csv" % tmpdir, 'test/data/exp.csv')

def cleanup():	
        shutil.rmtree(tmpdir)