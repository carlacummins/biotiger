import sys, os
sys.path.append(os.path.realpath('../modules'))
sys.path.append(os.path.realpath('./modules'))

from biotiger import *
import pytest, filecmp, tempfile, shutil

# create temp dir for writing tests
tmpdir = tempfile.mkdtemp(dir=os.getcwd())

fasta = ""
patterns = []
unique_patterns = []

def test_parse():
    global fasta
    input_file = "tests/data/sample.fa"
    exp_fasta = [['tax1', 'AAAAAAAAAA'],
    	         ['tax2', 'AAATTAAAAA'],
    	         ['tax3', 'AAATTCAAAT']]
    fasta = index.parse_fasta(input_file)
    
    assert fasta == exp_fasta

def test_pattern():
    site = ['A', 'A', 'A', 'G', 'G', 'T']
    exp = "0,1,2|3,4|5"
    got = index.site_pattern(site)

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
    patterns = index.patterns(fasta)

    assert patterns == exp

def test_pattern_counts():
    global unique_patterns
    exp = {
        '0,1,2' : {'count': 6, 'sites': [0, 1, 2, 6, 7, 8]},
        '0|1,2' : {'count': 2, 'sites': [3, 4]},
        '0,1|2' : {'count': 2, 'sites': [5, 9]}
    }
    unique_patterns = index.pattern_counts_sets(patterns)

    assert unique_patterns == exp

# def test_rates():
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
#    rates = rate.rate_list(unique_patterns, unique_patterns)

#    assert exp == rates

def test_set():
    exp = [set([0]), set([1,2])]
    got = rate.set_pattern('0|1,2')

    assert got == exp

def test_sort():
    patterns = ['A', 'B', 'C', 'D']
    rates = [0.4, 0.1, 0.3, 0.2]

    exp_patterns = ['A', 'C', 'D', 'B']
    exp_rates = [0.4, 0.3, 0.2, 0.1]

    [got_rates, got_patterns] = rate.sort(rates, patterns)

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