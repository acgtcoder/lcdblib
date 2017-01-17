"""Tests for odo backends."""

import os
import pytest
from math import ceil
from odo import odo, resource
import pandas as pd
from lcdblib.odo_backends.rseqc import *


@pytest.fixture
def infer_experiment_se(tmpdir):
    data = (
    'This is SingleEnd Data\n'
    'Fraction of reads failed to determine: 0.0964\n'
    'Fraction of reads explained by "++,--": 0.1195\n'
    'Fraction of reads explained by "+-,-+": 0.7841')
    fn = os.path.join(str(tmpdir), 'test_infer_experiment.txt')
    with open(fn, 'w') as fh:
        fh.write(data)
    return fn

@pytest.fixture
def infer_experiment_pe(tmpdir):
    data = (
        'This is PairEnd Data\n'
        'Fraction of reads failed to determine: 0.2530\n'
        'Fraction of reads explained by "1++,1--,2+-,2-+": 0.0126\n'
        'Fraction of reads explained by "1+-,1-+,2++,2--": 0.7345')
    fn = os.path.join(str(tmpdir), 'test_infer_experiment.txt')
    with open(fn, 'w') as fh:
        fh.write(data)
    return fn

def test_infer_experiment_se(infer_experiment_se):
    dat = resource(infer_experiment_se)
    assert dat[0] == 0.0964
    assert dat[1] == 0.1195
    assert dat[2] == 0.7841
    assert dat[3] == 'SingleEnd'

def test_infer_experiment_pe(infer_experiment_pe):
    dat = resource(infer_experiment_pe)
    assert dat[0] == 0.2530
    assert dat[1] == 0.0126
    assert dat[2] == 0.7345
    assert dat[3] == 'PairEnd'

@pytest.fixture
def geneBody_coverage(tmpdir):
    data = (
        'Percentile	1	2	3	4	5	6	7	8	9	10\n'
        'SRR9999999.prealn	46.0	70.0	115.0	154.0	224.0	246.0	275.0	343.0	418.0	457.0')
    fn = os.path.join(str(tmpdir), 'test_geneBody_coverage.txt')
    with open(fn, 'w') as fh:
        fh.write(data)
    return fn

def test_geneBody_coverage(geneBody_coverage):
    dat = resource(geneBody_coverage)
    assert dat.iloc[0, 0] == 46.0
    assert dat.columns[0] == '1'
    assert dat.index[0] == 'SRR9999999'
    assert dat.index.name == 'sample_id'


@pytest.fixture
def tin(tmpdir):
    data = (
        'Bam_file	TIN(mean)	TIN(median)	TIN(stdev)\n'
        'SRR9999999.prealn.bam	25.6425321019	22.2814373475	17.9599093927')
    fn = os.path.join(str(tmpdir), 'test_tin.txt')
    with open(fn, 'w') as fh:
        fh.write(data)
    return fn

def test_tin(tin):
    dat = resource(tin)
    assert ceil(dat.iloc[0, 0]) == 26
    assert dat.columns[0] == 'TIN(mean)'
    assert dat.index[0] == 'SRR9999999'
    assert dat.index.name == 'sample_id'


@pytest.fixture
def bam_stat(tmpdir):
    data = (
        'Load BAM file ...  Done\n'
        '\n'
        '#==================================================\n'
        '#All numbers are READ count\n'
        '#==================================================\n'
        '\n'
        'Total records:                          3759051\n'
        '\n'
        'QC failed:                              0\n'
        'Optical/PCR duplicate:                  0\n'
        'Non primary hits                        0\n'
        'Unmapped reads:                         0\n'
        'mapq < mapq_cut (non-unique):           0\n'
        '\n'
        'mapq >= mapq_cut (unique):              3759051\n'
        'Read-1:                                 1892634\n'
        'Read-2:                                 1866417\n'
        '''Reads map to '+':                       1880063\n'''
        '''Reads map to '-':                       1878988\n'''
        'Non-splice reads:                       3484515\n'
        'Splice reads:                           274536\n'
        'Reads mapped in proper pairs:           3474678\n'
        'Proper-paired reads map to different chrom:0')

    fn = os.path.join(str(tmpdir), 'test_bam_stat.txt')
    with open(fn, 'w') as fh:
        fh.write(data)
    return fn

def test_bam_stat(bam_stat):
    dat = resource(bam_stat)
    assert dat[0] == 3759051
    assert dat[-1] == 0
