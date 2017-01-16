"""Tests for odo backends."""

import os
import pytest
from odo import odo, resource
import pandas as pd
from lcdblib.odo_backends import rseqc


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
