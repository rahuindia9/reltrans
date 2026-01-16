import pytest
import numpy as np
from wrapper import DCP_Parameters


def test_basic_invocation(reltrans, assert_snapshot):
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.dcp(energy, DCP_Parameters())
    assert_snapshot(output)
