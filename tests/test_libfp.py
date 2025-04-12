#!/usr/bin/env python3

import unittest
import numpy as np
import os
import libfp

# Path to test files
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
EXAMPLES_DIR = os.path.join(os.path.dirname(TEST_DIR), 'examples')


def read_vasp(vp):
    """Read VASP structure file."""
    buff = []
    with open(vp) as f:
        for line in f:
            buff.append(line.split())

    lat = np.array(buff[2:5], float)
    try:
        typt = np.array(buff[5], int)
    except:
        del(buff[5])
        typt = np.array(buff[5], int)
    nat = sum(typt)
    pos = np.array(buff[7:7 + nat], float)
    types = []
    for i in range(len(typt)):
        types += [i+1]*typt[i]
    types = np.array(types, int)
    rxyz = np.dot(pos, lat)
    return lat, rxyz, types


class TestLibFP(unittest.TestCase):
    """Test cases for libfp functionality."""

    def setUp(self):
        """Set up test structures."""
        self.znucl = [14, 8]  # Si, O
        
        # Load example structures
        struct1_path = os.path.join(EXAMPLES_DIR, 'struct1.vasp')
        struct2_path = os.path.join(EXAMPLES_DIR, 'struct2.vasp')
        
        self.lat1, self.rxyz1, self.types1 = read_vasp(struct1_path)
        self.cell1 = (self.lat1, self.rxyz1, self.types1, self.znucl)
        
        if os.path.exists(struct2_path):
            self.lat2, self.rxyz2, self.types2 = read_vasp(struct2_path)
            self.cell2 = (self.lat2, self.rxyz2, self.types2, self.znucl)
        else:
            # Create slight variation of struct1 for testing if struct2 doesn't exist
            self.lat2 = self.lat1.copy()
            self.rxyz2 = self.rxyz1.copy()
            # Add small perturbation
            self.rxyz2 += np.random.normal(0, 0.01, self.rxyz2.shape)
            self.types2 = self.types1.copy()
            self.cell2 = (self.lat2, self.rxyz2, self.types2, self.znucl)

    def test_version(self):
        """Test version retrieval."""
        version = libfp.get_version()
        self.assertIsInstance(version, tuple)
        self.assertEqual(len(version), 3)
        for v in version:
            self.assertIsInstance(v, int)

    def test_get_lfp(self):
        """Test local fingerprint calculation."""
        fp1 = libfp.get_lfp(self.cell1, cutoff=5.0, log=True, orbital='s')
        self.assertIsInstance(fp1, np.ndarray)
        self.assertGreater(len(fp1), 0)
        
        # Test with sp orbital
        fp2 = libfp.get_lfp(self.cell1, cutoff=5.0, log=True, orbital='sp')
        self.assertIsInstance(fp2, np.ndarray)
        self.assertGreater(len(fp2), 0)
        
        # sp orbital fingerprint should be longer than s orbital
        self.assertGreater(len(fp2), len(fp1))

    def test_get_sfp(self):
        """Test structural fingerprint calculation."""
        fp = libfp.get_sfp(self.cell1, cutoff=5.0, log=True, orbital='s')
        self.assertIsInstance(fp, np.ndarray)
        self.assertGreater(len(fp), 0)

    def test_get_dfp(self):
        """Test fingerprint derivative calculation."""
        fp, dfp = libfp.get_dfp(self.cell1, cutoff=5.0, log=True, orbital='s')
        self.assertIsInstance(fp, np.ndarray)
        self.assertIsInstance(dfp, np.ndarray)
        self.assertGreater(len(fp), 0)
        self.assertGreater(len(dfp), 0)

    def test_fp_distance(self):
        """Test fingerprint distance calculation."""
        fp1 = libfp.get_lfp(self.cell1, cutoff=5.0, log=True, orbital='s')
        fp2 = libfp.get_lfp(self.cell2, cutoff=5.0, log=True, orbital='s')
        
        dist = libfp.get_fp_dist(fp1, fp2, self.types1)
        self.assertIsInstance(dist, float)
        self.assertGreaterEqual(dist, 0.0)
        
        # Test with assignment
        dist, assignment = libfp.get_fp_dist(fp1, fp2, self.types1, assignment=True)
        self.assertIsInstance(dist, float)
        self.assertIsInstance(assignment, np.ndarray)

    def test_input_validation(self):
        """Test input validation."""
        # Test with invalid cell
        invalid_cell = (self.lat1, self.rxyz1, self.types1[:5], self.znucl)
        with self.assertRaises(ValueError):
            libfp.get_lfp(invalid_cell)


if __name__ == '__main__':
    unittest.main()