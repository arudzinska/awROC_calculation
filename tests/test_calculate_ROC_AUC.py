#!/usr/bin/env python
#
# Command to run all the tests (from tests/ directory): $ python -m unittest discover -v -s ..
#

import unittest
import os
import calculate_ROC_AUC as roc

class  Test_calculate_ROC_AUC_TestCase(unittest.TestCase):
    
    def setUp(self):
        self.actives = open(os.path.join(rootdir,"test_awROC_awAUC-files/actives-test.sdf"),'r').readlines()
        self.molecs_in_clusters_result = {'Monica': 2, 'horse': 3, 'Peter': 1, 'Jenny': 2, 'Veronica': 2, 'Pauline': 2, 'Michael': 1, 'Lucy': 2, 'dog': 3, 'cat': 3, 'Andy': 1, 'monkey': 3, 'elephant': 3, 'Paul': 1, 'John': 1, 'mouse': 3, 'Sophie': 2, 'Katie': 2}
        self.struct_in_clusters_result = {1: 5, 2: 7, 3: 6}

    def test_process_actives_and_clusters_from_input_file(self):
        self.assertEqual(roc.process_actives(self.actives), (self.molecs_in_clusters_result, self.struct_in_clusters_result), "Reading the active molecules and their clusters did not succeed (function: process_actives).");

rootdir = os.getcwd()

if __name__ == '__main__':
    unittest.main()
