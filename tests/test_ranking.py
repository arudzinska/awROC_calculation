#!/usr/bin/env python
#
# Command to run all the tests (from tests/ directory): $ python -m unittest discover -v -s ..
#

import unittest
import ranking

class Test_ranking_TestCase(unittest.TestCase):
    
    def test_ranking_initial_attributes(self):
        with open('test_awroc_vs_1st/input_indices.txt') as f:
                indices = f.read().splitlines()
        with open('test_awroc_vs_1st/input_names.txt') as f2:
                names = f2.read().splitlines()
        self.rank = ranking.Ranking('test_awroc_vs_1st/rank_class.txt')
        self.assertEqual(self.rank.total_mols, 1528, "self.total_mols attribute of the ranking is wrong.")
        self.assertEqual(self.rank.indices, indices, "self.indices attribute of the ranking is wrong.")
        self.assertEqual(self.rank.names, names, "self.names attribute of the ranking is wrong.")
        
    def test_calculating_ROC_and_awROC(self):
        self.rank = ranking.Ranking('test_awroc_vs_1st/rank_class.txt')
        ROCE = [43.94117647058824, 23.435294117647057, 14.647058823529413, 7.030588235294118]
        AUC = 0.70810693755413956
        awROCE = [32.95588235294118, 17.576470588235296, 12.694117647058823, 5.412436974789917]
        awAUC = 0.64401021452757845
        first_dict = {'ZINC03832011': 10, 'ZINC00005722': 3, 'ZINC00005723': 1, 'ZINC04617885': 8, 'ZINC00598073': 10, 'ZINC01536588': 6, 'ZINC01492540': 10, 'ZINC03797553': 2, 'ZINC02004030': 6, 'ZINC03833870': 1, 'ZINC01625281': 2, 'ZINC03580965': 15, 'ZINC01905250': 6, 'ZINC00005563': 13, 'ZINC01624808': 2, 'ZINC00004778': 5, 'ZINC00018097': 7, 'ZINC00598801': 16, 'ZINC03832001': 10, 'ZINC03832004': 11, 'ZINC03832005': 3, 'ZINC03832007': 10, 'ZINC03832008': 10, 'ZINC03832009': 10, 'ZINC00014547': 14, 'ZINC03833865': 7, 'ZINC03833867': 6, 'ZINC03833866': 7, 'ZINC03833863': 17, 'ZINC03833862': 12, 'ZINC03833869': 6, 'ZINC03833868': 6, 'ZINC02020188': 4, 'ZINC00006488': 9}
        second_dict = {1: 2, 2: 3, 3: 2, 4: 1, 5: 1, 6: 6, 7: 3, 8: 1, 9: 1, 10: 7, 11: 1, 12: 1, 13: 1, 14: 1, 15: 1, 16: 1, 17: 1}
        with open('test_awroc_vs_1st/actives.txt') as f:
                actives = f.read().splitlines()
        self.rank.ROC_curve(actives, log_scale=False)
        self.rank.awROC_curve(first_dict, second_dict, log_scale=False)
        # AUC and awAUC are rounded to 5 digits after the decimal point (floating point rounding problem)
        self.assertEqual(self.rank.ROCE, ROCE, "Calculated ROCE are not correct")
        self.assertEqual(round(self.rank.AUC, 5), round(AUC, 5), "Calculated AUC is not correct")
        self.assertEqual(self.rank.awROCE, awROCE, "Calculated awROCE are not correct")
        self.assertEqual(round(self.rank.awAUC, 5), round(awAUC, 5), "Calculated awAUC is not correct")
