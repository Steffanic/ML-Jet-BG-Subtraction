from GetAndPrepareData import *
from modelPreparation import *
from plotData import *
import  pandas as pd
import unittest
from numpy import dtype

class TestAll(unittest.TestCase):
    
    def setUp(self):
        self.dat = get_data("../Analysis_Code/Generator Output/merged-ML-output-LOWSTATS-Rparam-%1.1f-pThardmin-%d.0.csv"%(0.4, 20), 10000)
    
    def testTyping(self):
        self.dat_typed = type_data(self.dat)
        self.assertCountEqual(list(self.dat_typed.dtypes), [dtype('float32'), dtype('float32'), dtype('float32'), dtype('float32'), dtype('float32'), dtype('float32'), dtype('uint32'), dtype('float32'), dtype('float32'), dtype('float32'), dtype('float32'), dtype('float32'), dtype('float32'), dtype('float32'), dtype('float32'), dtype('float32'), dtype('float32'), dtype('float32'), dtype('float32'), dtype('float32'), dtype('float32'), dtype('float32')])

    def testAddFeats(self):
        self.dat_feat_added = add_feats(type_data(self.dat))
        self.assertTrue(self.dat_feat_added['pythia-mom'].min()>=0)

    def testLabeling(self):
        self.dat_labeled = label_data(add_feats(type_data(self.dat)))
        self.assertLessEqual(self.dat_labeled['Label'].all(), 3)
        self.assertGreaterEqual(self.dat_labeled['Label'].all(), 1)

    def testDrop(self):
        self.dat_drop = drop_feat(label_data(add_feats(type_data(self.dat))))
        self.assertLess(len(self.dat_drop.columns), len(label_data(add_feats(type_data(self.dat))).columns))

    def testBalancing(self):
        self.dat_balanced = balance_classes(drop_feat(label_data(add_feats(type_data(self.dat)))))
        self.assertEqual(len(self.dat_balanced.loc[self.dat_balanced['Label']==2]),len(self.dat_balanced.loc[self.dat_balanced['Label']==1]))
        self.assertEqual(len(self.dat_balanced.loc[self.dat_balanced['Label']==2]),len(self.dat_balanced.loc[self.dat_balanced['Label']==3]))

    def testFeatLabelSplitting(self):
        '''
        Test the split of train into X and Y. Assumes that sklearn.model_selection.train_test_split works perfectly
        '''
        train, _ = split_data(balance_classes(drop_feat(label_data(add_feats(type_data(self.dat))))))
        self.X, self.Y = split_feat_label(train)
        self.assertEqual(type(self.Y), pd.Series)
        self.assertEqual(len(self.X.columns), len(train.columns)-3)

'''
TODO: Stop being trash and write the rest of the tests.
'''

if __name__ == '__main__':
    unittest.main()