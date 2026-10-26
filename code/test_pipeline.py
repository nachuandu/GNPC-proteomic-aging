"""Synthetic-only tests; no restricted study data needed."""
import importlib.util
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from analysis_core import (BINS, bh, classify, wave_membership, read_matrix,
                           align_metadata, gam_r2)
from replicate_st16 import concordance
from trajectory_backend import percentile


class PipelineTests(unittest.TestCase):
    def test_bh_known_values_and_undefined_test(self):
        np.testing.assert_allclose(bh([.01,.04,.03,np.nan]), [.04,.0533333333,.0533333333,1])

    def test_classification_signed_difference_and_boundaries(self):
        labels = classify(np.array([.2,.3,.25,-.1]), np.array([.1001,-.11,.9,.1]))
        self.assertEqual(labels.tolist(), ['Class I','Class II','Unclassified','Unclassified'])
        # Two large opposite-signed correlations differ even when |rho| is equal.
        self.assertEqual(classify(np.array([.4]), np.array([-.8]))[0], 'Class II')

    def test_gene_or_and_all_patterns_retained(self):
        patterns = {'a':'1100','b':'0011','c':'0111','d':'0010','e':'0000'}
        rows = [{'feature':f,'age_bin':label,'R2_deviation':.2 if bit=='1' else .1,
                 'q_pearson':.049} for f, pattern in patterns.items()
                for label, bit in zip(BINS, pattern)]
        ann = pd.DataFrame({'feature':list(patterns),'gene_label':['G','G','H','J','K']})
        members, counts = wave_membership(pd.DataFrame(rows), ann)
        by = members.set_index('gene_label')
        self.assertEqual(by.loc['G','pattern'], '1111')
        self.assertEqual(by.loc['H','wave'], 'Wave 2')
        self.assertEqual(by.loc['J','wave'], 'Other selected pattern')
        self.assertEqual(by.loc['K','wave'], 'Not selected')
        self.assertEqual(counts.n_gene_labels.sum(), 4)

    def test_missing_bin_rejected(self):
        rows = pd.DataFrame({'feature':['a']*3,'age_bin':BINS[:3],
                             'R2_deviation':[.8]*3,'q_pearson':[.01]*3})
        with self.assertRaises(ValueError):
            wave_membership(rows, pd.DataFrame({'feature':['a'],'gene_label':['G']}))

    def test_mapping_ambiguity_not_silently_replicated(self):
        g = pd.DataFrame({'feature':['a','b','c','d'],'delta_rho':[.2]*4,
                          'class':['Class I']*4})
        s = pd.DataFrame({'feature':['x','y'],'delta_rho':[.3,-.2]})
        m = pd.DataFrame({'gnpc_feature':['a','b','c','d'], 'st16_feature':['x','y','z','z']})
        detail, summary = concordance(g,s,m)
        row = summary.iloc[0]
        self.assertEqual(row.n_concordant,1)
        self.assertEqual(row.n_opposite,1)
        self.assertEqual(row.n_unmatched_or_ambiguous,2)
        self.assertEqual(row.percent_concordant,50.)

    def test_ecdf_ties_and_reference_weights(self):
        p,n = percentile(np.array([60.,60.,60.,60.]), np.array([0.,1.,1.,2.]),
                         np.array([60.]), np.array([1.]))
        np.testing.assert_allclose(p,[50.])
        np.testing.assert_allclose(n,[4.])

    def test_explicit_id_alignment_preserves_leading_zeros(self):
        with tempfile.TemporaryDirectory() as folder:
            folder=Path(folder)
            (folder/'x.csv').write_text('sample_id,a\n001,1\n002,2\n')
            (folder/'m.csv').write_text('sample_id,age\n002,80\n001,50\n')
            X=read_matrix(folder/'x.csv')
            m=align_metadata(X,folder/'m.csv',['age'])
            self.assertEqual(m.age.tolist(),[50,80])
            self.assertEqual(X.index.tolist(),['001','002'])

    @unittest.skipUnless(importlib.util.find_spec('pygam'), 'pygam not installed')
    def test_gam_has_separate_r2_and_pearson_test(self):
        x=np.linspace(-2,2,60)
        y=x*x
        self.assertGreater(gam_r2(x,y),.5)
        from analysis_core import correlation
        r,p,_=correlation(x,y)
        self.assertAlmostEqual(r,0,places=10)
        self.assertGreater(p,.9)


if __name__ == '__main__':
    unittest.main(verbosity=2)
