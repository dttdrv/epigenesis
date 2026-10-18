import math
from fractions import Fraction as F
import unittest

import scoring as check


def exact_beta_binomial(a,b,n):
    return [F(math.comb(n,k)*math.prod(range(a,a+k))*math.prod(range(b,b+n-k)),
              math.prod(range(a+b,a+b+n))) for k in range(n+1)]


class ReplayTests(unittest.TestCase):
    def test_scalar_baseline_and_pair_score(self):
        for a,b,n in ((1,1,1),(1,4,9),(5,1,8),(3,7,17),(100,300,17)):
            masses = exact_beta_binomial(a,b,n)
            self.assertEqual(sum(masses),1)
            actual = check.pooled_mass([(b-1,a-1)],n)
            check.same(actual,list(map(float,masses)),'exact beta-binomial',1e-13)
            for observed in (0,n//2,n):
                first = sum(p*F(abs(k-observed),n) for k,p in enumerate(masses))
                second = sum(p*q*F(abs(k-j),2*n)
                             for k,p in enumerate(masses) for j,q in enumerate(masses))
                self.assertAlmostEqual(check.pair_loss(actual,observed),float(first-second),places=13)
                cdf = list(map(float,check.prefix(masses)[:-1]))
                bounds = check.enclosure_loss(cdf,cdf,observed)
                check.same(bounds,[float(first-second)]*2,'CRPS grid endpoints',1e-13)

    def test_raw_population_moments(self):
        weights = [F(0)]*200+[F(1,4)]*300+[F(1,2)]*250+[F(1)]*250
        sums, squares = [75.,125.,0.,250.], [18.75,62.5,0.,250.]
        result = check.summarize(sums,squares,1000,1.,1,jobs=12,delta=.01,tolerance=0.)
        self.assertEqual(result['raw_particles'],1000)
        self.assertAlmostEqual(result['mean_weight'],.45,places=14)
        average = sum(weights)/len(weights)
        variance = sum((x-average)**2 for x in weights)/(len(weights)-1)
        rounding = 4*1000*check.sys.float_info.epsilon
        self.assertAlmostEqual(result['variance'][0],float(variance)+rounding,places=14)
        level = math.log(4*12*4/.01)
        radius = math.sqrt(2*(float(variance)+rounding)*level/1000)+7*level/(3*999)+rounding
        check.same(result['denominator_interval'],[.45-radius,.45+radius],'raw N interval',1e-14)
        self.assertLess(result['effective_particles'],1000)
        loss0 = check.summarize(sums,squares,1000,1.,0,tolerance=0.)
        self.assertEqual(result['cdf'],loss0['cdf'])
        self.assertNotEqual(result['loss'],loss0['loss'])

    def test_zero_and_invalid_moments(self):
        result = check.summarize([0.,0.],[0.,0.],1000,1.,0,tolerance=0.)
        self.assertIsNone(result['cdf'])
        self.assertIsNone(result['loss'])
        self.assertEqual(result['cdf_lower'],[0.])
        self.assertEqual(result['cdf_upper'],[1.])
        self.assertEqual(result['loss_interval'],[0.,1.])
        for sums,squares in (([1.,1.],[0.,0.]),([1.,0.,1.],[0.,.5,.5]),([-1.,0.],[1.,0.])):
            with self.assertRaises(ValueError):
                check.summarize(sums,squares,10,1.,0)

    def test_signs_and_comparison(self):
        pooled,bird,interval = .04,.03,[.02,.05]
        gain, bounds = pooled-bird,[pooled-interval[1],pooled-interval[0]]
        self.assertGreater(gain,0)
        self.assertLess(bounds[0],0)
        self.assertGreater(bounds[1],0)
        with self.assertRaises(AssertionError):
            check.same([.1,.2],[.1,.2001],'deliberate result change')


if __name__ == '__main__':
    unittest.main()
