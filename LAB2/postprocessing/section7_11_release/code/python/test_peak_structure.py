"""Synthetic tests are numerical QA, never experimental observations."""
import unittest
import numpy as np
from peak_structure import analyse,same_pair,sidelobe_conflict

class TopologyTests(unittest.TestCase):
    def setUp(self):
        self.x=np.linspace(-.5,.5,101);self.X,self.Y=np.meshgrid(self.x,self.x)
    def gauss(self,x,y,s=.04):return np.exp(-((self.X-x)**2+(self.Y-y)**2)/(2*s*s))
    def test_single_and_merged(self):
        for p in [self.gauss(0,0),self.gauss(-.02,0,.15)+self.gauss(.02,0,.15)]:
            self.assertIsNone(analyse(p,self.x,[-.2,.2])['pair'])
    def test_two_and_saddle(self):
        p=self.gauss(-.2,0,.1)+self.gauss(.2,0,.1)
        a=analyse(p,self.x,[-.2,.2]);self.assertIsNotNone(a['pair'])
        expected=10*np.log10(min(p[50,30],p[50,70])/p[50,50])
        self.assertAlmostEqual(a['saddle_depth_dB'],expected,places=8)
    def test_plateau_and_boundary(self):
        p=np.ones((101,101));p[40:43,40:43]=3
        a=analyse(p,self.x);interior=[z for z in a['peaks'] if not z['boundary']]
        self.assertEqual(len(interior),1);self.assertEqual(interior[0]['plateau_cells'],9)
        self.assertTrue(analyse(self.gauss(-.5,0),self.x)['peaks'][0]['boundary'])
    def test_sidelobe_control(self):
        p=self.gauss(-.2,0)+.2*self.gauss(.2,0)
        a=analyse(p,self.x,[-.2,.2]);c=analyse(p,self.x)
        self.assertTrue(sidelobe_conflict(a,c,.01))
    def test_no_robust_block_support(self):
        a=analyse(self.gauss(-.2,0)+self.gauss(.2,0),self.x,[-.2,.2])
        b=analyse(self.gauss(0,0,.2),self.x,[-.2,.2])
        self.assertFalse(all(same_pair(a,z) for z in [a,a,b,a,a]))
    def test_saddle_uses_two_dimensions(self):
        # A higher connecting path away from the horizontal cut reduces saddle depth.
        p=self.gauss(-.2,0,.08)+self.gauss(.2,0,.08)
        p+=.3*np.exp(-((self.Y-.1)/.07)**2)*np.exp(-(self.X/.3)**8)
        a=analyse(p,self.x,[-.2,.2]);q=a['pair']
        line_depth=min(z['level_dB'] for z in q)-10*np.log10(np.min(p[50,30:71])/p.max())
        self.assertLess(a['saddle_depth_dB'],line_depth)

if __name__=='__main__':unittest.main()
