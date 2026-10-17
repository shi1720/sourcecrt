import copy
import itertools as it
import unittest
from sourcecrt.core import *
from sourcecrt.verify import verify


def brute(ms, rs, bounds, supports):
    return {x for x in it.product(*(range(l,u+1) for l,u in bounds))
            if any(all(all(x[j] % m == r[j] % m for j in range(len(bounds)))
                       for i,(m,r) in enumerate(zip(ms,rs)) if i not in F) for F in supports)}


class ExactTests(unittest.TestCase):
    def test_exhaustive(self):
        self.checked = 0
        for ms in ([2,3,5], [3,4,5]):
            supports = grouped_supports([{0},{1},{2}], 1)
            bounds = [(-2,3),(-1,2)]
            for word in it.product(*(list(it.product(range(m),repeat=2)) for m in ms)):
                actual = brute(ms,word,bounds,supports)
                status,x = classify(analyze(ms,word,bounds,supports))
                expected = 'inconsistent' if not actual else 'unique' if len(actual)==1 else 'ambiguous'
                self.assertEqual(status,expected)
                if x is not None:
                    self.assertEqual(actual,{x})
                    c=certificate(ms,word,bounds,supports)
                    self.assertTrue(verify(ms,word,bounds,supports,c))
                self.checked += 1
        print('Exhaustive received arrays:',self.checked)

    def test_joint_separation_family(self):
        for k in range(1,101):
            ms=[5,5*k+2,5*k+3]
            rs=[[3,2],[0,0],[0,0]]
            b=[(0,5*k+3)]*2
            F=grouped_supports([{0},{1},{2}],1)
            self.assertEqual(classify(analyze(ms,rs,b,F)),('unique',(0,0)))
            self.assertEqual(scalar_status(ms,rs,b,F),['ambiguous','ambiguous'])
            self.assertEqual(separation_modulus(ms,F),5)
            self.assertTrue(verify(ms,rs,b,F,certificate(ms,rs,b,F)))

    def test_mutations(self):
        ms=[5,7,11];rs=[[1,2],[0,0],[0,0]];b=[(0,11)]*2;F=[(0,),(1,),(2,)]
        c=certificate(ms,rs,b,F)
        self.assertTrue(verify(ms,rs,b,F,c))
        for j in range(2):
            z=copy.deepcopy(c);z['candidate'][j]+=1
            self.assertFalse(verify(ms,rs,b,F,z))
        z=copy.deepcopy(c);z['records'].pop()
        self.assertFalse(verify(ms,rs,b,F,z))
        z=copy.deepcopy(c);z['records'][1]['exclude_coordinate']=0
        self.assertFalse(verify(ms,rs,b,F,z))

    def test_overlap_and_threshold(self):
        ms=[2,3,5,7]
        # An overlapping source model, not a partition of prime positions.
        F=grouped_supports([{0,1},{1,2},{3}],1)
        S=separation_modulus(ms,F)
        b=[(-3,-3+S)]
        found=False
        for A in F:
            for B in F:
                good=[i for i in range(len(ms)) if i not in set(A)|set(B)]
                if prod(ms[i] for i in good)!=S:continue
                x,y=-3,-3+S
                r=[[y%m if i in A else x%m] for i,m in enumerate(ms)]
                self.assertTrue({(x,),(y,)} <= brute(ms,r,b,F));found=True
        self.assertTrue(found)
        # Every received word at width S-1 is injective.
        for rs0 in it.product(*(range(m) for m in ms)):
            pts=brute(ms,[[a] for a in rs0],[(0,S-1)],F)
            self.assertLessEqual(len(pts),1)

    def test_invalid(self):
        with self.assertRaises(ValueError):analyze([4,6],[[0],[0]],[(0,1)],[()])
        self.assertEqual(classify(analyze([],[],[(0,1)],[()])),('ambiguous',None))
        self.assertEqual(classify(analyze([],[],[(1,1)],[()])),('unique',(1,)))

if __name__=='__main__':unittest.main()
