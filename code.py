import numpy as np
import math
from makeFinger import makeFinger
from unseen_1 import unseen
from entropy_1 import entropy_est
#Generate a sample of size 10,000 from the uniform distribution of support 100,000
n=100000
k=10000
n=1000
k=250
sample =  np.random.random_integers(low=1,high=n,size=k)
print "Sample:"
print sample

f = makeFinger(sample)
print "Fingerprint:"
print f

histx,x = unseen(f)
print "H(x) and x:"
print histx,x

print "True Entropy", math.log(n)

temp = np.array([float(i) for i in range(1,len(f)+1)])
temp /= k
print "Emperical Entropy", -np.dot(f.flatten(),[z*math.log(z) for z in temp])

print "Entropy of Recovered Histogram", -np.dot(histx,[z*math.log(z) for z in x])

print "Entropy using Code", entropy_est(f)