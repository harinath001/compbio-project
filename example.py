import numpy as np
import math
from makeFinger import makeFinger
from unseen import unseen
from entropy_estC import entropy_est
#Generate a sample of size 10,000 from the uniform distribution of support 100,000
n=100000
k=10000

sample =  np.random.random_integers(low=1,high=n,size=k)
# print "Sample:"
# print sample

f = makeFinger(sample)
print("Fingerprint:")
print f

histx,x = unseen(f)

print "The probability and Histogram of the probability"
for i in range(len(x)):
	print x[i], histx[i]


print("True Entropy", math.log(n))

temp = np.array([float(i) for i in range(1,len(f)+1)])
temp /= k
print("Emperical Entropy", -np.dot(f.flatten(),[z*math.log(z) for z in temp]))

print("Entropy of Recovered Histogram", -np.dot(histx,[z*math.log(z) for z in x]))

print("Entropy using Code", entropy_est(f))