import numpy as np
from utilities.utilities import find,greater
import math
import math.ceil as ceil
import math.sqrt as sqrt


def unseen(f):
    # Input: fingerprint f, where f(i) represents number of elements that
    # appear i times in a sample.  Thus sum_i i*f(i) = sample size.
    # File makeFinger.m transforms a sample into the associated fingerprint.
    #
    # Output: approximation of 'histogram' of true distribution.  Specifically,
    # histx(i) represents the number of domain elements that occur with
    # probability x(i).   Thus sum_i x(i)*histx(i) = 1, as distributions have
    # total probability mass 1.
    #
    # An approximation of the entropy of the true distribution can be computed
    # as:    Entropy = (-1)*sum(histx.*x.*log(x))

    f = f.T.flatten()
    k = np.matmul(f, np.array(range(1, f.shape[0]+1)).T)

    gridFactor = 1.05
    alpha = 0.5
    xLPmin = 1/(k*max(10, k))
    min_i = min(find(greater(f, 0)))
    if min_i > 0:
        xLPmin = min_i/k
    maxLPIters = 1000


    x = 0
    histx = 0
    fLP = [0]*len(f)
    for i in range(0, len(f)):
        if f[i]>0:
            wind = [ max(0, i-math.ceil(math.sqrt(i))), min(i+math.ceil(math.sqrt(i)), len(f))]
            if sum(f[wind[0]:wind[1]+1]) < sqrt(i):
                x = [x, i/k]
                histx = [histx, f[i]]
                fLP[i] = 0
            else:
                fLP[i] = f[i]

    fmax = max(find(greater(fLP, 0)))
    if

n = np.array([[1,2,3], [4,5,6]])
unseen(n)