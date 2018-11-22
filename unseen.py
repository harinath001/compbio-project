import numpy as np
from utilities.utilities import find



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

    print("k ", k)
    gridFactor = 1.05
    alpha = 0.5
    xLPmin = 1/(k*max(10, k))
    # min_i = min()

n = np.array([[1,2,3], [4,5,6]])
unseen(n)