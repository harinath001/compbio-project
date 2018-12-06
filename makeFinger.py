import numpy as np

def makeFinger(v):
    # Input:  vector of integers, v
    # Output: vector of fingerprints, f where f(i) = |{j: |{k:v(k)=j}|=i }|
    #         i.e. f(i) is the number of elements that occur exactly i times 
    #         in the vector v

    # "v" is a array
    mini = np.min(v)
    maxi = np.max(v)
    h1 = np.histogram(v, range=(mini, maxi), bins=maxi-mini+1)[0]
    f = np.histogram(h1, range=(0, max(h1)), bins=max(h1)+1)[0]
    f = f[1:].reshape(len(f)-1, 1)
    return f