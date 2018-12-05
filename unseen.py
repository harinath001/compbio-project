import numpy as np
from utilities.utilities import find,greater
import math
import math
from scipy.stats import poisson
from scipy.optimize import linprog
import pdb

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
    xLPmin = 1.00/float((k*max(10, k)))
    min_i = min(find(greater(f, 0)))
    if min_i > 0:
        xLPmin = float(min_i)/float(k)
    maxLPIters = 1000


    x = 0
    histx = 0
    fLP = [0]*len(f)
    for i in range(0, len(f)):
        if f[i]>0:
            wind = [ max(0, i-math.ceil(math.sqrt(i))), min(i+math.ceil(math.sqrt(i)), len(f))]
            if sum(f[int(wind[0]):int(wind[1])+1]) < math.sqrt(i):
                x = [x, i/k]
                histx = [histx, f[i]]
                fLP[i] = 0
            else:
                fLP[i] = f[i]

    fmax = max(find(greater(fLP, 0)))
    if not fmax:
        x = x[1:]
        histx = histx[1:]
        return
    LPmass = 1 - np.matmul(np.array([x]), np.array([histx]).T)
    fLP = fLP[0:fmax+1] + [0]*int(math.ceil(math.sqrt(fmax))) # assuming list append here
    szLPf = len(fLP)

    xLPmax = float(fmax)/float(k)
    xLP = xLPmin*(np.array([gridFactor**x for x in range(0, int(math.ceil(math.log(float(xLPmax)/float(xLPmin))/math.log(gridFactor))+1 ))]))
    # convert xLP round to 4 as default in matlab
    xLP = np.array([round(x, 4) for x in xLP])
    szLPx = len(xLP)

    objf = np.zeros(szLPx+2*szLPf).reshape(szLPx+2*szLPf, 1)
    sqrt_fLP = np.array([1.0/math.sqrt(x+1) for x in fLP])
    sqrt_fLP = sqrt_fLP.reshape(len(sqrt_fLP), 1)

    # objf[szLPx:] = sqrt_fLP*len(objf[szLPx:]) # converting to
    objf[szLPx::2] = sqrt_fLP
    objf[szLPx+1::2] = sqrt_fLP

    A = np.zeros(2*szLPf*(szLPx+2*szLPf)).reshape(2*szLPf, szLPx+2*szLPf)
    b = np.zeros(2*szLPf).reshape(2*szLPf, 1)

    for i in range(0, szLPf):
        A[2*i][:szLPx] = poisson.pmf(i+1, k*xLP) # [pmf]*times
        A[2*i+1][:szLPx] = (-1)*A[2*i][:szLPx]
        A[2*i][szLPx+2*i] = -1
        A[2*i+1][szLPx+2*i+1] = -1
        b[2*i] = fLP[i]
        b[2*i+1] = 0-fLP[i]

    Aeq = np.zeros(szLPx+2*szLPf).reshape(1, szLPx+2*szLPf)
    Aeq[0][:szLPx] = xLP.reshape(1, len(xLP)) # Aeq is row matrix so replace the first row first szLPx elements to xLP elements
    beq = LPmass

    for i in range(szLPx):
        A[:,i] /= xLP[i]
        Aeq[0][i] /= xLP[i]

    options = {"maxiter": maxLPIters, "disp": False}
    # pdb.set_trace()
    result1 = linprog(np.squeeze(objf), A, b, Aeq, beq, options=options)
    exitflag = result1["status"]
    fval = result1["fun"]
    if exitflag == 1:
        print 'maximum number of iterations reached--try increasing maxLPIters'
    elif exitflag > 1:
        print 'LP1 solution was not found, still solving LP2 anyway...', exitflag



    #106
    objf2 = np.zeros(objf.shape)
    objf2[0:szLPx] = 1
    A2 = np.vstack((A, objf.T))
    b2 = np.vstack((b,fval+alpha))
    for i in range(szLPx):
        objf2[i] /= xLP[i]

    #116 sol2,fval2 = linprog?
    # [sol2, fval2, exitflag2, output] = linprog(objf2, A2, b2, Aeq, beq, zeros(szLPx+2*szLPf,1), Inf*ones(szLPx+2*szLPf,1),[], options);
    result2 = linprog(np.squeeze(objf2), A2, b2, Aeq, beq, options=options)
    fval2 = result2["fun"]
    exitflag2 = result2["status"]
    sol2 = result2["x"]

    if exitflag2 != 0:
        print "LP2 solution was not found", exitflag2

    sol2[0:szLPx] = np.divide(sol2[0:szLPx], xLP)
    x = np.hstack((x,xLP))
    histx = np.hstack((histx,sol2.T))
    x.sort()
    ind = x.argsort()
    histx = histx[ind]
    ind = np.where(histx > 0)
    x = x[ind]
    histx = histx[ind]
    return [histx, x]

n = np.array([[8], [1]])
temp = unseen(n)
print(temp)
print(sum(temp[0]*temp[1]))