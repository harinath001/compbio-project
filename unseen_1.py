import numpy as np
import math
from scipy.stats import poisson
from scipy.optimize import linprog

def unseen(f):
	#Calculation of Sample Size
	f = np.transpose(f).flatten()
	k = np.array(range(1,len(f)+1))
	k = np.dot(f,k)

	#Variable definitions
	gridFactor = 1.05
	alpha = 0.5
	xLPmin = 1.0/float((k*max(10,k)))
	min_i = min([i for i in range(len(f)) if f[i] > 0])
	if min_i > 0:
		xLPmin = float(min_i/(k*1.0))
	maxLPIters = 1000

	x= list()
	histx = list()
	fLP = np.zeros(len(f))

	for i in range(len(f)):
		if f[i] > 0:
			window = [max(0,i-math.ceil(math.sqrt(i))),min(i+math.ceil(math.sqrt(i)),len(f))]
			window = [int(window[0]),int(window[1])]

			if np.sum(f[window[0]:window[1]]) < math.sqrt(i):
				x.append(float(i/(k*1.0)))
				histx.append(float(f[i]))
				fLP[i] = 0
			else:
				fLP[i] = f[i]

	fmax = [i for i in range(len(fLP)) if fLP[i] > 0]
	if len(fmax) == 0:
		return x,histx
	else:
		fmax = max(fmax)+1

	LPmass = 1.0 - np.dot(x,histx)

	fLP = np.concatenate((fLP[:fmax],np.zeros(int(math.ceil(math.sqrt(fmax))))))
	szLPf = len(fLP)

	xLPmax = float(fmax/(k*1.0))
	temp = [gridFactor**i for i in range(0,int(math.ceil(math.log((xLPmax*1.0)/xLPmin)/math.log(gridFactor))+1))]
	xLP = np.array(temp)
	xLP = xLPmin*xLP
	szLPx = len(xLP)

	objf = np.zeros((szLPx + 2*szLPf,1))
	for i in range(szLPx,len(objf),2):
		objf[i] = 1.0/math.sqrt(fLP[(i-szLPx)/2]+1)
		objf[i+1] = 1.0/math.sqrt(fLP[(i-szLPx)/2]+1)

	A = np.zeros((2*szLPf,szLPx+2*szLPf))
	b = np.zeros((2*szLPf,1))

	for i in range(szLPf):
		A[2*i][:szLPx] = poisson.pmf(i+1,k*xLP)
		A[(2*i)+1][:szLPx] = (-1)*A[2*i][:szLPx]
		A[2*i][szLPx+(2*i)] = -1
		A[(2*i)+1][szLPx+(2*i)+1] = -1
		b[2*i] = fLP[i]
		b[(2*i)+1] = (-1)*fLP[i]

	Aeq = np.zeros((1,szLPx+2*szLPf))
	#Aeq[0:szLPx][0] = xLP
	for i in range(szLPx):
		Aeq[0][i] = xLP[i]
	beq = np.zeros((1,1))
	beq[0][0] = LPmass

	for i in range(szLPx):
		Aeq[0][i] = Aeq[0][i]/xLP[i]
		A[:,i] = A[:,i]/xLP[i]

	# print objf.shape
	# print A.shape
	# print b.shape
	# print Aeq.shape
	# print beq.shape
	# print szLPx,szLPf
	options = {"maxiter": maxLPIters, "disp": False}
	
	A = np.array([[round(temp,4) for temp in y] for y in A])
	result1 = linprog(np.squeeze(objf), A, b, Aeq, beq, options=options)

	exitflag = result1["status"]
	fval = result1["fun"]
	if exitflag == 1:
		print('maximum number of iterations reached--try increasing maxLPIters')
	elif exitflag > 1:
		print('LP1 solution was not found, still solving LP2 anyway...', exitflag)
	
	print(result1['message'])

	objf2 = np.zeros((szLPx + 2*szLPf,1))
	objf2[:szLPx] = 1
	A2 = np.vstack((A,np.transpose(objf)))
	b2 = np.vstack((b,np.array(fval+alpha)))

	for i in range(szLPx):
		objf2[i] = objf2[i]/xLP[i]

	result2 = linprog(np.squeeze(objf2), A2, b2, Aeq, beq, options=options)

	exitflag = result2["status"]
	if exitflag > 1:
		print('LP2 solution was not found', exitflag)

	print(result2['message'])

	sol2 = result2['x']
	print(sol2)
	sol2[0:szLPx] = np.divide(sol2[0:szLPx],xLP)

	x = np.concatenate((x,xLP))
	histx = np.concatenate((histx,sol2))

	index = np.argsort(x)
	x = np.sort(x)
	histx = histx[index]

	index = np.where(histx > 0)
	x = x[index]
	histx = histx[index]
	return histx,x