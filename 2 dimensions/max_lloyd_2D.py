
# Max-Lloyd algorithm for finding the optimal quantizer
# in dimension 2

import math
import random
import scipy
import matplotlib.pyplot as plt
from scipy import integrate



def uniform(x,y):
    if x<=1 and x>=-1 and y<=1 and y>=-1:
        return 0.25
    else:
        return 0

def gaussian(x,y):
    return math.exp(-x**2/2-y**2/2)

# function studied
def f(x,y):
    return gaussian(x,y)

# same as f, but the returned value is zero if (x,y) is out of the boundaries
def f_with_constraints(x,y,boundaries):
    for b in boundaries:
        if b[0]*x + b[1]*y > b[2] :
            return 0
    return f(x,y)

# distribution studied
def random_distrib():
    #return [random.uniform(-1,1),random.uniform(-1,1)]
    return [random.gauss(0,1),random.gauss(0,1)]

## not 2-d adapted after this
# computes MSE inside one region
def region_MSE(x,germs,boundaries):
    return integrate.dblquad(lambda x,y: ((t - x)**2) * f(t), t1, t2)[0]# to complete

# computes mean square error on R
def MSE(germs,boundaries):
    s = 0
    for i in xrange(len(germs)):
        s = s + interval_MSE(x[i], t[i-1], t[i])
    return s/len(germs)

# boundaries is an array containing the constraints on the function
def centroid(boundaries):
    if integrate.quad(f, t1, t2)[0] == 0 or t1 == t2:#to complete: if  area == 0
        return 0
    else:
        return integrate.quad(lambda x,y:t*f(t), t1, t2)[0] / integrate.quad(f, t1, t2)[0]#to complete

def adjust_boundaries(germs,k):
    bounds = []
    main_germ = germs[i]
    for i in xrange(germs[i]):
        bounds.append() # to complete
    return bounds

# boundaries is an array containing the initial decision boundaries (array of n-1 constraints)
# constraints' model is [a1,a2,b] parameters of the equation a1*x + a2*y <= b
# germs is an array containing the germ for each region, on the model [x,y]
# error_threshold is the threshold to reach for the algorithm to stop
def maxlloyd(boundaries,germs,error_threshold):
    e = MSE(boundaries,germs)
    error = [e] # store the evolution of MSE through the loop
    c = 0 # counts the number of executions of the loop
    while e > error_threshold and c < 300:
        c = c+1
        if c%2 == 1:
            # adjust boundaries
            for i in xrange(len(boundaries)):
                boundaries[i] = adjust_boundaries(germs,i)
        else:
            # adjust germs
            for i in xrange(len(germs)):
                germs[i] = centroid(boundaries[i])
        e = MSE(boundaries,germs)
        error.append(e)
    return germs,boundaries,error


# Test of maxlloyd function
def test_maxlloyd():
    b = [-0.5,0,0.5]
    g = [[0,0],[1,1],[-1,-1]]
    x2,t2,error = maxlloyd(b,l,0.1)
    print x2,t2
    plt.plot(error)
    plt.show()
    return x2,t2

# test_maxlloyd()

def estimate(x,t,value):
    for i in xrange(len(t)):
        if t[i] > value:
            return x[i]
    return x[-1]

# Plot of average error
def plot_avg_error(N):
    x,t = test_maxlloyd()
    avg_E = []
    realizations = []
    square_error = []
    for i in xrange(N):
        realizations.append(random_distrib())
        square_error.append((realizations[-1] - estimate(x,t,realizations[-1]))**2)
        avg_E.append(sum(square_error)/len(square_error))
    plt.figure(2)
    plt.plot(avg_E)
    plt.show()

#plot_avg_error(20000)

































