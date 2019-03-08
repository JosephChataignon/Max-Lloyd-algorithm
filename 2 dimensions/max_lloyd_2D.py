
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
def f(t):
    return gaussian(t)

def random_distrib():
    #return [random.uniform(-1,1),random.uniform(-1,1)]
    return [random.gauss(0,1),random.gauss(0,1)]

##no 2d adapted after this
# computes MSE between 2 adjacent decision thresholds
def interval_MSE(x,t1,t2):
    return integrate.quad(lambda t: ((t - x)**2) * f(t), t1, t2)[0]

# computes mean square error on R
def MSE(t,x):
    s = interval_MSE(x[0], -float('Inf'), t[0]) + interval_MSE(x[-1], t[-1], float('Inf'))
    for i in xrange(1,len(x)-1):
        s = s + interval_MSE(x[i], t[i-1], t[i])
    return s

# t1 and t2 are the boundaries of the interval on which the centroid is calculated
def centroid(t1,t2):
    if integrate.quad(f, t1, t2)[0] == 0 or t1 == t2:
        return 0
    else:
        return integrate.quad(lambda t:t*f(t), t1, t2)[0] / integrate.quad(f, t1, t2)[0]

# t is an array containing the initial decision boundaries
# x is an array containing the representation levels
# error_threshold is the threshold to reach for the algorithm to stop
def maxlloyd(t,x,error_threshold):
    e = MSE(t,x)
    error = [e]
    c = 0
    while e > error_threshold and c < 300:
        c = c+1
        if c%2 == 1:
            # adjust boundaries
            for i in xrange(len(t)):
                t[i] = 0.5 * ( x[i] + x[i+1] )
        else:
            # adjust levels
            x[0] = centroid(-float('Inf'), t[0])
            x[-1] = centroid(t[-1], float('Inf'))
            for i in xrange(1,len(x)-1):
                x[i] = centroid(t[i-1], t[i])
        e = MSE(t,x)
        error.append(e)
    return x,t,error


# Test of maxlloyd function
def test_maxlloyd():
    t = [-0.5,0,0.5]
    x = [-1,0,1,1.5]
    x2,t2,error = maxlloyd(t,x,0.1)
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

plot_avg_error(20000)

































