def f_quadratic(x, a, b, c):
    return a * x * x + b * x + c

from scipy.optimize import curve_fit
import numpy

data = numpy.loadtxt("ofs_experiment_6.txt")
print(data)
x = data[:,0]
y = data[:,1]

input_params = [0, 0, 0]
output_params, covar = curve_fit(f_quadratic, x, y, p0=input_params)
print(output_params)

import matplotlib.pyplot as plt

a = output_params[0]
b = output_params[1]
c = output_params[2]
ymodel = f_quadratic(x, a, b, c)

# origional parameters
a2 = 238.44006955
b2 = -208.4151086
c2 = 379.22456916
# converted parameters
alpha = a2
x0 = b2 / (-2.0 * alpha)
beta = c - alpha * numpy.power(x0, 2.0)
print(alpha, x0, beta)
# adjusted parameters
x0 = 0.4
beta = 333
# converted back parameters
a2 = alpha
b2 = -2.0 * alpha * x0
c2 = alpha * numpy.power(x0, 2.0) + beta

ymodel2 = f_quadratic(x, a2, b2, c2)

plt.plot(x, y, 'g')
plt.plot(x, ymodel, 'b')
plt.plot(x, ymodel2, 'r')
plt.savefig('fig.png')
