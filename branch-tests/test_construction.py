import simple_pbox as spb
import numpy as np
import pandas as pd
import scipy.stats as sts

import matplotlib.pyplot as plt

#%% Intro
'''This file demonstrates the contruction of probability boxes from the pbox
class.'''

#%%
'''P-boxes are a generalisation of a deterministic number. The Pbox class
supports the construction of p-box-like scalars numbers.'''

R = pb.pbox(2) #This should return a scalar number - 2
print(f'R = {R} with width, w(R) = {R.width}')
R.plot();

#%% #The next four should return an interval
'''P-boxes also generalise intervals...

All four methods below can be used to construct an interval between 1 and 2 ->
two arguments - one for the left side and one for the right side and a single
array-like argument (list, array, or tuple).
'''

I = []
I.append(pb.pbox(1,2))
I.append(pb.pbox([1,2]))
I.append(pb.pbox(np.array([1,2])))
I.append(pb.pbox((1,2)))

for i in I:
    print(f'I = {i}')
    i.plot()
    plt.show()
    
#%% Comparison to intervals
'''To save on computation and to easily enforce information about epistemic 
uncertainty use the Interval class. The Pbox class will return an interval-
looking p-box object, whereas the Interval class will return a mathematical 
interval (a pair of bounds)'''

I_ival = pb.ival.I(1,2)

print(f'Interval-looking p-box {I[0]}')
print(f'Interval-looking p-box mean: {I[0].mean} and variance: {I[0].var}')
print(f'Mathametical interval: {I_ival}')
print(f'Mathametical interval mean: {I_ival.mean} and variance: {I_ival.var}') #This returns an error!

#%% Distributions
'''P-boxes also generalise distributions...
The two methods below can be used to construct a distribution from an array-like of
data points.

You can control the weight of each data point (step) in the EDF which will be
reflected in both the visual appearence of the uncertain number, but also in
its moments.
'''

D = []
D.append(pb.pbox([1,2,3,4])) #This should return a crisp distribution with equiprobable steps
D.append(pb.pbox([1,2,3,4], p_left=[0.05,0.35,0.82])) #This should return a a crisp distribution with step heights specified by p_left
D.append(pb.pbox([1,2,3,4], p_left=[1/17,1/3,1/1.5])) #This should return a a crisp distribution with step heights specified by p_left

for d in D:
    print(f'D = {d}')
    d.plot()
    plt.show()
    
#%% The next six should return a pbox with equiprobable steps
'''Finally the previous methods can be combined to construct a data-based p-box
with equiprobable steps. Any of the seven (nine) methods below will work.'''
# It is messy to work with imprecision in both x and p
import Interval as ival

Pe = []
Pe.append(pb.pbox([1,2,3,4], [5,6,7,8])) #Specify bounds as lists
Pe.append(pb.pbox(np.array([1,2,3,4]), np.array([5,6,7,8]))) #Specify bounds as arrays
Pe.append(pb.pbox([[1,5],[2,6],[3,7],[4,8]])) #Specify focal elements as a list of lists
Pe.append(pb.pbox(np.array([[1,5],[2,6],[3,7],[4,8]]))) #Specify focal elements as an array of arrays
Pe.append(pb.pbox([(1,5),(2,6),(3,7),(4,8)])) #Specify focal elements as a list of tuples
Pe.append(pb.pbox(np.array([(1,5),(2,6),(3,7),(4,8)]))) #Specify focal elements as an array of arrays, constructed with tuples
Pe.append(pb.pbox(np.array([
                            ival.I(1,5),
                            ival.I(2,6),
                            ival.I(3,7),
                            ival.I(4,8)]))) #Specify focal elements as an array of Interval objects
                                            #This can obviously be expanded to a list too

for pe in Pe:
    print(f'P = {pe}')
    pe.plot()
    plt.show()

#%% The next six should return a pbox with step heights, for both bounds, specified by p_left
'''Just like for distribution, any of the seven (nine) methods below will also
produce p-boxes with specified weights for the focal elements.'''

Pl = []
Pl.append(pb.pbox([1,2,3,4], [5,6,7,8], p_left=[0.05,0.35,0.82]))
Pl.append(pb.pbox(np.array([1,2,3,4]), np.array([5,6,7,8]), p_left=[0.05,0.35,0.82]) )
Pl.append(pb.pbox([[1,5],[2,6],[3,7],[4,8]], p_left=[0.05,0.35,0.82]))
Pl.append(pb.pbox(np.array([[1,5],[2,6],[3,7],[4,8]]), p_left=[0.05,0.35,0.82]))
Pl.append(pb.pbox([(1,5),(2,6),(3,7),(4,8)], p_left=[0.05,0.35,0.82]))
Pl.append(pb.pbox(np.array([(1,5),(2,6),(3,7),(4,8)]), p_left=[0.05,0.35,0.82]))
Pl.append(pb.pbox(np.array([
                            ival.I(1,5),
                            ival.I(2,6),
                            ival.I(3,7),
                            ival.I(4,8)]), p_left=[0.05,0.35,0.82]))

for pl in Pl:
    print(f'P = {pl}')
    pl.plot()
    plt.show()
    
#%% The next six should return a pbox with step heights specified by individually for the left and right bounds
'''Similar logic can be used to construct p-boxes with different weights for 
the left and right bounds. The sizes of probabilities (n) and quantiles (n+1)
must correspond.'''

Plr = []
Plr.append(pb.pbox([1,2,3,4], [5,6,7,8], p_left=[0.05,0.35,0.82], p_right=[0.1,0.25,0.95]))
Plr.append(pb.pbox(np.array([1,2,3,4]), np.array([5,6,7,8]), p_left=[0.05,0.35,0.82], p_right=[0.1,0.25,0.95]) )
Plr.append(pb.pbox([[1,5],[2,6],[3,7],[4,8]], p_left=[0.05,0.35,0.82], p_right=[0.1,0.25,0.95]))
Plr.append(pb.pbox(np.array([[1,5],[2,6],[3,7],[4,8]]), p_left=[0.05,0.35,0.82], p_right=[0.1,0.25,0.95]))
Plr.append(pb.pbox([(1,5),(2,6),(3,7),(4,8)], p_left=[0.05,0.35,0.82], p_right=[0.1,0.25,0.95]))
Plr.append(pb.pbox(np.array([(1,5),(2,6),(3,7),(4,8)]), p_left=[0.05,0.35,0.82], p_right=[0.1,0.25,0.95]))
Plr.append(pb.pbox(np.array([
                            ival.I(1,5),
                            ival.I(2,6),
                            ival.I(3,7),
                            ival.I(4,8)]), p_left=[0.05,0.35,0.82],
                                           p_right=[0.1,0.25,0.95]))

for plr in Plr:
    print(f'P = {plr}')
    plr.plot()
    plt.show()

#%% Discretisation
'''By default, when constructing a p-box from data, the Pbox class will create
objects with max(num_data_points_left, num_data_points_right, 100) steps for
each bound.
The default number of steps for construction (100) can be changed by passing 
the desired value to the n_step argument of the pbox constructor.

Attention must be paid if an n_step different to the default one is set. If the
number of data points is smaller than n_step, but n_step is not divisible by
the number of pointsthe result may and most
likely will be different to what is expexted (and thus probably incorrect).
Here is an example using the default, indivisible n_step and divisible n_step. 
Compare the plots - black is hidden by blue.'''

pb_ns = pb.pbox([1,2], [3,4,5])
pb_ns.plot()
plt.title('P-box with differently sized bounds', fontsize=16)
plt.show()


_, ax = plt.subplots()

pb_ns = pb.pbox([1,2,3,4], [5,6,7,8])
pb_ns.plot(ax=ax)

pb_ns = pb.pbox([1,2,3,4], [5,6,7,8], n_step=10)
pb_ns.plot(ax=ax, c='r')

pb_ns = pb.pbox([1,2,3,4], [5,6,7,8], n_step=12)
pb_ns.plot(ax=ax, c='b')

ax.grid(c=[0.9]*3)
ax.set_axisbelow(True)
ax.set_title('Incorrect p-boxes, when $n_{step}$ is incorrect', fontsize=16)

#%%
'''Similarly, if the probabilities supplied to p_left/p_right do not return an 
integer value when multiplied by n_step the result will be incorrect. Observe
that none of the coloured p-boxes matches the correct one (black).'''
_, ax = plt.subplots()

pb_ns = pb.pbox([1,2,3,4], [5,6,7,8], p_left=[0.05,0.35,0.82],
                p_right=[0.1,0.25,0.95])
pb_ns.plot(ax=ax)

pb_ns = pb.pbox([1,2,3,4], [5,6,7,8],  p_left=[0.05,0.35,0.82],
                p_right=[0.1,0.25,0.95], n_step=10)
pb_ns.plot(ax=ax, c='r')

pb_ns = pb.pbox([1,2,3,4], [5,6,7,8], p_left=[0.05,0.35,0.82],
                p_right=[0.1,0.25,0.95], n_step=12)
pb_ns.plot(ax=ax, c='b')

ax.grid(c=[0.9]*3)
ax.set_axisbelow(True)
ax.set_title('Incorrect p-boxes, when $n_{step}$ is incorrect', fontsize=16)


#%% Condensation

#%%
pb = reload(pb)

# pb_vp_100 = pb.pbox(np.linspace(1,2,100), np.linspace(2,3,100))
# pb_vp_1000 = pb.pbox(np.linspace(1,2,1000), np.linspace(2,3,1000))
pb_vp_100 = pb.pbox(1+np.random.rand(100,1), 2+np.random.rand(100,1))
pb_vp_1000 = pb.pbox(1+np.random.rand(1000,1), 2+np.random.rand(1000,1))


_, ax = plt.subplots()
pb_vp_100.plot(ax=ax, c=['r', 'k'])
pb_vp_1000.plot(ax=ax, c=['b', 'g'])
#%% Test the different moment constructors

#%% There is no capability to mix both - get from pbox.py

#%% Interval
# I = pb.pbox([1,1],[2,2],[0,1],[0,1]) #This is the native format and it will fail
I = pb.pbox([1,1],[2,2],[1],[1]) #But this works - check if first element of probability is 0
I.plot()

I = pb.pbox([1,1,1],[2,2,2],[0,1],[0,1]) #This also works - duplicate first element of qunatiles
I.plot()

#%% Distribution - e.g. Juan
# x = [-0.131, -0.032, 0.0917, 0.156, 0.215, 0.277, 0.339, 0.463, 0.587, 0.612] #This does not work
x = [-0.131, -0.131, -0.032, 0.0917, 0.156, 0.215, 0.277, 0.339, 0.463, 0.587, 0.612] #But this does
p = [0.001235245, 0.021613397, 0.158158633, 0.308400096, 0.5, 0.691623391, 0.841507078, 0.977440663, 0.998877725, 0.999543196]

D = pb.pbox(x, p_left=p) 
D.plot()

#%% P-box - e.g. DLR's as it is huge
data = pd.read_excel('../AIAAUQDGchallengeProblem2DataSubmission-v7-DLR_Stradtner.xlsx',
                      sheet_name='Alpha 0, Flap 0-CL', header=None)
total = data.loc[8:, 0:3].dropna(how='all').values.astype(np.float64)

#%%
# P = pb.pbox(total[:,0], total[:,2], total[:,1], total[:,3]) #This does not work
P = pb.pbox(total[:,0], total[:,2], total[1:,1], total[1:,3], verbatim_prob=True) #This does not work
P.plot()

#%% Try p-box condensation
p_large = np.linspace(0.001,0.999,10000)
p_alarge = np.arange(0.001, 1, 1/10000)

x = sts.norm(0,1).ppf(p_large)
xa = sts.norm(0,1).ppf(p_alarge)
#%%
import simple_pbox_dev as pb_rt

#%% Test generating p-boxes with or without external probabilities
pb_large_with_p = pb_rt.pbox(x, p_left=p_large[1:])
pb_large_wo_p = pb_rt.pbox(xa)

res = sts.ecdf(xa)
#%%
_, ax = plt.subplots(2,1,sharex=True, sharey=True)
ax[0].stairs(p_large[1:], x, baseline=None)
# ax[1].stairs(p_alarge[1:], xa, baseline=None)
res.cdf.plot(ax=ax[1])
pb_large_with_p.plot(ax=ax[0]) #This produces the correct behaviour, with the right bound wider than needed in some places
pb_large_wo_p.plot(ax=ax[1]) #This is incorrect - multiple crossings and inflations

ax[0].set_title('Probabilities passed')
ax[1].set_title('Probabilities generated')
