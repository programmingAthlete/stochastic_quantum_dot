import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import argparse

parser = argparse.ArgumentParser(description='a')
parser.add_argument('-c', '--compile', type=int, help='compile the C file')
parser.add_argument('-N', '--iterations', type=int, help='number of iterations')
parser.add_argument('-i', '--steps', type=int ,help='number of time steps')
parser.add_argument('-MAXT', '--maxTime', type=float ,help='maximim time fo iteration')
parser.add_argument('-dt', '--dt', type=float ,help='step')
parser.add_argument('-T', '--temperature', type=float ,help='temperature of the reservoir')
parser.add_argument('-mu', '--mu', type=float ,help='chemical poential of the reservoir')
parser.add_argument('-k', '--coef', type=float, help='reservoir coefficient')
parser.add_argument('-t', '--type', type=str,help='type of process')
args = parser.parse_args()

N = 10000
MAXT = 10.0
MAXSTEPS = 1000
dt = MAXT / MAXSTEPS
T = 1.
mu = 1.0
k = 1.0
type = "down"
if args.type:
    type = args.type
if args.temperature:
    T = args.temperature
if args.mu:
    T = args.mu
if args.iterations :
    N = args.iterations
if args.steps and not args.dt and not args.maxTime:
    MAXSTEPS = args.steps
    dt = MAXT/MAXSTEPS
if args.maxTime and args.steps and not args.dt:
    MAXT = args.MAXT
    MAXSTEPS = args.steps
    dt = MAXT/MAXSTEPS
if args.maxTime and args.steps and args.dt:
    print('[-] Arguments not valid')
    sys.exit()
if args.maxTime and args.dt:
    dt = args.dt
    MAXT = args.maxTime
    MAXSTEPS = MAXT / dt
if args.steps and args.dt:
    MAXSTEPS = args.steps
    dt = args.dt
    MAXT = dt * MAXSTEPS
if args.coef:
    k = args.coef
if args.compile == 1:
    y = os.system('gcc simulation.c -o simulation -lm')
    if y != 0:
        sys.exit()
os.system('./simulation %s %s %s %s %s %s %s %s' % (N, MAXT, MAXSTEPS, dt, T, mu, k, type))
data = pd.read_csv('data.txt', sep=',')

time = np.array(data['time'].values)
p = np.array(data['p'].to_numpy())
x = np.array(data['x'].values)
dQ = np.array(data['dQ'].values)
dQ_0 =  np.array(data['dQ_0'].values)
dw = data['dw'].values
dw_0 = np.array(data['dw_0'].values)
e = np.array(data['e'].values)
p_theory = np.array(data['p_theory'].values)
sigma = np.array(data['sigma'].values)
fermi = k / (np.exp(e/T) + 1)

evar = r'$x(t)$'
if type != "up" and type != "down":
    evar = r'$\epsilon(t)$'
fig, axs = plt.subplots(3)
fig.suptitle('N = %d ,  resolution = %d, dt = %.2f, process %s' % (N, MAXSTEPS,dt, type))
axs[0].set_title(r"Evolution of the occupation probability $(T ={:2.1f}, \mu ={:2.1f},  k = {:2.1f}$)".format(T, mu, k))
axs[0].plot(time, p, color="blue", label=r'$p(t)$')
axs[0].plot(time, x, ls='--',color="black", linewidth=0.8,label=r'$n(t)$')
axs[0].set_xlabel(r'$t$')
axs[0].legend(loc='upper right')

axs[1].plot(time, dQ, color='red', label=r'$\langle \dot{Q}(t) \rangle$')
axs[1].plot(time, dw, color='blue', label=r'$\langle \dot{w}(t) \rangle$')
axs[1].set_xlabel(r'$t$')
axs[1].legend(loc='upper right')

axs[2].plot(e, p,linewidth=1,label=r'$p(x(t))$')
axs[2].plot(e, fermi,ls='--',linewidth=1,label=r'$p_{Q-S}(t)$')
axs[2].set_xlabel(r'$x(t)$')
axs[2].legend(loc='upper right')

plt.subplots_adjust(hspace=0.5)

plt.show()
