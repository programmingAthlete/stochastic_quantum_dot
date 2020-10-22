import numpy as np
import matplotlib.pyplot as plt
import random
import sys

'''
    Reservoir Class
'''
class Reservoir:
    '''
        params: float:  T: reservoir's temperature
        params: float:  mu: reservoir's chemical potential
        params: float:  C: reservoir's-dot tunneling constant
        params: float:  e: energy of the dot in contact with the
                            reservoir
        return: Reservoir_Object
    '''
    def __init__(self, T,mu, C,e):
        self.T = T
        self.mu = mu
        self.C = C
        self.e = e
        self.up, self.down = self.rates(e)

    def fermi(self, e):
        return self.C /(1. + np.exp((e-self.mu)/self.T))

    def rates(self, e):
        return self.fermi(e), 1-self.fermi(e)


'''
    Calculates the miscrostate trajectories

    params: float:            MAXTIME:  maximum time of the
                                            simulation
    params  np_array:         times:    array of discretised time
    params: Reservoir_Object: res1:     reservoir 1
    params: Reservoir_Object: res2:     reservoir 2
    params: Reservoir_Object: res3:     reservoir 3
    reurn: tuple: (p_t, I, Q, jumps_t): occupation probability
                                          (2-D list),
                                         Currents (3-D list),
                                         Heats (3-D list),
                                         ransition times
                                          (MAXTIME-D list)
'''
def trajectory(MAXTIME, times, res1, res2, res3):
    t = 0
    p = [0,0]
    resolution = len(times)
    p_t = [[0,0] for _ in range(resolution)]
    jumps_t = [0]
    I_tau = [0,0,0]
    Q_tau = [0,0,0]
    I = {1 : [0], 2: [0], 3:[0]}
    Q = {1 : [0], 2: [0], 3:[0]}
    i = 0
    while i < resolution:
        pold = p
        while t <= times[i]:
            rate_uu = (1-p[0])* res1.up
            rate_ud = p[0] * res1.down
            rate_u = rate_uu + rate_ud
            rate_du = (1 - p[1]) * (res2.up + res3.up)
            rate_dd = p[1] * (res2.down + res3.down)
            rate_d = rate_du + rate_dd
            tau_u = - (1.0 / rate_u) * np.log(random.random())
            tau_d = - (1.0 / rate_d) * np.log(random.random())
            r = random.random()
            if tau_u < tau_d:
                t += tau_u
                if p[0] == 0:
                    I_tau = [1,0,0]
                    Q_tau = [res1.e - res1.mu, 0, 0]
                else:
                    I_tau = [-1,0,0]
                    Q_tau = [-res1.e + res1.mu, 0, 0]
                p[0] = 1- p[0]
            else:
                t += tau_d
                if p[1] == 0:
                    if r < res2.up / rate_du:
                        I_tau = [0,1,0]
                        Q_tau = [0, res2.e - res2.mu, 0]
                    else:
                        I_tau = [0,0,1]
                        Q_tau = [0, 0, res3.e - res3.mu]
                else:
                    if r < res3.down / rate_dd:
                        I_tau = [0,0,-1]
                        Q_tau = [0, 0, - res3.e + res3.mu]
                    else:
                        I_tau = [0,-1,0]
                        Q_tau = [0, -res2.e + res2.mu, 0]
                p[1] = 1 - p[1]
            jumps_t.append(t)
            for k in range(1,4):
                I[k].append(I_tau[k-1])
                Q[k].append(Q_tau[k-1])


        time_idx = np.searchsorted(times > t, True)
        corr_idx = min(resolution - 1, time_idx)
        for j in range(i, corr_idx - 1):
            p_t[j] = pold
        i = time_idx
    return p_t, I, Q, jumps_t


N = 500
MAXTIME = 3000
resolution =120000
e = 2.5
## Reservoirs
# Setup in such a way that
#res1.up < res2.up < res3.up <=> res1.down > res2.down > res3.down
res1 = Reservoir(T=0.9,mu=1.0,C=1.0,e=2.5)
res2 = Reservoir(T=1.0,mu=1.0,C=1.0, e=1.5)
res3 = Reservoir(T=5.2,mu=1.0,C=1.0, e=1.5)

print('res1.up = %s\nres1.down = %s' % (res1.up, res1.down))
print('\nres2.up = %s\nres2.down = %s' % (res2.up, res2.down))
print('\nres3.up = %s\nres3.down = %s' % (res3.up, res3.down))
if not res1.up < res2.up < res3.up and not res1.T == res2.T == res3.T:
    print('\nthe rates do not satisfy the conditions')
    sys.exit()


## Initiation of the variables
p_avg = [[0,0] for _ in range(resolution)]
jumped_t = []
times = np.linspace(0,MAXTIME, resolution)
sample_rate = resolution / MAXTIME
I = {1 : [0] , 2: [0], 3: [0]}
Q = {1 : [0] , 2: [0], 3: [0]}

dist_flag = False
for i in range(N):
    p_t, I_t, Q_t, jumps_t = trajectory(MAXTIME, times, res1, res2, res3)

    if dist_flag == False:
        jumped_t = jumps_t
        # Empirical comulative fucntion
        for i in range(1,len(I_t[1])):
            Q[1].append(sum(Q_t[1][:i]) / jumps_t[i])
            Q[2].append(sum(Q_t[2][:i]) / jumps_t[i])
            Q[3].append(sum(Q_t[3][:i]) / jumps_t[i])
            I[1].append(sum(I_t[1][:i]) / jumps_t[i])
            I[2].append(sum(I_t[2][:i]) / jumps_t[i])
            I[3].append(sum(I_t[3][:i]) / jumps_t[i])
        dist_flag = True
    ## MEAN VAlues
    for t in range(resolution):
        p_avg[t][0] += p_t[t][0] / N
        p_avg[t][1] += p_t[t][1] / N

p = {'u' : [item[0] for item in p_avg], 'd' : [item[1] for item in p_avg]}

#### PLOTS
fig, axs = plt.subplots(3)
plt.suptitle('N = %d ,  resolution = %.2d, sampling rate = %.1f' % (N, resolution, sample_rate))

axs[0].set_title(r"$T_1 = %.1f,\  \mu_1 = %.1f,\  T_2 = %.2f,\  \mu_2 = %.1f,\ T_3 = %.1f,\ \mu_3 = %.1f,\  e_u = %.1f, e_d = %.1f$" %(res1.T, res1.mu, res2.T, res2.mu, res3.T, res3.mu, res1.e, res2.e))
axs[0].plot(times, p['u'], label=r"$p_u$", color="blue")
axs[0].plot(times, p['d'], label=r'$p_d$', color="green")
axs[2].set_xlabel('time')
axs[0].legend(loc='upper right')

axs[1].plot(jumped_t[1:], I[1][1:], label=r"$I_{1}$", color="blue")
axs[1].plot(jumped_t[1:], I[2][1:], label=r"$I_{2}$", color="red")
axs[1].plot(jumped_t[1:], I[3][1:], label=r"$I_{3}$", color="green")
axs[2].set_xlabel('time')
axs[1].legend(loc='upper right')

axs[2].plot(jumped_t[1:], Q[1][1:], label=r"$Q_{1}$", color="blue")
axs[2].plot(jumped_t[1:], Q[2][1:], label=r"$Q_{2}$", color="red")
axs[2].plot(jumped_t[1:], Q[3][1:], label=r"$Q_{3}$", color="green")
axs[2].set_xlabel('time')
axs[2].legend(loc='upper right')

plt.show()
