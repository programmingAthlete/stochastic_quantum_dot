import numpy as np
import matplotlib.pyplot as plt
import random
import sys

def empiric_dist(p, t):
    resistence_t = np.array([p[i] * (t[i+1] - t[i]) for i in range( len(p)-1)])
    p_dist = np.array([0.0 for _ in range(len(p))])
    for i in range(len(p_dist)):
        if t[i] == 0.0:
            p_dist[i] = 0.
        else:
            p_dist[i] = sum(resistence_t[:i]) / t[i]
    return p_dist, t

def qd_trajectory(MAXTIME, times, rates):
    """
    Stochastic simulation algorithm with fixed time resolution
    """

    rate_u1 = rates['rate_up1']
    rate_u2 = rates['rate_up2']
    rate_d1 = rates['rate_down1']
    rate_d2 = rates['rate_down2']

    t = 0
    p = 0 # initial population: initially void quantum dot
    jumps_t = [0]
    resolution = len(times)
    p_t = [0 for i in range(resolution)]

    de1 = 0
    de2 = 0
    I_1 = 0
    I_2 = 0
    de1_t = [0]
    de2_t = [0]
    I_1_t =[0]
    I_2_t =[0]
    p_dist = []

    i = 0
    while i < resolution:
        pold = p

        while t <= times[i]:
            rate_u = (1 - p) * (rate_u1 + rate_u2)
            rate_d = p * (rate_d1 + rate_d2)
            rate = rate_u + rate_d
            tau = - (1.0 / rate) * np.log(random.random())
            t += tau
            jumps_t.append(t)
            r = random.random()
            if p == 0:
                if r < rate_u1/rate_u:
                    de1 = e-mu1
                    de2 = 0
                    I_1 = 1
                    I_2 = 0
                else:
                    de2 = e -mu2
                    de1 = 0
                    I_1 = 0
                    I_2 = 1

            elif p == 1:
                if r < rate_d2/rate_d:
                    de2 = -e + mu2
                    de1 = 0
                    I_2 = -1
                    I_1 = 0
                else:
                    de1 = -e + mu1
                    de2 = 0
                    I_2 = 0
                    I_1 = -1
            p = 1 - p
            de1_t.append(de1)
            de2_t.append(de2)
            I_1_t.append(I_1)
            I_2_t.append(I_2)
            p_dist.append(p)

        time_idx = np.searchsorted(times > t, True)
        corr_idx = min(resolution - 1, time_idx)

        for j in range(i, corr_idx - 1):
            p_t[j] = pold
        i = time_idx
    return p_t, de1_t, de2_t, I_1_t, I_2_t, jumps_t, p_dist

def ensemble():
    # kinetic constants
    global T1,T2,mu1,mu2, k1,k2,e
    k1 = 1.
    k2 = 1.
    T1 = 1.
    T2 = 5.2
    mu1 = 1
    mu2 = 1
    # time interval
    N=500
    MAXTIME = 3000
    t = 0
    # Thermodynamic parameters

    flag = ''
    if mu1*T1 == mu2*T2:
        flag = 'equilibrium'
    e = 2.5
    f1 = 1./(1. + np.exp((e-mu1)/T1))
    f2 = 1./(1. + np.exp((e-mu2)/T2))

    rate_u1 = k1 * f1
    rate_u2 = k2 * f2
    rate_d1 = k1 * (1 - f1)
    rate_d2 = k2 * (1 - f2)

    rates = {'rate_up1' : rate_u1,
                'rate_down1' : rate_d1 ,
                'rate_up2' : rate_u2,
                'rate_down2': rate_d2
                }

    res = 120000
    avgp_t = [0 for i in range(res)]
    ddQ1 = [0 for i in range(res)]
    ddQ2 = [0 for i in range(res)]
    dQ1 = [0]
    dQ2 = [0]
    I1 = [0]
    I2 = [0]
    jumped_t = [0]

    second_moment = [0 for i in range(res)]
    times = np.linspace(0, MAXTIME, res)
    # Check for undersampling
    srate = res / MAXTIME
    timescale = max(rate_d1 + rate_d2, rate_u1 + rate_u2)
    print("[+] Sampling rate: {:.2f}, expected up/down rates {:.2f}, {:.2f}".format(srate,
                                                                        rate_u1 + rate_u2,
                                                                        rate_d1 + rate_d2))
    if  srate < 2 * timescale:
        print("[!] Warning: Possible undersampling resolving SSA")

    dist_flag = False
    for j in range(N):
        p_t, dq1_t, dq2_t, n_1_t, n_2_t, jumps_t, p_0 = qd_trajectory(MAXTIME, times, rates)

        if dist_flag == False:

            dddQ1 = np.array([0.0 for _ in range(len(jumps_t))])
            dddQ2 = np.array([0.0 for _ in range(len(jumps_t))])

            ## Deirvate Q1 and Q2
            dt = 0.001
            for i in range(len(jumps_t)):
                jumped_t.append(jumps_t[i])
                if i != 0 and i != len(times) - 1:
                    dddQ1[i] = dq1_t[i]
                    dddQ2[i] = dq2_t[i]
                else:
                    dddQ1[i] = dq1_t[i]
                    dddQ2[i] = dq2_t[i]

            ddQ1 = np.array([0.0 for _ in range(len(jumps_t))])
            ddQ2 = np.array([0.0 for _ in range(len(jumps_t))])
            dI1 = np.array([0.0 for _ in range(len(jumps_t))])
            dI2 = np.array([0.0 for _ in range(len(jumps_t))])
            for i in range(len(ddQ1)):
                if times[i] == 0.0:
                    ddQ1[i] = dddQ1[i]
                    ddQ2[i] = dddQ2[i]
                    dI1[i] = n_1_t[i]
                    dI2[i] = n_2_t[i]
                else:
                    ddQ1[i] = sum(dddQ1[:i]) /jumps_t[i]
                    ddQ2[i] = sum(dddQ2[:i]) /jumps_t[i]
                    dI1[i] = sum(n_1_t[:i]) / jumps_t[i]
                    dI2[i] = sum(n_2_t[:i]) / jumps_t[i]
                dQ1.append(ddQ1[i])
                dQ2.append(ddQ2[i])
                I1.append(dI1[i])
                I2.append(dI2[i])
            dist_flag = True


        ## Mean on N trajectoreis
        for t in range(res):
            avgp_t[t] += p_t[t] / N
            second_moment[t] += p_t[t]**2 / N

    stddev = [np.sqrt((second_moment[t] - avgp_t[t]**2) / N) for t in range(res)]
    I1 = np.array(I1)
    I2 = np.array(I2)
    if mu1*T1 >= mu2*T2:
        F_m = ( mu1 / T1 - mu2 /T2)
        F_e = e* (1 / T1 - 1 / T2)
        sigma =  F_m * I1 - I1 * F_e
    elif mu2*T2 > mu1*T1:
        F_m = ( mu2 / T2 - mu1 /T1)
        F_e =  e*(1 / T2 - 1 / T1)
        sigma = F_m  * I2 - F_e * I2

    #####################################
    #### Theoretical calculations #######
    f1 = 1./(1. + np.exp((e-mu1)/T1))
    f2 = 1./(1. + np.exp((e-mu2)/T2))
    p_1 = f1 / 2.0
    p_2 = f2 / 2.0
    p_1Array = np.array([p_1 for _ in range(len(times))])
    p_2Array = np.array([p_2 for _ in range(len(times))])
    k = k1*k2 / (k1+k2)
    I1_theory = np.array([ k *(f1 - f2) for _ in range(len(times))])
    I2_theory = np.array([k *(f2 - f1) for _ in range(len(times))])
    dQ1_theory = ((e - mu1) * I1_theory)
    dQ2_theory = ((e - mu2) * I2_theory)
    sigma_theory = k * ((e-mu1) / T1 -  (e-mu2) /T2)* (f2 - f1)
    signma_Theory = np.array([sigma_theory for _ in range(len(times))])
    # plot
    fig, axs = plt.subplots(3)
    #fig, (ax0) = plt.subplots(nrows=1, figsize=(9,9))
    fig.suptitle('N = %d ,  resolution = %.2d, sampling rate = %.1f' % (N, res, srate))
    # empirical probability
    axs[0].plot(times,
             avgp_t,
             label=r'$p(t)$')
    zalpha = 1.96
    upperbar = [avgp_t[t] + zalpha * stddev[t] for t in range(res)]
    lowerbar = [avgp_t[t] - zalpha * stddev[t] for t in range(res)]
    axs[0].fill_between(times, upperbar, lowerbar, alpha=0.3,
                     color='darkorange', label=r'$ Gaussian CI $')
    # stationary probability
    pss = (rate_u1 + rate_u2) / (rate_u1 + rate_u2 + rate_d1 + rate_d2)
    peqleft = rate_u1 / (rate_u1 + rate_d1)
    peqright = rate_u2 / (rate_u2 + rate_d2)
    ptransient = [ avgp_t[0] * np.exp(-ti * (rate_u1 + rate_u2 + rate_d1 + rate_d2))
                   + pss * (1 - np.exp(-ti * (rate_u1 + rate_u2 + rate_d1 + rate_d2)))
                   for ti in times ]
    #axs[0].plot(t_dist[:len(t_dist)-1], p_dist, color='green', label=r'$p_{ECDF}(t)$')
    axs[0].plot(times, ptransient, 'y--', label=r'$p(t)$ (Theory)')
    #axs[0].plot(times, [peqright for i in times], 'r--', label=r'$p_{eq}^{L}$')
    #axs[0].plot(times, [peqleft for i in times], 'b--', label=r'$p_{eq}^{R}$')

    axs[0].set_ylabel(r'$p(n_{dot}=1)$')
    axs[0].set_xlabel(r'$t$')
    axs[0].set_title(r"Evolution of the occupation probability $(\mu_L T_L={:2.1f}, \mu_R T_R={:2.1f}, e={:2.1f})$".format(T1 * mu1, T2 * mu2, e))
    axs[0].legend(loc='upper right')

    axs[1].plot(jumped_t, dQ1,  color='blue',linewidth=1, label=r'$\dot{Q}^L(t)$')
    axs[1].plot(times, dQ1_theory, ls='--', color='blue', linewidth=1,label=r'$\dot{Q}^L(t)$ (Theory)')
    axs[1].plot(jumped_t, dQ2,  color='red', linewidth=1,label=r'$\dot{Q}^R(t)$')
    axs[1].plot(times, dQ2_theory, ls='--', color='red', linewidth=1,label=r'$\dot{Q}^R(t)$ (Theory)')
    #axs[1].set_ylim([-0.5,1.25])
    axs[1].legend(loc='upper right')


    axs[2].plot(jumped_t, I1, color='blue',linewidth=1, label=r'$I_L(t)$')
    axs[2].plot(times, I1_theory, ls="--", color="blue", linewidth=1, label=r'$I_L(t)$ (Theory)')
    axs[2].plot(jumped_t, I2, color='red',linewidth=1, label=r'$I_R(t)$')
    axs[2].plot(times, I2_theory, ls="--", color="red", linewidth=1, label=r'$I_R(t)$ (Theory)')
    axs[2].plot(jumped_t, sigma, color='green', linewidth=1, label=r'$\sigma (t)$')
    axs[2].plot(times, signma_Theory, ls="--",color='green', linewidth=1, label=r'$\sigma (t)$ (Theory)')
    #axs[2].set_ylim([-0.5,1.00])
    axs[2].legend(loc='upper right')

    j = 0
    for ax in axs.flat:
        if j == 0:
            ax.set(xlabel='time', ylabel=r'$p^{ss}(t)$')
        else:
            ax.set(xlabel="time")
        j += 1


    plt.show()

ensemble()
