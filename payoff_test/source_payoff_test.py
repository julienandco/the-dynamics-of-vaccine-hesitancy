from random import seed
import random
import math

def sig(x):
    return (1/(1+math.exp(-x)))

def u_pro (omega, theta, i, epsilon):
    bias_theta = 0.0
    bias_omega = 0.0
    #je näher omega am höchsten wert für risiko ist (1/1000), desto fetter muss das negative werden -> 1000*omega
    return ((1-(1000*omega)) * i * epsilon * (theta + bias_theta) - (1000*omega - bias_omega)*(1-epsilon)*(1-theta)*(1-i)) * 10

def u_con (omega, theta, i, epsilon):
    #noch sehr whack
    bias_theta = 0.0
    bias_omega = 0.0
    bias_epsilon = 0.0
    return ((theta-bias_theta)*(epsilon-bias_epsilon)*(1-(1000*omega)+bias_omega)-1000*omega*(1-epsilon)*(1-theta)*(1-i)) * 10

def single_u(omega, theta, i, epsilon):
    return epsilon*theta*i*(1-omega)

def counter_u_pro(omega, theta, i, epsilon):
    omega_payoff = -1
    theta_payoff = -1
    i_payoff = -1
    epsilon_payoff = -1
    if i > 0.4:
        i_payoff = 1
    if omega < 1/100000:
        omega_payoff = 1
    if theta > 0.25:
        theta_payoff = 1
    if epsilon > 0.8:
        epsilon_payoff = 1
    
    return 2.6 * i_payoff + 0.2 * omega_payoff + 0.4 * theta_payoff + 0.8 * epsilon_payoff


def counter_u_con(omega, theta, i, epsilon):
    omega_payoff = -1
    theta_payoff = -1
    i_payoff = -1
    epsilon_payoff = -1
    if i > 0.75:
        i_payoff = 1
    if omega < 1/1000000:
        omega_payoff = 1
    if theta > 0.6:
        theta_payoff = 1
    if epsilon > 0.95:
        epsilon_payoff = 1
    
    return 0.2 * i_payoff + 3.2 * omega_payoff + 0.28 * theta_payoff + 0.32 * epsilon_payoff

def apply(samples, fun):
    omega = samples[0]
    theta = samples[1]
    i = samples[2]
    epsilon = samples[3]
    return fun(omega, theta, i , epsilon)

def fun_to_string(fun):
    if fun == u_pro:
        return '$u_{pro}$'
    elif fun == u_con:
        return '$u_{con}$'
    elif fun == single_u:
        return '$single_{u}$'
    elif fun == counter_u_con:
        return '$counter_{u_{con}}$'
    elif fun == counter_u_pro:
        return '$counter_{u_{pro}}$'

seed(1)

def buildRandomSamples():
    samples = []
    #vaccine risk ist maximal 1/1000, vaccine efficacity ist mind. 75% sonst nicht zugelassen
    for i in range(1000):
        samples.append([random.uniform(0,1/10000), random.uniform(0,1), random.uniform(0,1), random.uniform(0.75, 1)])
    return samples

def compute(samples, fun_pro, fun_con, app_sigma):
    results = []
    for i in range(len(samples)):
        if app_sigma:
            results.append([i, round(samples[i][0], 6), round(samples[i][1], 6), round(samples[i][2], 6), round(samples[i][3], 6), round(sig(-apply(samples[i], fun_pro)), 6), round(sig(apply(samples[i], fun_con)), 6)])
        else:
            results.append([i, round(samples[i][0], 6), round(samples[i][1], 6), round(samples[i][2], 6), round(samples[i][3], 6), round(apply(samples[i], fun_pro), 6), round(apply(samples[i], fun_con), 6)])
    return (fun_to_string(fun_pro), fun_to_string(fun_con), results)

def outputToLatex(res, sig):
    txt = '\\begin{table}\n\\caption{Pro:' + res[0] + ', Con:' + res[1] + '$\\mathrm{sig}:' + str(sig) + '$}\n\\begin{tabular}{c|c|c|c|c|c|c}\n\\# & $\\omega$ & $\\theta$ & $I$ & $\\epsilon$ & Pro & Con' + "\\" + "\\\n" + '\\hline\n'
    counter = 0
    results = res[2]
    for i in range(len(results)):
        txt += str(results[i][0]) + ' & ' + str(results[i][1]) + ' & ' + str(results[i][2]) + ' & ' + str(results[i][3]) +  ' & ' + str(results[i][4]) + ' & ' + str(results[i][5]) + ' & ' + str(results[i][6]) + '\\\\\n'
        if i == counter + 45:
            counter = i
            txt += '\\end{tabular}\n\\end{table}\n\\newpage\n\\begin{table}\n\\begin{tabular}{c|c|c|c|c|c|c}\n\\# & $\\omega$ & $\\theta$ & $I$ & $\\epsilon$ & Pro & Con' + "\\" + "\\\n" + '\\hline\n'
    return (txt + '\\end{tabular}\n\\end{table}')

sam = buildRandomSamples()
sig = False
res= compute(sam, counter_u_pro, counter_u_con, sig)

with open('table_2_payoff_test.tex', 'w') as file:
    file.write(outputToLatex(res, sig))
