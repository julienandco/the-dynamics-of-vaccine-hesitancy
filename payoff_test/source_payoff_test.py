from random import seed
import random


def u_pro (omega, theta, i, epsilon):
    bias_theta = 0.0
    bias_omega = 0.0
    return (1-omega) * i * epsilon * (theta + bias_theta) - (omega - bias_omega)*(1-epsilon)*(1-theta)*(1-i)

def u_con (omega, theta, i, epsilon):
    bias_theta = 0.0
    bias_omega = 0.0
    bias_epsilon = 0.0
    return (theta-bias_theta)*(epsilon-bias_epsilon)*(omega+bias_omega)-omega*(1-epsilon)*theta*i

def new_u(omega, theta, i, epsilon):
    return epsilon*theta*i*(1-omega)

def apply(samples, fun):
    omega = samples[0]
    theta = samples[1]
    i = samples[2]
    epsilon = samples[3]
    return fun(omega, theta, i , epsilon)

seed(1)

def buildRandomSamples():
    samples = []
    for i in range(1000):
        samples.append([random.uniform(0,1/1000), random.uniform(0,1), random.uniform(0,1), random.uniform(0.7, 1)])
    return samples

def compute(samples):
    results = []
    for i in range(len(samples)):
        results.append([i, round(samples[i][0], 6), round(samples[i][1], 6), round(samples[i][2], 6), round(samples[i][3], 6), round(apply(samples[i], u_con), 6), round(apply(samples[i], u_pro), 6)])
    return results

def outputToLatex(results):
    txt = '\\begin{table}\n\\begin{tabular}{c|c|c|c|c|c|c}\n\\# & $\\omega$ & $\\theta$ & $I$ & $\\epsilon$ & $u_{con}$ & $u_{pro}$' + "\\" + "\\\n" + '\\hline\n'
    counter = 0
    for i in range(len(results)):
        txt += str(results[i][0]) + ' & ' + str(results[i][1]) + ' & ' + str(results[i][2]) + ' & ' + str(results[i][3]) +  ' & ' + str(results[i][4]) + ' & ' + str(results[i][5]) + ' & ' + str(results[i][6]) + '\\\\\n'
        if i == counter + 50:
            counter = i
            txt += '\\end{tabular}\n\\end{table}\n\\newpage\n\\begin{table}\n\\begin{tabular}{c|c|c|c|c|c|c}\n\\# & $\\omega$ & $\\theta$ & $I$ & $\\epsilon$ & $u_{con}$ & $u_{pro}$' + "\\" + "\\\n" + '\\hline\n'
    return (txt + '\\end{tabular}\n\\end{table}')

sam = buildRandomSamples()
res= compute(sam)

with open('table_2_payoff_test.tex', 'w') as file:
    file.write(outputToLatex(res))