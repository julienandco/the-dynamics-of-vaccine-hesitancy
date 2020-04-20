from random import seed
import random
import math

def sig(x):
    return (1/(1+math.exp(-x)))

def u_pro (omega, theta, i, epsilon):
    bias_theta = 0.0
    bias_omega = 0.0
    #je näher omega am höchsten wert für risiko ist (1/1000), desto fetter muss das negative werden -> 1-(1/1000 - omega) = 999/1000 + omega
    return ((1-omega) * i * epsilon * (theta + bias_theta) - (1000*omega - bias_omega)*(1-epsilon)*(1-theta)*(1-i)) * 10

def u_con (omega, theta, i, epsilon):
    bias_theta = 0.0
    bias_omega = 0.0
    bias_epsilon = 0.0
    return ((theta-bias_theta)*(epsilon-bias_epsilon)*(omega+bias_omega)-omega*(1-epsilon)*theta*i) * 10

def single_u(omega, theta, i, epsilon):
    return epsilon*theta*i*(1-omega)

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

seed(1)

def buildRandomSamples():
    samples = []
    #vaccine risk ist maximal 1/1000, vaccine efficacity ist mind. 70% sonst nicht zugelassen
    for i in range(1000):
        samples.append([random.uniform(0,1/1000), random.uniform(0,1), random.uniform(0,1), random.uniform(0.7, 1)])
    return samples

def compute(samples, fun_pro, fun_con):
    results = []
    for i in range(len(samples)):
        results.append([i, round(samples[i][0], 6), round(samples[i][1], 6), round(samples[i][2], 6), round(samples[i][3], 6), round(sig(-apply(samples[i], fun_pro)), 6), round(sig(apply(samples[i], fun_con)), 6)])
    return (fun_to_string(fun_pro), fun_to_string(fun_con), results)

def outputToLatex(res):
    txt = '\\begin{table}\n\\caption{Pro:' + res[0] + ', Con:' + res[1] + '}\n\\begin{tabular}{c|c|c|c|c|c|c}\n\\# & $\\omega$ & $\\theta$ & $I$ & $\\epsilon$ & Pro & Con' + "\\" + "\\\n" + '\\hline\n'
    counter = 0
    results = res[2]
    for i in range(len(results)):
        txt += str(results[i][0]) + ' & ' + str(results[i][1]) + ' & ' + str(results[i][2]) + ' & ' + str(results[i][3]) +  ' & ' + str(results[i][4]) + ' & ' + str(results[i][5]) + ' & ' + str(results[i][6]) + '\\\\\n'
        if i == counter + 50:
            counter = i
            txt += '\\end{tabular}\n\\end{table}\n\\newpage\n\\begin{table}\n\\begin{tabular}{c|c|c|c|c|c|c}\n\\# & $\\omega$ & $\\theta$ & $I$ & $\\epsilon$ & Pro & Con' + "\\" + "\\\n" + '\\hline\n'
    return (txt + '\\end{tabular}\n\\end{table}')

sam = buildRandomSamples()
res= compute(sam, u_pro, u_con)

with open('table_2_payoff_test.tex', 'w') as file:
    file.write(outputToLatex(res))
