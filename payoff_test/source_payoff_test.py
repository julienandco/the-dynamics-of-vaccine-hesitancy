from random import seed
import random
import math
import numpy as np
import matplotlib.pyplot as plt

#######################################################################################################
# first try to find payoff function (20 - 26 April 2020)
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

#######################################################################################################
# second try to find payoff function (27 - 30 April 2020)

def u_con_mark2(omega, theta, i, epsilon, kappa, X, p_meet_Z1, p_meet_Z0):
    return (X + p_meet_Z1)*(theta + i) - (2 - X - p_meet_Z1)*(1-epsilon+kappa*omega)

def u_pro_mark2(omega, theta, i, epsilon, kappa, X, p_meet_Z1, p_meet_Z0):
    return (2-1+X-p_meet_Z0) * (theta + i) - (1-X + p_meet_Z0) * (1-epsilon + kappa*omega)


#######################################################################################################
# third try (31 April - 7 May 2020) this time: payoff included in transition rates!

def transition_pro(omega, theta, i, epsilon, kappa, N, X, N_X, N_Y, theta_X, theta_Y, mu):
    return (mu * (N-X)/N * ((theta_X * (X/N + N_X/N + theta + i)) / (theta_X * (X/N + N_X/N + theta + i) + ((N-X)/N + N_Y/N + 1 - epsilon + kappa * omega))))

def transition_con(omega, theta, i, epsilon, kappa, N, X, N_X, N_Y, theta_X, theta_Y, mu):
    return (mu * X/N * ((theta_Y * ((N-X)/N + N_Y/N + 1 - epsilon + kappa * omega)) / ((X/N + N_X/N + theta + i) + theta_Y * ((N-X)/N + N_Y/N + 1 - epsilon + kappa * omega))))

#######################################################################################################

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
    elif fun == u_con_mark2:
        return '$u_{con} (mark 2)$'
    elif fun == u_pro_mark2:
        return '$u_{pro} (mark 2)$'
    elif fun == transition_con:
        return 'rate for $X_t \\rightarrow X_t - 1$'
    elif fun == transition_pro:
        return 'rate for $X_t \\rightarrow X_t + 1$'
    else:
        return 'You need to define a toString for this function'

def apply(samples, fun_pro, fun_con, **kwargs):
    args = kwargs.get('args',[])
    arg_list = samples
    for arg in args:
        arg_list.append(args[arg])
    res_pro = fun_pro(*arg_list)
    res_con = fun_con(*arg_list)
    return (res_pro,res_con)

def buildRandomSamples(numberOfSamples, kappa, allDynamic, maxElementsPerPage):
    samples = []
    if allDynamic:
        for j in range(numberOfSamples):
            omega = random.uniform(0, 1/kappa)
            theta = random.uniform(0, 1)
            i = random.uniform(0, 1)
            epsilon = random.uniform(0.75,1)
            samples.append([omega,theta,i,epsilon])
    else:
        omega = random.uniform(0, 1/kappa)
        theta = random.uniform(0, 1)
        epsilon = random.uniform(0.75,1)
        counter = 0
        for j in range(numberOfSamples):
            i = random.uniform(0, 1)
            samples.append([omega,theta,i,epsilon])
            if j == counter + maxElementsPerPage:
                counter = j
                omega = random.uniform(0, 1/kappa)
                theta = random.uniform(0, 1)
                epsilon = random.uniform(0.75,1)
    return samples

def buildRandomArguments(fun,kappa,**kwargs):
    args = {'\\kappa' : kappa}
    if fun == u_pro_mark2:
        args['X_t'] = round(random.uniform(0,1),6)
        args['p_{meetZ_1}'] = round(random.uniform(0,1),6)
        args['p_{meetZ_0}'] = round(random.uniform(0,1),6)
    elif fun == transition_pro:
        N = kwargs.get('N',0)
        max_zel_X = kwargs.get('maxZelX',0)
        max_zel_Y = kwargs.get('maxZelY',0)
        args['N'] = N
        args['X_t'] = random.randint(0,N)
        args['N_X']  = random.randint(0,max_Zel_X)
        args['N_Y'] = random.randint(0,max_Zel_Y)
        args['\\theta_X'] = round(random.uniform(0,1),6)
        args['\\theta_Y'] = round(random.uniform(0,1),6)
        args['\\mu'] = round(random.uniform(0,1),6)
    return args

def compute(samples,arguments, fun_pro, fun_con, **kwargs):
    results = []
    app_sigmoid = kwargs.get('sig',False)

    for i in range(len(samples)):
        current_res = [i]
        #all sampled parameters needed for computation, rounded to 6 digits
        for j in range(len(samples[i])):
            current_res.append(round(samples[i][j],6))
        #the two computation results (pro and con)
        res_pro,res_con = apply(samples[i], fun_pro, fun_con,args = arguments)

        #do we want the sigmoid function to transform the results?
        if app_sigmoid:
            current_res.append(round(sig(-res_pro),6))
            current_res.append(round(sig(res_con),6))
        else:
            current_res.append(round(res_pro,6))
            current_res.append(round(res_con,6))
        #push everything in results
        results.append(current_res)
    
    return (fun_to_string(fun_pro), fun_to_string(fun_con), str(app_sigmoid), results)

def getLatexIntro(args, skipKappa):
    intro = 'Following fixed values are considered:\n\\begin{itemize}\n'
    for itemName in args:
        if not (itemName == '\\kappa' and skipKappa):
           intro += '\\item $' + itemName + ' = ' + str(args[itemName]) + '$\n'
    intro += '\\end{itemize}\n'
    return intro

def getLatexSampleExplanation(counter,allDynamic):
    if allDynamic:
        expl = 'Following values are completely dynamic and (randomly) generated every step:\n\\begin{itemize}\n'
        expl += '\\item $\\omega$\n\\item $\\theta$\n\\item $\\epsilon$\n\\item $I$\n\\end{itemize}\n'
    else:
        expl = 'Following values are newly (randomly) generated all ' + str(counter) + ' steps:\n\\begin{itemize}\n'
        expl += '\\item $\\omega$\n\\item $\\theta$\n\\item $\\epsilon$\n\\end{itemize}\nThe value $I$ is completely dynamic and determined randomly at every step.'
    return expl

def getLatexCaption(res):
    fun_pro = res[0]
    fun_con = res[1]
    sig = res[2]
    results = res[3]
    if fun_pro == fun_to_string(transition_pro):
        caption = '\\begin{table}\n\\caption{Pro: ' + fun_pro + ', Con: ' + fun_con + '}\n\\begin{tabular*}{\\linewidth}{'
    else:
        caption = '\\begin{table}\n\\caption{Pro: ' + fun_pro + ', Con: ' + fun_con + ', $\\mathrm{sigmoid}: ' + sig + '$}\n\\begin{tabular*}{\\linewidth}{'
    for i in range(len(results[0])):
        if i == len(results[0]) - 1:
            caption += 'c}\n'
        else:
            caption += 'c|'
    caption += '\\# & $\\omega$ & $\\theta$ & $I$ & $\\epsilon$ & Pro & Con \\\\\n\\hline\n'
    return caption

def getLatexLine(results):
    line = ''
    for i in range(len(results)):
        if i == len(results) - 1:
            line += str(results[i]) + '\\\\\n'
        else:
            line += str(results[i]) + ' & '
    return line

def outputToLatex(res,args,skipKappa,maxElementsPerPage,allDynamic):
    txt = ''
    txt += getLatexIntro(args, skipKappa)
    expl = getLatexSampleExplanation(maxElementsPerPage,allDynamic)
    txt += expl
    caption = getLatexCaption(res)
    txt += caption
    results = res[3]
    counter = 0
    for i in range(len(results)):
        txt += getLatexLine(results[i])
        if i == counter + maxElementsPerPage:
            counter = i
            txt += '\\end{tabular*}\n\\end{table}\n\\newpage\n' + caption

    return (txt + '\\end{tabular*}\n\\end{table}')





#######################################################################################################
# set all parameters here
kappa = 100000
pop = 1000000
max_Zel_X = 0.15 * pop
max_Zel_Y = 0.05 * pop
maxElementsPerPage = 45
numberOfSamples = 1000
allDynamic = False
skipKappa = False

# reminder: if you are using transition_pro and transition_con, useSigmoid should always be set to False
#           (True makes no sense and: you get a math overflow error lel)
useSigmoid = False

fun_pro = transition_pro
fun_con = transition_con

#name of the Latex file (in this directory) you are writing:
fileName = 'table_5_payoff_test.tex'

# leave commented if you want unrepeated random numbers, uncomment if you want to have the same numbers every compilation
#seed(1)








#######################################################################################################
# do not touch!

samples = buildRandomSamples(numberOfSamples,kappa, allDynamic, maxElementsPerPage)
#note: we always pass all args, because of the if-else in buildRandomArgs: if we want args for u_pro_mark2 f.ex., he ain't gonna look at maxZel etcx anyhoobs!
args = buildRandomArguments(transition_pro,kappa,N=pop,maxZelX = max_Zel_X, maxZelY = max_Zel_Y)
res = compute(samples,args, fun_pro, fun_con, sig=useSigmoid)

with open(fileName, 'w+') as file:
    file.write(outputToLatex(res,args,skipKappa,maxElementsPerPage,allDynamic))