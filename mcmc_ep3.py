

# Import relevant modules
import numpy as np
import pandas as pd
import math
from scipy.stats import triang, norm
from runenergyplus import prepare_job_folders, run_eplus, generate_markup_value_pairs


#global markup_values_pairs, markup_value_pairs, path, totalarea

def tri_logdensity(x, min, max, mode):
    loc = min
    scale = max - min
    c = (mode - min)/scale
    log_density = triang.logpdf(x = x, c = c, loc = loc, scale = scale)
    return log_density

# Priors 
def prior(param):
    ROOF = param[0]
    WALL = param[1]
    WIN = param[2]
    SHGC = param[3]
    EPD = param[4]
    LPD = param[5]
    HSP = param[6]
    CSP = param[7]
    OCC = param[8]
    INF = param[9]
    Boiler = param[10]
    COP = param[11]

    ROOF_prior = tri_logdensity(x = ROOF, min = 0.01, max = 0.25, mode = 0.09667)
    WALL_prior = tri_logdensity(x = WALL, min = 0.01, max = 0.25, mode = 0.055)
    WIN_prior = tri_logdensity(x = WIN, min = 1, max = 5, mode = 2.792)
    SHGC_prior = tri_logdensity(x = SHGC, min = 0.1, max = 0.9, mode = 0.5)
    EPD_prior = tri_logdensity(x = EPD, min = 1, max = 70, mode = 22.8)
    LPD_prior = tri_logdensity(x = LPD, min = 1, max = 50, mode = 14.5)
    HSP_prior = tri_logdensity(x = HSP, min = 17, max = 28, mode = 21)
    CSP_prior = tri_logdensity(x = CSP, min = 17, max = 28, mode = 24)
    OCC_prior = tri_logdensity(x = OCC, min = 1, max = 50, mode = 23.4)
    INF_prior = tri_logdensity(x = INF, min = 0.1, max = 2, mode = 0.675)
    Boiler_prior = tri_logdensity(x = Boiler, min = 0.5, max = 0.99, mode = 0.72)
    COP_prior = tri_logdensity(x = COP, min = 2, max = 4, mode = 2.65)
    return ROOF_prior+WALL_prior+WIN_prior+SHGC_prior+EPD_prior+LPD_prior\
        +HSP_prior+CSP_prior+OCC_prior+INF_prior+Boiler_prior+COP_prior



#Likelihood function
def likelihood(param):
    ROOF = param[0]
    WALL = param[1]
    WIN = param[2]
    SHGC = param[3]
    EPD = param[4]
    LPD = param[5]
    HSP = param[6]
    CSP = param[7]
    OCC = param[8]
    INF = param[9]
    Boiler = param[10]
    COP = param[11]
    
    chain = make_chainlist(param)
    
    markup_values_pairs = dict(zip(['@@ROOF@@','@@WALL@@','@@WIN@@',
                                '@@SHGC@@','@@EPD@@','@@LPD@@',
                                '@@HSP@@','@@CSP@@','@@OCC@@',
                                '@@INF@@','@@Boiler@@','@@COP@@'], chain))
    prediction = param[12]
        
    singlelikelihoods = norm.logpdf(x=y, loc=prediction, scale=sd)
    #sumll = sum(singlelikelihoods)
    return singlelikelihoods


#Posterior 
def posterior(param):
    return likelihood(param) + prior(param)


#Proposal function- .tolist:to change numpy.ndarray to list

def proposalfunction(param):
    return abs(np.random.normal(loc=param, scale=[0.02,0.015,0.43,0.1,7,4,1.3,1.3,6.75,0.19,0.06,0.167])).tolist()



#Run metro-polis MCMC ######################################TEST
def run_metropolis_MCMC(startvalue, iterations, output_folder, template_idf_path, eplus_basic_folder, 
                        path, totalarea):
    chain = [[0 for x in xrange(12)] for x in xrange(iterations)]
    chain[0][0] = startvalue[0]
    chain[0][1] = startvalue[1]
    chain[0][2] = startvalue[2]
    chain[0][3] = startvalue[3]
    chain[0][4] = startvalue[4]
    chain[0][5] = startvalue[5]
    chain[0][6] = startvalue[6]
    chain[0][7] = startvalue[7]
    chain[0][8] = startvalue[8]
    chain[0][9] = startvalue[9]
    chain[0][10] = startvalue[10]
    chain[0][11] = startvalue[11]
    
    prediction = run_eplus(path, totalarea)
    chain[0].append(prediction)
    
    for i in range(1,iterations):
        proposal = proposalfunction(chain[i-1][:-1])
        chainlist = make_chainlist(proposal)
        markup_values_pairs = dict(zip(['@@ROOF@@','@@WALL@@','@@WIN@@',
                                '@@SHGC@@','@@EPD@@','@@LPD@@',
                                '@@HSP@@','@@CSP@@','@@OCC@@',
                                '@@INF@@','@@Boiler@@','@@COP@@'], chainlist))
        # prepares jobs
        markup_value_pairs = generate_markup_value_pairs(markup_values_pairs, count)
        path = prepare_job_folders(output_folder, template_idf_path, eplus_basic_folder, markup_value_pairs)
        prediction = run_eplus(path, totalarea)
        proposal.append(prediction) 
        
        probab = min(1, math.exp(posterior(proposal) - posterior(chain[i-1])))
        if np.random.uniform() < probab :
            chain[i] = list(proposal)
        else:
            chain[i] = list(chain[i-1])
    return chain


# To make chain list
def make_chainlist(list):
    chain = []
    for item in list:
        temp =[item]
        chain.append(temp)
    return chain



#Run
# Observation value
y = 1.619888
sd = 0.43
startvalue = [0.09667, 0.055, 2.792, 0.5, 22.8, 14.5, 21, 24, 23.4, 0.675, 0.72, 2.65]
iterations = 10
# initial values
template_idf_path = "test/campusbuilding.idf"
eplus_basic_folder = "test/basic_files"
output_folder = "test/out"
#sys.argv[1]
totalarea = 10336.99
count = 1

chain = make_chainlist(startvalue)
markup_values_pairs = dict(zip(['@@ROOF@@','@@WALL@@','@@WIN@@',
                                '@@SHGC@@','@@EPD@@','@@LPD@@',
                                '@@HSP@@','@@CSP@@','@@OCC@@',
                                '@@INF@@','@@Boiler@@','@@COP@@'], chain))

markup_value_pairs = generate_markup_value_pairs(markup_values_pairs, count)
path = prepare_job_folders(output_folder, template_idf_path, eplus_basic_folder, markup_value_pairs)

import time
start_time = time.time()
chain = run_metropolis_MCMC(startvalue, iterations, output_folder, template_idf_path, eplus_basic_folder, path, totalarea)
print time.time() - start_time, "seconds"
print chain

import csv
import os
    
with open(os.path.join("test", 'chain_ep.csv'), 'wb') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames = ['ROOF','WALL','WIN','SHGC','EPD','LPD',
                    'HSP','CSP','OCC','INF','Boiler','COP','EUI'], delimiter = ',')
    writer.writeheader()
    writer = csv.writer(csvfile)
    writer.writerows(chain)
