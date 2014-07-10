# -------------------------------------------------------------------
# Markov Chain Monte Calro for EnergyPlus
# 
# code by Hyunwoo Lim (thanks to sangheestyle)
# 07/07/2014
# -------------------------------------------------------------------
# UM Campus 
# This code is for 10 parameters. (exclude WIN, SHGC)
# Parameter ranges are changed. 
# You need to modify the idf file to insert mark. (e.g. @@WALL@@, @@WIN@@)


# Import relevant modules
import numpy as np
import pandas as pd
import math
import os
from scipy.stats import triang, norm, truncnorm
from runenergyplus import prepare_job_folders, run_eplus
import csv

# Make CSV file 
def make_csv(chain):
    with open(os.path.join(result_CSV_path), 'wb') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames = ['ROOF','WALL','EPD','LPD',
                    'HSP','CSP','OCC','INF','Boiler','COP','EUI'], delimiter = ',')
        writer.writeheader()
        writer = csv.writer(csvfile)
        writer.writerows([chain])
# Add
def add_chain(chain):
    writer = csv.writer(open(os.path.join(result_CSV_path), "ab"), delimiter=',', quotechar='|')
    writer.writerow(chain)


def tri_logdensity(x, min, max, mode):
    loc = min
    scale = max - min
    c = (mode - min)/scale
    log_density = triang.logpdf(x = x, c = c, loc = loc, scale = scale)
    return log_density

def truncated_normal(myclip_min, myclip_max, mu, sigma, n=1):
    a = (myclip_min - mu) / sigma
    b = (myclip_max - mu) / sigma
    c = truncnorm.rvs(a, b, loc=mu, scale=sigma,size=n)
    return c

# Priors 
def prior(param):
    ROOF = param[0]
    WALL = param[1]
    #WIN = param[2]
    #SHGC = param[3]
    EPD = param[2]
    LPD = param[3]
    HSP = param[4]
    CSP = param[5]
    OCC = param[6]
    INF = param[7]
    Boiler = param[8]
    COP = param[9]

    ROOF_prior = tri_logdensity(x = ROOF, min = 0.01, max = 0.3, mode = 0.116)
    WALL_prior = tri_logdensity(x = WALL, min = 0.01, max = 0.3, mode = 0.046)
    #WIN_prior = tri_logdensity(x = WIN, min = 1, max = 5, mode = 2.792)
    #SHGC_prior = tri_logdensity(x = SHGC, min = 0.1, max = 0.9, mode = 0.5)
    EPD_prior = tri_logdensity(x = EPD, min = 1, max = 70, mode = 11.67)
    LPD_prior = tri_logdensity(x = LPD, min = 1, max = 30, mode = 12.43)
    HSP_prior = tri_logdensity(x = HSP, min = 17, max = 28, mode = 21)
    CSP_prior = tri_logdensity(x = CSP, min = 17, max = 28, mode = 24)
    OCC_prior = tri_logdensity(x = OCC, min = 2, max = 56.7, mode = 14.37)
    INF_prior = tri_logdensity(x = INF, min = 0.1, max = 1.5, mode = 0.56)
    Boiler_prior = tri_logdensity(x = Boiler, min = 0.5, max = 0.95, mode = 0.72)
    COP_prior = tri_logdensity(x = COP, min = 2, max = 4, mode = 2.65)
    return ROOF_prior+WALL_prior+LPD_prior+LPD_prior+ \
            HSP_prior+CSP_prior+OCC_prior+INF_prior+Boiler_prior+COP_prior



#Likelihood function
def likelihood(param):
    prediction = param[10]
    singlelikelihoods = norm.logpdf(x=y, loc=prediction, scale=sd)
    #sumll = sum(singlelikelihoods)
    return singlelikelihoods


#Posterior 
def posterior(param):
    return likelihood(param) + prior(param)


#Proposal function- .tolist:to change numpy.ndarray to list
def proposalfunction(param):
    return np.random.normal(loc=param, 
        scale=[0.02, 0.015, 3.5, 3.5, 1.3, 1.3, 4, 0.15, 0.06, 0.15]).tolist()
 
"""
#Proposal function using truncated_normal
def proposalfunction(low_limit, upper_limit, mean, proposal_sd):
    proposal = []
    for i in range(0,len(mean)):
        proposal.extend(truncated_normal(low_limit[i], upper_limit[i], mean[i], proposal_sd[i]).tolist())
    return proposal
"""        

def generate_markup_value_pairs(markup, chain):
    markup_value_pairs = []
    markup_value_pairs.append(dict(zip(markup, chain)))
    return markup_value_pairs
"""
[{'@@Boiler@@': 0.72,
  '@@COP@@': 2.65,
  '@@CSP@@': 24,
  '@@EPD@@': 22.8,
  '@@HSP@@': 21,
  '@@INF@@': 0.675,
  '@@LPD@@': 14.5,
  '@@OCC@@': 23.4,
  '@@ROOF@@': 0.09667,
  '@@SHGC@@': 0.5,
  '@@WALL@@': 0.055,
  '@@WIN@@': 2.792}]
"""

#Run metro-polis MCMC ######################################TEST

def run_metropolis_MCMC(path):
    chain = [[0 for x in xrange(10)] for x in xrange(iterations)]
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
    
    prediction = run_eplus(path, totalarea)
    chain[0].append(prediction)
    make_csv(chain[0]) 

    for i in range(1,iterations):
        while True:
            try:
                current_dir = os.getcwd()
                #proposal = proposalfunction(low_limit, upper_limit, chain[i-1][:-1], proposal_sd)
                proposal = proposalfunction(chain[i-1][:-1])

                markup_value_pairs = generate_markup_value_pairs(markup, proposal)
                path = prepare_job_folders(output_folder, template_idf_path, 
                    eplus_basic_folder, markup_value_pairs)
                print "___________________________",i, "th iteration____________________________________"
                prediction = run_eplus(path, totalarea)
                proposal.append(prediction) 
                break
            except:
                os.chdir(current_dir)
                pass

        
        probab = min(1, math.exp(posterior(proposal) - posterior(chain[i-1])))
        if np.random.uniform() < probab :
            chain[i] = list(proposal)
            print chain[i]
            add_chain(chain[i])
        else:
            chain[i] = list(chain[i-1])
            print chain[i]
            add_chain(chain[i])
    return chain


 


# Initial values __________________________________________________________________________
    # Observation value
y = 1.619888
sd = 0.43
    # Set pathes amd folder
template_idf_path = "test/campusbuilding_p10.idf"
eplus_basic_folder = "test/basic_files"
output_folder = "test/out"
result_CSV_path = 'test/chain_ep_p10.csv'
#sys.argv[1]

startvalue = [0.116, 0.046, 11.67, 12.43, 21, 24, 14.37, 0.56, 0.72, 2.65]
    # Set the range for proposal function
#low_limit = [0.01, 0.01, 0.1, 0.1, 1, 1, 17, 17, 1, 0.1, 0.5, 2]
#upper_limit = [0.5, 0.5, 8, 1, 100, 80, 28, 28, 50, 4, 0.99, 5]
#proposal_sd = [0.02, 0.015, 3.5, 3.5, 1.3, 1.3, 4, 0.15, 0.06, 0.15]

iterations = 5
totalarea = 10336.99
count = 1
markup = ['@@ROOF@@','@@WALL@@','@@EPD@@','@@LPD@@',
            '@@HSP@@','@@CSP@@','@@OCC@@','@@INF@@','@@Boiler@@','@@COP@@']


markup_value_pairs = generate_markup_value_pairs(markup,startvalue)

path = prepare_job_folders(output_folder, template_idf_path, eplus_basic_folder, markup_value_pairs)

# RUN!!
import time
start_time = time.time()
chain = run_metropolis_MCMC(path)
print "Simulation time:",time.time() - start_time, "seconds"

# Acceptance
import itertools
print "Acceptance: ",len(list(chain for chian,_ in itertools.groupby(chain))) / float(iterations)

# make CSV file
#make_csv(chain)
# ____________________________________________________________________________________
