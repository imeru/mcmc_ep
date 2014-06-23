

# Import relevant modules
import numpy as np
import pandas as pd
import math
from scipy.stats import triang, norm
import time

# For Driver
import sys
import os
import shutil
import csv
from subprocess import call


# initial values
template_idf_path = "test/campusbuilding.idf"
eplus_basic_folder = "test/basic_files"
output_folder = "test/out"
#sys.argv[1]
count = 1
totalarea = 10336.99

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
    
    """
    #---------------------------------
    #to do: get the EUI from energyplus
    prediction = -0.2177929 - 0.2458677 * ROOF - 0.3842544 * WALL - 0.0049753 * WIN  \
        + 0.0176093 * SHGC + 0.014322 * EPD + 0.0125891 * LPD + 0.0518334 * HSP \
        +0.0026691 * CSP + 0.0003318 * OCC + 0.5664563* INF - 0.6163375 * Boiler \
        - 0.0425599 * COP
    #---------------------------------
    """
    chain = []
    for item in param:
        temp =[item]
        chain.append(temp)
    
    markup_values_pairs = dict(zip(['@@ROOF@@','@@WALL@@','@@WIN@@',
                                '@@SHGC@@','@@EPD@@','@@LPD@@',
                                '@@HSP@@','@@CSP@@','@@OCC@@',
                                '@@INF@@','@@Boiler@@','@@COP@@'], chain))

    # prepares jobs
    markup_value_pairs = generate_markup_value_pairs(markup_values_pairs, count)
    path = prepare_job_folders(output_folder, template_idf_path, eplus_basic_folder, markup_value_pairs)
    
    # Run energyplus and get eui value
    prediction = run_eplus(path)
         
    singlelikelihoods = norm.logpdf(x=y, loc=prediction, scale=sd)
    #sumll = sum(singlelikelihoods)
    return singlelikelihoods

#Posterior 
def posterior(param):
    return likelihood(param) + prior(param)



#Proposal function
def proposalfunction(param):
    return abs(np.random.normal(loc=param, scale=[0.02,0.015,0.43,0.1,7,4,1.3,1.3,6.75,0.19,0.06,0.167]))

#Metro-Polis MCMC
def run_metropolis_MCMC(startvalue, iterations):
    chain = [[0 for x in xrange(12)] for x in xrange(iterations)]
    chain[0] = startvalue
    for i in range(1,iterations):
        proposal = proposalfunction(chain[i-1])
        probab = min(1, math.exp(posterior(proposal) - posterior(chain[i-1])))
        if np.random.uniform() < probab:
            chain[i] = list(proposal)
        else:
            chain[i] = list(chain[i-1])
    return chain



#Run
# Observation value
y = 1.7
sd = 0.1
startvalue = [0.09667, 0.055, 2.792, 0.5, 22.8, 14.5, 21, 24, 23.4, 0.675, 0.72, 2.65]
iterations = 5


import time
start_time = time.time()
chain = run_metropolis_MCMC(startvalue, iterations)
print time.time() - start_time, "seconds"
print chain



# Execute Energyplys 


def replace_markup(line, markup_value_pairs):
    for markup in markup_value_pairs.keys():
        line = line.replace(markup, str(markup_value_pairs[markup]))
    return line


def generate_markup_value_pairs(markup_values_pairs, count):
    markup_value_pairs  = []
    for index in range(count):
        markup_value = {}
        for key in markup_values_pairs:
            markup_value[key] = markup_values_pairs[key].pop()
        markup_value_pairs.append(markup_value)
    return markup_value_pairs


def write_idf(template_path, output_path, markup_value_pairs):
    origin = open(template_path, 'r')
    new = open(output_path, 'w')
    for line in origin:
        replaced_line = replace_markup(line, markup_value_pairs)
        new.write(replaced_line)
    origin.close()
    new.close()

           
def copy_files(orig, dest):
    files = os.listdir(orig)
    for file_name in files:
        file_path = os.path.join(orig, file_name)
        shutil.copy(file_path, dest)

def prepare_job_folders(output_folder, template_idf_path,
                        eplus_basic_folder, markup_value_pairs):
    # check output path
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    pathes = []
    for index, markup_value_pair in enumerate(markup_value_pairs):
        path_to_write = output_folder + "/" + str(index)
        pathes.append(path_to_write)
        output_path = path_to_write + "/" + "in.idf"
        os.makedirs(path_to_write)
        copy_files(eplus_basic_folder, path_to_write)
        write_idf(template_idf_path, output_path, markup_value_pair)
        path = pathes[0]
    return path

    

def run_eplus(path):
    current_dir = os.getcwd()
    os.chdir(path)
    call(["EnergyPlus"])
    call(["ReadVarsESO", "my_results.rvi"])
    eui = get_eui(totalarea)
    os.chdir(current_dir)
    return eui 


def get_eui(totalarea):
    with open('eplusout.csv', 'rb') as f:
        reader = csv.reader(f)
        reader = list(reader)
        gas = float(reader[1][1])
        elec = float(reader[1][2])
        eui = (gas + elec) / totalarea / 10e8
    return eui

