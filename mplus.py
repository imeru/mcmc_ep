from multiprocessing import Pool, cpu_count
from subprocess import call
import os
import shutil

def run_eplus(path):
    current_dir = os.getcwd()
    os.chdir(path)
    call(["EnergyPlus"])
    call(["ReadVarsESO", "my_results.rvi"])
    os.chdir(current_dir)

def run_eplus_multi(pathes):
    p = Pool(cpu_count())
    result = p.map(run_eplus, pathes)
    return result