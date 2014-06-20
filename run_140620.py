#Run
# Observation value
startvalues = [0.09667, 0.055, 2.792, 0.5, 22.8, 14.5, 21, 24, 23.4, 0.675, 0.72, 2.65]

# Driver
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
totalarea = 10336.99


# check output path
if os.path.exists(output_folder):
    shutil.rmtree(output_folder)

# prepares jobs
chain = []
for item in startvalues:
    temp =[item]
    chain.append(temp)

markup_values_pairs = dict(zip(['@@ROOF@@','@@WALL@@','@@WIN@@',
                                '@@SHGC@@','@@EPD@@','@@LPD@@',
                                '@@HSP@@','@@CSP@@','@@OCC@@',
                                '@@INF@@','@@Boiler@@','@@COP@@'], chain))
count = 1


# Execute Energyplys 


def replace_markup(line, markup_value_pairs):
    for markup in markup_value_pairs.keys():
        line = line.replace(markup, str(markup_value_pairs[markup]))
    return line


def generate_markup_value_pairs(markup_values_pairs, count) :
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

# prepares jobs
markup_value_pairs = generate_markup_value_pairs(markup_values_pairs, count)
path = prepare_job_folders(output_folder, template_idf_path, eplus_basic_folder, markup_value_pairs)


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


predict = run_eplus(path)
print predict