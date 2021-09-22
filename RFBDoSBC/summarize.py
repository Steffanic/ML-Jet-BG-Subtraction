import os
from utility import log
import shutil as shu
import glob
import numpy as np
import pickle

BASE_DIR = ""

NICE_NAMES = {'random_state':'Random State',
'oob_score':'OOB Score',
'n_jobs':'Number of Jobs',
'n_estimators':'Number of Estimators',
'min_weight_fraction_leaf':'Minimum Weight Fraction Leaf',
'min_samples_split':'Minimum Samples Split',
'min_samples_leaf':'Minimum Samples Leaf',
'min_impurity_decrease':'Minimum Impurity Decrease',
'max_samples':'Max Samples',
'max_leaf_nodes':'Max Leaf Nodes',
'max_features':'Max Features',
'max_depth':'Max Depth',
'criterion':'Criterion',
'class_weight':'Class Weight'}

def read_info(filename):
    info_file = open(filename, "r")
    ...

def read_best_params(rad, ptm):
    best_params = pickle.load(open(BASE_DIR+f"/Objects/best_params_0.{rad}_{ptm}0.pickle", 'rb'))
    return best_params

def make_best_params_table(best_params, rad, ptm):
    columns = best_params.keys()
    columns = sorted(columns, reverse=True)
    result = [best_params[key] for key in columns]
    tab_frag = ''''''
    for i in range(len(columns)):
        tab_frag = tab_frag + f'''{NICE_NAMES[columns[i]]} & {result[i]} \\\ \n'''

    table = f'''Best Parameters for R=0.{rad} and pTmin={ptm}0\n\\begin{{tabular}}{{c  c}}\hline\nParameter Name & Value \\\ \hline \n''' +tab_frag + f'''\\end{{tabular}}\n\n\n'''
    print(table)
    return table
    
def latex_header(filename):
    with open(filename, 'w') as sfile:
        sfile.write(f'\\documentclass{{article}} \n\\usepackage[utf8]{{inputenc}} \n\\begin{{document}} \n\n')
        
def write_info():
    latex_header("Summary.latex")
    with open("Summary.latex", 'a') as sfile:
        for i in range(2, 7):
            for j in range(1,5):
                sfile.write(make_best_params_table(read_best_params(i,j),i,j))

def cd_summary_if_not_there():
    if(not os.path.isdir(BASE_DIR+"/Summary")):
        log("Please run build_file_structure before doing anything else!")
    if(os.getcwd()==BASE_DIR+"/Summary"):
        pass
    else:
        os.chdir(BASE_DIR+"/Summary")

def build_file_structure(base_dir):
    if(os.path.isdir(base_dir)):
        global BASE_DIR
        BASE_DIR=base_dir
        if(not os.path.isdir(BASE_DIR+"/Summary")):
            os.mkdir("Summary")
        os.chdir(BASE_DIR+"/Summary")
        for i in range(2,7):
            if(not os.path.isdir(BASE_DIR+f"/Summary/R=0.{i}")):
                os.mkdir(f"R=0.{i}")
            os.chdir(f"R=0.{i}")
            for j in range(1,5):
                if(not os.path.isdir(BASE_DIR+f"/Summary/R=0.{i}/pTmin={j}0")):
                    os.mkdir(f"pTmin={j}0")
            os.chdir("..")

def move_plots(rad, ptm):
    cd_summary_if_not_there()
    for fi in glob.glob(BASE_DIR+f"/Plots/R=0.{rad}/*{ptm}0.png"):

        fi_nm = fi.removeprefix(BASE_DIR+f"/Plots/R=0.{rad}/").removesuffix(f"_pTmin{ptm}0.png").removesuffix(f"_pthard={ptm}0.png")+".png"
        dest = BASE_DIR+f"/Summary/R=0.{rad}/pTmin={ptm}0/{fi_nm}"
        shu.copy2(fi, dest)

def move_oracles(rad, ptm):
    cd_summary_if_not_there()
    shu.copy2(BASE_DIR+f"/Oracles/Oracle_0.{rad}_{ptm}0.png",BASE_DIR+f"/Summary/R=0.{rad}/pTmin={ptm}0/Oracle.png")

if __name__=="__main__":
    build_file_structure("/mnt/d/Research/MLJetBGSub/ML-Jet-BG-Subtraction")
    for i in range(2, 7):
        for j in range(1, 5):
            move_plots(i, j)
            move_oracles(i,j)
            write_info()