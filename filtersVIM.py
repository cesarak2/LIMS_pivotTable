# Author: Cesar Krischer, April 2021
# script to open chemistry datasets from LIMS, transform (pivot table) and 
# plot historical values for a given chemistrydatetime A combination of a date and a time. Attributes: ()

# Paths, alloy compositions and results are hidden on purpose. 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
import time
import csv
import os # changing folders and select files
import re # for seaching files in the system
import fnmatch # for finding files to import

from paths import * # separated files containing all directories (LIMSFolderPath, save_files_to_path)
#LIMSFolderPath
#save_files_to_path

def find(pattern, path):
    '''
    looks for all files in a folder that contain a certain patter;
    use a blank list to append all files to be read and worked with
    input: pattern as RE and path to scan files
    '''
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


def pivot_205_table(alloy_code205, element_list):
    '''
    Function to create a Pivot Table from the 205 report by Melt ID
    Input: 205 report, list of elements of interest 
    Output: data frame pivot table containing Melt ID, heatNo, each element Beg and End, 
    '''
    # list files containing last 50+ chemistries heats
    raw205 = pd.read_csv(alloy_code205) #reads the file, filter for element and position
    alloyCode = raw205['U Alloy Codeid'][0]
    currentVAlloyCodeElement = raw205[raw205.Paramid.isin(element_list)] # filters the 205 report for the desired elements
    currentVAlloyCodeElement = currentVAlloyCodeElement[(currentVAlloyCodeElement['Location Id'] == 'Beg') |
                                              (currentVAlloyCodeElement['Location Id'] == 'End')] # select BEG and END values only
    currentVAlloyCodeElement = currentVAlloyCodeElement.sort_values(by=['Paramid'])
    currentVAlloyCodeElement.reset_index()
    pivotElement = currentVAlloyCodeElement.pivot(index = 'Melt Id', columns = ('Paramid', 'Location Id'), values=['Conc %'])
    pivotElement = pivotElement['Conc %'] # drops the upper most columns layer names 
    pivotElement['Melt_id'] = np.arrayindex=list(pivotElement.index)
    pivotElement['HeatNo'] = pivotElement.apply(lambda row: re.search('(?<=[ABCWLabcwl])\d+', str(row['Melt_id'])).group(0), axis=1)

    return pivotElement.sort_values(by='HeatNo', ascending = True)


def mark_target_heats(pivot_table, target_heats_list, alloy_code, label_if_true='trial_heats', save_PT=False): #RR for revert & reactives
    '''
    Adds to a data frame containing a list of heat numbers a column corresponding to the matching heats
    Input: Data frame from pivot_205_table, target heats for matching and a label (to be reused in diff projects)
    Output: Data frame 205 with new column for matching heats
    '''
    pivot_table[label_if_true] = 0
    try:
        pivot_table.loc[(pivot_table.Melt_id.isin(target_heats_list)),label_if_true] = 1
    except:
        print('error while getting Heat No')
    if save_PT == True:
        try:
            pivot_table.to_csv(alloy_code + '_' + label_if_true + '.csv')
        except:
            print('cannot override file')
    return(pivot_table).reset_index() #not sorting to preserve chemistry trends


def print_chart_heats(pivot_table_with_marked_heats, element, alloy_code, save_fig = False, label_if_true='trial_heats'):
    '''
    Function to plot the data from a data frame containing the beg and end chemistries from elements.
    Input: data frame from function mark_target_heats, element to plot, alloy code (for figure title), 
        bool for saving or not the charts, a label to add to the title.
    Output: scatter plot for a certain element vs. heats.
    '''
    try: #some elements are only analyzed at the beg or end.
        plt.scatter(x = pivot_table_with_marked_heats[pivot_table_with_marked_heats[label_if_true] == 0].index, 
                    y = pivot_table_with_marked_heats[pivot_table_with_marked_heats[label_if_true] == 0][element]['Beg'], c='black')
        plt.scatter(x = pivot_table_with_marked_heats[pivot_table_with_marked_heats[label_if_true] == 1].index, 
                    y = pivot_table_with_marked_heats[pivot_table_with_marked_heats[label_if_true] == 1][element]['Beg'], c='blue')
    except:
        print('no beg')
    try: #some elements are only analyzed at the beg or end.
        plt.scatter(x = pivot_table_with_marked_heats[pivot_table_with_marked_heats[label_if_true] == 0].index, 
                    y = pivot_table_with_marked_heats[pivot_table_with_marked_heats[label_if_true] == 0][element]['End'], c='gray')
        plt.scatter(x = pivot_table_with_marked_heats[pivot_table_with_marked_heats[label_if_true] == 1].index, 
                    y = pivot_table_with_marked_heats[pivot_table_with_marked_heats[label_if_true] == 1][element]['End'], c='red')    
    except:
        print('no end')
    plt.title(element + ', ' + alloy_code + '. Black|Blue = beg, Gray|Red = end')
    plt.tight_layout()
    if save_fig == True:
        plt.savefig(label_if_true + ', ' + alloy_code + ', ' + element +'.png', dpi = 200, transparent = False)
        plt.close()
    else:
        plt.show()
        plt.close()
        

def wrap_up(alloy_code205, element_list,
            target_heats_list, alloy_code, label_if_true,
            element, save_fig=False, save_PT = False):
    '''
    Function to aggregate all other functions.
    Input: all of other function inputs 'pivot_205_table', 'mark_target_heats', 'print_chart_heats'
    Outputs: same of other functions.
    '''
    temp205 = pivot_205_table(alloy_code205, element_list)
    temp_pivot_table_with_marked_heats = mark_target_heats(temp205, target_heats_list, alloy_code, label_if_true, save_PT)
    print_chart_heats(temp_pivot_table_with_marked_heats, element, alloy_code, save_fig, label_if_true)
       

main_elements = ['Al', 'B', 'C', 'Co', 'Cr', 'Cu', 'Fe', 'Mn', 'Mo', 'Nb', 'Ni',
                 'O', 'P', 'S', 'Si', 'Ta', 'Te', 'Ti', 'V', 'W', 'Zr']
RR_heats = ['140B62873','196B62874','186B62875','065C62890','062C62803',
            '085C62804','677B62970','563B62971','677B62976','563B62977',
            '813B62972','2005B62973','814B62974','681B62975','085C62928',
            '094C62961','094C62962','094C62963','224W62910','224W62965',
            '186W62914','140W62951']
alloy_codes_ref = ['2005', '814','813','681','677','563','224','196','186','140','094','085','065','062']


# = = = = = = = = = = =
#   GET DIRECTORIES   =
# = = = = = = = = = = =

rootFolderProject = os.getcwd() # save original folder
#os.chdir(rootFolderProject)
#os.chdir('..')     #os.getcwd()    #os.listdir()
#print('Root files access: ' + LIMSFolderPath)
#print('Saving files to: ' + save_files_to_path)

# = = = = = = = = = = =
#     TREAT DATA      =
# = = = = = = = = = = =

'''
os.chdir('062V')
for element in main_elements:
    wrap_up(LIMSFolderPath + '\\062V.csv', main_elements,
            RR_heats, '062', 'RR',
            element, save_fig=True, save_PT = True)
os.chdir(rootFolderProject)
'''


target = '2005' # saving the alloy code in a variable for recycling 
os.chdir(rootFolderProject)
os.chdir(target + 'V')
for element in main_elements:
    wrap_up(LIMSFolderPath + '\\' + target + 'V.csv', main_elements,
            RR_heats, target, 'RR',
            element, save_fig=True, save_PT = True)
os.chdir(rootFolderProject)




''' Doing manually for each element for each alloy code
alloy_140PT = pivot_205_table(LIMSFolderPath + '\\140V.csv', main_elements)
pivot_table_with_marked_heats = mark_target_heats(alloy_140PT, RR_heats, 'RR')
print_chart_heats(pivot_table_with_marked_heats, 'Al', '140', False, 'RR')


wrap_up(LIMSFolderPath + '\\140V.csv', main_elements,
        RR_heats, '140', 'RR', 
            'Al', save_fig=True, save_PT = True)
'''


'''
#https://stackoverflow.com/users/97828/nadia-alramli (find function)
'''