# 1/09/2017
# Adam Rauff

# the purpose of this script is to contain the variable declerations
# for the automatic pdf data extraction 

import numpy as np
import re
import csv

preNam = [] #list of patient names *** Not actually recorded in this script
preMRN = []
preCathDate = []
preFileNum = []
preDOB = []
preGender = []
preHeight = []
preWeight = []
preBSA = []

preLists = [preNam, preMRN, preFileNum, preDOB, preCathDate, preGender, \
    preHeight, preWeight, preBSA]

# RA = Room Air
# Hemodynamic values
preQpRA = [] # pulmonary flow normalized with BSA (L/min/m^2)
preQsRA = [] # systemic flow normalized (L/min/m^2)
preRpRA = [] # pulmonary resistance
preRsRA = [] # Systemic Resistance
preQp_Qs = [] # ratio
preRp_Rs = [] # res ratio

preHRRA = [] # Heart rate
preVO2 = []
preHgRA = []
preCITher = [] # Cardiac Index from thermodilution

# "Valused used in Fick's method calculation from cath report
preMVPO2RA = []
prePAPO2RA = []
prePVPO2RA = []
preSAPO2RA = []
preMV_satRA = []
prePAsatRA = []
preMMPA_RA = []
preMSys_RA = []
preHB_RA = []
prePVsat_RA = []
preSAsat_RA = []
preWedge_RA = []
preMRA_RA = []

# compose list of lists for more concise notation
preRA = [preQpRA, preQsRA, preRpRA, preRsRA, preQp_Qs, preRp_Rs, preHRRA, \
    preVO2, preHgRA, preCITher, preMV_satRA, prePAsatRA, preMMPA_RA, \
    preMSys_RA, preHB_RA, prePVsat_RA, preSAsat_RA, preWedge_RA, preMRA_RA, \
    preMVPO2RA, prePAPO2RA, prePVPO2RA, preSAPO2RA]

totPreList = preLists + preRA
# recall python is object oriented, so the above line is simply a pointer 
# that points to a bin. Thus, preQpRA.append('x') is equivalent to 
# preRA[0].append('x')

chkFeatList = ['Height' ,'Weight', 'BSA', 'Qp', 'Qs', 'Rp', 'Rs', \
    'Qp/Qs', 'Rp/Rs', 'HR', 'VO2', 'Hg', 'CITher', 'MV sat', 'PA sat', \
    'MMPA', 'MSys', 'HB', 'PV sat', 'SA sat', 'Wedge', 'MRA', \
    'MV PO2', 'PA PO2', 'PV PO2', 'SA PO2']

col = list(range(6,32))

featList = ['Name', 'MRN', 'CathNum', 'DOB', 'CathDate', 'Gender'] + chkFeatList

# ----------------------------------------------------------
