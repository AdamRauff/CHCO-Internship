# Adam Rauff
# 3/25/2017
# This script is written in ordr to match VVCR values for patients in the VVMM database 
# the script looks for the MRN and other demographics of a patient in the VVMM database,
# determines if that patient's pressure data has already been analyzed, and outputs 
# another csv file with the patient identifier and VVCR value

import csv

# read in the VVMM database from the 2 csv files
DBname = 'VVMM_Database.csv' 
demoName = 'VVMM_demographics.csv'

# list that hold to the patient numbers from the VVMM Database csv file
PatNum = []

# list that holds onto the dates of birth from the VVMM database csv
DOB = []

# read in the DOBs and PatNums of the VVMM database csv file
with open(DBname, 'r') as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=',', quotechar='|')

  # store all data of csv in a list of lists
  Data = list(csv_reader)

for i, row in enumerate(Data):
  for j, column in enumerate(row):
    if j == 0 and i > 0:
      PatNum.append(Data[i][j])
    elif j == 1 and i > 0:
      DOB.append(Data[i][j])

  # print('PatNum: ',PatNum)
  # print('DOB: ', DOB)
VVMM_DB_List = [PatNum, DOB]

# declare lists of the columns of the demographics
demPtNo = []
demFstNam = []
demLstNam = []
demGender = []
demMRN = []
demDOB = []

# read in the info of the VVMM dmeographics csv file
with open(demoName, 'r') as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=',', quotechar='|')

  # store all data of csv in a list of lists
  Data = list(csv_reader)

for i, row in enumerate(Data):
  for j, column in enumerate(row):
    if j == 0 and i > 0:
      demPtNo.append(Data[i][j])
    if j == 1 and i > 0:
      demFstNam.append(Data[i][j])
    if j == 3 and i > 0:
      demLstNam.append(Data[i][j])
    if j == 4 and i > 0:
      demGender.append(Data[i][j])
    if j == 5 and i > 0:
      demMRN.append(Data[i][j])
    elif j == 6 and i > 0:
      demDOB.append(Data[i][j])

VVMM_Demo_List = [demPtNo, demFstNam, demLstNam, demGender, demMRN, demDOB]
