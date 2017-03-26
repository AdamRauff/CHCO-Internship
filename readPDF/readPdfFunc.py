# Adam Rauff
# 12/28/2016
# pulmonary hypertension research - EPIC Outcomes

# the functions defined in this section allow
# the ReadCath code to be conscise and easy to read
# also, this script contains varaible declerations

import PyPDF2
import os
import re
import numpy as np
import csv
from varDec import *
# turn true to display printing statements used for debugging
debugger = False

# ---------------------------------------------------------
# Variable Declerations
# Get all the files in the cwd
dirs = os.listdir(".")

csvExists = False
# check if csv file exists in directory
for files in dirs:
  filesL = files.lower()
  if filesL.endswith('.csv'):
    if debugger:
      print('A csv file exists')
    csvExists = True
    csvFileName = files
    # OverWriteCSV('Outcomes.csv', 0, 0, 0, 0, 0)

# loop through csv file and record patient info to ensure 
# identical files are not analyzed twice
if csvExists:
  with open(csvFileName) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',', quotechar='|')

    # begin by updating the information of previous patient list
    for i, row in enumerate(csv_reader):
      if i > 0:
        for ColVal, preList in zip(row, totPreList):
          preList.append(str(ColVal))
      if debugger:
        print('row: ',row)
        print('i: ', i)

# ----------------------------------------------------------
# Date of Birth (DOB)
def findDate(pyStrText, Type):
  if Type == 'Birth':
    DOBind = pyStrText.find("Birth Date:") 
          
    if not DOBind == -1:
            
      if debugger:
        print('DOB index: ',str(DOBind))
        print('Index call: ', pyStrText[DOBind:(DOBind+25)])

      # establish a string that includes the Date of Birth
      DOBstr = pyStrText[DOBind:(DOBind+25)]

      # construct regualr expression object to detect DOB as mm/dd/yyyy
      DOBre = re.compile('\d+/\d+/\d+')
      
      # employ regx on string
      DOB = DOBre.findall(DOBstr)
      
      if debugger:
        if len(DOB) > 1:
          print('More than one DOB was found!')
        else:
          print('DOB is: ', DOB[0])

      return DOB[0]
    else:
      print('DOB was not obtained!')
      DOB = ' '
      return DOB

  elif Type == 'Cath':
    Cathind = pyStrText.find("Cath Date:") 
          
    if not Cathind == -1:
            
      if debugger:
        print('Cath index: ',str(Cathind))
        print('Index call: ', pyStrText[Cathind:(Cathind+25)])

      # establish a string that includes the Date of Birth
      Cathstr = pyStrText[Cathind:(Cathind+25)]

      # construct regualr expression object to detect DOB as mm/dd/yyyy
      Cathre = re.compile('\d+/\d+/\d+')
      
      # employ regx on string
      Cath = Cathre.findall(Cathstr)
      
      if debugger:
        if len(Cath) > 1:
          print('More than one Cath Date  was found!')
        else:
          print('Cath Date: ', Cath[0])

      return Cath[0]
    else:
      print('Cath Date was not obtained!')
      Cath = ' '
      return Cath

# ----------------------------------------------------------------
def findMRN(pyStrText):
  MRNind = pyStrText.find('MRN:')

  if not MRNind == -1:
          
    if debugger:
      print('MRN index: ',str(MRNind))
      print('Index call: ', pyStrText[MRNind:(MRNind+15)])

    # establish a string that includes the Date of Birth
    MRNstr = pyStrText[MRNind:(MRNind+15)]

    # construct regualr expression object to detect DOB as mm/dd/yyyy
    # MRNre = re.compile('\s+')
    
    # employ regx on string
    
    # MRN = MRNre.split(MRNstr)
    
    # MRN = MRN[1]

    # if debugger:
    #   print('MRN after split regexp: ',MRN)

    MRNre = re.compile('\d+')
    
    MRN = MRNre.findall(MRNstr)

    if debugger:
      print('MRN afer find numerics regexp: ',MRN)
    if len(MRN[0]) > 3:
      return MRN[0] 
    else:
      print('MRN  was not obtained!')
      return ' '
  else:
    print('MRN  was not obtained!')
    return ' ' 

# -----------------------------------------------------------
# Cath file number
# e.g Cath #: HA000014
def findCathNum(pyStrText):
  CathNumind = pyStrText.find('Cath #:')

  if not CathNumind == -1:
          
    if debugger:
      print('Cath # index: ',str(CathNumind))
      print('Index call: ', pyStrText[CathNumind:(CathNumind+20)])

    # establish a string that includes the Date of Birth
    CathNumstr = pyStrText[CathNumind:(CathNumind+20)]

    # construct regualr expression object to detect DOB as mm/dd/yyyy
    CathNumre = re.compile('\s+')
    
    # employ regx on string
    CathN = CathNumre.split(CathNumstr)
    
    CathN = CathN[2]

    if debugger:
      print('Cath # after split regexp: ',CathN)

    CathNumre = re.compile('\w\w\d+')
    
    CathN = CathNumre.findall(CathN)

    if debugger:
      print('CathN afer find numerics regexp: ',CathN)

    return CathN[0] 
  else:
    print('CathN  was not obtained!')
    CathN = ' '
    return CathN

# ----------------------------------------------------------
# Gender
# e.g Male or Female
def findGender(pyStrText):
  Gendind = pyStrText.find('Gender:')

  if not Gendind == -1:
          
    if debugger:
      print('Gender index: ',str(Gendind))
      print('Index call: ', pyStrText[Gendind:(Gendind+17)])

    # establish a string that includes the Date of Birth
    Gendstr = pyStrText[Gendind:(Gendind+17)]
    
    # if male:
    Male = Gendstr.find('Male')

    # if female:
    Female = Gendstr.find('Female')
    
    # 
    if Female == -1 and not Male == -1:
      return 'Male'
    elif Male == -1 and not Female == -1:
      return 'Female'
    # otherwise
    else:
      if debugger:
        print('The gender was not obtained!')
      return ' '
# ----------------------------------------------------------
# Height, Weight, BSA

def findHWB(pyStrText):
  Hind = pyStrText.find('Height:')
  if not Hind == -1:
    Hre = re.compile('\d+\.?\d*')
    Hstr = pyStrText[Hind:Hind+15]
    if debugger:
      print('Hstr: ',Hstr)
    Height = Hre.findall(Hstr)
    if debugger:
      print('Height after regexp: ',Height)
    if len(Height) == 0:
      Height = ' '
    else:
      Height = Height[0]
      chk = Height.lstrip('0')
      chk = chk.rstrip('0')
      if chk[0] == '.' or len(chk) == 0:
        Height = ' '
  else:
    Height = ' '

  Wind = pyStrText.find('Weight:')
  if not Wind == -1:
    Wre = re.compile('\d+\.?\d*')
    Wstr = pyStrText[Wind:Wind+15]
    Weight = Wre.findall(Wstr)
    if debugger:
      print('Weight after regexp: ',Weight)
    if len(Weight) == 0:
      Weight = ' '
    else:
      Weight = Weight[0]
      chk = Weight.lstrip('0')
      chk = chk.rstrip('0')
      if chk[0] == '.' or len(chk) < 1:
        Weight = ' '
  else:
    Weight = ' '

  Bind = pyStrText.find('BSA')
  if not Bind == -1:
    Bre = re.compile('\d+\.?\d*')
    Bstr = pyStrText[Bind:Bind+12]
    BSA = Bre.findall(Bstr)
    if debugger:
      print('BSA after regexp: ',BSA)
    if len(BSA) == 0:
      BSA = ' '
    else:
      BSA = BSA[0]
      chk = BSA.lstrip('0')
      chk = chk.rstrip('0')
      if chk[0] == '.' and len(chk) < 1:
        BSA = ' '
  else:
    BSA = ' '

  return Height, Weight, BSA

def Null():
  return ' '

def MasterFind(pyString, regexpStr, numIdxs, findStr):
  FirstInd = pyString.find(findStr)
  # print('MF First Ind: ', FirstInd)
  if not FirstInd == -1:
    rexp = re.compile(regexpStr)
    if numIdxs == -1:
      myStr = pyString
    elif FirstInd+numIdxs < len(pyString):
      myStr = pyString[FirstInd:(FirstInd+numIdxs)]
    else:
      myStr = pyString[FirstInd:]
    
    # print('MF myStr: ',myStr)
    if len(myStr) > 0:
      Value = rexp.findall(myStr)
      if len(Value) > 0:
        # print('MF after regexp: ',Value)
        return Value[0]
      else:
        return ' '
    else:
      return ' '
  else:
    return ' '

def OverWriteCSV(filename, Indx, col, pdfVal, CSVval, LofL):
  with open(filename, 'r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',', quotechar='|')

    # store all data of csv in a list of lists
    Data = list(csv_reader)

  for i, row in enumerate(Data):
    for j, column in enumerate(row):
      print('i, j: ', i, j)
      print('column: ', column)
      if i > 0 and (Indx+1) == i:
        if j == col and Data[i][j] == CSVval:
          print('col: ',col)
          print('pdfVal: ', pdfVal)
          # print(len(LofL))
          # print(len(LofL[i]))
          # print(len(LofL[i-1]))
          # print('LofL: ',LofL)
          # print('oldpreVal: ',LofL[j][i-1])
          Data[i][j] = pdfVal
          LofL[j][i-1] = pdfVal

  with open(filename, 'w') as csv_out:
    csv_writer = csv.writer(csv_out, delimiter=',', quotechar='|', \
        quoting=csv.QUOTE_MINIMAL)

    csv_writer.writerows(Data)
  #update Pre-lists
  return LofL
def AppendCSV(filename, CompRow):
  # open file in appending mode
  with open(filename, 'a') as csv_file:

    # instantiate writer object
    csv_writer = csv.writer(csv_file, delimiter=',', quotechar='|', \
        quoting=csv.QUOTE_MINIMAL)
    
    csv_writer.writerow(CompRow)

def findHemo(pyList, challStr):
  # convert pyStr to a string rather than a list
  pyStr = ''
  for i in pyList:
    pyStr = pyStr + str(i)
  
  # finding resistance and flow measurements -----------------------
  # find index of inspired O2: xx%
  ind = pyStr.find('Inspired O2: '+challStr)
  # print('ind: ', ind)
  
  if not ind == -1: 
    # initialize variables
    Qp = ' '
    Qs = ' '
    Rp = ' '
    Rs =  ' '
    QpQs = ' '
    RpRs = ' '
    HR = ' '
    VO2 = ' '
    Hg = ' '
    CITher = ' '
    MV_sat = ' '
    PA_sat = ' '
    MMPA = ' '
    MSys = ' '
    HB = ' '
    PV_sat = ' '
    SA_sat = ' '
    Wedge = ' '
    MRA = ' '
    MVPO2 = ' '
    PAPO2 = ' '
    PVPO2 = ' '
    SAPO2 = ' '

    # ----------------------------------------------------------------
    # find all instances where 'Qp =' is found
    QpL = np.array([m.start() for m in re.finditer('Qp =', pyStr)])
    # print('QpL: ',QpL)
    
    # find the PO2 values and values used in the calculations
    O2CL = np.array([m.start() for m in re.finditer('O2 consumption =', pyStr)])
    # print('O2CL: ',O2CL)

    if  len(QpL) > 0: 
    # find instance that is closest to ind and smaller
      diffAr = ind - QpL
      # print('diffAr: ',diffAr)
      
      # obtain positive number only
      diffAr2 = diffAr[diffAr > 0.0]

      if len(diffAr2) > 0:
          
        Qpind = np.argmin(diffAr2)
        # print('Qpind: ',Qpind)

        HemStr = pyStr[QpL[Qpind]:ind]

        # print(HemStr)
        # use MasterFind function to obtain hemodynamic values
        Qp = MasterFind(HemStr, '(?<=\s\()\d+\.\d+(?=\sL/min/m)', 33, 'Qp =')
        Qs = MasterFind(HemStr, '(?<=\s\()\d+\.\d+(?=\sL/min/m)', 33, 'Qs =')
        Rp = MasterFind(HemStr, '(?<=\s\()\d+\.\d+(?=\sunits\sx)', 34, 'Rp =')
        Rs =  MasterFind(HemStr, '(?<=\s\()\d+\.\d+(?=\sunits\sx)', 34, 'Rs =')
        QpQs = MasterFind(HemStr, '(?<=\s)\d+\.\d+\s+\:\s+\d(?=\s)', 22, 'Qp/Qs =')
        # strip whitespace from string
        if len(QpQs) > 1:
          QpQs = QpQs.replace(' ', '')
        RpRs = MasterFind(HemStr, '(?<=\=\s)\d+\.?\d*', 15, 'Rp/Rs =')
        HR = MasterFind(HemStr, '(?<=\s)\d+(?=\s+bpm)', 23, 'Heart Rate: ')
        VO2 = MasterFind(HemStr, '(?<=\s)\d+\.?\d*(?=\sml/min/m)', 23, 'VO2:')
        Hg = MasterFind(HemStr, '(?<=\s)\d+\.?\d*(?=\s\w+/\w+)', 25, 'Hemoglobin:')

    if  len(O2CL) > 0 and not ind == -1: 
      # find instance that is closest to ind and smaller
      diffO = O2CL - ind
      # print('diffO: ',diffO)
      
      # obtain positive number only
      diffO2 = diffO[diffO > 0.0]

      if len(diffO2) > 0:
          
        O2ind = np.argmin(diffO2)
        # print('Qpind: ',Qpind)

        O2Str = pyStr[ind:O2CL[O2ind]]

        # print(O2Str)
        # use MasterFind function to obtain hemodynamic values
        CITher = MasterFind(O2Str, '(?<=\:\s)\d+\.?\d*(?=\sL/min/m)', 25, 'CI:')
        MV_sat = MasterFind(O2Str, '(?<=MV\ssat\s\=\s)\d+\.?\d*(?=\s)', -1, 'MV sat')
        PA_sat = MasterFind(O2Str, '(?<=PA\ssat\s\=\s)\d+\.?\d*(?=\s)', -1, 'PA sat')
        MMPA = MasterFind(O2Str, '(?<=Mean\sMPA\s\=\s)\d+\.?\d*(?=\s)', -1, 'Mean MPA')
        MSys = MasterFind(O2Str, '(?<=Mean\sSys\s\=\s)\d+\.?\d*(?=\s)', -1, 'Mean Sys')
        HB  = MasterFind(O2Str, '(?<=HB\s\=\s)\d+\.?\d*(?=\s)', -1, 'HB =')
        PV_sat = MasterFind(O2Str, '(?<=PV\ssat\s\=\s)\d+\.?\d*(?=\s)', -1, 'PV sat')
        SA_sat = MasterFind(O2Str, '(?<=SA\ssat\s\=\s)\d+\.?\d*(?=\s)', -1, 'SA sat')
        Wedge = MasterFind(O2Str, '(?<=Wedge\s\=\s)\d+\.?\d*(?=\s)', -1, 'Wedge')
        MRA = MasterFind(O2Str, '(?<=Mean\sRA\s\=\s)\d+\.?\d*(?=\s)', -1, 'Mean RA')
        MVPO2 = MasterFind(O2Str, '(?<=MV\sPO2\s\=\s)\d+\.?\d*(?=\s)', -1, 'MV PO2')
        PAPO2 = MasterFind(O2Str, '(?<=PA\sPO2\s\=\s)\d+\.?\d*(?=\s)', -1, 'PA PO2')
        PVPO2 = MasterFind(O2Str, '(?<=PV\sPO2\s\=\s)\d+\.?\d*(?=\s)', -1, 'PV PO2')
        SAPO2 = MasterFind(O2Str, '(?<=SA\sPO2\s\=\s)\d+\.?\d*', -1, 'SA PO2')
    
    CurHem = [Qp, Qs, Rp, Rs, QpQs, RpRs, HR, VO2, Hg, CITher, \
        MV_sat, PA_sat, MMPA, MSys, HB, PV_sat, SA_sat, Wedge, \
        MRA, MVPO2, PAPO2, PVPO2, SAPO2]
    for i, val in enumerate(CurHem):
      CurHem[i] = val.replace(' ','')
    # print('CurHem: ',CurHem)

    return CurHem
  else:
    CurHem = []
    for i in range(23):
      CurHem.append('')
    return CurHem
