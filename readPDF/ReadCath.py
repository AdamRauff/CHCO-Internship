#!/bin/python
# Dr. Kendal Hunter
# Adam Rauff
# 12/27/2016

# The purpose of this script is to automatically read specific values
# from cath reports of PH patients for research purposes

# pyPDF is no longer supported, but PyPDF2 
# was forked from the original project and is supported
import PyPDF2
import os
import re # regular expressions
from readPdfFunc import * 
# The above is a written python script that hold function
# called on in this code
import csv
from varDec import *

# turn true to display printing statements used for debugging
debugger2 = False
# -----------------------------------------------------
# may be a good idea to add argparse to give a working directory 
# that contains the pdf files when called from command line. Or, can call
# on some file navigation GUI (like windows explorer)
# -------------------------------------------------------
# Get all the files in the cwd
dirs = os.listdir(".")

# loop through all pdf files present in the cwd
# and all subdirectories
for path, dirnames, filenames in os.walk("."):
  for files in filenames:
    
    # make files all lower case
    filesL = files.lower()

    # Only process pdf files 
    if filesL.endswith(".pdf"):
    
      print('pdf file name: ',files)
    
      # open the pdf file
      pdf_file= open(os.path.join(path, files),'rb')

      # instantiate a pdf file reader object (of the PyPDF2 module)
      pdf_reader = PyPDF2.PdfFileReader(pdf_file)

      # obtain number of pages
      num_pages = pdf_reader.numPages

      if debugger:  
        print('# pages: ',num_pages)

      # txt_file = open(files.replace('pdf','txt'),'wb')

      if num_pages > 0:
        # for page_no in range(num_pages):
        StrText = []
        for i in range(num_pages):
          # store the text from each page as an entry in a list
          # this way the contents of each page can be called and 
          # processed separatly
          StrText.append(pdf_reader.getPage(i).extractText())

        # declare variable for string of first page
        pyStrText = StrText[0]

        # find birthdate
        PatNam = MasterFind(pyStrText, '(?<=Laboratory)\w+\W?\s?\w*\W?\s?\w*(?=MRN)' \
            , -1, 'Laboratory')
        PatNam = PatNam.replace(',', ' - ')
        print('Name: ', PatNam)
        DOB = str(findDate(pyStrText, 'Birth'))

        # find MRN
        MRN = str(findMRN(pyStrText))

        # find Cath date
        CathDate = str(findDate(pyStrText, 'Cath'))

        # find cath # (filename, i.e HA0000XX)
        CathNum = str(findCathNum(pyStrText))

        # find gender
        Gender = findGender(pyStrText)

        # find the Height, Weight, and BSA
        Height, Weight, BSA = findHWB(pyStrText)

        # find RA Hemodynamics
        RA  = findHemo(StrText,'21%')

        # print('RA: ',RA)
        
        # create easy-to-use lists
        pdfValues = [Height, Weight, BSA]
        prePdfV = [preHeight, preWeight, preBSA]
        PatVals = [PatNam, MRN, CathNum, DOB, CathDate, Gender]
        TotRowVals = [PatNam, MRN, CathNum, DOB, CathDate, Gender, Height, Weight, BSA]
        CompRowVals = TotRowVals + RA

        # for feat, val in zip (featList, CompRowVals):
        #   print(feat+': '+val)
        # store text in bits so it can be easilt written to a text file
        # text = pyStrText.encode('UTF-8')

        # Write entries to the text files (example)
        # txt_file.write(text)
        # ----------------------------------------------- 
        # check to see if patient pdf information has already been analyzed
        if csvExists:
          
          # intialize all flags
          MRNFlag = False
          FileNumFlag = False
          CathDateFlag = False
          DOBFlag = False
          MRNinds = []
          FileNinds = []
          CathDinds = []
          DOBinds = []

          # determine if current patient has alrady been analyzed
          # compare MRNs with all previous MRNs
          if MRN in preMRN:
            MRNFlag = True
            MRNinds = [i for i,s in enumerate(preMRN) if MRN == preMRN[i]]
            # print('MRNinds: ', MRNinds)

          # compare file number (HA000XXX)
          if CathNum in preFileNum:
            FileNumFlag = True
            FileNinds = [i for i,s in enumerate(preFileNum) if \
                CathNum == preFileNum[i]]
            # print('FileNinds: ', FileNinds)

          # compare Cath date with all previous Cath dates
          if CathDate in preCathDate:
            CathDateFlag = True
            CathDinds = [i for i,s in enumerate(preCathDate) \
                if CathDate == preCathDate[i]]
            # print('CathDinds: ',CathDinds)

          # compare date of birth
          if DOB in preDOB:
            DOBFlag = True
            DOBinds = [i for i,s in enumerate(preDOB) if DOB == preDOB[i]]
            # print('DOBinds: ',DOBinds)

          # if patient has already been analyzed
          if MRNFlag == True and FileNumFlag == True and CathDateFlag == True \
              and DOBFlag == True:
            # find the index of all XXXinds that match
            SameIndx = [i for i in DOBinds if i in CathDinds and i in FileNinds \
                and i in MRNinds]
            if len(SameIndx) > 0:
              print('SameIndx: ', SameIndx)
              print('Patient of this file has already been analyzed!')
             
              # store the latest index in an integer rather than a list
              SameIndx = SameIndx[-1]
              tempGen = preGender[SameIndx]
              ChkVars = pdfValues + RA
              preLVals2Chk = prePdfV + preRA 
              
              tempVars = []
              for s in preLVals2Chk:
                tempVars.append(s[SameIndx])

              # print('tempVars: ', tempVars)

              if not tempGen == Gender:
                print('Difference in files: Gender!!!!!!!') 
                print(' ')
                if len(tempGen) == 1 and len(Gender) > 1:
                  # write over current value in csv file
                  # identify row of patient SameIndx
                  # write over it
                  AppendCSV(csvFileName, CompRowVals)
              # print('len(tempVars): ',len(tempVars))
              # print('len(ChkVars): ',len(ChkVars))
              # print('len(chkFeatList): ',len(chkFeatList))
              # print('len(col): ',len(col))
              # print('len(preLVals2Chk): ',len(preLVals2Chk))

              for itr, (CS, PD, F, C, L) in enumerate(zip(tempVars, ChkVars, \
                  chkFeatList, col, preLVals2Chk)):
               
                # if the stored value in the csv file does NOT 
                # match value obtained from current pdf
                if CS.replace(' ', '') == '' and not CS == PD:

                  print('filename: ', files)
                  print('Difference in files: ', F)
                  print('itr: ',itr)
                  print('CSV: ',CS)
                  print('Current PDF: ',PD)
                  # print(tempVars)
                  # print(ChkVars)
                  print('preVLS2: ',preLVals2Chk)
                  print('totPreL: ',totPreList)
                  # update pre-list
                  L[SameIndx] = PD
                  
                  print('Overwrite Current Patient CSV')
                  totPreList = OverWriteCSV(csvFileName, SameIndx, \
                      C, PD, CS, totPreList)
                  # if the stored value is not blank, and the new value
                  # is also not blank, and they are different, then add new patient
                elif len(CS.replace(' ','')) > 0 and len(PD.replace(' ','')) \
                    > 0 and not CS == PD:
                  # append a new patient with current information
                  # allow user to decipher when visualizing csv
                  print('Appending Difference')
                  AppendCSV(csvFileName, CompRowVals)
                
                  for feat, List in zip(CompRowVals, totPreList):
                    List.append(feat)
                  break
                    
                  # print(totPreList)
            # if patient has not already been analyzed
          else:
            print('Appending New Patient')  
            # append row to csv
            AppendCSV(csvFileName, CompRowVals)

            # add the current patient to the list of previous patients info 
            for feat, List in zip(CompRowVals, totPreList):
              List.append(feat)

            # print(totPreList)
            csv_file.close()
        else:
          # create a csv file named outcomes
          print('Creating New CSV')
          csv_file = open('Outcomes.csv', 'w')
          csv_writer = csv.writer(csv_file, delimiter=',', quotechar='|', \
              quoting=csv.QUOTE_MINIMAL)

          # first row is the titles
          csv_writer.writerow(['Patient Name', 'MRN', 'Cath #', \
              'DOB (mm/dd/yyyy)', 'Cath Date (mm/dd/yyyy)', 'Gender', \
              'Height (cm)' ,'Weight (Kg)', 'BSA (m^2)', 'Qp (L/min/m^2)', \
              'Qs (L/min/m^2)', 'Rp (units x m^2)', 'Rs (units x m^2)', \
              'Qp/Qs', 'Rp/Rs', 'HR (bpm)', 'VO2 (ml/min/m^2)', 'Hg (gm/dL)', \
              'Thermo CO' ,'MV sat', 'PA sat' ,'Mean MPA','Mean Sys', \
              'HB', 'PV sat', 'SA sat', 'Wedge', 'Mean RA' ,'MV PO2', \
              'PA PO2', 'PV PO2', 'SA PO2'])

          # prepare row to write  
          row2Writ = []
          for i in CompRowVals:
            row2Writ.append(i)
          csv_writer.writerow(row2Writ)

          # add patient info to list of previous patient info recorder
          for feat, List in zip(CompRowVals, totPreList):
            List.append(feat)
           
          # print(totPreList)
          # close csv writer
          csv_file.close()
          csvExists = True
          csvFileName = 'Outcomes.csv'

        pdf_file.close()
