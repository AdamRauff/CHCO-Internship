from readPdfFunc import *
import os
import numpy as np
# 
# x.append(3)
# print('x: ', x)
# import os
# 
# print('blah')
# for dirpath, dirnames, files in os.walk("."):
#   for filename in files:
#     print(os.path.join(dirpath, filename))

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
      # open the pdf file
      pdf_file= open(os.path.join(path, files),'rb')

      # instantiate a pdf file reader object (of the PyPDF2 module)
      pdf_reader = PyPDF2.PdfFileReader(pdf_file)

      # obtain number of pages
      num_pages = pdf_reader.numPages
      

      if num_pages > 0:
        # for page_no in range(num_pages):
        StrText = []
        with open(files.replace('pdf','txt'),'a') as txt_file:
          for i in range(num_pages):
            # store the text from each page as an entry in a list
            # this way the contents of each page can be called and 
            # processed separatly
            StrText.append(pdf_reader.getPage(i).extractText())
            txt_file.write(StrText[i])



