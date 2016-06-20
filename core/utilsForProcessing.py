# -*- coding: utf-8 -*-
"""
Created on Tue Aug 05 14:27:06 2014

@author: tcolomb
"""

import os
from Tkinter import *
import csv
from guidata.qt.QtGui import QFileDialog, QMessageBox
import numpy as np
from os import listdir
from os.path import isfile, join



#%% MessageBox to take password entry
class takeInput(object):
    def __init__(self, requestMessage, boolTextOrNumber, defaultText, hideText):
        self.root = Tk()
        self.string = ''
        self.frame = Frame(self.root)
        self.frame.pack()        
        self.acceptInput(requestMessage, defaultText, hideText)
        

    def acceptInput(self, requestMessage, defaultText, hideText):
        r = self.frame
        k = Label(r, text=requestMessage)
        k.pack(side='left')
        self.e = Entry(r, text='Name')
        if hideText:
            self.e["show"] = "*"
        self.e.pack(side='left')
        self.e.insert(0, defaultText)
        self.e.focus_set()
        b = Button(r, text='okay', command=self.gettext)
        b.pack(side='right')

    def gettext(self):
        self.string = self.e.get()
        self.root.destroy()

    def getString(self):
        return self.string

    def waitForInput(self):
        self.root.mainloop()


def getEntry(requestMessage, boolTextOrNumber, defaultText, hideText):
    msgBox = takeInput(requestMessage, boolTextOrNumber, defaultText, hideText)
    msgBox.waitForInput()
    if boolTextOrNumber: #True=text, False=Number
        return msgBox.getString()
    else:     
        return int(float(msgBox.getString()))
def getPassword(requestMessage, defaultText):
    msgBox = takeInput(requestMessage, True, defaultText, True)
    msgBox.waitForInput()

    return msgBox.getString()

  
## Directory and file
def OpenTxtFile(text, path):
    filename = QFileDialog.getOpenFileName(None, text, path, filter="txt (*.txt *.)")
    return filename

    
def CreateDirectory(directoryPath, directoryName):
    if not os.path.exists(directoryPath+'\\'+directoryName):
        os.makedirs(directoryPath+'\\'+directoryName)
    return directoryPath + '\\' + directoryName
        
def DeleteAllFilesInDirectory(directoryPath):
    filelist = [f for f in os.listdir(directoryPath) if f.endswith(".bin")]
    for f in filelist:
        os.remove(directoryPath+'\\'+f)

def FileExists(fname, extension):
    return os.path.isfile(fname) and fname.endswith(extension)
    
def FindFileInDirectory(directory, extension):
    onlyfiles = [f for f in listdir(directory) if isfile(join(directory, f))]
    return onlyfiles
def ErrorMessage(message):
    QMessageBox.warning(None, 'Error', message)

def Log(message):
    print message

def DisplayImage(img):
    plt.imshow(img)
    plt.gray()
    plt.show()

def SaveParamsFile(ColumnTitles, Values, fname):
    if np.size(ColumnTitles) != np.size(Values):
        ErrorMessage("Not possible to save because ColumnTiltes and Values are not the same size")
    else:
        f = open(fname, 'w')
        columntext = ''
        valuetext = ''
        for k in range(np.size(ColumnTitles)-1):
            columntext += ColumnTitles[k]+'\t'
            valuetext += str(Values[k])+'\t'
        columntext += ColumnTitles[-1]+'\n'
        valuetext += str(Values[-1])+'\n'
        f.writelines(columntext)
        f.writelines(valuetext)
        f.close()
#    
#def ReadParmsFile(fname):
#    with open(fname,'r') as f:
#        reader=csv.reader(f)
#        Values=[]
#        for row in islice(reader,1,None):
#            line=row[0].split()
#            for v in line:
#                Values.append((float)(v))
#    return Values           
        
def SaveParamsFileByLine(Data, fname):
    fname = QFileDialog.getSaveFileName(None, "Save file", fname)
    f = open(fname, 'w')
    for info in Data:
        f.writelines(info[0]+'\t'+str(info[1])+'\t'+info[2]+'\n')
    f.close()
    
def ReadParamsFromLine(fname):
    Data = []
    with open(fname, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            #line=row[0].split()
            Data.append((row[0], row[1], row[2]))
    return Data