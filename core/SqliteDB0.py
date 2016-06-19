# -*- coding: utf-8 -*-
"""
Created on Wed May 27 16:27:36 2015

@author: tcolomb
"""

import os
#from Tkinter import *
from guidata.qt.QtGui import QFileDialog
import sqlite3 as lite#import numpy as np
from tracking3d.core import utilsForProcessing
#from itertools import islice
#import csv

 
#Class to read and write database ObjectiveDB.db3
class DBReader:
    def __init__(self):
        self.dBpath = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))+'\\ObjectiveDB.db3'
        self.dBKoalapath = r'C:\ProgramData\Koala\Koala2L.db3'
        #self.dBKoalaVersion = None
        self.MagnificationName = None
        self.NA = None
        self.ProductionNo = None
        self.configType = None
        self.choiceItems = None
        self.con = lite.connect(self.dBpath)
        self.ObjectiveNo = None
        
        self.DefaultDirectory = None
        self.TransmissionDHM = None
        self.CCDPixelSizeUM = None
        self.IPKoala = None
        self.dbKoalaExist = True
        self.ConfigROI = [0, 0, 1024, 1024]
        self.conKoala = None
        #connection to database and import parameters from the database
        with self.con:
            self.con.row_factory = lite.Row
            cur = self.con.cursor()
            cur.execute("SELECT * FROM ObjectiveTable")
            rows = cur.fetchall()
            self.choiceItems = []
            for row in rows:
                config = row["ProductionNo"]
                name = str(row["MagnificationName"])
                self.choiceItems.append((config, name))
            cur.execute("SELECT Value FROM DefaultValueTable WHERE Name='DefaultDirectory'")
            data = cur.fetchone()
            self.DefaultDirectory = data[0]
            
            cur.execute("SELECT Value FROM DefaultValueTable WHERE Name='TransmissionDHM'")
            data = cur.fetchone()
            self.TransmissionDHM = (int)(data[0])
            
            cur.execute("SELECT Value FROM DefaultValueTable WHERE Name='CCDPixelSizeUM'")
            data = cur.fetchone()
            self.CCDPixelSizeUM = (float)(data[0])
            
            cur.execute("SELECT Value FROM DefaultValueTable WHERE Name='IPKoala'")
            data = cur.fetchone()
            self.IPKoala = data[0]
            

#        Exemple code if we want to extract some data from koala database Koala2L.db3
#        dbKoalaExist=True        
        try:
            open(self.dBKoalapath)
        except:
            self.dbKoalaExist = False
            utilsForProcessing.ErrorMessage('Koala DB is unaccessible')
        if self.dbKoalaExist:
#            conKoala = None
            self.conKoala = lite.connect(self.dBKoalapath)
#            #version of Koala
#            with conKoala:
#                conKoala.row_factory = lite.Row
#                curKoala = conKoala.cursor()
#                curKoala.execute("SELECT * FROM VersionTable")
#                data = curKoala.fetchone()
#                self.dBKoalaVersion = data["Version"]
#
#                curKoala.execute("SELECT * FROM PeripheralTable WHERE PeripheralCodeKoala=1")
#                data = curKoala.fetchone()
#                maxtext = "MaxLimitSoftQC"
#                mintext = "MinLimitSoftQC"
#                if self.dBKoalaVersion>29:
#                    maxtext = "MaxCoderPos"
#                    mintext = "MinCoderPos"
#                self.OPLParams = [data["PhysDim"],data[maxtext],data[mintext]]

    def choiceItemDB(self): #return the list of available configuration to define Magnification button in main panel
        return self.choiceItems

    def ParamFromChoice(self):#import the parameters associated to a configuration
        with self.con:
            self.con.row_factory = lite.Row
            cur = self.con.cursor()
            textsearch = "SELECT * FROM ObjectiveTable WHERE ProductionNo = "+str(self.ProductionNo)
            cur.execute(textsearch)
            data = cur.fetchone()
            self.configType = data["ConfigType"]
            self.NA = data["NA"]
            self.ObjectiveNo = data["ObjectiveNo"]
        if self.dbKoalaExist and self.conKoala is not None:
            self.conKoala.row_factory = lite.Row
            curKoala = self.conKoala.cursor()
            textsearch = "SELECT * FROM ProductionTable WHERE ProductionNo = "+str(self.ProductionNo)
            curKoala.execute(textsearch)
            data = curKoala.fetchone()
            AlgoNo = (int)(data["AlgoNo"])
            textsearch = "SELECT * FROM AlgoTable WHERE AlgoNo = "+str(AlgoNo)
            curKoala.execute(textsearch)
            data = curKoala.fetchone()
            RoiX = data["RoiX"]
            RoiY = data["RoiY"]
            RoiWidth = data["RoiWidth"]
            RoiHeight = data["RoiHeight"]
            self.ConfigROI = [RoiX, RoiY, RoiWidth, RoiHeight]
    def DefineDefaultDirectory(self):#import Default directory
        path = str(QFileDialog.getExistingDirectory(None, "Choose Default Directory"))
        with self.con:
            self.con.row_factory = lite.Row
            cur = self.con.cursor()
            insert_path = "UPDATE DefaultValueTable SET Value=? WHERE Name=?"
            cur.execute(insert_path,(path, 'DefaultDirectory'))
            self.con.commit()