# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 21:47:31 2016

@author: zerodeku
"""

import os
from guidata.qt.QtGui import QFileDialog
import sqlite3 as lite
from tracking3d.core import utilsForProcessing

class DBReader:
    def __init__(self):
        # init parameters
        self.dBpath = os.path.abspath(os.path.join(os.path.dirname(__file__), 
                                os.pardir)) + '\\ObjectiveDB.db3'
        self.dBKoalapath = r'C:\ProgramData\Koala\Koala2L.db3'
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
        self.testDirectory = r'C:\Users\zerodeku\Documents\DHM data\16.06.16 - E. coli M1 - 63x\2016.06.16 15-45'        
        
        # connection to db and import params
        with self.con:
            self.con.row_factory = lite.Row
            cur = self.con.cursor()
            cur.execute("SELECT * FROM ObjectiveTable")
            rows = cur.fetchall()
            self.choiceItems = []
            for row in rows:
                config = row["productionNo"]
                name = str(row["MagnificationName"])
                self.choiceItems.append((config, name))
                
            cur.execute("SELECT Value FROM DefaultValueTable" + 
                        " WHERE Name='DefaultDirectory'")
            data = cur.fetchone()
            self.DefaultDirectory = data[0]
            self.DefaultDirectory = self.testDirectory
            
            cur.execute("SELECT Value FROM DefaultValueTable" + 
                        " WHERE Name='TransmissionDHM'")
            data = cur.fetchone()
            self.TransmissionDHM = (int)(data[0])
        
            cur.execute("SELECT Value FROM DefaultValueTable" + 
                        " WHERE Name='CCDPixelSizeUM'")
            data = cur.fetchone()
            self.CCDPixelSizeUM = (float)(data[0])
            
            cur.execute("SELECT Value FROM DefaultValueTable" + 
                        " WHERE Name='IPKoala'")
            data = cur.fetchone()
            self.IPKoala = data[0]
        try:
            open(self.dBKoalapath)
        except:
            self.dbKoalaExist = False
            utilsForProcessing.ErrorMessage('Koala DB is unaccessible')
        if self.dbKoalaExist:
            self.conKoala = lite.connect(self.dBKoalapath)
        
    def choiceItemDB(self):
        return self.choiceItems
        
    def ParamFromChoice(self):
        with self.con:
            self.con.row_factory = lite.Row
            cur = self.con.cursor()
            query = "SELECT * FROM ObjectiveTable WHERE ProductionNo="
            query += str(self.ProductionNo)
            cur.execute(query)
            data = cur.fetchone()
            self.configType = data["ConfigType"]
            self.NA = data["NA"]
            self.ObjectiveNo = data["ObjectiveNo"]
            
        if self.dbKoalaExist and self.conKoala is not None:
            self.conKoala.row_factory = lite.Row
            curKoala = self.conKoala.cursor()
            query = "SELECT * FROM ProductionTable WHERE ProductionNo="
            query += str(self.ProductionNo)
            curKoala.execute(query)
            data = curKoala.fetchone()
            AlgoNo = (int)(data["AlgoNo"])
            query = "SELECT * FROM AlgoTable WHERE AlgoNo="
            query += str(AlgoNo)
            curKoala.execute(query)
            data = curKoala.fetchone()
            RoiX = data["RoiX"]
            RoiY = data["RoiY"]
            RoiWidth = data["RoiWidth"]
            RoiHeight = data["RoiHeight"]
            self.ConfigRoi = [RoiX, RoiY, RoiWidth, RoiHeight]
            
    def DefineDefaultDirectory(self):
        path = str(QFileDialog.getExistingDirectory(None, 
                "Choose Default Direcotry"))
        with self.con:
            self.con.row_factory = lite.Row
            cur = self.con.cursor()
            insert_path = "UPDATE DefaultValueTable SET Value=? WHERE Name=?"
            cur.execute(insert_path, (path, 'DefaultDirectory'))
            self.con.commit()

# debug
if __name__ == "__main__":
    db = DBReader()
    db.ProductionNo = 131
    db.ParamFromChoice()