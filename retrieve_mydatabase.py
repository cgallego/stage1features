# -*- coding: utf-8 -*-
"""
Process all pipeline for a record and
Send a record to database

Created on Fri Jul 25 15:46:19 2014

@ author (C) Cristina Gallego, University of Toronto
----------------------------------------------------------------------
 """
import os, os.path
import sys
from sys import argv, stderr, exit
import numpy as np
import dicom
import psycopg2
import pandas as pd

from query_mydatabase import *
import processDicoms
import datetime

from inputs_init import *
from display import *
from features_dynamic import *
from features_morphology import *
from features_texture import *
from features_T2 import *
from segment import *
import pylab      
import annot
 
from sqlalchemy import Column, Integer, String
from sqlalchemy.orm import sessionmaker


class Retrieve(object):
    """
    USAGE:
    =============
    RetrieveDB = Retrieve()
    """
    def __init__(self): 
        self.dataInfo = []
        self.queryDatalocal = QuerymyDatabase() 
        self.lesion = []
        self.casesFrame={}         
        
        
    def querymyDatabase(self, lesion_id):
        """ Querying with condition (e.g: mass, non-mass)"""
        #############################
        ###### 1) Querying Research database for clinical, pathology, radiology data
        #############################
        print "Executing SQL local..."
        # perform query
        self.lesion = self.queryDatalocal.queryby_lesionid(lesion_id)           
                
        ## append collection of cases
        self.casesFrame = pd.DataFrame(columns=self.lesion.__dict__.keys())
        self.casesFrame = self.casesFrame.append(self.lesion.__dict__, ignore_index=True) # 20    
        print self.casesFrame
        
        # ask for info about lesion row data from query
        if(self.casesFrame['exam_find_mri_mass_yn']):
            massrecord = self.lesion.mass_lesion[0]
            ## append collection of cases 
            self.massframe = pd.DataFrame(columns=massrecord.__dict__.keys()[1:])
            self.massframe = self.massframe.append(massrecord.__dict__, ignore_index=True)
            self.lesionFrame = self.massframe
            
        if(self.casesFrame['exam_find_mri_nonmass_yn']):
            nonmassrecord = self.lesion.nonmass_lesion[0]
            ## append collection of cases 
            self.nonmassframe = pd.DataFrame(columns=nonmassrecord.__dict__.keys()[1:])
            self.nonmassframe = self.nonmassframe.append(nonmassrecord.__dict__, ignore_index=True)
            self.lesionFrame = self.nonmassframe
        
        print self.lesionFrame
        return 
        

   