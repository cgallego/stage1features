# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 16:48:07 2014

@ author (C) Cristina Gallego, University of Toronto
"""
import sys, os
import string
import datetime
from numpy import *
import datetime
from base import Base, engine
import database

from sqlalchemy.orm import sessionmaker

#!/usr/bin/env python
class QuerymyDatabase(object):
    """
    USAGE:
    =============
    localdata = QuerymyDatabase()
    """        
        
    def __init__(self): 
        """ initialize QueryDatabase """
        
                              
    def queryby_lesionid(self, lesion_id):
        """
        run : Query by Lesion_id on local database
        
        Inputs
        ======
        lesion_id : (int)    My CADStudy lesion_id        
        StudyID : (int)    CAD StudyID
        redateID : (int)  CAD StudyID Data of exam (format yyyy-mm-dd)
        
        Output
        ======
        """               
        # Create the database: the Session. 
        self.Session = sessionmaker()
        self.Session.configure(bind=engine)  # once engine is available
        session = self.Session() #instantiate a Session
        
        # for cad_case in session.query(Cad_record).order_by(Cad_record.pt_id): 
        #     print cad_case.pt_id, cad_case.cad_pt_no_txt, cad_case.latest_mutation_status_int    
        for lesion in session.query(database.Lesion_record).filter(database.Lesion_record.lesion_id == str(lesion_id)):
            # print results
            if not lesion:
                print "lesion is empty"

        return lesion 
                   

           
           

