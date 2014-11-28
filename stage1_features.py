# -*- coding: utf-8 -*-
"""
Master python script to features engineer stage1 features

Arguments:
============
sys.argv[1] = input text file with one case per line in the following format:

lesion_id	BenignNMaligNAnt	DynSeries_id	T2Series_id	cad_pt_no_txt	lesionfile	exam_img_dicom_txt	exam_a_number_txt	exam_dt_datetime	exam_find_mri_mass_yn	exam_find_mri_nonmass_yn	exam_find_side_int	proc_pt_procedure_id	proc_proc_dt_datetime	proc_proc_side_int	lesion_label	lesion_diagnosis


Created on Thu Oct 30 11:24:12 2014
@ author (C) Cristina Gallego, University of Toronto
----------------------------------------------------------------------
 """
import os, os.path
import sys
from sys import argv

from retrieve_mydatabase import *
from send2_mydatabase import *
from newFeatures import *

if __name__ == '__main__':    
    # Get Root folder ( the directory of the script being run)
    path_rootFolder = os.path.dirname(os.path.abspath(__file__)) 
    print path_rootFolder
       
    # Open filename list
    print sys.argv[1]
    file_ids = open(sys.argv[1],"r")
    file_ids.seek(0)
    line = file_ids.readline()
    
    while ( line ) : 
        # Get the line: 
        fileline = line.split()
        lesion_id = int(fileline[0] )
        BenignNMaligNAnt = fileline[1]
        SeriesID = fileline[2]
        T2SeriesID = fileline[3]
        StudyID = fileline[4]
        Lesionfile = fileline[5]
        DicomExamNumber = fileline[6]
        AccessionNumber = fileline[7]
        dateID = fileline[8]
        ismass = bool(fileline[9])
        isnonmass = bool(fileline[10])
        cond = fileline[15] 
        Diagnosis = fileline[16]     
        sideBreast = fileline[14]
        
        if(lesion_id>272):
            DicomExamNumber = AccessionNumber
        
        #############################
        ###### 1) Querying Research database for clinical, pathology, radiology data
        #############################
        RetrieveDB = Retrieve()
        RetrieveDB.querymyDatabase(lesion_id)
        
        #############################
        ###### 2) Get T2 info
        #############################
        img_folder ='Z:/Cristina/MassNonmass'+os.sep+cond[:-1]+'/'
        Send2DB = Send()
        [path_T2Series, T2SeriesID] = Send2DB.processT2(T2SeriesID, img_folder, StudyID, DicomExamNumber)
    
        #############################                  
        ###### 3) Load segmentation and display
        #############################
        #[series_path, phases_series, lesion3D] = Send2DB.loadSegment(path_rootFolder, cond, StudyID, DicomExamNumber, SeriesID, Lesionfile, T2SeriesID, path_T2Series, lesion_id, sideBreast)
        [series_path, phases_series, lesion3D, chgSeg] = Send2DB.checkSegment(path_rootFolder, cond, StudyID, DicomExamNumber, SeriesID, Lesionfile, T2SeriesID, path_T2Series, lesion_id, sideBreast)

        #############################
        # 4) Extract Lesion and Muscle Major pectoralies signal                                   
        ############################# 
        if T2SeriesID != 'NONE':
            try:
                T2info = RetrieveDB.lesion.f_T2[0] 
                line_muscleVOI = T2info.bounds_muscleSI
            except:
                T2info=[]
                line_muscleVOI='[] []'
            
            [T2_muscleSI, muscle_scalar_range, bounds_muscleSI, T2_lesionSI, lesion_scalar_range, LMSIR, morphoT2features, textureT2features] = Send2DB.T2_extract(T2SeriesID, path_T2Series, line_muscleVOI, lesion3D)
              
        #############################
        ###### 4) Extract new features from each DCE-T1 and from T2 using segmented lesion
        #############################
        newfeatures = newFeatures(Send2DB.load, Send2DB.loadDisplay)
        [deltaS, t_delta, centerijk] = newfeatures.extract_MRIsamp(series_path, phases_series, lesion3D, T2SeriesID)
        
        # generate nodes from segmantation 
        [nnodes, curveT, earlySE, dce2SE, dce3SE, lateSE, ave_T2, prop] = newfeatures.generateNodesfromKmeans(deltaS['i0'], deltaS['j0'], deltaS['k0'], deltaS, centerijk, T2SeriesID)    
        [kmeansnodes, d_euclideanNodes] = prop
        
        # pass nodes to lesion graph
        G = newfeatures.createGraph(nnodes, curveT, prop)                   
    
        [degreeC, closenessC, betweennessC, no_triangles, no_con_comp] = newfeatures.analyzeGraph(G)        
        network_measures = [degreeC, closenessC, betweennessC, no_triangles, no_con_comp]
        
        #Send2DB.loadDisplay.iren1.Start()
        
        if chgSeg:
            #############################
            ###### 4) Extract Dynamic features
            #############################
            [dyn_inside, dyn_contour] = Send2DB.extract_dyn(series_path, phases_series, Send2DB.lesion3D)
            
            #############################
            ###### 5) Extract Morphology features
            #############################
            morphofeatures = Send2DB.extract_morph(series_path, phases_series, Send2DB.lesion3D)
            
            #############################        
            ###### 6) Extract Texture features
            #############################
            texturefeatures = Send2DB.extract_text(series_path, phases_series, Send2DB.lesion3D)       
        
            #############################
            ###### 5) End and Send record to DB
            #############################
            Send2DB.addRecordDB_features(5555, dyn_inside, dyn_contour, morphofeatures, texturefeatures)   
        
        if T2SeriesID != 'NONE':
            #############################
            # 4) Extract Lesion and Muscle Major pectoralies signal                                   
            ############################# 
            [T2_muscleSI, muscle_scalar_range, bounds_muscleSI, T2_lesionSI, lesion_scalar_range, LMSIR, morphoT2features, textureT2features] = Send2DB.T2_extract(T2SeriesID, path_T2Series, line_muscleVOI, lesion3D)       
            dataCase={}
            dataCase['finding.t2_signal_int']="NA"
            Send2DB.addRecordDB_T2(lesion_id, T2SeriesID, dataCase, morphoT2features, textureT2features, T2_muscleSI, muscle_scalar_range, bounds_muscleSI, T2_lesionSI, lesion_scalar_range, LMSIR)

        newfeatures.addRecordDB_stage1(lesion_id, d_euclideanNodes, earlySE, dce2SE, dce3SE, lateSE, ave_T2, network_measures)
                    
        ## continue to next case
        line = file_ids.readline()
        print line

        
    file_ids.close()