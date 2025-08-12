Trabont
2025
The following details the process for running this project's scripts as well as the location of data.
As the data has already been run, this file also includes documentation of output locations.

------------------------
Data Location: 
------------------------
Data can be found in and will run from: nfs/masi/trabont/Lymph (henceforth referred to as BaseDir)
- KEY:  
    '*'   = WM, GM, LYMPH, CSF
    '**'  = PSR, kba, T2b, T2R1
    '#'   = subject id folder(s)
    '##'  = files of various names identifiable with subject id number

Within BaseDir, the following files and folders should exist:
'/PreProcessed'
     |-- '/HC'
           |-- '/#' 
     |-- '/MS'
           |-- '/#' 
'/Processed'
     |-- '/HC'
           |-- '/MaskOverlays'
           |-- '/#' 
     |-- '/MS'
           |-- '/MaskOverlays'
           |-- '/#'
     |-- '/FullFit_HC'
           |-- 'FF_HC_all_tissue.mat'
           |-- 'FF_HC_all_SL_*.mat'
           |-- 'FF_HC_all_L_*.mat'
           |-- 'FF_HC_all_G_*.mat'
           |-- 'FF_HC_all_SL_**_MAP.png'
           |-- 'FF_HC_all_L_**_MAP.png'
           |-- 'FF_HC_all_G_**_MAP.png'
           |-- '/#'
     |-- '/FullFit_MS'
           |-- 'FF_MS_all_tissue.mat'
           |-- 'FF_MS_all_SL_*.mat'
           |-- 'FF_MS_all_L_*.mat'
           |-- 'FF_MS_all_G_*.mat'
           |-- 'FF_MS_all_SL_**_MAP.png'
           |-- 'FF_MS_all_L_**_MAP.png'
           |-- 'FF_MS_all_G_**_MAP.png'
           |-- '/#'
     |-- '/SinglePT_HC'
           |-- 'SP_HC_all_tissue.mat'
           |-- 'SP_HC_all_SL_*.mat'
           |-- 'SP_HC_all_L_*.mat'
           |-- 'SP_HC_all_G_*.mat'
           |-- 'SP_HC_all_SL_**_MAP.png'
           |-- 'SP_HC_all_L_**_MAP.png'
           |-- 'SP_HC_all_G_**_MAP.png'
           |-- '/#'
     |-- '/SinglePT_MS'
           |-- 'SP_MS_all_tissue.mat'
           |-- 'SP_MS_all_SL_*.mat'
           |-- 'SP_MS_all_L_*.mat'
           |-- 'SP_MS_all_G_*.mat'
           |-- 'SP_MS_all_SL_**_MAP.png'
           |-- 'SP_MS_all_L_**_MAP.png'
           |-- 'SP_MS_all_G_**_MAP.png'
           |-- '/#'
'/functions/'
     |-- 'yarnykh_pulseMT.m'
     |-- 'Analysis_Yarnykh_Full_Fit.m'
     |-- 'Analysis_Yarnykh_1pt.m'
     |-- 'fit_SSPulseMT_yarnykh_Full_Fit.m'
     |-- 'fit_SSPulseMT_yarnykh_1pt.m'
     |-- 'philipsRFpulse_FA.m'
     |-- 'CWEqMTPulse.m'
     |-- 'absorptionLineShape.m'
     |-- 'FullFit.m'
     |-- 'SinglePT.m'
     |-- 'CombineFullFitMaps.m'
     |-- 'CombineSinglePTMaps.m'
     |-- 'FullFit_Figures.m'
     |-- 'SinglePT_Figures.m'
     |-- 'FF_Excel_ROI_Analysis.m'
     |-- 'SP_Excel_ROI_Analysis.m'
'/MS_List.txt'
'/HC_List.txt'
'/sct_final_registration.sh'
'/ants_final_registration.sh'
'/file_check_all.m'
'/med_dice_centr.m'
'/file_checking.m'
'/Loop_FF.m'
'/Loop_1pt.m'
'/FF_Analysis.m'
'/SP_Analysis.m'

------------------------
Data Information Detailed:
------------------------
'BaseDir/MS_List' and 'BaseDir/HC_List' each contain the list of subject identification numbers found in SMITH/SMITH_HC in XNAT
'BaseDir/functions' contains all vital MATLAB functions required to perform the voxel-based Full Fit and Single Point Analysis
'BaseDir/PreProcessed/MS/#' and 'BaseDir/PreProcessed/HC/#' (where # is a folder for each subject) contain nine files each.
   Those files are:
             (1) *_mFFE_0.65_14slice_e1.nii.gz
             (2) *_B0_Map_In-phase_e2_real.nii.gz |or| *_B0_Map_In_phase_e2_real.nii.gz
             (3) *_WIP_B1_Map_Yarnykh_e2_real.nii.gz
             (4) *_T2W_DRIVE_CLEAR.nii.gz
             (5) *_PulseMT_16_dyn.nii.gz
             (6) *_Clinical_T1W_TSE.nii.gz
             (7) *_Clinical_T2W_TSE.nii.gz
             (8) *_MFA_T120170823.nii.gz
             (9) *_Clinical_T2W_TSE.01_ubMask.nii.gz
'BaseDir/Processed/MS/#' and 'BaseDir/Processed/HC/#' (where # is a folder for each subject), multiple files exist.
   The most pertinent files after registration of SCT and ANTs are:
             (1) #_registered_MT_1.nii.gz
             (2) #_registered_MFA_1.nii.gz
             (3) #_registered_B0_1.nii.gz
             (4) #_registered_B1_1.nii.gz
             (5) #_registered_T2AX_1.nii.gz
             (6) #_registered_T1_1.nii.gz
             (7) #_registered_GM.nii.gz
             (8) #_registered_WM.nii.gz
             (9) #_registered_CSF.nii.gz
            (10) #_registered_LYMPH_ANTS.nii.gz
            (11) #_registered_MT_ANTS.nii.gz
            (12) #_registered_MFA_ANTS.nii.gz
            (13) #_registered_B0_ANTS.nii.gz
            (14) #_registered_B1_ANTS.nii.gz
            (15) #_registered_T2AX_ANTS.nii.gz
            (16) #_registered_T1_ANTS.nii.gz


------------------------
Process to Run Project from BaseDir:
------------------------
  (1) 'sct_final_registration.sh'
  (2) 'ants_final_registration.sh'
  (3) 'file_check_all.m
  (5) 'Loop_FF.m'
  (6) 'FF_Analysis.sh'
  (7) 'Loop_1pt.sh'
  (8) 'SP_Analysis.sh'


------------------------
Software Required
------------------------
SCT 7.0 + FSL
AntsRegistration
MATLAB R2024b or newer
