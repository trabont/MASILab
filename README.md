# Lymph Node qMT Project
**Author:** Trabont  
**Year:** 2024 - 2025  

This document describes the process for running this project's scripts, the location of data, and the structure of outputs.  
As the data has already been processed, it also serves as documentation of output locations.

---

## 📍 Data Location

Data can be found in and will run from:

```bash
nfs/masi/trabont/Lymph
````

*(henceforth referred to as **BaseDir**)*

---
## 📂 Directory Structure

<ins> Key <ins>
* `*`   = WM, GM, LYMPH, CSF
* `**`  = PSR, kba, T2b, T2R1, etc.
* `***` = SL, L, G
* `#`   = subject ID folder(s)
* `##`  = files of various names identifiable with subject ID number  


> From _**BaseDir**_  

```plaintext  
├ /PreProcessed
    ├ /HC
      ├ /#
    ├ /MS
      ├ /#

├ /Processed
    ├ /HC
      ├ /MaskOverlays
      ├ /#
    ├ /MS
      ├ /MaskOverlays
      ├ /#

    ├ /FullFit_HC
      ├ All_tissue.mat
      ├ All_SL_*.mat
      ├ All_L_*.mat
      ├ All_G_*.mat
      ├ All_SL_**_MAP.png
      ├ All_L_**_MAP.png
      ├ All_G_**_MAP.png
      ├ /#

    ├ /FullFit_MS
      ├ All_tissue.mat
      ├ All_SL_*.mat
      ├ All_L_*.mat
      ├ All_G_*.mat
      ├ All_SL_**_MAP.png
      ├ All_L_**_MAP.png
      ├ All_G_**_MAP.png
      ├ /#

    ├ SinglePT_HC
      ├ All_tissue.mat
      ├ All_SL_*.mat
      ├ All_L_*.mat
      ├ All_G_*.mat
      ├ All_SL_**_MAP.png
      ├ All_L_**_MAP.png
      ├ All_G_**_MAP.png
      ├ /#

    ├ SinglePT_MS
      ├ All_tissue.mat
      ├ All_SL_*.mat
      ├ All_L_*.mat
      ├ All_G_*.mat
      ├ All_SL_**_MAP.png
      ├ All_L_**_MAP.png
      ├ All_G_**_MAP.png
      ├ /#

    ├ /FullFit_Analysis
      ├ Whisker_OUTLINE_**_***.png
      ├ Whisker_PTS_**_***.png
      ├ Box_**.png
      ├ ***_*_**_Hist.png
      ├ bin_widths_data_span.tsv
      ├ fit_metric_bounds.tsv
      ├ /#

├ /functions
    ├ yarnykh_pulseMT.m
    ├ Analysis_Yarnykh_Full_Fit.m
    ├ Analysis_Yarnykh_1pt.m
    ├ fit_SSPulseMT_yarnykh_Full_Fit.m
    ├ fit_SSPulseMT_yarnykh_1pt.m
    ├ philipsRFpulse_FA.m
    ├ CWEqMTPulse.m
    ├ absorptionLineShape.m

├ MS_List.txt
├ HC_List.txt
├ sct_final_registration.sh
├ ants_final_registration.sh
├ file_check_all.m
├ med_dice_centr.m
├ file_checking.m
├ Loop_FF.m
├ Loop_1pt.m
├ FF_Analysis.m
├ SP_Analysis.m
├ fullFit.m
├ singlePTFit.m
├ histBoxWhisker.m
├ combine.m
├ maskOverlay.m
├ fit_metrics.m
```

---

## 📝 Information Details
**DESCRIPTION OF FUNCTIONS AND FILES**
* [`/functions`](./functions) contains all vital MATLAB functions required to perform the voxel-based **Full Fit** and **Single Point** analysis.
  
* [`MS_List`](./MS_List.txt) and [`HC_List`](./HC_List.txt) each contain the list of subject identification numbers found in **SMITH** / **SMITH\_HC** in XNAT.
  
* [`file_check_all.m`](./file_check_all.m) will check registered file details and generate a subject exclusion list (files missing). All outputs from file_check_all.m **should be reviewed _before_ running Loop_FF.m**.
    - Calls [`maskOverlay.m`](./maskOverlay.m)
      - produces+saves images of ROI masks on MT and T2 images in <ins>/MaskOverlays/#_maskOverlay.jpg'</ins>
          - Note: this will only produce one overlay per subject
    - Calls [`file_checking.m`](./file_checking.m)
      - produces registered subject file dimensions ('NA' if file does not exist) in <ins>'FILE_CHECK.xlsx'</ins>
    - Calls [`med_dice_centr.m`](./med_dice_centr.m)
      - produces a registered subject DICE and Centroid Difference of spinal cord segmentations Excel file
      - produces a median ROI check of the Full Fit Analysis Excel file
        - Note: The median ROI Full Fit Analysis is _not_ a voxelwise analysis
        - Note: The Excel shows one fit per lineshape using a median value of each subject's ROI to validate [`/functions`](./functions) runability
      - produces fitted vs normalized data z-spectrum figures
        - Note: med_dice_centr.m will only use files that are NOT on the exclusion list created by file_checking.m
        - Note: Figures show one fit per lineshape using a median value of each subject's ROI to validate [`/functions`](./functions) runability
          
* [`Loop_FF.m`](./Loop_FF.m) and [`Loop_1pt.m`](./Loop_1pt.m) will loop through all groups and subjects and fit accordingly.
    - Calls [`fullFit.m`](./fullFit.m) and [`singlePTFit.m`](./singlePTFit.m)
      - produces each subject's parameter maps (.mat) for designated slices using fit designated functions (see /functions)
    - Calls [`combine.m`](./combine.m)
      - produces a combined parameter
      - produces a parameter subplot figure of all subject-slices for each parameter for each group
        
* [`FF_Analysis.m`](./FF_Analysis.m) and [`SP_Analysis.m`](./SP_Analysis.m) will create Excel files and figures necessary to analyze fit results. 
    - Calls [`histBoxWhisker.m`](./histBoxWhisker.m)
      - produces parameter, tissue (ROI), lineshape, and group specified Histograms and Box and Whisker Plots figures
      - produces a .mat file of the Median HC and MS parameters of each tissue per lineshape
      - produces an Excel file to view parameter, tissue (ROI), and lineshape specified histogram information
    - Calls [`fit_metrics.m`](./fit_metrics.m)
      - produces an Excel file to view fit metrics across all groups
 
**EXAMPLES OF OUTPUTS**
<details>
  <summary> 🧩 file_check_all.m (click to expand)</summary>

  ---
  
  ### (1) Mask Overlay (./maskOverlay.m)  
  
  Path to Output:
  
```ruby
BaseDir/Processed/MS/MaskOverlays/137929_maskOverlay.jpg
BaseDir/Processed/HC/MaskOverlays/14329_maskOverlay.jpg
```

Overlay images are used to verify ROI placement.
  
  ![137929_maskOverlay](https://github.com/user-attachments/assets/b6a54f3e-002b-4a45-961e-733c093aa862)

  Yellow = CSF\
  Blue = WM\
  Red = GM\
  Green = LYMPH

  ---

  ### (2) File Check (./file_checking.m)  
  
  Path to Output:
  
```ruby
BaseDir/Processed/FullFit_MS/FILE_CHECK.xlsx
BaseDir/Processed/FullFit_HC/FILE_CHECK.xlsx
```

File assists in checking file existence and dimensions.\
If all of a subject's files exist, the dimensions should be the following:
  
  #### (2.1) FILE_CHECK.xlsx : Sheet 'DIMS_SCT'

<img width="186" height="165" alt="image" src="https://github.com/user-attachments/assets/1f68455f-225d-4d6c-b11c-c526bb472455" />

  #### (2.2) FILE_CHECK.xlsx : Sheet 'DIMS_ANTS'

<img width="186" height="122" alt="image" src="https://github.com/user-attachments/assets/aaee806e-f5e7-49ff-b85a-1a02baa88bb2" />

  ---
  ### (3) DICE + Centroid Difference + Median Full Fit (./med_dice_centr.m)

  Path to Output:
  
```ruby
BaseDir/Processed/FullFit_MS/FF_MED_DICE_CENTROID.xlsx
BaseDir/Processed/FullFit_HC/FF_MED_DICE_CENTROID.xlsx
```

  #### (3.1) FF_MED_DICE_CENTROID.xlsx : Sheet 'SL'

<img width="625" height="141" alt="image" src="https://github.com/user-attachments/assets/7c51d901-26a4-4f94-9820-eae52d910dbc" />

  Similar details on Sheets: 'L' and 'G'

  #### (3.2) FF_MED_DICE_CENTROID.xlsx : Sheet 'DICE_CENTROID'

<img width="281" height="169" alt="image" src="https://github.com/user-attachments/assets/6b52afd6-f144-474f-ac5e-3cacd6e0e6ed" />

  ---
  ### (4) Fit Figures (./med_dice_centr.m)
  
  Path to Output:
  
```ruby
BaseDir/Processed/FullFit_MS/137929/137929_slice*_*_*_FitZ.png
BaseDir/Processed/FullFit_HC/14329/14329_slice*_*_*_FitZ.png
```
  These image names are differentiated by slice#_lineshape_ROI.\
  These images will have three figures (one for each lineshape) per tissue and slice.\
  [One Slice = 3 GM, 3 WM, 3 LYMPH, 3 CSF]
  
  These images are used to check the functionality of the fit functions.

  <img width="1313" height="914" alt="137929_slice08_WM_L_FitZ" src="https://github.com/user-attachments/assets/94171d49-386b-42af-ab94-ce9e8690a227" />

  The image above is the MS Median Full Fit of 137929's WM using a Lorentzian lineshape on slice 8.
  
---

</details>



<details>
  <summary> 🧩 Loop_FF.m (click to expand)</summary>
  
  ---

  ##### Outputs 1.1-1.3 are generated for each subject in each group.
  ##### Outputs 2.1-2.2 are generated using all subjects for each group.

  --- 
  
  ### (1) Full Fit Outputs (fullFit.m)

  #### (1.1) Dynamic Contrast Image

  Path to Output:
  
```ruby
BaseDir/Processed/FullFit_MS/137929/137929_slice08_MTdynamics.png
BaseDir/Processed/FullFit_HC/14329/14329_slice08_MTdynamics.png
```

These images check that there is a difference in contrast between Power and Offset Frequency.

<img width="2975" height="798" alt="137929_slice08_MTdynamics" src="https://github.com/user-attachments/assets/e241f7c0-bec6-4b1a-aa66-ed8ea639043f" />

  In this image, we see MS subject 137929's slice 8 having a good contrast difference between frequency offsets and powers. \
  In addition, we can see that the image contains no significant artifacts that depreciate the signal intensity.

 --- 

#### (1.2) Subject Full Fit Parameter Matrix

Path to Output:

```ruby
BaseDir/Processed/FullFit_MS/137929/137929_slice*_*_*.mat
BaseDir/Processed/FullFit_HC/14329/14329_slice*_*_*.mat
```
  Parameter matrix (.mat) files differentiated by slice_lineshape_ROI
  
<img width="280" height="519" alt="image" src="https://github.com/user-attachments/assets/46789b32-96d4-48fa-a529-d4e3bd910129" />

Each .mat will be a 256x256xNum_Parameters array.
 
--- 

#### (1.3) PSR Overlay Image

Path to Output:

```ruby
BaseDir/Processed/FullFit_MS/137929/137929_slice*_PSR_overlay.png
BaseDir/Processed/FullFit_HC/14329/14329_slice*_PSR_overlay.png
```

These overlay images are meant to ensure voxel-wise analysis location and value distribution.

<img width="256" height="197" alt="137929_slice08_PSR_overlay" src="https://github.com/user-attachments/assets/2f4d703c-ffc6-4bec-bc62-413b4753740f" />

Note: The color bar is inaccurate.


---

### (2) Combining all Subjects per Group (combined.m)

#### (2.1) Combined Tissues by Parameter Images

Path to Output:

```ruby
BaseDir/Processed/FullFit_MS/All_*_*_MAP.png
BaseDir/Processed/FullFit_HC/All_*_*_MAP.png
```
Combined Parameter Map images are differentiated by lineshape_parameter.

These images are meant to verify the merge of parameter-specific maps (LYMPH may not appear).

<img width="3300" height="2194" alt="MS_all_SL_PSR_MAP" src="https://github.com/user-attachments/assets/3fb072a9-31fc-4a3b-87b7-edd08f9d348d" />

As you can see in this image, some maps are empty, which requires either a re-run of combine.m or manual deletion.

--- 

#### (2.2) Combined Tissues by Parameter Matrix

Path to Output:

```ruby
BaseDir/Processed/FullFit_MS/All_*_*.mat
BaseDir/Processed/FullFit_HC/All_*_*.mat
```
Combined Parameter .mat files are differentiated by lineshape_parameter.

<img width="280" height="236" alt="image" src="https://github.com/user-attachments/assets/651f2a38-4781-4096-aa0f-cdca1819473e" />


All_*_*.mat files are merged into MS_all_tissues.mat having dimensions: [256x256x3xnSxnP]\
3 = number of lineshapes (SL, L, G)\
nS = number of slices for all subjects that have lymph node segmentations\
nP = number of parameter maps (PSR, kba, R1obs, etc.)

--- 
  
</details>
<details>
  <summary> 🧩 FF_Analysis.m (click to expand)</summary>

---

  ##### All outputs following are generated with HistBoxWhisker.m

  --- 

  ### (1) Box and Whiskers

  Path to Output:

  ```ruby
  BaseDir/Processed/FullFit_Analysis/Whisker_OUTLINE_**_***.png
  BaseDir/Processed/FullFit_Analysis/Whisker_PTS_**_***.png
  BaseDir/Processed/FullFit_Analysis/Box_**.png
  ```
  Images are differentiated by *** = lineshape and ** = parameter.

  #### (1.1) Box and Whisker (Lineshape Exclusivity) - Basic

  The following image is the basic outline of the box and whisker plot comparing Healthy Controls and MS Subjects across all tissues.

  <img width="1979" height="1250" alt="Whisker_OUTLINE_kba_SL" src="https://github.com/user-attachments/assets/98df36f8-9b4f-4edd-bdc9-e12d448f6fbe" />

  In this image, we see all tissue's kba box plots generated from a Super-Lorentzian lineshape.\
  There is one image for every parameter and lineshape. \
  Total Basic Box Plots = N_param x N_lineshapes.

--- 

  #### (1.2) Box and Whisker (Lineshape Exclusivity) - Filled

The following image is the whole box and whisker plot comparing Healthy Controls and MS Subjects across all tissues. In addition, each voxel (color coordinated by group type) is also plotted to show voxel value distribution.
  
<img width="1875" height="1250" alt="Whisker_PTS_kba_SL" src="https://github.com/user-attachments/assets/58637b3e-f314-4c77-a0b6-9afaa43b25e2" />

  In this image, we see all tissue's kba box plots generated from a Super-Lorentzian lineshape.\
  There is one image for every parameter and lineshape. \
  Total Basic Box Plots = N_param x N_lineshapes.

--- 

  #### (1.3) Box and Whisker (All Lineshapes) - Filled

The following image is the whole box and whisker plot comparing lineshape contribution to parameter values across all tissues and groups. In addition, each voxel (color coordinated by group type) is also plotted to show voxel value distribution.

<img width="3900" height="2025" alt="Box_kba_UnifiedMask" src="https://github.com/user-attachments/assets/a1c514c6-7458-4faa-b18b-0897d9184e97" />

 In this image, we see kba box plots generated from all lineshapes.\
  There is one image for every parameter.\
 Note: The bolded value is the median of the combined group data.

--- 

  ### (2) Histogram + KDE

  ```ruby
  BaseDir/Processed/FullFit_Analysis/***_*_**_Hist.png
  ```
  Images are differentiated by *** = lineshape, ** = parameter, * = tissue/ROI.

  These histograms are meant to identify the distribution of parameter values for each group.\

  <img width="3063" height="2374" alt="SL_WM_kba_Hist" src="https://github.com/user-attachments/assets/05edaeee-a9b5-4c25-bc5a-42c137565e0a" />

In this image, we see the White Matter (WM) kba distribution fitted by the Super-Lorentzian lineshape (SL).\
  There is one image for every parameter and lineshape.

--- 

  ### (3) Excel

  ```ruby
  BaseDir/Processed/FullFit_Analysis/bin_widths_data_span.xlsx
  BaseDir/Processed/FullFit_Analysis/fit_metric_bounds.xlsx
  ```

  #### (3.1) Histogram Analysis/Detail

  The following file details histogram information (values pertain to combined groups).

<img width="774" height="351" alt="image" src="https://github.com/user-attachments/assets/a6773a7e-3431-4389-adbd-8c7ceeca07c3" />

Note: This image has been cropped to only show kba.\
Note: The appearance of this sheet may be different.

--- 

  #### (3.2) Fitting Analysis (Chi, Chi2p, and Resn)

  The following file details each group's fitting quality for each lineshape.

<img width="891" height="544" alt="image" src="https://github.com/user-attachments/assets/5dd0a8ea-2636-4053-98ba-6565c2c6e063" />

Note: This image has been cropped to only show chi2.\
Note: The appearance of this sheet may be different.

--- 

</details>



---

## 📑 Pre-Processed Data Details

Path Directory to Pre-Processed Data:
```bash
BaseDir/PreProcessed/MS/#
BaseDir/PreProcessed/HC/#
```
> (`#` is a folder for each subject)


Each subject folder contains **nine files**:

1. `*_mFFE_0.65_14slice_e1.nii.gz`
2. `*_B0_Map_In-phase_e2_real.nii.gz` **or** `*_B0_Map_In_phase_e2_real.nii.gz`
3. `*_WIP_B1_Map_Yarnykh_e2_real.nii.gz`
4. `*_T2W_DRIVE_CLEAR.nii.gz`
5. `*_PulseMT_16_dyn.nii.gz`
6. `*_Clinical_T1W_TSE.nii.gz`
7. `*_Clinical_T2W_TSE.nii.gz`
8. `*_MFA_T120170823.nii.gz`
9. `*_Clinical_T2W_TSE.01_ubMask.nii.gz`

<details>
  <summary> 🧩 Example of MS Subject 141108 (click to expand)</summary>

Subject 141108 Pre-processed Path and Contained Files:
  
```bash
 BaseDir/PreProcessed/MS/141108
 
    141108_mFFE_0.65_14slice_e1.nii.gz
    141108_B0_Map_In-phase_e2_real.nii.gz
    141108_WIP_B1_Map_Yarnykh_e2_real.nii.gz
    141108_T2W_DRIVE_CLEAR.nii.gz
    141108_PulseMT_16_dyn.nii.gz
    141108_Clinical_T1W_TSE.nii.gz
    141108_Clinical_T2W_TSE.nii.gz
    141108_MFA_T120170823.nii.gz
    141108_Clinical_T2W_TSE.01_ubMask.nii.gz

```
</details>

---

## ⏳ Processed Data Details

Located in:

```bash
BaseDir/Processed/MS/#
BaseDir/Processed/HC/#
```

(where `#` is a folder for each subject)

**💾 Most pertinent files _after_ SCT and ANTs registration:**

Required to run DICE, Centroid, Full Fit, and Single Point Analysis

1. `#_registered_MT_1.nii.gz`
2. `#_registered_MFA_1.nii.gz`
3. `#_registered_B0_1.nii.gz`
4. `#_registered_B1_1.nii.gz`
5. `#_registered_T2AX_1.nii.gz`
6. `#_registered_T1_1.nii.gz`
7. `#_registered_GM.nii.gz`
8. `#_registered_WM.nii.gz`
9. `#_registered_CSF.nii.gz`
10. `#_registered_LYMPH_ANTS.nii.gz`
11. `#_registered_MT_ANTS.nii.gz`
12. `#_registered_MFA_ANTS.nii.gz`
13. `#_registered_B0_ANTS.nii.gz`
14. `#_registered_B1_ANTS.nii.gz`
15. `#_registered_T2AX_ANTS.nii.gz`
16. `#_registered_T1_ANTS.nii.gz`  
17. `#_registered_mFFE_1_sc.nii.gz`
18. `#_registered_T1_1_sc.nii.gz`
19. `#_registered_T2AX_1_sc.nii.gz`
20. `#_MFA_dyn_0000_in_MT_sc.nii.gz`
21. `#_MT_dyn_0000_reg_sc.nii.gz`
22. `#_MT_dyn_0008_reg_sc.nii.gz`

---

## ▶️ Process to Run Project from `BaseDir`

Run the following in order:

1. `./sct_final_registration.sh`
2. `./ants_final_registration.sh`
3. `file_check_all.m`
4. `Loop_FF.m`
5. `FF_Analysis.sh`
6. `Loop_1pt.sh`
7. `SP_Analysis.sh`

---

## 💻 Software Required

* **[SCT](https://spinalcordtoolbox.com/)** 7.0 + **FSLEyes**
* **[ANTs](http://stnava.github.io/ANTs/)** Registration
* **MATLAB** R2024b or newer
