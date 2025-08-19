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


From _**BaseDir**_  

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

## 📑 Description of Functions and Files

| File/Folder                         | Notes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| ----------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **file\_check\_all.m**              | **Purpose:** Checks registered file details and generates subject exclusion list when files are missing. <br><br> **Calls:** <br> • `maskOverlay.m` → ROI overlay images (`/MaskOverlays/#_maskOverlay.jpg`) <br> • `file_checking.m` → registered subject file dimensions (`FILE_CHECK.xlsx`) <br> • `med_dice_centr.m` → DICE + Centroid Excel, median ROI check, Z-spectrum figs <br><br> **Notes:** <br> • Only produces one overlay per subject <br> • `med_dice_centr.m` skips excluded subjects <br> • Median ROI fit ≠ voxelwise analysis (validation only) <br><br> [🔗 See Example](#example-file-check) |
| **Loop\_FF.m / Loop\_1pt.m**        | **Purpose:** Iterates through subjects in groups (HC/MS) to perform either full fit (`Loop_FF.m`) or single-point fit (`Loop_1pt.m`). <br><br> **Calls:** <br> • `fullFit.m` → voxelwise multi-parametric fitting <br> • `singlePTFit.m` → single-point parameter estimation <br><br> **Outputs:** <br> • Subject parameter maps (.mat) <br> • PSR overlay images <br> • Combined parameter maps/matrices (group-level) <br><br> **Notes:** <br> • Loops require PreProcessed subjects in both MS and HC <br> • Creates large outputs; watch disk usage <br><br> [🔗 See Example](#example-loop-ff)                |
| **FF\_Analysis.m / SP\_Analysis.m** | **Purpose:** Creates group-level analysis products. <br><br> **Calls:** <br> • `histBoxWhisker.m` → histogram + KDE + boxplots <br> • `fit_metrics.m` → chi², chi²p, resn analysis <br><br> **Outputs:** <br> • Boxplots (outline + filled + all lineshapes) <br> • Histograms + KDE overlays <br> • Excel files (bin widths, fit metrics) <br><br> **Notes:** <br> • Exclusion masks applied before figure generation <br> • Excel files are group-specific (HC/MS) <br><br> [🔗 See Example](#example-ff-analysis)                                                                                               |

---

## 🎯 Example of Outputs

<a name="example-file-check"></a>

<details>
  <summary> 🧩 file_check_all.m (click to expand)</summary>

---

### (1) Mask Overlay (`maskOverlay.m`)

**Path to Output:**

```ruby
BaseDir/Processed/MS/MaskOverlays/137929_maskOverlay.jpg
BaseDir/Processed/HC/MaskOverlays/14329_maskOverlay.jpg
```

Overlay images used to verify ROI placement.

![137929\_maskOverlay](https://github.com/user-attachments/assets/b6a54f3e-002b-4a45-961e-733c093aa862)

**Legend:** Yellow = CSF • Blue = WM • Red = GM • Green = LYMPH

---

### (2) File Check (`file_checking.m`)

**Path to Output:**

```ruby
BaseDir/Processed/FullFit_MS/FILE_CHECK.xlsx
BaseDir/Processed/FullFit_HC/FILE_CHECK.xlsx
```

Checks existence + dimensions.

**Excel Sheets:**

* **DIMS\_SCT:**

  <img width="186" height="165" alt="image" src="https://github.com/user-attachments/assets/1f68455f-225d-4d6c-b11c-c526bb472455" />  
* **DIMS\_ANTS:**

  <img width="186" height="122" alt="image" src="https://github.com/user-attachments/assets/aaee806e-f5e7-49ff-b85a-1a02baa88bb2" />

---

### (3) DICE + Centroid + Median Full Fit (`med_dice_centr.m`)

**Path to Output:**

```ruby
BaseDir/Processed/FullFit_MS/FF_MED_DICE_CENTROID.xlsx
BaseDir/Processed/FullFit_HC/FF_MED_DICE_CENTROID.xlsx
```

* **SL sheet:**

  <img width="625" height="141" alt="image" src="https://github.com/user-attachments/assets/7c51d901-26a4-4f94-9820-eae52d910dbc" />  
* **DICE\_CENTROID sheet:**

  <img width="281" height="169" alt="image" src="https://github.com/user-attachments/assets/6b52afd6-f144-474f-ac5e-3cacd6e0e6ed" />

---

### (4) Fit Figures (Median ROI Fits)

**Path to Output:**

```ruby
BaseDir/Processed/FullFit_MS/137929/137929_slice*_*_*_FitZ.png
BaseDir/Processed/FullFit_HC/14329/14329_slice*_*_*_FitZ.png
```

Figures: one per lineshape × tissue × slice.

![137929\_slice08\_WM\_L\_FitZ](https://github.com/user-attachments/assets/94171d49-386b-42af-ab94-ce9e8690a227)

*Example: MS subject 137929, WM, Lorentzian lineshape, slice 8.*

---

</details>

<a name="example-loop-ff"></a>

<details>
  <summary> 🧩 Loop_FF.m (click to expand)</summary>

---

#### (1) Full Fit Outputs (`fullFit.m`)

**(1.1) Dynamic Contrast Image**

```ruby
BaseDir/Processed/FullFit_MS/137929/137929_slice08_MTdynamics.png
BaseDir/Processed/FullFit_HC/14329/14329_slice08_MTdynamics.png
```

![137929\_slice08\_MTdynamics](https://github.com/user-attachments/assets/e241f7c0-bec6-4b1a-aa66-ed8ea639043f)

---

**(1.2) Subject Parameter Matrices (.mat)**

```ruby
BaseDir/Processed/FullFit_MS/137929/137929_slice*_*_*.mat
BaseDir/Processed/FullFit_HC/14329/14329_slice*_*_*.mat
```

![mat file screenshot](https://github.com/user-attachments/assets/46789b32-96d4-48fa-a529-d4e3bd910129)

---

**(1.3) PSR Overlay Image**

```ruby
BaseDir/Processed/FullFit_MS/137929/137929_slice*_PSR_overlay.png
BaseDir/Processed/FullFit_HC/14329/14329_slice*_PSR_overlay.png
```

![137929\_slice08\_PSR\_overlay](https://github.com/user-attachments/assets/256-psr-overlay.png)

---

#### (2) Combined Outputs (`combined.m`)

**(2.1) Combined Tissues by Parameter Images**

```ruby
BaseDir/Processed/FullFit_MS/All_*_*_MAP.png
BaseDir/Processed/FullFit_HC/All_*_*_MAP.png
```

![MS\_all\_SL\_PSR\_MAP](https://github.com/user-attachments/assets/3fb072a9-31fc-4a3b-87b7-edd08f9d348d)

---

**(2.2) Combined Tissues by Parameter Matrices**

```ruby
BaseDir/Processed/FullFit_MS/All_*_*.mat
BaseDir/Processed/FullFit_HC/All_*_*.mat
```

![mat merged screenshot](https://github.com/user-attachments/assets/651f2a38-4781-4096-aa0f-cdca1819473e)

---

</details>

<a name="example-ff-analysis"></a>

<details>
  <summary> 🧩 FF_Analysis.m (click to expand)</summary>

---

### (1) Box and Whiskers (`histBoxWhisker.m`)

**Paths:**

```ruby
BaseDir/Processed/FullFit_Analysis/Whisker_OUTLINE_**_***.png
BaseDir/Processed/FullFit_Analysis/Whisker_PTS_**_***.png
BaseDir/Processed/FullFit_Analysis/Box_**.png
```

Examples:

* Outline (per param × lineshape):
  ![Whisker\_OUTLINE\_kba\_SL](https://github.com/user-attachments/assets/98df36f8-9b4f-4edd-bdc9-e12d448f6fbe)
* Filled with points:
  ![Whisker\_PTS\_kba\_SL](https://github.com/user-attachments/assets/58637b3e-f314-4c77-a0b6-9afaa43b25e2)
* All lineshapes combined:
  ![Box\_kba\_UnifiedMask](https://github.com/user-attachments/assets/a1c514c6-7458-4faa-b18b-0897d9184e97)

---

### (2) Histograms + KDE

**Path:**

```ruby
BaseDir/Processed/FullFit_Analysis/***_*_**_Hist.png
```

![SL\_WM\_kba\_Hist](https://github.com/user-attachments/assets/05edaeee-a9b5-4c25-bc5a-42c137565e0a)

---

### (3) Excel

**Paths:**

```ruby
BaseDir/Processed/FullFit_Analysis/bin_widths_data_span.xlsx
BaseDir/Processed/FullFit_Analysis/fit_metric_bounds.xlsx
```

* Histogram analysis:
  ![bin\_widths\_excel](https://github.com/user-attachments/assets/a6773a7e-3431-4389-adbd-8c7ceeca07c3)
* Fit metrics:
  ![chi2\_excel](https://github.com/user-attachments/assets/5dd0a8ea-2636-4053-98ba-6565c2c6e063)

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
