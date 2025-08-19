## 📝 Information Details
**DESCRIPTION OF FUNCTIONS AND FILES**

##### [`/functions`](./functions) 
Contains all vital MATLAB functions required to perform the voxel-wise **Full Fit** and **Single Point** analysis

---

##### [`MS_List.txt`](./MS_List.txt), [`HC_List.txt`](./HC_List.txt)
Each file contains a list of Subject Identification Numbers for MS / HC groups found in **SMITH** / **SMITH_HC** in XNAT

---

##### [`file_check_all.m`](./file_check_all.m)
Checks registered file details and generates a subject exclusion lists (files missing).\
**All outputs should be reviewed *before* Loop\_FF.m.**

Calls:
  * [`maskOverlay.m`](./maskOverlay.m)
    * Produces ROI overlays on MT and T2 images saved in `Processed/*/MaskOverlays/`
  * [`file_checking.m`](./file_checking.m)
    * Produces registered subject file dimensions ('NA' if file does not exist) Excel
  * [`med_dice_centr.m`](./med_dice_centr.m)
    * Produces a registered subject spinal cord segmentation DICE and Centroid Difference Excel file
    * Produces a **median ROI Full Fit Analysis Excel File**
    * Produces median z-spectrum fit figures (one per lineshape per tissue/slice)

 ⚠️ **Notes:**
* One ROI overlay is generated per subject
* DICE + Centroid Difference + Median Z-spectrum uses only subjects **not excluded** in file_checking.m
* The median ROI Full Fit Analysis Excel File is not a voxel-wise analysis

**Example Outputs**

<details>
<summary>Show Examples</summary>

* Mask Overlay (`*_maskOverlay.jpg`)
* File Checks (`FILE_CHECK.xlsx`)
* Dice + Centroid + Median Fits (`FF_MED_DICE_CENTROID.xlsx`)
* Z-Spectrum Fit Plots (`*_FitZ.png`)

</details>

---

##### [`Loop_FF.m`](./Loop_FF.m) / [`Loop_1pt.m`](./Loop_1pt.m)
Loops through all groups and subjects → performs voxelwise **Full Fit** or **Single Point**

Calls:
  * [`fullFit.m`](./fullFit.m), [`singlePTFit.m`](./singlePTFit.m)
    * Produce per-subject parameter maps (`.mat`) for designated slices for each lineshape fit
    * Produce a figure of the PSR overlay on the MT slice 
  * [`combine.m`](./combine.m)
    * Produces combined parameter maps (`.mat`) using all per-subject parameter maps for each group
    * Produces a figure subplot (each plot designated by ID and slice) for each lineshape/parameter/group

**Example Outputs**

<details>
<summary>Show Examples</summary>

* Dynamics QC images (`*_MTdynamics.png`)
* Subject fit parameter matrices (`*_slice*_*.mat`)
* Parameter overlays (`*_PSR_overlay.png`)
* Group combined maps (`All_*_MAP.png`)
* Group combined matrices (`All_*_*.mat`)

</details>

---

##### [`FF_Analysis.m`](./FF_Analysis.m) / [`SP_Analysis.m`](./SP_Analysis.m)
Generates **figures and Excel outputs** for group comparisons.

Calls:
  * [`histBoxWhisker.m`](./histBoxWhisker.m)
    * Produces histograms and box/whisker plots (by parameter, tissue, lineshape, group).
    * Produces a `.mat` with median HC + MS parameters per tissue/lineshape.
    * Produces an Excel file with histogram statistics
  * [`fit_metrics.m`](./fit_metrics.m)
    * Produces an Excel file with fit metrics (chi², chi²p, resn) for each group/lineshape.

**Example Outputs**

<details>
<summary>Show Examples</summary>

* Box + Whisker plots (`Whisker_*.png`, `Box_*.png`)
* Histograms + KDE (`*_Hist.png`)
* Histogram details (`bin_widths_data_span.xlsx`)
* Fit metrics (`fit_metric_bounds.xlsx`)

</details>

---
