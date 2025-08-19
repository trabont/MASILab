## 📝 Information Details

#### [`/functions`](./functions) 
Core MATLAB functions for **Full Fit** and **Single Point** analysis

---

#### [`MS_List.txt`](./MS_List.txt), [`HC_List.txt`](./HC_List.txt)
Subject IDs for MS / HC groups from **SMITH** in XNAT.

---

#### [`file_check_all.m`](./file_check_all.m)
Checks registered files and generates exclusion lists.\
**Run this *before* Loop\_FF.m.**

Calls:
  * [`maskOverlay.m`](./maskOverlay.m)
    * Produces ROI overlays saved in `Processed/*/MaskOverlays/` (one per subject).
    * Used to verify ROI placement.
  * [`file_checking.m`](./file_checking.m)
    * Produces `FILE_CHECK.xlsx` (dimensions; “NA” = missing file).
    * If all subject files exist, dimensions should match expected SCT/ANTs values.
  * [`med_dice_centr.m`](./med_dice_centr.m)
    * Produces `FF_MED_DICE_CENTROID.xlsx` with Dice + centroid difference.
    * Generates a **median ROI Full Fit check** (not voxelwise).
    * Uses only subjects **not excluded** in `FILE_CHECK.xlsx`.
    * Produces median z-spectrum fit figures (one per lineshape per tissue/slice).
    * These are meant to validate [`/functions`](./functions) runability, not provide final analysis.

**Example Outputs**

<details>
<summary>Show Examples</summary>

* Mask Overlay (`*_maskOverlay.jpg`)
* File Checks (`FILE_CHECK.xlsx`)
* Dice + Centroid + Median Fits (`FF_MED_DICE_CENTROID.xlsx`)
* Z-Spectrum Fit Plots (`*_FitZ.png`)

</details>

---

#### [`Loop_FF.m`](./Loop_FF.m) / [`Loop_1pt.m`](./Loop_1pt.m)
Loops through all groups and subjects → performs voxelwise **Full Fit** or **Single Point** fitting.

* Calls:

  * [`fullFit.m`](./fullFit.m), [`singlePTFit.m`](./singlePTFit.m)

    * Produce per-subject parameter maps (`.mat`) for each slice/ROI.
  * [`combine.m`](./combine.m)

    * Produces combined parameter maps across all subjects.
    * Also generates parameter subplots per lineshape/parameter/group.

⚠️ **Notes:**

* Figures are **QC only**, meant to verify fit function execution.
* Combined group `.mat` files follow shape `[256 x 256 x 3 x nSlices x nParams]`:

  * `3` = number of lineshapes (SL, L, G)
  * `nSlices` = total slices with lymph node segmentations
  * `nParams` = number of parameter maps (PSR, kba, R1obs, etc.)

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

#### [`FF_Analysis.m`](./FF_Analysis.m) / [`SP_Analysis.m`](./SP_Analysis.m)
Generates **figures and Excel outputs** for group comparisons.

* Calls:

  * [`histBoxWhisker.m`](./histBoxWhisker.m)

    * Produces histograms and box/whisker plots (by parameter, tissue, lineshape, group).
    * Exports `.mat` with median HC + MS parameters per tissue/lineshape.
    * Exports Excel with histogram statistics.
  * [`fit_metrics.m`](./fit_metrics.m)

    * Produces Excel files with fit metrics (chi², chi²p, resn) for each group/lineshape.

⚠️ **Notes:**

* Median ROI checks are **not voxelwise** — they are validation summaries.
* Histograms and KDE plots show parameter distributions per ROI/group.
* Box/whisker plots can be generated per lineshape or across all lineshapes.
* Excel outputs provide **bin widths, histogram details, and fit quality ranges**.

**Example Outputs**

<details>
<summary>Show Examples</summary>

* Box + Whisker plots (`Whisker_*.png`, `Box_*.png`)
* Histograms + KDE (`*_Hist.png`)
* Histogram details (`bin_widths_data_span.xlsx`)
* Fit metrics (`fit_metric_bounds.xlsx`)

</details>

---
