#!/usr/bin/env bash
set -euo pipefail # Exit on error, undefined variable, or failed command in a pipeline
get_latest_file() {
  find "$1" -type f -iname "$2" \
    | sort -t/ -k5,5n -r \
    | head -n1 || true
}


# =============================================================================
# Batch registration pipeline for qMT study (Aligned to PAM50 via mFFE)
#   - Loops through each subject folder in the current directory
#   - Finds MT, MFA, B0, B1, T1, T2 (axial), T2 (sagittal), and mFFE files
#   - Orients everything to RPI
#   - T2SAG → mFFE for vertebral labels
#   - mFFE → PAM50 is the master warp
#   - All single-volumes and dynamics registered to mFFE → then to PAM50
#   - Chain‐registers 4D MT and 4D MFA within native space, then brings them into T2 sagittal space, then PAM50
#   - Cleans up temporary warp/inv files after each step
#
# Usage:
#   cd /path/to/parent_of_subject_folders
# =============================================================================

BASE_DIR="$(pwd)"                         # Directory where all subject folders live
PROCESSED_DIR="$BASE_DIR/Processed"       # Parent for processed outputs
QC_DIR="$PROCESSED_DIR/qc_singleSubj"     # QC outputs under Processed
mkdir -p "$QC_DIR"                        # Create the QC folder if it doesn’t exist

# Initialize counter for numbering subjects
counter=1

# Loop over each subject directory (e.g., 137773, 142019, etc.)
for subj_path in "$BASE_DIR"/*/; do
  # Skip if not a directory
  [[ "$(basename "$subj_path")" == "Processed" ]] && continue
  [ -d "$subj_path" ] || continue

  subj="$(basename "$subj_path")"
  id="$counter"

  echo "================================================================"
  echo "Processing subject: $subj  (ID = $id)"
  echo "----------------------------------------------------------------"

  # Move into the subject folder
  cd "$subj_path"

  # ---------------------------------------------------------------------------
  # 0) Locate input files for this subject (may vary in filenames, so we use `find`)
  # ---------------------------------------------------------------------------
  mfa_file=$(get_latest_file "$PWD" "*_MFA_*.nii*")
  mt_file=$(get_latest_file "$PWD" "*_PulseMT_16_dyn.nii*")
  b0_file=$(get_latest_file "$PWD" "*_B0_Map_In-phase_e2_real.nii*")
  b1_file=$(get_latest_file "$PWD" "*_WIP_B1_Map_Yarnykh_CLEAR_e2.nii*")
  t1_file=$(get_latest_file "$PWD" "*_Clinical_T1W_TSE.nii*")
  t2ax_file=$(get_latest_file "$PWD" "*_Clinical_T2W_TSE.nii*")
  sagittal_file=$(get_latest_file "$PWD" "*_DRIVE_CLEAR.nii*")
  mffe_file=$(get_latest_file "$PWD" "*_mFFE_0.65_14slice_e1.nii*")


  # Verify that required files exist; skip subject if any mandatory file is missing
  if [[ -z "$mfa_file" || -z "$mt_file" || -z "$b0_file" || -z "$b1_file" || \
        -z "$t1_file" || -z "$t2ax_file" || -z "$sagittal_file" || -z "$mffe_file" ]]; then
    echo "ERROR: One or more required files are missing for subject $subj. Skipping..."
    cd "$BASE_DIR"
    continue
  fi

  echo "  Found MFA file:      $mfa_file"
  echo "  Found MT file:       $mt_file"
  echo "  Found B0 file:       $b0_file"
  echo "  Found B1 file:       $b1_file"
  echo "  Found T1 file:       $t1_file"
  echo "  Found T2 axial file: $t2ax_file"
  echo "  Found T2 sag file:   $sagittal_file"
  echo "  Found mFFE file:     $mffe_file"
  echo "  Found all input files. Reorienting to RPI and copying to Processed folder..."
  echo

  mkdir -p "$PROCESSED_DIR/$id"
  cd "$PROCESSED_DIR/$id"

  # ----------------------------------------------------------
  # 1) Reorient all images to RPI (ensures consistent orientation)
  # ----------------------------------------------------------
  declare -A INPUTS=(
    [MFA]="$mfa_file"
    [MT]="$mt_file"
    [B0]="$b0_file"
    [B1]="$b1_file"
    [T1]="$t1_file"
    [T2AX]="$t2ax_file"
    [T2SAG]="$sagittal_file"
    [MFFE]="$mffe_file"
  )

  declare -A ORIENTED

  for key in "${!INPUTS[@]}"; do
    src="${INPUTS[$key]}"
    base="${id}_${key}.nii.gz"
    sct_image -i "$src" -setorient RPI -o "$base"
    ORIENTED[$key]="$base"
  done

  # ----------------------------------------------------------
  # 2) Segment spinal cord on T2 sagittal & label vertebrae
  # ----------------------------------------------------------
  sct_deepseg_sc -i "${ORIENTED[MFFE]}" -c t2 -o "${id}_mFFE_sc.nii.gz" -qc "$QC_DIR"
  sct_create_mask -i "${ORIENTED[MFFE]}" -p centerline,"${id}_mFFE_sc.nii.gz" -size 35mm -o "${id}_mFFE_mask.nii.gz"
  sct_deepseg_sc -i "${ORIENTED[T2SAG]}" -c t2 -o "${id}_T2SAG_sc.nii.gz" -qc "$QC_DIR"
  sct_label_vertebrae -i "${ORIENTED[T2SAG]}" \
                      -s "${id}_T2SAG_sc.nii.gz" \
                      -c t2 -qc "$QC_DIR"

  sct_label_utils -i "${id}_T2SAG_sc_labeled.nii.gz" -vert-body 2,8 -o ${id}_labels_vert.nii.gz

  # ----------------------------------------------------------
  # 3) Register T2sag → PAM50 template (then rename to numbered)
  #    (Also generating WM GM CSF binary masks)
  # ----------------------------------------------------------
  echo "  Registering T2SAG to PAM50..."
  sct_register_to_template -i "${ORIENTED[T2SAG]}" \
                           -s "${id}_T2SAG_sc.nii.gz" \
                           -l "${id}_labels_vert.nii.gz" \
                           -c t2 -qc "$QC_DIR" \

  cp "warp_template2anat.nii.gz" "${id}_warp_template2anat_T2SAG.nii.gz"
  cp "warp_anat2template.nii.gz" "${id}_warp_anat2template_T2SAG.nii.gz"

  sct_register_multimodal -i $SCT_DIR/data/PAM50/template/PAM50_t2.nii.gz \
                          -iseg $SCT_DIR/data/PAM50/template/PAM50_cord.nii.gz \
                          -d "${ORIENTED[MFFE]}" \
                          -dseg "${id}_mFFE_sc.nii.gz" \
                          -m "${id}_mFFE_mask.nii.gz" \
                          -initwarp "${id}_warp_template2anat_T2SAG.nii.gz" \
                          -initwarpinv "${id}_warp_anat2template_T2SAG.nii.gz" \
                          -param step=1,type=seg,algo=centermass:step=2,type=seg,algo=bsplinesyn,slicewise=1,iter=3 \
                          -owarp "${id}_warp_template2anat_native.nii.gz" \
                          -owarpinv "${id}_warp_anat2template_native.nii.gz"

  sct_warp_template -d "${ORIENTED[MFFE]}" -w "${id}_warp_template2anat_native.nii.gz" -a 1 -o ${id}_Labels

  # Register B0, B1, T1, and T2AX to mFFE, then to PAM50
  for key in B0 B1 T1 T2AX; do
    sct_register_multimodal -i "${ORIENTED[$key]}" \
                            -d "${ORIENTED[MFFE]}" \
                            -dseg "${id}_mFFE_sc.nii.gz" \
                            -param step=1,type=im,algo=rigid,metric=MI:step=2,type=im,algo=slicereg,metric=MI \
                            -x spline \
                            -o "${id}_${key}_to_mFFE.nii.gz"

    sct_apply_transfo -i "${id}_${key}_to_mFFE.nii.gz" \
                      -d "$SCT_DIR/data/PAM50/template/PAM50_t2.nii.gz" \
                      -w "${id}_warp_anat2template_native.nii.gz" \
                      -x spline \
                      -o "${id}_${key}_to_PAM.nii.gz"
  done

  # ----------------------------------------------------------
  # 4) Chain‐register 4D MFA → co-reg native → mFFE & PAM50
  # ----------------------------------------------------------
  echo "  Splitting and registering MFA to mFFE..."
  n_dyn=$(fslhd 1_MFA.nii.gz | grep "^dim4" | awk '{print $2}')
  fslsplit "${ORIENTED[MFA]}" "${id}_MFA_dyn_" -t
  mfa_ref="${id}_MFA_dyn_0000.nii.gz"
  sct_deepseg_sc -i "$mfa_ref" -c t2 -o "${mfa_ref%.*}_sc.nii.gz"

  last_index=$(printf "%04g" $((n_dyn - 1)))
  echo "Number of dynamics to loop: $last_index"

  sct_register_multimodal -i "$mfa_ref" \
                             -d "${ORIENTED[MFFE]}" \
                             -iseg "${mfa_ref%.*}_sc.nii.gz" \
                             -dseg "${id}_mFFE_sc.nii.gz" \
                             -param step=1,type=im,algo=rigid,metric=MI:step=2,type=seg,algo=centermass,metric=MeanSquares \
                             -x spline \
                             -o "${id}_MFA_dyn_0000_in_mFFE.nii.gz" \
                             -owarp "${id}_warp_MFAref2MFFE.nii.gz" \
                             -owarpinv "${id}_warp_MFFE2MFAref.nii.gz"
  
  sct_apply_transfo -i "${id}_MFA_dyn_0000_in_mFFE.nii.gz" \
                    -d "$SCT_DIR/data/PAM50/template/PAM50_t2.nii.gz" \
                    -w "${id}_warp_anat2template_native.nii.gz" \
                    -x spline \
                    -o "${id}_MFA_dyn_0000_to_PAM.nii.gz"

  for i in $(seq -f "%04g" 1 "$last_index"); do
    in_dyn="${id}_MFA_dyn_${i}.nii.gz"
    out_dyn="${id}_MFA_dyn_${i}_reg.nii.gz"
    sct_register_multimodal -i "$in_dyn" \
                             -d "$mfa_ref" \
                             -dseg "${mfa_ref%.*}_sc.nii.gz" \
                             -param step=1,type=im,algo=rigid,metric=MI \
                             -x spline \
                             -o "$out_dyn"

    sct_apply_transfo -i "$out_dyn" \
                      -d "${ORIENTED[MFFE]}" \
                      -w "${id}_warp_MFAref2MFFE.nii.gz" \
                      -x spline \
                      -o "${id}_MFA_dyn_${i}_in_mFFE.nii.gz"

    sct_apply_transfo -i "${id}_MFA_dyn_${i}_in_mFFE.nii.gz" \
                      -d "$SCT_DIR/data/PAM50/template/PAM50_t2.nii.gz" \
                      -w "${id}_warp_anat2template_native.nii.gz" \
                      -x spline \
                      -o "${id}_MFA_dyn_${i}_to_PAM.nii.gz"
  done

  rm -f warp_${id}_MFA_dyn_*2${id}_MFA_dyn_*.nii.gz \
      ${id}_MFA_dyn_*_reg_inv.nii.gz

  fslmerge -t "${id}_registered_MFA_native.nii.gz" \
  $(for i in $(seq -f "%04g" 0 "$last_index"); do echo "${id}_MFA_dyn_${i}_in_mFFE.nii.gz"; done)
  fslmerge -t "${id}_registered_MFA_PAM.nii.gz" \
  $(for i in $(seq -f "%04g" 0 "$last_index"); do echo "${id}_MFA_dyn_${i}_to_PAM.nii.gz"; done)

  # ----------------------------------------------------------
  # 5) Chain‐register 4D MT → co-reg native → mFFE → PAM50
  # ----------------------------------------------------------
  fslsplit "${ORIENTED[MT]}" "${id}_MT_dyn_"
  dyn_ref="${id}_MT_dyn_0008.nii.gz"
  sct_deepseg_sc -i "$dyn_ref" -c t2 -o "${dyn_ref%.*}_sc.nii.gz"
  cp "$dyn_ref" "${id}_MT_dyn_0008_reg.nii.gz"

  dyn_ref_sc="${dyn_ref%.*}_sc.nii.gz"

  sct_register_multimodal -i "$dyn_ref" \
                          -d "${ORIENTED[MFFE]}" \
                          -iseg "${dyn_ref%.*}_sc.nii.gz" \
                          -dseg "${id}_mFFE_sc.nii.gz" \
                          -param step=1,type=im,algo=rigid,metric=MI:step=2,type=im,algo=slicereg,metric=MI \
                          -x spline \
                          -o "${id}_MT_dyn_0008_to_mFFE.nii.gz" \
                         -owarp "${id}_warp_MTref2MFFE.nii.gz" \
                         -owarpinv "${id}_warp_MFFE2MTref.nii.gz"
  
  sct_apply_transfo -i "${id}_MT_dyn_0008_to_mFFE.nii.gz" \
                    -d "$SCT_DIR/data/PAM50/template/PAM50_t2.nii.gz" \
                    -w "${id}_warp_anat2template_native.nii.gz" \
                    -x spline \
                    -o "${id}_MT_dyn_0008_to_PAM.nii.gz"

  # Chain-register dynamics 0000–0007 (ascending) to 0008:
  prev="$dyn_ref"
  for i in $(seq -f "%04g" 0 7); do
    this_dyn="${id}_MT_dyn_${i}.nii.gz"
    out_reg="${id}_MT_dyn_${i}_reg.nii.gz"
    echo "    Registering $this_dyn → $prev ..."
    sct_register_multimodal -i "$this_dyn" \
                            -d "$prev" \
                            -dseg "$dyn_ref_sc" \
                            -param step=1,type=im,algo=rigid,metric=MI:step=2,type=im,algo=slicereg,metric=MI \
                            -x spline \
                            -o "$out_reg"
    rm -f warp_"${this_dyn%.*}"2"${prev%.*}".nii.gz \
          "${out_reg%.*}_inv.nii.gz"
    prev="$out_reg"

    sct_apply_transfo -i "$out_reg" \
                      -d "${ORIENTED[MFFE]}" \
                      -w "${id}_warp_MTref2MFFE.nii.gz" \
                      -x spline \
                      -o "${id}_MT_dyn_${i}_to_mFFE.nii.gz"

    sct_apply_transfo -i "${id}_MT_dyn_${i}_to_mFFE.nii.gz" \
                      -d "$SCT_DIR/data/PAM50/template/PAM50_t2.nii.gz" \
                      -w "${id}_warp_anat2template_native.nii.gz" \
                      -x spline \
                      -o "${id}_MT_dyn_${i}_to_PAM.nii.gz"
  done

  # Chain-register dynamics 0009–0015 (ascending) to 0008:
  prev="$dyn_ref"
  for i in $(seq -f "%04g" 9 15); do
    this_dyn="${id}_MT_dyn_${i}.nii.gz"
    out_reg="${id}_MT_dyn_${i}_reg.nii.gz"
    echo "    Registering $this_dyn → $prev ..."
    sct_register_multimodal -i "$this_dyn" \
                            -d "$prev" \
                            -dseg "$dyn_ref_sc" \
                            -param step=1,type=im,algo=rigid,metric=MI:step=2,type=im,algo=slicereg,metric=MI \
                            -x spline \
                            -o "$out_reg"
    prev="$out_reg"

    sct_apply_transfo -i "$out_reg" \
                      -d "${ORIENTED[MFFE]}" \
                      -w "${id}_warp_MTref2MFFE.nii.gz" \
                      -x spline \
                      -o "${id}_MT_dyn_${i}_to_mFFE.nii.gz"

    sct_apply_transfo -i "${id}_MT_dyn_${i}_to_mFFE.nii.gz" \
                      -d "$SCT_DIR/data/PAM50/template/PAM50_t2.nii.gz" \
                      -w "${id}_warp_anat2template_native.nii.gz" \
                      -x spline \
                      -o "${id}_MT_dyn_${i}_to_PAM.nii.gz"
  done

  rm *_inv.nii.gz

  # Step 6: Merge to final MT 4D file in PAM50 space
  fslmerge -t "${id}_registered_MT_native.nii.gz" $(for i in $(seq -f "%04g" 0 15); do echo "${id}_MT_dyn_${i}_to_mFFE.nii.gz"; done)
  fslmerge -t "${id}_registered_MT_PAM.nii.gz" $(for i in $(seq -f "%04g" 0 15); do echo "${id}_MT_dyn_${i}_to_PAM.nii.gz"; done)

  echo "Finished subject $subj (ID=$id). Files saved in $PROCESSED_DIR/$id"
  counter=$((counter + 1))
  cd "$BASE_DIR"
done

echo "All subjects processed."
