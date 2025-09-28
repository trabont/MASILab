#!/usr/bin/env bash
# sct_final_registration.sh
# Retrieves all PreProcessed Subject's files are registers them into MT Space
set -euo pipefail

# ------------------------------------------------------------------
# 0.  FOLDER SET UP (process BOTH MS and HC)
# ------------------------------------------------------------------
BASE="$(pwd)"
PREPROCESSED_BASE="$BASE/PreProcessed"
PROCESSED_BASE="$BASE/Processed"
GROUPS=(MS HC)

# -- DEFINE WHICH IDS TO SKIP --------------------------------------
SKIP_IDS=(1) # change as needed.

# ---- QC root (single place for all subjects) ----------------------
QC_ROOT="$BASE/qc"
export QC_ROOT
mkdir -p "$QC_ROOT"

# ---- Choose vertebral levels for cord-specific analyses --------------
choose_levels() {
  local lbl="$1"
  read -r minL maxL < <(fslstats "$lbl" -R | awk '{printf "%d %d\n", int($1+0.5), int($2+0.5)}')
  (( minL < 1 )) && minL=1
  local pref=(2 3 4 5 6 7 8)
  local avail=()
  for v in "${pref[@]}"; do
    (( v >= minL && v <= maxL )) && avail+=("$v")
  done
  if (( ${#avail[@]} == 0 )) && (( maxL >= minL )); then
    avail=("$minL"); (( maxL > minL )) && avail+=("$maxL")
  fi
  local pick=(); local n=${#avail[@]}
  if   (( n >= 3 )); then pick=("${avail[0]}" "${avail[n/2]}" "${avail[n-1]}")
  elif (( n == 2 )); then pick=("${avail[0]}" "${avail[1]}")
  elif (( n == 1 )); then pick=("${avail[0]}"); fi
  (IFS=,; echo "${pick[*]}")
}

# ------------------------------------------------------------------
# 1.  PER SUBJECT FUNCTION
# ------------------------------------------------------------------
registration_of_subjects() {
    subj_dir="$1"        # e.g., /.../PreProcessed/<GROUP>/<ID>
    PROCESSED_DIR="$2"   # e.g., /.../Processed/<GROUP>
    id="$(basename "$subj_dir")"

    # -- SKIP CHECK ------------------------------------------------
    for skip in "${SKIP_IDS[@]}"; do
        [[ $id == "$skip" ]] && {
            echo "=== Skipping subject $id ==="
            return
        }
    done

    echo "=== Processing subject $id ==="
    mkdir -p "$PROCESSED_DIR/$id"
    cd "$PROCESSED_DIR/$id"

    # ---------- helper to grab the newest match -------------------
    find_file() {
        local pat="$1"
        find "$subj_dir" -type f -iname "$pat" | sort -V | tail -n1
    }

    MT_FILE=$(find_file "*_PulseMT_16_dyn.nii*")
    B0_FILE=$(find_file "*_B0_Map_*_e2_real.nii*")
    B1_FILE=$(find_file "*_WIP_B1_Map_Yarnykh_CLEAR_e2.nii*")
    MFA_FILE=$(find_file "*_MFA_*.nii*")
    T1_FILE=$(find_file "*_Clinical_T1W_TSE.nii*")
    T2AX_FILE=$(find_file "*_Clinical_T2W_TSE.nii*")
    T2SAG_FILE=$(find_file "*_DRIVE_CLEAR.nii*")
    MFFE_FILE=$(find_file "*_mFFE_0.65_14slice_e1.nii*")
    LYMPH_RAW=$(find_file "*ubMask.nii")

    #--------------------- warn if missing -------------------------
    for var in MT_FILE B0_FILE B1_FILE MFA_FILE T1_FILE T2AX_FILE T2SAG_FILE MFFE_FILE LYMPH_RAW; do
        [[ -n "${!var}" ]] || echo "WARNING: $var not found"
    done

    #--------------------- orient everything RPI -------------------
    declare -A INPUTS=(
      [MFA]="$MFA_FILE"   [MT]="$MT_FILE"     [B0]="$B0_FILE"
      [B1]="$B1_FILE"     [T1]="$T1_FILE"     [T2AX]="$T2AX_FILE"
      [T2SAG]="$T2SAG_FILE" [MFFE]="$MFFE_FILE"
    )

    declare -A ORIENTED
    for key in "${!INPUTS[@]}"; do
        src="${INPUTS[$key]}"
        [[ -n $src ]] || { echo "skip $key (not found)"; continue; }
        out="${id}_${key}.nii.gz"
        sct_image -i "$src" -setorient RPI -o "$out"
        ORIENTED[$key]="$out"
    done

    declare -A INPUTS=(
      [LYMPH]="$LYMPH_RAW"
    )
    [[ -n "${INPUTS[LYMPH]}" ]] && cp "${INPUTS[LYMPH]}" "${id}_LYMPH.nii"

    # ------------------------------------------------------------------
    # 2.  Registration using SCT for cord specific images
    # ------------------------------------------------------------------
    # ========================= mFFE → MT Ref ==========================
    # Split MT dynamics, pick middle as ref
    fslsplit "${ORIENTED[MT]}" "${id}_MT_dyn_"
    dyn_ref="${id}_MT_dyn_0008.nii.gz"
    cp "$dyn_ref" "${id}_MT_dyn_0008_reg.nii.gz"
    # Segment MT ref
    sct_deepseg spinalcord -i "$dyn_ref" -c t2 -o "${id}_MT_dyn_0008_sc.nii.gz"
    cp "${id}_MT_dyn_0008_sc.nii.gz" "${id}_MT_dyn_0008_reg_sc.nii.gz"
    # QC: dyn_ref seg
    sct_qc -i "$dyn_ref" -s "${id}_MT_dyn_0008_sc.nii.gz" -p sct_deepseg_sc -qc "$QC_ROOT" -qc-subject $id

    # Create mask around segmentation
    sct_create_mask -i "$dyn_ref" -p centerline,"${id}_MT_dyn_0008_sc.nii.gz" -size 35mm -o "${id}_MT_dyn_0008_mask.nii.gz"

    # mFFE → MT (dyn_ref)
    sct_register_multimodal -i "${ORIENTED[MFFE]}" -d "$dyn_ref" \
        -param step=1,type=im,algo=rigid,metric=MI,iter=80:step=2,type=im,algo=slicereg,metric=MI,iter=60,poly=2 \
        -m "${id}_MT_dyn_0008_mask.nii.gz" \
        -x spline \
        -o "${id}_registered_mFFE_1.nii.gz" \
        -owarp "${id}_warp_mFFE2MT.nii.gz" \
        -owarpinv "${id}_warp_MT2mFFE.nii.gz"
    # QC: mFFE→MT reg (use dest seg)
    sct_qc -i "$dyn_ref" -d "${id}_registered_mFFE_1.nii.gz" -s "${id}_MT_dyn_0008_sc.nii.gz" -p sct_register_multimodal -qc "$QC_ROOT" -qc-subject $id

    # Segmentations of mFFE in MT space
    sct_deepseg spinalcord -i "${id}_registered_mFFE_1.nii.gz" -c t2 -o "${id}_registered_mFFE_1_sc.nii.gz"
    # Segmentation/Mask of mFFE in native space
    sct_deepseg spinalcord -i "${ORIENTED[MFFE]}" -c t2 -o "${id}_mFFE_sc.nii.gz"
    sct_create_mask -i "${ORIENTED[MFFE]}" -p centerline,"${id}_mFFE_sc.nii.gz" -size 35mm -o "${id}_mFFE_mask.nii.gz"
    fslmaths "${id}_mFFE_mask.nii.gz" -thr 0.5 -bin "${id}_mFFE_mask.nii.gz"
    # QC: mFFE seg
    sct_qc -i "${ORIENTED[MFFE]}" -s "${id}_mFFE_sc.nii.gz" -p sct_deepseg_sc -qc "$QC_ROOT" -qc-subject $id

    # ==================== T2SAG labels (native) ====================
    # Segment T2SAG & label vertebrae
    sct_deepseg spinalcord -i "${ORIENTED[T2SAG]}" -c t2 -o "${id}_T2SAG_sc.nii.gz"
    sct_label_vertebrae -i "${ORIENTED[T2SAG]}" -s "${id}_T2SAG_sc.nii.gz" -c t2
    # QC: vertebral labels on T2SAG
    sct_qc -i "${ORIENTED[T2SAG]}" -s "${id}_T2SAG_sc_labeled.nii.gz" -p sct_label_vertebrae -qc "$QC_ROOT" -qc-subject $id

    levels="$(choose_levels "${id}_T2SAG_sc_labeled.nii.gz")"
    sct_label_utils -i "${id}_T2SAG_sc_labeled.nii.gz" -vert-body "$levels"

    # ==================== T2SAG → mFFE (rigid → slicereg) ====================
    sct_register_multimodal \
      -i "${ORIENTED[T2SAG]}" -d "${ORIENTED[MFFE]}" \
      -param step=1,type=im,algo=rigid,metric=MI,iter=80 \
      -x linear \
      -o       "${id}_T2SAG_rigid2mFFE.nii.gz" \
      -owarp   "${id}_warp_T2SAGrigid2mFFE.nii.gz" \
      -owarpinv "${id}_warp_mFFE2T2SAGrigid.nii.gz"

    sct_register_multimodal \
      -i "${ORIENTED[T2SAG]}" -d "${ORIENTED[MFFE]}" -m "${id}_mFFE_mask.nii.gz" \
      -initwarp    "${id}_warp_T2SAGrigid2mFFE.nii.gz" \
      -initwarpinv "${id}_warp_mFFE2T2SAGrigid.nii.gz" \
      -param step=1,type=im,algo=slicereg,metric=MI,iter=60,poly=2 \
      -x linear \
      -o       "${id}_T2SAG_mFFE.nii.gz" \
      -owarp   "${id}_warp_T2SAG2mFFE.nii.gz" \
      -owarpinv "${id}_warp_mFFE2T2SAG.nii.gz"

    sct_qc -i "${ORIENTED[MFFE]}" -d "${id}_T2SAG_mFFE.nii.gz" -s "${id}_mFFE_sc.nii.gz" -p sct_register_multimodal -qc "$QC_ROOT" -qc-subject $id

    # Segment registered T2SAG in mFFE space
    sct_deepseg spinalcord -i "${id}_T2SAG_mFFE.nii.gz" -c t2 -o "${id}_T2SAG_mFFE_sc.nii.gz"

    # bring labeled T2SAG cord into mFFE space
    sct_apply_transfo -i "${id}_T2SAG_sc_labeled.nii.gz" -d "${ORIENTED[MFFE]}" -w "${id}_warp_T2SAG2mFFE.nii.gz" -o "${id}_mFFE_sc_labeled.nii.gz" -x nn
    sct_qc -i "${ORIENTED[MFFE]}" -s "${id}_mFFE_sc_labeled.nii.gz" -p sct_label_utils -qc "$QC_ROOT" -qc-subject $id
    levels="$(choose_levels "${id}_mFFE_sc_labeled.nii.gz")"
    sct_label_utils -i "${id}_mFFE_sc_labeled.nii.gz" -vert-body "$levels" -o "${id}_mFFE_labels.nii.gz"

    # ==================== mFFE → PAM50 (nonlinear) ====================
    sct_register_to_template -i "${ORIENTED[MFFE]}" -s "${id}_mFFE_sc.nii.gz" -l "${id}_mFFE_labels.nii.gz" -c t2 \
        -param step=1,type=seg,algo=rigid,iter=80:step=2,type=seg,algo=affine,iter=80:step=3,type=seg,algo=syn,iter=60,shrink=2,smooth=1
    # QC: register_to_template (mFFE)
    sct_qc -i "${ORIENTED[MFFE]}" -s "${id}_mFFE_sc.nii.gz" -p sct_register_to_template -qc "$QC_ROOT" -qc-subject $id

    # Invert warp (template → mFFE)
    sct_warp_template -d "${ORIENTED[MFFE]}" -w warp_template2anat.nii.gz
    # QC: warp_template (mFFE as dest)
    sct_qc -i "${ORIENTED[MFFE]}" -p sct_warp_template -qc "$QC_ROOT" -qc-subject $id

    # ================== Tissues (PAM50) in MT space ====================
    sct_apply_transfo -i label/template/PAM50_wm.nii.gz  -d "$dyn_ref" -w "${id}_warp_mFFE2MT.nii.gz" -o temp_registered_wm.nii.gz  -x linear
    sct_apply_transfo -i label/template/PAM50_gm.nii.gz  -d "$dyn_ref" -w "${id}_warp_mFFE2MT.nii.gz" -o temp_registered_gm.nii.gz  -x linear
    sct_apply_transfo -i label/template/PAM50_csf.nii.gz -d "$dyn_ref" -w "${id}_warp_mFFE2MT.nii.gz" -o temp_registered_csf.nii.gz -x linear
    sct_maths -i temp_registered_gm.nii.gz  -bin 0.5 -o "${id}_registered_GM.nii.gz"
    sct_maths -i temp_registered_csf.nii.gz -bin 0.5 -o "${id}_registered_CSF.nii.gz"
    fslmaths temp_registered_wm.nii.gz -sub temp_registered_gm.nii.gz -thr 0 -sub temp_registered_csf.nii.gz -thr 0 -bin "${id}_registered_WM.nii.gz"
    # QC: tissues in MT space (labels)
    sct_qc -i "$dyn_ref" -s "${id}_registered_GM.nii.gz"  -p sct_label_utils -qc "$QC_ROOT" -qc-subject $id
    sct_qc -i "$dyn_ref" -s "${id}_registered_CSF.nii.gz" -p sct_label_utils -qc "$QC_ROOT" -qc-subject $id
    sct_qc -i "$dyn_ref" -s "${id}_registered_WM.nii.gz"  -p sct_label_utils -qc "$QC_ROOT" -qc-subject $id

    # ==================== T1 & T2AX → mFFE ====================
    # T2AX → mFFE (rigid → slicereg)
    sct_register_multimodal \
      -i "${ORIENTED[T2AX]}" -d "${ORIENTED[MFFE]}" \
      -param step=1,type=im,algo=rigid,metric=MI,iter=80 \
      -x spline \
      -o       "${id}_T2AX_rigid2mFFE.nii.gz" \
      -owarp   "${id}_warp_T2AXrigid2mFFE.nii.gz" \
      -owarpinv "${id}_warp_mFFE2T2AXrigid.nii.gz"

    # Slicereg (use mask, init with fwd+inv; save fwd+inv)
    sct_register_multimodal \
      -i "${ORIENTED[T2AX]}" -d "${ORIENTED[MFFE]}" -m "${id}_mFFE_mask.nii.gz" \
      -initwarp    "${id}_warp_T2AXrigid2mFFE.nii.gz" \
      -initwarpinv "${id}_warp_mFFE2T2AXrigid.nii.gz" \
      -param step=1,type=im,algo=slicereg,metric=MI,iter=60,poly=2 \
      -x spline \
      -o       "${id}_T2AX_mFFE.nii.gz" \
      -owarp   "${id}_warp_T2AX2mFFE.nii.gz" \
      -owarpinv "${id}_warp_mFFE2T2AX.nii.gz"

    # QC for T2AX registration (needs fixed-image seg)
    sct_qc -i "${ORIENTED[MFFE]}" -d "${id}_T2AX_mFFE.nii.gz" -s "${id}_mFFE_sc.nii.gz" -p sct_register_multimodal -qc "$QC_ROOT" -qc-subject $id

    # T1 → mFFE (rigid → slicereg)
    sct_register_multimodal \
      -i "${ORIENTED[T1]}" -d "${ORIENTED[MFFE]}" \
      -param step=1,type=im,algo=rigid,metric=MI,iter=80 \
      -x spline \
      -o       "${id}_T1_rigid2mFFE.nii.gz" \
      -owarp   "${id}_warp_T1rigid2mFFE.nii.gz" \
      -owarpinv "${id}_warp_mFFE2T1rigid.nii.gz"

    sct_register_multimodal \
      -i "${ORIENTED[T1]}" -d "${ORIENTED[MFFE]}" -m "${id}_mFFE_mask.nii.gz" \
      -initwarp    "${id}_warp_T1rigid2mFFE.nii.gz" \
      -initwarpinv "${id}_warp_mFFE2T1rigid.nii.gz" \
      -param step=1,type=im,algo=slicereg,metric=MI,iter=60,poly=2 \
      -x spline \
      -o       "${id}_T1_mFFE.nii.gz" \
      -owarp   "${id}_warp_T12mFFE.nii.gz" \
      -owarpinv "${id}_warp_mFFE2T1.nii.gz"

    sct_qc -i "${ORIENTED[MFFE]}" -d "${id}_T1_mFFE.nii.gz" -s "${id}_mFFE_sc.nii.gz" -p sct_register_multimodal -qc "$QC_ROOT" -qc-subject $id

    # ==================== T1 & T2AX → MT ====================
    for MOD in T1 T2AX; do
      IMG="${id}_${MOD}_mFFE.nii.gz"
      [[ -f "$IMG" ]] || { echo "[${id}] Missing $IMG"; return 1; }
      OUT="${id}_registered_${MOD}_1.nii.gz"
      OUT_SC="${id}_registered_${MOD}_1_sc.nii.gz"
      sct_apply_transfo -i "$IMG" -d "$dyn_ref" -w "${id}_warp_mFFE2MT.nii.gz" -o "$OUT" -x linear
      sct_deepseg spinalcord -i "$OUT" -c t2 -o "$OUT_SC"
      # QC: modality moved into MT space (needs DEST seg → DYN_SC)
      sct_qc -i "$dyn_ref" -d "$OUT" -s "${id}_MT_dyn_0008_sc.nii.gz" -p sct_register_multimodal -qc "$QC_ROOT" -qc-subject $id
    done

    # ==================== CoRegister MT ====================
    # Chain registration of MT dynamics to middle dynamic (dyn_ref)
    prev="$dyn_ref"
    for i in $(seq -f "%04g" 0 7); do
        in="${id}_MT_dyn_${i}.nii.gz"; out="${id}_MT_dyn_${i}_reg.nii.gz"
        sct_register_multimodal -i "$in" -d "$prev" -m "${id}_MT_dyn_0008_mask.nii.gz" \
        -param step=1,type=im,algo=rigid,metric=MI:step=2,type=im,algo=slicereg,metric=MI \
        -x spline -o "$out"
        #QC: MT dyn chain reg
        sct_qc -i "$prev" -d "$out" -s "${id}_MT_dyn_0008_sc.nii.gz" -p sct_register_multimodal -qc "$QC_ROOT" -qc-subject $id
        prev="$out"
    done
    prev="$dyn_ref"
    for i in $(seq -f "%04g" 9 15); do
        in="${id}_MT_dyn_${i}.nii.gz"; out="${id}_MT_dyn_${i}_reg.nii.gz"
        sct_register_multimodal -i "$in" -d "$prev" -m "${id}_MT_dyn_0008_mask.nii.gz" \
        -param step=1,type=im,algo=rigid,metric=MI:step=2,type=im,algo=slicereg,metric=MI \
        -x spline -o "$out"
        #QC: MT dyn chain reg
        sct_qc -i "$prev" -d "$out" -s "${id}_MT_dyn_0008_sc.nii.gz" -p sct_register_multimodal -qc "$QC_ROOT" -qc-subject $id
        prev="$out"
    done
    fslmerge -t "${id}_registered_MT_1.nii.gz" ${id}_MT_dyn_*_reg.nii.gz
    sct_deepseg spinalcord -i "${id}_MT_dyn_0000_reg.nii.gz" -c t2 -o "${id}_MT_dyn_0000_reg_sc.nii.gz"

    # ==================== B0 & B1 → MT ====================
    # Register B0 and B1 to MT (no mask)
    for k in B0 B1; do
        [[ -n "${ORIENTED[$k]:-}" ]] || continue
        sct_register_multimodal -i "${ORIENTED[$k]}" -d "$dyn_ref" \
        -param step=1,type=im,algo=rigid,metric=MI:step=2,type=im,algo=slicereg,metric=MI \
        -x spline -o "${id}_registered_${k}_1.nii.gz"
        # QC: k→MT reg
        sct_qc -i "$dyn_ref" -d "${id}_registered_${k}_1.nii.gz" -s "${id}_MT_dyn_0008_sc.nii.gz" -p sct_register_multimodal -qc "$QC_ROOT" -qc-subject $id
    done

    # ==================== MFA → MT ====================
    # Register MFA reference to MT, then apply to others
    if [[ -n "${ORIENTED[MFA]:-}" ]]; then
      n_dyn=$(fslhd "${ORIENTED[MFA]}" | awk '/^dim4/ {print $2}')
      fslsplit "${ORIENTED[MFA]}" "${id}_MFA_dyn_" -t
      mfa_ref="${id}_MFA_dyn_0000.nii.gz"
      cp "$mfa_ref" "${id}_MFA_dyn_0000_reg.nii.gz"
      sct_deepseg spinalcord -i "${id}_MFA_dyn_0000.nii.gz" -c t1 -o "${id}_MFA_0000_sc.nii.gz"
      sct_create_mask -i "$mfa_ref" -p centerline,"${id}_MFA_0000_sc.nii.gz" -size 35mm -o "${id}_MFA_dyn_0000_mask.nii.gz"

      sct_register_multimodal -i "$mfa_ref" -d "$dyn_ref" \
          -m "${id}_MT_dyn_0008_mask.nii.gz" \
          -param step=1,type=im,algo=rigid,metric=MI:step=2,type=im,algo=slicereg,metric=MI \
          -x spline \
          -o "${id}_MFA_dyn_0000_in_MT.nii.gz" \
          -owarp "${id}_MFA_1_warp.nii.gz"
      # QC: MFA ref → MT
      sct_qc -i "$dyn_ref" -d "${id}_MFA_dyn_0000_in_MT.nii.gz" -s "${id}_MT_dyn_0008_sc.nii.gz" -p sct_register_multimodal -qc "$QC_ROOT" -qc-subject $id

      for i in $(seq -f "%04g" 1 $((n_dyn-1))); do
          in="${id}_MFA_dyn_${i}.nii.gz"; out="${id}_MFA_dyn_${i}_reg.nii.gz"
          sct_register_multimodal -i "$in" -d "$mfa_ref" -m "${id}_MFA_dyn_0000_mask.nii.gz" \
          -param step=1,type=im,algo=rigid,metric=MI:step=2,type=im,algo=slicereg,metric=MI \
          -x spline -o "$out"
          #QC: MFA dyn chain reg (needs fixed-image seg)
          sct_qc -i "$mfa_ref" -d "$out" -s "${id}_MFA_0000_sc.nii.gz" -p sct_register_multimodal -qc "$QC_ROOT" -qc-subject $id
          sct_apply_transfo -i "$out" -d "$dyn_ref" -w "${id}_MFA_1_warp.nii.gz" -x spline -o "${id}_MFA_dyn_${i}_in_MT.nii.gz"
          # QC: MFA dyn → MT
          sct_qc -i "$dyn_ref" -d "${id}_MFA_dyn_${i}_in_MT.nii.gz" -s "${id}_MT_dyn_0008_sc.nii.gz" -p sct_register_multimodal -qc "$QC_ROOT" -qc-subject $id
      done
      sct_deepseg spinalcord -i "${id}_MFA_dyn_0000_in_MT.nii.gz" -c t1 -o "${id}_MFA_dyn_0000_in_MT_sc.nii.gz"
      sct_qc -i "${id}_MFA_dyn_0000_in_MT.nii.gz" -s "${id}_MFA_dyn_0000_in_MT_sc.nii.gz" -p sct_deepseg_sc -qc "$QC_ROOT" -qc-subject $id
      fslmerge -t "${id}_registered_MFA_1.nii.gz" ${id}_MFA_dyn_*_in_MT.nii.gz
    fi

    echo "Finished processing subject $id"
    rm -f *_inv.nii.gz *.json || true
}

export -f registration_of_subjects

# ------------------------------------------------------------------
# 4.  LAUNCH (MS and HC)
# ------------------------------------------------------------------
for GROUP in "${GROUPS[@]}"; do
    ROOT="$PREPROCESSED_BASE/$GROUP"
    OUTDIR="$PROCESSED_BASE/$GROUP"
    mkdir -p "$OUTDIR"
    [[ -d "$ROOT" ]] || continue

    for subj_path in "$ROOT"/*; do
        [[ -d $subj_path ]] || continue
        b=$(basename "$subj_path")
        [[ $b =~ ^[0-9]+$ ]] || continue
        registration_of_subjects "$subj_path" "$OUTDIR"
    done
done
