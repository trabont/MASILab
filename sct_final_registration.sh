# sct_final_registration.sh
# Retrieves all PreProcessed Subject's files are registers them into MT Space

#!/usr/bin/env bash
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

# ------------------------------------------------------------------
# 1.  PER SUBJECT FUNCTION
# ------------------------------------------------------------------
registration_of_subjects() {
    subj_dir="$1"        # e.g., /.../PreProcessed/<GROUP>/<ID>
    PROCESSED_DIR="$2"   # e.g., /.../Processed/<GROUP>
    id="$(basename "$subj_dir")"

    # -- SKIP CHECK ------------------------------------------------
    for skip in "${SKIP_IDS[@]}"; do
        [[ "$id" == "$skip" ]] && {
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
    cp "${INPUTS[LYMPH]}" "${id}_LYMPH.nii"

    # ------------------------------------------------------------------
    # 2.  Registration using SCT for cord specific images
    # ------------------------------------------------------------------
    # First register mFFE to MT (mid-dynamic as reference)
    fslsplit "${ORIENTED[MT]}" "${id}_MT_dyn_"
    dyn_ref="${id}_MT_dyn_0008.nii.gz"
    cp "$dyn_ref" "${id}_MT_dyn_0008_reg.nii.gz"
    sct_deepseg spinalcord -i "$dyn_ref" -c t2 -o "${id}_MT_dyn_0008_sc.nii.gz"
    cp "${id}_MT_dyn_0008_sc.nii.gz" "${id}_MT_dyn_0008_reg_sc.nii.gz"
    sct_create_mask -i "$dyn_ref" -p centerline,"${id}_MT_dyn_0008_sc.nii.gz" -size 35mm -o "${id}_MT_dyn_0008_mask.nii.gz"
    sct_register_multimodal -i "${ORIENTED[MFFE]}" -d "$dyn_ref" \
        -param step=1,type=im,algo=rigid,metric=MI:step=2,type=im,algo=slicereg,metric=MI \
        -m "${id}_MT_dyn_0008_mask.nii.gz" \
        -x spline \
        -o "${id}_registered_mFFE_1.nii.gz" \
        -owarp "${id}_warp_mFFE2MT.nii.gz" \
        -owarpinv "${id}_warp_MT2mFFE.nii.gz"
    # Register PAM50 template to mFFE space to get tissue segmentations
    sct_deepseg spinalcord -i "${id}_registered_mFFE_1.nii.gz" -c t2 -o "${id}_registered_mFFE_1_sc.nii.gz"
    sct_deepseg spinalcord -i "${ORIENTED[MFFE]}" -c t2 -o "${id}_mFFE_sc.nii.gz"
    sct_create_mask -i "${ORIENTED[MFFE]}" -p centerline,"${id}_mFFE_sc.nii.gz" -size 35mm -o "${id}_mFFE_mask.nii.gz"
    sct_deepseg spinalcord -i "${ORIENTED[T2SAG]}" -c t2 -o "${id}_T2SAG_sc.nii.gz"
    sct_label_vertebrae -i "${ORIENTED[T2SAG]}" -s "${id}_T2SAG_sc.nii.gz" -c t2
    sct_label_utils -i "${id}_T2SAG_sc_labeled.nii.gz" -vert-body 2,8
    sct_register_to_template -i "${ORIENTED[MFFE]}" -s "${id}_mFFE_sc.nii.gz" -l "${id}_T2SAG_sc_labeled_discs.nii.gz" -c t2 

    sct_warp_template -d "${ORIENTED[MFFE]}" -w warp_template2anat.nii.gz
    sct_apply_transfo -i label/template/PAM50_wm.nii.gz  -d "$dyn_ref" -w "${id}_warp_mFFE2MT.nii.gz" -o temp_registered_wm.nii.gz  -x nn
    sct_apply_transfo -i label/template/PAM50_gm.nii.gz  -d "$dyn_ref" -w "${id}_warp_mFFE2MT.nii.gz" -o temp_registered_gm.nii.gz  -x nn
    sct_apply_transfo -i label/template/PAM50_csf.nii.gz -d "$dyn_ref" -w "${id}_warp_mFFE2MT.nii.gz" -o temp_registered_csf.nii.gz -x nn
    sct_maths -i temp_registered_gm.nii.gz  -bin 0.5 -o "${id}_registered_GM.nii.gz"
    sct_maths -i temp_registered_csf.nii.gz -bin 0.5 -o "${id}_registered_CSF.nii.gz"
    fslmaths temp_registered_wm.nii.gz -sub "${id}_registered_CSF.nii.gz" -thr 0 -bin tmp_wm1.nii.gz
    fslmaths tmp_wm1.nii.gz -sub "${id}_registered_GM.nii.gz"  -thr 0 -bin "${id}_registered_WM.nii.gz"
    rm tmp_wm1.nii.gz

    # Then register all other modalities to MT space (using same mid-dynamic as reference)
    prev="$dyn_ref"
    for i in $(seq -f "%04g" 0 7); do
        in="${id}_MT_dyn_${i}.nii.gz"; out="${id}_MT_dyn_${i}_reg.nii.gz"
        sct_register_multimodal -i "$in" -d "$prev" -m "${id}_MT_dyn_0008_mask.nii.gz" \
        -param step=1,type=im,algo=rigid,metric=MI:step=2,type=im,algo=slicereg,metric=MI \
        -x spline -o "$out"
        prev="$out"
    done
    prev="$dyn_ref"
    for i in $(seq -f "%04g" 9 15); do
        in="${id}_MT_dyn_${i}.nii.gz"; out="${id}_MT_dyn_${i}_reg.nii.gz"
        sct_register_multimodal -i "$in" -d "$prev" -m "${id}_MT_dyn_0008_mask.nii.gz" \
        -param step=1,type=im,algo=rigid,metric=MI:step=2,type=im,algo=slicereg,metric=MI \
        -x spline -o "$out"
        prev="$out"
    done
    fslmerge -t "${id}_registered_MT_1.nii.gz" ${id}_MT_dyn_*_reg.nii.gz
    sct_deepseg spinalcord -i "${id}_MT_dyn_0000_reg.nii.gz" -c t2 -o "${id}_MT_dyn_0000_reg_sc.nii.gz"
    
    # Register T1 and T2AX to MT
    for k in T1 T2AX; do
        sct_register_multimodal -i "${ORIENTED[$k]}" -d "$dyn_ref" \
        -m "${id}_MT_dyn_0008_mask.nii.gz" \
        -param step=1,type=im,algo=rigid,metric=MI:step=2,type=im,algo=slicereg,metric=MI \
        -x spline -o "${id}_registered_${k}_1.nii.gz"
    done
    sct_deepseg spinalcord -i "${id}_registered_T2AX_1.nii.gz" -c t2 -o "${id}_registered_T2AX_1_sc.nii.gz"
    sct_deepseg spinalcord -i "${id}_registered_T1_1.nii.gz" -c t1 -o "${id}_registered_T1_1_sc.nii.gz"

    # Register B0 and B1 to MT (no mask)
    for k in B0 B1; do
        sct_register_multimodal -i "${ORIENTED[$k]}" -d "$dyn_ref" \
        -param step=1,type=im,algo=rigid,metric=MI:step=2,type=im,algo=slicereg,metric=MI \
        -x spline -o "${id}_registered_${k}_1.nii.gz"
    done

    # Register MFA reference to MT
    n_dyn=$(fslhd "${ORIENTED[MFA]}" | awk '/^dim4/ {print $2}')
    fslsplit "${ORIENTED[MFA]}" "${id}_MFA_dyn_" -t
    mfa_ref="${id}_MFA_dyn_0000.nii.gz"
    cp "$mfa_ref" "${id}_MFA_dyn_0000_reg.nii.gz"
    sct_deepseg spinalcord -i "${id}_MFA_dyn_0000.nii.gz" -c t1 -o "${id}_MFA_0000_sc.nii.gz
    sct_create_mask -i "$mfa_ref" -p centerline,"${id}_MFA_0000_sc.nii.gz" -size 35mm -o "${id}_MFA_dyn_0000_mask.nii.gz"
    # Co-register MFA reference and apply the same warp to all other dynamics to MT space
    sct_register_multimodal -i "$mfa_ref" -d "$dyn_ref" \
        -m "${id}_MT_dyn_0008_mask.nii.gz" \
        -param step=1,type=im,algo=rigid,metric=MI:step=2,type=im,algo=slicereg,metric=MI \
        -x spline \
        -o "${id}_MFA_dyn_0000_in_MT.nii.gz" \
        -owarp "${id}_MFA_1_warp.nii.gz"
    for i in $(seq -f "%04g" 1 $((n_dyn-1))); do
        in="${id}_MFA_dyn_${i}.nii.gz"; out="${id}_MFA_dyn_${i}_reg.nii.gz"
        sct_register_multimodal -i "$in" -d "$mfa_ref" -m "${id}_MFA_dyn_0000_mask.nii.gz" \
        -param step=1,type=im,algo=rigid,metric=MI:step=2,type=im,algo=slicereg,metric=MI \
        -x spline -o "$out"
        sct_apply_transfo -i "$out" -d "$dyn_ref" -w "${id}_MFA_1_warp.nii.gz" -x spline -o "${id}_MFA_dyn_${i}_in_MT.nii.gz"
    done
    sct_deepseg spinalcord -i "${id}_MFA_dyn_0000_in_MT.nii.gz" -c t1 -o "${id}_MFA_dyn_0000_in_MT_sc.nii.gz"
    fslmerge -t "${id}_registered_MFA_1.nii.gz" ${id}_MFA_dyn_*_in_MT.nii.gz

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
