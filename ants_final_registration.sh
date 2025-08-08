#!/usr/bin/env bash
# ants_final_registration.sh
#
# Run after sct_final_registration.sh.
# Processes BOTH Processed/MS and Processed/HC.
# For each subject folder <ID>, generates:
#   <ID>_registered_{MT,MFA,T2AX,T1,B0,B1,LYMPH}_ANTS.nii.gz

# The input images are assumed to be **RPI-oriented** already, except
# the raw lymph mask (<ID>_LYMPH.nii) which the script flips on Y.


set -euo pipefail

THREADS=4                # cores for antsRegistrationSyNQuick
MT_REF_IDX=0008          # MT dynamic (0-based) chosen as reference
ROOT="$(pwd)"            # run from the Lymph folder
MS_DIR="$ROOT/Processed/MS"
HC_DIR="$ROOT/Processed/HC"

###############################################################################
process_subject () {
  local BASEDIR="$1"
  local ID="$2"
  local GROUP
  GROUP="$(basename "$BASEDIR")"   # "MS" or "HC"

  echo -e "\n================  $GROUP / $ID  ================"

  cd "$BASEDIR/$ID" || { echo "→ skip ($GROUP): folder missing"; return; }

  # ---------- inputs ----------
  local MT="${ID}_MT.nii.gz"
  local MFA="${ID}_MFA.nii.gz"
  local T1="${ID}_T1.nii.gz"
  local T2AX="${ID}_T2AX.nii.gz"
  local B0="${ID}_B0.nii.gz"
  local B1="${ID}_B1.nii.gz"
  local LYM_RAW="${ID}_LYMPH.nii"

  for f in "$MT" "$MFA" "$T1" "$T2AX" "$B0" "$B1" "$LYM_RAW"; do
    [[ -f $f ]] || { echo "→ skip ($GROUP/$ID): missing $f"; return; }
  done
  [[ $(fslval "$MT" dim4) -eq 16 ]] || { echo "→ skip ($GROUP/$ID): MT not 16 vols"; return; }

  # ---------- orient lymph mask ----------
  fslswapdim "$LYM_RAW" x -y z "${ID}_LYMPH_tmp.nii.gz"
  fslcpgeom  "$T2AX"    "${ID}_LYMPH_tmp.nii.gz"
  mv "${ID}_LYMPH_tmp.nii.gz" "${ID}_LYMPH_oriented_ANTS.nii.gz"

  # ---------- MT intra-dyn rigid ----------
  local MT_SPLIT="${ID}_MTdynANTS_"
  rm -f ${MT_SPLIT}* 2>/dev/null || true
  fslsplit "$MT" "$MT_SPLIT" -t

  local MT_REF_VOL="${ID}_MT_ref_ANTS.nii.gz"
  cp "${MT_SPLIT}${MT_REF_IDX}.nii.gz" "$MT_REF_VOL"

  for idx in $(seq -f "%04g" 0 15); do
    local mov="${MT_SPLIT}${idx}.nii.gz"
    local out="${MT_SPLIT}${idx}_reg_ANTS.nii.gz"
    if [[ $mov == "${MT_SPLIT}${MT_REF_IDX}.nii.gz" ]]; then
      cp "$mov" "$out"
      continue
    fi

    antsRegistrationSyNQuick.sh -d 3 -f "$MT_REF_VOL" -m "$mov" \
        -o "${MT_SPLIT}${idx}_ANTS_" -t r -n "$THREADS"

    antsApplyTransforms -d 3 -i "$mov" -r "$MT_REF_VOL" \
        -o "$out" -t "${MT_SPLIT}${idx}_ANTS_0GenericAffine.mat"
  done
  fslmerge -t "${ID}_registered_MT_ANTS.nii.gz" ${MT_SPLIT}*_reg_ANTS.nii.gz
  rm -f ${MT_SPLIT}* || true

  # ---------- MFA intra-dyn rigid ----------
  local MFA_SPLIT="${ID}_MFAdynANTS_"
  rm -f ${MFA_SPLIT}* 2>/dev/null || true
  fslsplit "$MFA" "$MFA_SPLIT" -t

  for mov in ${MFA_SPLIT}*.nii.gz; do
    [[ $mov == "${MFA_SPLIT}0000.nii.gz" ]] && continue
    antsRegistrationSyNQuick.sh -d 3 -f "${MFA_SPLIT}0000.nii.gz" -m "$mov" \
        -o "${mov%.nii.gz}_ANTS_" -t r -n "$THREADS"

    antsApplyTransforms -d 3 -i "$mov" -r "${MFA_SPLIT}0000.nii.gz" \
        -o "${mov%.nii.gz}_reg_ANTS.nii.gz" \
        -t "${mov%.nii.gz}_ANTS_0GenericAffine.mat"
  done
  cp "${MFA_SPLIT}0000.nii.gz" "${MFA_SPLIT}0000_reg_ANTS.nii.gz"
  fslmerge -t "${ID}_registered_MFA_ANTS.nii.gz" ${MFA_SPLIT}*_reg_ANTS.nii.gz
  rm -f ${MFA_SPLIT}* || true

  # ---------- single-volume mods → MT ----------
  for MOD in T2AX T1 B0 B1; do
    local SRC
    SRC=$(eval echo \$$MOD)
    antsRegistrationSyNQuick.sh -d 3 -f "$MT_REF_VOL" -m "$SRC" \
        -o "${ID}_${MOD}_to_MT_ANTS_" -t a -n "$THREADS"

    antsApplyTransforms -d 3 -i "$SRC" -r "$MT_REF_VOL" \
        -o "${ID}_registered_${MOD}_ANTS.nii.gz" \
        -t "${ID}_${MOD}_to_MT_ANTS_0GenericAffine.mat"
  done

  # ---------- warp lymph mask ----------
  antsApplyTransforms -d 3 \
      -i "${ID}_LYMPH_oriented_ANTS.nii.gz" -r "$MT_REF_VOL" \
      -o "${ID}_registered_LYMPH_ANTS.nii.gz" \
      -t "${ID}_T2AX_to_MT_ANTS_0GenericAffine.mat" \
      -n NearestNeighbor

  rm -f "$MT_REF_VOL"
  echo "✔ finished $GROUP / $ID"
}

###############################################################################
# Loop through every numeric folder in Processed/MS and Processed/HC
###############################################################################
for BASE in "$MS_DIR" "$HC_DIR"; do
  [[ -d "$BASE" ]] || continue
  for path in "$BASE"/*/; do
    [[ -d "$path" ]] || continue
    subj="$(basename "$path")"
    [[ $subj =~ ^[0-9]+$ ]] || continue
    process_subject "$BASE" "$subj"
  done
done
