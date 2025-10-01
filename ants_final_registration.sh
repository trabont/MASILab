#!/usr/bin/env bash
# ants_final_registration.sh
#
# Run after sct_final_registration.sh.
# Processes BOTH Processed/MS and Processed/HC.
# For each subject folder <ID>, generates:
#   <ID>_registered_{MT,MFA,MFFE,T2AX,T1,B0,B1,LYMPH}_ANTS.nii.gz

# The input images are assumed to be **RPI-oriented** already, except
# the raw lymph mask (<ID>_LYMPH.nii) which the script flips on Y.


set -euo pipefail

GROUPS=("HC" "MS")

# -------- helpers --------
need() {  # need <file>...
  local miss=0
  for f in "$@"; do
    [[ -f "$f" ]] || { echo "  MISSING: $f"; miss=1; }
  done
  return $miss
}

cleanup_outputs() {
  local ID="$1"
  # Remove only products we (re)create; keep raw inputs intact.
  rm -f \
    "${ID}_LYMPH_tmp.nii.gz" \
    "${ID}_LYMPH_oriented_ANTS.nii.gz" \
    "${ID}_registered_MT_ANTS.nii.gz" \
    "${ID}_registered_MFA_ANTS.nii.gz" \
    "${ID}_registered_mFFE_ANTS.nii.gz" \
    "${ID}_registered_T1_ANTS.nii.gz" \
    "${ID}_registered_T2AX_ANTS.nii.gz" \
    "${ID}_registered_B0_ANTS.nii.gz" \
    "${ID}_registered_B1_ANTS.nii.gz" \
    "${ID}_MT_ref_ANTS.nii.gz" \
    "${ID}_MFFE_to_MT_ANTS_"* \
    "${ID}_T1_to_MT_ANTS_"* \
    "${ID}_T2AX_to_MT_ANTS_"* \
    "${ID}_B0_to_MT_ANTS_"* \
    "${ID}_B1_to_MT_ANTS_"* \
    2>/dev/null || true

  # Split/rigid stacks (MT/MFA) and their intermediates
  rm -f "${ID}_MTdynANTS_"* 2>/dev/null || true
  rm -f "${ID}_MFAdynANTS_"* 2>/dev/null || true

  # Generic ANTs “*Warped*.nii.gz” that come from the prefixes above
  rm -f "${ID}"_*_ANTS_Warped.nii.gz "${ID}"_*_ANTS_InverseWarped.nii.gz 2>/dev/null || true

  # Misc json/inv files
  rm -f *.json *_inv.nii.gz 2>/dev/null || true
}

# -------- per-subject worker (subshell) --------
run_one() (
  set -euo pipefail
  shopt -s nullglob

  local subj_dir="$1"
  local ID; ID="$(basename "$subj_dir")"
  echo "=== [$ID] start ==="
  cd "$subj_dir"

  # ---------- inputs ----------
  local MT_REF_IDX="0008"           # 0-based index for ref vol

  local MT="${ID}_MT.nii.gz"
  local MFA="${ID}_MFA.nii.gz"
  local T1="${ID}_T1.nii.gz"
  local T2AX="${ID}_T2AX.nii.gz"
  local B0="${ID}_B0.nii.gz"
  local B1="${ID}_B1.nii.gz"
  local LYM_RAW="${ID}_LYMPH.nii"
  local MFFE="${ID}_MFFE.nii.gz"

  need "$MT" "$MFA" "$T1" "$T2AX" "$B0" "$B1" "$MFFE" "$LYM_RAW" || { echo "[$ID] ERROR: required inputs missing"; exit 1; }

  # Clean old outputs we’ll recreate
  cleanup_outputs "$ID"

  # ---------- orient lymph mask ----------
  fslswapdim "$LYM_RAW" x -y z "${ID}_LYMPH_tmp.nii.gz"
  fslcpgeom  "$T2AX"    "${ID}_LYMPH_tmp.nii.gz"
  mv "${ID}_LYMPH_tmp.nii.gz" "${ID}_LYMPH_oriented_ANTS.nii.gz"

  # ---------- MT intra-dyn rigid ----------
  MT_SPLIT="${ID}_MTdynANTS_"
  rm -f ${MT_SPLIT}* 2>/dev/null || true
  fslsplit "$MT" "$MT_SPLIT" -t

  MT_REF_VOL="${ID}_MT_ref_ANTS.nii.gz"
  cp "${MT_SPLIT}${MT_REF_IDX}.nii.gz" "$MT_REF_VOL"

  for idx in $(seq -f "%04g" 0 15); do
    mov="${MT_SPLIT}${idx}.nii.gz"
    out="${MT_SPLIT}${idx}_reg_ANTS.nii.gz"
    if [[ $mov == "${MT_SPLIT}${MT_REF_IDX}.nii.gz" ]]; then
      cp "$mov" "$out"
      continue
    fi

    antsRegistrationSyNQuick.sh -d 3 -f "$MT_REF_VOL" -m "$mov" \
      -o "${MT_SPLIT}${idx}_ANTS_" -t r -n 4

    antsApplyTransforms -d 3 -i "$mov" -r "$MT_REF_VOL" \
      -o "$out" \
      -t "${MT_SPLIT}${idx}_1Warp.nii.gz" \
      -t "${MT_SPLIT}${idx}_ANTS_0GenericAffine.mat"
  done
  fslmerge -t "${ID}_registered_MT_ANTS.nii.gz" ${MT_SPLIT}*_reg_ANTS.nii.gz
  rm -f ${MT_SPLIT}* || true

  # ---------- MFA intra-dyn rigid ----------
  MFA_SPLIT="${ID}_MFAdynANTS_"
  rm -f ${MFA_SPLIT}* 2>/dev/null || true
  fslsplit "$MFA" "$MFA_SPLIT" -t

  for mov in ${MFA_SPLIT}*.nii.gz; do
    [[ $mov == "${MFA_SPLIT}0000.nii.gz" ]] && continue
    antsRegistrationSyNQuick.sh -d 3 -f "${MFA_SPLIT}0000.nii.gz" -m "$mov" \
      -o "${mov%.nii.gz}_ANTS_" -t r -n 4

    antsApplyTransforms -d 3 -i "$mov" -r "${MFA_SPLIT}0000.nii.gz" \
      -o "${mov%.nii.gz}_reg_ANTS.nii.gz" \
      -t "${mov%.nii.gz}_1Warp.nii.gz" \
      -t "${mov%.nii.gz}_ANTS_0GenericAffine.mat"
  done
  cp "${MFA_SPLIT}0000.nii.gz" "${MFA_SPLIT}0000_reg_ANTS.nii.gz"
  fslmerge -t "${ID}_registered_MFA_ANTS.nii.gz" ${MFA_SPLIT}*_reg_ANTS.nii.gz
  rm -f ${MFA_SPLIT}* || true

  # ---------- mFFE to MT ----------
  antsRegistrationSyNQuick.sh -d 3 -f "$MT_REF_VOL" -m "$MFFE" \
    -o "${ID}_MFFE_to_MT_ANTS_" -t a -n 4

  antsApplyTransforms -d 3 -i "$MFFE" -r "$MT_REF_VOL" \
    -o "${ID}_registered_mFFE_ANTS.nii.gz" \
    -t "${ID}_MFFE_to_MT_ANTS_1Warp.nii.gz" \
    -t "${ID}_MFFE_to_MT_ANTS_0GenericAffine.mat"

  # ---------- single-volume mods (T1 & T2AX) → mFFE → MT ----------
  for MOD in T1 T2AX; do
    SRC="$(eval echo \$$MOD)"
    antsRegistrationSyNQuick.sh -d 3 -f "$MFFE" -m "$SRC" \
      -o "${ID}_${MOD}_to_MFFE_ANTS_" -t a -n 4

    antsApplyTransforms -d 3 -i "$SRC" -r "$MT_REF_VOL" \
      -o "${ID}_registered_${MOD}_ANTS.nii.gz" \
      -t "${ID}_${MOD}_to_MFFE_ANTS_1Warp.nii.gz" \
      -t "${ID}_${MOD}_to_MFFE_ANTS_0GenericAffine.mat" \
      -t "${ID}_MFFE_to_MT_ANTS_1Warp.nii.gz" \
      -t "${ID}_MFFE_to_MT_ANTS_0GenericAffine.mat"
  done

  # ---------- single-volume mods (B0 & B1) → MT ----------
  for MOD in B0 B1; do
    SRC="$(eval echo \$$MOD)"
    antsRegistrationSyNQuick.sh -d 3 -f "$MT_REF_VOL" -m "$SRC" \
      -o "${ID}_${MOD}_to_MT_ANTS_" -t a -n 4

    antsApplyTransforms -d 3 -i "$SRC" -r "$MT_REF_VOL" \
      -o "${ID}_registered_${MOD}_ANTS.nii.gz" \
      -t "${ID}_${MOD}_to_MT_ANTS_1Warp.nii.gz" \
      -t "${ID}_${MOD}_to_MT_ANTS_0GenericAffine.mat"
  done

  # ---------- warp lymph mask (via T2AX→mFFE→MT) ----------
  antsApplyTransforms -d 3 \
    -i "${ID}_LYMPH_oriented_ANTS.nii.gz" -r "$MT_REF_VOL" \
    -o "${ID}_registered_LYMPH_ANTS.nii.gz" \
    -t "${ID}_T2AX_to_MFFE_ANTS_1Warp.nii.gz" \
    -t "${ID}_T2AX_to_MFFE_ANTS_0GenericAffine.mat" \
    -t "${ID}_MFFE_to_MT_ANTS_1Warp.nii.gz" \
    -t "${ID}_MFFE_to_MT_ANTS_0GenericAffine.mat" \
    -n NearestNeighbor

  # Final tidy
  rm -f *.json *_inv.nii.gz 2>/dev/null || true

  echo "=== [$ID] done ==="
)

# -------- main --------
[[ -d "HC" && -d "MS" ]] || { echo "Run from Processed/ (must contain HC and MS folders)."; exit 1; }

declare -a OKS=()
declare -a FAILS=()

for G in "${GROUPS[@]}"; do
  [[ -d "$G" ]] || continue
  echo "--- Group: $G ---"
  for d in "$G"/*/ ; do
    [[ -d "$d" ]] || continue
    b="$(basename "$d")"
    [[ "$b" =~ ^[0-9]+$ ]] || { echo "skip non-numeric: $b"; continue; }
    if run_one "$d"; then
      OKS+=("$G/$b")
    else
      echo "!!! [$G/$b] FAILED — continuing"
      FAILS+=("$G/$b")
    fi
  done
done

echo
echo "================ SUMMARY ================"
echo "Succeeded: ${#OKS[@]}"
[[ ${#OKS[@]} -gt 0 ]] && printf '  %s\n' "${OKS[@]}"
echo "Failed:    ${#FAILS[@]}"
[[ ${#FAILS[@]} -gt 0 ]] && printf '  %s\n' "${FAILS[@]}"
echo "========================================="

# Exit nonzero if any failures (so CI can catch)
exit $(( ${#FAILS[@]} > 0 ? 2 : 0 ))

