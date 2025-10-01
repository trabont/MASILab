#!/usr/bin/env bash
# ants_final_registration.sh
# Run from: Processed/   (must contain HC/ and MS/)
# Nonlinear (SyN) everywhere -> always use 1Warp + 0GenericAffine
# Interp: Linear for intensities, NearestNeighbor for masks (LYMPH)
# Continues on errors; prints summary.

set -euo pipefail

COHORTS=("HC" "MS")

# ---------- helpers ----------
resolve() {  # resolve <basepath-without-ext>  -> echoes .nii.gz or .nii
  [[ -f "$1.nii.gz" ]] && { echo "$1.nii.gz"; return 0; }
  [[ -f "$1.nii"    ]] && { echo "$1.nii";    return 0; }
  return 1
}

cleanup_outputs() {
  local ID="$1"
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
    "${ID}_T1_to_MFFE_ANTS_"* \
    "${ID}_T2AX_to_MFFE_ANTS_"* \
    "${ID}_B0_to_MT_ANTS_"* \
    "${ID}_B1_to_MT_ANTS_"* \
    2>/dev/null || true
  rm -f "${ID}_MTdynANTS_"* "${ID}_MFAdynANTS_"* 2>/dev/null || true
  rm -f "${ID}"_*_ANTS_Warped.nii.gz "${ID}"_*_ANTS_InverseWarped.nii.gz 2>/dev/null || true
  rm -f *.json *_inv.nii.gz 2>/dev/null || true
}

run_one() (
  set -euo pipefail
  shopt -s nullglob
  local subj_dir="$1"
  local ID; ID="$(basename "$subj_dir")"
  echo "=== [$ID] start ==="
  cd "$subj_dir"

  # inputs (accept .nii or .nii.gz)
  local MT="$(resolve "${ID}_MT")"       || { echo "[$ID] missing MT"; exit 1; }
  local MFA="$(resolve "${ID}_MFA")"     || { echo "[$ID] missing MFA"; exit 1; }
  local T1="$(resolve "${ID}_T1")"       || { echo "[$ID] missing T1"; exit 1; }
  local T2AX="$(resolve "${ID}_T2AX")"   || { echo "[$ID] missing T2AX"; exit 1; }
  local B0="$(resolve "${ID}_B0")"       || { echo "[$ID] missing B0"; exit 1; }
  local B1="$(resolve "${ID}_B1")"       || { echo "[$ID] missing B1"; exit 1; }
  local MFFE="$(resolve "${ID}_MFFE")"   || { echo "[$ID] missing MFFE"; exit 1; }
  local LYM_RAW
  if   [[ -f "${ID}_LYMPH.nii.gz" ]]; then LYM_RAW="${ID}_LYMPH.nii.gz"
  elif [[ -f "${ID}_LYMPH.nii"    ]]; then LYM_RAW="${ID}_LYMPH.nii"
  else echo "[$ID] missing LYMPH"; exit 1; fi

  cleanup_outputs "$ID"

  # ---------- orient lymph mask (copy geometry from T2AX) ----------
  fslswapdim "$LYM_RAW" x -y z "${ID}_LYMPH_tmp.nii.gz"
  fslcpgeom  "$T2AX"    "${ID}_LYMPH_tmp.nii.gz"
  mv "${ID}_LYMPH_tmp.nii.gz" "${ID}_LYMPH_oriented_ANTS.nii.gz"

  # ---------- MT intra-dyn (SyN to ensure 1Warp exists) ----------
  local MT_REF_IDX="0008"      # 0-based index string produced by fslsplit
  local MT_SPLIT="${ID}_MTdynANTS_"
  rm -f ${MT_SPLIT}* 2>/dev/null || true
  fslsplit "$MT" "$MT_SPLIT" -t
  local MT_REF_VOL="${ID}_MT_ref_ANTS.nii.gz"
  cp "${MT_SPLIT}${MT_REF_IDX}.nii.gz" "$MT_REF_VOL"

  for idx in $(seq -f "%04g" 0 15); do
    mov="${MT_SPLIT}${idx}.nii.gz"
    out="${MT_SPLIT}${idx}_reg_ANTS.nii.gz"
    if [[ $mov == "${MT_SPLIT}${MT_REF_IDX}.nii.gz" ]]; then
      cp "$mov" "$out"
      continue
    fi
    # SyN (nonlinear) -> produces 1Warp + 0GenericAffine
    antsRegistrationSyNQuick.sh -d 3 -f "$MT_REF_VOL" -m "$mov" \
      -o "${MT_SPLIT}${idx}_ANTS_" -t s -n 4
    antsApplyTransforms -d 3 -i "$mov" -r "$MT_REF_VOL" \
      -o "$out" \
      -t "${MT_SPLIT}${idx}_ANTS_1Warp.nii.gz" \
      -t "${MT_SPLIT}${idx}_ANTS_0GenericAffine.mat" \
      -n Linear
  done
  fslmerge -t "${ID}_registered_MT_ANTS.nii.gz" ${MT_SPLIT}*_reg_ANTS.nii.gz
  rm -f ${MT_SPLIT}* || true

  # ---------- MFA intra-dyn (SyN as well) ----------
  local MFA_SPLIT="${ID}_MFAdynANTS_"
  rm -f ${MFA_SPLIT}* 2>/dev/null || true
  fslsplit "$MFA" "$MFA_SPLIT" -t
  for mov in ${MFA_SPLIT}*.nii.gz; do
    [[ $mov == "${MFA_SPLIT}0000.nii.gz" ]] && continue
    antsRegistrationSyNQuick.sh -d 3 -f "${MFA_SPLIT}0000.nii.gz" -m "$mov" \
      -o "${mov%.nii.gz}_ANTS_" -t s -n 4
    antsApplyTransforms -d 3 -i "$mov" -r "${MFA_SPLIT}0000.nii.gz" \
      -o "${mov%.nii.gz}_reg_ANTS.nii.gz" \
      -t "${mov%.nii.gz}_ANTS_1Warp.nii.gz" \
      -t "${mov%.nii.gz}_ANTS_0GenericAffine.mat" \
      -n Linear
  done
  cp "${MFA_SPLIT}0000.nii.gz" "${MFA_SPLIT}0000_reg_ANTS.nii.gz"
  fslmerge -t "${ID}_registered_MFA_ANTS.nii.gz" ${MFA_SPLIT}*_reg_ANTS.nii.gz
  rm -f ${MFA_SPLIT}* || true

  # ---------- mFFE → MT (SyN) ----------
  antsRegistrationSyNQuick.sh -d 3 -f "$MT_REF_VOL" -m "$MFFE" \
    -o "${ID}_MFFE_to_MT_ANTS_" -t s -n 4
  antsApplyTransforms -d 3 -i "$MFFE" -r "$MT_REF_VOL" \
    -o "${ID}_registered_mFFE_ANTS.nii.gz" \
    -t "${ID}_MFFE_to_MT_ANTS_1Warp.nii.gz" \
    -t "${ID}_MFFE_to_MT_ANTS_0GenericAffine.mat" \
    -n Linear

  # ---------- T1/T2AX → MFFE (SyN) → MT (SyN) ----------
  for MOD in T1 T2AX; do
    SRC="$(eval echo \$$MOD)"
    antsRegistrationSyNQuick.sh -d 3 -f "$MFFE" -m "$SRC" \
      -o "${ID}_${MOD}_to_MFFE_ANTS_" -t s -n 4
    antsApplyTransforms -d 3 -i "$SRC" -r "$MT_REF_VOL" \
      -o "${ID}_registered_${MOD}_ANTS.nii.gz" \
      -t "${ID}_MFFE_to_MT_ANTS_1Warp.nii.gz" \
      -t "${ID}_MFFE_to_MT_ANTS_0GenericAffine.mat" \
      -t "${ID}_${MOD}_to_MFFE_ANTS_1Warp.nii.gz" \
      -t "${ID}_${MOD}_to_MFFE_ANTS_0GenericAffine.mat" \
      -n Linear
  done

  # ---------- B0/B1 → MT (SyN) ----------
  for MOD in B0 B1; do
    SRC="$(eval echo \$$MOD)"
    antsRegistrationSyNQuick.sh -d 3 -f "$MT_REF_VOL" -m "$SRC" \
      -o "${ID}_${MOD}_to_MT_ANTS_" -t s -n 4
    antsApplyTransforms -d 3 -i "$SRC" -r "$MT_REF_VOL" \
      -o "${ID}_registered_${MOD}_ANTS.nii.gz" \
      -t "${ID}_${MOD}_to_MT_ANTS_1Warp.nii.gz" \
      -t "${ID}_${MOD}_to_MT_ANTS_0GenericAffine.mat" \
      -n Linear
  done

  # ---------- LYMPH mask: T2AX→MFFE (SyN) then MFFE→MT (SyN); NN interpolation ----------
  antsApplyTransforms -d 3 \
    -i "${ID}_LYMPH_oriented_ANTS.nii.gz" -r "$MT_REF_VOL" \
    -o "${ID}_registered_LYMPH_ANTS.nii.gz" \
    -t "${ID}_MFFE_to_MT_ANTS_1Warp.nii.gz" \
    -t "${ID}_MFFE_to_MT_ANTS_0GenericAffine.mat" \
    -t "${ID}_T2AX_to_MFFE_ANTS_1Warp.nii.gz" \
    -t "${ID}_T2AX_to_MFFE_ANTS_0GenericAffine.mat" \
    -n NearestNeighbor

  rm -f *.json *_inv.nii.gz 2>/dev/null || true
  echo "=== [$ID] done ==="
)

# ---------- main ----------
if [[ ! -d "HC" || ! -d "MS" ]]; then
  echo "ERROR: Run this from the Processed/ directory (must contain HC/ and MS/)."
  exit 1
fi

declare -a OKS=()
declare -a FAILS=()

shopt -s nullglob

for C in "${COHORTS[@]}"; do
  echo "--- Cohort: $C ---"
  declare -a SUBJS=()
  for d in "$C"/*; do
    [[ -d "$d" ]] || continue
    b="$(basename "$d")"
    [[ "$b" =~ ^[0-9]+$ ]] && SUBJS+=("$d")
  done

  echo "Found ${#SUBJS[@]} numeric subject dir(s) in $C."
  if [[ ${#SUBJS[@]} -eq 0 ]]; then
    echo "  Tip: expected $C/12345/, $C/00001/, etc."
    echo "  Contents sample:"; ls -1 "$C" | head -n 10 || true
    continue
  fi

  for d in "${SUBJS[@]}"; do
    if run_one "$d"; then
      OKS+=("${d}")
    else
      echo "!!! [${d}] FAILED — continuing"
      FAILS+=("${d}")
    fi
  done
done

echo
echo "================ SUMMARY ================"
echo "Succeeded: ${#OKS[@]}"
[[ ${#OKS[@]} -gt 0 ]] && printf '  %s\n' "${OKS[@]#./}"
echo "Failed:    ${#FAILS[@]}"
[[ ${#FAILS[@]} -gt 0 ]] && printf '  %s\n' "${FAILS[@]#./}"
echo "========================================="

exit $(( ${#FAILS[@]} > 0 ? 2 : 0 ))
