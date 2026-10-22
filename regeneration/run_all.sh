#!/usr/bin/env bash
# Re-run every regeneration script in this package, serially, from the deposited raw data.
#
# Serial execution is deliberate: the scripts build multiprocessing pools, and on Windows
# two overlapping pools can collide on a concurrent open of a geometry file
# (see REPRODUCTION.md, "One transient failure"). Running one at a time avoids it.
#
# Usage:  bash regeneration/run_all.sh          (from the package root)
# Writes: regeneration/logs2/<name>.log  and  regeneration/run_all_status.txt
set -u
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LOGS="$ROOT/regeneration/logs2"
STATUS="$ROOT/regeneration/run_all_status.txt"
mkdir -p "$LOGS"
: > "$STATUS"

run () {  # run <dir> <script> [arg]
  local dir="$1" script="$2" arg="${3:-}"
  local name="${script%.py}${arg:+_$arg}"
  local t0=$(date +%s)
  ( cd "$ROOT/multi_geometry_campaign/code/$dir" && python "$script" $arg ) \
      > "$LOGS/$name.log" 2>&1
  local rc=$?
  local t1=$(date +%s)
  printf '%-42s exit=%d  %ds\n' "$dir/$script $arg" "$rc" "$((t1-t0))" | tee -a "$STATUS"
}

run su5_ten_bound        verify.py            Z2xZ4_1_1
run su5_ten_bound        verify.py            Z3xZ6_1_1
run su5_ten_bound        pattern_census.py    Z2xZ4_1_1
run su5_ten_bound        pattern_census.py    Z3xZ6_1_1
run gate_su5_bound       d1_mechanism.py      Z2xZ4_1_1
run gate_su5_bound       d1_mechanism.py      Z3xZ6_1_1
run gate_su5_bound       d2_pshdup.py         Z2xZ4_1_1
run gate_su5_bound       d2_pshdup.py         Z3xZ6_1_1
run gate_su5_bound       d2_osc10.py
run gate_su5_bound       d5_excluded.py
run gate_su5_bound       d5_keycount.py
run gate_su5_bound       d6_pattern_recount.py Z2xZ4_1_1
run gate_su5_bound       d6_pattern_recount.py Z3xZ6_1_1
run gate_su5_bound       d8_alltwists.py
run gate_su5_bound       d8_classwide.py
run gate_su5_bound       d9_comarks.py
run gate_su5_construction g2_classes.py
run gate_su5_construction g4_lift.py
run gate_su5_construction g8c_doubled_smvalid.py
run su5_construction     s0_invariants.py
run su5_construction     s2_sizing.py

echo "--- two-class re-enumeration (referee standard) ---" | tee -a "$STATUS"
t0=$(date +%s)
( cd "$ROOT/two_classes" && python paperA_two_classes_enum.py ) \
    > "$LOGS/paperA_two_classes_enum.log" 2>&1
rc=$?; t1=$(date +%s)
printf '%-42s exit=%d  %ds\n' "two_classes/paperA_two_classes_enum.py" "$rc" "$((t1-t0))" | tee -a "$STATUS"
t0=$(date +%s)
( cd "$ROOT/two_classes" && python check_bfs_certificate.py ) \
    > "$LOGS/check_bfs_certificate.log" 2>&1
rc=$?; t1=$(date +%s)
printf '%-42s exit=%d  %ds\n' "two_classes/check_bfs_certificate.py" "$rc" "$((t1-t0))" | tee -a "$STATUS"

echo "--- root_lift (signature, scoped) ---" | tee -a "$STATUS"
for s in rl_a_centralizer.py rl_b_data.py rl_c_samefp.py rl_d_projection.py rl_e_summary.py; do
  t0=$(date +%s)
  ( cd "$ROOT/multi_geometry_campaign/code/root_lift" && python "$s" ) \
      > "$LOGS/${s%.py}.log" 2>&1
  rc=$?; t1=$(date +%s)
  printf '%-42s exit=%d  %ds\n' "root_lift/$s" "$rc" "$((t1-t0))" | tee -a "$STATUS"
done

echo "--- uniform-convention census (Table 1 of the paper) ---" | tee -a "$STATUS"
run uniform_census       uniform_census_e6.py Z2xZ4_1_1
run uniform_census       uniform_census_e6.py Z3xZ6_1_1
run uniform_census       uniform_totals.py
run uniform_census       rl_b_uniform.py     Z2xZ4_1_1
run uniform_census       rl_b_uniform.py     Z3xZ6_1_1
run uniform_census       sig_premise_counts.py

echo "--- lattice representatives (v2.2; analysis of the deposited orbifolder output) ---" | tee -a "$STATUS"
run representatives      b4_triality_nets.py
run representatives      b3_compare.py
run representatives      b11b_compare.py
run representatives      b10_brother_analysis.py
run representatives      b8_allbins.py
run representatives      b7_census_rank.py   Z2xZ4_1_1
run representatives      b7_census_rank.py   Z3xZ6_1_1

echo "--- comparison against results_deposited/ ---" | tee -a "$STATUS"
( cd "$ROOT/regeneration" && python compare_to_deposited.py ) > "$LOGS/compare.log" 2>&1
printf '%-42s exit=%d\n' "regeneration/compare_to_deposited.py" "$?" | tee -a "$STATUS"

echo "--- claim map ---" | tee -a "$STATUS"
( cd "$ROOT" && python regeneration/claims_check.py ) > "$LOGS/claims_check.log" 2>&1
printf '%-42s exit=%d\n' "regeneration/claims_check.py" "$?" | tee -a "$STATUS"

echo "DONE" | tee -a "$STATUS"
