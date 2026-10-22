#!/bin/bash
# z3xz6_sweep -- count `W = ` lines against distinct models in one raw output file.
# A count of ~2x the model number is the signature of the forked-child duplicate print.
f=$1
echo "W-lines:            $(grep -c '^  W = ' "$f")"
echo "distinct couplings prompts: $(grep -o '/MSSM[0-9]*/couplings>' "$f" | sort -u | wc -l)"
echo "distinct spectrum prompts:  $(grep -o '/MSSM[0-9]*/spectrum>' "$f" | sort -u | wc -l)"
echo "waiting done:       $(grep -c 'waiting done' "$f")"
echo "models in list:     $(grep -c '^Label:' "$(dirname "$f")/models.txt")"
