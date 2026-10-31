#!/bin/zsh
# A1 local test (variant L, 2026-10-31): single-call reference, then N=1,2,4 local nodes (workers=1 each),
# then distributed_summary.py with the expected body SHA-256 of the VEP chr1-22 output (repeat 2).
# Local processes share one machine: correctness and overhead only, not distribution.
set -u
D=~/workspace/data_vepyr
SC=${0:A:h}
PY=$D/.venv-native/bin/python3
O=$D/runs/a1_local_20261031
IN=$D/input/HG002_normalized.vcf.gz
$PY $SC/distributed_annotate.py single $IN n1 $O/nodes1.json $O/single > $O/single.out 2>&1 || { echo "single failed" >> $O/chain.log; exit 1; }
echo "single ok" >> $O/chain.log
for n in 1 2 4; do
  rm -rf $O/nodes
  $PY $SC/distributed_annotate.py plan $IN $O/nodes$n.json $O/N$n > $O/N$n.plan.out 2>&1 &&
  $PY $SC/distributed_annotate.py run $O/N$n > $O/N$n.run.out 2>&1
  echo "N$n rc=$?" >> $O/chain.log
  rm -rf $O/N$n/shards $O/N$n/parts
done
rm -rf $O/nodes

$PY -I $SC/distributed_summary.py $O $O/summary.json 5fee6224bf135d80075712141cc06194e0730fa253ff8b131a9d25d661647ee9 > $O/summary.out 2>&1
echo "summary rc=$?" >> $O/chain.log
