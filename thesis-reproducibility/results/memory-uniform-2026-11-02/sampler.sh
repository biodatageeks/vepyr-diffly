
cg=/sys/fs/cgroup/docker/$1
i=0; while [ ! -d "$cg" ] && [ $i -lt 100 ]; do sleep 0.05; i=$((i+1)); done
maxr=0; maxp=0; peak=0; n=0
while [ -d "$cg" ]; do
  set -- $(for p in $(cat "$cg/cgroup.procs" 2>/dev/null); do cat /proc/$p/smaps_rollup 2>/dev/null; done |
           awk '/^Rss:/{r+=$2} /^Pss:/{p+=$2} END{print r+0, p+0}') "$1"
  [ "$1" -gt "$maxr" ] && maxr=$1
  [ "$2" -gt "$maxp" ] && maxp=$2
  v=$(cat "$cg/memory.peak" 2>/dev/null) && peak=$v
  n=$((n+1)); sleep 0.1
done
echo "$maxr $maxp $peak $n"
