cd /home/user/3body-qc-v4
for t in A B; do for cut in 0.335 0.345; do for rep in 1 2; do
 fix=""; [ $t = B ] && fix="--fix KB KE"
 python3 -m qc3.fit --cutoffs $cut --form q0_consistent $fix --no-verify --out $1/${t}_${cut}_$rep.json > $1/${t}_${cut}_$rep.log 2>&1
done; done; done
echo done > $1/DONE
