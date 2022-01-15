rm PcalMqchi.dat
lncalMstep=92
for iM in `seq 0 $lncalMstep`
do
    ./PcalMqchi $iM $lncalMstep
done
