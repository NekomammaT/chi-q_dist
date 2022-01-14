rm PcalMqchi.dat
lncalMstep=9
for iM in `seq 0 $lncalMstep`
do
    ./PcalMqchi $iM $lncalMstep
done
