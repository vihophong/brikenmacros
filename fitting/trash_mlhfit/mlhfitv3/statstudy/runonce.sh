for run in `seq 1 $2`
do
  /home/phong/briken17/brikenmacrosupdate/brikenmacros/simulation/v2/simulatedecay.sh /home/phong/briken17/brikenmacrosupdate/brikenmacros/simulation/v2/decayparms/parmsIn134_exp.txt /home/phong/briken17/brikenmacrosupdate/brikenmacros/simulation/v2/simparms/simparmsCd130_highin_x0.05.txt temp/outhist$1.root
  ../offlinefit.sh fit$1 temp/outhist$1.root parmsIn134var.txt temp/temp$1.root
done
