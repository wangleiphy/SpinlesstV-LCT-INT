#python /mnt/lnec/lewang/spinlessctint/analysis/remaining_time.py -f `ls -tr ../jobs/*.log | tail -n $1`


for f in `ls -tr ../jobs/*.log | tail -n $1`
do

echo " "
echo "<<<<<<<<<<<<" $f ">>>>>>>>>>>>>"
head -n5 $f
grep Completed  $f | tail -n 5 
#tail -n5 $f 

echo " "

done
