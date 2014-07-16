echo "delete all jobs (y/n)?"
read a
if [[ $a == "Y" || $a == "y" ]]; then

#get all the jobids
#bjobs | grep lewang | cut -c1-8 > job2kill
#squeue -u lewang  | cut -c 13-18 > job2kill 

n=`cat job2kill | wc -l`

for ((i=1; i<=$n; i++))
do
jobid=`head -n $i job2kill | tail -n 1`
#echo $jobid
#bkill  $jobid
scancel  $jobid
done 
fi
