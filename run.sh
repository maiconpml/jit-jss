yourfilenames=`ls ./$1`
echo path lowerb makes penalties earlPenal tarPenalties
for eachfile in $yourfilenames 
do
#echo $eachfile	
./build/main $1/$eachfile lalala --searchType tabuSearch
done

