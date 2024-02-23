cdyourfilenames=`ls ./$1`
echo path lowerB makes penaltys earlPenaltys tarPenaltys  
for eachfile in $yourfilenames
do
    #echo $eachfile
    ./build/main $1$eachfile lalala --searchType tabuSearch
done