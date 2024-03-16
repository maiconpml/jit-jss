filenames=`ls ./$1`
echo path lowerb makes penalties earlPenal tarPenal
for eachfile in $filenames
do
	./build/main $1/$eachfile lala --searchType tabuSearch --maxSecs 10
done
