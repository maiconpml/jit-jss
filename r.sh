filenames=`ls ./$1`
echo path lowerb makes penalties earPenal tarPenal
for eachfile in $filenames
do
	./$2/main $1/$eachfile lala --maxSecs 30
done
