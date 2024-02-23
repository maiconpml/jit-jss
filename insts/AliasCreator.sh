
#!/usr/bin/env bash


##Copy inside all folder - enable exec

##fromFolder=../JSP/
##fromFolder=../OSP/
fromFolder=../DemirkolJsp/
instances=(`ls -1 ${fromFolder}`)


for inst in "${instances[@]}"
do
##ln -s ../TestFrom/test.txt ./
    ln -s "${fromFolder}${inst}" ./
done


