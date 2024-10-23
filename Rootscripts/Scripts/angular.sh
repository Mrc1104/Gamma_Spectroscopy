# !root 'angular_study.C("../Rootfiles/angular/tag040X.root", "0", "./angstudy.root")'
#!/bin/bash


declare rfile
declare rscript=.
declare dfile

angles=(0,5,10,15,20,30,-5,-10,-15,0,-20)
duration=(10,10,10,10,10,10,10,10,10,5,5);

while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
		-d | --data )
		shift; dfile=$1 ;;
		-r | --rootfile )
		shift; rfile=$1 ;;
		-R | --rootscript )
		shift; rscript=$1 ;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi


allfile=
for file in ${dfile}/tag*X.root; do
	echo ${file}	
	allfile=${allfile}${file}'\n'
done
allangle=
for val in "${angles[@]}"; do
	allangle=${allangle}${val}
done
echo ${allangle}
allduration=
for val in "${duration}"; do
	allduration=${allduration}${val}
done
echo ${allduration}

allangle=\"${allangle}\"
allfile=\"${allfile}\"
allduration=\"${allduration}\"
rfile=\"${rfile}\"

root "${rscript}/angular_study.C(${allfile}, ${allangle}, ${allduration}, ${rfile})"

