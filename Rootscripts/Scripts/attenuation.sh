#!/bin/bash

rscript=.
rfilepath=../Rootfiles/attenuation
filelist=
pdata=data.txt
prefix=tag
suffix=.root
choice=
rfile=./default.root


while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
		-e | --elem )
		shift; choice=$1 ;;
		-r | --rootfile )
		shift; rfile=$1 ;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi
if [[ $choice -eq 0 ]]; then
	echo "Selected Al"
	filelist=('060X' '061X' '062X' '063X' '067X')
	prefix=Al/${prefix}
	pdata=${rfilepath}/Al/${pdata}
elif [[ $choice -eq 1 ]]; then
	echo "Selected Cu"
	filelist=('060X' '064X' '065X' '066X')
	prefix=Cu/${prefix}
	pdata=${rfilepath}/Cu/${pdata}
else
	echo "INVALID INPUT"
	echo "INPUT 0 (Al) or 1 (Cu) with -e"
	exit
fi

files=
for val in "${filelist[@]}"; do
	files=${files}${rfilepath}/${prefix}${val}${suffix}'\n'
done
files=\"${files}\"
pdata=\"${pdata}\"
rfile=\"${rfile}\"
bkgFull=\"../Rootfiles/background/tag027X${suffix}\"

root "${rscript}/attenuation.C(${files}, ${pdata}, ${bkgFull}, ${rfile}, ${choice})"
