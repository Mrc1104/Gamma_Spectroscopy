#! /bin/bash
# example
# ./Scripts/makerootfiles.sh -d ../data/energy_calibration -r ../Rootfiles/energy_calibration
nargs=$#
args="$@"

runtype=('ec' 'ang' 'bkg' 'sc' 'hwc')
rundescr=('energy_calibration' 'angular' 'background' 'Scope' 'HW_coincidence')


declare rfile
declare rscript=.
declare dfile

while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
		-d | --data )
		shift; dfile=$1 ;;
		-r | --rootfile )
		shift; rfile=$1 ;;
		-R | --rootscript )
		shift; rscript=$1 ;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi


for file in ${dfile}/tag*0.xy; do
	path=\"${file%/tag*}\"
	TAG=${file#*tag}
	TAG=${TAG%.*}
	savefile=\"${rfile}/"tag"${TAG%?}"X.root"\"
	TAG=\"${TAG}\"
	root -q "${rscript}/make_rootfiles.C(${path}, ${TAG}, ${savefile})"
done


