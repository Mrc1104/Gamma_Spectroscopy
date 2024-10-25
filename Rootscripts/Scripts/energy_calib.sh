#!/bin/bash

rscript=.
rfilepath=../Rootfiles/energy_calibration
prefix=tag
suffix=.root
NA=('022X' '023X')
NADUR=(5 5)
CS=('025X' '026X')
CSDUR=(5 5)

fileNa=
for val in "${NA[@]}"; do
	fileNa=${fileNa}${rfilepath}/${prefix}${val}${suffix}'\n'
done
fileCs=
for val in "${CS[@]}"; do
	fileCs=${fileCs}${rfilepath}/${prefix}${val}${suffix}'\n'
done
NADURATION=
for val in "${NADUR[@]}"; do
	NADURATION=${NADURATION}${val},
done
CSDURATION=
for val in "${CSDUR[@]}"; do
	CSDURATION=${CSDURATION}${val},
done

NaFull=\"${fileNa}\"
CsFull=\"${fileCs}\"
savefile=\"${rfilepath}/energy_calib${suffix}\"
echo ${NaFull}
echo ${CsFull}
echo ${savefile}

root "${rscript}/energy_calibration.C(${NaFull},\"${NADURATION}\", ${CsFull}, \"${CSDURATION}\", ${savefile})"
