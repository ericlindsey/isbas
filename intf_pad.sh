#!/bin/bash

if [[ $# -lt 1 ]]; then
  fname="unwrap_mask_ll.grd"
else 
  if [[ "$1" == "-h" ]]; then
    echo "usage: $0 [grdfilename] [bounds_file]"
    echo "computes the exterior bounds matching all files intf/*/[name] and pads them all to match"
    echo "optional second argument gets the bounds from that grid file instead of computing the maximum"
    echo "default: unwrap_mask_ll.grd"
    exit 1
  fi
  fname=$1
fi
echo "Begin intf_pad.sh"

padfile=$(basename $fname .grd)_pad.grd
echo "reading files intf/*/$fname"
echo "writing files intf/*/$padfile"

if [[ $# -eq 2 ]]; then
  bounds=`gmt grdinfo -I- $2`
  echo using input bounds: $bounds
else
  bounds=`gmt grdinfo -C intf/*/$fname |cut -f 2-5 |gmt gmtinfo -C | cut -f 1,4,5,8 |awk '{printf("-R%.12f/%.12f/%.12f/%.12f",$1,$2,$3,$4)}'`
  echo "computed maximum exterior bounds for all files: $bounds"
fi

echo "renaming old intf.in, and creating a new one from logs_intf/"
time=`date +%Y_%m_%d-%H_%M_%S`
mv intf.in intf.in.bak.$time

echo "padding files..."
for dir in  `ls -d intf/*`
do
  if [[ -f $dir/$fname ]]; then
    echo padding $dir/$fname
    #gmt grdcut $dir/$fname $bounds -Gtemp.grd
    incr=`gmt grdinfo -I $dir/$fname`
    gmt grd2xyz $dir/$fname -bo | gmt xyz2grd -bi $bounds $incr -r -Gtemp_$fname
    nccopy -k classic temp_$fname $dir/$padfile
    cat logs_$dir.in >> intf.in
  else
    echo "$dir/$fname not found."
  fi
done
rm temp_$fname

echo "End intf_pad.sh"
