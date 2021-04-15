#!/bin/bash

# collect and compress the necessary files for ISBAS, in case you want to run them on a different computer.
# Eric Lindsey, April 2021
#

# first, check that the padding has been done
numfiles=`ls intf/*/unwrap_mask_ll_pad.grd | wc -l`
if [[ $numfiles -eq 0 ]]; then
  echo "Error: no files named unwrap_mask_ll_pad.grd found. Ensure you are in the right folder and have run 'intf_pad.sh' first."
  exit 1
fi
echo Creating minimum-size tar.gz file named gmtsar_output.zip containing $numfiles geocoded interferograms and metadata.
tar -cvzf gmtsar_output.zip topo/master.PRM raw/baseline_table.dat raw/data.in intf.in intf/*/unwrap_mask_ll_pad.grd 

