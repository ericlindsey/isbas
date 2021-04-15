#!/bin/bash

wavelen=0.236057

#create a gaussian shape
gmt grdmath -r -R103.5/104/1/1.5 -I0.01 X 103.75 SUB SQR Y 1.25 SUB SQR ADD NEG 100 MUL EXP  = temp.grd

#add an unwrapping error to one interferogram
gmt grdmath -Rtemp.grd X 103.6 LT X 103.55 GT MUL Y 1.2 GT Y 1.25 LT MUL MUL PI MUL = box.grd

#create interferograms with the gaussian shape uplifting at a constant rate, and convert from cm/yr to radians/day
gmt grdmath 539 355 SUB 365.25 DIV 4 MUL PI MUL 100 DIV NEG $wavelen DIV temp.grd MUL = intf/2006355_2007174/unwrap_mask_ll_pad.grd
gmt grdmath 585 355 SUB 365.25 DIV 4 MUL PI MUL 100 DIV NEG $wavelen DIV temp.grd MUL 17 ADD = intf/2006355_2007220/unwrap_mask_ll_pad.grd
gmt grdmath 631 355 SUB 365.25 DIV 4 MUL PI MUL 100 DIV NEG $wavelen DIV temp.grd MUL = intf/2006355_2007266/unwrap_mask_ll_pad.grd
gmt grdmath 585 539 SUB 365.25 DIV 4 MUL PI MUL 100 DIV NEG $wavelen DIV temp.grd MUL = intf/2007174_2007220/unwrap_mask_ll_pad.grd
gmt grdmath 631 539 SUB 365.25 DIV 4 MUL PI MUL 100 DIV NEG $wavelen DIV temp.grd MUL box.grd ADD = intf/2007174_2007266/unwrap_mask_ll_pad.grd
gmt grdmath 631 585 SUB 365.25 DIV 4 MUL PI MUL 100 DIV NEG $wavelen DIV temp.grd MUL = intf/2007220_2007266/unwrap_mask_ll_pad.grd

# clean up
rm -f box.grd temp.grd

