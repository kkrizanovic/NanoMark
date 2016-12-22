#! /bin/sh

./run-grpaper-ovlwin.sh
./run-dnadiff-generic.sh MiniasmRacon-grpaper-ovlwin MiniasmRacon-grpaper
./collect-grpaper-generic.sh MiniasmRacon-grpaper-ovlwin

./run-grpaper-normal.sh
./run-dnadiff-generic.sh MiniasmRacon-grpaper-normal MiniasmRacon-grpaper
./collect-grpaper-generic.sh MiniasmRacon-grpaper-normal
