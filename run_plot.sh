#!/bin/bash

# Base name of the data file
BASE_NAME="BWx022t_0.200308"
base_mesh="8x8_lv5"
slimiter="WENO5"
rsolver="RUSA"

datafile="/data/${BASE_NAME}.dat"

echo "Data file to be used: ${datafile}"

# Run Gnuplot with the script
gnuplot << EOF
# Construct the full data file name in Gnuplot
datafile="./data/${BASE_NAME}.dat"

# Common settings
set palette maxcolors 14
set style fill solid 1 border lc black
set size ratio -1
set size 1.0, 1.0

# Density plot
set term png size 8000, 8000
set output "./plots/${BASE_NAME}_${base_mesh}_${slimiter}_${rsolver}_density.png"
plot datafile using 1:2:3:4:5 with boxxyerror fc palette title "Density"

# U plot
set term png size 8000, 8000
set output "./plots/${BASE_NAME}_${base_mesh}_${slimiter}_${rsolver}_U.png"
plot datafile using 1:2:3:4:6 with boxxyerror fc palette title "U"

# V plot
set term png size 8000, 8000
set output "./plots/${BASE_NAME}_${base_mesh}_${slimiter}_${rsolver}_V.png"
plot datafile using 1:2:3:4:7 with boxxyerror fc palette title "V"

# Pressure plot
set term png size 8000, 8000
set output "./plots/${BASE_NAME}_${base_mesh}_${slimiter}_${rsolver}_pressure.png"
plot datafile using 1:2:3:4:9 with boxxyerror fc palette title "Pressure"

# Bx plot
set term png size 8000, 8000
set output "./plots/${BASE_NAME}_${base_mesh}_${slimiter}_${rsolver}_Bx.png"
plot datafile using 1:2:3:4:10 with boxxyerror fc palette title "Bx"

# Pressure plot
set term png size 8000, 8000
set output "./plots/${BASE_NAME}_${base_mesh}_${slimiter}_${rsolver}_By.png"
plot datafile using 1:2:3:4:11 with boxxyerror fc palette title "By"

# Add more plots as needed for other variables
EOF

echo "Plotting complete."

