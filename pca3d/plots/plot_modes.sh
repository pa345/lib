#!/bin/sh

mode_min="1"
mode_max="20"

for m in $(seq -w ${mode_min} ${mode_max}); do
  echo "Plotting mode ${m}..."
  gnuplot -c ./plot_modes.gp ${m}
done
