#!/bin/bash
#
# Script to plot benchmark results using gnuplot
#
# J.P. Clark
# August 14, 2016
# GNUPLOT V5.0 patch 0
########################


# Go to directory
cd results/venus

# Generate plots of the Venutian results (mu0 = 0.5).
echo '> Venus: Fourier mode k=0'
gnuplot -c venus.gp '0' 'mishchenko_0pt5'
paste -d' ' venus_refl_0.dat \
      mishchenko_0pt5/venus_refl_0.txt > err_0pt5/venus_err_0.txt
gnuplot -c venus_error.gp '0' 'err_0pt5'

echo '> Venus: Fourier mode k=1'
gnuplot -c venus.gp '1' 'mishchenko_0pt5'
paste -d' ' venus_refl_1.dat \
      mishchenko_0pt5/venus_refl_1.txt > err_0pt5/venus_err_1.txt
gnuplot -c venus_error.gp '1' 'err_0pt5'


echo '> Venus: Fourier mode k=2'
gnuplot -c venus.gp '2' 'mishchenko_0pt5'
paste -d' ' venus_refl_2.dat \
      mishchenko_0pt5/venus_refl_2.txt > err_0pt5/venus_err_2.txt
gnuplot -c venus_error.gp '2' 'err_0pt5'


echo '> Venus: Fourier mode k=3'
gnuplot -c venus.gp '3' 'mishchenko_0pt5'
paste -d' ' venus_refl_3.dat \
      mishchenko_0pt5/venus_refl_3.txt > err_0pt5/venus_err_3.txt
gnuplot -c venus_error.gp '3' 'err_0pt5'
