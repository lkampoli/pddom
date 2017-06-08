################################################################################
# GNUPLOT: Reflectance matrix
# Plot by J.P. Clark using template from:
#
#
# Motivation and original template:
#   http://www.sciencetronics.com/greenphotons/?p=570
#
# Explanation in future version.
#
# Modified by J.P. Dark
# August 11, 2016
# Tested with gnuplot version 5.0 patch 0
################################################################################

# User defined parameters.
nrows = 4 # Number of rows.
ncol = 4 # Number of columns.

# Outer margins (measurements in inches).
mtop = 0.6 # Outer top margin, title goes here
mbottom = 1 # Outer bottom margin, x label goes here
mleft = 0.6 # Outer left margin, y label goes here
mright = 0.6 # Outer right margin, y2 label goes here

# Define inner margins (measured in inches).
mdx = 1 # Inner horizontal spacing between plots.
mdy = 0.2 # Inner vertical spacing between plots.

# Define plot size (measured in inches).
plt_height = 1.5 # Inch height of individual plots
plt_width = 2.0 # Inch width of individual plots

# Calculate full dimensions given margins.
xsize = mleft + mright + (plt_width * ncol) + (ncol - 1) * mdx
ysize = mtop + mbottom + (nrows*plt_height) + (nrows-1) * mdy


################################################################################
# MACRO: Set margins.
#
# Modified from the following source:
#   http://www.sciencetronics.com/greenphotons/?p=570
#
bot(n) = (mbottom + (nrows - n) * plt_height + (nrows - n) * mdy) / ysize
top(n) = 1 - ((mtop + (n-1) * (plt_height + mdy)) / ysize)
left(n) = (mleft + (n-1) * plt_width + (n-1) * mdx) / xsize
right(n) = 1 - ((mright + (ncol-n) * (plt_width+mdx)) / xsize)

first_subplot = "row=1; col=1; \
    set tmargin at screen top(row); \
    set bmargin at screen bot(row); \
    set lmargin at screen left(col); \
    set rmargin at screen right(col); "

next_subplot = " \
    if (col < ncol) {col = col + 1;} else {col = 1; row = row+1;} \
    set tmargin at screen top(row); \
    set bmargin at screen bot(row); \
    set lmargin at screen left(col); \
    set rmargin at screen right(col)"
################################################################################

# # Begin plot
mode = ARG1
fileout = 'venus'.mode.'.pdf'
# fileout = 'venus'.mode.'tex'
file1 = 'venus_refl_'.mode.'.dat'
file2 = ARG2.'/venus_refl_'.mode.'.txt'

#set terminal postscript eps enhanced color dl 2.0 size xsize, ysize 28
set terminal pdf size xsize, ysize
set encoding iso_8859_1
set tics scale 1.5

set output fileout


# Set some properties for the plot
set offsets
set autoscale fix
set size 1,1
set nokey


# Format x-axis: formatting is the same for all subplots.
set xrange [0:1]
set mxtics 3
set xlabel '$\mu$';

# Format y-axis: some formatting is the same for all subplots.
# set format y "%-2.1f"
set ylabel '';
set format y '%G';
set mytics 2

unset grid

set multiplot

#  R11
@first_subplot
set xlabel ''; set format x ''
#set yrange [-1.2:1.2]
plot file1 u 1:2 title '' with lines lt 3 lw 2 lc 6, \
     file2 u 1:2 title '' with points lt 3 lw 2 lc 6


# R12
@next_subplot
set xlabel ''; set format x ''
plot file1 u 1:3 title '' with lines lt 3 lw 2 lc 6, \
     file2 u 1:3 title '' with points lt 3 lw 2 lc 6

# R13
@next_subplot
set xlabel ''; set format x ''
plot file1 u 1:4 title '' with lines lt 3 lw 2 lc 6, \
     file2 u 1:4 title '' with points lt 3 lw 2 lc 6

# R14
@next_subplot
set xlabel ''; set format x ''
plot file1 u 1:5 title '' with lines lt 3 lw 2 lc 6, \
     file2 u 1:5 title '' with points lt 3 lw 2 lc 6


# ROW 2
@next_subplot
set title ''
set xlabel ''; set format x ''
plot file1 u 1:6 title '' with lines lt 3 lw 2 lc 6, \
     file2 u 1:6 title '' with points lt 3 lw 2 lc 6

@next_subplot
set title ''
set xlabel ''; set format x ''
plot file1 u 1:7 title '' with lines lt 3 lw 2 lc 6, \
     file2 u 1:7 title '' with points lt 3 lw 2 lc 6

@next_subplot
set title ''
set xlabel ''; set format x ''
plot file1 u 1:8 title '' with lines lt 3 lw 2 lc 6, \
     file2 u 1:8 title '' with points lt 3 lw 2 lc 6

@next_subplot
set title ''
set xlabel ''; set format x ''
plot file1 u 1:9 title '' with lines lt 3 lw 2 lc 6, \
     file2 u 1:9 title '' with points lt 3 lw 2 lc 6

# ROW 3
@next_subplot
set title ''
set xlabel ''; set format x ''
plot file1 u 1:10 title '' with lines lt 3 lw 2 lc 6, \
     file2 u 1:10 title '' with points lt 3 lw 2 lc 6

@next_subplot
set title ''
set xlabel ''; set format x ''
plot file1 u 1:11 title '' with lines lt 3 lw 2 lc 6, \
     file2 u 1:11 title '' with points lt 3 lw 2 lc 6

@next_subplot
set title ''
set xlabel ''; set format x ''
plot file1 u 1:12 title '' with lines lt 3 lw 2 lc 6, \
     file2 u 1:12 title '' with points lt 3 lw 2 lc 6

@next_subplot
set title ''
set xlabel ''; set format x ''
plot file1 u 1:13 title '' with lines lt 3 lw 2 lc 6, \
     file2 u 1:13 title '' with points lt 3 lw 2 lc 6

# ROW 4
@next_subplot
set title ''
set xlabel '$\mu$'; set format x '%g'
plot file1 u 1:14 title '' with lines lt 3 lw 2 lc 6, \
     file2 u 1:14 title '' with points lt 3 lw 2 lc 6

@next_subplot
set title ''
set xlabel '$\mu$'; set format x '%g'
plot file1 u 1:15 title '' with lines lt 3 lw 2 lc 6, \
     file2 u 1:15 title '' with points lt 3 lw 2 lc 6

@next_subplot
set title ''
set xlabel '$\mu$'; set format x '%g'
plot file1 u 1:16 title '' with lines lt 3 lw 2 lc 6, \
     file2 u 1:16 title '' with points lt 3 lw 2 lc 6

@next_subplot
set title ''
set xlabel '$\mu$'; set format x '%g'
plot file1 u 1:17 title '' with lines lt 3 lw 2 lc 6, \
     file2 u 1:17 title '' with points lt 3 lw 2 lc 6




unset multiplot
