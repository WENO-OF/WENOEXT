#!/bin/bash

# ---   System Commands   ---
[ -d TEX ] || mkdir TEX
[ -d EPS ] || mkdir EPS
# --- END System Commands ---


gnuplot <<-EOFMarker
    # Select if a epslatex or an eps file
    # is created.
    DIR = "TEX";


    if (DIR eq "TEX") \
        ENDING = "tex"; \
    else \
        ENDING = "eps";

    if (DIR eq "TEX") \
        set term epslatex; \
    else \
        set term postscript eps enhanced color font 'Helvetica,16'

    set style line 1 lt 1 lc rgb "black"   lw 5  ps 2 pt 4 
    set style line 2 lt 1 lc rgb "red"     lw 5  ps 2 pt 5
    set style line 3 lt 1 lc rgb "blue"    lw 5  ps 2 pt 6
    set style line 5 lt 2 lc rgb "black"   lw 6  dt 2
    set style line 6 lt 2 lc rgb "red"     lw 6  dt 2
    set style line 7 lt 4 lc rgb "black"   lw 6  dt 2
    set style line 8 lt 2 lc rgb "blue"    lw 6  dt 2


    set xtics font ", 18"
    set ytics font ", 18"

    set grid

    # -------------------------------------------------------------------
    #                           Plots
    # -------------------------------------------------------------------

    set logscale y

    set output DIR."/divergenceSchemeAccuracy.".ENDING
    #set title 'Volume Fraction'
    if (DIR eq "TEX") \
        set format y '\huge %g'; \
        set format x '\huge %g'; \
        set xlabel '\Huge Mesh Size';\
        set ylabel '\Huge MeanError';\
    else \
        set xlabel 'Mesh Size'; \
        set ylabel 'MeanError';

    plot "data.dat"  using 1:3  title '\large limitedLinear' ls 3 with linespoint,\
         "data.dat"  using 1:2  title '\large linear' ls 1 with linespoint,\
         "data.dat"  using 1:4  title '\large WENO' ls 2 with linespoint

EOFMarker
