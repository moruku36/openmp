#set terminal jpeg  enhanced font "Times" 20
#set tics font 'Times,18'
#set rmargin 3
#set lmargin 6
set nokey

if (exist("n")==0 || n<1) n=1

file0(n) = sprintf("%d.dat",n) 
title(n) = sprintf("time step = %d",n) 

unset label 
set label title(n) at 0.5 , 1.005

plot [0:1][0.99:1.01]  file0(n) with points pointsize 4 pt 7
replot

pause 0.05

if (n < 300)  n=n+1; reread
