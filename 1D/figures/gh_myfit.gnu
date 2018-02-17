set term post eps enh 25 col solid
set out 'gh_myfit.eps'

g(x)=((x/a)**q/(1+(x/b)**p))*exp(-x/h)

a=5
b=6
q=4
p=6
h=5

set log y
set yr [1e-5:1]

set ytics ("10^{-5}" 1e-5,"10^{-4}" 1e-4,"10^{-3}" 1e-3,"10^{-2}" 1e-2,"10^{-1}" 1e-1,"1" 1)

set xlabel 'h'
set ylabel 'g(h)'

p 'pdf_ti300000-tf300000.dat' w lp pt 7 ps 2 t '{/=18 sim of our model}',g(x) w l lt -1 lw 3 t '{/=18 my fit}'

reset
set term x11
