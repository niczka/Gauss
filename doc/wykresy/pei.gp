set terminal postscript enhanced 
set output 'pei_plot.eps'
set key inside right bottom vertical Right noreverse enhanced 
set samples 10000
set xlabel "rozmiar macierzy"
set ylabel "||b-Ax||"
set grid
set title "Eliminacja Gaussa - macierz Pei"
set offsets 0
set multiplot
plot [0:100][0:0.000000000001] 'pei_without' ti "Bez wyboru" with lines,\
'pei_from_row' ti "Wybor z wiersza" with lines,\
'pei_from_col' ti "Wybor z kolumny" with lines,\
'pei_full' ti "Pelny wybor" with lines
