set terminal postscript enhanced 
set output 'dominating_plot.eps'
set key inside right bottom vertical Right noreverse enhanced 
set samples 10000
set ylabel "||b-Ax||"
set xlabel "rozmiar macierzy"
set grid
set title "Eliminacja Gaussa - macierz z dominujaca przekatna"
set offsets 0
set multiplot
plot [0:100][0:0.000000000001] 'dominating_without' ti "Bez wyboru" with lines,\
'dominating_from_row' ti "Wybor z wiersza" with lines,\
'dominating_from_col' ti "Wybor z kolumny" with lines,\
'dominating_full' ti "Pelny wybor" with lines
