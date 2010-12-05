set terminal postscript enhanced 
set output 'hilbert_plot.eps'
set key inside right bottom vertical Right noreverse enhanced 
set samples 10000
set xlabel "rozmiar macierzy"
set ylabel "||b-Ax||"
set grid
set title "Eliminacja Gaussa - macierz Hilberta"
set offsets 0
set multiplot
plot [0:100][0:100000] 'hilbert_without' ti "Bez wyboru" with lines,\
'hilbert_from_row' ti "Wybor z wiersza" with lines,\
'hilbert_from_col' ti "Wybor z kolumny" with lines,\
'hilbert_full' ti "Pelny wybor" with lines
