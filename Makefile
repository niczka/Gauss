all: sprawozdanie.pdf program

sprawozdanie.pdf: sprawozdanie.tex
	texi2pdf sprawozdanie.tex
program: program.cpp
	g++ -o program program.cpp
