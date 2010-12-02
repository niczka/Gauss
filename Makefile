all: doc prog

doc:
	texi2pdf sprawozdanie.tex
prog:
	g++ program.cpp
