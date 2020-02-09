all: compile run 
compile: 
	gfortran -fdefault-real-8  lib/*.f90 main.f90 -o md
run:
	./md