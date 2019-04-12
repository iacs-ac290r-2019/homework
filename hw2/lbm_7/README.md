Compile
make

Problem 1.1
./lbm_pbc

Problem 1.2
./lbm_mbc


plot.py is deprecated, use viz.ipynb to plot

Input parameter
	input()

Initialization
	init_hydro()
	equilibrium()
	init_pop()

Main
	bc() Receiver point of view, hence to bc first
	stream()
	hydro()
	equilibrium() CORRECT
	collide()

	force()
	obstacle()