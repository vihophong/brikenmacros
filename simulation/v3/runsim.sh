for i in {1..1}
do
	root -b -q  "simabinitio.C()"
	root -b -q  "checksimulation.C()"
	#./runfit.sh outhist.root parmsfitex.txt 0.62 0.062 0 0 500 -10 10
done