for i in {1..9000}
do
	root -b -q  "simabinitio.C()"
	root -b -q  "checksimulation.C()"
	root -b -q  "mlhfitsim.C()"	
done
