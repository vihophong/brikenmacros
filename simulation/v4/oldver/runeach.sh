for i in {1..1}
do
	root -b -q  "simabinitio.C()"
	root -b -q  "checksimulation.C()"
done
