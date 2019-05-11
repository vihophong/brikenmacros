for i in {1..9000}
do
	root -b -q  "gensim.C()"
	root -b -q  "mlhfitsim.C()"
done
