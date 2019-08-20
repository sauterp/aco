for i in {1..10}; do
	printf "\nrun %s\n" $i
	go run acotsp.go -tsp ../../testdata/benchmarks/Oliver30.tsp -log log$i.csv -seed $i;
done