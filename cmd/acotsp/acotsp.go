// solves a TSP with ACO as described in Dorigo et al. 96
package main

import (
	"flag"
	"fmt"
	"os"

	"bitbucket.org/baobabsoluciones/aco"
	"bitbucket.org/baobabsoluciones/aco/tsplib"
)

var (
	tsp = flag.String("tsp", "", "TSP problem file in TSPLIB format")
)

func main() {
	flag.Parse()
	if *tsp == "" {
		fmt.Println("fatal error: require flag -tsp")
		flag.Usage()
		os.Exit(1)
	} else {
		// TODO
		g, err := tsplib.ParseTSPLIBProblem(*tsp)
		if err != nil {
			// TODO
			panic(err)
		}
		s, stagBehv, err := aco.ASBestParams(g)
		if err != nil {
			// TODO
			panic(err)
		}
		fmt.Printf("solution: %v", s)
		if stagBehv {
			fmt.Println("AS terminated with stagnation behaviour")
		}
	}
}
