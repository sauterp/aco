// solves a TSP with ACO as described in Dorigo et al. 96
package main

import (
	"flag"
	"fmt"
	"io/ioutil"
	"io"
	"os"
	"time"

	"bitbucket.org/baobabsoluciones/aco"
	"bitbucket.org/baobabsoluciones/aco/tsplib"
	"github.com/pkg/errors"
)

var (
	tsp = flag.String("tsp", "", "TSP problem file in TSPLIB format")
	log = flag.String("log", "", "(optional)file where progress should be logged")
	seed = flag.Int64("seed", 0, "random seed for the algorithm; if 0 or unset a random seed will be generated; same seed and TSP problem will result in the same solution")
	alpha = flag.Float64("alpha", 1, "alpha and beta control the relative importance of trail versus visibility.")
	beta = flag.Float64("beta", 5, "alpha and beta control the relative importance of trail versus visibility.")
	rho = flag.Float64("rho", 0.5, "rho is a coefficient such that (1 - rho) represents the evaporation of trail between time t and t+n")
	Q = flag.Float64("Q", 100, "Q is a constant used in the LayTrail function to determine the amount of trail to be spread by each ant.")
	NCmax = flag.Int("NCmax", 5000, "NCmax is the maximum number of ant cycles.")
)

func main() {
	flag.Parse()

	// TODO change this, it can be misleading if the user really wants to set seed 0
	if *seed == 0 {
		*seed = time.Now().UTC().UnixNano()
	}

	var logWriter io.Writer

	if *log == "" {
		logWriter = ioutil.Discard
	} else {
		logFile, err := os.Create(*log)
		if err != nil {
			fmt.Printf("failed to create log file %s\nAborting...", *log)
			os.Exit(1)
		}

		logWriter = logFile

		// TODO handle possible error of Close()
		defer func() {
			err := logFile.Close()
			if err != nil {
				// TODO return this errors together with the other return values, wrap the other errors with this one
				err = errors.Wrapf(err, "failed to close solver progress log file: %s", *log)
				fmt.Println(err)
			}
		}()
	}

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

		var nAnts int = len(g.Vertices)
		trailUpdateFunc := aco.LayTrailAntCycle

		s, _, err := aco.AntSystemAlgorithm(g, nAnts, *NCmax, *Q, *rho, *alpha, *beta, trailUpdateFunc, *seed, logWriter)
		if err != nil {
			// TODO
			panic(err)
		}

		fmt.Printf("%d,%f,%f,%f,%f,%f\n", *NCmax, *Q, *rho, *alpha, *beta, aco.CompTotLength(g, s))

		//TODO find a good ouptut format
		/*
		fmt.Printf("solution: %v\n", s)
		fmt.Printf("solution length: %f\n", aco.CompTotLength(g, s))
		if stagBehv {
			fmt.Println("AS terminated with stagnation behaviour")
		}
		*/
	}
}
