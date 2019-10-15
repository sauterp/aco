// solves a TSP with ACO as described in Dorigo et al. 96
package main

import (
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"time"

	"bitbucket.org/baobabsoluciones/aco"
	"bitbucket.org/baobabsoluciones/aco/tsplib"
	"github.com/pkg/errors"
)


var (
	tsp  = flag.String("tsp", "", "TSP problem file in TSPLIB format") 
	log  = flag.String("log", "", "file where progress should be logged")
	seed = flag.Int64("seed", 0, "random seed for the algorithm; if 0 or unset a random seed will be generated; same seed and TSP problem will result in the same solution")
)

func main() {
	flag.Parse()

	if *seed == 0 {
		*seed = time.Now().UTC().UnixNano()
	}

	var logWriter io.Writer

	if *log == "" {
		fmt.Println("-log flag not specified; solver progress will not be logged")
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
		s, stagBehv, err := aco.ASBestParams(g, *seed, logWriter)
		if err != nil {
			// TODO
			panic(err)
		}

		fmt.Printf("solution: \n")
		for i := 0; i < len(s); i++ {
			fmt.Println(*(s[i]))
		}

		fmt.Printf("solution length: %f\n", aco.CompTotLength(g, s))
		if stagBehv {
			fmt.Println("AS terminated with stagnation behaviour")
		}
	}
}
