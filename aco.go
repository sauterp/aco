// package aco provides the Ant System algorithm according to Dorigo, Maniezzo and Colorni 1996:
// The Ant System: Optimization by a colony of cooperating agents
// Marco Dorigo, Vittorio Maniezzo and Alberto Colorni
// IEEE Transactions on Systems, Man, and Cybernetics–Part B, Vol.26, No.1, 1996, pp.1-13
package aco

import (
	"fmt"
	"math"
)

// Graph holds an undirected and fully connected graph represented as a triangular adjacency matrix.
type Graph struct {
	Vertices []Vertex
	// n := len(Vertices)
	// In a fully connected graph that is used to look up the edges connected to a Vertex very frequently, the most efficient way to store the Edges is by means of an triangular adjacency matrix.
	Edges [][]Edge
}

// Vertex holds the integer index used to address the rows in the Graph.Edges matrix and a string Label.
type Vertex struct {
	Index int
	// Label has no meaning to the Ant System algorithm it is only used to output the result.
	Label string
}

// Edge holds one edge of a Graph with Length information
type Edge struct {
	Length float64
	// Visibility increases the chance that an Ant takes this Edge. Usually Visibility = 1.0 / Length
	Visibility     float64 // Will be overwritten by Ant System Algorithm
	TrailIntensity float64 // Will be overwritten by Ant System Algorithm
}

// CompEuclid2dDist computes the euclidean distance between two 2 dimensional points a and b.
func CompEuclid2dDist(aX, aY, bX, bY float64) float64 {
	abXdist := aX - bX
	abYdist := aY - bY
	dist := math.Sqrt(abXdist*abXdist + abYdist*abYdist)
	return dist
}

type Tour []*Vertex

// Ant holds the state of one ant, the basic agent of this algorithm
type Ant struct {
	Position *Vertex
	Tour     Tour
	// An ant has to visit all n cities
	// If an ant has already visited a city we add it to the TabuList and cannot visit it again
	TabuList []*Vertex
}

// TODO implement func NewAnt()
// initialize TabuList with make([]*Vertex, n) or make([]*Vertex, 0, n)

// CheckFullyConnected will return non-nil error if graph is not fully connected
func CheckFullyConnected(graph Graph) error {
	errMsg := ""
	nVertices := len(graph.Vertices)

	if len(graph.Edges) != nVertices {
		errMsg += fmt.Sprintf("graph.Edges does not contain enough columns to act as an adjacency matrix for all nVertices = %d", nVertices)
	}

	for i := range graph.Edges {
		col := graph.Edges[i]
		if len(col) != i {
			errMsg += fmt.Sprintf("too many edges in column %d; want len(col) == %d; got len(col) == %d", i, i, len(col))
		}
		for j := range col {
			cell := col[j]
			if cell.Length < 0 {
				errMsg += fmt.Sprintf("Edges[%d][%d].Length < 0; Edges need to have non-negative Length", i, j)
			}
			if cell.Length == 0 {
				errMsg += fmt.Sprintf("Edges[%d][%d].Length == 0; all Vertices need to be connected", i, j)
			}
		}
	}

	if errMsg != "" {
		return fmt.Errorf(errMsg)
	} else {
		return nil
	}
}

// AntSystemAlgorithm is the main method for initiating the Ant System algorithm
func AntSystemAlgorithm(
	problemGraph Graph,
	nAnts int,
	NCmax int,
	// Q is a constant used in the LayTrail function to determine the amount of trail to be spread by each ant.
	Q float64,
	// rho is a coefficient such that (1 - rho) represents the evaporation of trail between time t and t+n
	// The coefficient rho must be set to a value < 1 to avoid unlimited accumulation of trail
	// In Dorigo96 a <<small>> constant c was used. rho = c
	rho float64,
	// alpha and beta control the relative importance of trail versus visibility. (autocataclytic process)
	alpha, beta float64,
	trailUpdateFunc *func(Graph, []Ant),
) (Tour, error) {
	err := CheckFullyConnected(problemGraph)
	if err != nil {
		return nil, err
	}
	return Tour{}, nil
}
