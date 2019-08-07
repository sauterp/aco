// package aco provides the Ant System algorithm according to Dorigo, Maniezzo and Colorni 1996:
// The Ant System: Optimization by a colony of cooperating agents
// Marco Dorigo, Vittorio Maniezzo and Alberto Colorni
// IEEE Transactions on Systems, Man, and Cybernetics–Part B, Vol.26, No.1, 1996, pp.1-13
package aco

import "math"

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
	// Visibility increases the chance that an Ant takes this Edge. It needs to be precomputed for each Vertex before the Ant System Algorithm is started. Usually Visibility = 1.0 / Length
	Visibility float64
	TrailIntensity float64
}

// CompEuclid2dDist computes the euclidean distance between two 2 dimensional points a and b.
func CompEuclid2dDist(aX, aY, bX, bY float64) float64 {
	abXdist := aX - bX
	abYdist := aY - bY
	dist := math.Sqrt(abXdist*abXdist + abYdist*abYdist)
	return dist
}
