// package aco provides the Ant System algorithm according to Dorigo, Maniezzo and Colorni 1996:
// The Ant System: Optimization by a colony of cooperating agents
// Marco Dorigo, Vittorio Maniezzo and Alberto Colorni
// IEEE Transactions on Systems, Man, and Cybernetics–Part B, Vol.26, No.1, 1996, pp.1-13
package aco

import (
	"fmt"
	"math"
	"math/rand"
)

// Graph holds an undirected and fully connected graph represented as a triangular adjacency matrix.
type Graph struct {
	Vertices []Vertex
	// n := len(Vertices)
	// In a fully connected graph that is used to look up the edges connected to a Vertex very frequently, the most efficient way to store the Edges is by means of an triangular adjacency matrix.
	Edges [][]Edge
}

// GetEdge retrieves the edge (vi, vj) from graph.Edges where vi and vj are Vertex.Index.
// GetEdge(vi, vj) == GetEdge(vj, vi)
func (graph *Graph) GetEdge(vi, vj int) (*Edge, error) {
	if vi == vj {
		return nil, fmt.Errorf("vi == vj == %d; graph does not have circular Edges", vi)
	}
	if vj < vi {
		return graph.Edges[vi][vj], nil
	}
	if vj > vi {
		return graph.Edges[vj][vi], nil
	}
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
	// An ant has to visit all n cities
	// If an ant has already visited a city we add it to the TabuList and cannot visit it again
	TabuList Tour
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

// MoveToNextVertex chooses the town to go to with a probability that is a function of the town distance and of the amount of trail present on the connecting edge
// TODO [#A]
func (ant *Ant) MoveToNextVertex() {
	// TODO many of the following computations are parallelizable and many of them can be precomputed
	normFact := 0.0
	pos := ant.Position
	//TODO find a better way to do this
	for _, v := range graph.Vertices {
		isInTabuList := false
		for tv := range ant.TabuList {
			if tv == v {
				isInTabuList = true
				break
			}
		}
		if !isInTabuList {
			edge = problemGraph.GetEdge(pos, v)
			normFact += pow(edge.TrailIntensity, alpha) * pow(edge.Visibility, beta)
		}
	}
	transProbs := make([]float64, len(problemGraph.Vertices))
	for tpi not in ant.TabuList {
		edge = problemGraph.GetEdge(pos, v)
		transProb[tpi] := pow(edge.TrailIntensity, alpha) * pow(edge.Visibility, beta) * normFact
	}
	// the ant cannot stay in its current position
	// TODO find a more elegant solution for this
	transProbs[pos] := 0.0
	new_pos := // random selection based on transProbs
	ant.Tour = append(ant.Tour, problemGraph.GetEdge(pos, new_pos))
	ant.TabuList = append(ant.TabuList, new_pos)
	ant.Position = new_pos
}

// CompTotTourLen computes the total length of this ant's tour
func CompTotLength(tour Tour) float64 {
	totLength := 0.0
	for i := range tour {
		totLength += tour.Length
	}
	return totLength
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
	// TODO is a check for rho > 0 necessary?
	if rho >= 1 {
		return nil, fmt.Errorf("rho >= 1.0")
	}

	err := CheckFullyConnected(problemGraph)
	if err != nil {
		return nil, err
	}

	// compute Visibility of all Edges
	// set all TrailIntensities to rho
	for i := range problemGraph.Edges {
		for j := range problemGraph.Edges[i] {
			edge := &problemGraph.Edges[i][j]
			edge.Visibility = 1.0 / edge.Length
			edge.TrailIntensity = rho
		}
	}


	// begin the Ant Cycle algorithm

	t := 0

	// set ants on Vertices using a uniform random distribution
	ants := make([]Ant, 0, nAnts)
	nVertices := len(problemGraph.Vertices)
	for i := 0; i < nAnts; i++ {
		ants = append(ants, Ant{
			TabuList: make(Tour, 0, nVertices),
		})
		newAnt := &ants[i]
		firstPos := &problemGraph.Vertices[rand.Intn(nVertices)]
		newAnt.Position = firstPos
		newAnt.TabuList = append(newAnt.TabuList, firstPos)
	}

	for _ := range NCmax {
		// Before computing the trails to be excreted before the next ant-cycle, we need to wait until all goroutines have finished their work.
		var wg sync.WaitGroup
		for ai := range ants {
			wg.Add(1)
			go func() {
				defer wg.Done()
				for _ := range problemGraph.Vertices {
					ants[ai].MoveToNextVertex()
				}
				// one ant has finished a tour
			} ()
		}

		// Wait until all ants have completed a tour
        		wg.Wait()

		// TODO [#A]
		trailUpdateFunc(problemGraph, ants)

		// save the shortestTour found by the ants
		shortestTour := ants[0].TabuList
		shortestLength := CompTotLength(shortestTour)
		for i := 1; i < len(ants); i++ {
			totLength := CompTotLength(ants[i].TabuList)
			if totLength < shortestLength {
				shortestTour := ants[i].TabuList
			}
		}

		for i := range ants {
			// TODO [#A]
			ants[i].EmptyTabuList()
		}

		// if all ants did the same tour, abort, stagnation behaviour, no alternative solutions will be explored
		// to test, whether all ants have found the same tour, we need to make C = (m * (m - 1)) / 2 comparisons between the tours.
		// the definition of a triangular number is (t * (t + 1)) / 2 now substitute t = m - 1 since we don't need to compare a tour with itself and you will obtain the above number
		// allocate a buffered channel of C boolean values
		C := (nAnts * (nAnts - 1)) / 2
		quit := make(chan struct{})
		comps := make(chan bool, C)
		// TODO see whether this loop can be optimized by not creating any goroutines as soon as <-quit
		for i := 0; i < nAnts; i++ {
			go func() {
				for j := i+1; j < nAnts; j++ {
					go func() {
                				select {
                				case <-quit:
                				        return
                				default:
							// if two ants have differing tours, all other comparisons can be aborted
							// TODO ensure that this comparison really compares both Tours elementwise and not the slice pointers
							if ants[i].Tour != ants[j].Tour {
								quit <- struct{}{}
							} else {
								comps <- true
							}
                				}
				}
			} ()
		}

		// TODO think about whether it makes sense to return information about when stangnation behaviour started
		stagnationBehaviour := true
		for comp := comps {
			if !comp {
				stagnationBehaviour = false
				break
			}
		}
		if stagnationBehaviour {
			return shortestTour, stagnationBehaviour, nil
		}

		// TODO do we really need to close these channels?
		close(comps)
		close(quit)
	}
	stagnationBehaviour := false
	return shortestTour, stagnationBehaviour, nil
}
