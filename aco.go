// package aco provides the Ant System algorithm according to Dorigo, Maniezzo and Colorni 1996:
// The Ant System: Optimization by a colony of cooperating agents
// Marco Dorigo, Vittorio Maniezzo and Alberto Colorni
// IEEE Transactions on Systems, Man, and Cybernetics–Part B, Vol.26, No.1, 1996, pp.1-13
package aco

import (
	"fmt"
	"math"
	"math/rand"
	"sync"
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
	switch {
	case vi == vj:
		return nil, fmt.Errorf("vi == vj == %d; graph should not have circular Edges", vi)
	case vj < vi:
		return &graph.Edges[vi][vj], nil
	case vj > vi:
		return &graph.Edges[vj][vi], nil
	default:
		return nil, fmt.Errorf("this error should be impossible; none of the cases vi == vj, vj < vi or vj > vi have evaluated to true")
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

// EqualTour returns true if both tours contain the same vertices, in the same order and have the same total number of vertices
func EqualTour(a, b Tour) bool {
	if len(a) != len(b) {
		return false
	}

	// check whether a and b contain the same vertices and that they are in the same order
	// find the first two elements that point to the same Vertex
	var eqInd int
	// indicates whether an element pointing to the same Vertex has ben found in each of the list.
	eqElFound := false
	for i := 0; i < len(a); i++ {
		if a[0] == b[i] {
			eqInd = i
			eqElFound = true
			break
		}
	}
	if eqElFound {
		toursForwardEqual := true
		bOffset := eqInd
		for i := 0; i < len(a); i++ {
			// if the last element of b has been reached, continue comparing from the first element onwards
			if i+bOffset >= len(b) {
				bOffset = -i
			}
			if a[i] != b[i+bOffset] {
				toursForwardEqual = false
				break
			}
		}
		if toursForwardEqual {
			return true
		} else {
			// same procedure in reverse
			toursReverseEqual := true
			bOffset = eqInd
			for i := 0; i < len(a); i++ {
				// if the last element of b has been reached, continue comparing from the first element onwards
				if bOffset-i < 0 {
					bOffset = len(a) + bOffset
				}
				if a[i] != b[bOffset-i] {
					toursReverseEqual = false
					break
				}
			}
			return toursReverseEqual
		}
	}
	return false
}

// Ant holds the state of one ant, the basic agent of this algorithm
type Ant struct {
	Position *Vertex
	// An ant has to visit all n cities
	// If an ant has already visited a city we add it to the TabuList and cannot visit it again
	// TODO should Tour be implemented via a (circular) List data structure?
		// What about loops? how do we terminate a loop looping through a circular list? Elements are unique.
	TabuList Tour
}

// TODO implement func NewAnt()
// initialize TabuList with make([]*Vertex, n) or make([]*Vertex, 0, n)

// TODO generalize together with initialization
// TODO you can probably delete this
func (ant *Ant) EmptyTabuList() {
	// TODO ant.TabuList = make(Tour, 0, nVertices)
	ant.TabuList = make(Tour, 0)
}

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
				errMsg += fmt.Sprintf("Edges[%d][%d].Length == 0; all Vertices need to be connected; distance zero between vertices not allowed", i, j)
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
// TODO [#A] move into AntSystemAlgorithm func and remove graph parameter
func (ant *Ant) MoveToNextVertex(alpha, beta float64, graph Graph) error {
	// TODO many of the following computations are parallelizable and many of them can be precomputed
	// TODO find a better way to do this
	nVertices := len(graph.Vertices)
	possVerts := make([]*Vertex, 0, nVertices-len(ant.TabuList))
	for vi := range graph.Vertices {
		v := &graph.Vertices[vi]
		isInTabuList := false
		for _, tv := range ant.TabuList {
			if tv == v {
				isInTabuList = true
				break
			}
		}
		if !isInTabuList {
			possVerts = append(possVerts, v)
		}
	}

	if len(possVerts) == 0 {
		return fmt.Errorf("all graph.Vertices in TabuList; No new position assigned")
	}

	var newPos *Vertex = nil
	if len(possVerts) == 1 {
		newPos = possVerts[0]
	} else {
		// Generate transition probability dist. for possVerts based on Edge.TrailIntensity and Edge.Visibility.
		probs := make([]float64, len(possVerts))

		// compute normalization factor of formula (4) on page 6 of Dorigo et al. 96
		var normFact float64 = 0
		for pi := range probs {
			v := possVerts[pi]
			edge, err := graph.GetEdge(ant.Position.Index, v.Index)
			if err != nil {
				// TODO
				panic(err)
			}
			normFact += math.Pow(edge.TrailIntensity, alpha) * math.Pow(edge.Visibility, beta)
		}

		// compute transition probabilities
		for pi := range probs {
			v := possVerts[pi]
			edge, err := graph.GetEdge(ant.Position.Index, v.Index)
			if err != nil {
				// TODO
				panic(err)
			}
			probs[pi] = math.Pow(edge.TrailIntensity, alpha) * math.Pow(edge.Visibility, beta) * normFact
		}

		// sum each element of probs with all preceding elements to get a "ladder"
		for pi := 1; pi < len(probs); pi++ {
			probs[pi] += probs[pi-1]
		}

		// Select random next position based on prob. dist.
		r := rand.Float64()
		for pi := 0; pi < len(probs)-1; pi++ {
			if r < probs[pi] {
				newPos = possVerts[pi]
			}
		}
		if newPos == nil {
			newPos = possVerts[len(possVerts)-1]
		}
	}

	ant.TabuList = append(ant.TabuList, newPos)
	ant.Position = newPos

	return nil
}

// CompTotTourLen computes the total length of this ant's tour
func CompTotLength(graph Graph, tour Tour) float64 {
	totLength := 0.0
	for i := 0; i < len(tour)-1; i++ {
		edge, err := graph.GetEdge(tour[i].Index, tour[i+1].Index)
		if err != nil {
			// TODO
			panic(err)
		}
		totLength += edge.Length
	}

	// The last edge is not explicitly stored in the tour, because it leads back to the first Vertex of the tour.
	edge, err := graph.GetEdge(tour[len(tour)-1].Index, tour[0].Index)
	if err != nil {
		// TODO
		panic(err)
	}
	totLength += edge.Length

	return totLength
}

// LayTrailAntCycle when ant completes a tour, it lays a substance called trail on each edge visited.
// This is the main procedure used in the publication, referred to as ant-cycle, but they also proposed two alternatives LayTrailAntDensity and LayTrailAntQuantity on page 8
// TODO [#B] This computation needs to be done concurrently without race conditions
func LayTrailAntCycle(Q float64, graph Graph, ant Ant) {
	L_k := CompTotLength(graph, ant.TabuList)
	for i := 0; i < len(ant.TabuList)-1; i++ {
		v := ant.TabuList[i].Index
		vNext := ant.TabuList[i+1].Index
		edge, err := graph.GetEdge(v, vNext)
		if err != nil {
			// TODO
			panic(err)
		}
		edge.TrailIntensity += Q / L_k
	}
	firstV := ant.TabuList[0].Index
	lastV := ant.TabuList[len(ant.TabuList)-1].Index
	edge, err := graph.GetEdge(firstV, lastV)
	if err != nil {
		// TODO
		panic(err)
	}
	edge.TrailIntensity += Q / L_k
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
	// TODO convert func to a type
	trailUpdateFunc func(float64, Graph, Ant),
) (shortestTour Tour, stagnationBehaviour bool, err error) {
	if rho >= 1 || rho < 0 {
		return nil, false, fmt.Errorf("requirement: 0 <= rho < 1")
	}

	if alpha < 0 || beta < 0 {
		return nil, false, fmt.Errorf("alpha and beta need to be larger or equal to 0")
	}

	err = CheckFullyConnected(problemGraph)
	if err != nil {
		return nil, false, err
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

	// TODO [A#]
	// t := 0

	for nc := 0; nc < NCmax; nc++ {
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

		// Before computing the trails to be excreted before the next ant-cycle, we need to wait until all goroutines have finished their work.
		var wg sync.WaitGroup
		for ai := range ants {
			wg.Add(1)
			go func(ant *Ant) {
				defer wg.Done()
				for vi := 0; vi < nVertices-1; vi++ {
					err := ant.MoveToNextVertex(alpha, beta, problemGraph)
					if err != nil {
						// TODO
						panic(err)
					}
				}
				// one ant has finished a tour
			}(&ants[ai])
		}

		// Wait until all ants have completed a tour
		wg.Wait()

		// TODO [#A]
		// TODO ensure this is correct
		// evaporate trail
		for i := range problemGraph.Edges {
			for j := range problemGraph.Edges[i] {
				edge := &problemGraph.Edges[i][j]
				edge.TrailIntensity = rho * edge.TrailIntensity
			}
		}
		// lay new trail
		for i := range ants {
			// TODO clean
			trailUpdateFunc(Q, problemGraph, ants[i])
		}

		// save the shortestTour found by the ants
		shortestTour = ants[0].TabuList
		shortestLength := CompTotLength(problemGraph, shortestTour)
		for i := 1; i < len(ants); i++ {
			totLength := CompTotLength(problemGraph, ants[i].TabuList)
			if totLength < shortestLength {
				shortestTour = ants[i].TabuList
			}
		}

		// if all ants did the same tour, abort, stagnation behaviour, no alternative solutions will be explored
		// TODO see whether this loop can be optimized by not creating any goroutines as soon as <-quit
		// TODO introduce concurrency
		stagnationBehaviour := true
		for i := 0; i < nAnts; i++ {
			for j := i + 1; j < nAnts; j++ {
				// if two ants have differing tours, all other comparisons can be aborted
				if !EqualTour(ants[i].TabuList, ants[j].TabuList) {
					stagnationBehaviour = false
					break
				}
			}
			if !stagnationBehaviour {
				break
			}
		}

		if stagnationBehaviour {
			// TODO think about whether it makes sense to return information about when stangnation behaviour started
			return shortestTour, stagnationBehaviour, nil
		}
	}
	stagnationBehaviour = false
	return shortestTour, stagnationBehaviour, nil
}

// ASBestParams calls AntSystemAlgorithm with the best parameters found by Dorigo et al 96. page 9 for the ant-cycle algorithm
func ASBestParams(problemGraph Graph) (shortestTour Tour, stagnationBehaviour bool, err error) {
	var nAnts int = len(problemGraph.Vertices)
	var NCmax int = 5000
	var Q float64 = 100
	var rho float64 = 0.5
	var alpha float64 = 1
	var beta float64 = 5
	trailUpdateFunc := LayTrailAntCycle

	return AntSystemAlgorithm(problemGraph, nAnts, NCmax, Q, rho, alpha, beta, trailUpdateFunc)
}