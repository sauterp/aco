// Tests for the Ant System Algorithm described in Dorigo et al. 96
package aco

import (
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"testing"
	"time"
)

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func TestEqualTour(t *testing.T) {
	vp := func(ind int) *Vertex {
		v := new(Vertex)
		v.Index = ind
		return v
	}

	t.Run("SameTour", func(t *testing.T) {
		a := Tour{vp(0), vp(1), vp(2)}
		if !EqualTour(a, a) {
			t.Fail()
		}
	})

	t.Run("SameTourDifferentObject", func(t *testing.T) {
		a := Tour{vp(0), vp(1), vp(2)}
		b := Tour{a[0], a[1], a[2]}
		if !EqualTour(a, b) {
			t.Fail()
		}
	})

	t.Run("SameTourDifferentOrder", func(t *testing.T) {
		a := Tour{vp(0), vp(1), vp(2)}
		b := Tour{a[1], a[2], a[0]}
		if !EqualTour(a, b) {
			t.Fail()
		}
	})

	t.Run("DiffTourSameLen", func(t *testing.T) {
		a := Tour{vp(0), vp(1), vp(2)}
		b := Tour{a[1], vp(3), a[0]}
		if EqualTour(a, b) {
			t.Fail()
		}
	})

	t.Run("DiffTourDiffLen", func(t *testing.T) {
		a := Tour{vp(0), vp(1), vp(2), vp(4)}
		b := Tour{a[1], vp(3), a[0]}
		if EqualTour(a, b) {
			t.Fail()
		}
	})

	// Two Tours are equal if their order is in reverse, that is because our problem Graph is undirected.
	t.Run("SameTourReversed", func(t *testing.T) {
		a := Tour{vp(0), vp(1), vp(2), vp(3)}
		b := Tour{a[3], a[2], a[1], a[0]}
		if !EqualTour(a, b) {
			t.Fail()
		}
	})

	t.Run("SameTourReversedDiffOrder", func(t *testing.T) {
		a := Tour{vp(0), vp(1), vp(2), vp(3), vp(4), vp(5)}
		b := Tour{a[5], a[4], a[3], a[2], a[1], a[0]}
		if !EqualTour(a, b) {
			t.Fail()
		}
	})
}

// TODO
/*
func TestMoveToNextVertex(t *testing.T) {
	t.Error("TODO implement")
}
*/

// TODO find a way to test replicability

func TestCheckFullyConnected(t *testing.T) {
	t.Run("OneVertex", func(t *testing.T) {
		// create triangle graph
		oneVertexG := Graph{
			Vertices: []Vertex{
				{1, "1"},
			},
			Edges: [][]Edge{
				{},
			},
		}
		err := CheckFullyConnected(oneVertexG)
		if err != nil {
			t.Fatal(err)
		}
	})

	t.Run("TwoVerticesNoEdge", func(t *testing.T) {
		// create triangle graph
		twoVertG := Graph{
			Vertices: []Vertex{
				{0, "0"},
				{1, "1"},
			},
			Edges: [][]Edge{
				{},
				{},
			},
		}
		err := CheckFullyConnected(twoVertG)
		if err == nil {
			t.Fail()
		}
	})

	t.Run("TwoVerticesOneEdge", func(t *testing.T) {
		// create triangle graph
		twoVertG := Graph{
			Vertices: []Vertex{
				{0, "0"},
				{1, "1"},
			},
			Edges: [][]Edge{
				{},
				{{1, 1, 0}},
			},
		}
		err := CheckFullyConnected(twoVertG)
		if err != nil {
			t.Fail()
		}
	})

	t.Run("ThreeVerticesTwoEdges", func(t *testing.T) {
		// create triangle graph
		threeVertG := Graph{
			Vertices: []Vertex{
				{0, "0"},
				{1, "1"},
				{2, "2"},
			},
			Edges: [][]Edge{
				{},
				{{1, 1, 0}},
				{{1, 1, 0}},
			},
		}
		err := CheckFullyConnected(threeVertG)
		if err == nil {
			t.Fail()
		}
	})
}

// TestTriangle tests the AS on a triangle graph. The triangle is the most trivial TSP and all ants will find the exact same tour, therefore the AS must terminate due to stagnation behaviour after the first cycle.
func TestTriangle(t *testing.T) {
	// create triangle graph
	triangleGraph := Graph{
		Vertices: []Vertex{
			{0, "0"},
			{1, "1"},
			{2, "2"},
		},
		Edges: [][]Edge{
			{},
			{{1, 1, 0}},
			{{1, 1, 0}, {1, 1, 0}},
		},
	}

	// TODO determine parameters
	var NCmax int = 1 // check whether AS terminated with stagnation behaviour after exactly 1 cycle
	var Q float64 = 1
	var rho float64 = 0.5
	var alpha float64 = 1
	var beta float64 = 1
	var seed int64 = 0

	solution, stagnationBehaviour, err := AntSystemAlgorithm(
		triangleGraph,
		len(triangleGraph.Vertices),
		NCmax,
		Q,
		rho,
		alpha, beta,
		LayTrailAntCycle,
		seed,
		os.Stdout,
	)
	if err != nil {
		t.Error(err)
	}

	err = CheckSolutionValid(solution, triangleGraph)
	if err != nil {
		t.Fatal(err)
	}

	if !stagnationBehaviour {
		t.Error("AntSystemAlgorithm should have terminated with stagnationBehaviour == true")
	}
}

// TestSquare tests the AS on a square graph. The graph consists of four vertices arranged as a square and all of them are connected to eachother. There are exactly three possible tours for this TSP(without visiting any Vertex twice). Two tours involve both diagonals and one involves all boundaries of the square. The boundaries are the optimal solution, since they are shorter than the diagonals. The test expects that AS returns the square boundaries solution.
// Be aware that this test may fail, if all ants come up with the same non-border solution in a cycle. That would cause stagnation behaviour. Use a large number of ants to avoid make this event unlikely.
// TODO Convert this to a benchmark?
func TestSquare(t *testing.T) {
	// The diagonals of a unit square have length Sqrt(2)
	sqrt2 := math.Sqrt(2)
	invSqrt2 := 1 / sqrt2

	// create square graph
	squareGraph := Graph{
		Vertices: []Vertex{
			// It would have been more intuitive if the square vertices were label as this:
			// 0	1
			// 2	3
			// instead of
			// 0	1
			// 3	2
			{0, "0"},
			{1, "1"},
			{2, "3"},
			{3, "2"},
		},
		Edges: [][]Edge{ //This needs a propper explanation
			{},
			{{1, 1, 0}},
			{{sqrt2, invSqrt2, 0}, {1, 1, 0}},
			{{1, 1, 0}, {sqrt2, invSqrt2, 0}, {1, 1, 0}},
		},
	}
	fmt.Println(squareGraph)

	var seed int64 = 0

	// run AS
	solution, _, err := ASBestParams(
		squareGraph,
		seed,
		os.Stdout,
	)
	for i := 0; i < len(solution); i++ {
		fmt.Println(*solution[i])
	}

	if err != nil {
		t.Error(err)
	}

	err = CheckSolutionValid(solution, squareGraph)
	if err != nil {
		t.Fatal(err)
	}

	// check that AS returns the square boundaries solution.
	// check that the solution Edges are in the right order.
	// check that the solution Vertices are never visited twice.
	sv := squareGraph.Vertices
	expectedSol := Tour{&sv[0], &sv[1], &sv[2], &sv[3]}
	if !EqualTour(solution, expectedSol) {
		t.Errorf("want square border edges in solution:\n%v\ngot:\n%v\n", solution, expectedSol)
	}

	solTotLen := CompTotLength(squareGraph, solution)
	if solTotLen > 4.5 {
		t.Errorf("want total length ~= 4.0\ngot total length == %f\n", solTotLen)
	}
}

// This function generates a Graph representing a equidistant grid of nxn
func generateGridGraph(nGridNodesPerDim int, dist float64) Graph {
	nGridNodes := nGridNodesPerDim * nGridNodesPerDim
	g := Graph{
		Vertices: make([]Vertex, nGridNodes),
		Edges:    make([][]Edge, nGridNodes),
	}
	for vx := 0; vx < nGridNodesPerDim; vx++ {
		for vy := 0; vy < nGridNodesPerDim; vy++ {
			v := vx*nGridNodesPerDim + vy
			g.Vertices[v].Index = v
			g.Vertices[v].Label = fmt.Sprintf("%d", v)
			g.Edges[v] = make([]Edge, v)
			for e := 0; e < v; e++ {
				aX := float64(vx) * dist
				aY := float64(vy) * dist
				bX := float64(g.Vertices[e].Index/nGridNodesPerDim) * dist
				bY := float64(g.Vertices[e].Index%nGridNodesPerDim) * dist
				g.Edges[v][e].Length = CompEuclid2dDist(aX, aY, bX, bY)
			}
		}
	}
	return g
}

// TestCompTotLength checks whether the function CompTotLength returns the rigth length
func TestCompTotLength(t *testing.T) {
	nGridNodesPerDim := 4
	var dist float64 = 10
	// For a graph that is an equidistant grid of n x n fully connected vertices with distance d between
	// neighbouring vertices, where n is even, an optimal solution has length: d * n * n.

	g := generateGridGraph(nGridNodesPerDim, dist)

	gv := func(i int) *Vertex {
		return &g.Vertices[i]
	}
	optSol := Tour{
		gv(0), gv(4), gv(8), gv(12),
		gv(13), gv(9), gv(10), gv(14),
		gv(15), gv(11), gv(7), gv(3),
		gv(2), gv(6), gv(5), gv(1),
	}

	optSolLenInt := int(CompTotLength(g, optSol))
	if optSolLenInt != 160 {
		t.Errorf("want: %d\ngot: %d\n", 160, optSolLenInt)
	}
}

func BenchmarkGrid(b *testing.B) {
	// func TestGrid(b *testing.T) {
	// TODO generalize to parameter
	nGridNodesPerDim := 8

	// For a graph that is an equidistant grid of n x n fully connected vertices with distance d between
	// neighbouring vertices, where n is even, an optimal solution has length: d * n * n.
	var dist float64 = 1
	optLen := dist * float64(nGridNodesPerDim*nGridNodesPerDim)

	// generate grid
	g := generateGridGraph(nGridNodesPerDim, dist)

	seed := time.Now().UTC().UnixNano()
	var nAnts int = len(g.Vertices)
	var NCmax int = 10000
	var Q float64 = 100
	var rho float64 = 0.5
	var alpha float64 = 1
	var beta float64 = 100
	trailUpdateFunc := LayTrailAntCycle

	solution, _, err := AntSystemAlgorithm(g, nAnts, NCmax, Q, rho, alpha, beta, trailUpdateFunc, seed, ioutil.Discard)
	check(err)
	solLen := CompTotLength(g, solution)
	b.Logf("BestSol %f\n", solLen)
	b.Logf("OptSol %f\n", optLen)

	epsilon := 0.0000001
	if solLen < optLen-epsilon || solLen > optLen+epsilon { 
		// Why not if solLen != optLen
		b.Fail()
	}
}
