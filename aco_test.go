// Tests for the Ant System Algorithm described in Dorigo et al. 96
package aco

import (
	"fmt"
	"math"
	"os"
	"testing"
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
func TestMoveToNextVertex(t *testing.T) {
}

// CheckSolutionValid checks that all Vertices in proglemGraph are visited exactly once.
func CheckSolutionValid(solution Tour, proglemGraph Graph) error {
	errMsg := ""

	if len(solution) == 0 {
		errMsg += "solution is empty"
	}

	// Check that no pointer is nil
	solContainsNilPointer := false
	for vi, v := range solution {
		if v == nil {
			solContainsNilPointer = true
			errMsg += fmt.Sprintf("solution[%d] == nil\n", vi)
		}
	}

	if !solContainsNilPointer {
		// Check that every Vertex is visited exactly once.
		// Check that all Vertices in solution are really in the problemGraph
		for vi, v := range solution {
			for vj := vi + 1; vj < len(solution); vj++ {
				if v == solution[vj] {
					errMsg += fmt.Sprintf("Vertex %v appears multiple times in solution\n", *v)
					break
				}
			}
			vFoundInProblemGraph := false
			for pgvi := 0; pgvi < len(proglemGraph.Vertices); pgvi++ {
				if v == &proglemGraph.Vertices[pgvi] {
					vFoundInProblemGraph = true
					break
				}
			}
			if !vFoundInProblemGraph {
				errMsg += fmt.Sprintf("Vertex %v is not in problemGraph\n", *v)
			}
		}
	}

	if errMsg == "" {
		return nil
	} else {
		return fmt.Errorf(errMsg)
	}
}

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

	solution, stagnationBehaviour, err := AntSystemAlgorithm(
		triangleGraph,
		len(triangleGraph.Vertices),
		NCmax,
		Q,
		rho,
		alpha, beta,
		LayTrailAntCycle,
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
			{0, "0"},
			{1, "1"},
			{2, "2"},
			{3, "2"},
		},
		Edges: [][]Edge{
			{},
			{{1, 1, 0}},
			{{sqrt2, invSqrt2, 0}, {1, 1, 0}},
			{{1, 1, 0}, {sqrt2, invSqrt2, 0}, {1, 1, 0}},
		},
	}

	// run AS
	solution, _, err := ASBestParams(
		squareGraph,
		os.Stdout,
	)
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
