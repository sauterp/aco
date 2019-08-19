// Tests for the Ant System Algorithm described in Dorigo et al. 96
package aco

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"io"
	"math"
	"os"
	"strconv"
	"strings"
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

	// TODO determine parameters
	var NCmax int = 1000
	var Q float64 = 1
	var rho float64 = 0.5
	var alpha float64 = 1
	var beta float64 = 1
	var nAnts int = 1000
	trailUpdateFunc := func(float64, Graph, Ant) {}

	// run AS
	solution, _, err := AntSystemAlgorithm(
		squareGraph,
		nAnts,
		NCmax,
		Q,
		rho,
		alpha, beta,
		trailUpdateFunc,
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

// BenchmarkOliver30 benchmarks AS with the prolem Oliver30 from Dorigo et al. 96 and compares the performance to that stated in the paper.
func BenchmarkOliver30(b *testing.B) {
	type City struct {
		X     float64
		Y     float64
		Label string
	}
	cities := make([]City, 0, 30)

	var oliver30Graph Graph
	oliver30Graph.Vertices = make([]Vertex, 0, 30)

	f, err := os.Open("testdata/benchmarks/Oliver30.csv")
	check(err)
	r := csv.NewReader(bufio.NewReader(f))
	for {
		record, err := r.Read()
		if err == io.EOF {
			break
		}
		check(err)

		newX, err := strconv.ParseFloat(record[1], 64)
		check(err)
		newY, err := strconv.ParseFloat(record[2], 64)
		check(err)
		newCity := City{
			X:     newX,
			Y:     newY,
			Label: record[0],
		}

		// The first column of the oliver30Graph.Edges 2d array will have length 0 which is correct, since it is an upper triangular matrix and no Vertex has a circular Edge.
		oliver30Graph.Edges = append(oliver30Graph.Edges, make([]Edge, len(cities)))
		newVertex := Vertex{
			Index: len(oliver30Graph.Vertices),
			Label: newCity.Label,
		}
		oliver30Graph.Vertices = append(oliver30Graph.Vertices, newVertex)
		lastEdgeIndex := len(oliver30Graph.Edges) - 1
		for i := 0; i < len(cities); i++ {
			oliver30Graph.Edges[lastEdgeIndex][i].Length = CompEuclid2dDist(cities[i].X, cities[i].Y, newCity.X, newCity.Y)
		}

		cities = append(cities, newCity)
	}

	// TODO determine parameters
	var NCmax int = 1000
	var Q float64 = 1
	var rho float64 = 0.5
	var alpha float64 = 1
	var beta float64 = 1
	trailUpdateFunc := LayTrailAntCycle

	// run AS
	solution, _, err := AntSystemAlgorithm(
		oliver30Graph,
		len(oliver30Graph.Vertices),
		NCmax,
		Q,
		rho,
		alpha, beta,
		trailUpdateFunc,
	)
	if err != nil {
		b.Error(err)
	}

	err = CheckSolutionValid(solution, oliver30Graph)
	if err != nil {
		b.Fatal(err)
	}

	b.Logf("total length: %f", CompTotLength(oliver30Graph, solution))
}

// TSPLIB
// official documentation
// [TSPLIBdoc] https://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/tsp95.pdf

// TODO Geographical distance
// TODO Distance for crystallography problems

// nint is the function used by TSPLIB to round to the nearest integer
func nint(x float64) int {
	return int(x + 0.5)
}

// ceilInt returns the least integer value greater than or equal to x.
func ceilInt(f float64) int {
	ni := nint(f)
	var ci int
	if float64(ni) < f {
		ci = ni + 1
	} else {
		ci = ni
	}
	return ci
}

// ATTdist computes the pseudo-Euclidean distance function according to [TSPLIBdoc] section 2.5
// The edge weight type ATT corresponds to a special “pseudo-Euclidean” distance function.
// Let x[i] and y[i] be the coordinates of node i. The distance between two points i and j
// is computed as follows:
//     xd = x[i] - x[j];
//     yd = y[i] - y[j];
//     rij = sqrt( (xd*xd + yd*yd) / 10.0 );
//     tij = nint( rij );
//     if (tij<rij) dij = tij + 1;
//     else dij = tij;
func ATTdist(aX, aY, bX, bY float64) int {
	xd := aX - bX
	yd := aY - bY
	rab := math.Sqrt((xd*xd + yd*yd) / 10)
	dab := ceilInt(rab)
	return dab
}

// EUC2Ddist computes the Euclidean distance function according to [TSPLIBdoc] section 2.1
// For edge weight type EUC_2D and EUC_3D, floating point coordinates must be specified for
// each node. Let x[i], y[i], and z[i] be the coordinates of node i.
// In the 2-dimensional case the distance between two points i and j is computed as follows:
//     xd = x[i] - x[j];
//     yd = y[i] - y[j];
//     dij = nint( sqrt( xd*xd + yd*yd) );
// ...
// where sqrt is the C square root function.
func EUC2Ddist(aX, aY, bX, bY float64) int {
	xd := aX - bX
	yd := aY - bY
	dab := nint(math.Sqrt(xd*xd + yd*yd))
	return dab
}

// EUC3Ddist computes the Euclidean distance function according to [TSPLIBdoc] section 2.1
// For edge weight type EUC_2D and EUC_3D, floating point coordinates must be specified for
// each node. Let x[i], y[i], and z[i] be the coordinates of node i.
// ...
// In the 3-dimensional case we have:
//     xd = x[i] - x[j];
//     yd = y[i] - y[j];
//     zd = z[i] - z[j];
//     dij = nint( sqrt( xd*xd + yd*yd + zd*zd) );
// where sqrt is the C square root function.
func EUC3Ddist(aX, aY, aZ, bX, bY, bZ float64) int {
	xd := aX - bX
	yd := aY - bY
	zd := aZ - bZ
	dab := nint(math.Sqrt(xd*xd + yd*yd + zd*zd))
	return dab
}

// CEIL2Ddist computes the ceiling of the Euclidean distance function according to [TSPLIBdoc] section 2.6
// The edge weight type CEIL_2D requires that the 2-dimensional Euclidean distances is
// rounded up to the next integer.
func CEIL2Ddist(aX, aY, bX, bY float64) int {
	xd := aX - bX
	yd := aY - bY
	dab := ceilInt(math.Sqrt(xd*xd + yd*yd))
	return dab
}

// MAN2Ddist computes the Manhattan distance (L1-metric) function according to [TSPLIBdoc] section 2.2
// Distances are given as Manhattan distances if the edge weight type is MAN_2D or MAN_3D.
// They are computed as follows.
// 2-dimensional case:
//     xd = abs( x[i] - x[j] );
//     yd = abs( y[i] - y[j] );
//     dij = nint( xd + yd );
func MAN2Ddist(aX, aY, bX, bY float64) int {
	xd := math.Abs(aX - bX)
	yd := math.Abs(aY - bY)
	dab := nint(xd + yd)
	return dab
}

// MAN3Ddist computes the Manhattan distance (L1-metric) function according to [TSPLIBdoc] section 2.2
// Distances are given as Manhattan distances if the edge weight type is MAN_2D or MAN_3D.
// They are computed as follows.
// ...
// 3-dimensional case:
//     xd = abs( x[i] - x[j] );
//     yd = abs( y[i] - y[j] );
//     zd = abs( z[i] - z[j] );
//     dij = nint( xd + yd + zd );
func MAN3Ddist(aX, aY, aZ, bX, bY, bZ float64) int {
	xd := math.Abs(aX - bX)
	yd := math.Abs(aY - bY)
	zd := math.Abs(aZ - bZ)
	dab := nint(xd + yd + zd)
	return dab
}

// MAX2Ddist computes the Maximum distance (L∞-metric) function according to [TSPLIBdoc] section 2.3
// Maximum distances are computed if the edge weight type is MAX_2D or MAX_3D.
// 2-dimensional case:
//     xd = abs( x[i] - x[j] );
//     yd = abs( y[i] - y[j] );
//     dij = max( nint( xd ), nint( yd ) ) );
func MAX2Ddist(aX, aY, bX, bY float64) int {
	xd := nint(math.Abs(aX - bX))
	yd := nint(math.Abs(aY - bY))
	if xd > yd {
		return xd
	} else {
		return yd
	}
}

// MAX3Ddist computes the Maximum distance (L∞-metric) function according to [TSPLIBdoc] section 2.3
// Maximum distances are computed if the edge weight type is MAX_2D or MAX_3D.
// ...
// 3-dimensional case:
//     xd = abs( x[i] - x[j] );
//     yd = abs( y[i] - y[j] );
//     zd = abs( z[i] - z[j] );
//     dij = max( nint( xd ), nint( yd ), nint( zd ) );
func MAX3Ddist(aX, aY, aZ, bX, bY, bZ float64) int {
	xd := nint(math.Abs(aX - bX))
	yd := nint(math.Abs(aY - bY))
	zd := nint(math.Abs(aZ - bZ))

	var max int
	if xd > yd {
		max = xd
	} else {
		max = yd
	}
	if zd > max {
		max = zd
	}
	return max
}

func parseTSPLIBProblem(probFilename string) Graph {
	f, err := os.Open(probFilename)
	check(err)

	var dimension int
	var dist2DFunc func(float64, float64, float64, float64) int
	var dist3DFunc func(float64, float64, float64, float64, float64, float64) int
	const (
		DistDim2D = iota
		DistDim3D = iota
	)
	var distDim int

	s := bufio.NewScanner(f)
	for s.Scan() {
		t := s.Text()
		if t == "NODE_COORD_SECTION" {
			break
		}

		ts := strings.Split(t, ":")
		fieldName := strings.TrimSpace(ts[0])
		switch fieldName {
		case "DIMENSION":
			dimension64, err := strconv.ParseInt(strings.TrimSpace(ts[1]), 10, 32)
			check(err)
			dimension = int(dimension64) // TSPLIB doc states that ints have word length 32
		case "EDGE_WEIGHT_TYPE":
			fieldValue := strings.TrimSpace(ts[1])
			switch fieldValue {
			case "EUC_2D":
				dist2DFunc = EUC2Ddist
				distDim = DistDim2D
			case "EUC_3D":
				dist3DFunc = EUC3Ddist
				distDim = DistDim3D
			case "CEIL_2D":
				dist2DFunc = CEIL2Ddist
				distDim = DistDim2D
			case "ATT":
				dist2DFunc = ATTdist
				distDim = DistDim2D
			case "MAN_2D":
				dist2DFunc = MAN2Ddist
				distDim = DistDim2D
			case "MAN_3D":
				dist3DFunc = MAN3Ddist
				distDim = DistDim3D
			case "MAX_2D":
				dist2DFunc = MAX2Ddist
				distDim = DistDim2D
			case "MAX_3D":
				dist3DFunc = MAX3Ddist
				distDim = DistDim3D
			}
		}
	}
	err = s.Err()
	check(err)

	type City struct {
		X     float64
		Y     float64
		Label string
	}
	cities := make([]City, 0, dimension)

	var g Graph
	g.Vertices = make([]Vertex, 0, dimension)

	for s.Scan() {
		t := s.Text()
		if t == "EOF" {
			break
		}
		record := strings.Fields(t)
		check(err)

		newX, err := strconv.ParseFloat(record[1], 64)
		check(err)
		newY, err := strconv.ParseFloat(record[2], 64)
		check(err)
		newCity := City{
			X:     newX,
			Y:     newY,
			Label: record[0],
		}

		// The first column of the g.Edges 2d array will have length 0 which is correct, since it is an upper triangular matrix and no Vertex has a circular Edge.
		g.Edges = append(g.Edges, make([]Edge, len(cities)))
		newVertex := Vertex{
			Index: len(g.Vertices),
			Label: newCity.Label,
		}
		g.Vertices = append(g.Vertices, newVertex)
		lastEdgeIndex := len(g.Edges) - 1
		for i := 0; i < len(cities); i++ {
			var newLength float64
			if distDim == DistDim2D {
				newLength = float64(dist2DFunc(cities[i].X, cities[i].Y, newCity.X, newCity.Y))
			} else if distDim == DistDim3D {
				// TODO
				newLength = float64(dist3DFunc(0, 0, 0, 0, 0, 0))
				panic("parsing of 3D data not implemented")
			} else {
				panic(fmt.Sprintf("unknown distDim = %d", distDim))
			}
			g.Edges[lastEdgeIndex][i].Length = newLength
		}

		cities = append(cities, newCity)
	}

	return g
}

// TestTSPLIB benchmarks the AntSystemAlgorithm with the publicly available dataset TSPLIB:
// https://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/index.html
// the data can be found in benchmarks/TSPLIB/problems/
//
// The problems in TSPLIB are specified in different formats, which all need to be converted into a format suitable for AS.
func BenchmarkTSPLIB(b *testing.B) {
	// TSP problems with vertices that are connected but have distance 0 are excluded from our benchmarks:
	// a280

	benchmarks := []struct {
		problemName string
		optimalValue int
	}{
		{ "att48", 10628 },
		// "att532", // takes too long
	}

	for _, bm := range benchmarks {
		b.Run(bm.problemName, func(b *testing.B) {
			g := parseTSPLIBProblem("testdata/benchmarks/TSPLIB/problems/" + bm.problemName + ".tsp")
			// TODO determine parameters
			var NCmax int = 1000
			var Q float64 = 1
			var rho float64 = 0.5
			var alpha float64 = 1
			var beta float64 = 1
			trailUpdateFunc := LayTrailAntCycle

			// run AS
			solution, _, err := AntSystemAlgorithm(
				g,
				len(g.Vertices),
				NCmax,
				Q,
				rho,
				alpha, beta,
				trailUpdateFunc,
			)
			if err != nil {
				b.Error(err)
			}

			err = CheckSolutionValid(solution, g)
			if err != nil {
				b.Fatal(err)
			}

			b.Logf("total length: %f", CompTotLength(g, solution))
			b.Logf("known optimal value: %d", bm.optimalValue)
		})
	}
}
