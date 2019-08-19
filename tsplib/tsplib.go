// package tsplib provides functionality for parsing TSP problems published as TSPLIB in its own format.
// official documentation
// [TSPLIBdoc] https://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/tsp95.pdf
package tsplib

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

func ParseTSPLIBProblem(probFilename string) Graph {
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