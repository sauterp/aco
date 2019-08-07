package aco

import (
	"bufio"
	"encoding/csv"
	"io"
	"os"
	"strconv"
	"testing"
)

func check(e error) {
	if e != nil {
		panic(e)
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
			Index: len(oliver30Graph.Edges),
			Label: newCity.Label,
		}
		oliver30Graph.Vertices = append(oliver30Graph.Vertices, newVertex)
		lastEdgeIndex := len(oliver30Graph.Edges) - 1
		for i := 0; i < len(cities); i++ {
			oliver30Graph.Edges[lastEdgeIndex][i].Length = CompEuclid2dDist(cities[i].X, cities[i].Y, newCity.X, newCity.Y)
		}

		cities = append(cities, newCity)
	}
	b.Log(oliver30Graph)
}
