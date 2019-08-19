// create a line plot for two columns of a CSV file if they contain numerical data that can be parsed as float64
package main

import (
	"bufio"
	"encoding/csv"
	"flag"
	"fmt"
	"image/color"
	"io"
	"log"
	"os"
	"strconv"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

var (
	X     = flag.Int("X", -1, "column to be plotted on X axis(first column is 0)")
	Y     = flag.Int("Y", -1, "column to be plotted on Y axis(first column is 0)")
	i = flag.String("i", "", "name of the input file containing CSV formated data, the columns you want to parse need to be convertible to float64")
	o     = flag.String("o", "", "name of the file where the PNG containing the plot should be saved")
	title = flag.String("title", "", "title of the plot (optional)")
)

func Abort() {
	flag.Usage()
	os.Exit(1)
}

// TODO
func check(err error) {
	if err != nil {
		panic(err)
	}
}

func main() {
	flag.Parse()
	if *X < 0 || *Y < 0 {
		fmt.Println("please specify non-negative integer flags -X and -Y")
		Abort()
	}
	if *o == "" {
		fmt.Println("flag -o required")
		Abort()
	}
	if *i == "" {
		fmt.Println("flag -i required")
		Abort()
	}

	p, err := plot.New()
	if err != nil {
		log.Fatal(err)
	}

	p.Title.Text = *title
	p.Y.Scale = plot.LinearScale{}
	p.Y.Tick.Marker = plot.DefaultTicks{}

	f, err := os.Open(*i)
	check(err)
	r := csv.NewReader(bufio.NewReader(f))

	record, err := r.Read()
	p.X.Label.Text = record[*X]
	p.Y.Label.Text = record[*Y]

	var pts plotter.XYs = make([]plotter.XY, 0)

	for {
		record, err := r.Read()
		if err == io.EOF {
			break
		}
		check(err)

		newX, err := strconv.ParseFloat(record[*X], 64)
		check(err)
		newY, err := strconv.ParseFloat(record[*Y], 64)
		check(err)

		pts = append(pts, plotter.XY{X: newX, Y: newY})
	}

	line, scatter, err := plotter.NewLinePoints(pts)
	if err != nil {
		log.Panic(err)
	}

	line.Color = color.RGBA{R: uint8(255 / 2), A: 255}
	scatter.GlyphStyle.Radius = 0.1
	p.Add(line, scatter)

	err = p.Save(10*vg.Centimeter, 10*vg.Centimeter, *o)
	if err != nil {
		log.Panic(err)
	}
}
