#include "graph.h"



int main()
{
	// creating directed graph
	graph g(true);

	// read city properties from txt file
	g.readSEHIR();

	// read csv files which stores each city's shortest distance
	// to other cities
	g.readCSV();

	// read map points for cities to plot path.
	g.readMAP();

	// create graph with these files.
	// graph stores city and its neighboor cities and distance between them.
	g.make_graph();

	// create start point and target points.
	// these points are plate number of cities.
	vector<int> targets = {6, 32, 35, 25};
	int start = 41;

	// find minimum path with genetic algorithm.
	g.minPath(start, targets);

	// print the map with path.
	g.print_map(0.7, true);


	system("pause");
	return 0;
}