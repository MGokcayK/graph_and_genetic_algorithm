# graph_and_genetic_algorithm
Learning graph and genetic algorithm with TSP (traveling salesman problem) problem.

Main focus is learning graph structure. TSP problem is an example of it. 
Genetic Algorithm (GA) solution is not the best solution for
TSP problem. To understand how GA works, implemented in the project. 

Example:

```cpp
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

    // Targets : Ankara, Isparta, Izmir, Erzurum
	vector<int> targets = {6, 32, 35, 25};

    // Start : Kocaeli
	int start = 41;

	// find minimum path with genetic algorithm.
	g.minPath(start, targets);

	// print the map with path.
	g.print_map(0.7, true);


	system("pause");
	return 0;
}
```
Result in cmd :
```
 GENERATION : 2

Kocaeli (41) >> Yalova (77) >> Bursa (16) >> Balikesir (10) >> Izmir (35) || FROM Kocaeli to Izmir => KM : 453
Izmir (35) >> Manisa (45) >> Usak (64) >> Afyonkarahisar (3) >> Isparta (32) || FROM Izmir to Isparta => KM : 382
Isparta (32) >> Afyonkarahisar (3) >> Eskisehir (26) >> Ankara (6) || FROM Isparta to Ankara => KM : 420
Ankara (6) >> Kirikkale (71) >> Yozgat (66) >> Sivas (58) >> Erzincan (24) >> Erzurum (25) || FROM Ankara to Erzurum => KM : 874
Erzurum (25) >> Erzincan (24) >> Sivas (58) >> Tokat (60) >> Amasya (5) >> Corum (19) >> Cankiri (18) >> Bolu (14) >> Duzce (81) >> Sakarya (54) >> Kocaeli (41) || FROM Erzurum to Kocaeli => KM : 1116

Total KM : 3245
```

Printed Map :
[[data/path.png]]