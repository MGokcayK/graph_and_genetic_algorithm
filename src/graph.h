#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <random>
#include <ctime>
#include <opencv2/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>

using namespace std;

// creating custom structure to store population data
struct gPop
{
	vector<int> dna;
	int dist;
	float fitness;
};

// definition of pair for graph
typedef pair<int, float> gTarget;

// range function to create vector from `start` to `end`. Incrementation is 1.
template <typename T>
vector<T> range(T start, T end) {
	size_t N = (int)floor(end - start);
	vector<T> vec;
	vec.resize(N);
	iota(vec.begin(), vec.end(), start);
	return vec;
}


// main class
class graph
{

private:

	// parameters for graph
	int maxNodeInd = 0;
	bool directed;
	vector<vector<gTarget>> adjList;
	vector<string> city;
	vector<int> city_plate;
	vector<vector<string>> neighboor;
	vector<vector<float>> distance;
	vector<int> sub_dest;
	vector<int> x_loc;
	vector<int> y_loc;
	float cost;

	// parameters for genetic algo.
	int popMember = 20;
	float pCross = 0.3;
	float pMut = 0.25;

	int b_start;
	int minVal = 10000000;

	vector<gPop> population;
	void gen_pop(int target_size);
	float get_fitness(int dist);
	int get_distance(int start, vector<int> target, vector<int> pop);
	gPop get_child(gPop p1, gPop p2, vector<int> targets);
	vector<int> crossover(gPop p1, gPop p2);
	vector<int> mutation(gPop child);

public:

	graph(bool directed);

	~graph();

	void addNode(int s_Ind, int d_Ind, float value);

	void readSEHIR();

	void readCSV();

	void readMAP();

	void make_graph();

	void print_graph();

	void print_map(float scale, bool save);

	void set_genetic_params(int popMember, float pCross, float pMut);

	vector<int> A_star(int s_Ind, int d_Int);

	void minPath(int start, vector<int> targets);



};

