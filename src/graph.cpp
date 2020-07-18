#include "graph.h"
#include <cassert>


// constructor. It set whether directed graph or not.
graph::graph(bool directed)
{
	this->directed = directed;
}

// destructor
graph::~graph()
{
}

// adding node to graph.
// s_Ind : start index
// d_Ind : destination index
// value : value of edge.

void graph::addNode(int s_Ind, int d_Ind, float value)
{
	// if index is bigger than maxNodeInd (maximum node index)
	// resize the adjancent list with maxNodeInd.
	if (d_Ind > this->maxNodeInd || s_Ind > this->maxNodeInd)
	{
		this->maxNodeInd = max(d_Ind,s_Ind);
		this->adjList.resize(this->maxNodeInd+1);
	}

	// add the pair
	this->adjList[s_Ind].push_back(make_pair(d_Ind, value));

	// if undirected graph, add same pair to destination index.
	if (this->directed == false)
	{
		this->adjList[d_Ind].push_back(make_pair(s_Ind, value));
	}
}


// printing graph
void graph::print_graph()
{
	for (int i = 0; i < this->adjList.size(); i++)
	{
		bool pr = false;
		for (gTarget t : this->adjList[i])
		{
			cout << "|| "<< i << " ---> " << t.first << " : " << t.second << ", ";
			pr = true;
		}

		if (pr) { cout << endl; }
	}
}

// reading city file which stores city plate number, its name and neighboorhood respectively.
void graph::readSEHIR()
{
	// load txt file
	ifstream txtFile;
	string line;
	txtFile.open("data/sehir.txt");
	vector<vector<string>> vecTXT;

	// read lines by storing vector of vector then got 2D vector.
	while (txtFile >> line)
	{
		stringstream ss(line);
		string token;
		vector<string> vecLine;
		while (getline(ss, token, ','))
		{
			vecLine.push_back(token);
		}

		vecTXT.push_back(vecLine);
	}


	// assign some properties to other parameters
	for (int i = 0; i < vecTXT.size(); i++)
	{

		vector<string> tmpNeigh;
		for (int j = 0; j < vecTXT[i].size(); j++)
		{
			if (j == 0)
			{
				string tmp = vecTXT[i][j];
				this->city_plate.push_back(stoi(tmp));
			}
			else if(j == 1)
			{
				this->city.push_back(vecTXT[i][j]);

			}
			else
			{
				tmpNeigh.push_back(vecTXT[i][j]);
			}
		}
		this->neighboor.push_back(tmpNeigh);
	}
}


// read map location file which store city location on image.
void graph::readMAP()
{
	// load txt files
	ifstream txtFile;
	string line;
	txtFile.open("data/map_loc.txt");
	vector<vector<string>> vecTXT;

	//read lines
	while (txtFile >> line)
	{
		stringstream ss(line);
		string token;
		vector<string> vecLine;
		while (getline(ss, token, ','))
		{
			vecLine.push_back(token);
		}

		vecTXT.push_back(vecLine);
	}

	// load locations separately
	for (int i = 0; i < vecTXT.size(); i++)
	{

		for (int j = 0; j < vecTXT[i].size(); j++)
		{
			if (j == 1)
			{
				string tmp = vecTXT[i][j];
				this->x_loc.push_back(stoi(tmp));
			}
			if (j == 2)
			{
				string tmp = vecTXT[i][j];
				this->y_loc.push_back(stoi(tmp));
			}
		}
	}
}


// read csv files which stores each city's shortest distance
void graph::readCSV()
{
	// load file
	ifstream csvFile;
	string line;
	csvFile.open("data/ilmesafe.csv");
	vector<vector<string>> vecTXT;

	// read lines then create vector of vector to have 2D vector.
	int lineCnt = 0;
	while (csvFile >> line)
	{
		int rowCnt = 0;
		stringstream ss(line);
		string token;
		vector<float> vecLine;
		while (getline(ss, token, ','))
		{
			if (lineCnt > 0)
			{
				if (rowCnt > 0)
				{
					if (token == "")
						token = "0.";
					vecLine.push_back(stof(token));
				}
			}
			rowCnt += 1;
		}
		if (lineCnt > 0)
		{
			this->distance.push_back(vecLine);
		}
		lineCnt += 1;
	}
}

// making graph with loaded files
void graph::make_graph()
{
	 //adding nodes
	for (int g = 0; g < this->city.size(); g++)
	{
		for (int h = 0; h < this->neighboor[g].size(); h++)
		{
			for (int k = 0; k < this->city.size(); k++)
			{
				if (this->neighboor[g][h] == this->city[k])
				{
					addNode(g, k, this->distance[g][k]);
				}
			}
		}
	}
}


// findin minimum path from start city to target cities.
void graph::minPath(int start, vector<int> targets)
{
	// substract 1 to each parameters bec of plate representation is not same as plate number.
	// Ex: Kocaeli's plate number is 41. Yet, 40 in data (because c++ start index is 0).
	for_each(targets.begin(), targets.end(), [](int& d) { d -= 1.0; });
	start -= 1;
	this->b_start = start;

	// generate population
	this->gen_pop(targets.size());

	// generation number
	int g = 1;

	// seach condition. If each generation has same shortest distance stop it.
	bool notsame = true;
	while(notsame)
	{
		g += 1;
		vector<gPop> nextGen;
		// loop for each population member
		for (int n = 1; n < this->popMember; n++)
		{

			for (int p = 0; p < this->population.size(); p++)
			{
				// get distances
				start = b_start;
				this->population[p].dist = this->get_distance(start, targets, this->population[p].dna);
				if (this->population[p].dist < this->minVal) { this->minVal = this->population[p].dist; }

				// get fitnesses
				this->population[p].fitness = this->get_fitness(this->population[p].dist);
			}

			// sort them w.r.t fitness
			sort(this->population.begin(), this->population.end(), [](const gPop& c1, const gPop& c2) {
				return c1.fitness > c2.fitness;
			});

			// get child from population
			gPop child;
			if (n == this->popMember - 1)
			{
				child = this->get_child(this->population[0], this->population[1], targets);
			}
			else
			{
				child = this->get_child(this->population[0], this->population[n], targets);
			}

			nextGen.push_back(child);
		}

		// make new population from child
		this->population = nextGen;

		// seach if each child has same shortest distance, stop seaching.
		auto s = std::adjacent_find(this->population.begin(), this->population.end(), [](const gPop& c1, const gPop& c2) {
			return c1.dist < c2.dist;
		});

		if (s == this->population.end())
		{
			notsame = false;
		}

	}

	cout << " GENERATION : " << g << endl << endl;

	if (this->population[0].dist > this->minVal)
	{
		cout << endl << " !!! Probably this is not the best path. You can try again." <<
			"If you want, you can set genetic algorithm parameters with `set_genetic_params` method." << endl << endl;
	}

	this->cost = this->population[0].dist;


	// Finding path to targets with A*. Target orders find by
	// genetic algo. Now, connect these point with graph by A* algo.

	// add start point to target at the end to return start city.
	this->population[0].dna.push_back(targets.size());
	targets.push_back(this->b_start);
	this->sub_dest.push_back(b_start);

	for (int g = 0; g < this->population[0].dna.size(); g++)
	{
		int tr = targets[this->population[0].dna[g]];
		int a = 0;
		vector<int> subDest = this->A_star(b_start, tr);

		for (int j = 0; j < subDest.size(); j++)
		{
			this->sub_dest.push_back(subDest[j]);
			cout <<this->city[subDest[j]] << " (" << this->city_plate[subDest[j]] <<
				")" << " >> ";
		}
		cout << this->city[tr] << " (" << this->city_plate[tr] <<
			")" << " || FROM " << this->city[b_start]
			<< " to " << this->city[tr] << " => KM : "<< this->distance[b_start][tr] ;
		b_start = tr;

		cout << endl;
	}

	this->sub_dest.push_back(b_start);

	cout << endl << "Total KM : " << cost << endl;
}

// A* implementation
vector<int> graph::A_star(int s_Ind, int d_Int)
{
	// initialize some params
	vector<int> A_dest;
	float g = 0, h = 0;
	int current = s_Ind;

	// loop
	while(current != d_Int)
	{
		// getting current node
		vector<gTarget> sGraph = this->adjList[current];

		// initialize some params
		float f = 9999999.;
		float g_pre = g;
		int tmp_ind=current;

		//cout << current << " " << endl;

		// search each neighboor nodes
		for (int i = 0; i < sGraph.size(); i++)
		{
			float temp = sGraph[i].second;
			h = this->distance[sGraph[i].first][d_Int];

			float f_temp = g_pre + temp + h;

			// if the node is shortest up to now, get properties
			if (f_temp < f)
			{
				g = g_pre + temp;
				tmp_ind = sGraph[i].first;
				f = f_temp;
			}
		}
		A_dest.push_back(current);
		current = tmp_ind;
	}

	return A_dest;
}


// print the map with path.
// scale mean scaling of map.
void graph::print_map(float scale, bool save)
{
	cv::Mat image;
	//load image
	image = cv::imread("data/map-tr.png");

	// print sub destinations on image
	for (int i = 1; i < this->sub_dest.size(); i++)
	{
		int t = this->sub_dest[i];
		int t1 = this->sub_dest[i-1];
		cv::Point2i point2(this->x_loc[t], this->y_loc[t]);
		cv::Point2i point1(this->x_loc[t1], this->y_loc[t1]);
		cv::arrowedLine(image, point1, point2, cv::Scalar(0, 0, 0),3);
	}

	// get size of image
	cv::Size im_size = image.size();

	// resize image with scale
	cv::resize(image, image, cv::Size(im_size.width * scale, im_size.height * scale));

	// if willing to save, save it
	cv::imwrite("data/path.png", image);

	// show image
	cv::imshow("Path", image);
	cv::waitKey(0);
}


// setting genetic algorithm parameters.
// popMember : population member number
// pCross : probability of crossover
// pMut : probability of mutation
void graph::set_genetic_params(int popMember, float pCross, float pMut)
{
	this->popMember = popMember;
	this->pCross = pCross;
	this->pMut = pMut;
}


// generate population
void graph::gen_pop(int target_size)
{
	// create population size
	for (int p = 0; p < this->popMember; p++)
	{
		// create range for target size then shuffle it
		vector<int> r = range(0, target_size);
		random_shuffle(r.begin(), r.end());
		gPop t;
		t.dna = r;
		t.dist = 0;
		t.fitness = 0.;
		this->population.push_back(t);
	}
}


// getting fitness
float graph::get_fitness(int dist)
{
	return 1./dist;
}

// getting child w.r.t parents p1, p2
gPop graph::get_child(gPop p1, gPop p2, vector<int> targets)
{
	gPop child;
	// seed to get random number
	srand(pow(static_cast<unsigned int>(std::time(nullptr)), 2));

	// get probabily to crossover
	float crossProb = (double)rand() / (RAND_MAX);

	// condition on crossover
	if (crossProb < this->pCross)
	{
		child.dna = this->crossover(p1, p2);
		child.dist = this->get_distance(this->b_start, targets, child.dna);
		child.fitness = this->get_fitness(child.dist);
	}
	else
	{
		child = p1;
	}

	// get probability to mutation
	float mutateProb = (double)rand() / (RAND_MAX);

	// condition on mutation
	if (mutateProb < this->pMut)
	{
		child.dna = this->mutation(child);
	}

	return child;
}


// crossovering
vector<int> graph::crossover(gPop p1, gPop p2)
{
	// getting index of parent genes with randomly
	int n1 = floor((double)rand() / (RAND_MAX)* p1.dna.size());
	int n2 = floor((double)rand() / (RAND_MAX)* p1.dna.size());

	int start = min(n1, n2);
	int end = max(n1, n2);

	// get p1's genes to child
	vector<int> child_pop = vector<int>(p1.dna.begin() + start, p1.dna.begin() + end);

	// check if genes is not from p2, add it.
	for (auto index : p2.dna)
	{
		vector<int>::iterator it = std::find(child_pop.begin(), child_pop.end(), index);
		if (it == child_pop.end())
		{
			child_pop.push_back(index);
		}

	}

	return child_pop;
}

// mutating the genes
vector<int> graph::mutation(gPop child)
{
	// getting index of child genes randomly
	int n1 = floor((double)rand() / (RAND_MAX)* child.dna.size());
	int n2 = floor((double)rand() / (RAND_MAX)* child.dna.size());

	int start = min(n1, n2);
	int end = max(n1, n2);

	vector<int> newGenesOrder;

	// mutate the genes.
	vector<int> mutate_genes = vector<int>(child.dna.begin() + start, child.dna.begin() + end);

	// shuffe the mutation genes
	random_shuffle(mutate_genes.begin(), mutate_genes.end());

	// check dna order
	bool addMutate = true;
	for (auto index : child.dna)
	{
		vector<int>::iterator it = std::find(mutate_genes.begin(), mutate_genes.end(), index);

		// if gen is not in mutated genes, add with same order.
		// else add mutated genes
		if (it == mutate_genes.end())
		{
			newGenesOrder.push_back(index);
		}
		else
		{
			if (addMutate)
			{
				for (auto m_genes : mutate_genes)
				{
					newGenesOrder.push_back(m_genes);
				}
				// added mutated genes.
				addMutate = false;
			}
		}
	}
	return mutate_genes;
}

// get distance of cities for population.
int graph::get_distance(int start, vector<int> targets, vector<int> pop)
{
	float tmp = 0;
	for (int t = 0; t < pop.size(); t++)
	{
		tmp += this->distance[start][targets[pop[t]]];
		start = targets[pop[t]];
	}
	// adding start location to distance
	tmp += this->distance[start][this->b_start];
	return tmp;
}