#include <iostream>
#include <iomanip>
#include <vector>
#include <time.h>
#include <chrono>
#include <algorithm>
#include <random>
#include <fstream>
#include <map>

using namespace std;
using namespace std::chrono;
typedef pair<int, int> Pair;

string DisplayPath(const vector<int> &vec)
{
    string temp = "";
    for (int i = 0; i < vec.size(); i++)
    {
        temp += to_string(vec[i]) + ' ';
        cout << vec[i] << " ";
    }
    cout << endl;
    return temp;
}

struct Node
{
    int data;
    bool visited;
    vector<pair<Node *, int>> edges;
    Node()
    {
        visited = false;
        data = 0;
    }
    Node(int d)
    {
        data = d;
    }
};

struct Graph
{
    vector<Node *> adjList;
    Node *addVertex(int s)
    {
        for (int i = 0; i < adjList.size(); i++)
        {
            if (adjList[i]->data == s)
            {
                cout << "Vertex already exists" << endl;
                return NULL;
            }
        }
        Node *temp = new Node(s);
        adjList.push_back(temp);
        return temp;
    }

    int numNodes()
    {
        return adjList.size();
    }

    void addEdge(int s, int d, int w)
    {
        bool s1 = false, s2 = false;
        Node *source;
        Node *dest;
        for (int i = 0; i < adjList.size(); i++)
        {
            if (adjList[i]->data == s)
            {
                s1 = true;
                source = adjList[i];
            }
            if (adjList[i]->data == d)
            {
                s2 = true;
                dest = adjList[i];
            }
        }
        if (!s1)
        {
            source = addVertex(s);
        }
        if (!s2)
        {
            dest = addVertex(d);
        }
        else
        {
            for (int i = 0; i < source->edges.size(); i++)
            {
                int e = source->edges[i].first->data;
                if (e == d)
                {
                    cout << "( " << s << " , " << d << " )Edge Already Exists" << endl;
                    return;
                }
            }
        }
        source->edges.push_back(make_pair(dest, w));
        dest->edges.push_back(make_pair(source, w));
    }

    vector<int> GeneratePath(int s, int d)
    {
        Node *source, *temp;
        for (int i = 0; i < adjList.size(); i++)
        {
            if (adjList[i]->data == s)
            {
                source = adjList[i];
                break;
            }
        }
        temp = source;
        int random;
        vector<int> Arr;
        Arr.push_back(s);
        int i = 1;
        do
        {
            //DisplayPath(Arr);
            if(temp->edges.size() == 0){
                // cout << temp->edges.size() << endl;
                vector<int> a;
                return a;
            }
            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<> distrib(0, temp->edges.size() - 1);
            random = d;
            bool check = true;
            int k = distrib(gen);
            for (; k < temp->edges.size(); k++)
            {
                random = k;
                check = false;
                for (int j = 0; j < i; j++)
                {
                    if (temp->edges[random].first->data == Arr[j])
                    {
                        check = true;
                        break;
                    }
                }
                if (check == false)
                {
                    break;
                }
                else if (k == temp->edges.size() - 1 && check == true)
                {
                    //DisplayPath(Arr);
                    vector<int> a;
                    return a;
                }
            }

            Arr.push_back(temp->edges[random].first->data);
            temp = temp->edges[random].first;
            i++;

        } while (Arr[i - 1] != d);
        //DisplayPath(Arr);
        return Arr;
    }

    vector<vector<int>> generateNPaths(int numofPath, int s, int d)
    {
        vector<vector<int>> paths;
        bool check = true;
        for (int i = 0; i < numofPath; i++)
        {
            cout << "Path Finding: " << i+1 << endl;
            check = true;
            while (check)
            {
                check = false;
                vector<int> temp = GeneratePath(s, d);
                if (temp.empty())
                {
                    check = true;
                    // cout << "A Cycle Detected" << endl;
                    continue;
                }
                if (paths.empty())
                {
                    paths.push_back(temp);
                    break;
                }
                else
                {
                    check = false;
                }
                if (check == false)
                {
                    paths.push_back(temp);
                    cout << "Paths found: " << paths.size() << endl;
                }
            }
        }
        return paths;
    }

    int findCost(const vector<int> &graphPath)
    {
        int s = graphPath[0];
        int cost = 0;
        Node *source;
        for (int i = 0; i < adjList.size(); i++)
        {
            if (graphPath[0] == adjList[i]->data)
            {
                source = adjList[i];
                break;
            }
        }
        for (int i = 1; i < graphPath.size(); i++)
        {
            int next = graphPath[i];
            for (int j = 0; j < source->edges.size(); j++)
            {
                if (next == source->edges[j].first->data)
                {
                    cost += source->edges[j].second;
                    source = source->edges[j].first;
                    break;
                }
            }
        }
        return cost;
    }
    void LoadData(string filename)
    {
        string line;
        ifstream myfile(filename);
        if (myfile.is_open())
        {
            while (getline(myfile, line))
            {
                int s, d, w;
                int count = 0;
                string temp = "";
                for (int i = 0; i < line.length(); i++)
                {
                    if (line[i] != ',')
                    {
                        temp += line[i];
                    }
                    else
                    {
                        if (count == 0)
                        {
                            s = stoi(temp);
                            count++;
                        }
                        else
                        {
                            d = stoi(temp);
                        }
                        temp = "";
                    }
                }
                w = stoi(temp);
                addEdge(s, d, w);
            }
            myfile.close();
            cout << " Data Loaded Successfully! " << endl;
        }

        else
            cout << "Unable to open file";
    }
    void Display()
    {
        for (int i = 0; i < adjList.size(); i++)
        {
            for (int j = 0; j < adjList[i]->edges.size(); j++)
            {
                int data = adjList[i]->edges[j].first->data;
                int weight = adjList[i]->edges[j].second;
                cout << "( " << adjList[i]->data << " , " << data << " , " << weight << " )" << endl;
            }
        }
    }
};

int Random(int max)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distrib(0, max - 1);
    return distrib(gen);
}

struct Path
{
    vector<int> route;
    int cost;
    double fitness;
    Path()
    {
        fitness = 0;
        cost = 0;
    }
    Path(vector<int> a, int c, double f)
    {
        route = a;
        cost = c;
        fitness = f;
    }
    void PathRoute(){
        DisplayPath(route);
    }
    bool operator==(const Path &n1)
    {
        if (n1.route == route)
        {
            return true;
        }
        return false;
    }
};

class ShortestPath
{
    vector<Path> Population;
    int PopulationSize;
    int Generation;
    double FitnessFunction(int cost)
    {
        return (1.0f / (double)cost);
    }
    void PopulationInitialization(Graph g, int s, int d)
    {
        //Generate Random Paths
        vector<vector<int>> paths = g.generateNPaths(PopulationSize, s, d);
        
        //Set the Population
        for (int i = 0; i < paths.size(); i++)
        {
            int c = g.findCost(paths[i]);
            double fit = FitnessFunction(c);
            Population.push_back({paths[i], c, fit});
        }
    }
    void DisplayPopulation()
    {
        for (int i = 0; i < Population.size(); i++)
        {
            cout << "Route: ";
            DisplayPath(Population[i].route);
            cout << "Cost: " << Population[i].cost << endl;
            cout << "Fitness: " << Population[i].fitness << endl;
        }

        cout << endl
             << endl;
    }

    void RepairFunction(vector<int> &vec)
    {
        //Write Repair
        map<int, int> count;
        int first = 0, sec = 0;
        bool cycle = false;
        for (int i = 0; i < vec.size(); i++)
        {
            if (count.find(vec[i]) == count.end())
            {
                count[vec[i]] = i;
            }
            else
            {
                cycle = true;
                first = count[vec[i]];
                sec = i;
                count[vec[i]]++;
            }
        }
        if (cycle)
            vec.erase(vec.begin() + first, vec.begin() + sec);
    }

    Path
    TournamentSelection()
    {
        int rd1 = Random(PopulationSize);
        Path R1 = Population[rd1];
        rd1 = Random(PopulationSize);
        Path R2 = Population[rd1];
        while (R1 == R2)
        {
            rd1 = Random(PopulationSize);
            R2 = Population[rd1];
        }
        if (R1.fitness > R2.fitness)
        {
            return R1;
        }
        
        return R2;
    }
    void CrossoverFunction(Graph g)
    {
        bool Converged = false;
        vector<vector<int>> crossOvers;
        for (int i = 0; i < Population.size(); i++)
        {
            //---------------------------------------------- Display
            // cout << "Tournament " << i + 1 << endl;
            // Tournament Selection
            Path Mother = TournamentSelection();
            Path Father = TournamentSelection();
            while (Mother == Father)
            {
                Father = TournamentSelection();
                // cout << "Father: ";
                // Father.PathRoute();
                // cout << "Mother: ";
                // Mother.PathRoute();

            }
            // if(Converged){
            //     break;
            // }
            //---------------------------------------------- Display
            // cout << "Mother: ";
            // DisplayPath(Mother.route);
            // cout << "Father: ";
            // DisplayPath(Father.route);
            
            // Perform Crossover
            vector<Pair> s;
            for (int j = 1; j < Mother.route.size() - 1; j++)
            {
                for (int k = 1; k < Father.route.size() - 1; k++)
                {
                    if (Mother.route[j] == Father.route[k])
                    {
                        s.push_back(make_pair(j, k));
                    }
                }
            }
            // ----------------------------------------------------------- Display
            // for (int j = 0; j < s.size(); j++)
            // {
            //     cout << "Crossing Site " << j + 1 << ": " << s[j].first << ", " << s[j].second << endl;
            // }
            //Perfrom Exchange between Mother and Father
            if (s.size() > 0)
            {
                random_device rd;
                mt19937 gen(rd());
                uniform_int_distribution<> distrib(0, s.size() - 1);
                int chosen = distrib(gen);
                //------------------------------------------------------------ Display
                // cout << "Chosen Site: " << chosen + 1 << endl;
                Pair site = s[chosen];
                vector<int> Offspring[2];
                for (int j = 0; j < site.first; j++)
                {
                    Offspring[0].push_back(Mother.route[j]);
                }
                for (int j = site.second; j < Father.route.size(); j++)
                {
                    Offspring[0].push_back(Father.route[j]);
                }
                for (int j = 0; j < site.second; j++)
                {
                    Offspring[1].push_back(Father.route[j]);
                }
                for (int j = site.first; j < Mother.route.size(); j++)
                {
                    Offspring[1].push_back(Mother.route[j]);
                }

                for (int j = 0; j < 2; j++)
                {
                    // Repair Offspring[j]
                    // Implement Repair Function
                    // Hint: Send Offspring[j] as a Reference
                    RepairFunction(Offspring[j]);
                    //--------------------------------------------------------- Display
                    // cout << "Offspring" << j + 1 << ": ";
                    // DisplayPath(Offspring[j]);
                    crossOvers.push_back(Offspring[j]);
                }
            }
            //---------------------------------------------- Display
            // cout << endl;
        }
        // Add the New Offsprings to the Existing Population
        for (int i = 0; i < crossOvers.size(); i++)
        {
            int c = g.findCost(crossOvers[i]);
            double fit = FitnessFunction(c);
            Population.push_back({crossOvers[i], c, fit});
        }
    }
    // Help Function to Compare Fitness in the Culling Stage
    static bool CompareFitness(Path p1, Path p2)
    {
        return (p1.fitness > p2.fitness);
    }
    // Keep the Fittest individuals
    void Culling(int n)
    {
        sort(Population.begin(), Population.end(), CompareFitness);
        Population.erase(Population.begin() + n, Population.end());
    }

    void MutationFunction(Graph &g, int i, int d)
    {
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> distrib(1, Population[i].route.size() - 2);
        // Generating a random index
        int random = distrib(gen);
        // Point to mutate from
        int s = Population[i].route[random];
        //---------------------------------------------- Display
        // cout << "Mutation Point: " << s << endl;
        vector<int> newPath = g.GeneratePath(s, d);
        while (newPath.empty())
        {
            newPath = g.GeneratePath(s, d);
        }
        Population[i].route.erase(Population[i].route.begin() + random, Population[i].route.end());
        for (int j = 0; j < newPath.size(); j++)
        {
            Population[i].route.push_back(newPath[j]);
        }
        // Repair Mutation
        RepairFunction(Population[i].route);
        int c = g.findCost(Population[i].route);
        Population[i].cost = c;
        Population[i].fitness = FitnessFunction(c);

        //---------------------------------------------- Display
        // cout << "Mutation Performed: ";
        // DisplayPath(Population[i].route);
    }

    float ConvergenceFunction()
    {
        int max = 0;
        map<double, int> count;
        for (int i = 0; i < Population.size(); i++)
        {
            if (count.find(Population[i].fitness) == count.end())
            {
                count[Population[i].fitness] = 1;
            }
            else
            {
                count[Population[i].fitness]++;
            }
        }
        for (auto i : count)
        {
            if (i.second > max)
                max = i.second;
        }

        return (float)max / (float)PopulationSize;
    }

    void OptimalPath(float dur)
    {
        ofstream myfile("LargeOutputs.txt", ios::app);
        if (myfile.is_open())
        {
            myfile << "--------------------- OUTPUT ---------------------" << endl;
            myfile << "Generation: " << Generation << endl;
            myfile << "Source: " << Population[0].route[0] << endl;
            myfile << "Destination: " << Population[0].route[Population[0].route.size() - 1] << endl;
            myfile << "Optimal Path: " << DisplayPath(Population[0].route) << endl;
            myfile << "Cost: " << Population[0].cost << endl;
            myfile << "Fitness: " << Population[0].fitness << endl;
            myfile << "Population Size: " << PopulationSize << endl;
            myfile << "Time taken by function: " << dur << " seconds" << endl;
            myfile << "----------------------- END ---------------------" << endl<<endl;
            myfile.close();
        }
        cout << "\nGeneration: " << Generation << endl;
        cout << "Source: " << Population[0].route[0] << endl;
        cout << "Destination: " << Population[0].route[Population[0].route.size()-1] << endl;
        cout << "Optimal Path: ";
        DisplayPath(Population[0].route);
        cout << "Cost: " << Population[0].cost << endl;
        cout << "Fitness: " << Population[0].fitness << endl;
        cout << "Population Size: " << PopulationSize<< endl;
        cout << "Time taken by function: " << dur << " seconds" << endl;
    }

public:
    void
    FindShortestPath(Graph g, int s, int d, float m, int p)
    {
        
        PopulationSize = p;
        Generation = 0;
        PopulationInitialization(g, s, d);
        // Get starting timepoint
        auto start = high_resolution_clock::now();
        bool Converged = false;
        while (!Converged)
        {
            cout << "Generation " << Generation << endl;
            // DisplayPopulation();
            CrossoverFunction(g);
            Culling(PopulationSize);
            // DisplayPopulation();
            for (int i = 0; i < PopulationSize; i++)
            {
                float r = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                if (r < m)
                {
                    // Implement this
                    MutationFunction(g, i, d);
                }
            }

            if (ConvergenceFunction() >= 0.5)
            {
                auto stop = high_resolution_clock::now();
                auto duration = duration_cast<microseconds>(stop - start);
                Converged = true;
                float dur = (float)duration.count() / (float)1000000;
                OptimalPath(dur);
            }
            Generation++;
        }
    }
};

int main()
{
    srand(time(NULL));
    string f = "Original.csv";
    int source = 1;
    int destination = 20;
    float mutationProb = 0.05;
  
    cout << "Enter File Name: ";
    cin >> f;
    cout << "Enter Source: ";
    cin >> source;
    cout << "Enter Destination: ";
    cin >> destination;
    cout << "Enter Muatation Prob: ";
    cin >> mutationProb;
  
    Graph graph;
    graph.LoadData(f);
    int psize = graph.numNodes();
    cout << "Number of Nodes in Graph: " << graph.numNodes() << endl;
    ShortestPath a;;
    a.FindShortestPath(graph, source, destination, mutationProb, psize);
}
