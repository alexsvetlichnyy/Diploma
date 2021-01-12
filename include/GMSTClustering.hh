#ifndef GMSTClustering_h
#define GMSTClustering_h 1

#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <vector>
#include "../TGlauber/TGlauberMC.hh"
#include "../TGlauber/TGlauNucleon.hh"
#include "TObject.h"
#include "TVector3.h"
using namespace std;

// Creating shortcut for an integer pair 
typedef  pair<G4int, G4int> iPair;

// Structure to represent a graph 
struct Graph
{
	// Vert and edges
    G4int V, E;
    vector< pair<G4double, iPair> > edges;

    // Pointer to an array containing adjacency lists 
    list<G4int>* adj;

    // Constructors 
    Graph(G4int V, G4int E);
    Graph();

    // Destructor
	~Graph();

    // Utility function to add an edge 
    void addEdge(G4int u, G4int v, G4double w);
    void DFSUtil(G4int v, G4bool visited[], vector < G4int > *clust_inside);
    // method to add a connection
    void addConn(G4int v, G4int w);
	// Method to print connected components in an 
	// undirected graph
    vector< vector <G4int> > connectedComponents();

    // Function to find MST using Kruskal's 
    // MST algorithm 
    vector< vector <G4int> > AdvancedKruskalMST(G4double CD_in);
};

// To represent Disjoint Sets
struct DisjointSets
{
    G4int* parent, * rnk;
    G4int n;

    // Constructor
    DisjointSets(G4int n);

    // Destructor
    ~DisjointSets();

    // Find the parent of a node 'u' 
    // Path Compression 
    G4int find(G4int u);

    // Union by rank 
    void merge(G4int x, G4int y);
};

class GMSTCluster{

	public: 
	GMSTCluster(G4int Z_in, G4int A_in);
	~GMSTCluster();

	public: 
	inline G4int GetZ() {return Z;};
	inline G4int GetA() {return A;};
	inline void SetZ(G4int Z_in) {Z = Z_in;}
	inline void SetA(G4int A_in) {A = A_in;}
	inline void PushBackCoordinateVector(TVector3 vec_in) {coord.push_back(vec_in);}
        inline std::vector<TVector3> GetCoordinates() {return coord;}

	private: 
	G4int Z;
	G4int A;
	std::vector<TVector3> coord;
};

typedef std::vector<GMSTCluster> GMSTClusterVector;

class GMSTClustering{

	public:
	GMSTClustering();
	GMSTClustering(G4double CD_in);
	~GMSTClustering();

	public:
	inline G4double SetCD(G4double CD_in) {CritDist = CD_in;}
	Graph ClusterToGraph(TObjArray* nucleons, G4double A);

	void GetClusters(TObjArray* nucleons, GMSTClusterVector& output_vector, GMSTClusterVector& output_vector_B);
	private:
	G4double CritDist;
};


#endif
