#include "../include/GMSTClustering.hh"

GMSTCluster::GMSTCluster(G4int Z_in, G4int A_in): 
A(A_in), Z(Z_in)
{
};

GMSTCluster::~GMSTCluster(){

};

GMSTClustering::GMSTClustering(){
CritDist = 0;
};

GMSTClustering::GMSTClustering(G4double CD_in){
CritDist = CD_in;
};

GMSTClustering::~GMSTClustering(){

};

Graph GMSTClustering::ClusterToGraph(TObjArray* nucleons, G4double A){
    Graph g(A, A*(A-1)/2);
    //  making full graph of nucleons
    for(G4int iArray = 0; iArray < nucleons->GetEntries(); iArray++){
    	TGlauNucleon *nucleon=(TGlauNucleon*)(nucleons->At(iArray));
    	for(G4int iArray_pairs = iArray + 1; iArray_pairs < nucleons->GetEntries(); iArray_pairs++){
    		TGlauNucleon *nucleon_pair=(TGlauNucleon*)(nucleons->At(iArray_pairs));
    		g.addEdge(iArray, iArray_pairs, std::sqrt(pow(nucleon->GetX() - nucleon_pair->GetX(),2) + pow(nucleon->GetY() - nucleon_pair->GetY(),2) + pow(nucleon->GetZ() - nucleon_pair->GetZ(),2)));
    	}
	}
    return g;
};

void GMSTClustering::GetClusters(TObjArray* nucleons_in, GMSTClusterVector& output_vector, GMSTClusterVector& output_vector_B){
	G4int A = 0;
	G4int Z = 0;
	G4int Ab = 0;
	G4int Zb = 0;
	TObjArray *nucleons = new TObjArray(1000); // nucleons from Side A
	TObjArray *nucleons_B = new TObjArray(1000); // nucleons from Side B
	for(G4int iArray = 0; iArray < nucleons_in->GetEntries(); iArray++){
		TGlauNucleon *nucleon=(TGlauNucleon*)(nucleons_in->At(iArray));
		if(nucleon->IsSpectator() && nucleon->IsInNucleusA()){
			A+=1;
			nucleons->AddLast((TObject*)(nucleons_in->At(iArray)));
		}
		if(nucleon->IsSpectator() && nucleon->IsInNucleusA() && nucleon->IsProton()){Z+=1;} 
		if(nucleon->IsSpectator() && nucleon->IsInNucleusB()){
			Ab+=1;
			nucleons_B->AddLast((TObject*)(nucleons_in->At(iArray)));
		}
		if(nucleon->IsSpectator() && nucleon->IsInNucleusB() && nucleon->IsProton()){Zb+=1;}
	}

	// Creating a graph of nucleons (side A)
	Graph g = this->ClusterToGraph(nucleons, A);

	// Creating a graph of nucleons (side B)
	Graph g_B = this->ClusterToGraph(nucleons_B, Ab);

	// Applying MST + critical distance cut + DFS clustering (side A)
	vector< vector <G4int> > clusters;
    clusters = g.AdvancedKruskalMST(this->CritDist);

    // Applying MST + critical distance cut + DFS clustering (side B)
	vector< vector <G4int> > clusters_B;
    clusters_B = g_B.AdvancedKruskalMST(this->CritDist);

    // Filling a Clister Vector (side A)
    for(G4int i = 0; i < clusters.size(); ++i) {
    	G4int Z_clust = 0;
    	G4int A_clust = 0;
	GMSTCluster outClust(0,0);
    	for(G4int j = 0; j < clusters[i].size(); ++j) {
        	TGlauNucleon *nucleon=(TGlauNucleon*)(nucleons->At((clusters[i])[j]));
        	if(nucleon->IsProton())
        	{
        		Z_clust += 1;
        	}
        	A_clust += 1;
		outClust.PushBackCoordinateVector(TVector3(nucleon->GetX(), nucleon->GetY(), nucleon->GetZ()));
        }
	outClust.SetZ(Z_clust);
	outClust.SetA(A_clust);
	output_vector.push_back(outClust);
        //output_vector.push_back(GMSTCluster(Z_clust, A_clust));
    }

    // Filling a Clister Vector (side B)
    for(G4int i = 0; i < clusters_B.size(); ++i) {
    	G4int Z_clust = 0;
    	G4int A_clust = 0;
	GMSTCluster outClust_B(0,0);
    	for(G4int j = 0; j < clusters_B[i].size(); ++j) {
        	TGlauNucleon *nucleon=(TGlauNucleon*)(nucleons_B->At((clusters_B[i])[j]));
                if(nucleon->IsProton())
        	{
        		Z_clust += 1;
        	}
        	A_clust += 1;
		outClust_B.PushBackCoordinateVector(TVector3(nucleon->GetX(), nucleon->GetY(), nucleon->GetZ()));
        }
	outClust_B.SetZ(Z_clust);
	outClust_B.SetA(A_clust);
        //output_vector_B.push_back(GMSTCluster(Z_clust, A_clust));
        output_vector_B.push_back(outClust_B);
    }

    delete nucleons;
    delete nucleons_B;
};

Graph::Graph(G4int V, G4int E)
{
    this->V = V;
    this->E = E;
    adj = new list<G4int>[V];
}

Graph::Graph()
{
    this->V = 0;
    this->E = 0;
}

Graph::~Graph()
{
    delete[] adj;
}

void Graph::addEdge(G4int u, G4int v, G4double w)
{
    edges.push_back({ w, {u, v} });
}

vector< vector <G4int> > Graph::AdvancedKruskalMST(G4double CD_in)
{
    // Sort edges in increasing order on basis of cost 
    sort(edges.begin(), edges.end());

    // Create disjoint sets 
    DisjointSets ds(V);

    // Iterate through all sorted edges 
    vector< pair<G4double, iPair> >::iterator it;
    for (it = edges.begin(); it != edges.end(); it++)
    {
        G4int u = it->second.first;
        G4int v = it->second.second;

        G4int set_u = ds.find(u);
        G4int set_v = ds.find(v);

        // Check if the selected edge is creating 
        // a cycle or not (Cycle is created if u 
        // and v belong to same set) 
        if (set_u != set_v)
        {
            // Current edge will be in the MST 
            if(it->first < CD_in){
            	this->addConn(u, v);
            }
            // Merge two sets 
            ds.merge(set_u, set_v);
        }
    }

    return this->connectedComponents();
}

void Graph::addConn(G4int v, G4int w)
{
    adj[v].push_back(w);
    adj[w].push_back(v);
}

vector< vector <G4int> > Graph::connectedComponents()
{
	vector< vector <G4int> > clusters;
    // Mark all the vertices as not visited 
    G4bool* visited = new G4bool[V];
    for (G4int v = 0; v < V; v++)
        visited[v] = false;
    for (G4int v = 0; v < V; v++)
    {
        if (visited[v] == false)
        {
        	// vector to fill with nucleon's numbers from one cluster
        	vector < G4int > clust_inside; 
            // print all reachable vertices from v 
            DFSUtil(v, visited, &clust_inside);
            // push cluster to vector of clusters
            clusters.push_back(clust_inside); 
        }
    }
    delete[] visited;
    return clusters;
}

void Graph::DFSUtil(G4int v, G4bool visited[], vector < G4int > *clust_inside)
{
    // Mark the current node as visited and print it 
    visited[v] = true;
    (*clust_inside).push_back(v); 

    // Recur for all the vertices adjacent to this vertex
    list<G4int>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
        if (!visited[*i])
            DFSUtil(*i, visited, &(*clust_inside));
}

DisjointSets::DisjointSets(G4int n)
{
    // Allocate memory 
    this->n = n;
    parent = new G4int[n + 1];
    rnk = new G4int[n + 1];

    // Initially, all vertices are in different sets and have rank 0.
    for (G4int i = 0; i <= n; i++)
    {
        rnk[i] = 0;
        //every element is parent of itself 
        parent[i] = i;
    }
}

DisjointSets::~DisjointSets()
{
	delete[] parent;
	delete[] rnk;
}

G4int DisjointSets::find(G4int u)
{
    /* Make the parent of the nodes in the path
       from u--> parent[u] point to parent[u] */
    if (u != parent[u])
        parent[u] = find(parent[u]);
    return parent[u];
}

void DisjointSets::merge(G4int x, G4int y)
{
    x = find(x), y = find(y);

    /* Make tree with smaller height
       a subtree of the other tree  */
    if (rnk[x] > rnk[y])
        parent[y] = x;
    else // If rnk[x] <= rnk[y] 
        parent[x] = y;

    if (rnk[x] == rnk[y])
        rnk[y]++;
}
