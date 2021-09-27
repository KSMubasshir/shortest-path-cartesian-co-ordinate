#include<bits/stdc++.h>
#include <cmath>
#define N 100
#define PI 3.14159265359
using namespace std;
//============disjoint set union==============
// a structure to represent an edge in graph
struct Edge{
    int src, dest;
};
// a structure to represent a graph
struct Graph{
    // V-> Number of vertices, E-> Number of edges
    int V, E;
    // graph is represented as an array of edges
    struct Edge* edge;
};
// Creates a graph with V vertices and E edges
struct Graph* createGraph(int V, int E){
    struct Graph* graph = (struct Graph*) malloc(sizeof(struct Graph));
    graph->V = V;
    graph->E = E;
    graph->edge = (struct Edge*) malloc( graph->E * sizeof(struct Edge));
    return graph;
}
// A utility function to find the subset of an element i
int find(int parent[], int i){
    if (parent[i] == -1)
        return i;
    return find(parent, parent[i]);
}
// A utility function to do union of two subsets
void Union(int parent[], int x, int y){
    int xset = find(parent, x);
    int yset = find(parent, y);
    if(xset!=yset){
        parent[xset] = yset;
    }
}
// The main function to check whether a given graph contains
// cycle or not
bool isCycle( struct Graph* graph ){
    // Allocate memory for creating V subsets
    int *parent = (int*) malloc( graph->V * sizeof(int) );
    // Initialize all subsets as single element sets
    memset(parent, -1, sizeof(int) * graph->V);
    // Iterate through all edges of graph, find subset of both
    // vertices of every edge, if both subsets are same, then
    // there is cycle in graph.
    for(int i = 0; i < graph->E; ++i){
        int x = find(parent, graph->edge[i].src);
        int y = find(parent, graph->edge[i].dest);
        if (x == y)
            return true;
        Union(parent, x, y);
    }
    return false;
}
//===========disjoint set ends here============
double costGeo(double lat1, double lon1, double lat2, double lon2);
double deg2rad(double deg);
bool sortbysec(const pair<int,double> &a, const pair<int,double> &b){
    return (a.second > b.second);
}
bool sortbysavings(const  pair < double,pair< int,int > >  &a, const pair < double,pair< int,int > > &b){
    return (a.first < b.first);
}
double costGeo(double lat1, double lon1, double lat2, double lon2){
    double R = 6371; // Radius of the earth in km
    double dLat = deg2rad(lat2 - lat1);
    double dLon = deg2rad(lon2 - lon1);
    double a = sin(dLat / 2) * sin(dLat / 2)
               + cos(deg2rad(lat1)) * cos(deg2rad(lat2))
               * sin(dLon / 2) * sin(dLon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    double d = R * c; // Distance in km
    return d;
}
double deg2rad(double deg){
    return deg * (PI / 180);
}
struct coord{
    double lat;
    double lon;
    coord() {}
    coord(double a,double b){
        lat=a;
        lon=b;
    }
    double distance(coord A){
        return costGeo(lat,lon,A.lat,A.lon);
    }
};
int n;
coord points[N];
bool visited[N];
vector<int>Tour;
int nearestUnvisted(int x){
    double Min=LLONG_MAX;
    int Idx=-1;
    vector< pair<int,int> > arr;
    for(int i=0; i<n; i++){
        if(visited[i]) continue;
        if(points[x].distance(points[i])<Min){
            Min=points[x].distance(points[i]);
            Idx=i;
        }
    }
    return Idx;
}
vector<int> nearest5Unvisted(int x){
    double Min=LLONG_MAX;
    int Idx=-1;
    vector< pair<int,double> > arr;
    for(int i=0; i<n; i++){
        if(visited[i]) continue;
        arr.push_back(make_pair(i,points[x].distance(points[i])));
    }
    sort(arr.begin(), arr.end(), sortbysec);
    vector<int> ret;
    int sizeOfVect=arr.size();
    int mySize=min(sizeOfVect,5);
    for(int i=0; i<mySize; i++){
        ret.push_back(arr[i].first);
    }
    return ret;
}
int nearest(){
    double Min=LLONG_MAX;
    int Idx=-1;
    for(int i=0; i<n; i++){
        if(!visited[i]) continue;
        int Nearest=nearestUnvisted(i);
        if(points[i].distance(points[Nearest])<Min) Min=points[i].distance(points[Nearest]), Idx=Nearest;
    }
    return Idx;
}
int nearestEdge(int node){
    double Min=LLONG_MAX;
    int Idx=-1;
    int Size=Tour.size();
    for(int edge=0; edge<Size; edge++){
        int x=Tour[edge];
        int y=Tour[(edge+1)%Size];
        double Dist=points[x].distance(points[node])+points[y].distance(points[node])-points[x].distance(points[y]);
        if(Dist<Min) Min=Dist, Idx=x;
    }
    return Idx;
}
double Cost(vector<int>Tour){
    double Ans=0;
    int Size=Tour.size();
    for(int i=0; i<Size; i++){
        int xx=Tour[i];
        int yy=Tour[(i+1)%Size];
        Ans+=points[xx].distance(points[yy]);
    }
    return Ans;
}
void printSolution(vector<int>Tour){
    printf("Solution Cost : %lf\n",Cost(Tour));
    printf("Tour Description : ");
    for(int i=0; i<Tour.size(); i++) printf("%d ",Tour[i]+1);
    printf("\n\n");
}
vector<int> nearestNeighbourHeuristic(){
    Tour.clear();
    memset(visited,0,sizeof(visited));
    int Start=rand()%n;
    int Last=Start;
    visited[Last]=1;
    Tour.push_back(Start);
    while(Tour.size()<n){
        int Idx=nearestUnvisted(Last);
        visited[Idx]=1;
        Tour.push_back(Idx);
    }
    return Tour;
}
vector<int> nearestNeighbourHeuristic2(int Start){
    Tour.clear();
    memset(visited,0,sizeof(visited));
    int Last=Start;
    visited[Last]=1;
    Tour.push_back(Start);
    while(Tour.size()<n){
        int Idx=nearestUnvisted(Last);
        visited[Idx]=1;
        Tour.push_back(Idx);
    }
    return Tour;
}
vector<int> nearest5NeighbourHeuristic(int Start){
    Tour.clear();
    memset(visited,0,sizeof(visited));
    vector<int> temp;
    int Last=Start;
    visited[Last]=1;
    Tour.push_back(Start);
    while(Tour.size()<n){
        temp=nearest5Unvisted(Last);
        int Idx=temp[rand()%temp.size()];
        visited[Idx]=1;
        Tour.push_back(Idx);
    }
    return Tour;
}
double savings(int h,int i,int j){
    return  costGeo(points[h].lat,points[h].lon,points[i].lat,points[i].lon)
            +   costGeo(points[h].lat,points[h].lon,points[j].lat,points[j].lon)
            -   costGeo(points[i].lat,points[i].lon,points[j].lat,points[j].lon);
}
void findAndErase(vector<int> vis,int elem){
    vector<int>::iterator it;
    it = find (vis.begin(), vis.end(), elem);
    vis.erase(it);
}
bool createsCycle(int **arr,int a,int b){
    int V = Tour.size()+2, E = V-1;
    int edges=0;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            //if()
        }
    }
    struct Graph* graph = createGraph(V, E);
    vector<int> temp=Tour;
    //temp.push_back(i);
    //temp.push_back(j);
    for(int k=0;k<E;k++){
        graph->edge[temp[k]].src = temp[k];
        graph->edge[temp[k]].dest = temp[k+1];
    }
    return isCycle(graph);
}
bool allDegreesLessThan2(int **arr,int a,int b){
    if(arr[a][b]==2||arr[b][a]==2)
        return false;
    return true;
}
void savingsHeuristics(){
    int arr[n][n];
    memset(arr,0,sizeof(arr));
    memset(visited,0,sizeof(visited));
    vector< pair < double,pair< int,int > > > savingsPair;
    Tour.clear();
    vector<int> VH;
    for(int i=0; i<n; i++){
        VH[i]=i;
    }
    int h=rand()%n;
    findAndErase(VH,h);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i!=j && i!=h && j!=h){
                //arr[i][j]=savingsHeuristics(h,i,j);
                savingsPair.push_back(make_pair(savings(h,i,j),make_pair(i,j)));
            }
        }
    }
    sort(savingsPair.begin(),savingsPair.end(),sortbysavings);
    int i,j;
    pair< int,int > temp;
    while(VH.size()>2){
        temp=savingsPair[0].second;
        i=temp.first;
        j=temp.second;
        if(!createsCycle((int**)arr,i,j)){
                //take edge(i,j)
        }
    }
}
void twoOptHeuristic(){
    Tour.clear();
    memset(visited,0,sizeof(visited));
    nearestNeighbourHeuristic();
    while(true){
        double Curr=Cost(Tour);
        bool Changed=false;
        for(int i=0; i<Tour.size(); i++){
            for(int j=i+2; j<Tour.size(); j++){
                reverse(Tour.begin()+i+1,Tour.begin()+j+1);
                double NewCost=Cost(Tour);
                if(NewCost<Curr){
                    Changed=true;
                    break;
                }
                reverse(Tour.begin()+i+1,Tour.begin()+j+1);
            }
            if(Changed) break;
        }
        if(!Changed) break;
    }
}
int main(){
    int j,idx,maxidx,minidx;
    scanf("%d",&n);
    for(int i=0; i<n; i++) scanf("%d %lf %lf",&j,&points[i].lat,&points[i].lon);
    vector<int> temp;
    double cost;
    double maxCost=-1;
    double minCost=9999999999;
    //=======Task1============
    for(int i=0; i<5; i++){
        temp.clear();
        temp=nearestNeighbourHeuristic();
        cost=Cost(temp);
        if(cost<minCost){
            minCost=cost;
            minidx=temp[0];
        }
        else if(cost>maxCost){
            maxCost=cost;
            maxidx=temp[0];
        }
    }
    cout<<"NN Heuristic"<<endl;
    cout<<"------------------"<<endl;
    cout<<"MinCost: "<<minCost<<" Minidx: "<<minidx+1<<endl;
    printSolution(nearestNeighbourHeuristic2(minidx));
    cout<<"MaxCost: "<<maxCost<<" Maxidx: "<<maxidx+1<<endl;
    printSolution(nearestNeighbourHeuristic2(maxidx));
    //=============Task2==================
    maxCost=-1;
    minCost=9999999999;
    int minIdx=minidx;
    int minTouridx,maxTourIdx;
    vector< vector<int> > tours;
    for(int i=0; i<10; i++){
        temp.clear();
        temp=nearest5NeighbourHeuristic(minIdx);
        tours.push_back(temp);
        cost=Cost(temp);
        if(cost<minCost){
            minCost=cost;
            minidx=temp[0];
            minTouridx=i;
        }
        else if(cost>maxCost){
            maxCost=cost;
            maxidx=temp[0];
            maxTourIdx=i;
        }
    }
    cout<<"NN5 Heuristic"<<endl;
    cout<<"------------------"<<endl;
    cout<<"MinCost: "<<minCost<<endl<<"Tour: ";
    for(int i=0; i<temp.size(); i++){
        cout<<tours[minTouridx][i]+1<< " ";
    }
    cout<<endl;
    cout<<"MaxCost: "<<maxCost<<endl<<"Tour: ";
    for(int i=0; i<temp.size(); i++){
        cout<<tours[maxTourIdx][i]+1<< " ";
    }
    cout<<endl;
}
