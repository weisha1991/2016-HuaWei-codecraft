#ifndef __ROUTE_H__
#define __ROUTE_H__
#include "lib_io.h"
#include "lib_record.h"
#include "lib_time.h"
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <ctime>
#include <time.h>
#include <fstream>
#include <random>
#include <iterator>
#include <limits.h>
#include <queue>
#include <functional>
#include <utility>
#include <unordered_map>
#include <string.h>
#include "macrodef.h"
using namespace std;


struct edgeNode {
	int linkid;
	int from;
	int to;
	int weight;
};
typedef edgeNode* edgeNode_ptr;

struct adjListNode {
	int to;
	int weight;
	int linkId;
	adjListNode *next;
};
typedef adjListNode* adjListNode_ptr;

struct vertexNode {
	int adjacentNum;
	int vertexType;
	adjListNode_ptr  adjList_head;
	adjListNode_ptr adjList_end;
};
typedef vertexNode * vertexNode_ptr;

struct Ant {
	int    weight;
	int mustPass_number;
	double   fitness_val;
	int  cur_pos;
	int*   pathNode;
	int cur_edge;
	int* pathEdge;
	int *nodeVisited;
	int*   tabuList;

};
typedef Ant* Ant_ptr;

struct graphInfo {
	int   startVertex;
	int   endVertex;
	int   subVertex1_num;
	int  subVertex2_num;
	int   node_num;
	int real_edge_num;
	int *subVertex1;
	int *subVertex2;

};
typedef graphInfo* graphInfo_ptr;

struct acoInfo {
	int   ant_num;
	int   alpha;
	int   beta;
	int   edge_maxWeight;
	double   best_fitval;
	int   best_pathLen;
	int pathNodes_num;
	int mustNum;
	int*  bestPathNode;
	int* bestPathEdge;
	double *pheromone;
	double *heuristicVal_vec;
};

struct recordStruct {
	vector<int> pathEdge;
	int weight;
};
typedef acoInfo* acoInfo_ptr;

struct Node {
	int to,  weight;
	Node(int t, int w) :to(t), weight(w) {}
	bool operator <(const Node &rhs)const { return weight > rhs.weight; }
};

int compute_sameArcs(const vector<int> &edge_path1, const vector<int> &edge_path2,int edge_num);

vertexNode_ptr createAdjList(graphInfo_ptr g, edgeNode_ptr e);

acoInfo_ptr acoInfo_init(graphInfo_ptr g, edgeNode_ptr edgeArray);
void reset_acoInfo(acoInfo_ptr acoPtr, graphInfo_ptr g, const vector<int> &edge_path1);

Ant_ptr createAntColony(acoInfo_ptr acoPtr, graphInfo_ptr g, int setAntNode);

void reset_antState(Ant_ptr antVec, acoInfo_ptr acoPtr, graphInfo_ptr g, int setAntNode);
void reset_vertexNodeType(vertexNode_ptr headPtr, graphInfo_ptr g, int choose_flag);


void ant_search(vertexNode_ptr head_ptr, Ant_ptr antVec,
	acoInfo_ptr mmasAlog, graphInfo_ptr graph, edgeNode_ptr edgeArray, int choosePath_flag);

void findBestAnt(Ant_ptr antVec,acoInfo_ptr acoPtr,graphInfo_ptr g,vertexNode_ptr head_ptr,int choosePath) ;

void calcFitval(Ant_ptr antVec, acoInfo_ptr acoPtr, graphInfo_ptr g, int choosePath);

void update_pheromone(Ant_ptr antVec,acoInfo_ptr acoPtr,graphInfo_ptr g,vertexNode_ptr head_ptr,int choosePath);

void DisplayBestWay(acoInfo_ptr acoPtr, graphInfo_ptr g, vertexNode_ptr headPtr);

void print_twoPathInfo(vector<int> edge_path1, vector<int> edge_path2, edgeNode_ptr edgeArray);

int compute_sameArcsV1(vector<int> edge_path1, vector<int> edge_path2);

void freeData(edgeNode_ptr edgeArray, graphInfo_ptr graph, acoInfo_ptr acoPtr, vertexNode_ptr headPtr, Ant_ptr antVec, int **adjMatrix);
	int parse_topo(edgeNode_ptr edgeArray,
	char *topo[MAX_EDGE_NUM],
	int edge_num);
void parse_demand(graphInfo_ptr graph, char *demand[], int demand_num);
void print_graphInfo(graphInfo_ptr g);
void print_topo(edgeNode_ptr e, graphInfo_ptr g, fstream&);
void count_degreeInfo(edgeNode_ptr edgeArray, int edge_num, int *indegree, int *outdegree);
void writeAdjList_toFile(vertexNode_ptr headPtr,int node_num,ofstream &os);

void search_route(char *graph[MAX_EDGE_NUM], int edge_num, char *condition[MAX_DEMAND_NUM], int demand_num);

#endif
