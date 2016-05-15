#ifndef __ROUTE_H__
#define __ROUTE_H__

#include "lib_io.h"
#include "lib_record.h"
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <stdio.h>
#include <stack>
#include <string.h>
#include <vector>
#include <queue>
#include <random>
#include <chrono>
#include <map>
#include <limits.h>
#include <iterator>
using namespace std;
typedef pair<int, int> edge;

void parse_file(char *graph[MAX_EDGE_NUM], int edge_num, char *condition[MAX_DEMAND_NUM], int demand_num);
void preprocess();
void build_graph();

int dijkstra(int node_num, int startVertex, int endVertex, const vector<vector<edge>> &adjList, vector<int> &path);
int dijkstra_New(int node_num, int startVertex, int endVertex, const vector<vector<edge>> &adjList, vector<int> &path);

pair<vector<int>,vector<int>> dijkstra_singleSource(int node_num, int startVertex,const vector<vector<edge>> &adjList);
int dijkstra_toSubset(int node_num, int startVertex, vector<int> subvertex, const vector<vector<edge>> &adjList, vector<int> &path);
int path_midNode(const vector<int> &midnode_seq, vector<int> &paths, int vertexs, const vector<vector<edge>> &adjList);

pair<int,vector<int>> forward_findPath(const vector<int> &pass_nodes);
pair<int,vector<int>> backward_findPath(const vector<int> &pass_nodes);

void print_path(const vector<int>&);
size_t compute_multiEdges(const vector<int>&,const vector<int>&);

void update_graphInfo(bool isReverse,const vector<int> &result_nodes);
void search_route(char *graph[MAX_EDGE_NUM], int edge_num, char *condition[MAX_DEMAND_NUM], int demand_num);

#endif
