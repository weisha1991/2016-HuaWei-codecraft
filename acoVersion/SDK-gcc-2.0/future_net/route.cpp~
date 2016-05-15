#include "route.h"
using namespace std;

// #define DEBUG


int *global_tabuList;
vector<double> transProb;
vector<bool> temp_tabu;
int *indegree;
int *outdegree;
int max_weight=0;
unordered_multimap<long, int> vertexTolinkId;
vector<vector<Node>> adjList_singleEdge;
vector<vector<int>> adjMatrix;

default_random_engine randSeed(time(NULL));
uniform_real_distribution<double> uniReal(0,100.0);
vector<int> vertex_id;
unordered_map<int,int> r_vertexId;
vector<int> edge_id;
unordered_map<int,int> r_edgeId;

void recordAntPath(int loop,Ant_ptr antVec,vertexNode_ptr headPtr,acoInfo_ptr acoPtr,graphInfo_ptr g,ofstream &outStream) {
	outStream<<"*******************************************************"<<endl;
	outStream<<"the "<<loop<<"th iteration"<<endl;
	for (int i = 0; i < acoPtr->ant_num; ++i)
	{
		// if(antVec[i].mustPass_number!=g->subVertex1_num)
		// 	continue;
		// if(antVec[i].pathNode[antVec[i].cur_pos]!=g->endVertex)
		// 	continue;
		outStream << "ant " << i << " path [weight=" << antVec[i].weight <<",must pass="<<antVec[i].mustPass_number<< "]:" << endl;
		for (int n = 0; n <= antVec[i].cur_pos; ++n) {
			int curNode = antVec[i].pathNode[n];
			if (headPtr[curNode].vertexType == NODE_TYPE_MID)
				outStream << "{" << curNode << "} ";
			else
				outStream << curNode << " ";
		}
		outStream<<endl;
		outStream<<"ant "<<i<<"'s tabuList:";
		for(int j=0;j<g->node_num;++j){
			if(antVec[i].nodeVisited[j]==ANT_TABU)
				outStream<<vertex_id[j]<<",";
		}
		outStream << endl;
	}
}

int dijkstra(int start,int end,vertexNode_ptr head_ptr,graphInfo_ptr g,vector<int> &path)
{
	priority_queue<Node> Q;
	vector<int> dist(g->node_num, MAX_DIST_VAL), dad(g->node_num, -1);
	Q.push(Node(start, 0));
	dist[start] = 0;

	while (!Q.empty())
	{
		auto temp = Q.top();

		if (temp.to == end)
			break;
		Q.pop();
		int curNode = temp.to;
		for (size_t i = 0; i < adjList_singleEdge[curNode].size(); ++i)
		{
			if (temp_tabu[adjList_singleEdge[curNode][i].to]) {
				continue;
			}
			if (dist[curNode] + adjList_singleEdge[curNode][i].weight < dist[adjList_singleEdge[curNode][i].to])
			{
				dist[adjList_singleEdge[curNode][i].to] = dist[curNode] + adjList_singleEdge[curNode][i].weight;
				dad[adjList_singleEdge[curNode][i].to] = curNode;
				Q.push(Node(adjList_singleEdge[curNode][i].to,dist[adjList_singleEdge[curNode][i].to]));
			}
		}
	}
	if (dist[end] < MAX_DIST_VAL) {
		for (int i = end; i != -1; i = dad[i])
			path.insert(path.begin(), i);
		return dist[end];
	}
	return -1;
}
vector<int> dijkstra_SingleSource(int start,vertexNode_ptr head_ptr,graphInfo_ptr g)
{
	priority_queue<Node> Q;
	vector<int> dist(g->node_num, MAX_DIST_VAL), dad(g->node_num, -1);
	Q.push(Node(start, 0));
	dist[start] = 0;

	while (!Q.empty())
	{
		auto temp = Q.top();
		Q.pop();
		int curNode = temp.to;
		for (size_t i = 0; i < adjList_singleEdge[curNode].size(); ++i)
		{

			if (dist[curNode] + adjList_singleEdge[curNode][i].weight < dist[adjList_singleEdge[curNode][i].to])
			{
				dist[adjList_singleEdge[curNode][i].to] = dist[curNode] + adjList_singleEdge[curNode][i].weight;
				dad[adjList_singleEdge[curNode][i].to] = curNode;
				Q.push(Node(adjList_singleEdge[curNode][i].to,dist[adjList_singleEdge[curNode][i].to]));
			}
		}
	}
	vector<int> nodes_notArrive;
	for(int i=0;i<g->node_num;++i){
		if(dist[i]==MAX_DIST_VAL)
			nodes_notArrive.push_back(i);
	}
	return nodes_notArrive;
}

pair<int,vector<int>> optimize_path1(const vector<int> &result_nodes1,
	graphInfo_ptr graph,
	vertexNode_ptr head_ptr) {

	temp_tabu.assign(graph->node_num, false);

	vector<int> node_seq;
	for (auto node : result_nodes1) {
		if (head_ptr[node].vertexType == NODE_TYPE_START)
			node_seq.push_back(node);
		if (head_ptr[node].vertexType == NODE_TYPE_END)
			node_seq.push_back(node);
		if (head_ptr[node].vertexType == NODE_TYPE_MID)
			node_seq.push_back(node);
	}
	vector<int> opt_path;
	int sum_weight = 0;
	for (size_t i = 0; i < node_seq.size()-1;++i)
	{
		for (auto node : result_nodes1)
			temp_tabu[node] = true;
		for (auto node : opt_path)
			temp_tabu[node] = true;
		temp_tabu[node_seq[i]] = false;
		temp_tabu[node_seq[i + 1]] = false;
		vector<int> path;
		int w = dijkstra(node_seq[i], node_seq[i + 1], head_ptr, graph, path);
		auto begin_iter = find(result_nodes1.begin(), result_nodes1.end(), node_seq[i]);
		auto end_iter = find(result_nodes1.begin(), result_nodes1.end(), node_seq[i + 1]);
		int temp_weight = 0;
		for (auto iter = begin_iter; iter != end_iter; ++iter)
			temp_weight += adjMatrix[*iter][*(iter + 1)];
		if (w != -1 && temp_weight>w) {
			for (auto n : path)
				opt_path.push_back(n);
			sum_weight += w;
		}
		else  {
			for (auto iter = begin_iter; iter != end_iter; ++iter)
				opt_path.push_back(*iter);
			opt_path.push_back(*end_iter);
			sum_weight += temp_weight;
		}
	}
	return make_pair(sum_weight, opt_path);
}

pair<int,vector<int>> optimize_path2(const vector<int> &result_nodes1,
	const  vector <int> &result_nodes2,
	graphInfo_ptr graph,
	vertexNode_ptr head_ptr) {
	temp_tabu.assign(graph->node_num, false);
	for (auto node : result_nodes1)
		temp_tabu[node] = true;
	vector<int> node_seq;
	for (auto node : result_nodes2) {
		if (head_ptr[node].vertexType == NODE_TYPE_START)
			node_seq.push_back(node);
		if (head_ptr[node].vertexType == NODE_TYPE_END)
			node_seq.push_back(node);
		if (head_ptr[node].vertexType == NODE_TYPE_MID)
			node_seq.push_back(node);
	}
	vector<int> optimizePath2;
	int sum_weight=0;
	for (size_t i = 0; i < node_seq.size()-1; ++i) {

		for (auto node : result_nodes2)
			temp_tabu[node] = true;
		for (auto node : optimizePath2)
			temp_tabu[node] = true;
		temp_tabu[node_seq[i]] = false;
		temp_tabu[node_seq[i + 1]] = false;
		vector<int> path;
		int w=dijkstra(node_seq[i], node_seq[i + 1], head_ptr, graph, path);
		auto begin_iter = find(result_nodes2.begin(), result_nodes2.end(), node_seq[i]);
		auto end_iter = find(result_nodes2.begin(), result_nodes2.end(), node_seq[i + 1]);
		int temp_weight = 0;
		for (auto iter = begin_iter; iter != end_iter; ++iter)
			temp_weight += adjMatrix[*iter][*(iter + 1)];
		if (w != -1&&temp_weight>w) {
			for (auto n : path)
				optimizePath2.push_back(n);
			sum_weight += w;
		}
		else{
			for (auto iter = begin_iter; iter != end_iter; ++iter)
				optimizePath2.push_back(*iter);
			optimizePath2.push_back(*end_iter);
			sum_weight += temp_weight;
		}
	}
	return make_pair(sum_weight,optimizePath2);
}

void search_route(char *topo[MAX_EDGE_NUM], int edge_num, char *demand[MAX_DEMAND_NUM], int demand_num)
{
#ifdef DEBUG
	clock_t begin_time, end_time;
	double totalTime=0.0;
	begin_time = clock();
#endif
	graphInfo_ptr graph = (graphInfo_ptr)malloc(sizeof(graphInfo));
	edgeNode_ptr edgeArray = (edgeNode_ptr)malloc(edge_num*sizeof(edgeNode));//存原图数据

	int ret=parse_topo(edgeArray, topo, edge_num);
	graph->node_num = ret;
	graph->real_edge_num = edge_num;
	parse_demand(graph, demand, demand_num);

	indegree = new int[graph->node_num];
	outdegree = new int[graph->node_num];
	count_degreeInfo(edgeArray, edge_num, indegree, outdegree);

	global_tabuList = (int*)calloc(graph->node_num, sizeof(int));

	adjMatrix.resize(graph->node_num);
	for (size_t i = 0; i < adjMatrix.size(); ++i) {
		adjMatrix[i].assign(graph->node_num, MAX_DIST_VAL);
	}
	for (int i = 0; i < graph->real_edge_num; ++i) {
		if(edgeArray[i].weight<adjMatrix[r_vertexId[edgeArray[i].from]][r_vertexId[edgeArray[i].to]])
			adjMatrix[r_vertexId[edgeArray[i].from]][r_vertexId[edgeArray[i].to]] = edgeArray[i].weight;
	}
	adjList_singleEdge.resize(graph->node_num);
	for (size_t i = 0; i < adjMatrix.size(); ++i)
	{
		for (size_t j = 0; j < adjMatrix[i].size(); ++j) {
			if (adjMatrix[i][j] < MAX_DIST_VAL) {
				adjList_singleEdge[i].push_back(Node(j, adjMatrix[i][j]));
			}
		}
	}

#ifdef DEBUG
	print_graphInfo(graph);
	cout << "node's indegree equal zero:";
#endif
	for (int i = 0; i < graph->node_num; ++i)
		if (indegree[i] == 0) {
#ifdef DEBUG
			cout << vertex_id[i] << ",";
#endif
			global_tabuList[i] = -1;
		}
#ifdef DEBUG
	cout << endl;
	cout << "node's outdegree equal zero:";
#endif
	for (int i = 0; i < graph->node_num; ++i)
		if (outdegree[i] == 0) {
#ifdef DEBUG
			cout << vertex_id[i] << ",";
#endif
			global_tabuList[i] = -1;
		}
#ifdef DEBUG
	cout << endl;
#endif
	global_tabuList[graph->startVertex] = 0;
	global_tabuList[graph->endVertex] = 0;

	vertexNode_ptr head_ptr = createAdjList(graph, edgeArray);
/*
	ofstream adjListFile("adjListCheck.txt");
	writeAdjList_toFile(head_ptr,graph->node_num,adjListFile);
*/
	// ofstream edgeArrayFile("edgeArray.txt");
	// for(int i=0;i<graph->real_edge_num;++i){
	// 	edgeArrayFile<<edgeArray[i].linkid<<","<<edgeArray[i].from<<","<<edgeArray[i].to
	// 		<<","<<edgeArray[i].weight<<endl;
	// }
	acoInfo_ptr mmasAlgo=acoInfo_init(graph,edgeArray);
	Ant_ptr Ant_vec = createAntColony(mmasAlgo,graph,graph->startVertex);
	int i = 0;
	transProb.assign(graph->node_num, 0.0);
#ifdef DEBUG
	end_time = clock();
	totalTime = (double)(end_time - begin_time) / CLOCKS_PER_SEC;
	cout << "read file and create data structre use " << totalTime << " s" << endl;
#endif
#ifdef DEBUG
	ofstream out("record.txt");
	begin_time = clock();
#endif

	vector<int> nodes_discard=dijkstra_SingleSource(graph->startVertex,head_ptr,graph);
#ifdef DEBUG
	cout<<"起点不可达的点("<<nodes_discard.size()<<" nodes):";
	for(auto node:nodes_discard)
		cout<<node<<",";
	cout<<endl;
#endif
	reset_vertexNodeType(head_ptr, graph, 2);
	while (i<10) {
		ant_search(head_ptr, Ant_vec, mmasAlgo, graph, edgeArray,2);
		findBestAnt(Ant_vec, mmasAlgo, graph,head_ptr,2);
		calcFitval(Ant_vec, mmasAlgo, graph,2);
		update_pheromone(Ant_vec, mmasAlgo, graph,head_ptr,2);
		// recordAntPath(i,Ant_vec, head_ptr, mmasAlgo, graph, out);
		++i;
	}


	vector<int> result_path1;
	vector<int> edge_path1;
	int w1 = mmasAlgo->best_pathLen;
	for (int i = 0; i < mmasAlgo->pathNodes_num; ++i) {
		result_path1.push_back(mmasAlgo->bestPathNode[i]);
	}
	for (int i = 0; i <mmasAlgo->pathNodes_num-1; ++i) {
		edge_path1.push_back(mmasAlgo->bestPathEdge[i]);
	}
	if (result_path1.empty()) {
		cout<<"cant not find path1"<<endl;
		return;
	}

#ifdef DEBUG

	end_time = clock();
	totalTime = (double)(end_time - begin_time) / CLOCKS_PER_SEC;
	cout << "find first path used " << totalTime << " s" << endl;
#endif

 #ifdef DEBUG
	begin_time = clock();
#endif
	reset_acoInfo(mmasAlgo, graph,edge_path1);
	reset_vertexNodeType(head_ptr, graph,1);
	i = 0;
	vector<recordStruct> save_iter_goodAnts;
	while (i++ < 10) {
		ant_search(head_ptr, Ant_vec, mmasAlgo, graph, edgeArray, 1);
		findBestAnt(Ant_vec, mmasAlgo, graph,head_ptr, 1);
		calcFitval(Ant_vec, mmasAlgo, graph, 1);
		update_pheromone(Ant_vec, mmasAlgo, graph, head_ptr,1);

		for (int i = 0; i < mmasAlgo->ant_num; ++i) {
			if (Ant_vec[i].pathNode[Ant_vec[i].cur_pos] == graph->endVertex
				&&Ant_vec[i].mustPass_number == graph->subVertex1_num)
			{
				recordStruct tmp;
				for (int j = 0; j <=Ant_vec[i].cur_edge; ++j) {
					tmp.pathEdge.push_back(Ant_vec[i].pathEdge[j]);
				}
				tmp.weight = Ant_vec[i].weight;
				save_iter_goodAnts.push_back(tmp);
			}

		}
	}


#ifdef DEBUG
	if(save_iter_goodAnts.empty())
	{
		cout<<"path2  not find"<<endl;
		return;
	}
#endif
	int minSameCnt = INT_MAX;
	int minLen = INT_MAX;
	vector<int> edge_path2;

	for (size_t i = 0; i < save_iter_goodAnts.size();++i){

		int cc = compute_sameArcs(edge_path1, save_iter_goodAnts[i].pathEdge,graph->real_edge_num);
		if ( cc<minSameCnt||(cc==minSameCnt&&save_iter_goodAnts[i].weight<minLen)) {
			minSameCnt = cc;
			minLen = save_iter_goodAnts[i].weight;
			edge_path2 = save_iter_goodAnts[i].pathEdge;
		}
	}
#ifdef DEBUG
	cout<<"***********************************"<<endl;
	cout<<"优化前"<<endl;
	cout << "path seq1:";
	for (auto e: edge_path1)
		cout << edgeArray[e].from << "->";
	cout << edgeArray[edge_path1.back()].to << endl;

	cout << "path seq2:";
	for (auto e : edge_path2)
		cout << edgeArray[e].from<< "->";

	cout <<edgeArray[edge_path2.back()].to<< endl;

	cout <<"重边数" <<minSameCnt << endl;
	cout << "path1 weight=" << w1 << endl;
	cout << "path2 weight=" << minLen << endl;
	cout << "权值和"<<w1 + minLen << endl;
#endif


	vector<int> result_path2;
	for(auto e:edge_path2)
	{
		result_path2.push_back(r_vertexId[edgeArray[e].from]);
	}
	result_path2.push_back(r_vertexId[edgeArray[edge_path2.back()].to]);
#ifdef DEBUG
	end_time = clock();
	totalTime = ((double)(end_time - begin_time) / CLOCKS_PER_SEC);
	cout << "find second path used " << totalTime << " s" << endl;
#endif
	// for(auto e:edge_path2)
	// 	record_result(WORK_PATH,edge_id[e]);
	// for(auto e:edge_path1)
	// 	record_result(BACK_PATH,edge_id[e]);

//***********************优化权值

	// optimize path 2
	reset_vertexNodeType(head_ptr, graph, 1);
	auto ret_val2 = optimize_path2(result_path1, result_path2, graph, head_ptr);
	auto endSame_iter = unique(ret_val2.second.begin(), ret_val2.second.end());
	ret_val2.second.erase(endSame_iter, ret_val2.second.end());
	vector<int> opt_edgePath2;
	for (auto iter = ret_val2.second.begin(); iter != ret_val2.second.end() - 1; ++iter)
	{
		auto temp = vertexTolinkId.find(hash_func(vertex_id[*iter], vertex_id[*(iter + 1)]));
		opt_edgePath2.push_back(r_edgeId[temp->second]);
	}
	//optimize path 1;
// 	reset_vertexNodeType(head_ptr, graph, 2);
// 	auto ret_val1 = optimize_path2(ret_val2.second, result_path1, graph, head_ptr);

// 	endSame_iter=unique(ret_val1.second.begin(), ret_val1.second.end());
// 	ret_val1.second.erase(endSame_iter, ret_val1.second.end());
// #ifdef DEBUG
// 	set<int> checkSet(ret_val1.second.begin(), ret_val1.second.end());
// 	if (checkSet.size() != ret_val1.second.size())
// 		cout << "path error"<<endl;
// #endif
// 	vector<int> opt_edgePath1;
// 	for (auto iter = ret_val1.second.begin(); iter != ret_val1.second.end() - 1; ++iter)
// 	{
// 		auto temp = vertexTolinkId.find(hash_func(vertex_id[*iter],vertex_id[*(iter+1)]));
// 		opt_edgePath1.push_back(r_edgeId[temp->second]);
// 	}

#ifdef DEBUG
	end_time = clock();
	totalTime = ((double)(end_time - begin_time) / CLOCKS_PER_SEC) - totalTime;
	cout << "optimize the path used " << totalTime << " s" << endl;
#endif
#ifdef DEBUG
	cout<<"***********************************"<<endl;
	cout<<"优化后"<<endl;
	cout<<"path seq1:";
	for (auto e: edge_path1)
		cout << edgeArray[e].from << "->";
	cout << edgeArray[edge_path1.back()].to << endl;

	cout << "path seq2:";
	for (auto e : opt_edgePath2)
		cout << edgeArray[e].from<< "->";
	cout << edgeArray[opt_edgePath2.back()].to << endl;

	cout<<"路径1优化权值："<<w1<<endl;
	cout<<"路径2优化权值:"<<ret_val2.first<<endl;
	cout<<"优化后的权重和:"<<w1+ret_val2.first<<endl;
	cout<<"重边数"<<compute_sameArcs(edge_path1,opt_edgePath2,graph->real_edge_num)<<endl;
	cout<<"重边数"<<compute_sameArcsV1(edge_path1,opt_edgePath2)<<endl;

#endif

	for(auto e:opt_edgePath2)
		record_result(WORK_PATH,edge_id[e]);
	for(auto e:edge_path1)
		record_result(BACK_PATH,edge_id[e]);
	//free(indegree)
	//delete[] indegree;
	//free(outdegree)
	//delete [] outdegree;
	//free(global_tabuList);
	//delete[] transProb;
	// freeData(edgeArray, graph, mmasAlgo, head_ptr, Ant_vec,adjMatrix);

	return;
}

int compute_sameArcsV1(vector<int> edge_path1, vector<int> edge_path2){
	sort(edge_path1.begin(),edge_path1.end());
	sort(edge_path2.begin(),edge_path2.end());
	vector<int> tmp;
	set_intersection(edge_path1.begin(),edge_path1.end(),edge_path2.begin(),edge_path2.end(),back_inserter(tmp));
	return tmp.size();
}

vector<int> edge_cnt;
int compute_sameArcs(const vector<int> &edge_path1, const vector<int> &edge_path2,int edge_num)
{
	edge_cnt.assign(edge_num,0);
	for(auto e:edge_path1)
		edge_cnt[e]++;
	for(auto e:edge_path2)
		edge_cnt[e]++;
	int ret_val=0;
	for(auto n:edge_cnt){
		if(n>1){
			ret_val++;
		}
	}

	return ret_val;
}

void print_twoPathInfo(vector<int> edge_path1,vector<int> edge_path2,edgeNode_ptr edgeArray){
	cout << "***********************************************" << endl;
	sort(edge_path1.begin(), edge_path1.end());
	sort(edge_path2.begin(), edge_path2.end());
	vector<int> mutiedge_set;
	set_intersection(edge_path1.begin(), edge_path1.end(), edge_path2.begin(), edge_path2.end(), back_inserter(mutiedge_set));
	cout << "same edge count=" << mutiedge_set.size() << endl;
	cout << "save edge Id:";
	for (auto e : mutiedge_set)
	cout << e << "(" << edgeArray[e].from << "->" << edgeArray[e].to << "),";
	cout << endl;
}

void print_topo(edgeNode_ptr e,graphInfo_ptr g, fstream &out) {
	for (int i = 0; i < g->real_edge_num; ++i) {
		out << e[i].linkid << "," << e[i].from << "," << e[i].to << "," << e[i].weight << endl;
	}
}
void print_graphInfo(graphInfo_ptr g) {
	cout << "***********************************************" << endl;
	cout << "startNode=" << vertex_id[g->startVertex] << ",endNode=" << vertex_id[g->endVertex] << endl;
	cout << "node nom=" << g->node_num << endl;
	cout <<"edge num=" << g->real_edge_num << endl;
	cout<<"max edge weight="<<max_weight<<endl;
	cout << "must pass set 1 (" << g->subVertex1_num << " nodes):" << endl;
	for (int i = 0; i < g->subVertex1_num; ++i) {
		cout <<vertex_id[g->subVertex1[i]] << ",";
	}
	cout << endl;
	cout << "must pass set 2 (" << g->subVertex2_num << " nodes):" << endl;
	for (int i = 0; i < g->subVertex2_num; ++i) {
		cout <<vertex_id[g->subVertex2[i]]<< ",";
	}
	cout << endl;
}

void freeData(edgeNode_ptr edgeArray, graphInfo_ptr graph,acoInfo_ptr acoPtr,vertexNode_ptr headPtr,Ant_ptr antVec,int **adjMatrix) {
	if (adjMatrix) {
		for (int i = 1; i < graph->node_num; ++i) {
			if (adjMatrix[i] != NULL)
				free(adjMatrix[i]);
		}
		free(adjMatrix);
	}
	free(edgeArray);
	if (headPtr != NULL) {
		//free the memory of adjList
		for (int i = 1; i < graph->node_num; ++i) {
			adjListNode_ptr curPtr = headPtr[i].adjList_head;
			if (curPtr == NULL)
				continue;
			adjListNode_ptr nextPtr = curPtr->next;
			while (curPtr != NULL) {
				free(curPtr);
				curPtr = nextPtr;
				if (curPtr == NULL)
					break;
				nextPtr = nextPtr->next;
			}
		}
		free(headPtr);
	}
	free(graph->subVertex1);
	free(graph->subVertex2);
	free(graph);

	//free the memory of antVec;
	if (antVec != NULL) {
		for (int i = 0; i < acoPtr->ant_num; ++i) {
			if (antVec[i].pathNode != NULL)
				delete[] antVec[i].pathNode;
			if (antVec[i].pathEdge != NULL)
				delete[] antVec[i].pathEdge;
			// if (antVec[i].tabuList != NULL)
			// 	delete[] antVec[i].tabuList;
			if (antVec[i].nodeVisited != NULL)
				delete[] antVec[i].nodeVisited;

		}
		free(antVec);
	}
	//free the memeory of acoInfo
	if (acoPtr->bestPathNode != NULL)
		delete[] acoPtr->bestPathNode;
	if (acoPtr->bestPathEdge != NULL)
		delete[] acoPtr->bestPathEdge;
	if (acoPtr->pheromone!=NULL)
		delete[] acoPtr->pheromone;
	if (acoPtr->heuristicVal_vec != NULL)
		delete[] acoPtr->heuristicVal_vec;
	free(acoPtr);



}

void count_degreeInfo(edgeNode_ptr edgeArray,int edge_num,int *indegree,int *outdegree)
{
	for (int i = 0; i < edge_num; ++i) {
		indegree[r_vertexId[edgeArray[i].to]] +=1 ;
		outdegree[r_vertexId[edgeArray[i].from]] += 1;
	}
}

int parse_topo(edgeNode_ptr edgeArray,char *topo[MAX_EDGE_NUM],  int edge_num)
{

	int j = 0;
	int max_Nodeid = 0;
	for (int i = 0; i<edge_num; ++i) {
		char *token = strtok(topo[i], ",");
		int temp[4];
		j = 0;
		while (token != NULL) {
			temp[j] = atoi(token);
			token = strtok(NULL, ",");
			++j;
		}
		if(temp[3]>max_weight)
			max_weight=temp[3];
		if (temp[1] > max_Nodeid)
			max_Nodeid = temp[1];
		if (temp[2] > max_Nodeid)
			max_Nodeid = temp[2];
		edgeArray[i].linkid=i;
		edgeArray[i].from = temp[1];
		edgeArray[i].to = temp[2];
		edgeArray[i].weight = temp[3];
		if (temp[1] == temp[2]) {
#ifdef DEBUG
			cout <<temp[1] <<"自环" << endl;
			exit (1);
#endif
		}
		vertex_id.push_back(temp[1]);
		vertex_id.push_back(temp[2]);
		edge_id.push_back(temp[0]);
		vertexTolinkId.insert(make_pair(hash_func(temp[1], temp[2]), temp[0]));
	}

	sort(vertex_id.begin(),vertex_id.end());
	auto end_iter=unique(vertex_id.begin(),vertex_id.end());
	vertex_id.erase(end_iter,vertex_id.end());
#ifdef DEBUG
	for(auto iter=vertex_id.begin();iter!=vertex_id.end()-1;++iter){
		if(*(iter+1)-*(iter)!=1){
			cout<<"顶点序号不连续"<<endl;
			break;
		}
	}
#endif
	for(size_t i=0;i<vertex_id.size();++i)
		r_vertexId.insert(make_pair(vertex_id[i],i));
	for(size_t i=0;i<edge_id.size();++i)
		r_edgeId.insert(make_pair(edge_id[i],i));
#ifdef DEBUG
	ofstream vertexIdMap("vertexIdmap.txt");
	for(size_t i=0;i<vertex_id.size();++i){
		vertexIdMap<<i<<"<->"<<vertex_id[i]<<endl;
	}
	ofstream edgeIdMap("edgeIdmap.txt");
	for(size_t i=0;i<edge_id.size();++i){
		edgeIdMap<<i<<"<->"<<edge_id[i]<<endl;
	}
#endif
	return vertex_id.size();
}

void parse_demand(graphInfo_ptr graph ,char *demand[],int demand_num) {
	
	string parse_line1(demand[0]);	
	string parse_line2(demand[1]);
	size_t endpos1=parse_line1.find("NA");
	size_t endpos2=parse_line2.find("NA");
	if(endpos1!=string::npos&&endpos2==string::npos){//subset1 is empty
#ifdef DEBUG
	cout<<"subset 1 is empty"<<endl;
#endif	
		vector<int> line1;
		char *token=strtok(demand[0],",");
		while(token!=NULL){
			line1.push_back(atoi(token));
			token=strtok(NULL,",");
		}
		graph->startVertex=r_vertexId[line1[1]];
		graph->endVertex=r_vertexId[line1[2]];
		graph->subVertex1_num=0;
		graph->subVertex1=NULL;
		line1.clear();
		token=strtok(demand[1],",|");
		while(token!=NULL){
			line1.push_back(atoi(token));
			token=strtok(NULL,",|");
		}
		graph->subVertex2_num=line1.size()-3;
		graph->subVertex2 = (int*)calloc(line1.size() - 3, sizeof(int));
		for(int i=0;i<graph->subVertex2_num;++i){
			graph->subVertex2[i] = r_vertexId[line1[i+3]];
		}
	}
	else if(endpos2!=string::npos&&endpos1==string::npos){//subset2 is empty
#ifdef DEBUG
	cout<<"subset 2 is empty"<<endl;
#endif
		vector<int> line2;
		char *token=strtok(demand[1],",");
		while(token!=NULL){
			line2.push_back(atoi(token));
			token=strtok(NULL,",");
		}
		graph->startVertex=r_vertexId[line2[1]];
		graph->endVertex=r_vertexId[line2[2]];
		graph->subVertex2_num=0;
		graph->subVertex2=NULL;
		line2.clear();
		token=strtok(demand[0],",|");
		while(token!=NULL){
			line2.push_back(atoi(token));
			token=strtok(NULL,",|");
		}
		graph->subVertex1_num=line2.size()-3;
		graph->subVertex1 = (int*)calloc(line2.size() - 3, sizeof(int));
		for(int i=0;i<graph->subVertex1_num;++i){
			graph->subVertex1[i] = r_vertexId[line2[i+3]];
		}

	}
	else if(endpos1==string::npos&&endpos2==string::npos){
		vector<vector<int>> temp(2);
		for (int i = 0; i < demand_num; ++i) {
			char *token = strtok(demand[i], ",|");
			while (token != NULL) {
				temp[i].push_back(atoi(token));
				token = strtok(NULL, ",|");
			}
		}
		graph->startVertex = r_vertexId[temp[0][1]];
		graph->endVertex =r_vertexId[temp[0][2]];
		graph->subVertex1_num = temp[0].size() - 3;
		graph->subVertex1 = (int*)calloc(temp[0].size() - 3, sizeof(int));
		graph->subVertex2_num = temp[1].size() - 3;
		graph->subVertex2 = (int*)calloc(temp[1].size() - 3, sizeof(int));
		for (int i = 0; i < graph->subVertex1_num; ++i)
			graph->subVertex1[i] = r_vertexId[temp[0][i + 3]];
		for (int i = 0; i < graph->subVertex2_num; ++i)
			graph->subVertex2[i] = r_vertexId[temp[1][i + 3]];
	}
}

vertexNode_ptr createAdjList(graphInfo_ptr g,edgeNode_ptr e) {
	vertexNode_ptr head = (vertexNode_ptr)malloc(sizeof(vertexNode)*g->node_num);
	for (int i = 0; i < g->node_num; ++i) {
		head[i].adjList_head = NULL;
		head[i].adjList_end = NULL;
		head[i].adjacentNum = 0;
		head[i].vertexType = NODE_TYPE_COMM;
	}
	(head + g->startVertex)->vertexType = NODE_TYPE_START;
	(head + g->endVertex)->vertexType = NODE_TYPE_END;
	adjListNode_ptr temp = NULL;
	for (int i = 0; i < g->real_edge_num; ++i)
	{
		temp= (adjListNode_ptr)malloc(sizeof(adjListNode));
		temp->to = r_vertexId[e[i].to];
		temp->linkId = e[i].linkid;
		temp->weight = e[i].weight;
		temp->next = NULL;
		if ((head + r_vertexId[e[i].from])->adjList_head == NULL) {
			(head + r_vertexId[e[i].from])->adjList_head = temp;
			head[r_vertexId[e[i].from]].adjList_end = temp;
			(head + r_vertexId[e[i].from])->adjacentNum++;

		}
		else {
			head[r_vertexId[e[i].from]].adjList_end->next = temp;
			head[r_vertexId[e[i].from]].adjList_end = temp;
			head[r_vertexId[e[i].from]].adjacentNum++;
		}

	}
	for (int i = 0; i < g->subVertex1_num; ++i) {
		head[g->subVertex1[i]].vertexType = NODE_TYPE_MID;
	}
	return head;
}

void reset_vertexNodeType(vertexNode_ptr headPtr, graphInfo_ptr g, int choose_flag) {
	for (int i = 0; i < g->node_num; ++i) {
		headPtr[i].vertexType = NODE_TYPE_COMM;
	}
	headPtr[g->startVertex].vertexType = NODE_TYPE_START;
	headPtr[g->endVertex].vertexType = NODE_TYPE_END;
	if (choose_flag == 1){
		for (int i = 0; i < g->subVertex1_num; ++i) {
			headPtr[g->subVertex1[i]].vertexType = NODE_TYPE_MID;
		}
	}
	else if (choose_flag == 2) {
		for (int i = 0; i < g->subVertex2_num; ++i) {
			headPtr[g->subVertex2[i]].vertexType = NODE_TYPE_MID;
		}
	}
}


acoInfo_ptr acoInfo_init(graphInfo_ptr g,edgeNode_ptr edgeArray)
{
	acoInfo_ptr acoInfo_obj = (acoInfo_ptr)malloc(sizeof(acoInfo));
	acoInfo_obj->alpha = ANT_ALPHA;
	acoInfo_obj->beta = ANT_BETA;
	acoInfo_obj->ant_num = ANT_NUMBER;
	acoInfo_obj->edge_maxWeight = max_weight;
	acoInfo_obj->best_pathLen = MAX_DIST_VAL;
	acoInfo_obj->pathNodes_num = 0;
	acoInfo_obj->best_fitval = MAX_DIST_VAL;
	acoInfo_obj->bestPathNode = new int[g->node_num];
	acoInfo_obj->bestPathEdge = new int[g->real_edge_num];
	acoInfo_obj->pheromone = new double[g->real_edge_num];
	acoInfo_obj->heuristicVal_vec = new double[g->real_edge_num];
	acoInfo_obj->mustNum = 0;
	for (int i = 0; i < g->node_num; ++i)
		acoInfo_obj->bestPathNode[i] = INVALID;
	for (int i = 0; i < g->real_edge_num; ++i) {
		acoInfo_obj->bestPathEdge[i] = INVALID;
		acoInfo_obj->pheromone[i] = PHEROMONE_INIT;
		acoInfo_obj->heuristicVal_vec[i] = pow(MAX_WEIGHT - edgeArray[i].weight, ANT_BETA);
	}
	return acoInfo_obj;
}

void reset_acoInfo(acoInfo_ptr acoPtr,graphInfo_ptr g,const vector<int> &edge_path1) {
	acoPtr->alpha = ANT_ALPHA;
	acoPtr->beta = ANT_BETA;
	acoPtr->ant_num = ANT_NUMBER;
	acoPtr->edge_maxWeight = max_weight;
	acoPtr->best_pathLen = MAX_DIST_VAL;
	acoPtr->pathNodes_num = 0;
	acoPtr->mustNum = 0;
	acoPtr->best_fitval = MAX_DIST_VAL;
	for(int i=0;i<g->node_num;++i)
		acoPtr->bestPathNode[i]=INVALID;
	for(int i=0;i<g->real_edge_num;++i)
		acoPtr->bestPathEdge[i]=INVALID;
	for (int i = 0; i < g->real_edge_num; ++i) {
		acoPtr->pheromone[i] = PHEROMONE_INIT;
	}

	for (size_t i = 0; i < edge_path1.size(); ++i)
		acoPtr->pheromone[edge_path1[i]] = MIN_PHEROMONE;

	for (size_t i = 0; i < edge_path1.size(); ++i)
		acoPtr->heuristicVal_vec[edge_path1[i]] =0.0001;
}

Ant_ptr createAntColony(acoInfo_ptr acoPtr,graphInfo_ptr g,int setAntNode) {
	Ant_ptr antVec = (Ant_ptr)malloc(sizeof(Ant)*acoPtr->ant_num);
	for (int i = 0; i < acoPtr->ant_num; ++i) {
		antVec[i].weight = 0;
		antVec[i].fitness_val = 0.0;
		antVec[i].mustPass_number = 0;
		antVec[i].pathNode = new int[g->node_num];
		for (int j = 1; j < g->node_num; ++j)
			antVec[i].pathNode[j] = INVALID;
		antVec[i].cur_pos = 0;
		antVec[i].cur_edge = -1;
		antVec[i].pathNode[0] = setAntNode;
		// antVec[i].tabuList = new int[g->node_num];
		// for (int j = 0; j < g->node_num; ++j)
		// 	antVec[i].tabuList[j] = INVALID;
		antVec[i].pathEdge = new int[g->real_edge_num];
		for (int j = 0; j < g->real_edge_num; ++j)
			antVec[i].pathEdge[j] = INVALID;
		antVec[i].nodeVisited = new int[g->node_num];
		for (int j = 0; j < g->node_num; ++j)
			antVec[i].nodeVisited[j] = ANT_UNGO;
		antVec[i].nodeVisited[setAntNode] = ANT_GOST;

	}
	return antVec;
}

void reset_antState(Ant_ptr antVec,acoInfo_ptr acoPtr,graphInfo_ptr g,int setAntNode) {
	if (acoPtr != NULL) {
		for (int i = 0; i < acoPtr->ant_num; ++i)
		{
			antVec[i].weight = 0;
			antVec[i].fitness_val = 0.0;
			antVec[i].mustPass_number = 0;

			for (int j = 0; j < g->node_num; ++j)
				antVec[i].pathNode[j] = INVALID;
			antVec[i].pathNode[0] = setAntNode;
			antVec[i].cur_pos = 0;
			antVec[i].cur_edge = -1;

			for (int j = 0; j < g->real_edge_num; ++j)
				antVec[i].pathEdge[j] = INVALID;

			for (int j = 0; j < g->node_num; ++j)
			{
				antVec[i].nodeVisited[j] = ANT_UNGO;
				// antVec[i].tabuList[j] = INVALID;
			}
			antVec[i].nodeVisited[setAntNode] = ANT_GOST;

		}
	}
}

//every time call this function,params trans_prob must be malloc memort at first
void ant_search(vertexNode_ptr head_ptr,Ant_ptr antVec,
	acoInfo_ptr mmasAlog,graphInfo_ptr graph, edgeNode_ptr edgeArray,int choosePath_flag) {

	int tabuCount=0;
	int mustPassCount = 0;
	double probTotal = 0.0;
	int subVertex_num = 0;
	if (choosePath_flag == 1)
		subVertex_num = graph->subVertex1_num;
	else if (choosePath_flag == 2)
		subVertex_num = graph->subVertex2_num;

	reset_antState(antVec, mmasAlog, graph,graph->startVertex);
	for (int i = 0; i < mmasAlog->ant_num; ++i)
	{
		int loops = 0;
		tabuCount = 0;
		mustPassCount = 0;

		while (antVec[i].pathNode[antVec[i].cur_pos] != graph->endVertex)
		{
			loops++;
			int curNode = antVec[i].pathNode[antVec[i].cur_pos];
			transProb.assign(graph->node_num, 0.0);

			probTotal = 0.0;

			//access the adjList
			adjListNode_ptr curPtr = head_ptr[curNode].adjList_head;
			for (int j = 0; j < head_ptr[curNode].adjacentNum; ++j)
			{
				if (antVec[i].nodeVisited[curPtr->to] == ANT_GOST) {
					curPtr = curPtr->next;
					continue;
				}
				if (antVec[i].nodeVisited[curPtr->to] == ANT_TABU) {
					curPtr = curPtr->next;
					continue;
				}
				if (curPtr->to == graph->endVertex&&mustPassCount != subVertex_num) {
					curPtr = curPtr->next;
					continue;
				}
				// if (global_tabuList[curPtr->to] == -1) {
				// 	curPtr = curPtr->next;
				// 	continue;
				// }

				if (head_ptr[curPtr->to].vertexType == NODE_TYPE_MID&&antVec[i].cur_pos>1) {
					mmasAlog->pheromone[curPtr->linkId] *= 1.2;
					mmasAlog->pheromone[antVec[i].pathEdge[antVec[i].cur_edge]] *= 1.2;
				}

				transProb[curPtr->to]= pow(mmasAlog->pheromone[curPtr->linkId], ANT_ALPHA)*mmasAlog->heuristicVal_vec[curPtr->linkId];
				probTotal += transProb[curPtr->to];
				curPtr = curPtr->next;
			}

			curPtr= head_ptr[curNode].adjList_head;
			for (int j = 0; j < head_ptr[curNode].adjacentNum; ++j)
			{
				if (antVec[i].nodeVisited[curPtr->to] == ANT_GOST) {
					curPtr = curPtr->next;
					continue;
				}
				//if (antVec[i].edgeVisited[curPtr->linkId] == ANT_GOST) {
				//	curPtr = curPtr->next;
				//	continue;
				//}
				if (antVec[i].nodeVisited[curPtr->to] == ANT_TABU) {
					curPtr = curPtr->next;
					continue;
				}
				if (curPtr->to == graph->endVertex&&mustPassCount != subVertex_num) {
					curPtr = curPtr->next;
					continue;
				}
				/*if (global_tabuList[curPtr->to] == -1) {
					curPtr = curPtr->next;
					continue;
				}*/
				transProb[curPtr->to] /= probTotal;
				curPtr = curPtr->next;
			}
			double min_prob = 0.0, max_prob = 0.0;
			int selectedNode = -1;
			int selectedEdge = -1;
			double roulette = uniReal(randSeed)/100.0;
			if (roulette <= 0)
			{
				double proTemp = 0.0;
				curPtr = head_ptr[curNode].adjList_head;
				for (int j = 0; j < head_ptr[curNode].adjacentNum; ++j)
				{
					if (antVec[i].nodeVisited[curPtr->to] == ANT_GOST) {
						curPtr = curPtr->next;
						continue;
					}
					if (antVec[i].nodeVisited[curPtr->to] == ANT_TABU) {
						curPtr = curPtr->next;
						continue;
					}
					if (curPtr->to == graph->endVertex&&mustPassCount !=subVertex_num) {
						curPtr = curPtr->next;
						continue;
					}
					// if (global_tabuList[curPtr->to] == -1) {
					// 	curPtr = curPtr->next;
					// 	continue;
					// }

					if (transProb[curPtr->to] > proTemp) {
						proTemp = transProb[curPtr->to];
						selectedNode = curPtr->to;
						selectedEdge = curPtr->linkId;
					}
					curPtr = curPtr->next;
				}
			}
			else
			{
				curPtr = head_ptr[curNode].adjList_head;
				for (int j = 0; j < head_ptr[curNode].adjacentNum; ++j)
				{
					if (antVec[i].nodeVisited[curPtr->to] == ANT_GOST) {
						curPtr = curPtr->next;
						continue;
					}
					if (antVec[i].nodeVisited[curPtr->to] == ANT_TABU) {
						curPtr = curPtr->next;
						continue;
					}
					if (curPtr->to == graph->endVertex&&mustPassCount !=subVertex_num) {
						curPtr = curPtr->next;
						continue;
					}
					// if (global_tabuList[curPtr->to] == -1) {
					// 	curPtr = curPtr->next;
					// 	continue;
					// }
					max_prob += transProb[curPtr->to];
					if (roulette >= min_prob&&roulette <= max_prob) {
						selectedNode = curPtr->to;
						selectedEdge = curPtr->linkId;
						break;
					}
					else
						min_prob = max_prob;
					curPtr = curPtr->next;
				}
			}
			if (selectedNode == -1 && selectedEdge == -1)//cant not find any node to move
			{
				if (antVec[i].cur_pos == 0)
					break;
				if (head_ptr[curNode].vertexType == NODE_TYPE_MID) {
					mustPassCount--;
					antVec[i].mustPass_number--;
					antVec[i].nodeVisited[curNode]=ANT_UNGO;
					// antVec[i].pathNode[antVec[i].cur_pos] = INVALID;
					// antVec[i].pathEdge[antVec[i].cur_edge] = ANT_UNGO;
					// antVec[i].tabuList[tabuCount] = INVALID;
				}
				else {
					// antVec[i].tabuList[tabuCount] = antVec[i].pathNode[antVec[i].cur_pos];
					if (indegree[curNode] == 1 && outdegree[curNode] == 1)
						antVec[i].nodeVisited[curNode] = ANT_TABU;
					//if (indegree[curNode] == 1 && outdegree[curNode] == 1)
					//	global_tabuList[curNode] = -1;
				}
				antVec[i].pathNode[antVec[i].cur_pos] = INVALID;
				--antVec[i].cur_pos;
				int edgeId = antVec[i].pathEdge[antVec[i].cur_edge];
				antVec[i].weight -= edgeArray[edgeId].weight;
				antVec[i].pathEdge[antVec[i].cur_edge] = INVALID;
				--antVec[i].cur_edge;
				tabuCount++;
				if (tabuCount > graph->node_num*0.5 || loops > graph->node_num)
					break;
			}
			else
			{
				antVec[i].cur_pos++;
				antVec[i].pathNode[antVec[i].cur_pos] = selectedNode;
				antVec[i].cur_edge++;
				antVec[i].pathEdge[antVec[i].cur_edge] = selectedEdge;
				antVec[i].weight += edgeArray[selectedEdge].weight;
				antVec[i].nodeVisited[selectedNode] = ANT_GOST;
				mmasAlog->pheromone[selectedEdge] *= 0.9;
				if (head_ptr[selectedNode].vertexType == NODE_TYPE_MID)
				{
					mustPassCount++;
					antVec[i].mustPass_number++;
				}
			}
			if (mustPassCount == subVertex_num&&antVec[i].pathNode[antVec[i].cur_pos] == graph->endVertex)
				break;
		}

	}

}

void findBestAnt(Ant_ptr antVec,acoInfo_ptr acoPtr,graphInfo_ptr g,vertexNode_ptr head_ptr,int choosePath) {
	int subVertex_num = 0;
	if (choosePath == 1)
		subVertex_num = g->subVertex1_num;
	else if (choosePath == 2)
		subVertex_num = g->subVertex2_num;

	int bestAntindex = -1;
	//提取每次迭代蚂蚁构造的可行解中的最优解
	for (int i = 0; i < acoPtr->ant_num; ++i) {
		if (antVec[i].mustPass_number !=subVertex_num)
			continue;
		if (antVec[i].weight < acoPtr->best_pathLen&&antVec[i].pathNode[antVec[i].cur_pos] == g->endVertex) {
			acoPtr->best_pathLen = antVec[i].weight;
			bestAntindex = i;
			acoPtr->mustNum = antVec[i].mustPass_number;
		}
	}
	
	int max_mustPass = acoPtr->mustNum;
	int min_weight=acoPtr->best_pathLen;
	if (bestAntindex == -1) {
		for (int i = 0; i < acoPtr->ant_num; ++i) {

			if (antVec[i].mustPass_number >max_mustPass||
				(antVec[i].mustPass_number==max_mustPass&&antVec[i].weight<min_weight)) {
				max_mustPass = antVec[i].mustPass_number;
				acoPtr->best_pathLen = antVec[i].weight;
				min_weight = antVec[i].weight;
				bestAntindex = i;
			}
		}
	}
	
	if(bestAntindex!=-1){
		if(antVec[bestAntindex].mustPass_number==subVertex_num
			&&antVec[bestAntindex].pathNode[antVec[bestAntindex].cur_pos]==g->endVertex){

			vector<int> temp;
			for(int i=0;i<=antVec[bestAntindex].cur_pos;++i){
				temp.push_back(antVec[bestAntindex].pathNode[i]);
			}

			auto opt_path=optimize_path1(temp,g,head_ptr);
#ifdef DEBUG
			cout<<"best_pathLen before opt:"<<acoPtr->best_pathLen<<endl;
			cout<<"best_pathLen adter opt:"<<opt_path.first<<endl;
#endif
			acoPtr->best_pathLen=opt_path.first;
			acoPtr->mustNum=subVertex_num;
			acoPtr->best_fitval=antVec[bestAntindex].fitness_val;
			acoPtr->pathNodes_num =temp.size();
			for(size_t i=0;i<temp.size();++i)
				acoPtr->bestPathNode[i]=temp[i];
			for(size_t i=temp.size();i<g->node_num;++i)
				acoPtr->bestPathNode[i]=INVALID;
			for(int i=0;i<g->real_edge_num;++i)
				acoPtr->bestPathEdge[i]=INVALID;
			for(size_t i=0;i<temp.size()-1;++i){
				auto iter=vertexTolinkId.find(hash_func(vertex_id[temp[i]],vertex_id[temp[i+1]]));
				acoPtr->bestPathEdge[i]=r_edgeId[iter->second];
			}
		}
		else{
			acoPtr->best_pathLen = antVec[bestAntindex].weight;
			acoPtr->best_fitval = antVec[bestAntindex].fitness_val;
			acoPtr->mustNum = antVec[bestAntindex].mustPass_number;

			acoPtr->pathNodes_num = antVec[bestAntindex].cur_pos + 1;
			for (int i = 0; i < g->node_num; ++i)
				acoPtr->bestPathNode[i] = antVec[bestAntindex].pathNode[i];
			for (int i = 0; i < g->real_edge_num; ++i) {
				acoPtr->bestPathEdge[i] = antVec[bestAntindex].pathEdge[i];
			}

		}
	}
#ifdef DEBUG
	cout<<acoPtr->mustNum<<endl;
#endif

}

void calcFitval(Ant_ptr antVec,acoInfo_ptr acoPtr,graphInfo_ptr g,int choosePath) {
	int subVertex_num = 0;
	if (choosePath == 1)
		subVertex_num = g->subVertex1_num;
	else if (choosePath == 2)
		subVertex_num = g->subVertex2_num;

	double fitnessVal=0.0;
	for (int i = 0; i < acoPtr->ant_num; ++i) {

		if (antVec[i].mustPass_number != subVertex_num)
			continue;
		if (antVec[i].pathNode[antVec[i].cur_pos] == g->endVertex) {
			fitnessVal = (double)antVec[i].weight / (double)acoPtr->best_pathLen;
			antVec[i].fitness_val = fitnessVal;
		}

	}
}
void update_pheromone(Ant_ptr antVec,acoInfo_ptr acoPtr,graphInfo_ptr g,vertexNode_ptr head_ptr,int choosePath) {
	int subVertex_num = 0;
	if (choosePath == 1)
		subVertex_num = g->subVertex1_num;
	else if (choosePath == 2)
		subVertex_num = g->subVertex2_num;
	for (int i = 0; i < g->real_edge_num;++i)
	{
		acoPtr->pheromone[i] = (1 - EVAPORATION_RATE)*acoPtr->pheromone[i];
	}
	double pheromoneTotal = 0.0;
	for (int i = 0; i < acoPtr->ant_num; ++i)
	{
		if (antVec[i].mustPass_number != subVertex_num)
			continue;
		if (antVec[i].pathNode[antVec[i].cur_pos] == g->endVertex) {
			if (antVec[i].fitness_val == 0.0)
				pheromoneTotal = 0.0;
			else
				pheromoneTotal = POSITIVE_CONTS / (antVec[i].fitness_val);
			for (int j = 0; j <=antVec[i].cur_edge; ++j) {
				acoPtr->pheromone[antVec[i].pathEdge[j]] += pheromoneTotal;
			}
		}
	}
	for (int i = 0; i < g->real_edge_num; ++i) {
		if (acoPtr->pheromone[i] > MAX_PHEROMONE)
			acoPtr->pheromone[i] = MAX_PHEROMONE;
		else if (acoPtr->pheromone[i] < MIN_PHEROMONE)
			acoPtr->pheromone[i] = MIN_PHEROMONE;
	}
}

void DisplayBestWay(acoInfo_ptr acoPtr,graphInfo_ptr g,vertexNode_ptr headPtr) {
	cout << "********************************************" << endl;
	cout << "weight=" << acoPtr->best_pathLen<<endl;
	cout << "fitness value=" << acoPtr->best_fitval << endl;
	cout << "contain " << acoPtr->pathNodes_num << " nodes" << endl;
	size_t cnt= 0;
	for (int i = 0; i < acoPtr->pathNodes_num; ++i) {
		if (acoPtr->bestPathNode[i] != INVALID) {
			if (headPtr[acoPtr->bestPathNode[i]].vertexType == NODE_TYPE_MID)
				cnt++;
		}
	}
	cout << "must pass cnt=" << cnt << endl;

	for (int i = 0; i <acoPtr->pathNodes_num; ++i) {
		if (headPtr[acoPtr->bestPathNode[i]].vertexType == NODE_TYPE_MID)
			cout <<"{"<< acoPtr->bestPathNode[i]<<"}" << "->";
		else
			cout << acoPtr->bestPathNode[i] << "->";
	}
	cout << endl;

}

void writeAdjList_toFile(vertexNode_ptr headPtr,int node_num,ofstream &os){
	for(int i=0;i<node_num;++i){
		if(headPtr[i].vertexType==NODE_TYPE_START){
			os<<vertex_id[i]<<" [startNode,adjacentNum="<<headPtr[i].adjacentNum<<"]:";
		}
		else if(headPtr[i].vertexType==NODE_TYPE_END)
			os<<vertex_id[i]<<" [endNode,adjacentNum="<<headPtr[i].adjacentNum<<"]:";
		else if(headPtr[i].vertexType==NODE_TYPE_MID)
			os<<vertex_id[i]<<" [midNode,adjacentNum="<<headPtr[i].adjacentNum<<"]:";
		else
			os<<vertex_id[i]<<" [commNode,adjacentNum="<<headPtr[i].adjacentNum<<"]:";
		auto curPtr=headPtr[i].adjList_head;
		for(int j=0;j<headPtr[i].adjacentNum;++j){
			os<<"("<<vertex_id[curPtr->to]<<","<<curPtr->weight<<","<<edge_id[curPtr->linkId]<<"), ";
			curPtr=curPtr->next;
		}
		os<<endl;
	}
}

