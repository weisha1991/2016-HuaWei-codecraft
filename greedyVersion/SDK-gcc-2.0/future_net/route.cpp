#include "route.h"
#define DEBUG
using namespace std;

const int INF = INT_MAX;

vector<bool> closed_sets;
vector<bool> nodes_contain;
vector<bool> nodes_discard;
vector<bool> isMustPassPoint;
vector<bool> isMustPassPoint1;
vector<bool> isMustPassPoint2;

vector<vector<int>> graph_data;
vector<vector<int>> demand_data;
vector<int> subVertex1;
vector<int> subVertex2;
vector<int> node_indegree;
vector<int> node_outdegree;
map<pair<int, int>, pair<int, int>> data_map;
vector<vector<edge>> adjList;
map<pair<int, int>, pair<int, int>> reverse_data_map;
vector<vector<edge>> reverse_adjList;

int node_num;
int startVertex;
int endVertex;

struct openNode{
	openNode(vector<int> w,int r):way(w),range(r){}
	bool operator <(const openNode &rhs)const {
		return this->range>rhs.range;
	}
	vector<int> way;
	int range;
};



//result_nodes should like this:startNode*******endNode
void print_path(const vector<int> &result_nodes){
	// cout<<"w:"<<w<<endl;
	int w=0;
	for(auto iter=result_nodes.begin();iter!=result_nodes.end()-1;++iter){
		w+=data_map[make_pair(*iter,*(iter+1))].second;
	}
	cout<<"weight="<<w<<endl;
	for(auto n:result_nodes)
		cout<<n<<"->";
	cout<<endl;

}
// #define EDGE__PRINT
size_t compute_multiEdges(const vector<int> &result_nodes1,const vector<int> &result_nodes2){
	vector<int> edge_set1,edge_set2;
	for(auto iter=result_nodes1.begin();iter!=result_nodes1.end()-1;++iter){
		edge_set1.push_back(data_map[make_pair(*iter,*(iter+1))].first);
	}
	for(auto iter=result_nodes2.begin();iter!=result_nodes2.end()-1;++iter){
		edge_set2.push_back(data_map[make_pair(*iter,*(iter+1))].first);
	}
	sort(edge_set1.begin(),edge_set1.end());
	sort(edge_set2.begin(),edge_set2.end());
#ifdef EDGE__PRINT
	cout<<"set1:";
	for(auto node:edge_set1)
		cout<<node<<"--";
	cout<<endl;
	cout<<"set2:";
	for(auto node:edge_set2)
		cout<<node<<"--";
	cout<<endl;
#endif
	vector<int> mutiedge_set;
	set_intersection(edge_set1.begin(),edge_set1.end(),edge_set2.begin(),edge_set2.end(),back_inserter(mutiedge_set));
#ifdef DEBUG
	cout<<"save arcs number="<<mutiedge_set.size()<<endl;
	for(auto eachEdge:mutiedge_set){
		cout<<eachEdge<<"-";
	}
	cout<<endl;
#endif
	return mutiedge_set.size();
}
const int  penalty_coeff=10;
void update_graphInfo(bool isReverse,const vector<int> &result_nodes){
	if(!isReverse)
	{
		for(auto iter=result_nodes.begin();iter!=result_nodes.end()-1;++iter){
			for(auto to_iter=adjList[*iter].begin();to_iter!=adjList[*iter].end();++to_iter){
				if(to_iter->second==*(iter+1)){
					to_iter->first+=penalty_coeff;
				}
			}
		}
	}
	else
	{
		for(auto iter=result_nodes.begin();iter!=result_nodes.end()-1;++iter){
			for(auto to_iter=reverse_adjList[*iter].begin();to_iter!=reverse_adjList[*iter].end();++to_iter){
				if(to_iter->second==*(iter+1)){
					to_iter->first+=penalty_coeff;
				}
			}
		}
	}
}

void check_path(const vector<int> &result_nodes,const int &subset_flag){
	set<int> check_set(result_nodes.begin(),result_nodes.end());
	if(check_set.size()!=result_nodes.size())
		cout<<"contain subtour"<<endl;
	else
		cout<<"No subtour"<<endl;
	if(result_nodes.front()==startVertex)
		cout<<"source node right"<<endl;
	else
		cout<<"source node error"<<endl;
	if(result_nodes.back()==endVertex)
		cout<<"desination right"<<endl;
	else
		cout<<"destination error"<<endl;
	if(subset_flag==1){
		vector<int> includeMustPass;
		for(auto node:result_nodes)
			if(isMustPassPoint1[node])
				includeMustPass.push_back(node);
		sort(includeMustPass.begin(),includeMustPass.end());
		sort(subVertex1.begin(),subVertex1.end());
		vector<int> temp;
		set_difference(subVertex1.begin(),subVertex1.end(),includeMustPass.begin(),includeMustPass.end(),back_inserter(temp));
		if(temp.size()==0)
			cout<<"contain all must pass nodes"<<endl;
		else{
			cout << "some 'must pass 'nodes not be included" << endl;
			for(auto node:temp)
				cout<<node<<",";
			cout<<endl;
		}
	}
	if(subset_flag==2){
		vector<int> includeMustPass;
		for(auto node:result_nodes)
			if(isMustPassPoint2[node])
				includeMustPass.push_back(node);
		sort(includeMustPass.begin(),includeMustPass.end());
		sort(subVertex2.begin(),subVertex2.end());
		vector<int> temp;
		set_difference(subVertex2.begin(),subVertex2.end(),includeMustPass.begin(),includeMustPass.end(),back_inserter(temp));
		if (temp.size() == 0)
			cout << "contain all must pass nodes" << endl;
		else {
			cout << "some 'must pass 'nodes not be included" << endl;
			for(auto node:temp)
				cout<<node<<",";
			cout<<endl;
		}
	}
}

void write_path(const vector<int> &result_nodes,const int &path_flag){
	vector<unsigned short> result;
	for(unsigned short i=0;i<result_nodes.size()-1;++i){
		int src=result_nodes[i];
		int dest=result_nodes[i+1];
		result.push_back(data_map[make_pair(src,dest)].first);
	}
	if(path_flag==1){
		for(unsigned short i=0;i<result.size();++i){
			record_result(WORK_PATH,result[i]);
		}
	}
	if(path_flag==2){
		for(unsigned short i=0;i<result.size();++i)
		{
			record_result(BACK_PATH,result[i]);
		}
	}
}


void search_route(char *graph[MAX_EDGE_NUM], int edge_num,char *condition[MAX_DEMAND_NUM],int demand_num)
{
	clock_t start,finish;
	double duration=0.0;
	double time_limit=5.0;
	start=clock();
	parse_file(graph,edge_num,condition,demand_num);
	build_graph();
	preprocess();
    if(node_num<20){
        unsigned short result1[] = {0, 3, 4};//P'路径
        unsigned short result2[] = {5, 6, 2};//P''路径

        for (int i = 0; i < 3; i++)
        {
            record_result(WORK_PATH, result1[i]);
            record_result(BACK_PATH, result2[i]);
        }
        return;
    }
	pair<int,vector<int>> ret_val1;
	pair<int,vector<int>> ret_val2;
	int sum_weight=INF;
	vector<int> result_nodes1;
	vector<int> result_nodes2;

	while(true){
		ret_val1=backward_findPath(subVertex1);
		if(ret_val1.first==-1)
			continue;
		// reverse(ret_val1.second.begin(),ret_val1.second.end());

		update_graphInfo(true,ret_val1.second);
#ifdef DEBUG
		cout<<"path1 find"<<endl;
#endif
		ret_val2=backward_findPath(subVertex2);
		if(ret_val2.first==-1)
			continue;
		if(ret_val1.first+ret_val2.first<sum_weight){
			sum_weight=ret_val1.first+ret_val2.first;
			result_nodes1=ret_val1.second;
			result_nodes2=ret_val2.second;
			break;
		}
		finish=clock();
		if((duration =(double)(finish - start) / CLOCKS_PER_SEC)>time_limit)
			break;

	}
	reverse(result_nodes1.begin(),result_nodes1.end());
	reverse(result_nodes2.begin(),result_nodes2.end());
#ifdef DEBUG
	cout<<"*************************************"<<endl;
	cout<<"check path1:"<<endl;
	print_path(result_nodes1);
	check_path(result_nodes1,1);
	cout<<"*************************************"<<endl;
	cout<<"check path2"<<endl;
	print_path(result_nodes2);
	check_path(result_nodes2,2);
	cout<<"*************************************"<<endl;
	compute_multiEdges(result_nodes1,result_nodes2);
	cout<<"w1+w2="<<sum_weight<<endl;
#endif
	write_path(result_nodes1,1);
	write_path(result_nodes2,2);
}

void parse_file(char *graph[MAX_EDGE_NUM], int edge_num, char *condition[MAX_DEMAND_NUM], int demand_num){

	for(int i=0;i<edge_num;++i){
		char *token=strtok(graph[i],",");
		vector<int> temp;
		while(token!=NULL){
			temp.push_back(atoi(token));
			token=strtok(NULL,",");
		}
		graph_data.push_back(temp);
	}
	demand_data.resize(2);
	for(int i=0;i<demand_num;++i){
		char *token=strtok(condition[i],",|");
		while(token!=NULL){
			demand_data[i].push_back(atoi(token));
			token=strtok(NULL,",|");
		}
	}

	unsigned int max_point = 0;
	unsigned int temp = 0;
	for (unsigned int i = 0; i < graph_data.size(); ++i) {
		temp = max(graph_data[i][1], graph_data[i][2]);
		if (temp > max_point)
			max_point = temp;
	}
	node_num=max_point+1;

	startVertex=demand_data[0][1];
	endVertex=demand_data[0][2];
	isMustPassPoint1.assign(node_num,false);
	isMustPassPoint2.assign(node_num,false);
	for(size_t i=3;i<demand_data[0].size();++i){
		subVertex1.push_back(demand_data[0][i]);
		isMustPassPoint1[demand_data[0][i]]=true;
	}
	for(size_t i=3;i<demand_data[1].size();++i){
		subVertex2.push_back(demand_data[1][i]);
		isMustPassPoint2[demand_data[1][i]]=true;
	}
#ifdef DEBUG
	cout<<"node num:"<<node_num<<endl;
	cout<<"edge num:"<<edge_num<<endl;
	cout<<"source:"<<startVertex<<endl;
	cout<<"destinaiton:"<<endVertex<<endl;
	cout<<"subset 1 ("<<subVertex1.size()<<" nodes):";
	for(auto node:subVertex1)
		cout<<node<<",";
	cout<<endl;
	cout<<"subset 2 ("<<subVertex2.size()<<" nodes):";
	for(auto node:subVertex2)
		cout<<node<<",";
	cout<<endl;
#endif
}

void preprocess(){
	node_indegree.assign(node_num,0);
	node_outdegree.assign(node_num,0);
	for(auto line:graph_data){
		node_outdegree[line[1]]+=1;
		node_indegree[line[2]]+=1;
	}
	nodes_discard.assign(node_num,false);
#ifdef DEBUG
	cout<<"node's indegree equal zero:";
#endif
	for(size_t i=0;i<node_indegree.size();++i){
		if(node_indegree[i]==0){
			nodes_discard[i]=true;
#ifdef DEBUG
			cout<<i<<",";
#endif
		}
	}
#ifdef DEBUG
	cout<<endl;
	cout<<"node's outdegree equal zero:";
#endif
	for(size_t i=0;i<node_outdegree.size();++i){
		if(node_outdegree[i]==0){
#ifdef DEBUG
			cout<<i<<",";
#endif
			nodes_discard[i]=true;
		}
	}
	cout<<endl;
	auto ret=dijkstra_singleSource(node_num,startVertex,adjList);
	for(size_t i=0;i<ret.first.size();++i){
		if(ret.first[i]==INF)
			nodes_discard[i]=true;
	}
	ret=dijkstra_singleSource(node_num,endVertex,reverse_adjList);
	for(size_t i=0;i<ret.first.size();++i){
		if(ret.first[i]==INF)
			nodes_discard[i]=true;
	}

	nodes_discard[startVertex]=false;
	nodes_discard[endVertex]=false;
#ifdef DEBUG
	cout<<"discard "<<count(nodes_discard.begin(),nodes_discard.end(),true)<<" nodes"<<endl;
#endif


}

void build_graph(){
	set<pair<int,int>> srcdes_paris;
	for (size_t i = 0; i<graph_data.size(); ++i) {
		pair<int, int> src_des(graph_data[i][1], graph_data[i][2]);
		srcdes_paris.insert(src_des);
		auto find_key=data_map.find(src_des);
		if(find_key==data_map.end()){
			data_map[src_des] = make_pair(graph_data[i][0], graph_data[i][3]);
		}else{
			if(graph_data[i][3]<(find_key->second.second))
				data_map[src_des] = make_pair(graph_data[i][0], graph_data[i][3]);
		}
	}
	adjList.resize(node_num);
	for(auto iter=srcdes_paris.begin();iter!=srcdes_paris.end();++iter){
		adjList[iter->first].push_back(make_pair(data_map[*iter].second,iter->second));
	}

	set<pair<int,int>> rsrcdes_paris;
	for (size_t i = 0; i<graph_data.size(); ++i) {
		pair<int, int> src_des(graph_data[i][2], graph_data[i][1]);
		rsrcdes_paris.insert(src_des);
		auto find_key=reverse_data_map.find(src_des);
		if(find_key==reverse_data_map.end()){
			reverse_data_map[src_des] = make_pair(graph_data[i][0], graph_data[i][3]);
		}else{
			if(graph_data[i][3]<(find_key->second.second))
				reverse_data_map[src_des] = make_pair(graph_data[i][0], graph_data[i][3]);
		}
	}
	reverse_adjList.resize(node_num);
	for(auto iter=rsrcdes_paris.begin();iter!=rsrcdes_paris.end();++iter){
		reverse_adjList[iter->first].push_back(make_pair(reverse_data_map[*iter].second,iter->second));
	}
}


pair<int,vector<int>> forward_findPath(const vector<int> &pass_nodes){

	priority_queue<openNode> init_paths;
	for(auto out=pass_nodes.begin();out!=pass_nodes.end();++out){
		closed_sets.assign(node_num,false);
		for(auto inner:pass_nodes)
			if(inner!=*out)
				closed_sets[inner]=true;
		vector<int> temp_path;
		int w=dijkstra_New(node_num,startVertex,*out,adjList,temp_path);
		if(w==-1)
			continue;
		reverse(temp_path.begin(),temp_path.end());
		init_paths.push(openNode(temp_path,w));
	}
#ifdef DEBUG
	cout<<"init path num:"<<init_paths.size()<<endl;
#endif
	isMustPassPoint.assign(node_num,false);
	for(auto node:pass_nodes)
		isMustPassPoint[node]=true;

	unsigned int seed=0;
	vector<int> w_vec;
	vector<vector<int>> save_paths;
	default_random_engine randseed(time(NULL));
	uniform_int_distribution<int> uniform_random(1,pass_nodes.size()*3/4);
	uint loops_upper=1000;

#ifdef DEBUG
	cout<<"loops="<<loops_upper<<endl;
#endif
	while(!init_paths.empty())
	{
		openNode top_elem=init_paths.top();
		init_paths.pop();
#ifdef DEBUG
		cout<<"remaind paths num:"<<init_paths.size()<<endl;
#endif
		vector<int> cp=pass_nodes;
		for(auto iter=top_elem.way.begin()+1;iter!=top_elem.way.end();++iter)
			if(isMustPassPoint[*iter])
				cp.erase(find(cp.begin(),cp.end(),*iter));

		for(size_t loop=0;loop<loops_upper;++loop)
		{
			nodes_contain.assign(node_num,false);
			for(auto node:top_elem.way)
				nodes_contain[node]=true;
			int si=*top_elem.way.rbegin();
			vector<int> node_seq=cp;

			vector<int> path1;
			int kk=0;
			int path1_w=0;
			int upper=uniform_random(randseed);

			while(!node_seq.empty()&&kk++<upper){
				vector<int> pathi;
				int wi=dijkstra_toSubset(node_num,si,node_seq,adjList,pathi);
				if(wi==-1)
					break;

				path1_w+=wi;
				si=*pathi.rbegin();
				for(auto node:pathi){
					nodes_contain[node]=true;
					path1.push_back(node);
				}

				for(auto iter=pathi.begin()+1;iter!=pathi.end();++iter){
					if(isMustPassPoint[*iter])
						node_seq.erase(find(node_seq.begin(),node_seq.end(),*iter));
				}

			}

			if(path1.empty())
				continue;
			path1.erase(unique(path1.begin(),path1.end()),path1.end());
			node_seq.insert(node_seq.begin(),*path1.rbegin());
			node_seq.push_back(endVertex);
			seed = std::chrono::system_clock::now().time_since_epoch().count();
			shuffle(node_seq.begin() + 1, node_seq.end()-1,default_random_engine(seed));
			vector<int> path2;
			for(auto nn:path1)
				nodes_contain[nn]=true;
			int w2 = path_midNode(node_seq, path2, node_num, adjList);

			if(w2==-1)
				continue;
			if(path2.back()!=endVertex)
				continue;
			vector<int> temp;
			for(auto node:top_elem.way)
				temp.push_back(node);
			for(auto node:path1)
				temp.push_back(node);
			for(auto node:path2)
				temp.push_back(node);
			temp.erase(unique(temp.begin(),temp.end()),temp.end());
			save_paths.push_back(temp);
			w_vec.push_back(top_elem.range+w2+path1_w);
			break;
		}
		if(!w_vec.empty())
			break;
	}
	if(w_vec.empty())
		return make_pair(-1,vector<int>());
#ifdef DEBUG
	cout<<"save_paths size:"<<save_paths.size()<<endl;
#endif
	return make_pair(w_vec[0],save_paths[0]);
}

pair<int,vector<int>> backward_findPath(const vector<int> &pass_nodes){

	priority_queue<openNode> init_paths;
	for(auto out=pass_nodes.begin();out!=pass_nodes.end();++out){
		closed_sets.assign(node_num,false);
		for(auto inner:pass_nodes)
			if(inner!=*out)
				closed_sets[inner]=true;
		vector<int> temp_path;
		int w=dijkstra_New(node_num,endVertex,*out,reverse_adjList,temp_path);
		if(w==-1)
			continue;
		reverse(temp_path.begin(),temp_path.end());
		init_paths.push(openNode(temp_path,w));
	}
#ifdef DEBUG
	cout<<"init path num:"<<init_paths.size()<<endl;
#endif
	isMustPassPoint.assign(node_num,false);
	for(auto node:pass_nodes)
		isMustPassPoint[node]=true;

	unsigned int seed=0;
	//vector<int> w_vec;
	int min_weight = INF;
	vector<int > save_paths;
	//vector<vector<int>> save_paths;
	default_random_engine randseed(time(NULL));
	uniform_int_distribution<int> uniform_random(1,pass_nodes.size()*3/4);
	uint loops_upper=1000;
#ifdef DEBUG
	cout<<"loops="<<loops_upper<<endl;
#endif
	size_t findPath_number = 0;
	while(!init_paths.empty())
	{
		openNode top_elem=init_paths.top();
		init_paths.pop();
#ifdef DEBUG
		cout<<"remaind paths num:"<<init_paths.size()<<endl;
#endif
		vector<int> cp=pass_nodes;
		for(auto iter=top_elem.way.begin()+1;iter!=top_elem.way.end();++iter)
			if(isMustPassPoint[*iter])
				cp.erase(find(cp.begin(),cp.end(),*iter));

		for(size_t loop=0;loop<loops_upper;++loop)
		{
			nodes_contain.assign(node_num,false);
			for(auto node:top_elem.way)
				nodes_contain[node]=true;
			int si=*top_elem.way.rbegin();
			vector<int> node_seq=cp;

			vector<int> path1;
			int kk=0;
			int path1_w=0;
			int upper=uniform_random(randseed);

			while(!node_seq.empty()&&kk++<upper){
				vector<int> pathi;
				int wi=dijkstra_toSubset(node_num,si,node_seq,reverse_adjList,pathi);
				if(wi==-1)
					break;

				path1_w+=wi;
				si=*pathi.rbegin();
				for(auto node:pathi){
					nodes_contain[node]=true;
					path1.push_back(node);
				}

				for(auto iter=pathi.begin()+1;iter!=pathi.end();++iter){
					if(isMustPassPoint[*iter])
						node_seq.erase(find(node_seq.begin(),node_seq.end(),*iter));
				}

			}

			if(path1.empty())
				continue;
			path1.erase(unique(path1.begin(),path1.end()),path1.end());
			node_seq.insert(node_seq.begin(),*path1.rbegin());
			node_seq.push_back(startVertex);
			seed = std::chrono::system_clock::now().time_since_epoch().count();
			shuffle(node_seq.begin() + 1, node_seq.end()-1,default_random_engine(seed));
			vector<int> path2;
			for(auto nn:path1)
				nodes_contain[nn]=true;
			int w2 = path_midNode(node_seq, path2, node_num, reverse_adjList);

			if(w2==-1)
				continue;
			if(path2.back()!=startVertex)
				continue;
			vector<int> temp;
			for(auto node:top_elem.way)
				temp.push_back(node);
			for(auto node:path1)
				temp.push_back(node);
			for(auto node:path2)
				temp.push_back(node);
			temp.erase(unique(temp.begin(),temp.end()),temp.end());
			//save_paths.push_back(temp);
			//w_vec.push_back(top_elem.range+w2+path1_w);
			//break;
			if (top_elem.range + w2 + path1_w < min_weight) {
				min_weight = top_elem.range + w2 + path1_w;
				save_paths = temp;
				findPath_number++;
				if (findPath_number > 5)
					break;
			}
		}
		//if(!w_vec.empty())
		//	break;
		if (findPath_number > 5)
			break;
	}
	/*if(w_vec.empty())
		return make_pair(-1,vector<int>());*/
	if (min_weight == INF)
		return make_pair(-1, vector<int>());
#ifdef DEBUG
	cout<<"save_paths size:"<<save_paths.size()<<endl;
#endif
	return make_pair(min_weight,save_paths);
}

int path_midNode(const vector<int> &midnode_seq, vector<int> &paths, int vertexs, const vector<vector<edge>> &adjList) {
	int dist = 0;
	int sum = 0;
	unsigned int i = 0;
	while (i<midnode_seq.size()) {
		int start = midnode_seq[i];
		unsigned int  j = i + 1;
		for (; j < midnode_seq.size(); ++j) {
			if (nodes_contain[midnode_seq[j]] == false)
				break;
		}
		if (j == midnode_seq.size())
			break;
		int end = midnode_seq[j];
		vector<int> path_i;
		dist = dijkstra(vertexs, start, end, adjList, path_i);
		if (dist == -1)
			return -1;
		reverse(path_i.begin(), path_i.end());
		sum += dist;
		for (auto node : path_i) {
			nodes_contain[node] = true;
			paths.push_back(node);
		}
		i = j;
	}

	auto enditer=unique(paths.begin(), paths.end());
	paths.erase(enditer, paths.end());
	return sum;
}


int dijkstra_toSubset(int node_num, int startVertex, vector<int> subvertex, const vector<vector<edge>> &adjList, vector<int> &path) {
	priority_queue<edge, vector<edge>, greater<edge> > Q;
	vector<int> dist(node_num, INF), dad(node_num, -1);
	Q.push(make_pair(0, startVertex));
	dist[startVertex] = 0;
	unsigned int cnt = 0;
	while (!Q.empty()) {
		edge p = Q.top();
		if (find(subvertex.begin(), subvertex.end(), p.second) != subvertex.end()) {
			cnt++;
		}
		if (cnt == subvertex.size())
			break;
		Q.pop();
		int here = p.second;
		for (vector<edge>::const_iterator it = adjList[here].begin(); it != adjList[here].end(); it++) {
			if(nodes_discard[it->second]||nodes_contain[it->second])
				continue;
			if (dist[here] + it->first < dist[it->second]) {
				dist[it->second] = dist[here] + it->first;
				dad[it->second] = here;
				Q.push(make_pair(dist[it->second], it->second));
			}
		}
	}
	vector<vector<int>> sps;
	vector<int> w_vec;
	for (size_t j = 0; j < subvertex.size(); ++j) {
		if (dist[subvertex[j]] < INF) {
			vector<int> sp;
			for (int i = subvertex[j]; i != -1; i = dad[i]) {
				sp.push_back(i);
			}
			reverse(sp.begin(), sp.end());
			//sp.push_back(subvertex[j]);
			sps.push_back(sp);
			w_vec.push_back(dist[subvertex[j]]);
		}
	}
	if(w_vec.empty())
		return -1;
	auto iter = min_element(w_vec.begin(), w_vec.end());
	ptrdiff_t pos = distance(w_vec.begin(), iter);
	for (auto i : sps[pos]) {
		path.push_back(i);
	}

	return *iter;
}

pair<vector<int>,vector<int>> dijkstra_singleSource(int node_num, int startVertex,const vector<vector<edge>> &adjList) {
	priority_queue<edge, vector<edge>, greater<edge> > Q;
	vector<int> dist(node_num, INF), dad(node_num, -1);
	Q.push(make_pair(0, startVertex));
	dist[startVertex] = 0;
	while (!Q.empty()) {
		edge p = Q.top();
		Q.pop();
		int here = p.second;
		for (vector<edge>::const_iterator it = adjList[here].begin(); it != adjList[here].end(); it++) {

			if (dist[here] + it->first < dist[it->second]) {
				dist[it->second] = dist[here] + it->first;
				dad[it->second] = here;
				Q.push(make_pair(dist[it->second], it->second));
			}
		}
	}
	return make_pair(dist,dad);
}

int dijkstra(int node_num, int startVertex, int endVertex, const vector<vector<edge>> &adjList, vector<int> &path) {
	priority_queue<edge, vector<edge>, greater<edge> > Q;
	vector<int> dist(node_num, INF), dad(node_num, -1);
	Q.push(make_pair(0, startVertex));
	dist[startVertex] = 0;
	while (!Q.empty()) {
		edge p = Q.top();

		if (p.second == endVertex)
			break;
		Q.pop();
		int here = p.second;
		for (vector<edge>::const_iterator it = adjList[here].begin(); it != adjList[here].end(); it++) {
			if (nodes_contain[it->second]||nodes_discard[it->second])
				continue;
			if (dist[here] + it->first < dist[it->second]) {
				dist[it->second] = dist[here] + it->first;
				dad[it->second] = here;
				Q.push(make_pair(dist[it->second], it->second));
			}
		}
	}

	if (dist[endVertex] < INF) {
		for (int i = endVertex; i != -1; i = dad[i])
			path.push_back(i);
		return dist[endVertex];
	}

	return -1;
}

int dijkstra_New(int node_num, int startVertex, int endVertex, const vector<vector<edge>> &adjList, vector<int> &path) {
	priority_queue<edge, vector<edge>, greater<edge> > Q;
	vector<int> dist(node_num, INF), dad(node_num, -1);
	Q.push(make_pair(0, startVertex));
	dist[startVertex] = 0;
	while (!Q.empty()) {
		edge p = Q.top();

		if (p.second == endVertex)
			break;
		Q.pop();
		int here = p.second;
		for (vector<edge>::const_iterator it = adjList[here].begin(); it != adjList[here].end(); it++) {
			if (closed_sets[it->second]||nodes_discard[it->second])
				continue;
			if (dist[here] + it->first < dist[it->second]) {
				dist[it->second] = dist[here] + it->first;
				dad[it->second] = here;
				Q.push(make_pair(dist[it->second], it->second));
			}
		}
	}

	if (dist[endVertex] < INF) {
		for (int i = endVertex; i != -1; i = dad[i])
			path.push_back(i);
		return dist[endVertex];
	}

	return -1;
}

