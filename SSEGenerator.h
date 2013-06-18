#ifndef _ssegenerator_h
#define _ssegenerator_h

typedef unordered_map<string,string> specequiv;
typedef graph_traits<Graph> vertex_descriptor;
typedef graph_traits<Graph> edge_descriptor;

enum colors { White, Gray, Black };

class SSEGenerator
{
	public: 
		SSGenerator(Graph G, Graph H);
		specequiv nextSSE();

		

	private:
		Graph* ggraph;
		Graph* hgraph;
		Graph ghproduct;
		Graph hgproduct;
		SSETree st;
};

class SSETree
{
public: 
	SSETree(Graph*,Graph*);
	void GrowBranch();


private:
	stack<string> createAssocationStack();
	Graph* ghprod;
	Graph* hgprod;
	SSENode* root;
	int maxLevel;
	vector<ED> edVector;
};

class SSENode
{
public:
	SSENode(SSENode*,int,string);
	const vector<SSENode*> children();
	bool NoChildren();
	void addChild(SSENode*);
	const int getLevel();
	SSENode* getParent();
	string getAssoc();

private:
	SSENode* parent_;
	vector<SSENode*> children_;
	int level_;
	string assoc_;
	colors color_;
};

#endif