#ifndef _ssegenerator_h
#define _ssegenerator_h

typedef unordered_map<string,string> specequiv;

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

private:
	SSENode* parent_;
	vector<SSENode*> children_;
	int level_;
	string assoc_;
	colors color_;
};

#endif