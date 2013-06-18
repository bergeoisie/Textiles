#include "SSEGenerator.h"
#include "textileHelper.h"
#include "output.h"

typedef unordered_map<string,string> specequiv;
typedef graph_traits<Graph> vertex_descriptor;
typedef graph_traits<Graph> edge_descriptor;
typedef property_map<Graph,edge_name_t>::type g_edge_name_type;



SSEGenerator::SSEGenerator(Graph G,Graph H) : ggraph(G),hgraph(H),ghproduct(ProductGraph(G,H)),hgproduct(ProductGraph(H,G))
{}

specequiv SSEGenerator::NextSSE()
{

}

SSETree::SSETree(Graph* gh,Graph* hg): ghprod(gh),hgprod(hg)
{
	SSENode* root = new SSENode(0,0,string(""));


} 

void SSETree::MakeSSE()
{

}

stack<string> SSETree::createAssociationStack(SSENode* node)
{
	set<string> seenEdges;
	stack<string> assocStack;

	int convertedLevel = node->getLevel() - 1;

	SSENode* parent =node.getParent();

	g_edge_name_type hg_e_name = get(edge_name,*hg);

	// 
	ED e = edVector[convertedLevel];
	VD tar = target(e,*gh);

	while(parent != null)
	{
		seenEdges.insert(parent->getAssoc());
		parent = parent->getParent();
	}

	// We want to iterate through the out edges of HG to 
	for(tie(oei,oei_end)=out_edges(source(e,*gh),*hg);oei!=oei_end;oei++)
	{
		if(target(*oei,*hg)==tar && seenEdges.find(hg_e_name(*oei)))
		{
			assocStack.push(hg_e_name(*oei));
		}
	}

	return assocStack;
}

void SSETree::GrowBranch()
{
	bool done=false;
	SSENode* currentNode = root;

	while(!done)
	{
		if(currentNode->NoChildren())
		{
			stack<string> assocStack = createAssociationStack(currentNode);

			while(!assocStack.empty()){
				SSENode* newChild = new SSENode(currentNode,level_+1,currAssoc);
				currentNode->addChild(newChild);
			}
		}
		else{
			auto childrenIt = currentNode->children.begin();
			while(childrenIt->IsBlack())
			{
				++childrenIt;
			}
			if(childrenIt->level() == maxLevel)
			{
				done = true;
			}
			else
			{
				currentNode = *childrenIt;
			}
		}
	}
}

void SSENode::addChild(SSENode* child)
{
	children_.push_back(child);
}

bool SSENode::NoChildren()
{
	return !children_.size();
}

SSENode::SSENode(SSENode* p,int l,string a) : parent_(p), level_(l), assoc_(a),color_(White)
{
	children_ = vector<SSENode*>();
}

const int SSENode::getLevel()
{
	return level_;
}

string SSENode::getAssoc()
{
	return assoc_;
}

SSENode* SSENode::getParent()
{
	return parent_;
}