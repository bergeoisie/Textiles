#include "SSEGenerator.h"
#include "textileHelper.h"
#include "output.h"

typedef unordered_map<string,string> specequiv;


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

	int convertedLevel = node->getLevel() - 1;


	// We want to iterate through the out edges of HG to 
	for(tie(oei,oei_end)=out_edges(source(e,*gh),*hg);oei!=oei_end;oei++)
	{
		if(target(*oei,*hg)==tar)
		{

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