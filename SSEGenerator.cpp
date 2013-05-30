#include "SSEGenerator.h"
#include "textileHelper.h"
#include "output.h"

typedef unordered_map<string,string> specequiv;

SSEGenerator::SSEGenerator(Graph G,Graph H) : ggraph(G),hgraph(H),ghproduct(ProductGraph(G,H)),hgproduct(ProductGraph(H,G))
{}

specequiv SSEGenerator::NextSSE()
{

}