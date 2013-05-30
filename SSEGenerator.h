typedef unordered_map<string,string> specequiv;

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

};