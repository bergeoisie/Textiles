#include "output.h"
#include "textileHelper.h"


#include <iostream>

#include <time.h>


#include <boost/config.hpp> // put this first to suppress some VC++ warnings

#include <iostream>
#include <iterator>
#include <algorithm>
#include <time.h>
#include <fstream>
#include <unordered_map>
#include <stdlib.h>
#include <queue>
#include <tuple>

#include <boost/utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/transpose_graph.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

// A string collection is a set of strings which are the names of applicable vertices
typedef vector<string> sColl;

typedef pair<int,string> sPair;

// Typedef an edge as a pair of integers.
typedef pair<int,int> Edge;

// Set up all the p,q homomorphism properties
typedef property<edge_name_t, string> Name;
typedef property<vertex_name_t,string> VName;
typedef property<edge_q_homom_t, string, Name> Q_Homom;
typedef property<edge_p_homom_t, string, Q_Homom> PQ_Homoms;
typedef property<vertex_q_vhomom_t, int, VName> Q_VHomom;
typedef property<vertex_p_vhomom_t, int, Q_VHomom> PQ_VHomoms;


// Typedef a graph, using the STL vector structure. This is more
// space efficient than the STL list structure.
typedef adjacency_list<vecS,vecS, bidirectionalS, VName, Name> Graph;

typedef adjacency_list<vecS,vecS, bidirectionalS, PQ_VHomoms,
PQ_Homoms> GammaGraph;

// Define two types of unordered_maps
typedef std::unordered_map<string,graph_traits<Graph>::vertex_descriptor> VertexMap;
typedef std::unordered_map<string,Edge> EdgeMap;

typedef pair<GammaGraph,Graph> Textile;

typedef vector<graph_traits<GammaGraph>::vertex_descriptor> VVec;
typedef graph_traits<GammaGraph>::vertex_iterator GammaVI;
typedef graph_traits<GammaGraph>::out_edge_iterator GammaOEI;
typedef graph_traits<Graph>::out_edge_iterator GOEI;
typedef graph_traits<GammaGraph>::edge_iterator GammaEI;
typedef graph_traits<Graph>::vertex_iterator GVI;
typedef graph_traits<Graph>::edge_iterator GEI;

typedef graph_traits<Graph>::vertex_descriptor GVD;
typedef graph_traits<GammaGraph>::vertex_descriptor VD;

// A vertex collection is a set of vertex descriptors
typedef set<graph_traits<GammaGraph>::vertex_descriptor> vColl;

typedef std::tuple<int,int,int> PQOEIElement;

void PrintRepMatrix(Textile T)
{
    bool found;
    
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end,wi,wi_end;
    graph_traits<GammaGraph>::edge_descriptor e;
    
    graph_traits<GammaGraph>::out_edge_iterator oei,oei_end;
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,T.first);
    
    
    
    
    for(tie(vi,vi_end)=vertices(T.first); vi!=vi_end; vi++)
    {
        map<graph_traits<GammaGraph>::vertex_descriptor,string> strStor;
        for(tie(oei,oei_end)=out_edges(*vi,T.first);oei!=oei_end;oei++)
        {
            graph_traits<GammaGraph>::vertex_descriptor t=target(*oei,T.first);
            if(strStor.find(t)==strStor.end()) {
                strStor[t]=p_homom(*oei) + "/" + q_homom(*oei);
            }
            else {
                strStor[t]=strStor[t]+ string(" + ") + p_homom(*oei) + string("/") + q_homom(*oei);
            }
            
        }
        
        for(tie(wi,wi_end)=vertices(T.first); wi!=wi_end; wi++)
        {
            
            if(strStor.find(*wi)!=strStor.end()) {
                cout << strStor[*wi] << " ";
            }
            else {
                cout << "0 ";
            }
        }
        cout << endl;
        
    }
    
}

void PrintFullTextileInfo(Textile T,ostream& os)
{
    bool found;
    
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end,wi,wi_end;
    graph_traits<GammaGraph>::edge_descriptor e;
    graph_traits<GammaGraph>::out_edge_iterator oei,oei_end;
    
    graph_traits<Graph>::vertex_iterator gvi,gvi_end,gwi,gwi_end;
    graph_traits<Graph>::out_edge_iterator goei,goei_end;
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,T.first);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    pvhom = get(vertex_p_vhomom,T.first);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    qvhom = get(vertex_q_vhomom,T.first);
    
    property_map<GammaGraph,edge_name_t>::type
    ename = get(edge_name,T.first);
    
    property_map<GammaGraph,vertex_name_t>::type
    vname = get(vertex_name,T.first);
    
    property_map<Graph,edge_name_t>::type
    gename = get(edge_name,T.second);
    
    property_map<Graph,vertex_name_t>::type
    gvname = get(vertex_name,T.second);
    
    os << "Printing Gamma Graph which has " << num_vertices(T.first) << " vertices and " << num_edges(T.first) << " edges." << endl;
    
    for(tie(vi,vi_end)=vertices(T.first); vi!=vi_end; vi++)
    {
        os << vname(*vi) << " ";
        map<graph_traits<GammaGraph>::vertex_descriptor,string> strStor;
        for(tie(oei,oei_end)=out_edges(*vi,T.first);oei!=oei_end;oei++)
        {
            graph_traits<GammaGraph>::vertex_descriptor t=target(*oei,T.first);
            if(strStor.find(t)==strStor.end()) {
                strStor[t]=ename(*oei);
            }
            else {
                strStor[t]=strStor[t]+ string(" + ") + ename(*oei);
            }
            
        }
        
        for(tie(wi,wi_end)=vertices(T.first); wi!=wi_end; wi++)
        {
            
            if(strStor.find(*wi)!=strStor.end()) {
                os << strStor[*wi] << " ";
            }
            else {
                os << "0 ";
            }
        }
        os << endl;
        
    }
    
    os << "Printing Representation Matrix" << endl;
    
    for(tie(vi,vi_end)=vertices(T.first); vi!=vi_end; vi++)
    {
        cout << vname(*vi) << " ";
        map<graph_traits<GammaGraph>::vertex_descriptor,string> strStor;
        for(tie(oei,oei_end)=out_edges(*vi,T.first);oei!=oei_end;oei++)
        {
            graph_traits<GammaGraph>::vertex_descriptor t=target(*oei,T.first);
            if(strStor.find(t)==strStor.end()) {
                strStor[t]=p_homom(*oei) + "/" + q_homom(*oei);
            }
            else {
                strStor[t]=strStor[t]+ string(" + ") + p_homom(*oei) + string("/") + q_homom(*oei);
            }
            
        }
        
        for(tie(wi,wi_end)=vertices(T.first); wi!=wi_end; wi++)
        {
            
            if(strStor.find(*wi)!=strStor.end()) {
                os << strStor[*wi] << " ";
            }
            else {
                os << "0 ";
            }
        }
        os << endl;
        
    }
    
    os << "Printing vertex homoms" << endl;
    
    for(tie(vi,vi_end)=vertices(T.first); vi!=vi_end; vi++)
    {
        os << vname(*vi) << " " << pvhom(*vi) << " " << qvhom(*vi) << endl;
        
    }
    
    
    os << "Printing G Graph which has " << num_vertices(T.second) << " vertices and " << num_edges(T.second) << " edges." << endl;
    
    
    for(tie(gvi,gvi_end)=vertices(T.second); gvi!=gvi_end; gvi++)
    {
        os << gvname(*gvi) << " ";
        map<graph_traits<Graph>::vertex_descriptor,string> strStor;
        for(tie(goei,goei_end)=out_edges(*gvi,T.second);goei!=goei_end;goei++)
        {
            graph_traits<Graph>::vertex_descriptor t=target(*goei,T.second);
            if(strStor.find(t)==strStor.end()) {
                strStor[t]=gename(*goei);
            }
            else {
                strStor[t]=strStor[t]+ string(" + ") + gename(*goei);
            }
            
        }
        
        for(tie(gwi,gwi_end)=vertices(T.second); gwi!=gwi_end; gwi++)
        {
            
            if(strStor.find(*gwi)!=strStor.end()) {
                os << strStor[*gwi] << " ";
            }
            else {
                os << "0 ";
            }
        }
        os << endl;
        
    }
    
}

void PrintBasicTextileInfo(Textile T,ostream& os)
{
	
	os << "Gamma Graph has " << num_vertices(T.first) << " vertices and " << num_edges(T.first) << " edges." << endl;
	os << "G Graph which has " << num_vertices(T.second) << " vertices and " << num_edges(T.second) << " edges." << endl;
    
}

void SmartPrintTextileInfo(Textile T, ostream& os)
{
	if(num_edges(T.first) > 500)
	{
		cout << "Textile is large, just printing the basic info." << endl;
		PrintBasicTextileInfo(T,os);
	}
	else{
		PrintFullTextileInfo(T,os);
	}
}

void printVVec(VVec a)
{
    VVec::iterator it;
    cout << "{";
    
    for(it = a.begin(); it < a.end()-1; it++)
    {
        cout << *it << " ";
    }
    cout << *it << "}";
}

void printVColl(vColl a)
{
    vColl::iterator it;
    cout << "{";
    
    for(it = a.begin(); it < a.end()-1; it++)
    {
        cout << *it << " ";
    }
    cout << *it << "}";
}

// Makes a VVec (a vector of vertex descriptors) into a stringstream
string ssVVec(VVec a)
{
    stringstream ss;
    VVec::iterator it;
    ss << "{";
    
    for(it = a.begin(); it < a.end()-1; it++)
    {
        ss << *it << " ";
    }
    ss << *it << "}";
    
    return ss.str();
}

// This is a helper function that produces shortened strings based on
// a given length. You can also change the starting point.
string Namer(int i, int len, int start)
{
    int div,j,curr;
    string name;
    char * temp = new char[len+1];
    
    for(j=0;j<len;j++)
    {
        div = i / pow(26.0,j);
        if(div > 0)
        {
            curr = div % 26;
        }
        else {
            curr = 0;
        }
        temp[j]= (char) (start + curr);
    }
	temp[len] = '\0';
    name.append(temp);
    
    return name;
}
// This function is really designed for outputting the G graph of a textile so that you input
// the matrix into Octave and 
void OctaveOutput(Textile T,string filename)
{
	ofstream os(filename.c_str(),ios_base::trunc);
	
	os << "[" ;
	
	GVI vi,vi_end,wi,wi_end;
	GOEI oei,oei_end;
	int i,N = num_vertices(T.second);
	
	for(tie(vi,vi_end)=vertices(T.second); vi!=vi_end; vi++)
	{
		// degStor stores a vertex of G and an integer representing the degree of the edges
		// between *vi and the GVD.
		map<GVD,int> degStor;
		for(tie(oei,oei_end)=out_edges(*vi,T.second); oei!=oei_end; oei++)
		{
			GVD t=target(*oei,T.second);
			if(degStor.find(t)==degStor.end()) // we didn't find it
			{
				degStor[t] = 1;
			}
			else
			{
				degStor[t] = degStor[t]+1;
			}
		}
		
		for(tie(wi,wi_end)=vertices(T.second), i=0; wi!=wi_end; wi++, i++)
		{
			if(degStor.find(*wi)==degStor.end())
			{
				os << "0";
			}
			else
			{
				os << degStor[*wi];
			}
			if(i!=N-1)
			{
				os << ", ";
			}
			else
			{
				os << ";" << endl;
			}
		}
		
	}
	
	os << "]";
}

void PrintGraph(Graph G,ostream& os)
{
    GVI gvi, gvi_end,gwi,gwi_end;
    GOEI goei,goei_end;


    property_map<Graph,edge_name_t>::type
    G_ename = get(edge_name,G);

    property_map<Graph,vertex_name_t>::type
    G_vname = get(vertex_name,G);    

    os << "Printing G Graph which has " << num_vertices(G) << " vertices and " << num_edges(G) << " edges." << endl;
    
    
    for(tie(gvi,gvi_end)=vertices(G); gvi!=gvi_end; gvi++)
    {
        os << G_vname(*gvi) << " ";
        map<graph_traits<Graph>::vertex_descriptor,string> strStor;
        for(tie(goei,goei_end)=out_edges(*gvi,G);goei!=goei_end;goei++)
        {
            graph_traits<Graph>::vertex_descriptor t=target(*goei,G);
            if(strStor.find(t)==strStor.end()) {
                strStor[t]=G_ename(*goei);
            }
            else {
                strStor[t]=strStor[t]+ string(" + ") + G_ename(*goei);
            }
            
        }
        
        for(tie(gwi,gwi_end)=vertices(G); gwi!=gwi_end; gwi++)
        {
            
            if(strStor.find(*gwi)!=strStor.end()) {
                os << strStor[*gwi] << " ";
            }
            else {
                os << "0 ";
            }
        }
        os << endl;
        
    }
}


vector<string> StringSplitter(string s,vector<int> sizes)
{
 
  int i,j=0,n=sizes.size(),length = s.length();
  vector<string> split(n);

//  cout << "Initial String: " << s << endl;

  for(i=0;i<n;i++)
    {
      split[i] = s.substr(j,sizes[i]);
 //     cout << "The "<< i << "th chunk is " << split[i] << " while sizes[" << i << "] is " << sizes[i] << endl; 
      j += sizes[i];
    }
  

    return split;
}
