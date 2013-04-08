//
//  textileHelper.h
//  
//
//  Created by Brendan Berg on 9/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef _textileHelper_h
#define _textileHelper_h


#include <boost/config.hpp> // put this first to suppress some VC++ warnings

#include <iostream>
#include <iterator>
#include <algorithm>
#include <time.h>
#include <fstream>
#include <tr1/unordered_map>
#include <stdlib.h>
#include <queue>

#include <boost/utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace boost;

namespace boost { 
    enum edge_p_homom_t { edge_p_homom };
    enum edge_q_homom_t { edge_q_homom };
    enum vertex_p_vhomom_t { vertex_p_vhomom };
    enum vertex_q_vhomom_t { vertex_q_vhomom };
    
    BOOST_INSTALL_PROPERTY(edge, p_homom);
    BOOST_INSTALL_PROPERTY(edge, q_homom);
    BOOST_INSTALL_PROPERTY(vertex, p_vhomom);
    BOOST_INSTALL_PROPERTY(vertex, q_vhomom);
}
 

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
typedef std::tr1::unordered_map<string,graph_traits<Graph>::vertex_descriptor> VertexMap;
typedef std::tr1::unordered_map<string,Edge> EdgeMap;

typedef pair<GammaGraph,Graph> Textile;

typedef vector<graph_traits<GammaGraph>::vertex_descriptor> VVec;
typedef graph_traits<GammaGraph>::vertex_iterator GammaVI;
typedef graph_traits<GammaGraph>::out_edge_iterator GammaOEI;
typedef graph_traits<GammaGraph>::edge_iterator GammaEI;
typedef graph_traits<Graph>::vertex_iterator GVI;
typedef graph_traits<Graph>::edge_iterator GEI;

typedef graph_traits<Graph>::vertex_descriptor GVD;
typedef graph_traits<GammaGraph>::vertex_descriptor VD;

// A vertex collection is a set of vertex descriptors
typedef set<graph_traits<GammaGraph>::vertex_descriptor> vColl;

typedef tuple<int,int,int> PQOEIElement;

enum colors { White, Gray, Black };


Textile CreateDual(Textile);
bool checkNondegen(Graph);
bool checkHomoms(Textile);
Textile Trim(Textile T);
Textile HigherBlock(Textile);
Textile CreateInverse(Textile);
Textile ProductTextile(Textile,Textile);
bool IsomLanguages(Textile);
bool IsLR(Textile);
bool IspRightResolving(Textile);
bool IspLeftResolving(Textile);
bool IsqRightResolving(Textile);
bool IsqLeftResolving(Textile);
int NDefinite(Textile);
int IsqRightDefinite(Textile);
int IspLeftDefinite(Textile);
int IspRightDefinite(Textile);
int IsqLeftDefinite(Textile);
void PrintRepMatrix(Textile);
void PrintFullTextileInfo(Textile,ostream& os = cout);
void PrintBasicTextileInfo(Textile, ostream& os = cout);
void SmartPrintTextileInfo(Textile, ostream& os = cout);
Textile RenameTextile(Textile,map<string,string>);
Textile HigherNBlock(Textile,int);
Textile Quotient(Textile,vector<VVec>);
Textile HomomComp(Textile, Graph,map<string,string>);
Textile CreateTranspose(Textile);
Textile CreateOneOneTextile(Textile);
VVec compatibleSet(Textile&,
                   bool,
                   VVec,
                   string);
int maxLblOutDegree(Textile,string);
Textile inducedRp(Textile);
Textile inducedLp(Textile);
Textile inducedRq(Textile);
Textile inducedLq(Textile);
bool is1to1(Textile);
bool isOneSided1to1(Textile);
void printVVec(VVec);
bool VVecSubset(VVec&,VVec&);
string ssVVec(VVec);
bool hasCycleHelper(Graph&,GVD,colors*);
void Analyzer(Textile);
bool NewIsomLanguages(Textile);
Textile DirectSum(Textile,Textile);
Textile DFAMinimization(Textile);
Textile Composition(Textile,Textile);
Textile LookForConjugacy(Textile,int, bool*, int start = 2);
Textile CompositionPower(Textile, int);
Textile CreateNMTextile(Textile,int,int);
Textile AutoHomom(Textile);
Textile AutoHomomLite(Textile);
Textile AutoRenamer(Textile);
string Namer(int,int,int=65);
double PFEigenvalue(Graph);
void OctaveOutput(Textile,string);
//Textile RandomTrim(Textile);
Textile ArrayTrim(Textile);
Textile NewInducedRp(Textile);
Textile FromSSE(Graph,Graph,std::tr1::unordered_map<string,string>);
Graph ProductGraph(Graph,Graph);
void PrintGraph(Graph,ostream& os = cout);

#endif
