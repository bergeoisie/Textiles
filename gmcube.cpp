//
//  gmcube.cpp
//  
//
//  Created by Brendan Berg on 9/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "textileHelper.h"

using namespace std;


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

typedef std::tuple<int,int,int> PQOEIElement;



int main(void)
{
    bool found=false;
    int i,j;
    cout << "Defining a GammaGraph" << endl;
    
    
    GammaGraph G(3);
    
    add_edge(0, 0, PQ_Homoms(string("a"),Q_Homom(string("c"),string("01"))),G);
    add_edge(0, 0, PQ_Homoms(string("c"),Q_Homom(string("a"),string("02"))),G);
    add_edge(0, 1, PQ_Homoms(string("a"),Q_Homom(string("e"),string("03"))),G);
    add_edge(0, 1, PQ_Homoms(string("c"),Q_Homom(string("d"),string("04"))),G);
    add_edge(0, 2, PQ_Homoms(string("d"),Q_Homom(string("b"),string("05"))),G);
    add_edge(1, 0, PQ_Homoms(string("b"),Q_Homom(string("g"),string("06"))),G);
    add_edge(1, 1, PQ_Homoms(string("b"),Q_Homom(string("h"),string("07"))),G);
    add_edge(1, 2, PQ_Homoms(string("e"),Q_Homom(string("f"),string("08"))),G);
    add_edge(2, 0, PQ_Homoms(string("f"),Q_Homom(string("c"),string("09"))),G);
    add_edge(2, 0, PQ_Homoms(string("g"),Q_Homom(string("a"),string("10"))),G);
    add_edge(2, 1, PQ_Homoms(string("f"),Q_Homom(string("e"),string("11"))),G);
    add_edge(2, 1, PQ_Homoms(string("g"),Q_Homom(string("d"),string("12"))),G);
    add_edge(2, 2, PQ_Homoms(string("h"),Q_Homom(string("b"),string("13"))),G);

    
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_vhom = get(vertex_p_vhomom,G);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_vhom = get(vertex_q_vhomom,G);
    
    property_map<GammaGraph,vertex_name_t>::type
    vname = get(vertex_name,G);
    
    property_map<GammaGraph,edge_name_t>::type
    name = get(edge_name,G);
    
    put(p_vhom,0,0);
    put(q_vhom,0,0);
    put(vname,0,string("u"));
    put(p_vhom,1,0);
    put(q_vhom,1,1);
    put(vname,1,string("v"));
    put(p_vhom,2,1);
    put(q_vhom,2,0);
    put(vname,2,string("w"));
    
    
    // Set up a file to put the graphviz data in.
    ofstream gviz("gviz.dat");
    write_graphviz(gviz, G, make_label_writer(vname),make_label_writer(name));
    
    cout << "Defining a Graph " << endl;
    
    
    // Define Gr, the target graph
    
    
    Graph Gr(2);
    
    property_map<Graph,vertex_name_t>::type
    grvname = get(vertex_name,Gr);
    
    add_edge(0,0,string("a"),Gr);
    add_edge(0,0,string("b"),Gr);
    add_edge(0,0,string("c"),Gr);
    add_edge(0,1,string("d"),Gr);
    add_edge(0,1,string("e"),Gr);
    add_edge(1,0,string("f"),Gr);
    add_edge(1,0,string("g"),Gr);
    add_edge(1,1,string("h"),Gr);
    
    put(grvname,0,string("D"));
    put(grvname,1,string("E"));
    
    Textile T(G,Gr);
    
    PrintFullTextileInfo(T);
    
    for(i=-1; i>-4;i--)
    {
        for(j=7; j<9; j++)
        {
            Textile Tij = AutoHomomLite(CreateNMTextile(T,i,j));
            Textile Tconj = LookForConjugacy(Tij,5,&found);
            ofstream os("conjugaciescube3212.txt",ios_base::app);
            os << "PRINTING IJ CONJUGACY FOR " << i << j << endl;
            PrintFullTextileInfo(Tij,os);
            PrintFullTextileInfo(Tconj,os);
        }
    }
    
    return 0;
}