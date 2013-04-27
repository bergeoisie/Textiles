//
//  gmcube.cpp
//  
//
//  Created by Brendan Berg on 9/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <algorithm>
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
typedef std::unordered_map<string,graph_traits<Graph>::vertex_descriptor> VertexMap;
typedef std::unordered_map<string,Edge> EdgeMap;

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
    int i;
	Graph GM(2),HM(2),GH,HG;
	std::unordered_map<string,string> sequiv;

    property_map<Graph,vertex_name_t>::type
    GM_vname = get(vertex_name,GM); 

    property_map<Graph,vertex_name_t>::type
    HM_vname = get(vertex_name,HM); 

    property_map<Graph,edge_name_t>::type
    GH_ename = get(edge_name,GH);

    property_map<Graph,edge_name_t>::type
    HG_ename = get(edge_name,HG);

    GEI ei,ei_end,fi,fi_end;

    add_edge(0,0,string("u"),GM);
    add_edge(0,1,string("v"),GM);
    add_edge(1,0,string("w"),GM);
    
    put(GM_vname,0,string("A"));
    put(GM_vname,1,string("B"));


    PrintGraph(GM);

    add_edge(0,0,string("x"),HM);
    add_edge(0,1,string("y"),HM);
    add_edge(1,0,string("z"),HM);
    
    put(HM_vname,0,string("C"));
    put(HM_vname,1,string("D"));


    GH = ProductGraph(GM,HM);

    HG = ProductGraph(HM,GM);

    PrintGraph(GH);

    PrintGraph(HG);

  /*  sequiv[string("xu")]=string("vz");
    sequiv[string("yw")]=string("ux");
    sequiv[string("xv")]=string("uy");
    sequiv[string("zu")]=string("wx");
    sequiv[string("zv")]=string("wy");
*/
  //  Textile T = FromSSE(GM,HM,sequiv);

//    PrintFullTextileInfo(T);

    vector<graph_traits<Graph>::edge_descriptor> permTest(5);

    for(tie(ei,ei_end)=edges(HG),i=0;ei!=ei_end;ei++,i++)
    {
        permTest[i]=*ei;
    } 

    auto start = permTest.begin(),finish = permTest.end();
do {
    for(tie(ei,ei_end)=edges(GH),tie(fi,fi_end)=edges(HG),i=0;ei!=ei_end;ei++,i++,fi++)
    {
        cout<< GH_ename(*ei) << " and " << HG_ename(*fi) << " and " << HG_ename(permTest[i]) << endl;
        sequiv[GH_ename(*ei)]=HG_ename(permTest[i]);
        cout << GH_ename(*ei) << " -> " << sequiv[GH_ename(*ei)] << endl;
    } 

    cout << endl;

    Textile T = FromSSE(GM,HM,sequiv);

    PrintFullTextileInfo(T);

    if(is1to1(T))
    {
        cout << "T is 1-1" << endl;
        Textile Td = CreateDual(T);
        cout << "Is Td 1-1?" << is1to1(Td) << endl;
    }
    else
    {
        cout << "T is not 1-1" << endl;
    }

} while(std::next_permutation(start,finish));


	return 0;
}
