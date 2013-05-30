//
//  gmsquare.cpp
//  
//
//  Created by Brendan Berg on 10/4/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
//
//  gmcube.cpp
//  
//
//  Created by Brendan Berg on 9/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <tuple>
#include <unordered_map>
#include "textileHelper.h"
#include "output.h"

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
    bool found;
    int i,j,k,l;
    
    GammaGraph A(6);
    
    
    ofstream os("Odense.txt");

    add_edge(0, 0, PQ_Homoms(string("a"),Q_Homom(string("a"),string("01"))),A); // aa -(a/a)-> aa 
    add_edge(0, 1, PQ_Homoms(string("a"),Q_Homom(string("b"),string("02"))),A); // aa -(a/b)-> ab
    add_edge(1, 2, PQ_Homoms(string("a"),Q_Homom(string("c"),string("03"))),A); // ab -(a/c)-> bc
    add_edge(1, 3, PQ_Homoms(string("a"),Q_Homom(string("d"),string("04"))),A); // ab -(a/d)-> bd

    add_edge(2, 4, PQ_Homoms(string("b"),Q_Homom(string("a"),string("05"))),A); // bc -(b/a)-> ca
    add_edge(3, 5, PQ_Homoms(string("b"),Q_Homom(string("e"),string("06"))),A); // bd -(b/e)-> cb
    add_edge(2, 6, PQ_Homoms(string("b"),Q_Homom(string("b"),string("07"))),A); // bc -(b/b)-> de
    add_edge(3, 7, PQ_Homoms(string("b"),Q_Homom(string("f"),string("08"))),A); // bd -(b/f)-> df

    add_edge(4, 0, PQ_Homoms(string("c"),Q_Homom(string("a"),string("09"))),A); // ca -(c/a)-> aa
    add_edge(7, 1, PQ_Homoms(string("c"),Q_Homom(string("e"),string("10"))),A); // df -(c/e)-> ab ****
    add_edge(5, 2, PQ_Homoms(string("c"),Q_Homom(string("c"),string("11"))),A); // cb -(c/c)-> bc
    add_edge(5, 3, PQ_Homoms(string("c"),Q_Homom(string("d"),string("12"))),A); // cb -(c/d)-> bd
    
    add_edge(6, 8, PQ_Homoms(string("d"),Q_Homom(string("c"),string("13"))),A); // de -(d/c)-> ec 
    add_edge(6, 9, PQ_Homoms(string("d"),Q_Homom(string("d"),string("14"))),A); // de -(d/d)-> ed  
    add_edge(4, 10, PQ_Homoms(string("d"),Q_Homom(string("b"),string("15"))),A); // ca -(d/b)-> fe ****
    add_edge(7, 11, PQ_Homoms(string("d"),Q_Homom(string("f"),string("16"))),A); // df -(d/f)-> ff

    add_edge(8, 4, PQ_Homoms(string("e"),Q_Homom(string("a"),string("17"))),A); // ec -(e/a)-> ca
    add_edge(9, 5, PQ_Homoms(string("e"),Q_Homom(string("e"),string("18"))),A); // ed -(e/e)-> cb
    add_edge(8, 6, PQ_Homoms(string("e"),Q_Homom(string("b"),string("19"))),A); // ec -(e/b)-> de
    add_edge(9, 7, PQ_Homoms(string("e"),Q_Homom(string("f"),string("20"))),A); // ed -(e/f)-> df

    add_edge(10, 8, PQ_Homoms(string("f"),Q_Homom(string("c"),string("21"))),A); // fe -(f/c)-> ec
    add_edge(10, 9, PQ_Homoms(string("f"),Q_Homom(string("d"),string("22"))),A); // fe -(f/d)-> ed
    add_edge(11,10, PQ_Homoms(string("f"),Q_Homom(string("e"),string("23"))),A); // ff -(f/e)-> fe
    add_edge(11,11, PQ_Homoms(string("f"),Q_Homom(string("f"),string("24"))),A); // ff -(f/f)-> ff
    
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    ap_vhom = get(vertex_p_vhomom,A);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    aq_vhom = get(vertex_q_vhomom,A);
    
    property_map<GammaGraph,vertex_name_t>::type
    avname = get(vertex_name,A);
    
    property_map<GammaGraph,edge_name_t>::type
    aname = get(edge_name,A);
    
    put(ap_vhom,0,0);
    put(aq_vhom,0,0);
    put(avname,0,string("aa"));
    put(ap_vhom,1,0);
    put(aq_vhom,1,1);
    put(avname,1,string("ab"));
    put(ap_vhom,2,0);
    put(aq_vhom,2,0);
    put(avname,2,string("bc"));
    put(ap_vhom,3,0);
    put(aq_vhom,3,2);
    put(avname,3,string("bd"));
    put(ap_vhom,4,1);
    put(aq_vhom,4,0);
    put(avname,4,string("ca"));
    put(ap_vhom,5,1);
    put(aq_vhom,5,1);
    put(avname,5,string("cb"));
    put(ap_vhom,6,1);
    put(aq_vhom,6,1);
    put(avname,6,string("de"));
    put(ap_vhom,7,1);
    put(aq_vhom,7,2);
    put(avname,7,string("df"));
    put(ap_vhom,8,2);
    put(aq_vhom,8,0);
    put(avname,8,string("ec"));
    put(ap_vhom,9,2);
    put(aq_vhom,9,2);
    put(avname,9,string("ed"));
    put(ap_vhom,10,2);
    put(aq_vhom,10,1);
    put(avname,10,string("fe"));
    put(ap_vhom,11,2);
    put(aq_vhom,11,2);
    put(avname,11,string("ff"));


    
    Graph B(3);
    
    add_edge(0, 0, string("a"),B);
    add_edge(0, 1, string("b"),B);
    add_edge(1, 0, string("c"),B);
    add_edge(1, 2, string("d"),B);
    add_edge(2, 1, string("e"),B);
    add_edge(2, 2, string("f"),B);
    
    property_map<Graph,vertex_name_t>::type
    bvname = get(vertex_name,B);
    
    put(bvname,0,string("Y"));
    put(bvname,1,string("Z"));
    
    
    Textile T(A,B);
    
    cout << "Our original T is" << endl;

    PrintFullTextileInfo(T);

    cout << "isLR(T) is " << IsLR(T) << endl;
    
    Textile TNTwoOne = CreateNMTextile(T,-2,1);

  //  Textile TNOneOne = CreateNMTextile(T,-1,1);

 //   PrintFullTextileInfo(TNTwoOne);


  //  LookForConjugacy(TNOneOne,5,&found);

    LookForConjugacy(TNTwoOne,5,&found);

 //   LookForConjugacy(T,5,&found);

    return 0;
}