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
    
    GammaGraph A(3);
    
    
    ofstream os("gms.txt");

    add_edge(0, 0, PQ_Homoms(string("a"),Q_Homom(string("b"),string("01"))),A);
    add_edge(0, 1, PQ_Homoms(string("a"),Q_Homom(string("c"),string("02"))),A);
    add_edge(1, 0, PQ_Homoms(string("b"),Q_Homom(string("d"),string("03"))),A);
    add_edge(1, 1, PQ_Homoms(string("b"),Q_Homom(string("e"),string("04"))),A);
    add_edge(0, 2, PQ_Homoms(string("c"),Q_Homom(string("a"),string("05"))),A);
    add_edge(2, 0, PQ_Homoms(string("d"),Q_Homom(string("b"),string("06"))),A);
    add_edge(2, 1, PQ_Homoms(string("d"),Q_Homom(string("c"),string("07"))),A);
    add_edge(2, 2, PQ_Homoms(string("e"),Q_Homom(string("a"),string("08"))),A);
    
    
    
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
    put(avname,0,string("u"));
    put(ap_vhom,1,0);
    put(aq_vhom,1,1);
    put(avname,1,string("v"));
    put(ap_vhom,2,1);
    put(aq_vhom,2,0);
    put(avname,2,string("w"));
    
    Graph B(2);
    
    add_edge(0, 0, string("a"),B);
    add_edge(0, 0, string("b"),B);
    add_edge(0, 1, string("c"),B);
    add_edge(1, 0, string("d"),B);
    add_edge(1, 1, string("e"),B);
    
    property_map<Graph,vertex_name_t>::type
    bvname = get(vertex_name,B);
    
    put(bvname,0,string("Y"));
    put(bvname,1,string("Z"));
    
    
    Textile T(A,B);
    
    PrintFullTextileInfo(T);
    

    Textile Tone = ArrayTrim(AutoHomom(CreateNMTextile(T,-1,2)));
    
    Textile Tdthree = Trim(HigherNBlock(CreateDual(Tone),3));


//    PrintFullTextileInfo(Tone);


//    Textile Tdthree = Trim(HigherNBlock(CreateDual(Tone),3));

//    PrintFullTextileInfo(Tdthree);

    Textile Ttwo = CreateDual(Tdthree);

    cout << is1to1(Ttwo) << endl;



  
    vector<vector<graph_traits<GammaGraph>::vertex_descriptor> > E(4);
    
    graph_traits<GammaGraph>::vertex_descriptor EAi[] = {0,5};
    graph_traits<GammaGraph>::vertex_descriptor EBi[] = {1,3,6};
    graph_traits<GammaGraph>::vertex_descriptor ECi[] = {2,7,4};
    graph_traits<GammaGraph>::vertex_descriptor EDi[] = {8,9,10};
   
    vector<graph_traits<GammaGraph>::vertex_descriptor> EA (EAi,EAi+sizeof(EAi)/sizeof(graph_traits<GammaGraph>::vertex_descriptor));
    vector<graph_traits<GammaGraph>::vertex_descriptor> EB (EBi,EBi+sizeof(EBi)/sizeof(graph_traits<GammaGraph>::vertex_descriptor));
    vector<graph_traits<GammaGraph>::vertex_descriptor> EC (ECi,ECi+sizeof(ECi)/sizeof(graph_traits<GammaGraph>::vertex_descriptor));
    vector<graph_traits<GammaGraph>::vertex_descriptor> ED (EDi,EDi+sizeof(EDi)/sizeof(graph_traits<GammaGraph>::vertex_descriptor));
   
    E[0]=EA;
    E[1]=EB;
    E[2]=EC;
    E[3]=ED;

    Textile Tq = Quotient(Tone,E);
    
    PrintFullTextileInfo(Tq);
    
    Textile Taq = AutoHomom(Tone);
    
    PrintFullTextileInfo(Taq);
    
    for(k=-1; k>-3;k--)
    {
        for(l=3; l<7; l++)
        {
            found = false;
          Textile Tkl = AutoHomom(CreateNMTextile(T,k,l));
              /*  if(IsLR(Tkl))
                {
                    grid << "We are LR at SE " << j << " where k = " << k << " and l = "  << l << endl;
                }
                else{
                    grid << "We are NOT LR at SE " << j << " where k = " << k << " and l = "  << l << endl;
                } */
               // cout << "CHECKING CONJUGACY FOR k = " << k << " and l = " << l << endl;
                Textile Tconj = LookForConjugacy(Tkl,4,&found);
                cout << "Found = " <<  found << endl;
                if(found){
                    os << "PRINTING IJ CONJUGACY FOR " << k << l << endl;
                    PrintFullTextileInfo(Tkl,os);
                    PrintFullTextileInfo(Tconj,os);
                }
        }
    }
    
    return 0;
}