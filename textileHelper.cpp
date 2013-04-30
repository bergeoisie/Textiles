//
//  textileHelper.cpp
//  
//
//  Created by Brendan Berg on 9/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "textileHelper.h"

//#include "tnt.h"
//#include "jama_eig.h"

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

using namespace std;
//using namespace TNT;
//using namespace JAMA;
using namespace boost;

// This allows us to create new properties later on and will
// be used to store the homomorphism data on the GammaGraph.
/*namespace boost { 
    enum edge_p_homom_t { edge_p_homom };
    enum edge_q_homom_t { edge_q_homom };
    enum vertex_p_vhomom_t { vertex_p_vhomom };
    enum vertex_q_vhomom_t { vertex_q_vhomom };
    
    BOOST_INSTALL_PROPERTY(edge, p_homom);
    BOOST_INSTALL_PROPERTY(edge, q_homom);
    BOOST_INSTALL_PROPERTY(vertex, p_vhomom);
    BOOST_INSTALL_PROPERTY(vertex, q_vhomom);
}*/

struct eqstr
{
    bool operator()(const char* s1, const char* s2) const
    {
        return strcmp(s1, s2) == 0;
    }
};

#define DEBUG 1




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

//enum colors { White, Gray, Black };

random::mt19937 gen(time(0));

class pqcomparison
{
public:
    bool operator() (const PQOEIElement& lhs, const PQOEIElement& rhs) const
    {
        return std::get<0>(rhs) > std::get<0>(lhs);
    }
};

/*
 * This function takes in a textile system T=(p,q:Gamma->G)
 * and outputs its dual system T*=(pT,qT:GammaT->GT)
 */
Textile CreateDual(Textile T)
{
    // int sizeOfVG = num_vertices(T.second);
    int sizeOfAG = num_edges(T.second);
    // int sizeOfVGamma = num_vertices(T.first);
    //int sizeOfAGamma = num_edges(T.first);
    int i;
    string p,q;
    
    VertexMap vmap;
    EdgeMap edgemap;
    
    // Initialize GT and GammaT
    Graph GT=T.second;
    GammaGraph GammaT(sizeOfAG);
    
    //  vector<bool> visited(sizeOfVG,false);
    
//	cout << "Creating Dual" << endl;
    
    // Set up access for the vertex names, indices, etc
    // for G
    property_map<Graph,vertex_name_t>::type
    g_vname = get(vertex_name,T.second);
    
    property_map<Graph,vertex_index_t>::type
    g_vindex = get(vertex_index,T.second);
    
    property_map<Graph,edge_name_t>::type
    g_name = get(edge_name,T.second);
    
    // The next few lines allow us to access the p_homoms
    // and q_homoms and enames of Gamma.
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,T.first);
    
    property_map<GammaGraph,edge_name_t>::type
    gamma_name = get(edge_name,T.first);
    
    property_map<GammaGraph,vertex_name_t>::type
    gamma_vname = get(vertex_name,T.first);
    
    property_map<GammaGraph,vertex_index_t>::type
    gamma_vindex = get(vertex_index,T.first);
    
    // And now we get access to the p_vhomoms and q_vhomoms
    // and vertex names of Gamma
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_vhom = get(vertex_p_vhomom,T.first);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_vhom = get(vertex_q_vhomom,T.first);
    
    // We now get access to the GT and GammaT vertex names
    property_map<Graph,vertex_name_t>::type
    gt_vname = get(vertex_name,GT);
    
    property_map<Graph,edge_name_t>::type
    gt_name = get(edge_name,GT);
    
    property_map<GammaGraph,edge_name_t>::type
    gammat_name = get(edge_name,GammaT);
    
    property_map<GammaGraph,vertex_name_t>::type
    gammat_vname = get(vertex_name,GammaT);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    gammat_pvhom = get(vertex_p_vhomom,GammaT);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    gammat_qvhom = get(vertex_q_vhomom,GammaT);
    
    // Let us clear the edges of GT just leaving the named vertices
    // which are there because we created GT by copying G
    graph_traits<Graph>::vertex_iterator wi, wi_end;
    for(tie(wi,wi_end)=vertices(GT); wi != wi_end; ++wi)
    {
        clear_vertex(*wi,GT);
    }
    
    
    // We add the edges to GT
//	cout << "Adding Edges to GT" << endl;
    graph_traits<GammaGraph>::vertex_iterator vi, vi_end;
    for(tie(vi,vi_end)=vertices(T.first); vi != vi_end; ++vi)
    {
        add_edge(p_vhom(*vi),q_vhom(*vi),Name(gamma_vname(*vi)),GT);
    }
    
    
    
    // Add the names to the vertices of GammaT
//	cout << "Adding names to verts of GammaT" << endl;
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for(tie(ei,ei_end)=edges(T.second), i=0; ei != ei_end; ++ei, i++)
    {
        put(gammat_vname,i,g_name(*ei));
        //   cout << "GammaT vname " << i << " is " << g_name(*ei) << endl;
        put(gammat_pvhom,i,source(*ei,T.second));
        put(gammat_qvhom,i,target(*ei,T.second));
        //    edgemap.insert(EdgeMap::value_type(g_name(*ei),
        //				 Edge(source(*ei,T.second),
        //				      target(*ei,T.second))));
        vmap.insert(VertexMap::value_type(g_name(*ei),i));
        
    }
    // Add the edges to GammaT
//	cout << "Adding Edges to GammaT" << endl;
    graph_traits<GammaGraph>::edge_iterator fi, fi_end;
    for(tie(fi,fi_end)=edges(T.first),i=0; fi != fi_end; ++fi, i++)
    {
	//	cout << i<< endl;
        p=p_homom(*fi);
        q=q_homom(*fi);
        
        //    EdgeMap::iterator emi=edgemap.find(p),fmi=edgemap.find(q);
        VertexMap::iterator vmi=vmap.find(p),wmi=vmap.find(q);
        
        
        add_edge(vmi->second,
                 wmi->second,
                 PQ_Homoms(gamma_vname(source(*fi,T.first)),
                           Q_Homom(gamma_vname(target(*fi,T.first)),
                                   gamma_name(*fi))),
                 GammaT);
        
    }
    
    cout << "GammaT vnames really are: ";
    graph_traits<GammaGraph>::vertex_iterator ci,ci_end;
    for(tie(ci,ci_end)=vertices(GammaT); ci!=ci_end; ++ci)
    {
        cout << gammat_vname(*ci) << " ";
    }
    cout << endl;
    
    return Textile(GammaT,GT);
}

Textile CreateInverse(Textile T)
{
    int i;
    GammaGraph A(num_vertices(T.first));
    Graph B=T.second;
    
    graph_traits<GammaGraph>::edge_iterator ei,ei_end;
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end,gi,gi_end;
    
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,T.first);
    
    property_map<GammaGraph,edge_name_t>::type
    ename = get(edge_name,T.first);
    
    property_map<GammaGraph,vertex_name_t>::type
    gammavname = get(vertex_name,T.first);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_avhom = get(vertex_p_vhomom,A);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_avhom = get(vertex_q_vhomom,A);  
    
    property_map<GammaGraph,vertex_name_t>::type
    avname = get(vertex_name,A);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_gammavhom = get(vertex_p_vhomom,T.first);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_gammavhom = get(vertex_q_vhomom,T.first);
    
    for(tie(ei,ei_end)=edges(T.first);ei!=ei_end;ei++)
    {
        add_edge(source(*ei,T.first),target(*ei,T.first),PQ_Homoms(q_homom(*ei),Q_Homom(p_homom(*ei),ename(*ei))),A);
    }
    
    for(tie(gi,gi_end)=vertices(T.first),tie(vi,vi_end)=vertices(A);gi!=gi_end;gi++,vi++)
    {
        put(p_avhom,*vi,q_gammavhom(*gi));
        put(q_avhom,*vi,p_gammavhom(*gi));
        put(avname,*vi,gammavname(*gi));
    }
    return Textile(A,B);
    
}

Textile ProductTextile(Textile T, Textile S)
{
    Graph GProd(num_vertices(T.second));
    GammaGraph GammaProd(num_vertices(T.first));
    string pword,qword,name;
    
    // P, Q Homom definitions  
    property_map<GammaGraph,edge_p_homom_t>::type
    p_Thomom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_Thomom = get(edge_q_homom,T.first);
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_Shomom = get(edge_p_homom,S.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_Shomom = get(edge_q_homom,S.first);  
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_Tvhom = get(vertex_p_vhomom,T.first);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_Tvhom = get(vertex_q_vhomom,T.first);  
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_Svhom = get(vertex_p_vhomom,S.first);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_Prodvhom = get(vertex_q_vhomom,GammaProd);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_Prodvhom = get(vertex_p_vhomom,GammaProd);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_Svhom = get(vertex_q_vhomom,S.first);
    
    // Edge and Vertex Names
    property_map<GammaGraph,edge_name_t>::type
    Tgamma_ename = get(edge_name,T.first);
    
    property_map<GammaGraph,vertex_name_t>::type
    Tgamma_vname = get(vertex_name,T.first);
    
    property_map<GammaGraph,edge_name_t>::type
    Sgamma_ename = get(edge_name,T.first);
    
    property_map<GammaGraph,vertex_name_t>::type
    Sgamma_vname = get(vertex_name,T.first);
    
    property_map<Graph,edge_name_t>::type
    Tg_ename = get(edge_name,T.second);
    
    property_map<Graph,vertex_name_t>::type
    Tg_vname = get(vertex_name,T.second);
    
    property_map<Graph,edge_name_t>::type
    Sg_ename = get(edge_name,T.second);
    
    property_map<Graph,vertex_name_t>::type
    Sg_vname = get(vertex_name,T.second);
    
    property_map<Graph,edge_name_t>::type
    Prodg_ename = get(edge_name,GProd);
    
    property_map<Graph,vertex_name_t>::type
    Prodg_vname = get(vertex_name,GProd);
    
    property_map<GammaGraph,edge_name_t>::type
    Prodgamma_ename = get(edge_name,GammaProd);
    
    property_map<GammaGraph,vertex_name_t>::type
    Prodgamma_vname = get(vertex_name,GammaProd);
    
    
    // Necessary vertex iterators
    graph_traits<GammaGraph>::vertex_iterator tgammavi,tgammavi_end,ProdGammavi,ProdGammavi_end;
    graph_traits<Graph>::vertex_iterator tgvi,tgvi_end,pgvi,pgvi_end,sgvi,sgvi_end;
    
    // Necessary out_edge iterators
    graph_traits<GammaGraph>::out_edge_iterator tgammaoei,tgammaoei_end,sgammaoei,sgammaoei_end;
    graph_traits<Graph>::out_edge_iterator tgoei,tgoei_end,sgoei,sgoei_end;
    
    // We need to make sure our graphs are compatible
    for(tie(tgvi,tgvi_end)=vertices(T.second);tgvi!=tgvi_end;tgvi++)
    {
        cout << "We're in the for loop" << endl;
        graph_traits<GammaGraph>::vertex_descriptor currV=*tgvi;
        if(num_vertices(T.first)!= num_vertices(S.first) || num_vertices(T.second) != num_vertices(S.second) 
           || p_Tvhom(currV) != p_Svhom(currV) || q_Tvhom(currV) != q_Svhom(currV) 
           ||  Tgamma_vname(currV) != Sgamma_vname(currV))
        {
            cout << "WARNING! THIS PRODUCT WILL NOT WORK" << endl;
            return Textile(GammaProd,GProd);
        }
        else {
            cout << "We can make a product! Yay!" << endl;
        }
    }
    
    
    // We first create G_1G_2
    
    // Add names
    for(tie(tgvi,tgvi_end)=vertices(T.second),tie(pgvi,pgvi_end)=vertices(GProd);
        tgvi!=tgvi_end; tgvi++, pgvi++) {
        put(Prodg_vname,*pgvi,Tg_vname(*tgvi));
    }
    
    // We now need to create the edges of G_1G_2
    for(tie(tgvi,tgvi_end)=vertices(T.second); tgvi!=tgvi_end;tgvi++)
    {
        for(tie(tgoei,tgoei_end)=out_edges(*tgvi,T.second);tgoei!=tgoei_end;tgoei++)
        {
            for(tie(sgoei,sgoei_end)=out_edges(target(*tgoei,T.second),S.second);sgoei!=sgoei_end;sgoei++)
            {
                add_edge(*tgvi,target(*sgoei,S.second),string( Tg_ename(*tgoei) + Sg_ename(*sgoei)),GProd);
            }
        }
    }
    // Let's set up the vertices of Gamma_1Gamma_2
    for(tie(tgammavi,tgammavi_end)=vertices(T.first);
        tgammavi!=tgammavi_end;tgammavi++)
    {
        put(Prodgamma_vname,*tgammavi,Tgamma_vname(*tgammavi));
        put(p_Prodvhom,*tgammavi,p_Tvhom(*tgammavi));
        put(q_Prodvhom,*tgammavi,q_Tvhom(*tgammavi));
    }
    
    // Now's the hard part. We have to set up the equivalence relation on Gamma_1Gamma_2
    for(tie(tgammavi,tgammavi_end)=vertices(T.first);tgammavi!=tgammavi_end;tgammavi++)
    {
        graph_traits<GammaGraph>::vertex_descriptor currV=*tgammavi;
        for(tie(tgammaoei,tgammaoei_end)=out_edges(currV,T.first);
            tgammaoei!=tgammaoei_end; tgammaoei++)
        {
            graph_traits<GammaGraph>::edge_descriptor currTEdge=*tgammaoei;
            for(tie(sgammaoei,sgammaoei_end)=out_edges(target(currTEdge,T.first),S.first);
                sgammaoei!=sgammaoei_end; sgammaoei++)
            {
                graph_traits<GammaGraph>::edge_descriptor currSEdge=*sgammaoei;
                pword=p_Thomom(currTEdge)+p_Shomom(currSEdge);
                qword=q_Thomom(currTEdge)+q_Shomom(currSEdge);
                name=Tgamma_ename(currTEdge)+Sgamma_ename(currSEdge);	     
                if(DEBUG)
                {
         //           cout << "We have created an edge from " << currV << " to " << target(currSEdge,S.first)
         //           << " named " << name << " where p and q are " << pword << " and " << qword
         //           << " respectively." << endl;
                }
                add_edge(currV,target(currSEdge,S.first),PQ_Homoms(pword,Q_Homom(qword,name)),GammaProd);
            }
        }
    }
    
    
    return Textile(GammaProd,GProd);
}

bool CheckNondegen(Graph G)
{
    unsigned int i;
    vector<bool> inits(num_vertices(G),false),tails(num_vertices(G),false);
    
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for(tie(ei,ei_end)=edges(G); ei != ei_end; ++ei) {
        inits[source(*ei,G)]=true;
        tails[target(*ei,G)]=true;
    }
    
    for(i=0;i<num_vertices(G);i++) {
        if(!inits[i] || !tails[i])
            cout << "WARNING, GRAPH IS DEGENERATE." << endl;
        return false;
    }
    
    return true;
}

bool checkHomoms(Textile T)
{
    int i,pi,qi,pt,qt;
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,T.first);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_vhom = get(vertex_p_vhomom,T.first);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_vhom = get(vertex_q_vhomom,T.first);
    
    property_map<Graph,edge_name_t>::type
    g_name = get(edge_name,T.second);
    
    
    // Create and populate our EdgeMap. This map associates edges in G with
    // their labels.
    EdgeMap edgemap;
    graph_traits<Graph>::edge_iterator fi, fi_end;
    for(tie(fi,fi_end)=edges(T.second), i=0; fi != fi_end; ++fi, i++)
    {
        const string tstr = g_name(*fi);
        const Edge tedge = Edge(source(*fi,T.second), target(*fi,T.second));
        edgemap.insert(EdgeMap::value_type(tstr,tedge));
    }
    
    
    graph_traits<GammaGraph>::edge_iterator ei, ei_end;
    for(tie(ei,ei_end)=edges(T.first); ei != ei_end; ++ei) 
    {
        EdgeMap::iterator p=edgemap.find(p_homom(*ei)),q=edgemap.find(q_homom(*ei));
        
        if(p != edgemap.end()) {
            pi=(p->second).first;
            pt=(p->second).second;
        }
        else {
            cout << "ERROR, CANNOT FIND P HOMOM EDGE" << endl;
        }
        
        if(q != edgemap.end()) {
            qi=(q->second).first;
            qt=(q->second).second;
        }
        else {
            cout << "ERROR, CANNOT FIND Q HOMOM EDGE" << endl;
        }
        if(pi!=p_vhom(source(*ei,T.first)) ||
           qi!=q_vhom(source(*ei,T.first)) ||
           pt!=p_vhom(target(*ei,T.first)) ||
           qt!=q_vhom(target(*ei,T.first))) 
        {
            cout << "There's a problem with " << "(" << source(*ei, T.first) 
            << "," << target(*ei, T.first) << ") " << endl;
            
            if(pi!=p_vhom(source(*ei,T.first))) {
                cout << "i_Gp_A = " << pi << " while "
                << "p_Vi_Gamma = " << p_vhom(source(*ei,T.first)) << endl;
            }
            
            if(qi!=q_vhom(source(*ei,T.first))) {
                cout << "i_Gq_A = " << qi << " while "
                << "q_Vi_Gamma = " << q_vhom(source(*ei,T.first)) << endl;
            }
            
            if(pt!=p_vhom(target(*ei,T.first))) {
                cout << "t_Gp_A = " << pt << " while "
                << "p_Vt_Gamma = " << p_vhom(target(*ei,T.first)) << endl;
            }
            
            if(qt!=q_vhom(target(*ei,T.first))) {
                cout << "t_Gq_A = " << qt << " while "
                << "q_Vt_Gamma = " << q_vhom(target(*ei,T.first)) << endl;
            }
            
            return false;
        }
        
    }
    
    return true;
}

Textile Trim(Textile T)
{
    Textile Trimmed=T;
    bool isEndIT=true,isEndTI=true,isEndPQ=true,isEndQP=true;
    unsigned int i,j=0,delcount=0,N;
    //  vector<graph_traits<GammaGraph>::edge_descriptor> toDelete;
    vector<int> vertToDelete;
    bool done,stable;
    time_t start,end;
    double dif,average=0;
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,Trimmed.first);
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,Trimmed.first);
    property_map<GammaGraph,edge_name_t>::type
    ename = get(edge_name,Trimmed.first);
    property_map<GammaGraph,vertex_name_t>::type
    vname = get(vertex_name,Trimmed.first);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_vhom = get(vertex_p_vhomom,Trimmed.first);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_vhom = get(vertex_q_vhomom,Trimmed.first);
    
    property_map<Graph,edge_name_t>::type
    gename = get(edge_name,Trimmed.second);
    
    graph_traits<GammaGraph>::edge_iterator ei, ei_end, fi, fi_end;
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end;
    graph_traits<GammaGraph>::out_edge_iterator oei,oei_end;
    graph_traits<GammaGraph>::in_edge_iterator iei,iei_end;
    graph_traits<Graph>::edge_iterator gei,gei_end;
    //  graph_traits<GammaGraph>::edge_descriptor e;
    GVI gvi,gvi_end;
	bool needToCheckVhoms=false,found=false,topStart=true;;
    
    stack<VD> toCheck;  

	cout << "Our starting graph has " << num_edges(Trimmed.first) << " edges and " << num_vertices(Trimmed.first) << " vertices." << endl; 

    time(&start);
    do{
        j++;
       if(j % 100 == 0)
       {
           time(&end);
           dif = difftime(end,start);
           time(&start);
           N = j/100;
           average = (((N-1)*average)+dif)/N;
           cout << "On round " << j << " deleted at least " << delcount << " edges so far. The last round took " 
                << dif << " seconds, and the average time is " << average << " seconds." <<  '\r';
           cout.flush();
       } 
        stable = true;
        done = false;
        
        // Check for sinks 
        for(tie(vi,vi_end)=vertices(Trimmed.first); vi!=vi_end && !done; vi++)
        {
            //	cout << "We are checking " << vname(*vi) << endl;
            if(out_degree(*vi,Trimmed.first)==0 || in_degree(*vi,Trimmed.first)==0)
            {
                //   cout << "The vertex " << vname(*vi) << "is a sink or a source" << endl;
                clear_vertex(*vi,Trimmed.first);
                remove_vertex(*vi,Trimmed.first);
                stable = false;
                done=true;
                break;
            }
            //	cout << "It's not." << endl;
        }
        if(!done)
        {
            //	cout << "Starting ends loop" << endl;
            for(tie(ei,ei_end)=edges(Trimmed.first); ei != ei_end && !done; ++ei) {
                //  cout << "Checking " << ename(*ei) << " for ends" << endl;
                for(tie(fi,fi_end)=edges(Trimmed.first); fi != fi_end; ++fi) {
                    
                    if(source(*ei,Trimmed.first)==target(*fi,Trimmed.first)) {
                        isEndIT=false;
                    }
                    
                    if(target(*ei,Trimmed.first)==source(*fi,Trimmed.first)) {
                        isEndTI=false;
                    }
                    
                    if(p_homom(*ei)==q_homom(*fi)) {
                        isEndPQ=false;
                    }
                    
                    if(q_homom(*ei)==p_homom(*fi)) {
                        isEndQP=false;
                    }
                }
                if(isEndIT || isEndTI || isEndPQ || isEndQP) {
                    VD s=source(*ei,Trimmed.first),t=target(*ei,Trimmed.first);
                    toCheck.push(s);
                    toCheck.push(t);
                  //            cout << "Removed Edge (" << s << ", "
                  //            << t << ")" << " named " << ename(*ei) << endl;
                    remove_edge(*ei,Trimmed.first);
                    delcount++;
                    done = true;
                    stable = false;
                    // the goal of this inner while loop is to check to see if we've just created a dead end for the surrounding vertices.
                    while(!toCheck.empty())
                    {
                        VD currv = toCheck.top();
                        toCheck.pop();
                        
                        if(out_degree(currv,Trimmed.first)==0 && in_degree(currv,Trimmed.first)!=0)
                        {
							for(tie(iei,iei_end)=in_edges(currv,Trimmed.first);iei!=iei_end;iei++)
							{
								VD sour = source(*iei,Trimmed.first);
								toCheck.push(sour);
							}
							delcount += in_degree(currv,Trimmed.first);
                            //                  cout << "Clearing " << currv << endl;
                            clear_vertex(currv,Trimmed.first);
                        }
						else if(out_degree(currv,Trimmed.first)!=0 && in_degree(currv,Trimmed.first)==0)
						{
								for(tie(oei,oei_end)=out_edges(currv,Trimmed.first);oei!=oei_end;oei++)
								{
									VD tar = target(*oei,Trimmed.first);
									toCheck.push(tar);
								}
								delcount += out_degree(currv,Trimmed.first);
	                            //                  cout << "Clearing " << currv << endl;
	                            clear_vertex(currv,Trimmed.first);
						}
                        
                    }
                }
                isEndIT=true;
                isEndTI=true;
                isEndPQ=true;
                isEndQP=true;
            }
            //	cout << "Done with the ends, checking for empty vertices" << endl;
            if(!done)
            {
                for(tie(vi,vi_end)=vertices(Trimmed.first); vi!=vi_end && !done; vi++)
                {
                    if(out_degree(*vi,Trimmed.first)==0 && in_degree(*vi,Trimmed.first)==0)
                    {
                        remove_vertex(*vi,Trimmed.first);
                        done = true;
                        stable = false;
                    }
                }
            }
        }
        
    } while(!stable); // do while statement (don't panic)
    
    // Now that we have trimmed the GammaGraph, we need to check the downstairs Graph
    // First, we check for dead edges
    stable = false;
    while(!stable)
    {
        for(tie(gei,gei_end)=edges(Trimmed.second);gei!=gei_end; gei++)
        {
            found = false;
            for(tie(ei,ei_end)=edges(Trimmed.first); ei!=ei_end; ei++)
            {
                if(p_homom(*ei) == gename(*gei) || q_homom(*ei) == gename(*gei))
                {
                    found = true;
                    break;
                } // if
            } // for ei
            if(!found)
            {
                remove_edge(*gei,Trimmed.second);
                break;
            } // if found 
        } // for gei
        if(gei==gei_end)
        {
            stable = true;
        }
    }
    
    // Now we check for dead vertices
    stack<GVD> toDel;
    
    for(tie(gvi,gvi_end)=vertices(Trimmed.second); gvi!=gvi_end; gvi++)
    {
        bool found=false;
        
        for(tie(vi,vi_end)=vertices(Trimmed.first); vi!=vi_end; vi++)
        {
            if(p_vhom(*vi) == *gvi || q_vhom(*vi) == *gvi)
            {
                found = true;
                break;
            }
        }
        if(!found)
        {
            toDel.push(*gvi);
        }
        
    }
    
    if(!toDel.empty()) {
        needToCheckVhoms = true;
    }
    
    while(!toDel.empty())
    {
        GVD u = toDel.top();
        toDel.pop();
        
        clear_vertex(u,Trimmed.second);
        remove_vertex(u,Trimmed.second);
        
     //   cout << "We have deleted " << u << " from G" << endl;
    }
    
    
    if(needToCheckVhoms)
    {
        for(tie(vi,vi_end)=vertices(Trimmed.first);vi!=vi_end;vi++)
        {
            tie(oei,oei_end)=out_edges(*vi,Trimmed.first);
            found = false;
            tie(gei,gei_end)=edges(Trimmed.second);
            while(!found)
            {
                if(gei==gei_end)
                {
                    cout << "WE DIDN'T FIND A MATCH FOR A HOMOM AND AN G EDGE. THIS IS A BIG PROBLEM" << endl;
                    return T;
                }
                
                if(p_homom(*oei)==gename(*gei))
                {
                    found = true;
                    put(p_vhom,*vi,source(*gei,Trimmed.second));
                }
                else {
                    gei++;
                }
            }
            
            found = false;
            
            tie(gei,gei_end)=edges(Trimmed.second);
            while(!found)
            {
                if(gei==gei_end)
                {
                    cout << "WE DIDN'T FIND A MATCH FOR A HOMOM AND AN G EDGE. THIS IS A BIG PROBLEM" << endl;
                    return T;
                }
                
                if(q_homom(*oei)==gename(*gei))
                {
                    found = true;
                    put(q_vhom,*vi,source(*gei,Trimmed.second));
                }
                else {
                    gei++;
                }
            }
            
        }
        
    }
    
    cout << "Original T.first had " << num_edges(T.first) << " edges and " << num_vertices(T.first)
    << " vertices, now has " << num_edges(Trimmed.first) << " edges and " << num_vertices(Trimmed.first) << endl;
    
    cout << "Original T.second had " << num_edges(T.second) << " edges and " << num_vertices(T.second)
    << " vertices, now has " << num_edges(Trimmed.second) << " edges and " << num_vertices(Trimmed.second) << endl;
    
    return Trimmed;
}

Textile HigherBlock(Textile T)
{
    int sizeOfAG = num_edges(T.second);
    int sizeOfAGamma = num_edges(T.first);
    int i,s,t;
    
    Graph GTwo(sizeOfAG);
    GammaGraph GammaTwo(sizeOfAGamma);
    
    VertexMap vmap,xmap;
    
    property_map<Graph,edge_name_t>::type
    g_name = get(edge_name,T.second);
    
    property_map<GammaGraph,edge_name_t>::type
    gamma_name = get(edge_name,T.first);
    
    property_map<Graph,vertex_name_t>::type
    gtwo_vname = get(vertex_name,GTwo);
    
    property_map<GammaGraph,vertex_name_t>::type
    gammatwo_vname = get(vertex_name,GammaTwo);
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,T.first);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    tp_vhom = get(vertex_p_vhomom,GammaTwo);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    tq_vhom = get(vertex_q_vhomom,GammaTwo);
    
    graph_traits<Graph>::vertex_iterator vi,vi_end;
    graph_traits<Graph>::edge_iterator ei,ei_end;
    graph_traits<Graph>::out_edge_iterator fi,fi_end,gi,gi_end;
    VertexMap::iterator vmi,wmi,xmi,ymi;
    
    graph_traits<GammaGraph>::vertex_iterator gvi,gvi_end;
    graph_traits<GammaGraph>::edge_iterator gei,gei_end;
    graph_traits<GammaGraph>::out_edge_iterator gfi,gfi_end,ggi,ggi_end;
    
    // Populate our VertexMap
    for(tie(ei,ei_end)=edges(T.second), i=0; ei!=ei_end; ++ei,i++)
    {
        vmap.insert(VertexMap::value_type(g_name(*ei),i));
        put(gtwo_vname,i,g_name(*ei));
    }
    
    // Fill the edges of GTwo
    for(tie(vi,vi_end)=vertices(T.second); vi!=vi_end; ++vi)
    {
        for(tie(fi,fi_end)=out_edges(*vi,T.second); fi!=fi_end;++fi)
        {
            for(tie(gi,gi_end)=out_edges(target(*fi,T.second),T.second);
                gi!=gi_end;
                ++gi)
            {
                vmi=vmap.find(g_name(*fi));
                wmi=vmap.find(g_name(*gi));
                
                add_edge(vmi->second,
                         wmi->second,
                         Name(g_name(*fi)+g_name(*gi)),
                         GTwo);
            }
            
        }
        
    }
    
    
    // Populate our second VertexMap and the Vertex Homoms
    for(tie(gei,gei_end)=edges(T.first), i=0; gei!=gei_end; ++gei,i++)
    {
        xmap.insert(VertexMap::value_type(gamma_name(*gei),i));
        vmi=vmap.find(p_homom(*gei));
        wmi=vmap.find(q_homom(*gei));
        
        put(gammatwo_vname,i,gamma_name(*gei));
        put(tp_vhom,i,vmi->second);
        put(tq_vhom,i,wmi->second);
    }
    
    
    // Fill the edges of GTwo
    for(tie(gvi,gvi_end)=vertices(T.first); gvi!=gvi_end; ++gvi)
    {
        for(tie(gfi,gfi_end)=out_edges(*gvi,T.first); gfi!=gfi_end;++gfi)
        {
            for(tie(ggi,ggi_end)=out_edges(target(*gfi,T.first),T.first);
                ggi!=ggi_end;
                ++ggi)
            {
                xmi=xmap.find(gamma_name(*gfi));
                ymi=xmap.find(gamma_name(*ggi));
                
                add_edge(xmi->second,
                         ymi->second,
                         PQ_Homoms(p_homom(*gfi)+p_homom(*ggi),
                                   Q_Homom(q_homom(*gfi)+q_homom(*ggi),
                                           Name(gamma_name(*gfi)+gamma_name(*ggi)))),
                         GammaTwo);
                //     cout << "Created edge: " << gamma_name(*gfi)+gamma_name(*ggi) << endl;
            }
            
        }
        
    }
    
    
    
    return Textile(GammaTwo,GTwo);
}

// This function takes in a textile system and determines
// whether or not p,q generate the same language
bool IsomLanguages(Textile T)
{
    // The way this algorithm works is essentially depth first search. We will iterate through
    // the vertices and at each one create every p and q string of length k starting there. We then store
    // all of these string in memory and compare them once finished. We delete strings that have a match
    // in an effort to save memory.
    
    
    string currentpWord="",currentqWord;
    
    // These are the lists that store the p and q strings.
    list<string> plist,qlist;
    
    int n=num_edges(T.first);
    
    // The OEI pairs allow us to keep both the current location of the iterator and its end.
    typedef std::tuple<graph_traits<GammaGraph>::out_edge_iterator,
    graph_traits<GammaGraph>::out_edge_iterator,
    string,string> OEITuple;
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,T.first);
    
    property_map<Graph,edge_name_t>::type
    g_name = get(edge_name,T.second);
    
    
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end;
    
    graph_traits<GammaGraph>::out_edge_iterator oi,oi_end,currenti,currenti_end;
    
    stack<OEITuple> OEIStack;
    
    bool done;
    
    for(tie(vi,vi_end)=vertices(T.first); vi!=vi_end; vi++)
    {
        tie(oi,oi_end)=out_edges(*vi,T.first);
        
        OEIStack.push(OEITuple(oi,oi_end,"",""));
        
        while(!OEIStack.empty())
        {
            OEITuple current = OEIStack.top();
            
            currenti = get<0>(current);
            currenti_end = get<1>(current);
            currentpWord = get<2>(current);
            currentqWord = get<3>(current);
            
            if(currentpWord.length() < n-1) // We want to go deeper
            {
                graph_traits<GammaGraph>::out_edge_iterator ni,ni_end;
                //boost::graph_traits<GammaGraph>::vertex_descriptor v;
                tie(ni,ni_end)=out_edges(target(*currenti,T.first),T.first);
                OEIStack.push(OEITuple(ni,ni_end,currentpWord+p_homom(*currenti),currentqWord+q_homom(*currenti)));
                
                cout << "Our new p and q words are " << currentpWord+p_homom(*currenti) << " and " << currentqWord+q_homom(*currenti) << endl;
                
            }
            else // We're the deepest we want to go. 
            {
                OEIStack.pop();
                plist.push_front(currentpWord+p_homom(*currenti));
                qlist.push_front(currentqWord+q_homom(*currenti));
                currenti++;
                if(currenti!=currenti_end) {
                    OEIStack.push(OEITuple(currenti,currenti_end,currentpWord,currentqWord));
                }
                else { // Our top iterator has finished. We need to go down a level (or two..)
                    done = false;
                    while(!done && !OEIStack.empty()) {
                        current = OEIStack.top();
                        OEIStack.pop();
                        currenti = get<0>(current);
                        currenti_end = get<1>(current);
                        currentpWord = get<2>(current);
                        currentqWord = get<3>(current);
                        currenti++;
                        if(currenti!=currenti_end)
                        {
                            done=true;
                            OEIStack.push(OEITuple(currenti,currenti_end,currentpWord,currentqWord));
                        }
                        
                    }
                }
            }
            
        }
        while(plist.empty())
        {
            string first = plist.front();
            plist.pop_front();
            // WE NEED TO ITERATE THROUGH FIX THIS !!!!!
        }
    }
}


bool IsLR(Textile T)
{
    int matches=0;
    
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end;
    graph_traits<GammaGraph>::in_edge_iterator gammaiei, gammaiei_end;
    graph_traits<Graph>::in_edge_iterator giei, giei_end;
    graph_traits<GammaGraph>::out_edge_iterator gammaoei, gammaoei_end;
    graph_traits<Graph>::out_edge_iterator goei, goei_end;
    
    property_map<Graph,edge_name_t>::type
    g_name = get(edge_name,T.second);
    
    property_map<GammaGraph,vertex_name_t>::type
    gamma_vname = get(vertex_name,T.first);
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,T.first);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_vhom = get(vertex_p_vhomom,T.first);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_vhom = get(vertex_q_vhomom,T.first);
    
    
    for(tie(vi,vi_end)=vertices(T.first); vi!=vi_end; vi++)
    {
        cout << "Looking at vertex " << gamma_vname(*vi) << endl;
        
        for(tie(giei,giei_end)=in_edges(p_vhom(*vi),T.second);
            giei!=giei_end; giei++)
        {
            cout << "Looking at G edge " << g_name(*giei) << endl;
            for(tie(gammaiei,gammaiei_end)=in_edges(*vi,T.first);
                gammaiei!=gammaiei_end; gammaiei++) 
            {
                if(p_homom(*gammaiei)==g_name(*giei))
                {
                    cout << "Found a match!" << endl;
                    matches++;
                }
            }
            if(matches!=1) {
                cout << "This G edge has " << matches << " matches" << endl;
                return false;
            }
            matches=0;
        }
        
        cout << "Done checking in edges. Moving on to out edges" << endl;
        
        matches=0;
        for(tie(goei,goei_end)=out_edges(q_vhom(*vi),T.second);
            goei!=goei_end; goei++)
        {
            cout << "Looking at G edge " << g_name(*goei) << endl;
            for(tie(gammaoei,gammaoei_end)=out_edges(*vi,T.first);
                gammaoei!=gammaoei_end; gammaoei++)
            {
                if(q_homom(*gammaoei)==g_name(*goei))
                {
                    cout << "Found a match!" << endl;
                    matches++;
                }
            }
            if(matches!=1) {
                cout << "The G edge " << g_name(*goei) << " has " << matches << " matches" << endl;
                return false;
            }
            matches=0;
        }
    }
    
    cout << "We're LR!" << endl;
    return true;
}

bool IspRightResolving(Textile T)
{
    int matches=0;
    
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end;
    graph_traits<GammaGraph>::in_edge_iterator gammaiei, gammaiei_end;
    graph_traits<Graph>::in_edge_iterator giei, giei_end;
    graph_traits<GammaGraph>::out_edge_iterator gammaoei, gammaoei_end;
    graph_traits<Graph>::out_edge_iterator goei, goei_end;
    
    property_map<Graph,edge_name_t>::type
    g_name = get(edge_name,T.second);
    
    property_map<GammaGraph,vertex_name_t>::type
    gamma_vname = get(vertex_name,T.first);
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_vhom = get(vertex_p_vhomom,T.first);
    
    
    for(tie(vi,vi_end)=vertices(T.first); vi!=vi_end; vi++)
    {
        //      cout << "Looking at vertex " << gamma_vname(*vi) << endl;
        
        matches=0;
        for(tie(goei,goei_end)=out_edges(p_vhom(*vi),T.second);
            goei!=goei_end; goei++)
        {
            for(tie(gammaoei,gammaoei_end)=out_edges(*vi,T.first);
                gammaoei!=gammaoei_end; gammaoei++)
            {
                if(p_homom(*gammaoei)==g_name(*goei))
                {
                    //		  cout << "Found a match!" << endl;
                    matches++;
                }
            }
            if(matches!=1) {
                cout << "The G edge " << g_name(*goei) << " has " << matches << " matches" << endl;
                return false;
            }
            matches=0;
        }
    }
    
    cout << "p is right resolving!" << endl;
    return true;
}

bool IspLeftResolving(Textile T)
{
    int matches=0;
    
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end;
    graph_traits<GammaGraph>::in_edge_iterator gammaiei, gammaiei_end;
    graph_traits<Graph>::in_edge_iterator giei, giei_end;
    graph_traits<GammaGraph>::out_edge_iterator gammaoei, gammaoei_end;
    graph_traits<Graph>::out_edge_iterator goei, goei_end;
    
    property_map<Graph,edge_name_t>::type
    g_name = get(edge_name,T.second);
    
    property_map<GammaGraph,vertex_name_t>::type
    gamma_vname = get(vertex_name,T.first);
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_vhom = get(vertex_p_vhomom,T.first);
    
    
    for(tie(vi,vi_end)=vertices(T.first); vi!=vi_end; vi++)
    {
        //      cout << "Looking at vertex " << gamma_vname(*vi) << endl;
        
        for(tie(giei,giei_end)=in_edges(p_vhom(*vi),T.second);
            giei!=giei_end; giei++)
        {
            //  cout << "Looking at G edge " << g_name(*giei) << endl;
            for(tie(gammaiei,gammaiei_end)=in_edges(*vi,T.first);
                gammaiei!=gammaiei_end; gammaiei++) 
            {
                if(p_homom(*gammaiei)==g_name(*giei))
                {
                    //  cout << "Found a match!" << endl;
                    matches++;
                }
            }
            if(matches!=1) {
                cout << "The G edge " << g_name(*giei) << " has " << matches << " matches" << endl;
                return false;
            }
            matches=0;
        }
        
    }
    cout << "p is left resolving" << endl;
    return true;
    
}

bool IsqRightResolving(Textile T)
{
    int matches=0;
    
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end;
    graph_traits<GammaGraph>::in_edge_iterator gammaiei, gammaiei_end;
    graph_traits<Graph>::in_edge_iterator giei, giei_end;
    graph_traits<GammaGraph>::out_edge_iterator gammaoei, gammaoei_end;
    graph_traits<Graph>::out_edge_iterator goei, goei_end;
    
    property_map<Graph,edge_name_t>::type
    g_name = get(edge_name,T.second);
    
    property_map<GammaGraph,vertex_name_t>::type
    gamma_vname = get(vertex_name,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,T.first);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_vhom = get(vertex_q_vhomom,T.first);
    
    
    for(tie(vi,vi_end)=vertices(T.first); vi!=vi_end; vi++)
    {
        //      cout << "Looking at vertex " << gamma_vname(*vi) << endl;
        
        matches=0;
        for(tie(goei,goei_end)=out_edges(q_vhom(*vi),T.second);
            goei!=goei_end; goei++)
        {
            for(tie(gammaoei,gammaoei_end)=out_edges(*vi,T.first);
                gammaoei!=gammaoei_end; gammaoei++)
            {
                if(q_homom(*gammaoei)==g_name(*goei))
                {
                    //	  cout << "Found a match!" << endl;
                    matches++;
                }
            }
            if(matches!=1) {
                cout << "The G edge " << g_name(*goei) << " has " << matches << " matches" << endl;
                return false;
            }
            matches=0;
        }
    }
    
    cout << "q is right resolving!" << endl;
    return true;
}

bool IsqLeftResolving(Textile T)
{
    int matches=0;
    
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end;
    graph_traits<GammaGraph>::in_edge_iterator gammaiei, gammaiei_end;
    graph_traits<Graph>::in_edge_iterator giei, giei_end;
    graph_traits<GammaGraph>::out_edge_iterator gammaoei, gammaoei_end;
    graph_traits<Graph>::out_edge_iterator goei, goei_end;
    
    property_map<Graph,edge_name_t>::type
    g_name = get(edge_name,T.second);
    
    property_map<GammaGraph,vertex_name_t>::type
    gamma_vname = get(vertex_name,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,T.first);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_vhom = get(vertex_q_vhomom,T.first);
    
    
    for(tie(vi,vi_end)=vertices(T.first); vi!=vi_end; vi++)
    {
        //      cout << "Looking at vertex " << gamma_vname(*vi) << endl;
        
        for(tie(giei,giei_end)=in_edges(q_vhom(*vi),T.second);
            giei!=giei_end; giei++)
        {
            //  cout << "Looking at G edge " << g_name(*giei) << endl;
            for(tie(gammaiei,gammaiei_end)=in_edges(*vi,T.first);
                gammaiei!=gammaiei_end; gammaiei++) 
            {
                if(q_homom(*gammaiei)==g_name(*giei))
                {
                    //  cout << "Found a match!" << endl;
                    matches++;
                }
            }
            if(matches!=1) {
                cout << "The G edge " << g_name(*giei) << " has " << matches << " matches" << endl;
                return false;
            }
            matches=0;
        }
        
    }
    cout << "q is left resolving" << endl;
    return true;
}


int IsqRightDefinite(Textile T)
{
    return IspRightDefinite(CreateInverse(T));
}

int IspRightDefinite(Textile T)
{
    Graph matrixHelper;
    vector<vColl*> graphRef;
    int i=0,j=0,k=0,m=0,newM=0,maxM=0,typeNumber=num_vertices(T.second);
    graph_traits<Graph>::vertex_iterator gvi,gvi_end;
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end;
    graph_traits<GammaGraph>::out_edge_iterator oei,oei_end;
    graph_traits<Graph>::out_edge_iterator goei,goei_end;
    bool done, found, firstfound;
    stack<sPair> prodStack;
    
    // The production rules depend only on the type. They can be generated by looking at the out edges
    // of the downstairs vertex the type represents.
    vector<sColl> productionRules;
    
    vector<string> types;
    vector<string>::iterator it;

    cout << "About to create type matrix" << endl;
    
    // The type matrix is indexed so that the 0 index is actually representing the types of
    // size M (where M is the largest vertex set).
    boost::numeric::ublas::matrix<vector<vColl* >* > typeMatrix;
    
    cout << "Created type matrix" << endl;

    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_vhomom = get(vertex_p_vhomom,T.first);
    
    property_map<Graph,vertex_name_t>::type
    mh_vname = get(vertex_name,matrixHelper);
    
    property_map<Graph,vertex_name_t>::type
    gvname = get(vertex_name,T.second);
    
    property_map<Graph,edge_name_t>::type
    gename = get(edge_name,T.second);
    
    // We first need to find initial types and find the maximum size of the inital blocks
    // After, we can create the multi-matrix of size/type. However, we also need to check that
    // each vertex above has been given homomorphism data because you are lazy about that.
    
    for(tie(gvi,gvi_end)=vertices(T.second); gvi!=gvi_end; gvi++)
    {
        sColl currSColl;
        
        types.push_back(gvname(*gvi));
        
	cout << "We have inserted " << gvname(*gvi) << " into types." << endl;
        
        for(tie(goei,goei_end)=out_edges(*gvi,T.second); goei!=goei_end; goei++)
	  {
            currSColl.push_back(gename(*goei));
            cout << "We have pushed " << gename(*goei) << " onto the current string collection" << endl;
	  } // inner for
        productionRules.push_back(currSColl);
    } // outer for
    
    
    // We have all the types and their production rules, now we iterate through and create 
    // the initial vColl for each type
    for(it = types.begin(), k=0; it != types.end(); it++, k++)
    {
        vColl *currColl = new vColl;
        vector<vColl*> *dummy = new vector<vColl*>(1);
        int M=0;
        
        cout << "Looking at type: " << *it << endl;
        
        
        // iterate through the upstairs vertices and find all those whose vertex homom matches
        // the current type
        for(tie(vi,vi_end)=vertices(T.first); vi!=vi_end; vi++)
        {
            if(gvname(p_vhomom(*vi))==*it)
            {
                // if it matches, we insert it and increase the M, the size of the initial vColl of this type
                currColl->insert(*vi);
                M++;
                
  //              cout << "Found a match, new m value is: " << M << endl;
            }
 //           cout << "We are done checking vertex: " << *vi << endl;
        } // inner for
        
  //      cout << "Out of for loop" << endl;
        if(M>maxM) {
            cout << "New maxM is " << M << endl;
            typeMatrix.resize(M,typeNumber,true);
            cout << "Resizing complete." << endl;
            for(i=maxM;i<M;i++)
            {
                for(j=0;j<typeNumber;j++)
                {
                    typeMatrix(i,j) = 0;		
                }
            }
            
            
            maxM = M;
        }
        
        (*dummy)[0]=currColl;
        
            cout << "Dummy[0] set and maxM-M is " << maxM-M << " while curr typeNumber is " << k << endl;
        
        typeMatrix(M-1,k)=dummy;
        
              cout << "successfully initialized dummy" << endl;
    } // outer for
    
    cout << "Moving on with a maxM of " << maxM << endl;
    
    // Now that the type matrix has been computed, we want to iterate though each level, apply
    // our production rules, and if our production rules lead us to a set of vertices of the same size,
    // we want to enter it into the graph of that level. If adding that edge gives us a cycle, then
    // we know that our graph is not definite.
    for(i = maxM-1; i != 0; i--)
    {
        cout << "Starting to look at vertex collections of size " << i+1 << endl;
        
        for(j = 0; j < typeNumber; j++)
        {
            if(typeMatrix(i,j)!=0)
            {
                cout << "Looking at type " << j << endl;
                for(m = 0; m < typeMatrix(i,j)->size(); m++)
                {
                    GVD oldv;
                    unsigned int p;
                    vector<vColl*>::iterator vvcit;
                    
                    vColl* currVColl = (*typeMatrix(i,j))[m];
		    cout << "currVColl.empty? " << currVColl->empty() << endl;
                    // Push all of the production rules for this type onto the stack
                    sColl currentSColl = productionRules[j];
                    bool grfound=false;
                    
                    for(vvcit=graphRef.begin(), p=0; vvcit!=graphRef.end(); vvcit++, p++)
                    {
                        if(*(*vvcit)== *currVColl)
                        {
			  cout << "Found a match for our vvcit: " << p << endl;
			  oldv=p;
			  grfound = true;
                        }
                    }
                    if(!grfound)
                    {
                        oldv = add_vertex(matrixHelper);  
                        graphRef.push_back(currVColl);
                        
                    }
                    
                    for(k = 0; k < currentSColl.size(); k++)
                    {
		      cout << "Pushing rule " << j << "," << currentSColl[k] << endl;
                        prodStack.push(sPair(j,currentSColl[k]));
                    }
                    
                    // Apply the the production rules until we're empty.
                    while(!prodStack.empty())
                    {
                        vColl::iterator vit;
                        
                        firstfound = false;
                        
                        sPair currPair = prodStack.top();
                        prodStack.pop();
                        
                        //	      cout << "Got pair " << currPair.first << "," << currPair.second << endl;
                        
                        // We now create the new vertex collection
                        vColl * newVColl = new vColl;
                        GVD newType;
                        
                        for(vit = currVColl->begin(); vit != currVColl->end(); vit++)
                        {
                            //found = false;
                            //		  cout << "Looking at element " << *vit << " of our current vertex collection" << endl;
                            
                            // iterate through the out_edges looking for a match
                            for(tie(oei,oei_end)=out_edges(*vit,T.first); oei != oei_end; oei++)
                            {
                                //		      cout << "Current phom is " << p_homom(*oei) << endl;
                                // we found a match for our rule
                                if(p_homom(*oei) == currPair.second)
                                {
                                    //		  cout << "We found a match: " << target(*oei,T.first) <<  endl;
                                    //found = true;
                                    newVColl->insert(target(*oei,T.first));
                                    if(!firstfound) {
                                        newType = p_vhomom(target(*oei,T.first));
                                        firstfound = true;
                                    } // firstfound
                                }
                            } // oei for
                            //} // while found
                        } // vit for
                        
                        // We have the new vColl, now we need to find its size and add it to the matrix and graph.
                        newM = newVColl->size();
                        //	      cout << "newM is " << newM << " and newType is " << newType << endl;
                        //	      cout << "typeMatrix is " << typeMatrix.size1() << " by " << typeMatrix.size2() << endl;
                        if(newM!=0){
                            if(typeMatrix(newM-1,newType)==0){
                                //		  cout << "Target is 0, putting in a new vector<vColl>" << endl;
                                typeMatrix(newM-1,newType)= new vector< vColl *>;
                            }
                            // We need to check to make sure that newVColl isn't already in there
                            bool tmfound= false;
                            
                            for(vvcit = (*typeMatrix(newM-1,newType)).begin(); vvcit != (*typeMatrix(newM-1,newType)).end(); vvcit++)
                            {
                                if(*(*vvcit) == *newVColl) {
                                    tmfound = true;
                                    break;
                                }
                            }
                            if(!tmfound) {
                                typeMatrix(newM-1,newType)->push_back(newVColl);
                            }
                            // if they are the same size, we need to create an edge
                            if(newM == i+1)
                            {
                                bool vfound=false;
                                //	    stringstream newName;
                                GVD newv;
                                //	    newName << newType << "," << typeMatrix(newM-1,newType)->size()-1;
                                
                                //         cout << newName.str() << endl;
                                // check to see if vertex already exists
                                for(tie(gvi,gvi_end) = vertices(matrixHelper); gvi != gvi_end; gvi++)
                                {
                                    if(*graphRef[*gvi] == *newVColl) {
                                        vfound = true;
                                        newv = *gvi;
                                    }
                                    
                                }
                                if(!vfound)
                                {
                                    //             cout << "Could not find the vertex, creating a new one" << endl;
                                    newv = add_vertex(matrixHelper);
                                    graphRef.push_back(newVColl);
                                }
                                add_edge(oldv,newv,currPair.second,matrixHelper);
                                //	    cout << "We have created an edge between " << oldv
                                //		 << " and " << newv << endl;
                            } // newM == i+1
                        } // newM != 0
                    } //while (all rules for a given vColl)
                } // m for (all vColls for a given type and size)
            } // if
            else {
           //     cout << "We are 0" << endl;
            } // else
        } // j for (all types of a given size)
        // We now need to check to see if any cycles developed in matrixHelper over the course of the algorithm
        //   cout << "Checking to see if we have cycles." << endl;
        //    cout << "MatrixHelper has " << num_vertices(matrixHelper) << " vertices and " << num_edges(matrixHelper)
        //   << " edges." << endl;
        vector<colors> color(num_vertices(matrixHelper),White);
        for(tie(gvi,gvi_end)=vertices(matrixHelper); gvi!=gvi_end; gvi++)
        {
            
            if(color[*gvi] == White) {
                cout << "Vertex " << *gvi << " is white" << endl;
                if(hasCycleHelper(matrixHelper,*gvi,&color[0]))
                    return 0;
                
            }
        }
        
    } // i for 
    
    cout << "We are definite" << endl;
    return 1;
    
}

bool hasCycleHelper(Graph& g, GVD u, colors * color)
{
    color[u]=Gray;
    graph_traits<Graph>::adjacency_iterator vi,vi_end;
    for(tie(vi,vi_end)=adjacent_vertices(u,g); vi!=vi_end; ++vi)
    {
        cout << "Looking at vertex " << *vi << " which is connected to " << u  << " and is colored " << color[*vi] << endl;
        if(color[*vi] == White) {
            if(hasCycleHelper(g,*vi,color)) {
                return true;
            }
        }
        else if(color[*vi]==Gray)
            return true; 
    }
    color[u] = Black;
    return false;
} 


int IspLeftDefinite(Textile T)
{
    Graph matrixHelper;
    vector<vColl*> graphRef;
    int i=0,j=0,k=0,m=0,newM=0,maxM=0,typeNumber=num_vertices(T.second);
    graph_traits<Graph>::vertex_iterator gvi,gvi_end;
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end;
    graph_traits<GammaGraph>::in_edge_iterator oei,oei_end;
    graph_traits<Graph>::in_edge_iterator goei,goei_end;
    bool done, found, firstfound;
    stack<sPair> prodStack;
    
    // The production rules depend only on the type. They can be generated by looking at the out edges
    // of the downstairs vertex the type represents.
    vector<sColl> productionRules;
    
    vector<string> types;
    vector<string>::iterator it;

    cout << "About to create type matrix" << endl;
    
    // The type matrix is indexed so that the 0 index is actually representing the types of
    // size M (where M is the largest vertex set).
    boost::numeric::ublas::matrix<vector<vColl* >* > typeMatrix;
    
    cout << "Created type matrix" << endl;

    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_vhomom = get(vertex_p_vhomom,T.first);
    
    property_map<Graph,vertex_name_t>::type
    mh_vname = get(vertex_name,matrixHelper);
    
    property_map<Graph,vertex_name_t>::type
    gvname = get(vertex_name,T.second);
    
    property_map<Graph,edge_name_t>::type
    gename = get(edge_name,T.second);
    
    // We first need to find initial types and find the maximum size of the inital blocks
    // After, we can create the multi-matrix of size/type. However, we also need to check that
    // each vertex above has been given homomorphism data because you are lazy about that.
    
    for(tie(gvi,gvi_end)=vertices(T.second); gvi!=gvi_end; gvi++)
    {
        sColl currSColl;
        
        types.push_back(gvname(*gvi));
        
        cout << "We have inserted " << gvname(*gvi) << " into types." << endl;
        
        for(tie(goei,goei_end)=in_edges(*gvi,T.second); goei!=goei_end; goei++)
        {
            currSColl.push_back(gename(*goei));
            cout << "We have pushed " << gename(*goei) << " onto the current string collection" << endl;
        } // inner for
        productionRules.push_back(currSColl);
    } // outer for
    
    
    // We have all the types and their production rules, now we iterate through and create 
    // the initial vColl for each type
    for(it = types.begin(), k=0; it != types.end(); it++, k++)
    {
        vColl *currColl = new vColl;
        vector<vColl*> *dummy = new vector<vColl*>(1);
        int M=0;
        
        cout << "Looking at type: " << *it << endl;
        
        
        // iterate through the upstairs vertices and find all those whose vertex homom matches
        // the current type
        for(tie(vi,vi_end)=vertices(T.first); vi!=vi_end; vi++)
        {
            if(gvname(p_vhomom(*vi))==*it)
            {
                // if it matches, we insert it and increase the M, the size of the initial vColl of this type
                currColl->insert(*vi);
                M++;
                
                cout << "Found a match, new m value is: " << M << endl;
            }
            cout << "We are done checking vertex: " << *vi << endl;
        } // inner for
        
   //     cout << "Out of for loop" << endl;
        if(M>maxM) {
            cout << "New maxM is " << M << endl;
            typeMatrix.resize(M,typeNumber,true);
            cout << "Resizing complete." << endl;
            for(i=maxM;i<M;i++)
            {
                for(j=0;j<typeNumber;j++)
                {
                    typeMatrix(i,j) = 0;		
                }
            }
            
            
            maxM = M;
        }
        
        (*dummy)[0]=currColl;
        
        cout << "Dummy[0] set and maxM-M is " << maxM-M << " while curr typeNumber is " << k << endl;
        
        typeMatrix(M-1,k)=dummy;
        
        cout << "successfully initialized dummy" << endl;
    } // outer for
    
    cout << "Moving on with a maxM of " << maxM << endl;
    
    // Now that the type matrix has been computed, we want to iterate though each level, apply
    // our production rules, and if our production rules lead us to a set of vertices of the same size,
    // we want to enter it into the graph of that level. If adding that edge gives us a cycle, then
    // we know that our graph is not definite.
    for(i = maxM-1; i != 0; i--)
    {
      cout << "Starting to look at vertex collections of size " << i+1 << endl;
      
      for(j = 0; j < typeNumber; j++)
        {
	  if(typeMatrix(i,j)!=0)
            {
	      cout << "Looking at type " << j << endl;
	      for(m = 0; m < typeMatrix(i,j)->size(); m++)
                {
		  GVD oldv;
		  unsigned int p;
		  vector<vColl*>::iterator vvcit;
                  
		  vColl* currVColl = (*typeMatrix(i,j))[m];
		  cout << "currVColl.empty? " << currVColl->empty() << endl;
		  // Push all of the production rules for this type onto the stack
		  sColl currentSColl = productionRules[j];
		  bool grfound=false;
                  
		  for(vvcit=graphRef.begin(), p=0; vvcit!=graphRef.end(); vvcit++, p++)
                    {
		      if(*(*vvcit)==*currVColl)
                        {
			  cout << "Found a match for our vvcit: " << p << endl;
			  oldv=p;
			  grfound = true;
                        }
                    }
		  if(!grfound)
                    {
		      oldv = add_vertex(matrixHelper);
		      graphRef.push_back(currVColl);
		      cout << "Did not find a match for vvcit, oldv set as " << oldv << endl;
                    }
		  
		  for(k = 0; k < currentSColl.size(); k++)
                    {
		      cout << "Pushing rule " << j << "," << currentSColl[k] << endl;
		      prodStack.push(sPair(j,currentSColl[k]));
                    }
		  
		  // Apply the the production rules until we're empty.
		  while(!prodStack.empty())
                    {
		      vColl::iterator vit;
                      
		      firstfound = false;
                      
		      sPair currPair = prodStack.top();
		      prodStack.pop();
                      
		      cout << "Got pair " << currPair.first << "," << currPair.second << endl;
                      
		      // We now create the new vertex collection
		      vColl * newVColl = new vColl;
		      GVD newType;
                      
		      for(vit = currVColl->begin(); vit != currVColl->end(); vit++)
			
                        {
			  //found = false;
			  cout << "Looking at element " << *vit 
			       << " of our current vertex collection" << endl;
			  
			  // iterate through the in_edges looking for a match
			  for(tie(oei,oei_end)=in_edges(*vit,T.first); oei != oei_end; oei++)
                            {
			      //	      cout << "Current phom is " << p_homom(*oei) << endl;
			      // we found a match for our rule
			      if(p_homom(*oei) == currPair.second)
				{
				  cout << "We found a match: " << target(*oei,T.first) <<  endl;
				  //found = true;
				  newVColl->insert(source(*oei,T.first));
				  if(!firstfound) {
				    newType = p_vhomom(source(*oei,T.first));
				    firstfound = true;
				  } // firstfound
				}
                            } // oei for
			  //} // while found
                        } // vit for
		      
		      // We have the new vColl, now we need to find its size 
		      // and add it to the matrix and graph.
		      newM = newVColl->size();
		      cout << "newM is " << newM << " and newType is " << newType << endl;
		      cout << "typeMatrix is " << typeMatrix.size1() 
			   << " by " << typeMatrix.size2() << endl;
		      if(newM!=0){
			if(typeMatrix(newM-1,newType)==0){
			  cout << "Target is 0, putting in a new vector<vColl>" << endl;
			  typeMatrix(newM-1,newType)= new vector< vColl *>;
                          
			}
                        
			// We need to check to make sure that newVColl isn't already in there
			bool tmfound= false;
                        
			for(vvcit = (*typeMatrix(newM-1,newType)).begin();
			    vvcit != (*typeMatrix(newM-1,newType)).end();
			    vvcit++)
			  {
			    if(*(*vvcit) == *newVColl) {
			      tmfound = true;
			      break;
			    }
			  }
			if(!tmfound) {
			  typeMatrix(newM-1,newType)->push_back(newVColl);
			}
			// if they are the same size, we need to create an edge
			if(newM == i+1)
			  {
			    bool vfound=false;
			    //	    stringstream newName;
			    GVD newv;
			    // check to see if vertex already exists
			    for(tie(gvi,gvi_end) = vertices(matrixHelper); gvi != gvi_end; gvi++)
			      {
				if(*graphRef[*gvi] == *newVColl) {
				  vfound = true;
				  newv = *gvi;
				}
                                
			      }
			    
			    if(!vfound)
			      {
				//		  cout << "Could not find the vertex, creating a new one" << endl;
				newv = add_vertex(matrixHelper);
				graphRef.push_back(newVColl);
			      }
			    add_edge(newv,oldv,currPair.second,matrixHelper);
			    cout << "We have created an edge between " << newv
				 << " and " << oldv << endl;
			  } // newM == i
		      } // newM != 0
                    } //while (all rules for a given vColl)
                } // m for (all vColls for a given type and size)
            } // if
            else {
        //        cout << "We are 0" << endl;
            } // else
        } // j for (all types of a given size)
        // We now need to check to see if any cycles developed in matrixHelper over the course of the algorithm
        //   cout << "Checking to see if we have cycles." << endl;
        //   cout << "MatrixHelper has " << num_vertices(matrixHelper) << " vertices and " << num_edges(matrixHelper)
        //   << " edges." << endl;
        vector<colors> color(num_vertices(matrixHelper),White);
        for(tie(gvi,gvi_end)=vertices(matrixHelper); gvi!=gvi_end; gvi++)
        {
            if(color[*gvi] == White)
                if(hasCycleHelper(matrixHelper,*gvi,&color[0]))
                    return -1;
        }
        
    } // i for 
    
    cout << "We are definite" << endl;
    return 1;
}

int IsqLeftDefinite(Textile T)
{
    return IspLeftDefinite(CreateInverse(T));
}



void PrintRepMatrix(Textile T)
{
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


// This will rename the edges of the G graph and the p,q maps in Gamma with whatever
// you fed in as a map of strings
Textile RenameTextile(Textile T,map<string,string> names)
{
    GammaGraph A(num_vertices(T.first));
    Graph B(num_vertices(T.second));
    
    graph_traits<GammaGraph>::edge_iterator ei,ei_end;
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end,gi,gi_end;
    graph_traits<Graph>::edge_iterator gei,gei_end;
    graph_traits<Graph>::vertex_iterator gvi,gvi_end,gwi,gwi_end;
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,T.first);
    
    property_map<GammaGraph,edge_name_t>::type
    ename = get(edge_name,T.first);
    
    property_map<GammaGraph,vertex_name_t>::type
    gammavname = get(vertex_name,T.first);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_avhom = get(vertex_p_vhomom,A);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_avhom = get(vertex_q_vhomom,A);  
    
    property_map<GammaGraph,vertex_name_t>::type
    avname = get(vertex_name,A);
    
    property_map<Graph,vertex_name_t>::type
    bvname = get(vertex_name,B);
    
    property_map<Graph,vertex_name_t>::type
    gvname = get(vertex_name,T.second);
    
    property_map<Graph,edge_name_t>::type
    gename = get(edge_name,T.second);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_gammavhom = get(vertex_p_vhomom,T.first);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_gammavhom = get(vertex_q_vhomom,T.first);
    
    for(tie(ei,ei_end)=edges(T.first);ei!=ei_end;ei++)
    {
        add_edge(source(*ei,T.first),target(*ei,T.first),PQ_Homoms(names[p_homom(*ei)],Q_Homom(names[q_homom(*ei)],ename(*ei))),A);
    }
    
    for(tie(gi,gi_end)=vertices(T.first),tie(vi,vi_end)=vertices(A);gi!=gi_end;gi++,vi++)
    {
        put(p_avhom,*vi,p_gammavhom(*gi));
        put(q_avhom,*vi,q_gammavhom(*gi));
        put(avname,*vi,gammavname(*gi));
    }
    
    for(tie(gei,gei_end)=edges(T.second);gei!=gei_end;gei++)
    {
        add_edge(source(*gei,T.second),target(*gei,T.second),names[gename(*gei)],B);
    }
    
    for(tie(gvi,gvi_end)=vertices(T.second), tie(gwi,gwi_end)=vertices(B);gvi!=gvi_end;gvi++,gwi++)
    {
        put(bvname,*gwi,gvname(*gvi));
    }
    
    return Textile(A,B);
    
    
}

Textile HigherNBlock(Textile T, int n)
{ 
    GammaGraph A;
    Graph B;
    
    int i=0;
    
    string currentpWord,currentqWord,currenteWord,first,last,newpWord,newqWord,neweWord;
    
    // This map stores n strings and their terminal vertices. Before inserting a new
    // string, we will check to see if the string is already in the map. If so, we check
    // their terminal vertices. If they match, we do not add the string, if not, return false
    map<string,graph_traits<GammaGraph>::vertex_descriptor> GammaVertexMap;
    map<string,graph_traits<Graph>::vertex_descriptor> GVertexMap;
    //  map<string,string>::iterator tmi;
    
    // The OEI pairs allow us to keep both the current location of the iterator, its end
    // and the string it is carrying
    typedef std::tuple<graph_traits<GammaGraph>::out_edge_iterator,
    graph_traits<GammaGraph>::out_edge_iterator,
    string,string,string> GammaHOEITuple;
    
    typedef std::tuple<graph_traits<Graph>::out_edge_iterator,
    graph_traits<Graph>::out_edge_iterator,
    string> GHOEITuple;
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,T.first);
    
    property_map<Graph,edge_name_t>::type
    g_name = get(edge_name,T.second);
    
    property_map<GammaGraph,vertex_name_t>::type
    gamma_vname = get(vertex_name,T.first);
    
    property_map<GammaGraph,vertex_name_t>::type
    a_vname = get(vertex_name,A);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    a_pvhom = get(vertex_p_vhomom,A);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    a_qvhom = get(vertex_q_vhomom,A);  
    
    property_map<Graph,vertex_name_t>::type
    b_vname = get(vertex_name,B);
    
    property_map<GammaGraph,edge_name_t>::type
    ename = get(edge_name,T.first);
    
    property_map<Graph,edge_name_t>::type
    gename = get(edge_name,T.second);
    
    property_map<Graph,vertex_name_t>::type
    g_vname = get(vertex_name,T.second);
    
    GammaVI vi,vi_end;
    graph_traits<Graph>::vertex_iterator gvi,gvi_end;
    
    GammaOEI oi,oi_end,currenti,currenti_end;
    graph_traits<Graph>::out_edge_iterator goi,goi_end,gcurrenti,gcurrenti_end;
    
    stack<GammaHOEITuple> GammaHOEIStack;
    stack<GHOEITuple> GHOEIStack;
    
    int ChunkSize,BlockSize,NameChunkSize,NameBlockSize,HomChunkSize,HomBlockSize;
    
    bool done= false;
    if(n <= 1)
    {
        cout << "n must be greater than 1, it is " << n << endl;
        return T;
    }
    
    for(tie(gvi,gvi_end)=vertices(T.second); gvi!=gvi_end; gvi++)
    {
        i=1;
        tie(goi,goi_end)=out_edges(*gvi,T.second);
        
        ChunkSize = (gename(*goi)).size();
        BlockSize = ChunkSize*(n-1);
        
        //    cout << "We've started looking at vertex " << g_vname(*gvi) << endl;
        
        GHOEIStack.push(GHOEITuple(goi,goi_end,""));
        
        while(!GHOEIStack.empty())
        {
            GHOEITuple current = GHOEIStack.top();
            
            gcurrenti = get<0>(current);
            gcurrenti_end = get<1>(current);
            currenteWord = get<2>(current);
            
            boost::graph_traits<Graph>::edge_descriptor e=*gcurrenti;
            
            if(i < n) // We want to go deeper
            {
                graph_traits<Graph>::out_edge_iterator ni,ni_end;
                
                tie(ni,ni_end)=out_edges(target(e,T.second),T.second);
                GHOEIStack.push(GHOEITuple(ni,ni_end,currenteWord+gename(e)));
                i++;
                //	      cout << "i value is " << i << endl;
                
                //	      cout << "Our new edge word is " << currenteWord+gename(e) << endl;
            }
            else // We're the deepest we want to go. 
            {
                graph_traits<Graph>::vertex_descriptor s,t;
                map<string,graph_traits<Graph>::vertex_descriptor>::iterator si,ti;
                GHOEIStack.pop();
                neweWord = currenteWord+gename(e);
                //    cout << "Our final length " << n << " word is " << neweWord << endl;
                first = string(neweWord,0,BlockSize);
                last = string(neweWord,ChunkSize,BlockSize);
                si = GVertexMap.find(first);
                
                if(si==GVertexMap.end())
                {
                    s = add_vertex(B);
                    put(b_vname,s,first);
                    GVertexMap[first]=s;
                    //   cout << "We have created the vertex " << first << endl;
                }
                else {
                    s = (*si).second; 
                }
                ti = GVertexMap.find(last);
                if(ti==GVertexMap.end())
                {
                    t = add_vertex(B);
                    put(b_vname,t,last);
                    GVertexMap[last]=t;
                    //  cout << "We have created the vertex " << last << endl;
                }
                else {
                    t = (*ti).second;
                }	      
                
                add_edge(s,t,neweWord, B);
                gcurrenti++;
                if(gcurrenti!=gcurrenti_end) {
                    GHOEIStack.push(GHOEITuple(gcurrenti,gcurrenti_end,currenteWord));
                }
                else { // Our top iterator has finished. We need to go down a level (or two..)
                    done = false;
                    while(!done && !GHOEIStack.empty()) {
                        current = GHOEIStack.top();
                        GHOEIStack.pop();
                        gcurrenti = get<0>(current);
                        gcurrenti_end = get<1>(current);
                        currenteWord = get<2>(current);
                        gcurrenti++;
                        i--;
                        if(gcurrenti!=gcurrenti_end)
                        {
                            done=true;
                            GHOEIStack.push(GHOEITuple(gcurrenti,gcurrenti_end,currenteWord));
                        }
                        
                    }
                }
            }
            
        }
    }
    
    cout << "WE ARE DONE WITH THE G GRAPH" << endl;
    
    
    
    for(tie(vi,vi_end)=vertices(T.first); vi!=vi_end; vi++)
    {
        i=1;
        tie(oi,oi_end)=out_edges(*vi,T.first);
        
        NameChunkSize = (ename(*oi)).size();
        NameBlockSize = NameChunkSize*(n-1);
        HomChunkSize = (p_homom(*oi)).size();
        HomBlockSize = HomChunkSize*(n-1);
        /*     
         cout << "name chunk size = " << NameChunkSize << ", name block size = " << NameBlockSize
         << ", hom chunk size = " << HomChunkSize << ", hom block size = " << HomBlockSize << endl;
         */   
     //   cout << "We've started looking at vertex " << gamma_vname(*vi) << endl;
        
        GammaHOEIStack.push(GammaHOEITuple(oi,oi_end,"","",""));
        
        while(!GammaHOEIStack.empty())
        {
            GammaHOEITuple current = GammaHOEIStack.top();
            
            currenti = get<0>(current);
            currenti_end = get<1>(current);
            currentpWord = get<2>(current);
            currentqWord = get<3>(current);
            currenteWord = get<4>(current);
            
            boost::graph_traits<GammaGraph>::edge_descriptor e=*currenti;
            
            if(i < n) // We want to go deeper
            {
                graph_traits<GammaGraph>::out_edge_iterator ni,ni_end;
                i++;
                tie(ni,ni_end)=out_edges(target(e,T.first),T.first);
                GammaHOEIStack.push(GammaHOEITuple(ni,ni_end,currentpWord+p_homom(e),currentqWord+q_homom(e),currenteWord+ename(e)));
                
                //	     cout << "Our new edge word is " << currenteWord+ename(e) << endl;
                //  cout << "Our new p word is " << currentpWord+p_homom(e) << endl;	      
                // cout << "Our new q word is " << currentqWord+q_homom(e) << endl;
                // cout << "i value is " << i << endl;
            }
            else // We're the deepest we want to go. 
            {
                VD s,t;
                map<string,VD>::iterator si,ti;
                GammaHOEIStack.pop();
                neweWord = currenteWord+ename(e);
                newpWord = currentpWord+p_homom(e);
                newqWord = currentqWord+q_homom(e);
                //            cout << "Our final length " << n << " word is " << neweWord << " with associated p and q as "
                //       << newpWord << " and " << newqWord << endl;
                first = string(neweWord,0,NameBlockSize);
                last = string(neweWord,NameChunkSize,NameBlockSize);
                si = GammaVertexMap.find(first);
                
                if(si==GammaVertexMap.end()) // We did not find first, so we need to create its vertex
                {
                    s = add_vertex(A);
                    put(a_vname,s,first);
                    map<string,GVD>::iterator ploc=GVertexMap.find(string(newpWord,0,HomBlockSize)),
                    qloc=GVertexMap.find(string(newqWord,0,HomBlockSize));
                    if(ploc==GVertexMap.end() || qloc==GVertexMap.end()) {
                                         cout << "Something very bad happened and we can't find a p or q word for "
                                         << string(newpWord,0,HomBlockSize) << " or " << string(newqWord,0,HomBlockSize) << endl;
                    }
                    put(a_pvhom,s,(*ploc).second);
                    put(a_qvhom,s,(*qloc).second);
              //      cout << "We have created the vertex " << first << endl;
                    GammaVertexMap[first]=s;
                }
                else {
                    s = (*si).second; 
                }
                ti = GammaVertexMap.find(last);
                if(ti==GammaVertexMap.end()) // We did not find last, so we need to create its vertex
                {
                    t = add_vertex(A);
                    put(a_vname,t,last);
                    map<string,graph_traits<Graph>::vertex_descriptor>::iterator ploc=GVertexMap.find(string(newpWord,HomChunkSize,HomBlockSize)),
                    qloc=GVertexMap.find(string(newqWord,HomChunkSize,HomBlockSize));
                    if(ploc==GVertexMap.end() || qloc==GVertexMap.end()) {
                                         cout << "Something very bad happened and we can't find a p or q word" << endl;
                    }
                    put(a_pvhom,t,(*ploc).second);
                    put(a_qvhom,t,(*qloc).second);
            //        cout << "We have created the vertex " << last << endl;
                    GammaVertexMap[last]=t;
                }
                else {
                    t = (*ti).second;
                }	      
                
                add_edge(s,t,PQ_Homoms(newpWord,Q_Homom(newqWord, neweWord)), A);
            //                cout << "We have created the edge from " << a_vname(s) << " to " << a_vname(t) << " with name " << neweWord 
            //                << " and p and q words " << newpWord << " and " << newqWord << endl;
                currenti++;
                if(currenti!=currenti_end) {
                    GammaHOEIStack.push(GammaHOEITuple(currenti,currenti_end,currentpWord,currentqWord,currenteWord));
                }
                else { // Our top iterator has finished. We need to go down a level (or two..)
                    done = false;
                    while(!done && !GammaHOEIStack.empty()) {
                        current = GammaHOEIStack.top();
                        GammaHOEIStack.pop();
                        currenti = get<0>(current);
                        currenti_end = get<1>(current);
                        currentpWord = get<2>(current);
                        currentqWord = get<3>(current);
                        currenteWord = get<4>(current);
                        currenti++;
                        i--;
                        if(currenti!=currenti_end)
                        {
                            done=true;
                            GammaHOEIStack.push(GammaHOEITuple(currenti,currenti_end,currentpWord,currentqWord,currenteWord));
                        }
                        
                    }
                }
            }
            
        }
    }
    
    cout << "Done with creating Higher N block" << endl;
    
    return Textile(A,B);
}

// NOTE: Quotient requires a Textile with named vertices.
Textile Quotient(Textile T, vector<vector<VD> > E)
{
    Graph G = T.second;
    GammaGraph Gamma(E.size());
    
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end;
    graph_traits<GammaGraph>::out_edge_iterator oei,oei_end,ofi,ofi_end;
    vector<graph_traits<GammaGraph>::vertex_descriptor> currE;
    VD va,vb;
    int i,j,k,m,l,ea,eb;
    bool done,founda,foundb,found;
    string name;
    
    property_map<GammaGraph,vertex_name_t>::type
    gamma_vname = get(vertex_name,Gamma);
    
    property_map<GammaGraph,vertex_name_t>::type
    Tf_vname = get(vertex_name,T.first);
    
    property_map<GammaGraph,edge_p_homom_t>::type
    Tf_p = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    Tf_q = get(edge_q_homom,T.first);
    
    property_map<GammaGraph,edge_name_t>::type
    Tf_e = get(edge_name,T.first);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_vhom= get(vertex_p_vhomom,Gamma);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_vhom = get(vertex_q_vhomom,Gamma);

    property_map<GammaGraph,vertex_p_vhomom_t>::type
    Tp_vhom= get(vertex_p_vhomom,T.first);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    Tq_vhom = get(vertex_q_vhomom,T.first);
    
    
    VertexMap vmap;
    
    cout << "Printing E, the state equivalences" << endl;
    for(int i=0; i<E.size();i++)
    {
        printVVec(E[i]);
        cout << endl;
    }
    
    // We need to check to make sure our equivalence relation E is actually a state equivalence    
     for(int i=0;i<E.size();i++)
     {
         if( (E[i]).size() > 1) // Our state equivalence actually contains more than one element. If it doesn't, we can move on
         {
             for(j=0; j< E[i].size()-1; j++)
             {
                 for(k=j; k<E[i].size(); k++)
                 {
                     // We want to check that we satisfy condition 1 of the definition of a state equiv.
                     for(tie(oei,oei_end)=out_edges(E[i][j],T.first); oei != oei_end; oei++)
                     {
                         found = false;
                         for(tie(ofi,ofi_end)=out_edges(E[i][k],T.first); ofi != ofi_end; ofi++)
                         {
                             if(Tf_p(*ofi)==Tf_p(*oei) && Tf_q(*ofi) == Tf_q(*oei))
                             {
                                 found = true;
                                 va = target(*oei,T.first);
                                 vb = target(*ofi,T.first);
                                 founda = false; foundb = false;
                                 // Iterate through the classes
                                 for(m=0;m<E.size() && !(founda && foundb); m++)
                                 {
                                    // Iterate within the classes
                                     for(l=0;l<E[m].size();l++)
                                     {
                                         if(E[m][l]==va) // We've found the class, set ea to it.
                                         {
                                             ea = m;
                                             founda = true;
                                         }
                                         if(E[m][l]==vb)
                                         {
                                             eb = m;
                                             foundb = true;
                                         }
                                     } // for l
                                 } // for m
                                 if(ea!=eb)
                                 {
                                     cerr << "NOT A STATE EQUIVALENCE" << endl;
                                     printVVec(E[i]);
                                     cout << endl;
                                     SmartPrintTextileInfo(T);
                                     return T;
                                 }
                             } // if
                         } // for ofi
                         if(!found)
                         {
                             cerr << "NOT A STATE EQUIVALENCE" << endl;
                             printVVec(E[i]);
                             cout << endl;
                             SmartPrintTextileInfo(T);
                             return T; 
                         } 
                     } // for oei

                     // We repeat the same check as above, but checking the other direction
                     for(tie(oei,oei_end)=out_edges(E[i][k],T.first); oei != oei_end; oei++)
                     {
                         found = false;
                         for(tie(ofi,ofi_end)=out_edges(E[i][j],T.first); ofi != ofi_end; ofi++)
                         {
                             if(Tf_p(*ofi)==Tf_p(*oei) && Tf_q(*ofi) == Tf_q(*oei))
                             {
                                 found = true;
                             } // if
                         } // for ofi
                         if(!found)
                         {
                             cerr << "NOT A STATE EQUIVALENCE" << endl;
                             printVVec(E[i]);
                             cout << endl;
                             SmartPrintTextileInfo(T);
                             return T; 
                         } 
                     } // for oei

                 } // for k
             } // for j
         } // if size > 1
     }
     
    
    for(tie(vi,vi_end)=vertices(Gamma),i=0; vi!=vi_end; vi++, i++)
    {
        currE=E[i];
   //     cout << "Naming the " << i << "th vertex" << endl;
        put(gamma_vname,*vi,Tf_vname(currE[0]));
        vmap.insert(VertexMap::value_type(gamma_vname(i),i));
    //    cout << "Assigning vhoms" << endl;
        put(p_vhom,*vi,Tp_vhom(currE[0]));
        put(q_vhom,*vi,Tq_vhom(currE[0]));
    }
    
    cout << "Made it through the naming process" << endl;
    
    // Now we can actually create our quotient system
    for(tie(vi,vi_end)=vertices(Gamma),i=0; vi!=vi_end; vi++, i++)
    {
        map<pair<graph_traits<GammaGraph>::vertex_descriptor,string>,bool> mmap;
        map<pair<graph_traits<GammaGraph>::vertex_descriptor,string>,bool>::iterator it;
        
        
    //    cout << "Starting the " << i << "th round" << endl;
        
        currE=E[i];
        
        for(j=0;j<currE.size();j++)
        {
   //         cout << "Looking at currE[" << j << "], which is " << currE[j] << endl;
            for(tie(oei,oei_end)=out_edges(currE[j],T.first);oei!=oei_end;oei++)
            {
                graph_traits<GammaGraph>::vertex_descriptor t=target(*oei,T.first),Et;
                done = false;
       //         cout << "We are trying to find the repn for " << t << endl;
                for(k=0;k<E.size() && !done;k++)
                {
                    vector<graph_traits<GammaGraph>::vertex_descriptor> searchE=E[k];
                    for(m=0;m<searchE.size() && !done;m++)
                    {
                        if(t==searchE[m])
                        {
          //                  cout << "We found it! It's in the " << k << "th equiv class" << endl;
                            Et=searchE[0]; // This gives us the equivalence class representative of the target
                            done = true;
                        } // if 
                    } // for m
                    
                } // for k
                // Have we already added this? 
                it=mmap.find(pair<VD,string>(vmap[Tf_vname(Et)],Tf_p(*oei)+"/"+Tf_q(*oei)));
                if(it == mmap.end()) 
                {
                    add_edge(*vi,vmap[Tf_vname(Et)],PQ_Homoms(Tf_p(*oei),Q_Homom(Tf_q(*oei),Tf_vname(*vi) + Tf_p(*oei) + "/" + Tf_q(*oei))),Gamma);
                    mmap[pair<VD,string>(vmap[Tf_vname(Et)],Tf_p(*oei)+"/"+Tf_q(*oei))]=true;
                } // if it
            } // for oei
        } // for j
        
    } // for vi
    return Textile(Gamma,G);
    
}


Textile HomomComp(Textile T,Graph G,map<string,string> h)
{
    GammaGraph Gamma=T.first;
    string p,q,np,nq;
    
    
    graph_traits<GammaGraph>::edge_iterator ei,ei_end;
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end;
    graph_traits<Graph>::vertex_descriptor vp,vq;
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,Gamma);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,Gamma);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_vhomom = get(vertex_p_vhomom,Gamma);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_vhomom = get(vertex_q_vhomom,Gamma);
    
    property_map<Graph,vertex_name_t>::type
    gvname = get(vertex_name,G);
    
    VertexMap vmap;
    
    for(tie(ei,ei_end)=edges(Gamma); ei!=ei_end; ei++)
    {
        p=p_homom(*ei);
        q=q_homom(*ei);
        np = h[p]; nq = h[q];
        
        put(p_homom,*ei,np);
        put(q_homom,*ei,nq);
    }
    
    for(int i=0;i<num_vertices(G);i++)
    {
        vmap.insert(VertexMap::value_type(gvname(i),i));
    }
    
    for(tie(vi,vi_end)=vertices(Gamma); vi!=vi_end; vi++)
    {
        p=p_vhomom(*vi);
        q=q_vhomom(*vi);
        
        np = h[p]; nq = h[q];
        
        vp = vmap[np]; vq = vmap[nq];
        
        put(p_vhomom,*vi,vp);
        put(q_vhomom,*vi,vq);
    }
    
    return Textile(Gamma,G);
    
}

Textile CreateTranspose(Textile T)
{
    int i;
    GammaGraph A(num_vertices(T.first));
    Graph B;
    transpose_graph(T.first,A);
    transpose_graph(T.second,B);
    
    graph_traits<GammaGraph>::edge_iterator ei,ei_end;
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end,gi,gi_end;
    graph_traits<GammaGraph>::out_edge_iterator aoei,aoei_end;
    graph_traits<Graph>::edge_iterator gei,gei_end;
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,T.first);
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_ahomom = get(edge_p_homom,A);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_ahomom = get(edge_q_homom,A);
    
    property_map<GammaGraph,edge_name_t>::type
    ename = get(edge_name,T.first);
    
    property_map<GammaGraph,vertex_name_t>::type
    gammavname = get(vertex_name,T.first);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_avhom = get(vertex_p_vhomom,A);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_avhom = get(vertex_q_vhomom,A);  
    
    property_map<GammaGraph,vertex_name_t>::type
    avname = get(vertex_name,A);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_gammavhom = get(vertex_p_vhomom,T.first);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_gammavhom = get(vertex_q_vhomom,T.first);
    
    
    property_map<Graph,edge_name_t>::type
    be_name = get(edge_name,B);
    
    cout << "About to add edges" << endl;
/**
    for(tie(ei,ei_end)=edges(T.first);ei!=ei_end;ei++)
    {
        cout << "Inserting an edge from " << target(*ei,T.first) << " to " << source(*ei,T.first) " with p = " 
             << p_homom(*ei) << 
        add_edge(target(*ei,T.first),source(*ei,T.first),PQ_Homoms(p_homom(*ei),Q_Homom(q_homom(*ei),ename(*ei))),A);
    }
    
    for(tie(vi,vi_end)=vertices(A);vi!=vi_end;++vi)
    {
        tie(aoei,aoei_end)=out_edges(*vi,A);
        
        for(tie(gei,gei_end)=edges(B);gei!=gei_end;gei++)
        {
            if(be_name(*gei)==p_ahomom(*aoei))
            {
                cout << "Putting " << source(*gei,B) << " as pvhom for " << *vi << endl;
                put(p_avhom,*vi,source(*gei,B));
                break;
            }
        }
        
        for(tie(gei,gei_end)=edges(B);gei!=gei_end;gei++)
        {
            if(be_name(*gei)==q_ahomom(*aoei))
            {
                cout << "Putting " << source(*gei,B) << " as qvhom for " << *vi << endl;
                put(q_avhom,*vi,source(*gei,B));
                break;
            }
        }
        
        
        
        put(avname,*vi,gammavname(*vi));
    }

    */
    return Textile(A,B);
    
    
}

Textile CreateOneOneTextile(Textile T)
{
    int blocksize;
    
    Textile S = HigherNBlock(T,2);
    
    GammaGraph A = S.first;
    Graph  B = T.second;
    
    GammaOEI oei,oei_end;
    GammaVI vi,vi_end;
    GEI gei,gei_end;
    
    property_map<GammaGraph,edge_p_homom_t>::type
    pa = get(edge_p_homom,A);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    qa = get(edge_q_homom,A);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    qv = get(vertex_q_vhomom,A);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    pv = get(vertex_p_vhomom,A);
    
    property_map<Graph,edge_name_t>::type
    bename = get(edge_name,B);
    
    
    
    
    graph_traits<GammaGraph>::edge_iterator ei,ei_end;
    
    tie(ei,ei_end)=edges(A);
    
    if(pa(*ei).length() % 2 != 0)
    {
        cout << "THE T[2] P LENGTH IS ODD, RETURNING T" << endl;
        return T;
    }
    else
    {
        blocksize = pa(*ei).length()/2;
    }
    
    
    cout << "Block size is "  << blocksize << endl;
    
    
    for(tie(ei,ei_end)=edges(A); ei!=ei_end; ei++)
    {
        string newp=string(pa(*ei),0,blocksize),newq=string(qa(*ei),blocksize,blocksize);  
        
        //    cout << "Old p is " << pa(*ei) << " and new p is " << newp << endl;
        //    cout << "Old q is " << qa(*ei) << " and new q is " << newq << endl;
        
        put(pa,*ei,newp);
        put(qa,*ei,newq);
        
    }
    
    for(tie(vi,vi_end)=vertices(A); vi!=vi_end; vi++)
    {
        tie(oei,oei_end) = out_edges(*vi,A);

        for(tie(gei,gei_end)=edges(B); gei!=gei_end; gei++)
        {
            if(bename(*gei) == pa(*oei))
            {
                put(pv,*vi,source(*gei,B));
                break;
            } // if 
        } // for gei

        for(tie(gei,gei_end)=edges(B); gei!=gei_end; gei++)
        {
            if(bename(*gei) == qa(*oei))
            {
                put(qv,*vi,source(*gei,B));
                break;
            } // if 
        } // for gei       
        
    } // for vi
    
    return Textile(A,B);
    
}


VVec compatibleSet(Textile & T, bool porq, VVec U, string s)
{
    VVec C;
    
    graph_traits<GammaGraph>::out_edge_iterator oei,oei_end;
    
    vector<graph_traits<GammaGraph>::vertex_descriptor>::iterator uit;
    
    property_map<GammaGraph,edge_p_homom_t>::type
    pt = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    qt = get(edge_q_homom,T.first);
    
    typedef std::tuple<graph_traits<GammaGraph>::vertex_descriptor,string> CSTuple;
    
    stack<CSTuple> CSStack;
    bool found;
    
    //cout << "We are computing S+(U," << s << ")" << endl;
    
    for(uit = U.begin(); uit < U.end(); uit++)
    {
        CSStack.push(CSTuple(*uit,s));
    }
    
    while(!CSStack.empty())
    {
        CSTuple current = CSStack.top();
        CSStack.pop();
        
        graph_traits<GammaGraph>::vertex_descriptor v = get<0>(current);
        string currs = get<1>(current);
        
        
        graph_traits<GammaGraph>::out_edge_iterator ni,ni_end;
        
        for(tie(ni,ni_end)=out_edges(v,T.first);ni!=ni_end;ni++)
        {
            int len = (pt(*ni)).length(),lens = currs.length();
            
            //	  cout << "p(edge) = " << pt(*ni) << "and our s is " << currs << endl;
            
            
            if(pt(*ni) == string(currs,0,len) && pt(*ni) != currs)
            {
                //   cout << "We are pushing (" << target(*ni,T.first) <<  "," << string(currs,len,lens) << ")" << endl;
                CSStack.push(CSTuple(target(*ni,T.first),string(currs,len,lens)));
            }
            else if(pt(*ni)==currs)
            {
                found = false;
                for(uit = C.begin(); uit < C.end(); uit++)
                {
                    if(*uit == target(*ni,T.first))
                    {
                        found = true;
                        
                    }
                }
                if(!found)
                {
                    C.push_back(target(*ni,T.first));
                    
                }
                
            }
        }
    }
    
    sort(C.begin(),C.end());
    
    return C;
    
}

/* This function generates the induced right resolving labelled graph of p
 * which will be stored as the p map in a new graph. It does not, however, have a q map
 * defined. One must proceed with caution. Also note that p might not be right resolving as
 * a homomorphism.
 */
Textile inducedRp(Textile T)
{
    GammaVI vi,vi_end;
    GammaOEI oei,oei_end;
    GEI gei,gei_end;
    
    GammaGraph Gamma;
    Graph G=T.second;
    
    property_map<GammaGraph,edge_p_homom_t>::type
    pte = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_p_homom_t>::type
    pe = get(edge_p_homom,Gamma);
    
    property_map<GammaGraph,edge_name_t>::type
    tfename = get(edge_name,T.first);
    
    property_map<GammaGraph,edge_name_t>::type
    gamma_ename = get(edge_name,Gamma);
    
    property_map<GammaGraph,vertex_name_t>::type
    gamma_vname = get(vertex_name,Gamma);
    
    property_map<Graph,edge_name_t>::type
    tsename = get(edge_name,T.second);
    
    property_map<GammaGraph,vertex_name_t>::type
    tfvname = get(vertex_name,T.first);
    
    property_map<Graph,edge_name_t>::type
    gename = get(edge_name,G);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    pvgamma = get(vertex_p_vhomom, Gamma);
    
    stack<VVec> toInv;
    
    set<VVec> seen;
    vector<string> codex;
    
	for(tie(vi,vi_end)=vertices(T.first);vi!=vi_end;vi++)
	{
		VVec v(1,*vi);
		toInv.push(v);
		seen.insert(v);
	}

	for(tie(gei,gei_end)=edges(T.second);gei!=gei_end;gei++)
	{
		string s = tsename(*gei);
		codex.push_back(s);
	}

	while(!toInv.empty())
	{
		vector<string>::iterator sit;
		VVec::iterator it;
		VVec v = toInv.top();
		toInv.pop();

	/*	cout << "Look at {";

		for(it = v.begin(); it < v.end(); it++)
		{
			cout << *it << " ";
		}

		cout << "}" << endl; */

		for(sit=codex.begin(); sit < codex.end(); sit++)
		{
			VVec S = compatibleSet(T,1,v,*sit);
			if(S.size() > 0)
			{
				cout << "S is of size " << S.size() << endl << "It is {";
				for(it = S.begin(); it < S.end(); it++)
				{
					cout << *it << " ";
				}
				cout << "}" << endl; 

				if(seen.find(S)==seen.end())
				{
					seen.insert(S);
					toInv.push(S);
				}
				else{
		//			cout << "We've already seen it, not inserting." << endl;
				}
			}
		}

	}

    
    
    // Before we create the actual graph, we need to trim our seen set. We need to check if we have non-maximal
    // VVecs in there and if so, get rid of them.
    set<VVec>::iterator si,ti;
    set<VVec> toDelete;
    set<VVec>::iterator tdi;
    bool done=false,fdone=false,found=false;
    VertexMap vmap;
    
    for(si=seen.begin(); si!=seen.end(); si++)
    {
        for(ti=seen.begin(); ti!=seen.end(); ti++)
        {
            VVec a=*si,b=*ti;
			if(VVecSubset(a,b) && toDelete.find(a) == toDelete.end())
            {
				
       //         cout << "We are planning to erase";
       //         printVVec(a); cout << endl;
       //         toDelete.insert(a);
            }
        }
    }

	cout << "We have " << toDelete.size() << " vectors to delete" << endl;
    
    for(tdi=toDelete.begin(); tdi!=toDelete.end(); tdi++)
    {
         //   cout << "We are actually erasing "; printVVec(*tdi); cout << endl;
        //seen.erase(*tdi);
    }
    
    
     cout << "Our final Seen consists of:";
    
    for(si=seen.begin(); si != seen.end(); si++)
    {
        
        graph_traits<GammaGraph>::vertex_descriptor v = add_vertex(Gamma);
        VVec::iterator it;
        string vname;
        VVec curr=*si;
        
               printVVec(curr);
        vname = ssVVec(curr);
        put(gamma_vname,v,vname);
        
        //We need to create a VMap so we can add edges later
        vmap.insert(VertexMap::value_type(vname,v));
    }
    cout << endl;
    
    for(si=seen.begin(); si != seen.end(); si++)
    {
        string s=ssVVec(*si),t;
        graph_traits<GammaGraph>::vertex_descriptor v=vmap[s],w;
        vector<string>::iterator cit;
        VVec S;
        
        for(cit = codex.begin(); cit < codex.end(); cit++)
        {
            S = compatibleSet(T,1,(*si),(*cit));
            if(S.size() > 0) {
                t = ssVVec(S);
                w = vmap[t];
                add_edge(v,w,PQ_Homoms((*cit),Q_Homom(string(""),string("(" + s + "," + (*cit) + ")"))),Gamma);
            }  
        }
        
        
    }
    
    /*  while(!toInv.empty())
     {
     VVec::iterator it;
     VVec v = toInv.top();
     for(it = v.begin(); it < v.end(); it++)
     {
     cout << *it;
     }
     cout << endl;
     toInv.pop();
     }
     */
    
    // We now need to fill in the vertex homomorphisms
    for(tie(vi,vi_end)=vertices(Gamma); vi!=vi_end; vi++)
    {
        string currPE;
        tie(oei,oei_end)=out_edges(*vi,Gamma);
        
        currPE = pe(*oei);
        tie(gei,gei_end)=edges(G);
        while(!found)
        {
            if(currPE == gename(*gei))
            {
                put(pvgamma,*vi,source(*gei,G));
                found = true;
            }
            gei++;
        } // while not found
        found = false;
    } // for vi
    
    return Textile(Gamma,G);
}

Textile inducedRq(Textile T)
{
    Textile Tinv = CreateInverse(T);
    Textile iirp = NewInducedRp(Tinv);
    cout << "IIRP created" << endl;
    PrintFullTextileInfo(iirp);
    return CreateInverse(iirp);
    
}

Textile inducedLp(Textile T)
{
    Textile Ttrans = CreateTranspose(T);
    cout << "Created transpose" << endl;
    Textile tirq = NewInducedRp(Ttrans);
    cout << "TIRQ created" << endl;
    PrintFullTextileInfo(tirq);
    return CreateTranspose(tirq);
}

Textile inducedLq(Textile T)
{
    Textile Tinv = CreateInverse(T);
    cout << "Created inverse" << endl;
    Textile iilq = inducedLp(Tinv);
    cout << "IILQ created" << endl;
    return CreateInverse(iilq);
}

// Checks if a Textile is 1-1 by looking at its induced left and right homomorphisms and checking if they are definite.
bool is1to1(Textile T)
{
    Textile irq=inducedRq(T),ilp=inducedLp(T),irp=NewInducedRp(T),ilq=inducedLq(T);
    
        bool rq,lp,rp,lq;
        cout << "In 1-1 function for textile " << endl;
        SmartPrintTextileInfo(T);
        
        cout << "Printing IRQ" << endl;
        SmartPrintTextileInfo(irq);
        rq = IsqRightDefinite(irq);
        
        cout << "Printing ILQ" << endl;
        SmartPrintTextileInfo(ilq);
        lq = IsqLeftDefinite(ilq);
        
        cout << "Printing IRP" << endl;
        SmartPrintTextileInfo(irp);
        rp = IspRightDefinite(irp);
        
        cout << "Printing ILP" << endl;
        SmartPrintTextileInfo(ilp);
        lp = IspLeftDefinite(ilp);
        
        return rq && lq && rp && lp;
        
  }

// Checks if a Textile is one sided 1-1 by looking at its induced left and right homomorphisms and checking if they are definite.
bool isOneSided1to1(Textile T)
{
    Textile ilp=inducedLp(T),irp=NewInducedRp(T);
    

      	bool lp,rp;
        cout << "In One sided 1-1 function for textile " << endl;
        SmartPrintTextileInfo(T);
        
        cout << "Printing IRP" << endl;
        SmartPrintTextileInfo(irp);
        rp = IspRightDefinite(irp);
        
        cout << "Printing ILP" << endl;
        SmartPrintTextileInfo(ilp);
        lp = IspLeftDefinite(ilp);
        
        return rp && lp;
        
  }


// Checks to see if one VVec (vector of vertex descriptors) is a subset of another
bool VVecSubset(VVec& a, VVec& b)
{
    VVec::iterator ait,bit;
    for(ait = a.begin(); ait < a.end(); ait++)
    {
        bool found=false;
        for(bit = b.begin(); bit < b.end(); bit++)
        {
            if(*ait == *bit)
            {
                found = true;
            }
        }
        if(!found)
        {
            // printVVec(a);
            //	  cout << " is NOT a subset of ";
            //	  printVVec(b);	cout << endl;
            return false;
        }
    }
    
    if(a.size() == b.size())
    {
        return false;
    }
    //  printVVec(a); cout << "is a subset of "; printVVec(b); cout << endl;
    return true;
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

/*
 * The analyzer function takes in a Textile system and tries to find LR components of the input textile system.
 * It will also try to find other TMS's in the system and then look for their LR components
 */

void Analyzer(Textile T)
{
    // Our first goal is to check to see if the phi associated with T is an LR map, we will use this as a starting 
    // point to find the LR cone of sigma.
	double maxTopAngle,minTopAngle,maxBotAngle,minBotAngle;
	bool topAngleKnown = false,botAngleKnown = false;

	// Set the currTextile to be T initially
	Textile currTextile = T;
	// T corresponds to looking at the (0,1) direction of the textile system
	double x = 0.0, y = 1.0;

	while(!topAngleKnown && !botAngleKnown)
	{
		if(IsLR(currTextile)) 
		{
			minTopAngle = atan(y/x);
			if(false)  //IsDefinite(currTextile))
			{
				topAngleKnown = true;
				maxTopAngle = atan(y/x);
			} 
		}
	}
}


/*
 * This is a new, better version of the checking isom languages program. It minimizes the automata given by
 * p and q graphs and then checks to see if they are the same.
 */ 
bool NewIsomLanguages(Textile T)
{
    cout << "Checking to see that T has isom languages" << endl;
    SmartPrintTextileInfo(T);
    // We first construct the induced right resolving representations of the two maps. This essentially
    // makes them DFAs that we can minimize then test their equality.
    Textile irp = NewInducedRp(T), irq = CreateInverse(inducedRq(T));
    
    cout << "Printing the induced p and q graphs" << endl;
  //  PrintFullTextileInfo(irp);
  //  PrintFullTextileInfo(irq);
    
    if(!IspRightResolving(irp) || !IspRightResolving(irq))
    {
        return false;
    }
    
    cout << "Creating minimized DFAs" << endl;
    
    Textile mirp = DFAMinimization(irp), mirq = DFAMinimization(irq);
    
    
    int psize = num_vertices(mirp.first), qsize = num_vertices(mirq.first),i,N;
    bool dble,found;
    GammaOEI oei,oei_end,ofi,ofi_end;
    
    Textile DS = DirectSum(mirp,mirq);
    
    set<pair<VD,VD> > seen;
    
  //  PrintFullTextileInfo(DS);
    

    
    N = psize + qsize - 2;
    
    property_map<GammaGraph,edge_p_homom_t>::type
    dsp_homom = get(edge_p_homom,DS.first);

    property_map<GammaGraph,vertex_name_t>::type
    dsp_vname = get(vertex_name,DS.first);

    property_map<GammaGraph,vertex_index_t>::type
    dsp_index = get(vertex_index,DS.first);

  //  ofstream gviz("NI.gv");
    
  //  write_graphviz(gviz, DS.first, make_label_writer(dsp_index),make_label_writer(dsp_homom));

    
    // We fix a start vertex in the mirp section of our direct sum, then iterate through the
    // mirq section looking for a vertex that is equivalent to our start vertex.
    for(i=psize; i < psize+qsize; i++)
    {
        dble = false;
        
        // Here we check to see if states 0 and i are equivalent
        cout << "Checking if " << 0 << " and " << i << " are equivalent. We are at "  << i-psize << " of " << qsize  << '\r';
        cout.flush();
        // Iterate through words of length k in a breadth first search
        
        priority_queue<PQOEIElement,vector<PQOEIElement>,pqcomparison> OEIPQ;
        
        tie(oei,oei_end)=out_edges(0,DS.first);   
        
    //       cout << "Pushing " << 0 << ", i = " << 0 << ", j = " << i << ", phomom(oei) = " << dsp_homom(*oei) << endl; 
        OEIPQ.push(PQOEIElement(0,0,i));
        
        while(!OEIPQ.empty() && !dble)
        {
            PQOEIElement current = OEIPQ.top();
            OEIPQ.pop();
    //        cout << "Popped: k = " << get<0>(current) << ", i = " << get<1>(current) << ", j = " << get<2>(current) << " stacksize is " << OEIPQ.size() << endl;
            
            // Make sure both vertices have the same out degree.
            if(out_degree(get<1>(current),DS.first) != out_degree(get<2>(current),DS.first))
            {
                           cout << "Different out degrees" << endl;
                dble = true;
            } // if
            else {
                for(tie(oei,oei_end)=out_edges(get<1>(current),DS.first); oei != oei_end; oei++)
                {
                    found = false;
                    for(tie(ofi,ofi_end)=out_edges(get<2>(current),DS.first); ofi != ofi_end; ofi++)
                    {
                        if(dsp_homom(*ofi) == dsp_homom(*oei))
                        {
                            found = true;
       //                                              cout << "Found a match for " << dsp_homom(*oei) <<endl;
                            if(get<0>(current) < N && target(*oei,DS.first) != target(*ofi,DS.first)  && target(*oei,DS.first) != source(*oei,DS.first) 
                               && seen.find(pair<int,int>(target(*oei,DS.first),target(*ofi,DS.first))) == seen.end()) { 
        //                                                      cout << "Pushing: k = " << get<0>(current)+1 << ", i = " << target(*oei,DS.first)
        //                                                      << ", j = " << target(*ofi,DS.first) << endl;
                                seen.insert(pair<int,int>(target(*oei,DS.first),target(*ofi,DS.first)));
                                OEIPQ.push(PQOEIElement(get<0>(current)+1,target(*oei,DS.first),target(*ofi,DS.first)));
                            }
                            break;
                        } // if p = p
                        
                    } // for ofi
                    if(!found) {
                                         cout << "Did not find a match. They are distinguisable" << endl;
                        dble = true;
                    }// if !found
                } // for oei
                // Now we have to find a match for i and j
            } // else
            
        } // while
        if(!dble)
        {
            cout << "We have found that " << i << " and 0 are indistinguishable" << endl;
            return true;
        }
        
    } // for i
    
    return false;
    
}

/*
 *
 */
Textile DirectSum(Textile S,Textile T)
{
    int sSize = num_vertices(S.first), tSize = num_vertices(T.first);
    
    GammaVI vi,vi_end;
    GammaOEI oei,oei_end;
    
    GammaGraph Sum(sSize + tSize);
    
    property_map<GammaGraph,edge_p_homom_t>::type
    sp_homom = get(edge_p_homom,S.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    sq_homom = get(edge_q_homom,S.first);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    sp_vhomom = get(vertex_p_vhomom,S.first);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    sq_vhomom = get(vertex_q_vhomom,S.first);    
    
    property_map<GammaGraph,edge_name_t>::type
    se_name = get(edge_name,S.first);
    
    property_map<GammaGraph,vertex_name_t>::type
    s_vname = get(vertex_name,S.first);
    
    property_map<GammaGraph,edge_p_homom_t>::type
    tp_homom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    tq_homom = get(edge_q_homom,T.first);
    
    property_map<GammaGraph,edge_name_t>::type
    te_name = get(edge_name,T.first);    
    
    property_map<GammaGraph,vertex_name_t>::type
    t_vname = get(vertex_name,T.first);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    tp_vhomom = get(vertex_p_vhomom,T.first);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    tq_vhomom = get(vertex_q_vhomom,T.first);   
    
    property_map<GammaGraph,vertex_name_t>::type
    sum_vname = get(vertex_name,Sum);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    sump_vhomom = get(vertex_p_vhomom,Sum);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    sumq_vhomom = get(vertex_q_vhomom,Sum);   
    
    for(tie(vi,vi_end)=vertices(S.first); vi != vi_end; vi++)
    {
        put(sum_vname,*vi,s_vname(*vi));
        put(sump_vhomom,*vi,sp_vhomom(*vi));
        put(sumq_vhomom,*vi,sq_vhomom(*vi));
        
        for(tie(oei,oei_end)=out_edges(*vi,S.first); oei != oei_end; oei++)
        {
            add_edge(*vi,target(*oei,S.first),PQ_Homoms(sp_homom(*oei),Q_Homom(sq_homom(*oei),se_name(*oei))),Sum);
        } // for oei
    } // for vi
    
    for(tie(vi,vi_end)=vertices(T.first); vi != vi_end; vi++)
    {
        put(sum_vname,(*vi)+sSize,t_vname(*vi));
        put(sump_vhomom,(*vi)+sSize,tp_vhomom(*vi));
        put(sumq_vhomom,(*vi)+sSize,tq_vhomom(*vi));
        
        for(tie(oei,oei_end)=out_edges(*vi,T.first); oei != oei_end; oei++)
        {
            add_edge((*vi)+sSize,target(*oei,T.first)+sSize,PQ_Homoms(tp_homom(*oei),Q_Homom(tq_homom(*oei),te_name(*oei))),Sum);
        } // for oei
    } // for vi    
    
    
    
    return Textile(Sum,S.second);
}

/*
 *  This function actual performs the DFA minimization on the p map, which is assumed to be right resolving.
 */
Textile DFAMinimization(Textile T)
{
    // N is the size of the longest word we need to check when looking at distinguishing states
    int N,i,j,numStates=num_vertices(T.first),blocksize;
    GammaEI ei,ei_end;
    GammaOEI oei,oei_end,ofi,ofi_end;
    bool dble,found,seenlast=false,unnec;
    vector<VD>::iterator it;
    vector< vector<VD> >::iterator vit;
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,T.first);
    
    
    if(!IspRightResolving(T)) {
        cout << "The given p map is not right resolving we have returned the original textile" << endl;
        SmartPrintTextileInfo(T);
        return T;
    }
    
    if(numStates-2 < 1) {
        N = 1;
    }
    else {
        N = numStates - 2;
    }
    
    tie(ei,ei_end) = edges(T.first);
    blocksize = p_homom(*ei).size();
    
    vector< vector<VD> > stateEquiv;
    
    
    for(i=0;i<numStates-1;i++) {
        unnec = false;
        for(vit=stateEquiv.begin(); vit < stateEquiv.end(); vit++)
        {
            for(it=(*vit).begin(); it < (*vit).end(); it++)
            {
                if(*it == i)
                {
         //           cout << *it << " is unneccesary, skipping" << endl;
                    unnec = true;
                }
            }
        }
        if(!unnec) {
            vector<VD> EquivForI;  
            EquivForI.push_back(i);
            
            for(j=i+1;j<numStates;j++) {
                dble = false;
                
                // Here we check to see if states i and j are equivalent
                //    cout << "Checking if " << i << " and " << j << " are equivalent." << endl;
                // Iterate through words of length k in a breadth first search
                
                priority_queue<PQOEIElement,vector<PQOEIElement>,pqcomparison> OEIPQ;
                
                tie(oei,oei_end)=out_edges(i,T.first);           
                //   tie(ofi,ofi_end)=out_edges(j,T.first);
                
                //   cout << "Pushing " << 0 << ", i = " << i << ", j = " << j << ", phomom(oei) = " << p_homom(*oei) << endl; 
                OEIPQ.push(PQOEIElement(0,i,j));
                
                while(!OEIPQ.empty() && !dble)
                {
                    PQOEIElement current = OEIPQ.top();
                    OEIPQ.pop();
                    //        cout << "Popped: k = " << get<0>(current) << ", i = " << get<1>(current) << ", j = " << get<2>(current) << endl;
                    
                    // Make sure both vertices have the same out degree.
                    if(out_degree(get<1>(current),T.first) != out_degree(get<2>(current),T.first))
                    {
                        //             cout << "Different out degrees" << endl;
                        dble = true;
                    } // if
                    else {
                        for(tie(oei,oei_end)=out_edges(get<1>(current),T.first); oei != oei_end; oei++)
                        {
                            found = false;
                            for(tie(ofi,ofi_end)=out_edges(get<2>(current),T.first); ofi != ofi_end; ofi++)
                            {
                                if(p_homom(*ofi) == p_homom(*oei))
                                {
                                    found = true;
                                    //                         cout << "Found a match for " << p_homom(*oei) <<endl;
                                    if(get<0>(current) < N && target(*oei,T.first) != target(*ofi,T.first)) { 
                                        //                             cout << "Pushing: k = " << get<0>(current)+1 << ", i = " << target(*oei,T.first)
                                        //                             << ", j = " << target(*ofi,T.first) << endl;
                                        OEIPQ.push(PQOEIElement(get<0>(current)+1,target(*oei,T.first),target(*ofi,T.first)));
                                    }
                                    break;
                                } // if p = p
                                
                            } // for ofi
                            if(!found) {
                                //                     cout << "Did not find a match. They are distinguisable" << endl;
                                dble = true;
                            }// if !found
                        } // for oei
                        // Now we have to find a match for i and j
                    } // else
                    
                } // while
                
                if(!dble)
                {
                    //             cout << "We have created an entry in stateEquiv for " << i << ", " << j << endl;
                    EquivForI.push_back(j);
                    if(j==numStates-1)
                    {
                        seenlast=true;
                    }
                }
                
            } // for j
            
     //       cout << "State equivalences for " << i << " are : ";
            for(it=EquivForI.begin(); it < EquivForI.end(); it++)
            {
      //          cout << *it << ", ";
            }
            
            stateEquiv.push_back(EquivForI);
            
            // Now we need to search through the previous entries
            
            
    //        cout << endl;
        } // if unnec
    } // for i
    
    if(!seenlast)
    {
        stateEquiv.push_back(vector<VD>(1,numStates-1));
    }
    
    for(vit=stateEquiv.begin(); vit < stateEquiv.end(); vit++)
    {
   //     cout << "State equiv contains: ";
        for(it = (*vit).begin(); it < (*vit).end(); it++)
        {
  //          cout << *it << ", ";
        }
 //       cout << endl;
    }
    
    return Quotient(T,stateEquiv);
}

Textile Composition(Textile S, Textile T)
{
    return CreateDual(ProductTextile(CreateDual(S),CreateDual(T)));
}

// Looks for a conjugacy to a given textile system. MAX indicates how large of higher block
// shifts to make. Making MAX larger than 4 or 5 will take a really long time.
Textile LookForConjugacy(Textile T, int MAX, bool* found, int start)
{
    bool TND,TOO,NND,NOO,DOO,DND;
    int i;
    // We first check that T is nondegenerate
    cout << "We are looking for conjugacy and so are going to print T" << endl;
    SmartPrintTextileInfo(T);
    
    cout << "Checking to see if T is nondegenerate" << endl;
    TND = NewIsomLanguages(T);
    if(TND)
    {
        cout << "T is nondegen, checking to see if it is 1-1" << endl;
        TOO = is1to1(T);
        
        
        
        if(TOO)
        {
            Textile TD = CreateDual(T);
            cout << " T is 1-1, created T* and now checking if it's nondegen after printing it" << endl;
            SmartPrintTextileInfo(TD);
            DND = NewIsomLanguages(TD);
            
            if(DND)
            {
                if(is1to1(TD))
                {
                    cout << "TD is good to go already" << endl;
                    return T;
                }
            }
            else{
                cout << "TD is not nondegenerate" << endl;
                for(i=start;i<MAX;i++)
                {
                    cout << "We are starting to look at trim((T*)^[" << i << "])*" << endl;
                    Textile TiD = AutoRenamer(AutoHomomLite(Trim(AutoHomomLite(HigherNBlock(TD,i)))));
                    cout << "TiD has been created, starting to check if nondegen" << endl;
                    NND = NewIsomLanguages(TiD);
                    if(NND)
                    {
                        cout << "TiD is nondegenerate, checking if it's 1-1" << endl;
                        NOO = is1to1(TiD);
                        if(NOO)
                        {
							*found = true;
                            cout << endl << endl << "FOUND CONJUGACY!!!!!!" << endl << endl;
                            cout << "TiD is also 1-1 at level " << i << endl;
			    //                cout << "P-F Eigenvalue is " << PFEigenvalue(TiD.second) << endl;
                            return TiD;
                        }
                        else
                        {
                            cout << "TiD is not 1-1" << endl;
                        }
                    }
                } // for
            }
        } // TOO
        
    } // TND
    cout << "Could not find a conjugacy, return the original" << endl;
	*found = false;
    return T;
    
}

// Recursively creates the composition power of a textile system
Textile CompositionPower(Textile T, int n)
{
    cout << "Composition power " << n << endl;
    if(n==2)
    {
        return Composition(T,T);
    }
    else
    {
        return Composition(T,CompositionPower(T,n-1));
    }
}

// This function creates the textile system who phi is (n,m).
// It changes the basis to use the (1,1) textile and the original system.
Textile CreateNMTextile(Textile T, int n, int m)
{
    int newn = n, newm = -n + m;
    
    cout << "New N is " << n << " and new M is " << newm << endl;
    
    Textile Toneone = CreateOneOneTextile(T);
  
 //   cout << "Printing 1,1 textile" << endl;
    
 //   PrintFullTextileInfo(Toneone);
    
    
    Textile N,M;
    
    if(newn<0 && newn!=-1)
    {
        N = CompositionPower(CreateInverse(Toneone),-newn);
        cout << "N has been created, it is: " << endl;
        SmartPrintTextileInfo(N);

    }
    else if(newn == 1){
        N = Toneone;
        cout << "N has been created, it is: " << endl;
        SmartPrintTextileInfo(N);

    }
    else if(newn == -1)
    {
        N = CreateInverse(Toneone);
        cout << "N has been created, it is: " << endl;
        SmartPrintTextileInfo(N);

    }
    else{
        N = CompositionPower(Toneone,newn);
        cout << "N has been created, it is: " << endl;
        SmartPrintTextileInfo(N);

    }
    
    if(newm<0 && newm!=-1)
    {
        M = CompositionPower(CreateInverse(T),-newm);
        cout << "Now, M has been created, it is: " << endl;
        SmartPrintTextileInfo(M);

    }
    else if(newm == 1){
        M = T;
        cout << "Now, M has been created, it is: " << endl;
        SmartPrintTextileInfo(M);

    }
    else if(newm == -1)
    {
        M = CreateInverse(T);
        cout << "Now, M has been created, it is: " << endl;
        SmartPrintTextileInfo(M);

    }
    else{
        M = CompositionPower(T,newm);
        cout << "Now, M has been created, it is: " << endl;
        SmartPrintTextileInfo(M);

    }

//	cout << "Using old , we get" << endl;
	
//	PrintFullTextileInfo(Trim(Composition(N,M)));
    
    return Trim(AutoHomomLite(Composition(N,M)));
    
}

// This function will go through a textile, look for identical rows, then add all of those rows to the same equiv class.
// The quotient function will check to make sure that this is actually a state equiv and if its not, it will just
// return the original textile.
Textile AutoHomom(Textile T)
{
    // We need to check to see if we've already seen a vertex
    vector<VD> seen;
    vector<VD>::iterator si;
    
    vector<vector<VD> > classes;
    
    vector<multiset<string> * > classData;
    vector<multiset<string> * >::iterator ci;
    
    GammaVI vi,vi_end;
    GammaOEI oei,oei_end;
    
    bool done,found;
    int i;
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_hom = get(edge_p_homom,T.first);
    
    property_map<GammaGraph,edge_q_homom_t>::type
    q_hom = get(edge_q_homom,T.first);
    
	cout << "We are going to AutoHomom the following textile:" << endl;
	SmartPrintTextileInfo(T);


    for(tie(vi,vi_end)=vertices(T.first); vi!=vi_end; vi++)
    {
        found = false;
        done = false;
        for(si = seen.begin(); si < seen.end(); si++)
        {
            if(*si == *vi)
                done = true;
        }
        if(!done)
        {
            multiset<string> * currData = new multiset<string>;
            for(tie(oei,oei_end)=out_edges(*vi,T.first); oei!=oei_end; oei++)
            {
				cout << "Checking for the edge "<< source(*oei,T.first) << "," << target(*oei,T.first) << endl;
                currData->insert(p_hom(*oei) + "/" + q_hom(*oei));
            }
            for(ci = classData.begin(), i=0; ci < classData.end() && !found; ci++, i++)
            {
                if(*(*ci)==(*currData))
                {
                    (classes[i]).push_back(*vi);
                    found = true;
                } // if
            } // for ci
            if(!found)
            {
                classData.push_back(currData);
                classes.push_back(vector<VD>(1,*vi));
            }
        }
    }
    
    return Quotient(T,classes);
}

// This function is for if Autohomom doesn't work. It is a little more demanding than autohomom.
// It cares about the target edge and the p/q homoms, as opposed to just the p/q homoms.
Textile AutoHomomLite(Textile T)
{
	property_map<GammaGraph,edge_p_homom_t>::type
		p_hom = get(edge_p_homom,T.first);

	property_map<GammaGraph,edge_q_homom_t>::type
		q_hom = get(edge_q_homom,T.first);
		
	property_map<GammaGraph,edge_name_t>::type
		ename = get(edge_name,T.first);

	property_map<GammaGraph,vertex_name_t>::type
		vname = get(vertex_name,T.first);

	GammaVI vi,vi_end,wi,wi_end;
	GammaOEI oei,oei_end,ofi,ofi_end;
	
	vector<vector<VD> > classes;
	vector<vector<VD> >::iterator ci;
	
	vector<VD>::iterator si;
	
	bool found,seen;

	cout << "We are going to AutoHomomLite(tm) the following textile:" << endl;
	SmartPrintTextileInfo(T);

	for(tie(vi,vi_end)=vertices(T.first);vi!=vi_end;vi++)
	{
		seen=false;
		for(ci=classes.begin(); ci!= classes.end(); ci++)
		{
			for(si=(*ci).begin(); si!=(*ci).end() && !seen; si++)
			{
				if(*si==*vi)
				{
					seen=true;
					cout << "We have seen " << *si << endl;
				}
			}
		} 

		if(!seen)
		{
			VD v = *vi;
			vector<VD> vivec = vector<VD>(1,v); 

			for(tie(wi,wi_end)=vertices(T.first),wi=vi,wi++;wi!=wi_end;wi++)
			{
	//			cout << "Checking vertex " << vname(*wi) << " for compatibility" << endl;
				if(out_degree(*vi,T.first)==out_degree(*wi,T.first))
				{
					for(tie(oei,oei_end)=out_edges(*vi,T.first);oei!=oei_end;oei++)
					{
		//				cout << "Looking at the edge from " << *vi << " to " << target(*oei,T.first) <<  " with p q " << p_hom(*oei) << " / " << q_hom(*oei) << endl;
						tie(ofi,ofi_end)=out_edges(*wi,T.first);
						found = false;
						while(!found && ofi!=ofi_end)
						{
				//			cout << "Comparing to edge from " << *wi << " to " << target(*ofi,T.first) <<  " with p q " << p_hom(*ofi) << " / " << q_hom(*ofi) << endl;
							if(target(*ofi,T.first)==target(*oei,T.first) && p_hom(*ofi)==p_hom(*oei) && q_hom(*ofi)==q_hom(*oei))
							{
								found = true;
							}
							else{
				//				cout << "Didn't find it" << endl;
								ofi++;
							}
						}
						if(ofi==ofi_end)
						{
				//			cout << "No match" << endl;
							break;
						}
					}
					if(oei==oei_end) // We made it all the way through the oei iterator, meaning all the out_edges found a match
					{
						VD w=*wi;
						vivec.push_back(w);
					}

				}
				else{
			//		cout << "Different out degrees" << endl;
				}

			} // wi

			classes.push_back(vivec);
		}
	} // vi 


	return Quotient(T,classes);
}

// Useful when making compositions and powers of textiles, this will shorten all
// of the names of the edges and p,q homoms. 
Textile AutoRenamer(Textile T)
{
    GammaEI ei,ei_end;
    GEI gei,gei_end;
    GammaVI vi,vi_end;

    int numGedges = num_edges(T.second), len,i,vlen;
    
    len = (int) ceil(log(numGedges)/log(26));
    vlen = (int) ceil(log(num_vertices(T.first))/log(26));
    
//    cout << "Length is " << len << " while numGedges is " << numGedges << endl;
    
    property_map<Graph,edge_name_t>::type
    gename = get(edge_name,T.second);
    
    property_map<GammaGraph,vertex_name_t>::type
    vname = get(vertex_name,T.first);
    
    map<string,string> GedgeMap;
    
    for(tie(vi,vi_end)=vertices(T.first), i=0; vi!=vi_end; vi++, i++)
    {
        put(vname,i,Namer(i,vlen,97));
    }
    
    // We first create the string map to feed into RenameTextile
    for(tie(gei,gei_end)=edges(T.second), i=0; gei!=gei_end; gei++, i++)
    {
        GedgeMap[gename(*gei)] = Namer(i,len);
    }
    
    return RenameTextile(T,GedgeMap);
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

// This function finds the Perron Frobenius Eigenvalue of the input matrix. We use JAMA for this as
// I would definitely screw it up.
/**double PFEigenvalue(Graph G)
{
    // We use the graph G to create an adjacency matrix, then we make a JAMA Eigen class for it
    int N = num_vertices(G),i;
    Array2D<double> Adj(N,N,0.0);
    double max = 0;
    GEI gei,gei_end;
    
    for(tie(gei,gei_end)=edges(G); gei!=gei_end; gei++)
    {
        Adj[source(*gei,G)][target(*gei,G)]+=1.0;
    }
    
    Eigenvalue<double> E(Adj);
    
    Array1D<double> Eig; 
    E.getRealEigenvalues(Eig);
    
    for(i=0; i<N; i++)
    {
        cout << Eig[i] << endl;
        if(Eig[i] > max)
        {
            max = Eig[i];
        }
    }
    
    return max;
    
}
*/

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
/*
Textile RandomTrim(Textile T)
{
    Textile Trimmed=T;
    bool isEndIT=true,isEndTI=true,isEndPQ=true,isEndQP=true;
    unsigned int i,j=0,delcount=0,N,rint;
    //  vector<graph_traits<GammaGraph>::edge_descriptor> toDelete;
    vector<int> vertToDelete;
    bool done,stable,needToCheckVhoms=false,found=false,rdone,rmax;
    time_t start,end;
    double dif,average=0;
    
    property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,Trimmed.first);
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,Trimmed.first);
    property_map<GammaGraph,edge_name_t>::type
    ename = get(edge_name,Trimmed.first);
    property_map<GammaGraph,vertex_name_t>::type
    vname = get(vertex_name,Trimmed.first);
    
    property_map<GammaGraph,vertex_p_vhomom_t>::type
    p_vhom = get(vertex_p_vhomom,Trimmed.first);
    
    property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_vhom = get(vertex_q_vhomom,Trimmed.first);
    
    property_map<Graph,edge_name_t>::type
    gename = get(edge_name,Trimmed.second);
    
    graph_traits<GammaGraph>::edge_iterator ei, ei_end, fi, fi_end;
    graph_traits<GammaGraph>::vertex_iterator vi,vi_end;
    graph_traits<GammaGraph>::out_edge_iterator oei,oei_end;
    graph_traits<Graph>::edge_iterator gei,gei_end;
    //  graph_traits<GammaGraph>::edge_descriptor e;
    GVI gvi,gvi_end;
	
    stack<VD> toCheck;  

	cout << "Our starting graph has " << num_edges(Trimmed.first) << " edges and " << num_vertices(Trimmed.first) << " vertices." << endl; 

    time(&start);
    
	// We will randomly check all of the edges (might change this to half). 
	rmax = num_edges(T.first);

	do{
        j++;
       if(j % 100 == 0)
       {
           time(&end);
           dif = difftime(end,start);
           time(&start);
           N = j/100;
           average = (((N-1)*average)+dif)/N;
           cout << "On round " << j << " deleted at least " << delcount << " edges so far. The last round took " 
                << dif << " seconds, and the average time is " << average << " seconds." <<  '\r';
           cout.flush();
       } 
        stable = true;
        done = false;
		if(rint >= rmax){ 
			rdone = true;
        }
		else{
			rdone = false;
		}

        // Check for sinks 
        for(tie(vi,vi_end)=vertices(Trimmed.first); vi!=vi_end && !done; vi++)
        {
            //	cout << "We are checking " << vname(*vi) << endl;
            if(out_degree(*vi,Trimmed.first)==0 || in_degree(*vi,Trimmed.first)==0)
            {
                //   cout << "The vertex " << vname(*vi) << "is a sink or a source" << endl;
                clear_vertex(*vi,Trimmed.first);
                remove_vertex(*vi,Trimmed.first);
                stable = false;
                done=true;
				rdone=true;
                break;
            }
            //	cout << "It's not." << endl;
        }

        if(!rdone)
        {
			rint++;
            //	cout << "Starting ends loop" << endl;
            for(tie(ei,ei_end)=edges(Trimmed.first); ei != ei_end && !rdone; ++ei) {
				random::uniform_int_distribution<> dist(1, num_edges(Trimmed.first));
				ei = ei+dist(gen);
                //  cout << "Checking " << ename(*ei) << " for ends" << endl;
                for(tie(fi,fi_end)=edges(Trimmed.first); fi != fi_end; ++fi) {
                    
                    if(source(*ei,Trimmed.first)==target(*fi,Trimmed.first)) {
                        isEndIT=false;
                    }
                    
                    if(target(*ei,Trimmed.first)==source(*fi,Trimmed.first)) {
                        isEndTI=false;
                    }
                    
                    if(p_homom(*ei)==q_homom(*fi)) {
                        isEndPQ=false;
                    }
                    
                    if(q_homom(*ei)==p_homom(*fi)) {
                        isEndQP=false;
                    }
                }
                if(isEndIT || isEndTI || isEndPQ || isEndQP) {
                    VD s=source(*ei,Trimmed.first),t=target(*ei,Trimmed.first);
                    toCheck.push(s);
                    toCheck.push(t);
                  //            cout << "Removed Edge (" << s << ", "
                  //            << t << ")" << " named " << ename(*ei) << endl;
                    remove_edge(*ei,Trimmed.first);
                    delcount++;
                    rdone = true;
					done = true;
                    stable = false;
                    // the goal of this inner while loop is to check to see if we've just created a dead end for the surrounding vertices.
                    while(!toCheck.empty())
                    {
                        VD currv = toCheck.top();
                        toCheck.pop();
                        
                        if(out_degree(currv,Trimmed.first)==0 ^ in_degree(currv,Trimmed.first)==0 )
                        {
                            //                  cout << "Clearing " << currv << endl;
                            clear_vertex(currv,Trimmed.first);
                        }
                        
                    }
                }
                isEndIT=true;
                isEndTI=true;
                isEndPQ=true;
                isEndQP=true;
            }
            //	cout << "Done with the ends, checking for empty vertices" << endl;
            if(!done)
            {
                for(tie(vi,vi_end)=vertices(Trimmed.first); vi!=vi_end && !done; vi++)
                {
                    if(out_degree(*vi,Trimmed.first)==0 && in_degree(*vi,Trimmed.first)==0)
                    {
                        remove_vertex(*vi,Trimmed.first);
                        done = true;
                        stable = false;
                    }
                } // for vi
            } // if not done
        } // if not rdone

       if(rdone) // if we're done with randomly checking, we want to do more passes to double check
        {
            //	cout << "Starting ends loop" << endl;
            for(tie(ei,ei_end)=edges(Trimmed.first); ei != ei_end && !done; ++ei) {
                //  cout << "Checking " << ename(*ei) << " for ends" << endl;
                for(tie(fi,fi_end)=edges(Trimmed.first); fi != fi_end; ++fi) {

                    if(source(*ei,Trimmed.first)==target(*fi,Trimmed.first)) {
                        isEndIT=false;
                    }

                    if(target(*ei,Trimmed.first)==source(*fi,Trimmed.first)) {
                        isEndTI=false;
                    }

                    if(p_homom(*ei)==q_homom(*fi)) {
                        isEndPQ=false;
                    }

                    if(q_homom(*ei)==p_homom(*fi)) {
                        isEndQP=false;
                    }
                }
                if(isEndIT || isEndTI || isEndPQ || isEndQP) {
                    VD s=source(*ei,Trimmed.first),t=target(*ei,Trimmed.first);
                    toCheck.push(s);
                    toCheck.push(t);
                  //            cout << "Removed Edge (" << s << ", "
                  //            << t << ")" << " named " << ename(*ei) << endl;
                    remove_edge(*ei,Trimmed.first);
                    delcount++;
                    done = true;
                    stable = false;
                    // the goal of this inner while loop is to check to see if we've just created a dead end for the surrounding vertices.
                    while(!toCheck.empty())
                    {
                        VD currv = toCheck.top();
                        toCheck.pop();

                        if(out_degree(currv,Trimmed.first)==0 ^ in_degree(currv,Trimmed.first)==0 )
                        {
                            //                  cout << "Clearing " << currv << endl;
                            clear_vertex(currv,Trimmed.first);
                        }

                    }
                }
                isEndIT=true;
                isEndTI=true;
                isEndPQ=true;
                isEndQP=true;
            }
            //	cout << "Done with the ends, checking for empty vertices" << endl;
            if(!done)
            {
                for(tie(vi,vi_end)=vertices(Trimmed.first); vi!=vi_end && !done; vi++)
                {
                    if(out_degree(*vi,Trimmed.first)==0 && in_degree(*vi,Trimmed.first)==0)
                    {
                        remove_vertex(*vi,Trimmed.first);
                        done = true;
                        stable = false;
                    }
                }
            }
        }		
        
    } while(!stable); // do while statement (don't panic)
    
    // Now that we have trimmed the GammaGraph, we need to check the downstairs Graph
    // First, we check for dead edges
    stable = false;
    while(!stable)
    {
        for(tie(gei,gei_end)=edges(Trimmed.second);gei!=gei_end; gei++)
        {
            found = false;
            for(tie(ei,ei_end)=edges(Trimmed.first); ei!=ei_end; ei++)
            {
                if(p_homom(*ei) == gename(*gei) || q_homom(*ei) == gename(*gei))
                {
                    found = true;
                    break;
                } // if
            } // for ei
            if(!found)
            {
                remove_edge(*gei,Trimmed.second);
                break;
            } // if found 
        } // for gei
        if(gei==gei_end)
        {
            stable = true;
        }
    }
    
    // Now we check for dead vertices
    stack<GVD> toDel;
    
    for(tie(gvi,gvi_end)=vertices(Trimmed.second); gvi!=gvi_end; gvi++)
    {
        bool found=false;
        
        for(tie(vi,vi_end)=vertices(Trimmed.first); vi!=vi_end; vi++)
        {
            if(p_vhom(*vi) == *gvi || q_vhom(*vi) == *gvi)
            {
                found = true;
                break;
            }
        }
        if(!found)
        {
            toDel.push(*gvi);
        }
        
    }
    
    if(!toDel.empty()) {
        needToCheckVhoms = true;
    }
    
    while(!toDel.empty())
    {
        GVD u = toDel.top();
        toDel.pop();
        
        clear_vertex(u,Trimmed.second);
        remove_vertex(u,Trimmed.second);
        
     //   cout << "We have deleted " << u << " from G" << endl;
    }
    
    
    if(needToCheckVhoms)
    {
        for(tie(vi,vi_end)=vertices(Trimmed.first);vi!=vi_end;vi++)
        {
            tie(oei,oei_end)=out_edges(*vi,Trimmed.first);
            found = false;
            tie(gei,gei_end)=edges(Trimmed.second);
            while(!found)
            {
                if(gei==gei_end)
                {
                    cout << "WE DIDN'T FIND A MATCH FOR A HOMOM AND AN G EDGE. THIS IS A BIG PROBLEM" << endl;
                    return T;
                }
                
                if(p_homom(*oei)==gename(*gei))
                {
                    found = true;
                    put(p_vhom,*vi,source(*gei,Trimmed.second));
                }
                else {
                    gei++;
                }
            }
            
            found = false;
            
            tie(gei,gei_end)=edges(Trimmed.second);
            while(!found)
            {
                if(gei==gei_end)
                {
                    cout << "WE DIDN'T FIND A MATCH FOR A HOMOM AND AN G EDGE. THIS IS A BIG PROBLEM" << endl;
                    return T;
                }
                
                if(q_homom(*oei)==gename(*gei))
                {
                    found = true;
                    put(q_vhom,*vi,source(*gei,Trimmed.second));
                }
                else {
                    gei++;
                }
            }
            
        }
        
    }
    
    cout << "Original T.first had " << num_edges(T.first) << " edges and " << num_vertices(T.first)
    << " vertices, now has " << num_edges(Trimmed.first) << " edges and " << num_vertices(Trimmed.first) << endl;
    
    cout << "Original T.second had " << num_edges(T.second) << " edges and " << num_vertices(T.second)
    << " vertices, now has " << num_edges(Trimmed.second) << " edges and " << num_vertices(Trimmed.second) << endl;
    
    return Trimmed;
}
*/


// This version of trim creates four sets for t(A_Gamma),i(A_Gamma),p(A_Gamma), and q(A_Gamma)
Textile ArrayTrim(Textile T)
{
	set<VD> iAGamma,tAGamma;
	set<string> pAGamma,qAGamma;
	Textile Trimmed = T;
	
	GammaEI ei,ei_end;
	GEI gei,gei_end;
	GVI gvi,gvi_end;
	GammaVI vi,vi_end;
	GammaOEI oei,oei_end;
	
	bool done=false,removed=true,needToCheckVhoms,stable=false,found;
	int roundcount = 0,vertcount=0;
	
	property_map<GammaGraph,edge_p_homom_t>::type
    p_homom = get(edge_p_homom,Trimmed.first);
    property_map<GammaGraph,edge_q_homom_t>::type
    q_homom = get(edge_q_homom,Trimmed.first);
    property_map<GammaGraph,edge_name_t>::type
    ename = get(edge_name,Trimmed.first);
    property_map<GammaGraph,vertex_name_t>::type
    vname = get(vertex_name,Trimmed.first);
	property_map<Graph,edge_name_t>::type
	gename = get(edge_name,Trimmed.second);
	property_map<GammaGraph,vertex_p_vhomom_t>::type
	p_vhom = get(vertex_p_vhomom,Trimmed.first);
  	property_map<GammaGraph,vertex_q_vhomom_t>::type
    q_vhom = get(vertex_q_vhomom,Trimmed.first);
	
	cout << "Original T.first had " << num_edges(T.first) << " edges and " << num_vertices(T.first) << endl;
	
	do {
		
		// checking for dead vertices
		do
		{
			vertcount=0;
			for(tie(vi,vi_end)=vertices(Trimmed.first); vi!=vi_end; vi++)
			{
				if(out_degree(*vi,Trimmed.first)==0 || in_degree(*vi,Trimmed.first)==0)
				{
		//					cout << "Deleting vertex " << *vi << endl;
					clear_vertex(*vi,Trimmed.first);
					remove_vertex(*vi,Trimmed.first);
					vertcount++;
					break;
				}
			} // for vi
		} while(vertcount>0);


		roundcount = 0;
	// First we'll populate the sets with the current version of Trimmed.	
		cout << "Creating the sets" << endl;
		for(tie(ei,ei_end)=edges(Trimmed.first); ei!=ei_end; ei++)
		{

			iAGamma.insert(source(*ei,Trimmed.first));
			tAGamma.insert(target(*ei,Trimmed.first));
			pAGamma.insert(p_homom(*ei));
			qAGamma.insert(q_homom(*ei));
		}
		
		cout << "iAGamma size is " << iAGamma.size() << endl;
		cout << "Sets created, moving on to removed while loop" << endl;	
		
	// Now that we've populated the sets, we can go through the edges and check their properties.
	while(removed) {
//		cout << "Back in removed while loop" << endl; 
		removed = false;
		for(tie(ei,ei_end)=edges(Trimmed.first); ei!=ei_end; ei++)
		{
			if(iAGamma.find(target(*ei,Trimmed.first))==iAGamma.end())
			{
//				cout << "Deleting edge from " << source(*ei,Trimmed.first) << " to " << target(*ei,Trimmed.first) << endl;
				roundcount++;
				removed = true;
				remove_edge(*ei,Trimmed.first);
				break;
			}
			if(tAGamma.find(source(*ei,Trimmed.first))==tAGamma.end())
			{
//				cout << "Deleting edge from " << source(*ei,Trimmed.first) << " to " << target(*ei,Trimmed.first) << endl;
				roundcount++;
				removed = true;
				remove_edge(*ei,Trimmed.first);
				break;
			}
			if(pAGamma.find(q_homom(*ei))==pAGamma.end())
			{
//				cout << "Deleting edge from " << source(*ei,Trimmed.first) << " to " << target(*ei,Trimmed.first) << endl;
				roundcount++;
				removed = true;
				remove_edge(*ei,Trimmed.first);
				break;
			}
			if(qAGamma.find(p_homom(*ei))==qAGamma.end())
			{
//				cout << "Deleting edge from " << source(*ei,Trimmed.first) << " to " << target(*ei,Trimmed.first) << endl;
				roundcount++;
				removed = true;
				remove_edge(*ei,Trimmed.first);
				break;
			}
		} // for ei
	} // while

	cout << "Done with removed while loop." << endl;

		// checking for dead vertices
		do
		{
			vertcount=0;
			for(tie(vi,vi_end)=vertices(Trimmed.first); vi!=vi_end; vi++)
			{
				if(out_degree(*vi,Trimmed.first)==0 || in_degree(*vi,Trimmed.first)==0)
				{
//					cout << "Deleting vertex " << *vi << endl;
					clear_vertex(*vi,Trimmed.first);
					remove_vertex(*vi,Trimmed.first);
					vertcount++;
					break;
				}
			} // for vi
		} while(vertcount>0);
		cout << "Roundcount is " << roundcount << endl;
	//	cout << "Clearing sets" << endl;
		
		iAGamma.clear();
		tAGamma.clear();
		pAGamma.clear();
		qAGamma.clear();
		
		
	} while(roundcount > 0);

   // Now that we have trimmed the GammaGraph, we need to check the downstairs Graph
    // First, we check for dead edges
    while(!stable)
    {
		stable = true;
        for(tie(gei,gei_end)=edges(Trimmed.second);gei!=gei_end; gei++)
        {
            found = false;
            for(tie(ei,ei_end)=edges(Trimmed.first); ei!=ei_end; ei++)
            {
                if(p_homom(*ei) == gename(*gei) || q_homom(*ei) == gename(*gei))
                {
                    found = true;
                    break;
                } // if
            } // for ei
            if(!found)
            {
				stable = false;
                remove_edge(*gei,Trimmed.second);
                break;
            } // if found 
        } // for gei
    }

	cout << "Done removing downstairs edges" << endl;
    
    // Now we check for dead vertices
    stack<GVD> toDel;
    
    for(tie(gvi,gvi_end)=vertices(Trimmed.second); gvi!=gvi_end; gvi++)
    {
        bool found=false;
        
        for(tie(vi,vi_end)=vertices(Trimmed.first); vi!=vi_end; vi++)
        {
            if(p_vhom(*vi) == *gvi || q_vhom(*vi) == *gvi)
            {
                found = true;
                break;
            } 
        } // for vi
        if(!found)
        {
            toDel.push(*gvi);
        } // if !found
    }
    
    if(!toDel.empty()) {
        needToCheckVhoms = true;
    }
    
    while(!toDel.empty())
    {
        GVD u = toDel.top();
        toDel.pop();
        
        clear_vertex(u,Trimmed.second);
        remove_vertex(u,Trimmed.second);
        
     //   cout << "We have deleted " << u << " from G" << endl;
    }
    
	cout << "Do we need to check vHoms? :" << needToCheckVhoms << endl;
	
	cout << "Original T.first had " << num_edges(T.first) << " edges and " << num_vertices(T.first)
    << " vertices, now has " << num_edges(Trimmed.first) << " edges and " << num_vertices(Trimmed.first) << endl;

    cout << "Original T.second had " << num_edges(T.second) << " edges and " << num_vertices(T.second)
    << " vertices, now has " << num_edges(Trimmed.second) << " edges and " << num_vertices(Trimmed.second) << endl;
   
	
    if(needToCheckVhoms)
    {
        for(tie(vi,vi_end)=vertices(Trimmed.first);vi!=vi_end;vi++)
        {
			cout << "Out degree of vi is " << out_degree(*vi,Trimmed.first) << " and in degree is " << in_degree(*vi,Trimmed.first) << endl;
            tie(oei,oei_end)=out_edges(*vi,Trimmed.first);
            found = false;
            tie(gei,gei_end)=edges(Trimmed.second);
            while(!found)
            {
                if(gei==gei_end)
                {
                    cout << "WE DIDN'T FIND A MATCH FOR A HOMOM AND AN G EDGE. THIS IS A BIG PROBLEM" << endl;
                    return T;
                }
                
                if(p_homom(*oei)==gename(*gei))
                {
                    found = true;
                    put(p_vhom,*vi,source(*gei,Trimmed.second));
					cout << "Assigned " << source(*gei,Trimmed.second) << "as the p_vhom for " << *vi << endl;
                }
                else {
                    gei++;
                }
            }
            
            found = false;
            
            tie(gei,gei_end)=edges(Trimmed.second);
            while(!found)
            {
                if(gei==gei_end)
                {
                    cout << "WE DIDN'T FIND A MATCH FOR A HOMOM AND AN G EDGE. THIS IS A BIG PROBLEM" << endl;
                    return T;
                }
                
                if(q_homom(*oei)==gename(*gei))
                {
                    found = true;
                    put(q_vhom,*vi,source(*gei,Trimmed.second));
					cout << "Assigned " << source(*gei,Trimmed.second) << "as the q_vhom for " << *vi << endl;

                }
                else {
                    gei++;
                }
            }
            
        }
        
    }
    
    cout << "Original T.first had " << num_edges(T.first) << " edges and " << num_vertices(T.first)
    << " vertices, now has " << num_edges(Trimmed.first) << " edges and " << num_vertices(Trimmed.first) << endl;
    
    cout << "Original T.second had " << num_edges(T.second) << " edges and " << num_vertices(T.second)
    << " vertices, now has " << num_edges(Trimmed.second) << " edges and " << num_vertices(Trimmed.second) << endl;
    
    return Trimmed;

}


Textile NewInducedRp(Textile T)
{
	GammaGraph Gamma;
    Graph G=T.second;

 
    property_map<GammaGraph,edge_p_homom_t>::type
    pe = get(edge_p_homom,Gamma);

   	property_map<GammaGraph,edge_name_t>::type
    gamma_ename = get(edge_name,Gamma);

    property_map<GammaGraph,vertex_name_t>::type
    gamma_vname = get(vertex_name,Gamma);

    property_map<Graph,edge_name_t>::type
    gename = get(edge_name,G);

    property_map<Graph,vertex_name_t>::type
    g_vname = get(vertex_name,G);

    property_map<GammaGraph,vertex_p_vhomom_t>::type
    pvgamma = get(vertex_p_vhomom, Gamma);
  
	
	GammaVI vi,vi_end;
	GEI gei,gei_end;
	GammaOEI oei,oei_end;
	
	set<VVec> cSets;
	set<VVec>::iterator li;
	
	stack<VVec> toExamine;
	
	vector<string> codex;
	vector<string>::iterator sit;
	
	property_map<Graph,edge_name_t>::type
    tsename = get(edge_name,T.second);

	VVec::iterator si;

	bool added,isSubset,found=false;
	

    // Insert all the singletons onto our stack of sets to examine and into the compatible sets
	for(tie(vi,vi_end)=vertices(T.first);vi!=vi_end;vi++)
	{
		VD v = *vi;
		VVec singleton(1,v);
		cout << "Pushing " << v << endl;
		toExamine.push(singleton);
		cSets.insert(singleton);
	}
	

    // Fix the codex with all of the G edge names
	for(tie(gei,gei_end)=edges(T.second);gei!=gei_end;gei++)
	{
		string s = tsename(*gei);
		codex.push_back(s);
	}
	

    // Start going through the stack of sets to examine
	while(!toExamine.empty())
	{
		VVec currv = toExamine.top();
		toExamine.pop();
		
		for(sit=codex.begin(); sit < codex.end(); sit++)
		{
			VVec S = compatibleSet(T,1,currv,*sit);
			if(S.size() > 0)
			{
				cout << "Looking at CS " << ssVVec(S) << endl;
				// We have a (sorted) compatible set, we want to check if its already seen

                if(cSets.find(S)==cSets.end())
                {
                    cSets.insert(S);
                    toExamine.push(S);
                }

				/*added = false;
				isSubset = false;
				for(li = cSets.begin(); !isSubset && li!=cSets.end();li++)
				{
					cout << "Comparing it to " << ssVVec(*li) << endl;
					if(VVecSubset(*li,S))
					{
						cout << "*li is a subset" << endl;
						li = cSets.erase(li);
						li--;
						if(!added)
						{
						cSets.push_back(S);
						added = true;
						}
					}
					else if(VVecSubset(S,*li))
					{
			//			cout << "S is a subset" << endl;
						isSubset = true;
					}
                    
				}*/
			}
		}
		
	}
	
	
	VertexMap vmap;
    
	
	cout << "Our final Seen consists of:";
    
    for(li=cSets.begin(); li != cSets.end(); li++)
    {
        
        graph_traits<GammaGraph>::vertex_descriptor v = add_vertex(Gamma);
        VVec::iterator it;
        string vname;
        VVec curr=*li;
        
        printVVec(curr);
        vname = ssVVec(curr);
        put(gamma_vname,v,vname);
        
        //We need to create a VMap so we can add edges later
        vmap.insert(VertexMap::value_type(vname,v));
    }
    cout << endl;
    
    for(li=cSets.begin(); li != cSets.end(); li++)
    {
        string s=ssVVec(*li),t;
        graph_traits<GammaGraph>::vertex_descriptor v=vmap[s],w;
        vector<string>::iterator cit;
        VVec S;
        
        for(cit = codex.begin(); cit < codex.end(); cit++)
        {
            S = compatibleSet(T,1,(*li),(*cit));
            if(S.size() > 0) {
                t = ssVVec(S);
                w = vmap[t];
                add_edge(v,w,PQ_Homoms((*cit),Q_Homom(string(""),string("(" + s + "," + (*cit) + ")"))),Gamma);
            }  
        }
        
        
    }

   // We now need to fill in the vertex homomorphisms
    for(tie(vi,vi_end)=vertices(Gamma); vi!=vi_end; vi++)
    {
        string currPE;
        tie(oei,oei_end)=out_edges(*vi,Gamma);
        
        currPE = pe(*oei);
        cout << "In order to place a vertex homom on " << gamma_vname(*vi) << " we are looking at " << currPE << endl;
        tie(gei,gei_end)=edges(G);
        while(!found)
        {
            if(currPE == gename(*gei))
            {
                cout << "We found a match between " << currPE << " and " << gename(*gei) << " which has source " << g_vname(source(*gei,G)) << endl;
                put(pvgamma,*vi,source(*gei,G));
                found = true;
            }
            gei++;
        } // while not found
        found = false;
    } // for vi
	
    cout << "Done with creating induced presentation" << endl;

	return Textile(Gamma,G);
	
}

/*
 * Given graphs G and H and a specified equivalence between M_G M_H and M_H M_G, we create the corresponding
 * textile system T = (p',q': Gamma -> G) for this specified equivalence. 
 */
Textile FromSSE(Graph G, Graph H, std::unordered_map<string,string> sequiv)
{
    bool found;
    vector<int> RHSsizes(2),LHSsizes(2);

    GammaGraph Gamma(num_edges(H));

    Graph Gprime = G;

    GammaVI vi,vi_end;
    GEI gei,gei_end;
    GammaOEI oei,oei_end;

    std::unordered_map<string,VD> AHNameToGammaVD;
    std::unordered_map<string,GVD> GLookup;

    property_map<GammaGraph,vertex_name_t>::type
    Gamma_vname = get(vertex_name,Gamma);

    property_map<Graph,edge_name_t>::type
    H_ename = get(edge_name,H);

    property_map<Graph,edge_name_t>::type
    G_ename = get(edge_name,G);

    property_map<GammaGraph,vertex_p_vhomom_t>::type
    Gamma_pvhom = get(vertex_p_vhomom,Gamma);

    property_map<GammaGraph,vertex_q_vhomom_t>::type
    Gamma_qvhom = get(vertex_q_vhomom,Gamma);

    property_map<Graph,vertex_name_t>::type
    G_vname = get(vertex_name,G);

    property_map<GammaGraph,edge_p_homom_t>::type
    Gamma_phom = get(edge_p_homom,Gamma);    

    property_map<GammaGraph,edge_q_homom_t>::type
    Gamma_qhom = get(edge_q_homom,Gamma);

    property_map<Graph,edge_name_t>::type
    Gprime_ename = get(edge_name,Gprime);

    tie(gei,gei_end)=edges(G);
    LHSsizes[0]=G_ename(*gei).length();
    RHSsizes[1]=G_ename(*gei).length();

    tie(gei,gei_end)=edges(H);
    LHSsizes[1]=H_ename(*gei).length();
    RHSsizes[0]=H_ename(*gei).length();

    // First let's add names to the Gamma Graph vertices, we will also create the map for the AHnames
    for(tie(vi,vi_end)=vertices(Gamma), tie(gei,gei_end)=edges(H); vi!=vi_end; vi++,gei++)
    {
        cout << "Adding " << H_ename(*gei) << endl;
        put(Gamma_vname,*vi,H_ename(*gei));
        AHNameToGammaVD[H_ename(*gei)] = *vi;
    }

    // Now we need to add the Gamma Graph edges and assign their properties
    for(auto it = sequiv.begin(); it != sequiv.end(); ++it)
    {
        string currLHS = (*it).first,currRHS = (*it).second,a,b,aprime,bprime,name;
        cout << "Attempting to split " << currLHS << " and " << currRHS << endl;
        vector<string> LHSSplit = StringSplitter(currLHS,LHSsizes);
        vector<string> RHSSplit = StringSplitter(currRHS,RHSsizes);

        a = LHSSplit[0];
        b = RHSSplit[0];
        aprime = RHSSplit[1];
        bprime = LHSSplit[1];
        name = a+b+aprime+bprime;
        cout << "a = " << a << ", b = " << b << ", aprime = " << aprime << ", bprime = " << bprime  << endl;
        add_edge(AHNameToGammaVD[b], AHNameToGammaVD[bprime], PQ_Homoms(a,Q_Homom(aprime,name)),Gamma);
    }

   // We now need to fill in the vertex homomorphisms
    for(tie(vi,vi_end)=vertices(Gamma); vi!=vi_end; vi++)
    {
        string currPE;
        tie(oei,oei_end)=out_edges(*vi,Gamma);
        cout << "Looking at vertex: " << Gamma_vname(*vi) << endl;
        currPE = Gamma_phom(*oei);
        tie(gei,gei_end)=edges(Gprime);
	cout << "Looking for " << currPE << endl;
        while(!found)
        {
	  cout << "Currently seeing " << Gprime_ename(*gei) << endl;
            if(currPE == Gprime_ename(*gei))
            {
                put(Gamma_pvhom,*vi,source(*gei,Gprime));
                found = true;
            }
            gei++;
        } // while not found
        found = false;
    } // for vi

    for(tie(vi,vi_end)=vertices(Gamma); vi!=vi_end; vi++)
    {
        string currQE;
        tie(oei,oei_end)=out_edges(*vi,Gamma);
        
        currQE = Gamma_qhom(*oei);
        tie(gei,gei_end)=edges(Gprime);
        while(!found)
        {
            if(currQE == Gprime_ename(*gei))
            {
                put(Gamma_qvhom,*vi,source(*gei,Gprime));
                found = true;
            }
            gei++;
        } // while not found
        found = false;
    } // for vi



    return Textile(Gamma,Gprime);
}
/**
Textile GenerateRandomTextileFromSSE(Graph G, Graph H)
{

    int i;
    std::unordered_map<string,string> SpecEq;

    for(i=0;i<n;i++)
    {

    }

    return T;
}**/

Graph ProductGraph(Graph G, Graph H)
{
    Graph GH(num_vertices(G));


    property_map<Graph,vertex_name_t>::type
    G_vname = get(vertex_name,G);

    property_map<Graph,vertex_name_t>::type
    GH_vname = get(vertex_name,GH);

    property_map<Graph,edge_name_t>::type
    G_ename = get(edge_name,G);

    property_map<Graph,edge_name_t>::type
    GH_ename = get(edge_name,GH);

    property_map<Graph,edge_name_t>::type
    H_ename = get(edge_name,H);


    GVI tgvi,tgvi_end,pgvi,pgvi_end;
    GOEI tgoei,tgoei_end,sgoei,sgoei_end;

       // Add names
    for(tie(tgvi,tgvi_end)=vertices(G),tie(pgvi,pgvi_end)=vertices(GH);
        tgvi!=tgvi_end; tgvi++, pgvi++) {
        put(GH_vname,*pgvi,G_vname(*tgvi));
    }
    
    // We now need to create the edges of GH
    for(tie(tgvi,tgvi_end)=vertices(G); tgvi!=tgvi_end;tgvi++)
    {
        for(tie(tgoei,tgoei_end)=out_edges(*tgvi,G);tgoei!=tgoei_end;tgoei++)
        {
            for(tie(sgoei,sgoei_end)=out_edges(target(*tgoei,G),H);sgoei!=sgoei_end;sgoei++)
            {
                add_edge(*tgvi,target(*sgoei,H),string( G_ename(*tgoei) + H_ename(*sgoei)),GH);
            }
        }
    }


    return GH;
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

bool SEquivChecker(Graph& GH, Graph& HG, std::unordered_map<string,string> sequiv)
{
  bool validse=true,found=false;

  property_map<Graph,vertex_name_t>::type
    GH_vname = get(vertex_name,GH);

  property_map<Graph,vertex_name_t>::type
    HG_vname = get(vertex_name,HG);

  property_map<Graph,edge_name_t>::type
    GH_ename = get(edge_name,GH);
  
  property_map<Graph,edge_name_t>::type
    HG_ename = get(edge_name,HG);

  GEI ei,ei_end,fi,fi_end;

  for(tie(ei,ei_end)=edges(GH);ei!=ei_end;ei++)
    {
      found = false;
      string currName = sequiv[GH_ename(*ei)];
      tie(fi,fi_end)=edges(HG);
      while(!found && validse && fi!=fi_end)
	{
	  if(currName == HG_ename(*fi))
	    {
	      found = true;
          cout << "Checking that " << GH_vname(source(*ei,GH)) << " matches " << HG_vname(source(*fi,HG)) << endl;
          cout << "Checking that " << GH_vname(target(*ei,GH)) << " matches " << HG_vname(target(*fi,HG)) << endl;
	      if(GH_vname(source(*ei,GH)) != HG_vname(source(*fi,HG)))
		{
		  /*
		  cout << "GH_vname is " << GH_vname(source(*ei,GH)) << " while HG_vname is "
		       << HG_vname(source(*fi,HG)) << endl;
		  */
		  validse = false;
		}
	      if(GH_vname(target(*ei,GH)) != HG_vname(target(*fi,HG)))
		{
		  /*
		  cout << "GH_vname is " << GH_vname(target(*ei,GH)) << " while HG_vname is "
		       << HG_vname(target(*fi,HG)) << endl;
		  */
		  validse = false;
		  }
	    }
	  fi++;
	}

    }

    return validse;
}
