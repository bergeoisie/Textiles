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
    int i,j=1,k,l;
    bool found=false;
    Graph GM(2),HM(2),GH,HG,PM,MP,M;
	std::unordered_map<string,string> sequiv;

    property_map<Graph,vertex_name_t>::type
    GM_vname = get(vertex_name,GM); 

    property_map<Graph,vertex_name_t>::type
    HM_vname = get(vertex_name,HM); 

    property_map<Graph,edge_name_t>::type
    GH_ename = get(edge_name,GH);

    property_map<Graph,edge_name_t>::type
    HG_ename = get(edge_name,HG);
    
    property_map<Graph,edge_name_t>::type
    PM_ename = get(edge_name,PM);

    property_map<Graph,edge_name_t>::type
    MP_ename = get(edge_name,MP);

    GEI ei,ei_end,fi,fi_end;

    ofstream grid("gridGMSSE.txt"),loSEs("listOfSEs.txt");

    add_edge(0,0,string("u"),GM);
    add_edge(0,1,string("v"),GM);
    add_edge(1,0,string("w"),GM);
    
    put(GM_vname,0,string("A"));
    put(GM_vname,1,string("B"));


    PrintGraph(GM);

    add_edge(0,0,string("x"),HM);
    add_edge(0,1,string("y"),HM);
    add_edge(1,0,string("z"),HM);
    
    put(HM_vname,0,string("A"));
    put(HM_vname,1,string("B"));


    M = ProductGraph(GM,HM);

    MP = ProductGraph(M,GM);

    PM = ProductGraph(GM,M);

    //    HG = ProductGraph(HM,GM);

    PrintGraph(M);

    PrintGraph(MP);

    PrintGraph(PM);

  /*  sequiv[string("xu")]=string("vz");
    sequiv[string("yw")]=string("ux");
    sequiv[string("xv")]=string("uy");
    sequiv[string("zu")]=string("wx");
    sequiv[string("zv")]=string("wy");
*/
  //  Textile T = FromSSE(GM,HM,sequiv);

//    PrintFullTextileInfo(T);

    vector<graph_traits<Graph>::edge_descriptor> permTest(num_edges(MP));

    for(tie(ei,ei_end)=edges(PM),i=0;ei!=ei_end;ei++,i++)
    {
        permTest[i]=*ei;
    } 

    auto start = permTest.begin(),finish = permTest.end();
    do {
        for(tie(ei,ei_end)=edges(MP),tie(fi,fi_end)=edges(PM),i=0;ei!=ei_end;ei++,i++,fi++)
        {
      //  cout<< MP_ename(*ei) << " and " << PM_ename(*fi) << " and " << PM_ename(permTest[i]) << endl;
            sequiv[MP_ename(*ei)]=PM_ename(permTest[i]);
	//	    cout << MP_ename(*ei) << " -> " << sequiv[MP_ename(*ei)] << endl;
        } 

    //  cout << endl;

        if(SEquivChecker(MP,PM,sequiv))
        {
           /* if(j==5046)
            { */
                cout << "Beginning printing the " << j << "th SE " << endl;
                for(tie(ei,ei_end)=edges(MP),i=0;ei!=ei_end;ei++,i++) {
                    cout << MP_ename(*ei) << " -> " << sequiv[MP_ename(*ei)] << endl;
                }
                cout << "Printing trimmed SSE textile" << endl;
                Textile T = Trim(FromSSE(M,GM,sequiv));

	
	if(is1to1(T))
	  {
	    cout << "T is 1-1" << endl;
	    Textile Td = CreateDual(T);

        PrintFullTextileInfo(Td);
	    cout << "Is Td 1-1?" << is1to1(Td) << endl;
	  }
	else
	  {
	    cout << "T is not 1-1" << endl;
	  }
        for(k=-3; k>-4;k--)
        {
            for(l=1; l<2; l++)
            {
                Textile Tkl = AutoHomom(CreateNMTextile(T,k,l));
                if(IsLR(Tkl))
                {
                    grid << "We are LR at SE " << j << " where k = " << k << " and l = "  << l << endl;
                }
                else{
                    grid << "We are NOT LR at SE " << j << " where k = " << k << " and l = "  << l << endl;
                }
                cout << "CHECKING CONJUGACY FOR k = " << k << " and l = " << l << endl;
                Textile Tconj = LookForConjugacy(Tkl,4,&found);
                cout << "Found = " <<  found << endl;
                if(found){
                    grid << "PRINTING IJ CONJUGACY FOR " << i << j << endl;
                    PrintFullTextileInfo(Tkl,grid);
                    PrintFullTextileInfo(Tconj,grid);
                }
            } // l for loop
        } // k for loop  
    //}
    }

j++;
} while(std::next_permutation(start,finish));


	return 0;
}
