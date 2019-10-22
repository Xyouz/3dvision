#include <Imagine/Images.h>
#include <iostream>
#include <fstream>

#include "maxflow/graph.h"

using namespace std;
using namespace Imagine;

// This section shows how to use the library to compute a minimum cut on the following graph:
//
//		        SOURCE
//		       /       \
//		     1/         \6
//		     /      4    \
//		   node0 -----> node1
//		     |   <-----   |
//		     |      3     |
//		     \            /
//		     5\          /1
//		       \        /
//		          SINK
//
///////////////////////////////////////////////////

void testGCuts()
{
	Graph<int,int,int> g(/*estimated # of nodes*/ 4, /*estimated # of edges*/ 4); 
	g.add_node(4); 
	
	g.add_tweights( 0,   /* capacities */  3, 0 );
	g.add_tweights( 1,   /* capacities */  2, 0 );
	g.add_tweights( 2,   /* capacities */  0, 6 );
	g.add_tweights( 3,   /* capacities */  0, 5 );
	
	g.add_edge( 0, 3,    /* capacities */  1, 0 );
	g.add_edge( 0, 2,    /* capacities */  1, 0 );
	g.add_edge( 1, 3,    /* capacities */  1, 0 );
	g.add_edge( 1, 2,    /* capacities */  1, 0 );

	int flow = g.maxflow();
	cout << "Flow = " << flow << endl;
	for (int i=0;i<4;i++)
		if (g.what_segment(i) == Graph<int,int,int>::SOURCE)
			cout << i << " is in the SOURCE set" << endl;
		else
			cout << i << " is in the SINK set" << endl;
}


int main() {
	testGCuts();
	char dummy; cin >> dummy;
	return 0;
}
