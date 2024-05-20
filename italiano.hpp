/**
 * \file italiano.hpp
 *   \author Marcus Ritt <mrpritt@inf.ufrgs.br> 
 *   \version $Id: emacs 4749 2013-08-16 11:47:04Z ritt $
 *   \date Time-stamp: <2013-11-27 15:09:45 ritt>
 *
 * A simple implementation of Italiano's (1986) semi-dynamic
 * (incremental) transitive closure.
 */
#pragma once

#include <vector>

typedef unsigned ItalianoNode;



template <unsigned Size>
struct ClosedGraph {
  struct Vertex {
	  Vertex(ItalianoNode _node)
  : node(_node), firstchild(0), nextsibling(0)
    {}
    ItalianoNode node;
    Vertex *firstchild, *nextsibling;
  };

  ClosedGraph(unsigned _n) : n(_n), cycle(false)
#if defined(ITALIANO_WITH_COUNTER)
	, a(0)
#endif
{
    for(unsigned i=0; i<n; i++)
      for(unsigned j=0; j<n; j++)
	r[i][j]=nullptr;
    for(unsigned i=0; i<n; i++)
      r[i][i]=new Vertex(i);
  }
  ~ClosedGraph() {
    for(unsigned i=0; i<n; i++)
      for(unsigned j=0; j<n; j++)
	if (r[i][j]!=nullptr)
	  delete r[i][j];
  }
  ClosedGraph(const ClosedGraph& other) : n(other.n), cycle(other.cycle) {
#if defined(ITALIANO_WITH_COUNTER)
	a=other.a;
#endif
    for(unsigned i=0; i<n; i++)
      for(unsigned j=0; j<n; j++)
	if (other.r[i][j]!=nullptr)
	  r[i][j]=new Vertex(other.r[i][j]->node);
	else
	  r[i][j]=nullptr;
    for(unsigned i=0; i<n; i++)
      for(unsigned j=0; j<n; j++)
	if (other.r[i][j]!=nullptr) {
	  Vertex *w = other.r[i][j]->firstchild;
	  while (w != nullptr) {
	    r[i][w->node]->nextsibling = r[i][j]->firstchild;
	    r[i][j]->firstchild = r[i][w->node];
	    w = w->nextsibling;
	  }
	}
  }

  ClosedGraph(ClosedGraph&& other) : ClosedGraph(other.n) {
    rip(*this,other);
  }

  ClosedGraph& operator=(ClosedGraph other) {
    rip(*this,other);
    return *this;
  }
  void rip(ClosedGraph& cga, ClosedGraph& cgb) {
    using std::swap;
    swap(cga.n,cgb.n);
    swap(cga.cycle,cgb.cycle);
    for(unsigned i=0; i<n; i++)
      for(unsigned j=0; j<n; j++)
	      swap(cga.r[i][j],cgb.r[i][j]);
	  swap(cga.a, cgb.a);
  }
  
  // is there a path from u to v?
  bool path(ItalianoNode u, ItalianoNode v) const {
    return r[u][v] != nullptr;
  }
  // add arc (u,v)
  void add(ItalianoNode u, ItalianoNode v) {
    if (r[u][v] != nullptr)
      return;
    // update all nodes
    for(ItalianoNode x=0; x<n; x++)
      if (r[x][u] != nullptr && r[x][v] == nullptr)
	meld(x,v,u,v);
    if (r[v][u]!=nullptr)
      cycle = true;
  }
  bool hasCycle() const {
    return cycle;
  }

#if defined(ITALIANO_WITH_COUNTER)
  bool isLinearOrder() const {
    return n*(n-1)==2*a;
  }
#endif

  //protected:
  void meld(ItalianoNode x, ItalianoNode y, ItalianoNode u, ItalianoNode v) {
    r[x][v] = new Vertex(v);
#if defined(ITALIANO_WITH_COUNTER)
    a++;
#endif
    r[x][v]->nextsibling = r[x][u]->firstchild;
    r[x][u]->firstchild = r[x][v];
    Vertex *w = r[y][v]->firstchild;
    while (w != nullptr) {
      if (r[x][w->node] == nullptr)
	meld(x,y,v,w->node);
      w = w->nextsibling;
    }
  }

  unsigned n;
  bool cycle;
  Vertex *r[Size][Size];
#if defined(ITALIANO_WITH_COUNTER)
  unsigned a; // number of arcs
#endif
};
