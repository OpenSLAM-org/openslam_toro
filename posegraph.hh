/**********************************************************************
 *
 * This source code is part of the Tree-based Network Optimizer (TORO)
 *
 * TORO Copyright (c) 2007 Giorgio Grisetti, Cyrill Stachniss, and
 * Wolfram Burgard
 *
 * TORO is licences under the Common Creative License,
 * Attribution-NonCommercial-ShareAlike 3.0
 *
 * You are free:
 *   - to Share - to copy, distribute and transmit the work
 *   - to Remix - to adapt the work
 *
 * Under the following conditions:
 *
 *   - Attribution. You must attribute the work in the manner specified
 *     by the author or licensor (but not in any way that suggests that
 *     they endorse you or your use of the work).
 *  
 *   - Noncommercial. You may not use this work for commercial purposes.
 *  
 *   - Share Alike. If you alter, transform, or build upon this work,
 *     you may distribute the resulting work only under the same or
 *     similar license to this one.
 *
 * Any of the above conditions can be waived if you get permission
 * from the copyright holder.  Nothing in this license impairs or
 * restricts the author's moral rights.
 *
 * TORO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  
 **********************************************************************/

/** \file posegraph.hh 
 * 
 * \brief The template class for the node parameters. The graph of
 * poses with support to tree construction functionalities.
 **/

#ifndef _TREEPOSEGRAPH_HXX_
#define _TREEPOSEGRAPH_HXX_

#include <iostream>
#include <assert.h>

#include <set>
#include <list>
#include <map>
#include <deque>
#include <vector>
#include <values.h>
#include <algorithm>

/** \brief The template class for representing an abstract tree
    without specifing the dimensionality of the exact parameterization
    of the nodes. This definition is passed in via the Operation (Ops)
    template class **/
template <class Ops>
struct TreePoseGraph{
  typedef typename Ops::PoseType Pose;
  typedef typename Ops::RotationType Rotation;
  typedef typename Ops::TranslationType Translation;
  typedef typename Ops::TransformationType Transformation;
  typedef typename Ops::CovarianceType Covariance;
  typedef typename Ops::InformationType Information;
  typedef typename Ops::ParametersType  Parameters;
  
  struct Vertex;

  /** \brief Definition of an edge in the graph based on the template
      input from Ops **/
    struct Edge{
    Vertex* v1;   /**< The constraint is defined between v1 and v2 **/
    Vertex* v2;   /**< The constraint is defined between v1 and v2 **/
    Vertex* top;  /**< The node with the smallest level in the path **/
    int length;   /**< Length of the path on the tree (number of vertieces involved) **/
    Transformation transformation;    /**< Transformation describing the constraint (relative mapping) **/
    Information    informationMatrix; /**< Uncertainty encoded in the information matrix **/
  };

  typedef typename std::list< Edge* > EdgeList;
  typedef typename std::map< int, Vertex* > VertexMap;
  typedef typename std::map< Edge*, Edge* > EdgeMap;

  /** \brief Definition of a vertex in the graph based on the
      template input from Ops **/
  struct Vertex {

    // Graph-related elements
    int id;         /**< Id of the vertex in the graph **/
    EdgeList edges; /**< The edges related to this vertex **/

    // Tree-related elements
    int level;         /**< level in the tree. It is the distance on the tree to the root **/
    Vertex* parent;    /**< Parent vertex **/
    Edge* parentEdge;  /**< Constraint between the parent and the current vertex in the tree **/
    EdgeList children;   /**< All constraints involving the children of this vertex **/

    // Parameterization-related elements
    Pose pose;              /**< The pose of the vertex **/
    Parameters parameters;  /**< The parameter representation **/
  };

  /** Returns the vertex with the given id **/
  Vertex* vertex(int id);

  /** Returns a const pointer to the vertex with the given id **/
  const Vertex* vertex (int id) const;

  /** Returns the edge between the two vertices **/
  Edge* edge(int id1, int id2);

  /** Returns a const pointer tothe edge between the two vertices **/
  const Edge* edge(int id1, int id2) const;

  /** Add a vertex to the graph **/
  Vertex* addVertex(int id, const Pose& pose);
  /** Remove a vertex from the graph **/
  Vertex* removeVertex (int id);

  /** Add an edge/constraint to the graph **/
  Edge* addEdge(Vertex* v1, Vertex* v2, const Transformation& t, const Information& i);

  /** Remove an edge/constraint from the graph **/
  Edge* removeEdge(Edge* eq);

  /** Function to perform a breadth-first visit of the nodes in the tree to carry out a specific action act**/
  template <class Action>
  void treeBreadthVisit(Action& act);

  /** Function to perform a depth-first visit of the nodes in the tree to carry out a specific action act **/
  template <class Action>
  void treeDepthVisit(Action& act, Vertex *v);

  /** Constructs the tree be computing a minimal spanning tree **/
  bool buildMST(int id);

  /** Constructs the incremental tree according to the input trajectory **/
  bool buildSimpleTree();

  /** Trun around an edge (used to ensure a certain oder on the vertexes) **/
  void revertEdge(Edge* e);

  /** Revert edge info. This function needs to be implemented by a subclass **/
  virtual void revertEdgeInfo(Edge* e) = 0;

  /** Delete all edges and vertices **/
  void clear();

  /** Destructor **/
  virtual ~TreePoseGraph();

  /** Sort constraints for correct processing order **/
  EdgeList* sortEdges();

  /** Determines the length of the longest path in the tree **/
  int maxPathLength();

  /** Determines the path length of all pathes in the tree **/
  int totalPathLength();

  /** remove gaps in the indices of the vertex ids **/
  void compressIndices();

  /** compute the highest index of an vertex **/
  int maxIndex();

  /** The root node of the tree **/
  Vertex* root;

  /** All vertices **/
  VertexMap vertices;

  /** All edges **/
  EdgeMap edges;

protected:
  void fillEdgeInfo(Edge* e);
  void fillEdgesInfo();


};

//include the template implementation part

#include "posegraph.hxx"

#endif
