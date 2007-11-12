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

/** \file posegraph3.hh
 *
 * \brief Defines the graph of 3D poses, with specific functionalities
 * such as loading, saving, merging constraints, and etc.
 **/

#ifndef _POSEGRAPH3_HH_
#define _POSEGRAPH3_HH_

#include "posegraph.hh"
#include "transformation3.hh"
#include <iostream>
#include <vector>



/** \brief The class (struct) that contains 2D graph related functions
    such as loading, saving, merging, etc. **/
struct TreePoseGraph3: public TreePoseGraph<Operations3D<double> >{
  typedef Operations3D<double> Ops;
  
  typedef Ops::PoseType           Pose;
  typedef Ops::RotationType       Rotation;
  typedef Ops::TranslationType    Translation;
  typedef Ops::TransformationType Transformation;
  typedef Ops::CovarianceType     CovarianceMatrix;
  typedef Ops::InformationType    InformationMatrix;


  /** Load a graph from a file ignoring the equivalence constraints
      @param filename the graph file
      @param overrideCovariances ignore the covariances from the file, and use identities instead
   **/
  bool load( const char* filename, bool overrideCovariances=false, bool twoDimensions=false); 

  /** Load only the equivalence constraints from a  graph file (call load before)  **/
  bool loadEquivalences( const char* filename); 

  /** Saves the graph in the graph-format**/
  bool save( const char* filename);

  /** Saved the graph for visualizing it using gnuplot **/
  bool saveGnuplot( const char* filename);

  /** Debug function  **/
  void printDepth( std::ostream& os );

  /** Debug function  **/
  void printWidth( std::ostream& os );

  /** Debug function  **/
  void printEdgesStat( std::ostream& os);

  void initializeOnTree();
  
  virtual void initializeFromParentEdge(Vertex* v);

  /** Turn around the edge (<i,j>  => <j,i>)  **/
  virtual void revertEdgeInfo(Edge* e);

  /** Function to compress a graph. Needed if, for example, equivalence
      constraints are used to build a graoh structure with indices
      without gaps.  **/
  virtual void collapseEdge(Edge* e);

  /** Specifies the verbose level for debugging **/  
  int verboseLevel;
};

#endif



