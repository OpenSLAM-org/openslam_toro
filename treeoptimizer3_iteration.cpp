#include "treeoptimizer3.hh"
#include <fstream>
#include <string>

using namespace std;

#define DEBUG(i) \
        if (verboseLevel>i) cerr


struct NodeInfo{
  TreeOptimizer3::Vertex* n;
  double translationalWeight;
  double rotationalWeight;
  int direction;
  NodeInfo(TreeOptimizer3::Vertex* v=0, double tw=0, double rw=0, int dir=0){
    n=v;
    translationalWeight=tw;
    rotationalWeight=rw;
    direction=dir;
  }
};

typedef std::vector<NodeInfo> NodeInfoVector;

/********************************** unpreconditioned error distribution ************************************/

void TreeOptimizer3::propagateErrors(){
  iteration++;
  int edgeCount=0;
  // this is the workspace for computing the paths without 
  // bothering too much the memory allocation
  static NodeInfoVector path;
  path.resize(edges.size()+1);
  static Rotation zero(0.,0.,0.);

  onIterationStart(iteration);
  for (EdgeSet::iterator it=sortedEdges->begin(); it!=sortedEdges->end(); it++){
    edgeCount++;
    if (! (edgeCount%1000))
      DEBUG(1) << "c";
    
    if (isDone())
      return;

    Edge* e=*it;
    Vertex* top=e->top;
    Vertex* v1=e->v1;
    Vertex* v2=e->v2;
    int l=e->length;

    onStepStart(*it);
    
    
    DEBUG(2) << "Edge: " << v1->id << " " << v2->id << ", top=" << top->id << ", length="<< l <<endl;

    //BEGIN: Path computation
    int pc=0;
    Vertex* aux=v1;
    while(aux!=top){
      path[pc++]=NodeInfo(aux,0.,0.,-1);
      aux=aux->parent;
    }
    int topIndex=pc;
    path[pc++]=NodeInfo(top,0.,0.,0);
    pc=l;
    aux=v2;
    while(aux!=top){
      path[pc--]=NodeInfo(aux,0.,0.,1);
      aux=aux->parent;
    }
    //END: Path computation
    

    Transformation topTransformation=top->transformation;
    Transformation topParameters=top->parameters;

    //BEGIN: Rotational Error
    Rotation r1=getRotation(v1, top);
    Rotation r2=getRotation(v2, top);
    Rotation re=e->transformation.rotation();
    Rotation rR=r2.inverse()*(r1*re);
    Translation angles=rR.toAngles();
    rR=Rotation(angles.roll(), angles.pitch(), angles.yaw());
    rR=rR.normalized();
    
    double rotationFactor=sqrt(double(l))*rotGain/(double)iteration;
    
    if (rotationFactor>1)
      rotationFactor=1;

    Rotation irR=rR.inverse();
    irR=irR.normalized();

    Translation axis=rR.axis();
    double angle=rR.angle();
    //cerr <<  "EE(" << v1->id << "," << v2->id <<")" << axis.x() << " " << axis.y() << " " << axis.z() << " alpha=" << angle << endl;
    Translation iAxis=irR.axis();
    double iAngle=irR.angle();
    
    double cw=0;
    Rotation Gamma_prev=zero;
    for (int i= topIndex-1; i>=0; i--){
      cw+=1./l;
      Vertex* n=path[i].n;
      Rotation R=n->parameters.rotation();
      Rotation B(iAxis, iAngle*cw*rotationFactor);
      R= Gamma_prev.inverse() * R * B;
      Gamma_prev=B;
      n->parameters.setRotation(R);
      n->transformation=n->parent->transformation*n->parameters;
    }
  
    cw=0;
    Gamma_prev=zero;
    for (int i= topIndex+1; i<=l; i++){
      cw+=1./l;
      Vertex* n=path[i].n;
      Rotation R=n->parameters.rotation();
      Rotation B(axis, angle*cw*rotationFactor);
      R= Gamma_prev.inverse() * R * B;
      Gamma_prev=B;
      n->parameters.setRotation(R);
      n->transformation=n->parent->transformation*n->parameters;
    }

    //END: Rotational Error

    //BEGIN: Translational Error
    Translation topTranslation=top->transformation.translation();
   
    Transformation tr12=v1->transformation*e->transformation;
    Translation tR=tr12.translation()-v2->transformation.translation();

    double translationFactor=trasGain*l/(double)iteration;
    if (translationFactor>1)
      translationFactor=1;
    Translation dt=tR*translationFactor;

    //left wing
    double lcum=0;
    for (int i=topIndex-1; i>=0; i--){
      Vertex* n=path[i].n;
      lcum-=1./l;
      double fraction=lcum;
      Translation offset= dt*fraction;
      Translation T=n->transformation.translation()+offset;
      n->transformation.setTranslation(T);
    }
    //right wing
    
    double rcum=0;
    for (int i=topIndex+1; i<=l; i++){
      Vertex* n=path[i].n;
      rcum+=1./l;
      double fraction=rcum;
      Translation offset= dt*fraction;
      Translation T=n->transformation.translation()+offset;
      n->transformation.setTranslation(T);
    }
    assert(fabs(lcum+rcum)-1<1e-6);

    recomputeParameters(v1, top);
    recomputeParameters(v2, top);
    //END: Translational Error

    onStepFinished(e);

    if (verboseLevel>2){
      Rotation newRotResidual=v2->transformation.rotation().inverse()*(v1->transformation.rotation()*re);
      Translation newRotResidualAxis=newRotResidual.axis();
      double      newRotResidualAngle=newRotResidual.angle();
      Translation rotResidualAxis=rR.axis();
      double      rotResidualAngle=rR.angle();
      Translation newTransResidual=(v1->transformation*e->transformation).translation()-v2->transformation.translation();

      cerr << "RotationalFraction:  " << rotationFactor << endl;
      cerr << "Rotational residual: " 
	   << " axis " << rotResidualAxis.x()  << "\t" << rotResidualAxis.y()  << "\t" << rotResidualAxis.z() << " --> "
	   << "   ->  " << newRotResidualAxis.x()  << "\t" << newRotResidualAxis.y()  << "\t" << newRotResidualAxis.z()  << endl;
      cerr << " angle " << rotResidualAngle  << "\t" << newRotResidualAngle  << endl;
      
      cerr << "Translational Fraction:  " << translationFactor << endl;
      cerr << "Translational Residual" << endl;
      cerr << "    " << tR.x()  << "\t" << tR.y()  << "\t" << tR.z()  << endl;
      cerr << "    " << newTransResidual.x()  << "\t" << newTransResidual.y()  << "\t" << newTransResidual.z()  << endl;
    }
    
    if (verboseLevel>101){
      char filename [1000];
      sprintf(filename, "po-%02d-%03d-%03d-.dat", iteration, v1->id, v2->id);
      recomputeAllTransformations();
      saveGnuplot(filename);
    }
  }
  onIterationFinished(iteration);
}


/********************************** Preconditioned error distribution ************************************/


inline double max3( const double& a, const double& b, const double& c){
  double m=a>b?a:b;
  return m>c?m:c;
}
inline double min3( const double& a, const double& b, const double& c){
  double m=a<b?a:b;
  return m<c?m:c;
}



void TreeOptimizer3::computePreconditioner(){
  for (uint i=0; i<M.size(); i++){
    M[i][0]=0;
    M[i][1]=0;
  }
  gamma[0] = gamma[1] = MAXDOUBLE;

  int edgeCount=0;
  for (EdgeSet::iterator it=sortedEdges->begin(); it!=sortedEdges->end(); it++){
    edgeCount++;
    if (! (edgeCount%1000))
      DEBUG(1) << "m";

    Edge* e=*it;
    Transformation t=e->transformation;
    InformationMatrix W=e->informationMatrix;
    
    Vertex* top=e->top;
    for (int dir=0; dir<2; dir++){
      Vertex* n = (dir==0)? e->v1 : e->v2;
      while (n!=top){
	uint i=n->id;
	double rW=min3(W[0][0], W[1][1], W[2][2]);
	double tW=min3(W[3][3], W[4][4], W[5][5]);
	M[i][0]+=rW;
	M[i][1]+=tW;
	gamma[0]=gamma[0]<rW?gamma[0]:rW;
	gamma[1]=gamma[1]<tW?gamma[1]:tW;
	n=n->parent;
      }
    }

  }
  
  if (verboseLevel>1){
    for (uint i=0; i<M.size(); i++){
      cerr << "M[" << i << "]=" << M[i][0] << " " << M[i][1] << endl;
    }
  }
}


void TreeOptimizer3::propagateErrorsPreconditioner(){
  iteration++;
  int edgeCount=0;
  // this is the workspace for computing the paths without 
  // bothering too much the memory allocation
  static NodeInfoVector path;
  path.resize(edges.size()+1);

  static Rotation zero(0.,0.,0.);
  onIterationStart(iteration);
  for (EdgeSet::iterator it=sortedEdges->begin(); it!=sortedEdges->end(); it++){
    edgeCount++;
    if (! (edgeCount%1000))
      DEBUG(1) << "c";

    if (isDone())
      return;
    
    Edge* e=*it;
    Vertex* top=e->top;
    Vertex* v1=e->v1;
    Vertex* v2=e->v2;
    int l=e->length;
    onStepStart(e);
    
    
    DEBUG(2) << "Edge: " << v1->id << " " << v2->id << ", top=" << top->id << ", length="<< l <<endl;

    //BEGIN: Path and weight computation 
    int pc=0;
    Vertex* aux=v1;
    double totTW=0, totRW=0;
    while(aux!=top){
      int index=aux->id;
      double tw=1./M[index][1], rw=1./M[index][1];
      totTW+=tw;
      totRW+=rw;
      path[pc++]=NodeInfo(aux,tw,rw,-1);
      aux=aux->parent;
    }
    int topIndex=pc;
    path[pc++]=NodeInfo(top,0.,0.,0);
    pc=l;
    aux=v2;
    while(aux!=top){
      int index=aux->id;
      double tw=1./M[index][1], rw=1./M[index][1];
      totTW+=tw;
      totRW+=rw;
      path[pc--]=NodeInfo(aux,tw,rw,1);
      aux=aux->parent;
    }
    //END: Path and weight computation
    

    Transformation topTransformation=top->transformation;
    Transformation topParameters=top->parameters;

    //BEGIN: Rotational Error
    Rotation r1=getRotation(v1, top);
    Rotation r2=getRotation(v2, top);
    Rotation re=e->transformation.rotation();
    Rotation rR=r2.inverse()*(r1*re);
    Translation angles=rR.toAngles();
    rR=Rotation(angles.roll(), angles.pitch(), angles.yaw());
    rR=rR.normalized();
    
    double rotationFactor=
      sqrt(double(l))*rotGain/
      ( gamma[0]* (double)iteration * min3(e->informationMatrix[0][0],
					   e->informationMatrix[1][1], 
					   e->informationMatrix[2][2]));
    if (rotationFactor>1)
      rotationFactor=1;

    Rotation irR=rR.inverse();
    irR=irR.normalized();

    Translation axis=rR.axis();
    double angle=rR.angle();
    Translation iAxis=irR.axis();
    double iAngle=irR.angle();
    
    double cw=0;
    Rotation Gamma_prev=zero;
    for (int i= topIndex-1; i>=0; i--){
      Vertex* n=path[i].n;
      cw+=path[i].rotationalWeight/totRW;
      Rotation R=n->parameters.rotation();
      Rotation B(iAxis, iAngle*cw*rotationFactor);
      R= Gamma_prev.inverse() * R * B;
      Gamma_prev=B;
      n->parameters.setRotation(R);
      n->transformation=n->parent->transformation*n->parameters;
    }
  
    cw=0;
    Gamma_prev=zero;
    for (int i= topIndex+1; i<=l; i++){
      cw+=path[i].rotationalWeight/totRW;
      Vertex* n=path[i].n;
      Rotation R=n->parameters.rotation();
      Rotation B(axis, angle*cw*rotationFactor);
      R= Gamma_prev.inverse() * R * B;
      Gamma_prev=B;
      n->parameters.setRotation(R);
      n->transformation=n->parent->transformation*n->parameters;
    }

    //END: Rotational Error

    //BEGIN: Translational Error
    Translation topTranslation=top->transformation.translation();
   
    Transformation tr12=v1->transformation*e->transformation;
    Translation tR=tr12.translation()-v2->transformation.translation();

    double translationFactor=trasGain*l/
	( gamma[1]* (double)iteration * min3(e->informationMatrix[3][3],
					     e->informationMatrix[4][4], 
					     e->informationMatrix[5][5]));
    if (translationFactor>1)
      translationFactor=1;
    Translation dt=tR*translationFactor;

    //left wing
    double lcum=0;
    for (int i=topIndex-1; i>=0; i--){
      Vertex* n=path[i].n;
      lcum-=path[i].translationalWeight/totTW;
      double fraction=lcum;
      Translation offset= dt*fraction;
      Translation T=n->transformation.translation()+offset;
      n->transformation.setTranslation(T);
    }
    //right wing
    
    double rcum=0;
    for (int i=topIndex+1; i<=l; i++){
      Vertex* n=path[i].n;
      rcum+=path[i].translationalWeight/totTW;;
      double fraction=rcum;
      Translation offset= dt*fraction;
      Translation T=n->transformation.translation()+offset;
      n->transformation.setTranslation(T);
    }
    assert(fabs(lcum+rcum)-1<1e-6);

    recomputeParameters(v1, top);
    recomputeParameters(v2, top);


    //END: Translational Error
    onStepFinished(e);

    if (verboseLevel>2){
      Rotation newRotResidual=v2->transformation.rotation().inverse()*(v1->transformation.rotation()*re);
      Translation newRotResidualAxis=newRotResidual.axis();
      double      newRotResidualAngle=newRotResidual.angle();
      Translation rotResidualAxis=rR.axis();
      double      rotResidualAngle=rR.angle();
      Translation newTransResidual=(v1->transformation*e->transformation).translation()-v2->transformation.translation();

      cerr << "RotationalFraction:  " << rotationFactor << endl;
      cerr << "Rotational residual: " 
	   << " axis " << rotResidualAxis.x()  << "\t" << rotResidualAxis.y()  << "\t" << rotResidualAxis.z() << " --> "
	   << "   ->  " << newRotResidualAxis.x()  << "\t" << newRotResidualAxis.y()  << "\t" << newRotResidualAxis.z()  << endl;
      cerr << " angle " << rotResidualAngle  << "\t" << newRotResidualAngle  << endl;
      
      cerr << "Translational Fraction:  " << translationFactor << endl;
      cerr << "Translational Residual" << endl;
      cerr << "    " << tR.x()  << "\t" << tR.y()  << "\t" << tR.z()  << endl;
      cerr << "    " << newTransResidual.x()  << "\t" << newTransResidual.y()  << "\t" << newTransResidual.z()  << endl;
    }
    
    if (verboseLevel>101){
      char filename [1000];
      sprintf(filename, "po-%02d-%03d-%03d-.dat", iteration, v1->id, v2->id);
      recomputeAllTransformations();
      saveGnuplot(filename);
    }
  }
  onIterationFinished(iteration);
}
