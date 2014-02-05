/**
 * nodegroups - TGroups header file
 *
 * @author  gw0 [http://gw.tnode.com/] <gw.2014@tnode.com>
 * @version 0.1
 */

#ifndef tgroups_h
#define tgroups_h

#include "Snap.h"

/////////////////////////////////////////////////

/** Results of node group extraction */
struct TGroup {
  double W;
  int N;
  int M;
  int NodesSLen;
  int NodesTLen;
  double LinksST;
  double Tau;
  double ModularityS;
  TIntV NodesS;
  TIntV NodesT;
};
typedef TVec<TGroup> TGroupV;

/////////////////////////////////////////////////

double GroupTau(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT);
double GroupLinks(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT, bool IfDelEdge, double& LinksST, double& LinksSInvT);
double GroupLinks(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT, bool IfDelEdge=false);
double GroupW(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT, double& WNorm, double& LinksST, double& LinksSInvT);
double GroupW(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT);
double GroupWFast(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT, bool IfSwapS, int SwapNId, double& WNorm, double& LinksST, double& LinksSInvT);
double GroupExtract(PUNGraph& Graph, int Steps, double& W, TIntV& NodesS, TIntV& NodesT);
double GroupExtractRerunner(PUNGraph& Graph, int Iters, int Steps, double& WBest, TIntV& NodesSBest, TIntV& NodesTBest);

/////////////////////////////////////////////////

#endif