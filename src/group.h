/**
 * nodegroups - Group structures header file
 *
 * @author  gw0 [http://gw.tnode.com/] <gw.2014@tnode.com>
 * @version 0.2
 */

#ifndef group_h
#define group_h
#define group_h_VERSION 0.2

#include "Snap.h"

/////////////////////////////////////////////////

#define DEF_OptRestarts 1000  /* 1000 */
#define DEF_OptMxSteps 100000  /* 100000 */
#define DEF_OptStopSteps 1000  /* 1000, disable: DEF_OptMxSteps */
#define DEF_OptInitSample 0  /* random: 0, add: >0, del: <0 */
#define DEF_RndRestarts 100  /* 100 */
#define DEF_RndRecompW 1.1  /* 1.1, disable: INFINITY */
#define DEF_RndStopW 1.01  /* 1.01 */

/** Results of node group extraction (into S, T) */
class TGroupST {
public:
  int N;  // number of nodes in graph
  int M;  // number of edges in graph

  int SubSN;  // number of nodes in subgraph on S
  int SubSM;  // number of edges in subgraph on S
  TIntV SubSNIdV;  // list of node IDs in group S

  int SubTN;  // number of nodes in subgraph on T
  int SubTM;  // number of edges in subgraph on T
  TIntV SubTNIdV;  // list of node IDs in group T

  int SubSTN;  // number of nodes in subgraph on S and T
  int SubSTM;  // number of edges in subgraph on S and T
  int LinksST;  // number of edges L(S,T) between groups S and T
  int LinksSTc;  // number of edges L(S,-T) between groups S and complement T

  double Tau;  // group type parameter Tau(S,T)
  double W;  // group critetion W(S,T)
  double ModularityS;  // modularity measure on group S
  double ModularityT;  // modularity measure on group T

  TStr GetStr(bool Verbose=true);
  void RecomputeAll(const PUNGraph& Graph, const TIntV& NewSubSNIdV, const TIntV& NewSubTNIdV);
};
typedef TVec<TGroupST> TGroupSTV;

/////////////////////////////////////////////////

double LinksCnt(int& LinksST, int& LinksSTc, const PUNGraph& Graph, const TIntV& SubSNIdV, const TIntV& SubTNIdV, bool DoDelEdges=false);
double LinksCnt(const PUNGraph& Graph, const TIntV& SubSNIdV, const TIntV& SubTNIdV, bool DoDelEdges=false);
double GroupTau(const PUNGraph& Graph, const TIntV& SubSNIdV, const TIntV& SubTNIdV);
double GroupW(int N, int SubSN, int SubTN, int LinksST, int LinksSTc);
double GroupW(TGroupST& G, const PUNGraph& Graph, const TIntV& SubSNIdV, const TIntV& SubTNIdV);
double GroupWFast(TGroupST& G, const PUNGraph& Graph, const TIntV& SubSNIdV, const TIntV& SubTNIdV, const TIntV& AddSNIdV, const TIntV& DelSNIdV, const TIntV& AddTNIdV, const TIntV& DelTNIdV, int LinksST, int LinksSTc);
double GroupExtractSingle(TGroupST& GBest, const PUNGraph& Graph, int OptMxSteps=DEF_OptMxSteps, int OptStopSteps=DEF_OptStopSteps, int OptInitSample=DEF_OptInitSample);
double GroupExtractRestarter(TGroupST& GBest, const PUNGraph& Graph, int OptRestarts=DEF_OptRestarts, int OptMxSteps=DEF_OptMxSteps, int OptStopSteps=DEF_OptStopSteps, int OptInitSample=DEF_OptInitSample);
double GroupExtractAvgRndGnm(TGroupST& RAvg, int N, int M, int RndRestarts=DEF_RndRestarts, int OptMxSteps=DEF_OptMxSteps, int OptStopSteps=DEF_OptStopSteps, int OptInitSample=DEF_OptInitSample);
int GroupExtractFramework(TGroupSTV& GroupV, PUNGraph& Graph, int OptRestarts=DEF_OptRestarts, int OptMxSteps=DEF_OptMxSteps, int OptStopSteps=DEF_OptStopSteps, int OptInitSample=DEF_OptInitSample, int RndRestarts=DEF_RndRestarts, double RndRecompW=DEF_RndRecompW, double RndStopW=DEF_RndStopW);

/////////////////////////////////////////////////

#endif