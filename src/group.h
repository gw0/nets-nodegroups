/**
 * nodegroups - Group structures header file
 *
 * @author  gw0 [http://gw.tnode.com/] <gw.2014@tnode.com>
 * @version 0.22
 */

#ifndef group_h
#define group_h
#define group_h_VERSION 0.22

#include "Snap.h"

/////////////////////////////////////////////////

#define DEF_OptRestarts 2000  /* 2000 */
#define DEF_OptMxSteps 100000  /* 100000 */
#define DEF_OptStopSteps 1000  /* 1000, disable: DEF_OptMxSteps */
#define DEF_OptInitSample 1  /* random: 0, specified: >0, all but specified: <0 */
#define DEF_RndGraphs (DEF_OptRestarts / 10)  /* 200 */
#define DEF_RndRestarts 1  /* 1 */
#define DEF_RndRecompW INFINITY  /* 1.1, disable: INFINITY */
#define DEF_RndStopW 1.1  /* 1.1 */

/** Results of node group extraction (into S, T) */
class TGroupST {
public:
  int N;  // N - number of nodes in graph
  int M;  // M - number of edges in graph

  int SubSN;  // N_S - number of nodes in subgraph on S
  int SubSM;  // M_S - number of edges in subgraph on S (not used)
  TIntV SubSNIdV;  // list of node IDs in group S

  int SubTN;  // N_T - number of nodes in subgraph on T
  int SubTM;  // M_T - number of edges in subgraph on T (not used)
  TIntV SubTNIdV;  // list of node IDs in group T

  int SubSTN;  // N_ST - number of nodes in subgraph on S and T intersection (not used)
  int SubSTM;  // M_ST - number of edges in subgraph on S and T intersection (not used)
  int LinksST;  // L_ST - number of edges L(S,T) between groups S and T
  int LinksSTc;  // L_STc - number of edges L(S,-T) between groups S and complement T

  double Tau;  // Tau - group type parameter Tau(S,T)
  double W;  // W - group critetion W(S,T)
  double ModularityS;  // Mod_S - modularity measure on group S (not used)
  double ModularityT;  // Mod_T - modularity measure on group T (not used)

  bool operator== (const TGroupST& B) const { return this->W == B.W; };
  bool operator< (const TGroupST& B) const { return this->W < B.W; };

  TStr GetStr(int Type=10);
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
double GroupExtractAvgRndGnm(TGroupST& RAvg, int N, int M, int RndGraphs=DEF_RndGraphs, int RndRestarts=DEF_RndRestarts, int OptMxSteps=DEF_OptMxSteps, int OptStopSteps=DEF_OptStopSteps, int OptInitSample=DEF_OptInitSample);
int GroupExtractFramework(TGroupSTV& GroupV, PUNGraph& Graph, int OptRestarts=DEF_OptRestarts, int OptMxSteps=DEF_OptMxSteps, int OptStopSteps=DEF_OptStopSteps, int OptInitSample=DEF_OptInitSample, int RndGraphs=DEF_RndGraphs, int RndRestarts=DEF_RndRestarts, double RndRecompW=DEF_RndRecompW, double RndStopW=DEF_RndStopW);

/////////////////////////////////////////////////

#endif