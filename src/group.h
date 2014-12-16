/**
 * nets-nodegroups - Group structure functions header file
 *
 * @author  gw0 [http://gw.tnode.com/] <gw.2014@tnode.com>
 * @version 0.4
 */

#ifndef group_h
#define group_h
#define group_h_VERSION 0.4

#include "Snap.h"

/////////////////////////////////////////////////

#define DEF_OptRestarts 2000  /* 2000 */
#define DEF_OptMxSteps 100000  /* 100000 */
#define DEF_OptStopSteps 1000  /* 1000, disable: DEF_OptMxSteps */
#define DEF_OptInitSample 1  /* random: 0, specified: >0, all but specified: <0 */
#define DEF_FinishCnt 0  /* 0, infinity: 0 */
#define DEF_FinishRndW 1.0  /* 1.0 for 1% */
#define DEF_RndGraphs 500  /* 500 */
#define DEF_RndRestarts (DEF_OptRestarts / 200)  /* 10 */
#define DEF_RndRecompW INFINITY  /* 1.1, disable: INFINITY */

/** ST-group extraction results */
class TGroupST {
public:
  int N;  // N - number of nodes left in graph
  int M;  // M - number of edges left in graph

  int SubSN;  // N_S - number of nodes in subgraph on group *S*
  int SubSM;  // M_S - number of edges in subgraph on group *S* (not used)
  TIntV SubSNIdV;  // list of node IDs in group S

  int SubTN;  // N_T - number of nodes in subgraph on linking pattern *T*
  int SubTM;  // M_T - number of edges in subgraph on linking pattern *T* (not used)
  TIntV SubTNIdV;  // list of node IDs in linking pattern T

  int SubSTN;  // N_ST - number of nodes in subgraph on intersection of *S* and *T* (not used)
  int SubSTM;  // M_ST - number of edges in subgraph on intersection of *S* and *T* (not used)
  int LinksST;  // L_ST - number of edges *L(S,T)* between groups *S* and *T*
  int LinksSTc;  // L_STc - number of edges *L(S,Tc)* between groups *S* and complement of *T*

  double Tau;  // Tau - group type parameter *Tau(S,T)*
  double W;  // W - group critetion *W(S,T)*
  double ModularityS;  // Mod_S - modularity measure on group *S* (not used)
  double ModularityT;  // Mod_T - modularity measure on group *T* (not used)

  bool operator== (const TGroupST& B) const { return this->W == B.W; };
  bool operator< (const TGroupST& B) const { return this->W < B.W; };

  TStr GetStr(int Type=20);
  void RecomputeAll(const PUNGraph& Graph, const TIntV& NewSubSNIdV, const TIntV& NewSubTNIdV);
};
typedef TVec<TGroupST> TGroupSTV;

/////////////////////////////////////////////////

double LinksCnt(int& LinksST, int& LinksSTc, const PUNGraph& Graph, const TIntV& SubSNIdV, const TIntV& SubTNIdV, bool DoDelEdges=false);
double LinksCnt(const PUNGraph& Graph, const TIntV& SubSNIdV, const TIntV& SubTNIdV, bool DoDelEdges=false);
double GroupTau(const TIntV& SubSNIdV, const TIntV& SubTNIdV);
TStr GroupName(const TIntV& SubSNIdV, const TIntV& SubTNIdV);
double GroupW(int N, int SubSN, int SubTN, int LinksST, int LinksSTc);
double GroupW(TGroupST& G, const PUNGraph& Graph, const TIntV& SubSNIdV, const TIntV& SubTNIdV);
double GroupWFast(TGroupST& G, const PUNGraph& Graph, const TIntV& SubSNIdV, const TIntV& SubTNIdV, const TIntV& AddSNIdV, const TIntV& DelSNIdV, const TIntV& AddTNIdV, const TIntV& DelTNIdV, int LinksST, int LinksSTc);
double GroupExtractSingle(TGroupST& GBest, const PUNGraph& Graph, int OptMxSteps=DEF_OptMxSteps, int OptStopSteps=DEF_OptStopSteps, int OptInitSample=DEF_OptInitSample);
double GroupExtractRestarter(TGroupST& GBest, const PUNGraph& Graph, int OptRestarts=DEF_OptRestarts, int OptMxSteps=DEF_OptMxSteps, int OptStopSteps=DEF_OptStopSteps, int OptInitSample=DEF_OptInitSample);
int GroupExtractRndGnms(TGroupSTV& GroupERV, int N, int M, int RndGraphs=DEF_RndGraphs, int RndRestarts=DEF_RndRestarts, int OptMxSteps=DEF_OptMxSteps, int OptStopSteps=DEF_OptStopSteps, int OptInitSample=DEF_OptInitSample);
int GroupExtractFramework(TGroupSTV& GroupV, PUNGraph& Graph, int OptRestarts=DEF_OptRestarts, int OptMxSteps=DEF_OptMxSteps, int OptStopSteps=DEF_OptStopSteps, int OptInitSample=DEF_OptInitSample, int FinishCnt=DEF_FinishCnt, double FinishRndW=DEF_FinishRndW, int RndGraphs=DEF_RndGraphs, int RndRestarts=DEF_RndRestarts, double RndRecompW=DEF_RndRecompW);

/////////////////////////////////////////////////

#endif
