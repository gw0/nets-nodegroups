/**
 * Node group structure
 *
 * Implementation of the *node group extraction framework* enables the exploration of node group structure of different networks, such as communities, modules, hubs and similar. Description of the algorithm can be found in:
 *
 * - L. Šubelj, N. Blagus, and M. Bajec, "Group extraction for real-world networks: The case of communities, modules, and hubs and spokes," in Proc. of NetSci '13, 2013, p. 152.
 * - L. Šubelj, S. Žitnik, N. Blagus, and M. Bajec, "Node mixing and group structure of complex software networks," Adv. Complex Syst., 2014.
 *
 * Author: gw0 [http://gw.tnode.com/] <gw.2014@tnode.com>
 * Version: 0.1
 */

#include "Snap.h"


/**
 * Group type parameter \tau(S,T)
 *
 * @param Graph Input graph
 * @param NodesS List of node IDs in subgraph S
 * @param NodesT List of node IDs in subgraph T
 * @return Value of group type parameter \tau
 */
double GroupTau(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT) {
  return (double)NodesS.IntrsLen(NodesT) / NodesS.UnionLen(NodesT);
}


/**
 * Count edges/links between subgraphs L(S,T)
 *
 * @param Graph Input graph
 * @param NodesS List of node IDs in subgraph S
 * @param NodesT List of node IDs in subgraph T
 * @param Lst Output number of links between S and T
 * @param Lsinvt Output number of links between S and inverse T
 * @return Number of links between S and T
 */
double GroupLinks(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT, double& Lst, double& Lsinvt) {
  Lst = 0.0;
  Lsinvt = 0.0;

  // iterate through S
  for (int s = 0; s < NodesS.Len(); ++s) {
    TUNGraph::TNodeI NI = Graph->GetNI(NodesS[s]);

    // iterate through its out-edges
    for (int e = 0; e < NI.GetOutDeg(); ++e) {
      // check if endpoint is inside T or not
      if (NodesT.IsIn(NI.GetOutNId(e))) {
        //printf("edge (%d %d)\n", NI.GetId(), NI.GetOutNId(e));
        Lst += 1.0;
      } else {
        Lsinvt += 1.0;
      }
    }
  }
  return Lst;
}
double GroupLinks(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT) {
  double Lst, Lsinvt;
  return GroupLinks(Graph, NodesS, NodesT, Lst, Lsinvt);
}


/**
 * Compute group criterion W(S,T)
 *
 * @param Graph Input graph
 * @param NodesS List of node IDs in subgraph S
 * @param NodesT List of node IDs in subgraph T
 * @param Wnorm Output W normalization factor
 * @param Lst Output number of links between S and T
 * @param Lsinvt Output number of links between S and inverse T
 * @return Computed value of group critetion W
 */
double GroupW(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT, double& Wnorm, double& Lst, double& Lsinvt) {
  int Len = Graph->GetEdges();
  int LenS = NodesS.Len();
  int LenT = NodesT.Len();

  // Normalization factor (geometric mean of s and t)
  double W2st = 2.0 * LenS * LenT;
  double Wnst = Len * (LenS + LenT);
  Wnorm = W2st * (Wnst - W2st) / Wnst;

  // Count edges/links between subgraphs L(S,T)
  GroupLinks(Graph, NodesS, NodesT, Lst, Lsinvt);

  // Group criterion
  return Wnorm * (Lst / LenT - Lsinvt / (Len - LenT)) / LenS;
}
double GroupW(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT) {
  double Wnorm, Lst, Lsinvt;
  return GroupW(Graph, NodesS, NodesT, Wnorm, Lst, Lsinvt);
}


/**
 * Fast estimate of group criterion W(S,T) after swapping nodes
 *
 * @param Graph Input graph
 * @param NodesS List of node IDs in subgraph S
 * @param NodesT List of node IDs in subgraph T
 * @param IfSwapS Estimate swap in S or T
 * @param DelNId Node ID to estimate deletion
 * @param AddNId Node ID to estimate addition
 * @param Wnorm Precomputed W normalization factor
 * @param Lst Output number of links between estimated S and T
 * @param Lsinvt Output number of links between estimated S and inverse T
 * @return Estimated value of group critetion W
 */
double GroupWFast(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT, bool IfSwapS, int DelNId, int AddNId, double& Wnorm, double& Lst, double& Lsinvt) {
  int Len = Graph->GetEdges();
  int LenS = NodesS.Len();
  int LenT = NodesT.Len();

  // Correction for removing node DelNId
  TIntV TmpNIdV;
  double LstTmp, LsinvtTmp;
  TmpNIdV.Clr();
  TmpNIdV.Add(DelNId);
  if (IfSwapS) {  // swap in S
    GroupLinks(Graph, TmpNIdV, NodesT, LstTmp, LsinvtTmp);
    Lst -= LstTmp;
    Lsinvt -= LsinvtTmp;

  } else {  // swap in T
    GroupLinks(Graph, TmpNIdV, NodesS, LstTmp, LsinvtTmp);
    Lst -= LstTmp;
    Lsinvt += LstTmp;
  }

  // Correction for adding node AddNId
  TmpNIdV.Clr();
  TmpNIdV.Add(AddNId);
  if (IfSwapS) {  // swap in S
    GroupLinks(Graph, TmpNIdV, NodesT, LstTmp, LsinvtTmp);
    Lst += LstTmp;
    Lsinvt += LsinvtTmp;

  } else {  // swap in T
    GroupLinks(Graph, TmpNIdV, NodesS, LstTmp, LsinvtTmp);
    Lst += LstTmp;
    Lsinvt -= LstTmp;      
  }

  // Corrected group criterion
  return Wnorm * (Lst / LenT - Lsinvt / (Len - LenT)) / LenS;
}


/**
 * Node group extraction algorithm
 *
 * @param Graph Input graph (edges will be removed)
 * @param LenS Size of subgraph S
 * @param LenT Size of subgraph T
 * @param Iters Number of optimization iterations
 * @param WBest Output best group criterion W
 * @return List of node IDs in extracted group
 */
TIntV GroupExtract(PUNGraph& Graph, int LenS, int LenT, int Iters, double& WBest) {
  // Initial random subgraph S and T
  TIntV NodesS;
  for (int i = 0; NodesS.Len() < LenS; ++i) {
    NodesS.AddMerged(Graph->GetRndNId());
  }
  TIntV NodesT;
  for (int i = 0; NodesT.Len() < LenT; ++i) {
    NodesT.AddMerged(Graph->GetRndNId());
  }

  // Initial group criterion
  double W, Wnorm, Lst, Lsinvt;
  WBest = GroupW(Graph, NodesS, NodesT, Wnorm, Lst, Lsinvt);
  TIntV NodesSBest = NodesS;
  TIntV NodesTBest = NodesT;

  // Inefficient random walk optimization
  /*for (int i = 0; i < Iters; ++i) {
    // Select nodes to swap in either S or T
    bool IfSwapS = TBool::GetRnd();
    TIntV& NodesX = (IfSwapS) ? NodesS : NodesT;
    int DelI = TUInt::GetRnd(NodesX.Len());
    int AddNId;
    do {
      AddNId = Graph->GetRndNId();
    } while(NodesX.IsIn(AddNId));

    // Swap and recompute group criterion
    int DelNId = NodesX[DelI];
    NodesX.Del(DelI);
    NodesX.AddMerged(AddNId);
    W = GroupW(Graph, NodesS, NodesT, Wnorm, Lst, Lsinvt);

    if (W > WBest) {
      //printf("%f %f %f %f\n", W, Wnorm, Lst, Lsinvt);
      WBest = W;
      NodesSBest = NodesS;
      NodesTBest = NodesT;
    //} else {  // gradient descent
    //  NodesX.DelAll(AddNId);
    //  NodesX.AddMerged(DelNId);
    }
  }*/

  // Greedy descent optimization with fast estimate of group criterion
  for (int i = 0; i < Iters; ++i) {
    // Select nodes to swap in either S or T
    bool IfSwapS = TBool::GetRnd();
    TIntV& NodesX = (IfSwapS) ? NodesS : NodesT;
    int DelI = TUInt::GetRnd(NodesX.Len());
    int AddNId;
    do {
      AddNId = Graph->GetRndNId();
    } while(NodesX.IsIn(AddNId));

    // Fast estimate of group criterion after swapping nodes
    double LstNew = Lst;
    double LsinvtNew = Lsinvt;
    W = GroupWFast(Graph, NodesS, NodesT, IfSwapS, NodesX[DelI], AddNId, Wnorm, LstNew, LsinvtNew);

    if (W > WBest) {
      NodesX.Del(DelI);
      NodesX.AddMerged(AddNId);
      Lst = LstNew;
      Lsinvt = LsinvtNew;

      //printf("%f %f %f %f\n", W, Wnorm, LstNew, LsinvtNew);
      WBest = W;
      NodesSBest = NodesS;
      NodesTBest = NodesT;
    }
  }

  printf("%f %d %d\n", WBest, NodesSBest.Len(), NodesTBest.Len());

  // Delete edges between S and T
  // iterate through S
  for (int s = 0; s < NodesSBest.Len(); ++s) {
    TUNGraph::TNodeI NI = Graph->GetNI(NodesSBest[s]);

    // iterate through its out-edges
    for (int e = 0; e < NI.GetOutDeg(); ++e) {
      // check if endpoint is inside T or not
      if (NodesTBest.IsIn(NI.GetOutNId(e))) {
        //printf("edge (%d %d)\n", NI.GetId(), NI.GetOutNId(e));
        Graph->DelEdge(NodesSBest[s], NI.GetOutNId(e));
      }
    }
  }

  return NodesSBest;
}


/**
 * Estimate average W in Erdos-Renyi random graph
 *
 * @param N Number of nodes in Erdos-Renyi random graph
 * @param M Number of edges in Erdos-Renyi random graph
 * @param LenS Size of subgraph S
 * @param LenT Size of subgraph T
 * @param Iters Number of iterations for average
 * @return Average value of group critetion W
 */
double RndWErdosRenyi(int N, int M, int LenS, int LenT, int Iters) {
  PUNGraph GraphER;
  TIntV NodesS, NodesT;
  double RndW = 0.0;

  // iterations for average
  for (int i = 0; i < Iters; ++i) {
    // Generate Erdos-Renyi random graph
    GraphER = TSnap::GenRndGnm<PUNGraph>(N, M, false);

    // Select random subgraphs S and T
    NodesS.Clr();
    for (int i = 0; NodesS.Len() < LenS; ++i) {
      NodesS.AddMerged(GraphER->GetRndNId());
    }
    NodesT.Clr();
    for (int i = 0; NodesT.Len() < LenT; ++i) {
      NodesT.AddMerged(GraphER->GetRndNId());
    }

    // Compute group criterion W
    RndW += GroupW(GraphER, NodesS, NodesT);
  }

  return RndW / Iters;
}


/**
 * Console application entry point
 */
int main(int argc, char* argv[]) {
  // Header
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("Node group structure. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try

  // Parameters
  const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "graph.edgelist", "Input graph (list of undirected edges)");
  const TStr LabelFNm = Env.GetIfArgPrefixStr("-l:", "graph.labels", "Optional input node labels (node ID, node label)");
  const TStr OutFNm = Env.GetIfArgPrefixStr("-o:", "graph.groups", "Output group assignments");

  // Input
  PUNGraph Graph = TSnap::LoadEdgeList<PUNGraph>(InFNm, false);
  //PUNGraph Graph = TSnap::GenRndGnm<PUNGraph>(5000, 10000); // generate a random graph
  TIntStrH NIDNameH;
  if (LabelFNm.Len() > 0) {
    TSsParser Ss(LabelFNm, ssfTabSep);
    while (Ss.Next()) {
      if (Ss.Len() > 0) {
        NIDNameH.AddDat(Ss.GetInt(0), Ss.GetFld(1));
      }
    }
  }

  // Run node group extraction framework
  double RndW = RndWErdosRenyi(Graph->GetNodes(), Graph->GetEdges(), 20, 20, 1000);
  double W;
  TCnComV GroupsV;
  do {
    GroupsV.Add(TCnCom(GroupExtract(Graph, 20, 20, 100000, W)));
  } while(W > RndW);

  // Output
  FILE *F = fopen(OutFNm.CStr(), "wt");
  fprintf(F, "# Input: %s\n", InFNm.CStr());
  fprintf(F, "# Nodes: %d    Edges: %d\n", Graph->GetNodes(), Graph->GetEdges());
  fprintf(F, "# Groups: %d\n", GroupsV.Len());
  fprintf(F, "# NId\tGroupId\tNLabel\n");
  for (int c = 0; c < GroupsV.Len(); c++) {
    for (int i = 0; i < GroupsV[c].Len(); i++) {
      int NId = GroupsV[c][i].Val;
      if (NIDNameH.IsKey(NId)) {
        fprintf(F, "%d\t%d\t%s\n", NId, c, NIDNameH.GetDat(NId).CStr());
      } else {
        fprintf(F, "%d\t%d\t-\n", NId, c);
      }
    }
  }
  fclose(F);

  // Footer
  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
