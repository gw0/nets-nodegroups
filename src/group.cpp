/**
 * nodegroups - Group structures functions
 *
 * @author  gw0 [http://gw.tnode.com/] <gw.2014@tnode.com>
 */

#include "group.h"


/**
 * Convert to string
 *
 * @param this Input TGroupST object
 * @param Verbose Switch output verbose level
 */
TStr TGroupST::GetStr(bool Verbose/*=true*/) {
  if (Verbose) {
    return TStr::Fmt("N=%-5d M=%-5d N_S=%-5d M_S=%-5d N_T=%-5d M_T=%-5d N_ST=%-5d M_ST=%-5d L_S,T=%-5d L_S,-T=%-5d W=%-12.6f Tau=%-8.6f Mod_S=%-9.6f Mod_T=%-9.6f", N, M, SubSN, SubSM, SubTN, SubTM, SubSTN, SubSTM, LinksST, LinksSTInv, W, Tau, ModularityS, ModularityT);

  } else {
    return TStr::Fmt("W=%-12.6f N=%-5d M=%-5d |S|=%-5d |T|=%-5d L(S,T)=%-5d Tau=%-8.6f Mod(S)=%-9.6f", W, N, M, SubSN, SubTN, LinksST, Tau, ModularityS);
  }
}


/**
 * Recompute all results of node group extraction (except W)
 *
 * @param[out] this Output TGroupST object
 * @param Graph Input graph
 * @param NewSubSNIdV List of node IDs in group S
 * @param NewSubTNIdV List of node IDs in group T
 */
void TGroupST::RecomputeAll(PUNGraph& Graph, TIntV& NewSubSNIdV, TIntV& NewSubTNIdV) {
  N = Graph->GetNodes();
  M = Graph->GetEdges();

  SubSNIdV = NewSubSNIdV;
  PUNGraph SubS = TSnap::GetSubGraph(Graph, SubSNIdV);
  SubSN = SubS->GetNodes();  // or SubSNIdV.Len()
  SubSM = SubS->GetEdges();

  SubTNIdV = NewSubTNIdV;
  PUNGraph SubT = TSnap::GetSubGraph(Graph, SubTNIdV);
  SubTN = SubT->GetNodes();  // or SubTNIdV.Len()
  SubTM = SubT->GetEdges();

  TIntV SubSTNIdV;
  SubSNIdV.Intrs(SubTNIdV, SubSTNIdV);
  PUNGraph SubST = TSnap::GetSubGraph(Graph, SubSTNIdV);
  SubSTN = SubST->GetNodes();  // or SubSNIdV.IntrsLen(SubTNIdV)
  SubSTM = SubST->GetEdges();
  LinksCnt(LinksST, LinksSTInv, Graph, SubSNIdV, SubTNIdV, false);

  Tau = GroupTau(Graph, SubSNIdV, SubTNIdV);
  W = GroupW(N, SubSN, SubTN, LinksST, LinksSTInv);
  ModularityS = TSnap::GetModularity(Graph, SubSNIdV, N);
  ModularityT = TSnap::GetModularity(Graph, SubTNIdV, N);
}


/**
 * Count edges L(S,T) between groups S and T
 *
 * @param[out] LinksST Number of edges L(S,T) between groups S and T
 * @param[out] LinksSTInv Number of edges L(S,T) between groups S and inverse T
 * @param Graph Input graph
 * @param SubSNIdV List of node IDs in group S
 * @param SubTNIdV List of node IDs in group T
 * @param DoDelEdges Delete found edges between S and T
 * @return Number of edges L(S,T) between groups S and T
 */
double LinksCnt(int& LinksST, int& LinksSTInv, PUNGraph& Graph, TIntV& SubSNIdV, TIntV& SubTNIdV, bool DoDelEdges/*=false*/) {
  LinksST = 0;
  LinksSTInv = 0;

  // iterate through S
  for (int s = 0; s < SubSNIdV.Len(); ++s) {
    TUNGraph::TNodeI NI = Graph->GetNI(SubSNIdV[s]);

    // iterate through its out-edges
    for (int e = 0; e < NI.GetOutDeg(); ++e) {
      if (SubTNIdV.IsIn(NI.GetOutNId(e))) {  // endpoint is inside T
        //printf("edge (%d %d)\n", NI.GetId(), NI.GetOutNId(e));
        LinksST += 1;
        if (DoDelEdges) {  // delete edge between S and T
          Graph->DelEdge(SubSNIdV[s], NI.GetOutNId(e));
        }

      } else {  // endpoint is not inside T
        LinksSTInv += 1;
      }
    }
  }

  return LinksST;
}
double LinksCnt(PUNGraph& Graph, TIntV& SubSNIdV, TIntV& SubTNIdV, bool DoDelEdges/*=false*/) {
  int LinksST, LinksSTInv;
  return LinksCnt(LinksST, LinksSTInv, Graph, SubSNIdV, SubTNIdV, DoDelEdges);
}


/**
 * Group type parameter Tau(S,T)
 *
 * @param Graph Input graph
 * @param SubSNIdV List of node IDs in group S
 * @param SubTNIdV List of node IDs in group T
 * @return Value of group type parameter Tau
 */
double GroupTau(PUNGraph& Graph, TIntV& SubSNIdV, TIntV& SubTNIdV) {
  return (double)SubSNIdV.IntrsLen(SubTNIdV) / SubSNIdV.UnionLen(SubTNIdV);
}


/**
 * Compute group criterion W(S,T)
 *
 * @param N Number of nodes in graph
 * @param SubSN Number of nodes in subgraph on S
 * @param SubSNIdV List of node IDs in group S
 * @param SubTN Number of nodes in subgraph on T
 * @param SubTNIdV List of node IDs in group T
 * @param LinksST Number of edges L(S,T) between groups S and T
 * @param LinksSTInv Number of edges L(S,T) between groups S and inverse T
 * @param Graph Input graph
 * @param[out] G Output results of node group extraction (into S, T)
 * @return Value of group critetion W
 */
double GroupW(int N, int SubSN, int SubTN, int LinksST, int LinksSTInv) {
  // Normalization factor (geometric mean of |S| and |T|)
  double W2st = 2.0 * SubSN * SubTN;
  double Wnst = N * (SubSN + SubTN);
  double WNorm = W2st * (Wnst - W2st) / Wnst;

  // Main group criterion W
  double WMain = ((double)LinksST / SubTN - (double)LinksSTInv / (N - SubTN)) / SubSN;
  return WNorm * WMain;
}
double GroupW(TGroupST& G, PUNGraph& Graph, TIntV& SubSNIdV, TIntV& SubTNIdV) {
  // Recompute essential fields
  G.N = Graph->GetNodes();
  G.SubSNIdV = SubSNIdV;
  G.SubSN = SubSNIdV.Len();
  G.SubTNIdV = SubTNIdV;
  G.SubTN = SubTNIdV.Len();
  LinksCnt(G.LinksST, G.LinksSTInv, Graph, SubSNIdV, SubTNIdV, false);

  // Group criterion W
  G.W = GroupW(G.N, G.SubSN, G.SubTN, G.LinksST, G.LinksSTInv);
  return G.W;
}


/**
 * Fast estimate of group criterion W after swapping
 *
 * @param[out] G Output estimate of node group extraction (into S, T)
 * @param SubSNIdV List of node IDs in group S
 * @param SubTNIdV List of node IDs in group T
 * @param AddSNIdV Proposed list of node IDs to add into S
 * @param DelSNIdV Proposed list of node IDs to delete from S
 * @param AddTNIdV Proposed list of node IDs to add into T
 * @param DelTNIdV Proposed list of node IDs to delete from T
 * @param LinksST Number of edges L(S,T) between groups S and T
 * @param LinksSTInv Number of edges L(S,T) between groups S and inverse T
 * @return Estimated value of group critetion W
 */
double GroupWFast(TGroupST& G, PUNGraph& Graph, TIntV& SubSNIdV, TIntV& SubTNIdV, TIntV& AddSNIdV, TIntV& DelSNIdV, TIntV& AddTNIdV, TIntV& DelTNIdV, int LinksST, int LinksSTInv) {
  // Recompute essential fields
  G.N = Graph->GetNodes();
  G.SubSN = SubSNIdV.Len();
  G.SubTN = SubTNIdV.Len();
  G.LinksST = LinksST;
  G.LinksSTInv = LinksSTInv;
  int LinksSTDiff, LinksSTInvDiff;

  // Proposed delete from S
  G.SubSN -= DelSNIdV.Len();
  LinksCnt(LinksSTDiff, LinksSTInvDiff, Graph, DelSNIdV, SubTNIdV, false);
  G.LinksST -= LinksSTDiff;
  G.LinksSTInv -= LinksSTInvDiff;

  // Proposed add into S
  G.SubSN += AddSNIdV.Len();
  LinksCnt(LinksSTDiff, LinksSTInvDiff, Graph, AddSNIdV, SubTNIdV, false);
  G.LinksST += LinksSTDiff;
  G.LinksSTInv += LinksSTInvDiff;

  // Proposed delete from T
  G.SubTN -= DelTNIdV.Len();
  LinksCnt(LinksSTDiff, LinksSTInvDiff, Graph, DelTNIdV, SubSNIdV, false);
  G.LinksST -= LinksSTDiff;
  G.LinksSTInv += LinksSTDiff;

  // Proposed add into T
  G.SubTN += AddTNIdV.Len();
  LinksCnt(LinksSTDiff, LinksSTInvDiff, Graph, AddTNIdV, SubSNIdV, false);
  G.LinksST += LinksSTDiff;
  G.LinksSTInv -= LinksSTDiff;

  // Group criterion W
  G.W = GroupW(G.N, G.SubSN, G.SubTN, G.LinksST, G.LinksSTInv);
  return G.W;
}


/**
 * Node group extraction algorithm (into S, T)
 *
 * @param[out] GBest Output results of node group extraction (into S, T)
 * @param Graph Input graph
 * @param OptRestarts Number of restarts of the optimization algorithm
 * @param OptMxSteps Maximal number of steps in each optimization run
 * @param OptStopSteps Stop optimization if no W improvement in steps
 * @return Best value of group critetion W
 */
double GroupExtractSingle(TGroupST& GBest, PUNGraph& Graph, int OptMxSteps/*=DEF_OptMxSteps*/, int OptStopSteps/*=DEF_OptStopSteps*/) {
  // Initialize random groups S and T
  TGroupST G;
  G.SubSNIdV.AddMerged(Graph->GetRndNId());
  G.SubTNIdV.AddMerged(Graph->GetRndNId());
  // uncomment to start with random nodes in S and T
  //for (int i = TInt::GetRnd(Graph->GetNodes()); i > 0; --i) {
  //  NodesS.AddMerged(Graph->GetRndNId());
  //}
  //for (int i = TInt::GetRnd(Graph->GetNodes()); i > 0; --i) {
  //  NodesT.AddMerged(Graph->GetRndNId());
  //}
  // uncomment to start with all but one node
  //int SkipSNId = TInt::GetRnd(Graph->GetNodes());
  //int SkipTNId = TInt::GetRnd(Graph->GetNodes());
  //for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
  //  if (NI.GetId() != SkipSNId) {
  //    NodesS.AddMerged(NI.GetId());
  //  }
  //  if (NI.GetId() != SkipTNId) {
  //    NodesT.AddMerged(NI.GetId());
  //  }
  //}
  GroupW(G, Graph, G.SubSNIdV, G.SubTNIdV);
  double WBestPrev1, WBestPrev2;
  WBestPrev1 = WBestPrev2 = -INFINITY;
  GBest = G;

  // Hill climbing optimization
  /*for (int i = 0; i < OptMxSteps && GBest.W > WBestPrev2; ++i) {
    if (i % OptStopSteps == 0) {  // for early stopping
      WBestPrev2 = WBestPrev1;
      WBestPrev1 = GBest.W;
    }

    // Select and swap single node in either S or T
    bool DoSwapS = TBool::GetRnd();
    TIntV& SubX = (DoSwapS) ? G.SubSNIdV : G.SubTNIdV;
    int SwapNId = Graph->GetRndNId();
    if (SubX.IsIn(SwapNId)) {  // delete node
      if (SubX.Len() > 1) {
        SubX.DelAll(SwapNId);
      }
    } else {  // add node
      if (SubX.Len() + 1 < Graph->GetNodes()) {
        SubX.AddMerged(SwapNId);
      }
    }

    // Recompute group criterion W
    GroupW(G, Graph, G.SubSNIdV, G.SubTNIdV);

    // Evaluate for best
    if (G.W > GBest.W) {
      GBest = G;

    } else if (true) {  // false to use blind random walk
      if (SubX.IsIn(SwapNId)) {  // delete node again
        SubX.DelAll(SwapNId);
      } else {  // add node again
        SubX.AddMerged(SwapNId);
      }
    }
  }*/

  // Hill climbing optimization with fast estimation
  for (int i = 0; i < OptMxSteps && GBest.W > WBestPrev2; ++i) {
    if (i % OptStopSteps == 0) {  // for early stopping
      WBestPrev2 = WBestPrev1;
      WBestPrev1 = GBest.W;
    }

    // Select single node for proposed swapping in either S or T
    TIntV AddSNIdV, DelSNIdV, AddTNIdV, DelTNIdV;
    int SwapNId = Graph->GetRndNId();
    if (TBool::GetRnd()) {  // single swap in S
      if (GBest.SubSNIdV.IsIn(SwapNId)) {  // delete node from S
        if (GBest.SubSN > 1) {
          DelSNIdV.Add(SwapNId);
        }
      } else {  // add node into S
        if (GBest.SubSN + 1 < GBest.N) {
          AddSNIdV.Add(SwapNId);
        }
      }
    } else {  // single swap in T
      if (GBest.SubTNIdV.IsIn(SwapNId)) {  // delete node from T
        if (GBest.SubTN > 1) {
          DelTNIdV.Add(SwapNId);
        }
      } else {  // add node into T
        if(GBest.SubTN + 1 < GBest.N) {
          AddTNIdV.Add(SwapNId);
        }
      }
    }

    // Fast estimate of group criterion W after swapping
    GroupWFast(G, Graph, GBest.SubSNIdV, GBest.SubTNIdV, AddSNIdV, DelSNIdV, AddTNIdV, DelTNIdV, GBest.LinksST, GBest.LinksSTInv);

    // Evaluate for best
    if (G.W > GBest.W) {
      // Actually perform swap
      G.SubSNIdV.Diff(DelSNIdV);
      G.SubSNIdV.Union(AddSNIdV);
      G.SubTNIdV.Diff(DelTNIdV);
      G.SubTNIdV.Union(AddTNIdV);

      GBest = G;
    }
  }

  // Populate results of node group extraction
  GBest.RecomputeAll(Graph, GBest.SubSNIdV, GBest.SubTNIdV);
  return GBest.W;
}


/**
 * Restarter for finding best node group extraction (into S, T)
 *
 * @param[out] GBest Output results of best node group extraction (into S, T)
 * @param Graph Input graph
 * @param OptRestarts Number of restarts of the optimization algorithm
 * @param OptMxSteps Maximal number of steps in each optimization run
 * @param OptStopSteps Stop optimization if no W improvement in steps
 * @return Best value of group critetion W
 */
double GroupExtractRestarter(TGroupST& GBest, PUNGraph& Graph, int OptRestarts/*=DEF_OptRestarts*/, int OptMxSteps/*=DEF_OptMxSteps*/, int OptStopSteps/*=DEF_OptStopSteps*/) {
  TGroupST G;
  GBest.W = -INFINITY;

  for (int i = 0; i < OptRestarts; ++i) {
    // Find node group extraction (into S, T)
    GroupExtractSingle(G, Graph, OptMxSteps, OptStopSteps);

    // Evaluate for best
    if (G.W > GBest.W) {
      GBest = G;
    }
  }

  return GBest.W;
}


/**
 * Estimate average W on Erdos-Renyi random graphs
 *
 * @param[out] RAvg Output results of average node group extraction (into S, T)
 * @param N Number of nodes in random graphs
 * @param M Number of edges in random graphs
 * @param RndRestarts Number of restarts on Erdos-Renyi random graphs
 * @param OptMxSteps Maximal number of steps in each optimization run
 * @param OptStopSteps Stop optimization if no W improvement in steps
 * @return Average value of group critetion W
 */
double GroupExtractAvgRndGnm(TGroupST& RAvg, int N, int M, int RndRestarts/*=DEF_RndRestarts*/, int OptMxSteps/*=DEF_OptMxSteps*/, int OptStopSteps/*=DEF_OptStopSteps*/) {
  PUNGraph GraphER;
  TGroupST R;
  double WSum = 0.0;
  double WDiff;
  double WAvgDiff = INFINITY;

  for (int i = 0; i < RndRestarts; ++i) {
    // Generate a Erdos-Renyi random graph
    GraphER = TSnap::GenRndGnm<PUNGraph>(N, M, false);

    // Find node group extraction (into S, T)
    GroupExtractSingle(R, GraphER, OptMxSteps, OptStopSteps);

    // Evaluate for average
    WSum += R.W;
    if (i > RndRestarts / 5) {  // when enough data is available
      WDiff = R.W - WSum / (i + 1);
      if (WDiff < 0.0)
        WDiff = -WDiff;
      if (WDiff < WAvgDiff) {  // better estimate of average
        WAvgDiff = WDiff;
        RAvg = R;
      }
    }
  }

  return WSum / RndRestarts;
}


/**
 * Group extraction framework (into S, T)
 *
 * @param[out] GroupV Output list of results of node group extraction (into S, T)
 * @param[in,out] Graph Input graph
 * @param OptRestarts Number of restarts of the optimization algorithm
 * @param OptMxSteps Maximal number of steps in each optimization run
 * @param OptStopSteps Stop optimization if no W improvement in steps
 * @param RndRestarts Number of restarts on Erdos-Renyi random graphs
 * @param RndRecompW Force recomputation on random graphs if relative W difference smaller
 * @param RndStopW Stop group extraction if relative W difference smaller
 * @return Number of extracted groups
 */
int GroupExtractFramework(TGroupSTV& GroupV, PUNGraph& Graph, int OptRestarts/*=DEF_OptRestarts*/, int OptMxSteps/*=DEF_OptMxSteps*/, int OptStopSteps/*=DEF_OptStopSteps*/, int RndRestarts/*=DEF_RndRestarts*/, double RndRecompW/*=DEF_RndRecompW*/, double RndStopW/*=DEF_RndStopW*/) {
  Graph = TSnap::GetMxWcc(Graph);  // largest weakly-connected component
  TGroupST G, R;
  R.W = INFINITY;
  R.LinksST = 0;

  do {
    // Find node group extraction (into S, T)
    GroupExtractRestarter(G, Graph, OptRestarts, OptMxSteps, OptStopSteps);

    // Recompute on corresponding Erdos-Renyi random graphs
    if (G.W < RndRecompW * R.W || G.LinksST > R.LinksST) {
      R.W = GroupExtractAvgRndGnm(R, G.N, G.M, RndRestarts, OptMxSteps, OptStopSteps);
      if (R.W < 0.0)
        R.W = 0.0;
    }

    // Print status
    printf("\n");
    printf("%-3d %s\n", GroupV.Len(), G.GetStr().CStr());
    printf("  r %s\n", R.GetStr().CStr());
    GroupV.Add(G);

    // Delete links between S and T
    LinksCnt(Graph, G.SubSNIdV, G.SubTNIdV, true);
    Graph = TSnap::GetMxWcc(Graph);  // largest weakly-connected component
  } while(G.W > RndStopW * R.W);

  return GroupV.Len();
}
