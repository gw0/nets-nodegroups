/**
 * nodegroups - Group structures functions
 *
 * @author  gw0 [http://gw.tnode.com/] <gw.2014@tnode.com>
 * @version (see group_h_VERSION)
 */

#include "group.h"


/**
 * Convert to string
 *
 * @param this Input TGroupST object
 * @param Type Select output format type
 */
TStr TGroupST::GetStr(int Type/*=20*/) {
  switch (Type) {
    case 0:
      return TStr::Fmt("W=%-10.4f N=%-5d M=%-5d |S|=%-5d |T|=%-5d L(S,T)=%-5d Tau=%-6.4f Mod(S)=%-7.4f", W, N, M, SubSN, SubTN, LinksST, Tau, ModularityS);
    case 10:
      return TStr::Fmt("N=%-5d M=%-5d N_S=%-5d M_S=%-5d N_T=%-5d M_T=%-5d N_ST=%-5d M_ST=%-5d L_ST=%-5d L_STc=%-5d W=%-10.4f Tau=%-6.4f Mod_S=%-7.4f Mod_T=%-7.4f", N, M, SubSN, SubSM, SubTN, SubTM, SubSTN, SubSTM, LinksST, LinksSTc, W, Tau, ModularityS, ModularityT);
    case 11:
      return TStr::Fmt("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f", N, M, SubSN, SubSM, SubTN, SubTM, SubSTN, SubSTM, LinksST, LinksSTc, W, Tau, ModularityS, ModularityT);
    case 20:
      return TStr::Fmt("N=%-5d M=%-5d N_S=%-5d M_S=%-5d N_T=%-5d M_T=%-5d N_ST=%-5d M_ST=%-5d L_ST=%-5d L_STc=%-5d W=%-10.4f Tau=%-6.4f Mod_S=%-7.4f Mod_T=%-7.4f Type=%s", N, M, SubSN, SubSM, SubTN, SubTM, SubSTN, SubSTM, LinksST, LinksSTc, W, Tau, ModularityS, ModularityT, GroupName(SubSNIdV, SubTNIdV).CStr());
    case 21:
      return TStr::Fmt("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%s", N, M, SubSN, SubSM, SubTN, SubTM, SubSTN, SubSTM, LinksST, LinksSTc, W, Tau, ModularityS, ModularityT, GroupName(SubSNIdV, SubTNIdV).CStr());
  }
  return "";
}


/**
 * Recompute all results of node group extraction (except W)
 *
 * @param[out] this Output TGroupST object
 * @param Graph Input graph
 * @param NewSubSNIdV List of node IDs in group S
 * @param NewSubTNIdV List of node IDs in group T
 */
void TGroupST::RecomputeAll(const PUNGraph& Graph, const TIntV& NewSubSNIdV, const TIntV& NewSubTNIdV) {
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
  LinksCnt(LinksST, LinksSTc, Graph, SubSNIdV, SubTNIdV, false);

  Tau = GroupTau(SubSNIdV, SubTNIdV);
  W = GroupW(N, SubSN, SubTN, LinksST, LinksSTc);
  ModularityS = TSnap::GetModularity(Graph, SubSNIdV, N);
  ModularityT = TSnap::GetModularity(Graph, SubTNIdV, N);
}


/**
 * Count edges L(S,T) between groups S and T
 *
 * @param[out] LinksST Number of edges L(S,T) between groups S and T
 * @param[out] LinksSTc Number of edges L(S,T) between groups S and complement T
 * @param Graph Input graph
 * @param SubSNIdV List of node IDs in group S
 * @param SubTNIdV List of node IDs in group T
 * @param DoDelEdges Delete found edges between S and T
 * @return Number of edges L(S,T) between groups S and T
 */
double LinksCnt(int& LinksST, int& LinksSTc, const PUNGraph& Graph, const TIntV& SubSNIdV, const TIntV& SubTNIdV, bool DoDelEdges/*=false*/) {
  int SNId, TNId;
  LinksST = 0;
  LinksSTc = 0;

  // iterate through S
  for (int s = 0; s < SubSNIdV.Len(); ++s) {
    SNId = SubSNIdV[s];
    TUNGraph::TNodeI NI = Graph->GetNI(SNId);

    // iterate through its out-edges
    for (int e = 0; e < NI.GetOutDeg(); ++e) {
      TNId = NI.GetOutNId(e);

      if (SubTNIdV.IsInBin(TNId)) {  // endpoint is inside T
        //printf("edge (%d %d)\n", SNId, TNId);
        LinksST += 1;
        if (DoDelEdges) {  // delete edge between S and T
          Graph->DelEdge(SNId, TNId);
        }

      } else {  // endpoint is not inside T
        LinksSTc += 1;
      }
    }
  }

  return LinksST;
}
double LinksCnt(const PUNGraph& Graph, const TIntV& SubSNIdV, const TIntV& SubTNIdV, bool DoDelEdges/*=false*/) {
  int LinksST, LinksSTc;
  return LinksCnt(LinksST, LinksSTc, Graph, SubSNIdV, SubTNIdV, DoDelEdges);
}


/**
 * Group type parameter Tau(S,T)
 *
 * @param SubSNIdV List of node IDs in group S
 * @param SubTNIdV List of node IDs in group T
 * @return Value of group type parameter Tau
 */
double GroupTau(const TIntV& SubSNIdV, const TIntV& SubTNIdV) {
  return (double)SubSNIdV.IntrsLen(SubTNIdV) / SubSNIdV.UnionLen(SubTNIdV);
}


/**
 * Group type name
 *
 * Values:
 *   - "COM": community (S = T)
 *   - "MOD": module (S inters T = 0)
 *   - "HSD": hub&spokes module (module and |T| = 1)
 *   - "MIX": mixture (else)
 *   - "CPX": core/periphery mixture (S subset T or T subset S)
 *
 * @param SubSNIdV List of node IDs in group S
 * @param SubTNIdV List of node IDs in group T
 * @return Name of group type
 */
TStr GroupName(const TIntV& SubSNIdV, const TIntV& SubTNIdV) {
  int IntrsLen = SubSNIdV.IntrsLen(SubTNIdV);

  if(SubSNIdV == SubTNIdV) {  // community
    return "COM";

  } else if(IntrsLen == 0) {  // module
    if(SubTNIdV.Len() == 1)  // hub&spokes module
      return "HSD";
    else
      return "MOD";

  } else {  // mixture
    if(IntrsLen == SubSNIdV.Len() || IntrsLen == SubTNIdV.Len())  // core/periphery mixture
      return "CPX";
    else
      return "MIX";
  }
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
 * @param LinksSTc Number of edges L(S,T) between groups S and complement T
 * @param Graph Input graph
 * @param[out] G Output results of node group extraction (into S, T)
 * @return Value of group critetion W
 */
double GroupW(int N, int SubSN, int SubTN, int LinksST, int LinksSTc) {
  // Normalization factor (geometric mean of |S| and |T|)
  double WMean = 2.0 * SubSN * SubTN / (SubSN + SubTN);

  // Main group criterion W
  double WMain = ((double)LinksST / SubTN - (double)LinksSTc / (N - SubTN)) / SubSN;
  return WMean * (N - WMean) * WMain;
}
double GroupW(TGroupST& G, const PUNGraph& Graph, const TIntV& SubSNIdV, const TIntV& SubTNIdV) {
  // Recompute essential fields
  G.N = Graph->GetNodes();
  G.SubSNIdV = SubSNIdV;
  G.SubSN = SubSNIdV.Len();
  G.SubTNIdV = SubTNIdV;
  G.SubTN = SubTNIdV.Len();
  LinksCnt(G.LinksST, G.LinksSTc, Graph, SubSNIdV, SubTNIdV, false);

  // Group criterion W
  G.W = GroupW(G.N, G.SubSN, G.SubTN, G.LinksST, G.LinksSTc);
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
 * @param LinksSTc Number of edges L(S,T) between groups S and complement T
 * @return Estimated value of group critetion W
 */
double GroupWFast(TGroupST& G, const PUNGraph& Graph, const TIntV& SubSNIdV, const TIntV& SubTNIdV, const TIntV& AddSNIdV, const TIntV& DelSNIdV, const TIntV& AddTNIdV, const TIntV& DelTNIdV, int LinksST, int LinksSTc) {
  // Recompute essential fields
  G.N = Graph->GetNodes();
  G.SubSN = SubSNIdV.Len();
  G.SubTN = SubTNIdV.Len();
  G.LinksST = LinksST;
  G.LinksSTc = LinksSTc;
  int LinksSTDiff, LinksSTcDiff;

  // Proposed add into S
  if (AddSNIdV.Len() > 0) {
    G.SubSN += AddSNIdV.Len();
    LinksCnt(LinksSTDiff, LinksSTcDiff, Graph, AddSNIdV, SubTNIdV, false);
    G.LinksST += LinksSTDiff;
    G.LinksSTc += LinksSTcDiff;
  }

  // Proposed delete from S
  if (DelSNIdV.Len() > 0) {
    G.SubSN -= DelSNIdV.Len();
    LinksCnt(LinksSTDiff, LinksSTcDiff, Graph, DelSNIdV, SubTNIdV, false);
    G.LinksST -= LinksSTDiff;
    G.LinksSTc -= LinksSTcDiff;
  }

  // Proposed add into T
  if (AddTNIdV.Len() > 0) {
    G.SubTN += AddTNIdV.Len();
    LinksCnt(LinksSTDiff, LinksSTcDiff, Graph, AddTNIdV, SubSNIdV, false);
    G.LinksST += LinksSTDiff;
    G.LinksSTc -= LinksSTDiff;
  }

  // Proposed delete from T
  if (DelTNIdV.Len() > 0) {
    G.SubTN -= DelTNIdV.Len();
    LinksCnt(LinksSTDiff, LinksSTcDiff, Graph, DelTNIdV, SubSNIdV, false);
    G.LinksST -= LinksSTDiff;
    G.LinksSTc += LinksSTDiff;
  }

  // Group criterion W
  G.W = GroupW(G.N, G.SubSN, G.SubTN, G.LinksST, G.LinksSTc);
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
 * @param OptInitSample Initial random-sample size of S ant T (0 for random)
 * @return Best value of group critetion W
 */
double GroupExtractSingle(TGroupST& GBest, const PUNGraph& Graph, int OptMxSteps/*=DEF_OptMxSteps*/, int OptStopSteps/*=DEF_OptStopSteps*/, int OptInitSample/*=DEF_OptInitSample*/) {
  // Initialize random groups S and T
  TGroupST G = {};
  G.N = Graph->GetNodes();
  if (G.N < 2 || OptInitSample > G.N - 1 || -OptInitSample > G.N - 1) {  // invalid case
    OptInitSample = 1;
  } else if (OptInitSample == 0) {  // random initial size
    OptInitSample = TInt::GetRnd(G.N - 1) + 1;
  }
  if (OptInitSample > 0) {  // initial size as specified
    for (int i = 0; i < OptInitSample; ++i) {
      G.SubSNIdV.AddMerged(Graph->GetRndNId());
      G.SubTNIdV.AddMerged(Graph->GetRndNId());
    }
  } else {  // initial size as all but specified
    TIntV SubSInvNIdV, SubTcNIdV;
    for (int i = 0; i < -OptInitSample; ++i) {
      SubSInvNIdV.AddMerged(Graph->GetRndNId());
      SubTcNIdV.AddMerged(Graph->GetRndNId());
    }
    for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
      if (!SubSInvNIdV.IsInBin(NI.GetId())) {
        G.SubSNIdV.AddMerged(NI.GetId());
      }
      if (!SubTcNIdV.IsInBin(NI.GetId())) {
        G.SubTNIdV.AddMerged(NI.GetId());
      }
    }
  }

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
    if (SubX.IsInBin(SwapNId)) {  // delete node
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
      if (SubX.IsInBin(SwapNId)) {  // delete node again
        SubX.DelAll(SwapNId);
      } else {  // add node again
        SubX.AddMerged(SwapNId);
      }
    }
  }*/

  // Hill climbing optimization with fast estimation
  TIntV AddSNIdV, DelSNIdV, AddTNIdV, DelTNIdV;
  for (int i = 0; i < OptMxSteps && GBest.W > WBestPrev2; ++i) {
    if (i % OptStopSteps == 0) {  // for early stopping
      WBestPrev2 = WBestPrev1;
      WBestPrev1 = GBest.W;
    }

    // Select single node for proposed swapping in either S or T
    AddSNIdV.Clr();
    DelSNIdV.Clr();
    AddTNIdV.Clr();
    DelTNIdV.Clr();
    int SwapNId = Graph->GetRndNId();
    if (TBool::GetRnd()) {  // single swap in S
      if (!GBest.SubSNIdV.IsInBin(SwapNId)) {  // add node into S
        if (GBest.SubSN + 1 < GBest.N) {
          AddSNIdV.Add(SwapNId);
        }
      } else {  // delete node from S
        if (GBest.SubSN > 1) {
          DelSNIdV.Add(SwapNId);
        }
      }
    } else {  // single swap in T
      if (!GBest.SubTNIdV.IsInBin(SwapNId)) {  // add node into T
        if(GBest.SubTN + 1 < GBest.N) {
          AddTNIdV.Add(SwapNId);
        }
      } else {  // delete node from T
        if (GBest.SubTN > 1) {
          DelTNIdV.Add(SwapNId);
        }
      }
    }

    // Fast estimate of group criterion W after swapping
    GroupWFast(G, Graph, GBest.SubSNIdV, GBest.SubTNIdV, AddSNIdV, DelSNIdV, AddTNIdV, DelTNIdV, GBest.LinksST, GBest.LinksSTc);

    // Evaluate for best
    if (G.W > GBest.W) {
      // Actually perform swap
      if (AddSNIdV.Len() > 0)
        G.SubSNIdV.Union(AddSNIdV);
      if (DelSNIdV.Len() > 0)
        G.SubSNIdV.Diff(DelSNIdV);
      if (AddTNIdV.Len() > 0)
        G.SubTNIdV.Union(AddTNIdV);
      if (DelTNIdV.Len() > 0)
        G.SubTNIdV.Diff(DelTNIdV);

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
 * @param OptInitSample Initial random-sample size of S ant T (0 for random)
 * @return Best value of group critetion W
 */
double GroupExtractRestarter(TGroupST& GBest, const PUNGraph& Graph, int OptRestarts/*=DEF_OptRestarts*/, int OptMxSteps/*=DEF_OptMxSteps*/, int OptStopSteps/*=DEF_OptStopSteps*/, int OptInitSample/*=DEF_OptInitSample*/) {
  TGroupST G = {};
  GBest.W = -INFINITY;

  for (int i = 0; i < OptRestarts; ++i) {
    // Find node group extraction (into S, T)
    GroupExtractSingle(G, Graph, OptMxSteps, OptStopSteps, OptInitSample);

    // Evaluate for best
    if (G.W > GBest.W) {
      GBest = G;
    }
  }

  return GBest.W;
}


/**
 * Estimate group extraction results on random Erdos-Renyi graphs
 *
 * @param[out] GroupERV List of results of node group extraction (into S, T)
 * @param N Number of nodes in random graphs
 * @param M Number of edges in random graphs
 * @param RndGraphs Number of different Erdos-Renyi random graphs
 * @param RndRestarts Number of restarts on each Erdos-Renyi random graph
 * @param OptMxSteps Maximal number of steps in each optimization run
 * @param OptStopSteps Stop optimization if no W improvement in steps
 * @param OptInitSample Initial random-sample size of S ant T (0 for random)
 * @return Number of estimated group extraction results
 */
int GroupExtractRndGnms(TGroupSTV& GroupERV, int N, int M, int RndGraphs, int RndRestarts, int OptMxSteps, int OptStopSteps, int OptInitSample) {
  PUNGraph GraphER;

  for (int i = 0; i < RndGraphs; ++i) {
    // Generate a Erdos-Renyi random graph
    GraphER = TSnap::GenRndGnm<PUNGraph>(N, M, false);
    GraphER = TSnap::GetMxWcc(GraphER);  // largest weakly-connected component

    // Find node group extraction (into S, T)
    TGroupST R = {};
    GroupExtractRestarter(R, GraphER, RndRestarts, OptMxSteps, OptStopSteps, OptInitSample);

    // Add to list
    GroupERV.AddSorted(R);
  }

  return GroupERV.Len();
}


/**
 * Group extraction framework (into S, T)
 *
 * @param[out] GroupV Output list of results of node group extraction (into S, T)
 * @param[in,out] Graph Input graph
 * @param OptRestarts Number of restarts of the optimization algorithm
 * @param OptMxSteps Maximal number of steps in each optimization run
 * @param OptStopSteps Stop optimization if no W improvement in steps
 * @param OptInitSample Initial random-sample size of S ant T (0 for random)
 * @param FinishCnt Finish after extracting so many groups
 * @param FinishRndW Finish if W smaller than top percentile on random graphs
 * @param RndGraphs Number of different Erdos-Renyi random graphs
 * @param RndRestarts Number of restarts on each Erdos-Renyi random graph
 * @param RndRecompW Force W recomputation on random graphs when relative difference smaller
 * @return Number of extracted groups
 */
int GroupExtractFramework(TGroupSTV& GroupV, PUNGraph& Graph, int OptRestarts/*=DEF_OptRestarts*/, int OptMxSteps/*=DEF_OptMxSteps*/, int OptStopSteps/*=DEF_OptStopSteps*/, int OptInitSample/*=DEF_OptInitSample*/, int FinishCnt/*=DEF_FinishCnt*/, double FinishRndW/*=DEF_FinishRndW*/, int RndGraphs/*=DEF_RndGraphs*/, int RndRestarts/*=DEF_RndRestarts*/, double RndRecompW/*=DEF_RndRecompW*/) {
  Graph = TSnap::GetMxWcc(Graph);  // largest weakly-connected component
  TGroupST G = {}, R = {};
  R.W = INFINITY;

  do {
    // Find node group extraction (into S, T)
    GroupExtractRestarter(G, Graph, OptRestarts, OptMxSteps, OptStopSteps, OptInitSample);

    // Recompute on corresponding Erdos-Renyi random graphs
    if (G.W < RndRecompW * R.W || G.LinksST > R.LinksST) {
      TGroupSTV GroupERV;
      GroupExtractRndGnms(GroupERV, G.N, G.M, RndGraphs, RndRestarts, OptMxSteps, OptStopSteps, OptInitSample);
      R = GroupERV[(int)((GroupERV.Len() - 1) * (100.0 - FinishRndW) / 100.0)];
      if (R.W < 0.0)
        R.W = 0.0;
    }

    // Print status
    printf("\n");
    printf("%-3d %s\n", GroupV.Len(), G.GetStr().CStr());
    printf("  r %s\n", R.GetStr().CStr());

    if(G.W < R.W)
      break;
    GroupV.Add(G);

    // Delete links between S and T
    LinksCnt(Graph, G.SubSNIdV, G.SubTNIdV, true);
    Graph = TSnap::GetMxWcc(Graph);  // largest weakly-connected component
  } while(FinishCnt == 0 || GroupV.Len() < FinishCnt);

  return GroupV.Len();
}
