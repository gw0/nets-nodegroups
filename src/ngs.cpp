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


/** Results of node group extraction */
struct TGroup {
  double W;
  int N;
  int M;
  int NodesSLen;
  int NodesTLen;
  double LinksST;
  double Tau;
  TIntV NodesS;
  TIntV NodesT;
};
typedef TVec<TGroup> TGroupV;


/**
 * Group type parameter Tau(S,T)
 *
 * @param Graph Input graph
 * @param NodesS List of node IDs in subgraph S
 * @param NodesT List of node IDs in subgraph T
 * @return Value of group type parameter Tau
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
 * @param IfDelEdge Delete visited edge between S and T
 * @param LinksST Output number of links between S and T
 * @param LinksSInvT Output number of links between S and inverse T
 * @return Number of links between S and T
 */
double GroupLinks(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT, bool IfDelEdge, double& LinksST, double& LinksSInvT) {
  LinksST = 0.0;
  LinksSInvT = 0.0;

  // iterate through S
  for (int s = 0; s < NodesS.Len(); ++s) {
    TUNGraph::TNodeI NI = Graph->GetNI(NodesS[s]);

    // iterate through its out-edges
    for (int e = 0; e < NI.GetOutDeg(); ++e) {
      // check if endpoint is inside T or not
      if (NodesT.IsIn(NI.GetOutNId(e))) {
        //printf("edge (%d %d)\n", NI.GetId(), NI.GetOutNId(e));
        LinksST += 1.0;
        if (IfDelEdge) {  // delete edge between S and T
          Graph->DelEdge(NodesS[s], NI.GetOutNId(e));
        }
      } else {
        LinksSInvT += 1.0;
      }
    }
  }
  return LinksST;
}
double GroupLinks(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT, bool IfDelEdge=false) {
  double LinksST, LinksSInvT;
  return GroupLinks(Graph, NodesS, NodesT, IfDelEdge, LinksST, LinksSInvT);
}


/**
 * Compute group criterion W(S,T)
 *
 * @param Graph Input graph
 * @param NodesS List of node IDs in subgraph S
 * @param NodesT List of node IDs in subgraph T
 * @param WNorm Output W normalization factor
 * @param LinksST Output number of links between S and T
 * @param LinksSInvT Output number of links between S and inverse T
 * @return Computed value of group critetion W
 */
double GroupW(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT, double& WNorm, double& LinksST, double& LinksSInvT) {
  int N = Graph->GetNodes();
  int NodesSLen = NodesS.Len();
  int NodesTLen = NodesT.Len();

  // Normalization factor (geometric mean of s and t)
  double W2st = 2.0 * NodesSLen * NodesTLen;
  double Wnst = N * (NodesSLen + NodesTLen);
  WNorm = W2st * (Wnst - W2st) / Wnst;

  // Count edges/links between subgraphs L(S,T)
  GroupLinks(Graph, NodesS, NodesT, false, LinksST, LinksSInvT);

  // Group criterion
  return WNorm * (LinksST / NodesTLen - LinksSInvT / (N - NodesTLen)) / NodesSLen;
}
double GroupW(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT) {
  double WNorm, LinksST, LinksSInvT;
  return GroupW(Graph, NodesS, NodesT, WNorm, LinksST, LinksSInvT);
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
 * @param WNorm Precomputed W normalization factor
 * @param LinksST Output number of links between estimated S and T
 * @param LinksSInvT Output number of links between estimated S and inverse T
 * @return Estimated value of group critetion W
 */
double GroupWFast(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT, bool IfSwapS, int DelNId, int AddNId, double& WNorm, double& LinksST, double& LinksSInvT) {
  int Len = Graph->GetNodes();
  int NodesSLen = NodesS.Len();
  int NodesTLen = NodesT.Len();

  // Correction for removing node DelNId
  TIntV TmpNIdV;
  double LinksSTTmp, LinksSInvTTmp;
  TmpNIdV.Clr();
  TmpNIdV.Add(DelNId);
  if (IfSwapS) {  // swap in S
    GroupLinks(Graph, TmpNIdV, NodesT, false, LinksSTTmp, LinksSInvTTmp);
    LinksST -= LinksSTTmp;
    LinksSInvT -= LinksSInvTTmp;

  } else {  // swap in T
    GroupLinks(Graph, TmpNIdV, NodesS, false, LinksSTTmp, LinksSInvTTmp);
    LinksST -= LinksSTTmp;
    LinksSInvT += LinksSTTmp;
  }

  // Correction for adding node AddNId
  TmpNIdV.Clr();
  TmpNIdV.Add(AddNId);
  if (IfSwapS) {  // swap in S
    GroupLinks(Graph, TmpNIdV, NodesT, false, LinksSTTmp, LinksSInvTTmp);
    LinksST += LinksSTTmp;
    LinksSInvT += LinksSInvTTmp;

  } else {  // swap in T
    GroupLinks(Graph, TmpNIdV, NodesS, false, LinksSTTmp, LinksSInvTTmp);
    LinksST += LinksSTTmp;
    LinksSInvT -= LinksSTTmp;      
  }

  // Corrected group criterion
  return WNorm * (LinksST / NodesTLen - LinksSInvT / (Len - NodesTLen)) / NodesSLen;
}


/**
 * Node group extraction algorithm
 *
 * @param Graph Input graph
 * @param Iters Number of optimization iterations
 * @param WBest Output best group criterion W
 * @param NodesSBest Output list of nodes IDs in best subgraph S
 * @param NodesTBest Output list of nodes IDs in best subgraph T
 * @return List of nodes IDs in best subgraph S
 */
TIntV GroupExtract(PUNGraph& Graph, int Iters, double& WBest, TIntV& NodesSBest, TIntV& NodesTBest) {
  // Initial random subgraph S and T (not empty or whole)
  int N = Graph->GetNodes();
  TIntV NodesS;
  for (int i = 1 + TInt::GetRnd(N - 2); i > 0; --i) {
    NodesS.AddMerged(Graph->GetRndNId());
  }
  TIntV NodesT;
  for (int i = 1 + TInt::GetRnd(N - 2); i > 0; --i) {
    NodesT.AddMerged(Graph->GetRndNId());
  }

  // Initial group criterion
  double W, WNorm, LinksST, LinksSInvT;
  WBest = GroupW(Graph, NodesS, NodesT, WNorm, LinksST, LinksSInvT);
  NodesSBest = NodesS;
  NodesTBest = NodesT;

  // Inefficient random walk optimization
  for (int i = 0; i < Iters; ++i) {
    // Select node to add or delete in either S or T
    bool IfSwapS = TBool::GetRnd();
    TIntV& NodesX = (IfSwapS) ? NodesS : NodesT;
    int SwapNId = Graph->GetRndNId();

    // Add or delete node and recompute group criterion
    if (NodesX.IsIn(SwapNId) && NodesX.Len() > 1) {  // delete node
      NodesX.DelAll(SwapNId);
    } else if(NodesX.Len() + 1 < N) {  // add node
      NodesX.AddMerged(SwapNId);
    }
    W = GroupW(Graph, NodesS, NodesT, WNorm, LinksST, LinksSInvT);

    if (W > WBest) {
      //printf("%f %f %f %f\n", W, WNorm, LinksST, LinksSInvT);
      WBest = W;
      NodesSBest = NodesS;
      NodesTBest = NodesT;

    } else if (false) {  // enable to use greedy descent
      if (NodesX.IsIn(SwapNId)) {  // delete node again
        NodesX.DelAll(SwapNId);
      } else {  // add node again
        NodesX.AddMerged(SwapNId);
      }
    }
  }

  // Greedy descent optimization with fast estimate of group criterion
  /*for (int i = 0; i < Iters; ++i) {
    // Select nodes to swap in either S or T
    bool IfSwapS = TBool::GetRnd();
    TIntV& NodesX = (IfSwapS) ? NodesS : NodesT;
    int DelI = TUInt::GetRnd(NodesX.Len());
    int AddNId;
    do {
      AddNId = Graph->GetRndNId();
    } while(NodesX.IsIn(AddNId));

    // Fast estimate of group criterion after swapping nodes
    double LinksSTNew = LinksST;
    double LinksSInvTNew = LinksSInvT;
    W = GroupWFast(Graph, NodesS, NodesT, IfSwapS, NodesX[DelI], AddNId, WNorm, LinksSTNew, LinksSInvTNew);

    if (W > WBest) {
      NodesX.Del(DelI);
      NodesX.AddMerged(AddNId);
      LinksST = LinksSTNew;
      LinksSInvT = LinksSInvTNew;

      //printf("%f %f %f %f\n", W, WNorm, LinksSTNew, LinksSInvTNew);
      WBest = W;
      NodesSBest = NodesS;
      NodesTBest = NodesT;
    }
  }*/

  return NodesSBest;
}


/**
 * Estimate average W in Erdos-Renyi random graph
 *
 * @param N Number of nodes in Erdos-Renyi random graph
 * @param M Number of edges in Erdos-Renyi random graph
 * @param NodesSLen Number of nodes in subgraph S
 * @param NodesTLen Number of nodes in subgraph T
 * @param Iters Number of iterations for average
 * @return Average value of group critetion W
 */
double WRndErdosRenyi(int N, int M, int NodesSLen, int NodesTLen, int Iters) {
  PUNGraph GraphER;
  TIntV NodesS, NodesT;
  double WRnd = 0.0;

  // iterations for average
  for (int i = 0; i < Iters; ++i) {
    // Generate Erdos-Renyi random graph
    GraphER = TSnap::GenRndGnm<PUNGraph>(N, M, false);

    // Select random subgraphs S and T
    NodesS.Clr();
    for (int i = 0; NodesS.Len() < NodesSLen; ++i) {
      NodesS.AddMerged(GraphER->GetRndNId());
    }
    NodesT.Clr();
    for (int i = 0; NodesT.Len() < NodesTLen; ++i) {
      NodesT.AddMerged(GraphER->GetRndNId());
    }

    // Compute group criterion W
    WRnd += GroupW(GraphER, NodesS, NodesT);
  }

  return WRnd / Iters;
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
  TGroupV GroupV;
  TGroup g;
  double AlphaW = 0.01;
  double WRnd = INFINITY;
  do {
    // Extract group criterion W and subgraphs S and T
    GroupExtract(Graph, 1000000, g.W, g.NodesS, g.NodesT);
    g.N = Graph->GetNodes();
    g.M = Graph->GetEdges();
    g.NodesSLen = g.NodesS.Len();
    g.NodesTLen = g.NodesT.Len();
    g.LinksST = GroupLinks(Graph, g.NodesS, g.NodesT);
    g.Tau = GroupTau(Graph, g.NodesS, g.NodesT);

    // Recompute corresponding Erdos-Renyi random graph
    if (g.W < WRnd && WRnd > 0.0) {
      WRnd = WRndErdosRenyi(g.N, g.M, g.NodesSLen, g.NodesTLen, 1000);
      if (WRnd < 0.0) {
        WRnd = 0.0;
      }
    }

    // Print status
    printf("%d  W=%-12.6f N=%-5d M=%-5d |S|=%-5d |T|=%-5d L(S,T)=%-5.0f Tau=%-8.6f ; WRnd=%-12.6f\n", GroupV.Len(), g.W, g.N, g.M, g.NodesSLen, g.NodesTLen, g.LinksST, g.Tau, WRnd);
    GroupV.Add(g);

    // Delete edges between S and T
    GroupLinks(Graph, g.NodesS, g.NodesT, true);
  } while(g.W * (1.0 - AlphaW) > WRnd);

  // Output status and groups
  FILE *F = fopen(OutFNm.CStr(), "wt");
  fprintf(F, "# Input: %s\n", InFNm.CStr());
  fprintf(F, "# Nodes: %d    Edges: %d\n", GroupV[0].N, GroupV[0].M);
  fprintf(F, "# Groups: %d\n", GroupV.Len());
  for (int j = 0; j < GroupV.Len(); ++j) {
    TGroup& g = GroupV[j];
    fprintf(F, "#   %d  W=%-12.6f N=%-5d M=%-5d |S|=%-5d |T|=%-5d L(S,T)=%-5.0f Tau=%-8.6f\n", j, g.W, g.N, g.M, g.NodesSLen, g.NodesTLen, g.LinksST, g.Tau);
  }
  fprintf(F, "# NId\tGroupId\tNLabel\n");
  for (int j = 0; j < GroupV.Len(); ++j) {
    for (int i = 0; i < GroupV[j].NodesS.Len(); ++i) {
      int NId = GroupV[j].NodesS[i].Val;
      if (NIDNameH.IsKey(NId)) {
        fprintf(F, "%d\t%d\t%s\n", NId, j, NIDNameH.GetDat(NId).CStr());
      } else {
        fprintf(F, "%d\t%d\t-\n", NId, j);
      }
    }
  }
  fclose(F);

  // Footer
  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
