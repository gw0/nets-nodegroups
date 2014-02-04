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

  // Count edges/links between subgraphs L(S,T)
  GroupLinks(Graph, NodesS, NodesT, false, LinksST, LinksSInvT);

  // Normalization factor (geometric mean of s and t)
  double W2st = 2.0 * NodesSLen * NodesTLen;
  double Wnst = N * (NodesSLen + NodesTLen);
  WNorm = W2st * (Wnst - W2st) / Wnst;

  // Group criterion
  return WNorm * (LinksST / NodesTLen - LinksSInvT / (N - NodesTLen)) / NodesSLen;
}
double GroupW(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT) {
  double WNorm, LinksST, LinksSInvT;
  return GroupW(Graph, NodesS, NodesT, WNorm, LinksST, LinksSInvT);
}


/**
 * Fast estimate of group criterion W after swapping a single node
 *
 * @param Graph Input graph
 * @param NodesS List of node IDs in subgraph S
 * @param NodesT List of node IDs in subgraph T
 * @param IfSwapS Estimate swap in S or T
 * @param SwapNId Node ID to add or delete from either S or T
 * @param WNorm Output W normalization factor
 * @param LinksST Output number of links between estimated S and T
 * @param LinksSInvT Output number of links between estimated S and inverse T
 * @return Estimated value of group critetion W
 */
double GroupWFast(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT, bool IfSwapS, int SwapNId, double& WNorm, double& LinksST, double& LinksSInvT) {
  int N = Graph->GetNodes();
  int NodesSLen = NodesS.Len();
  int NodesTLen = NodesT.Len();

  double LinksSTDiff, LinksSInvTDiff;
  TIntV TmpNIdV;
  TmpNIdV.Clr();
  TmpNIdV.Add(SwapNId);
  if (IfSwapS) {  // single swap in S
    if (NodesS.IsIn(SwapNId)) {  // delete node from S
      if (NodesSLen > 1) {
        NodesSLen -= 1;
        GroupLinks(Graph, TmpNIdV, NodesT, false, LinksSTDiff, LinksSInvTDiff);
        LinksST -= LinksSTDiff;
        LinksSInvT -= LinksSInvTDiff;
      }

    } else {  // add node into S
      if (NodesSLen + 1 < N) {
        NodesSLen += 1;
        GroupLinks(Graph, TmpNIdV, NodesT, false, LinksSTDiff, LinksSInvTDiff);
        LinksST += LinksSTDiff;
        LinksSInvT += LinksSInvTDiff;
      }
    }

  } else {  // single swap in T
    if (NodesT.IsIn(SwapNId)) {  // delete node from T
      if (NodesTLen > 1) {
        NodesTLen -= 1;
        GroupLinks(Graph, TmpNIdV, NodesS, false, LinksSTDiff, LinksSInvTDiff);
        LinksST -= LinksSTDiff;
        LinksSInvT += LinksSTDiff;
      }

    } else {  // add node into T
      if(NodesTLen + 1 < N) {
        NodesTLen += 1;
        GroupLinks(Graph, TmpNIdV, NodesS, false, LinksSTDiff, LinksSInvTDiff);
        LinksST += LinksSTDiff;
        LinksSInvT -= LinksSTDiff;
      }
    }
  }

  // Normalization factor (geometric mean of s and t)
  double W2st = 2.0 * NodesSLen * NodesTLen;
  double Wnst = N * (NodesSLen + NodesTLen);
  WNorm = W2st * (Wnst - W2st) / Wnst;

  // Corrected group criterion
  return WNorm * (LinksST / NodesTLen - LinksSInvT / (N - NodesTLen)) / NodesSLen;
}


/**
 * Node group extraction algorithm
 *
 * @param Graph Input graph
 * @param Steps Number of optimization steps
 * @param WBest Output best group criterion W
 * @param NodesSBest Output list of nodes IDs in best subgraph S
 * @param NodesTBest Output list of nodes IDs in best subgraph T
 * @return List of nodes IDs in best subgraph S
 */
TIntV GroupExtract(PUNGraph& Graph, int Steps, double& WBest, TIntV& NodesSBest, TIntV& NodesTBest) {
  // Initial random subgraph S and T (with one node)
  TIntV NodesS;
  NodesS.AddMerged(Graph->GetRndNId());
  TIntV NodesT;
  NodesT.AddMerged(Graph->GetRndNId());

  // Initial group criterion
  double W, WNorm, WNormBest, LinksST, LinksSTBest, LinksSInvT, LinksSInvTBest;
  WBest = GroupW(Graph, NodesS, NodesT, WNorm, LinksST, LinksSInvT);
  NodesSBest = NodesS;
  NodesTBest = NodesT;
  WNormBest = WNorm;
  LinksSTBest = LinksST;
  LinksSInvTBest = LinksSInvT;

  // Random walk optimization
  /*for (int i = 0; i < Steps; ++i) {
    // Select node to add or delete in either S or T
    bool IfSwapS = TBool::GetRnd();
    TIntV& NodesX = (IfSwapS) ? NodesS : NodesT;
    int SwapNId = Graph->GetRndNId();

    // Swap single node and recompute group criterion W
    if (NodesX.IsIn(SwapNId)) {  // delete node
      if (NodesX.Len() > 1) {
        NodesX.DelAll(SwapNId);
      }
    } else {  // add node
      if (NodesX.Len() + 1 < Graph->GetNodes()) {
        NodesX.AddMerged(SwapNId);
      }
    }
    W = GroupW(Graph, NodesS, NodesT, WNorm, LinksST, LinksSInvT);

    if (W > WBest) {
      //printf("%f %f %.0f %.0f\n", W, WNorm, LinksST, LinksSInvT);
      WBest = W;
      WNormBest = WNorm;
      LinksSTBest = LinksST;
      LinksSInvTBest = LinksSInvT;
      NodesSBest = NodesS;
      NodesTBest = NodesT;

    } else if (true) {  // enable to use greedy descent
      if (NodesX.IsIn(SwapNId)) {  // delete node again
        NodesX.DelAll(SwapNId);
      } else {  // add node again
        NodesX.AddMerged(SwapNId);
      }
    }
  }*/

  // Greedy descent optimization with fast estimation
  /*for (int i = 0; i < Steps; ++i) {
    // Select node to add or delete in either S or T
    bool IfSwapS = TBool::GetRnd();
    TIntV& NodesX = (IfSwapS) ? NodesS : NodesT;
    int SwapNId = Graph->GetRndNId();

    // Fast estimate of group criterion W after swapping a single node
    double WNormNew = WNorm;
    double LinksSTNew = LinksST;
    double LinksSInvTNew = LinksSInvT;
    W = GroupWFast(Graph, NodesS, NodesT, IfSwapS, SwapNId, WNormNew, LinksSTNew, LinksSInvTNew);

    if (W > WBest) {
      //printf("%f %f %.0f %.0f\n", W, WNorm, LinksST, LinksSInvT);
      WBest = W;
      WNormBest = WNormNew;
      LinksSTBest = LinksSTNew;
      LinksSInvTBest = LinksSInvTNew;

      // Actually swap best node
      if (NodesX.IsIn(SwapNId)) {  // delete node
        NodesX.DelAll(SwapNId);
      } else {  // add node
        NodesX.AddMerged(SwapNId);
      }
      NodesSBest = NodesS;
      NodesTBest = NodesT;
      WNorm = WNormBest;
      LinksST = LinksSTBest;
      LinksSInvT = LinksSInvTBest;
    }
  }*/

  // Steepest descent optimization with fast estimation
  for (int i = 0; i < (Steps / Graph->GetNodes() + 1); ++i) {
    // Select either S or T
    bool IfSwapS = TBool::GetRnd();
    TIntV& NodesX = (IfSwapS) ? NodesS : NodesT;
    int SwapNIdBest = -1;

    // iterate through all nodes in Graph
    for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
      int SwapNId = NI.GetId();

      // Fast estimate of group criterion W after swapping a single node
      double WNormNew = WNorm;
      double LinksSTNew = LinksST;
      double LinksSInvTNew = LinksSInvT;
      W = GroupWFast(Graph, NodesS, NodesT, IfSwapS, SwapNId, WNormNew, LinksSTNew, LinksSInvTNew);

      if (W > WBest) {
        SwapNIdBest = SwapNId;
        WBest = W;
        WNormBest = WNormNew;
        LinksSTBest = LinksSTNew;
        LinksSInvTBest = LinksSInvTNew;
      }
    }

    if (SwapNIdBest != -1) {
      //printf("%f %f %.0f %.0f\n", WBest, WNormBest, LinksSTBest, LinksSInvTBest);
      // Actually swap best node
      if (NodesX.IsIn(SwapNIdBest)) {  // delete node
        NodesX.DelAll(SwapNIdBest);
      } else {  // add node
        NodesX.AddMerged(SwapNIdBest);
      }
      NodesSBest = NodesS;
      NodesTBest = NodesT;
      WNorm = WNormBest;
      LinksST = LinksSTBest;
      LinksSInvT = LinksSInvTBest;
    }
  }

  return NodesSBest;
}


/**
 * Estimate maximum W in Erdos-Renyi random graph
 *
 * @param N Number of nodes in Erdos-Renyi random graph
 * @param M Number of edges in Erdos-Renyi random graph
 * @param NodesSLen Number of nodes in subgraph S
 * @param NodesTLen Number of nodes in subgraph T
 * @param Iters Number of iterations for finding maximum
 * @return Maximum value of group critetion W
 */
double WRndErdosRenyi(int N, int M, int NodesSLen, int NodesTLen, int Iters) {
  PUNGraph GraphER;
  TIntV NodesS, NodesT;
  double W;
  double WRnd = -INFINITY;

  // iterations for finding maximum
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
    W = GroupW(GraphER, NodesS, NodesT);
    if (W > WRnd) {
      WRnd = W;
    }
  }

  return WRnd;
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
  TGroup g, r;
  double AlphaW = 0.01;
  r.W = INFINITY;
  r.LinksST = 0;
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
    if (g.W < r.W || r.LinksST < g.LinksST) {
      r = g;  // remember when it was built
      r.W = WRndErdosRenyi(r.N, r.M, r.NodesSLen, r.NodesTLen, 10000);
      if (r.W < 0.0)
        r.W = 0.0;
    }

    // Print status
    printf("%d  W=%-12.6f N=%-5d M=%-5d |S|=%-5d |T|=%-5d L(S,T)=%-5.0f Tau=%-8.6f ; r.W=%-12.6f\n", GroupV.Len(), g.W, g.N, g.M, g.NodesSLen, g.NodesTLen, g.LinksST, g.Tau, r.W);
    GroupV.Add(g);

    // Delete edges between S and T
    GroupLinks(Graph, g.NodesS, g.NodesT, true);
  } while(g.W * (1.0 - AlphaW) > r.W);

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
