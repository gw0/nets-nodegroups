/**
 * nodegroups - TGroups functions
 *
 * @author  gw0 [http://gw.tnode.com/] <gw.2014@tnode.com>
 * @version 0.1
 */

#include "tgroups.h"


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
 * @param[out] LinksST Output number of links between S and T
 * @param[out] LinksSInvT Output number of links between S and inverse T
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
double GroupLinks(PUNGraph& Graph, TIntV& NodesS, TIntV& NodesT, bool IfDelEdge/*=false*/) {
  double LinksST, LinksSInvT;
  return GroupLinks(Graph, NodesS, NodesT, IfDelEdge, LinksST, LinksSInvT);
}


/**
 * Compute group criterion W(S,T)
 *
 * @param Graph Input graph
 * @param NodesS List of node IDs in subgraph S
 * @param NodesT List of node IDs in subgraph T
 * @param[out] WNorm Output W normalization factor
 * @param[out] LinksST Output number of links between S and T
 * @param[out] LinksSInvT Output number of links between S and inverse T
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
 * @param[in,out] WNorm Output W normalization factor
 * @param[in,out] LinksST Output number of links between estimated S and T
 * @param[in,out] LinksSInvT Output number of links between estimated S and inverse T
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
 * @param[out] W Output best group criterion W
 * @param[in,out] NodesS Initial and output list of nodes IDs in best subgraph S
 * @param[in,out] NodesT Initial and output list of nodes IDs in best subgraph T
 * @return Best group criterion W
 */
double GroupExtract(PUNGraph& Graph, int Steps, double& W, TIntV& NodesS, TIntV& NodesT) {
  // Initial random subgraph S and T should containt at least one node
  if (NodesS.Len() < 1)
    NodesS.AddMerged(Graph->GetRndNId());
  if (NodesT.Len() < 1)
    NodesT.AddMerged(Graph->GetRndNId());

  // Initial group criterion W
  double WNorm, LinksST, LinksSInvT;
  double WBest, WNormBest, LinksSTBest, LinksSInvTBest;
  TIntV NodesSBest, NodesTBest;
  WBest = GroupW(Graph, NodesS, NodesT, WNorm, LinksST, LinksSInvT);
  WNormBest = WNorm;
  LinksSTBest = LinksST;
  LinksSInvTBest = LinksSInvT;
  NodesSBest = NodesS;
  NodesTBest = NodesT;

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
  for (int i = 0; i < Steps; ++i) {
    // Select node to add or delete in either S or T
    bool IfSwapS = TBool::GetRnd();
    TIntV& NodesX = (IfSwapS) ? NodesS : NodesT;
    int SwapNId = Graph->GetRndNId();

    // Fast estimate of group criterion W after swapping a single node
    double WNormNew = WNorm;
    double LinksSTNew = LinksST;
    double LinksSInvTNew = LinksSInvT;
    double WNew = GroupWFast(Graph, NodesS, NodesT, IfSwapS, SwapNId, WNormNew, LinksSTNew, LinksSInvTNew);

    if (WNew > WBest) {
      //printf("%f %f %.0f %.0f\n", WNew, WNormNew, LinksSTNew, LinksSInvTNew);
      WBest = WNew;
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
      W = WBest;
      WNorm = WNormBest;
      LinksST = LinksSTBest;
      LinksSInvT = LinksSInvTBest;
    }
  }

  // Steepest descent optimization with fast estimation
  /*for (int i = 0; i < (Steps / Graph->GetNodes() + 1); ++i) {
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
      double WNew = GroupWFast(Graph, NodesS, NodesT, IfSwapS, SwapNId, WNormNew, LinksSTNew, LinksSInvTNew);

      if (WNew > WBest) {
        SwapNIdBest = SwapNId;
        WBest = WNew;
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
      W = WBest;
      WNorm = WNormBest;
      LinksST = LinksSTBest;
      LinksSInvT = LinksSInvTBest;
    }
  }*/

  W = WBest;
  NodesS = NodesSBest;
  NodesT = NodesTBest;
  return W;
}


/**
 * Rerunner of node group extraction algorithm
 *
 * @param Graph Input graph
 * @param Iters Number of reruns
 * @param Steps Number of optimization steps
 * @param[out] WBest Output best group criterion W
 * @param[out] NodesSBest Output list of nodes IDs in best subgraph S
 * @param[out] NodesTBest Output list of nodes IDs in best subgraph T
 * @return Best group criterion W
 */
double GroupExtractRerunner(PUNGraph& Graph, int Iters, int Steps, double& WBest, TIntV& NodesSBest, TIntV& NodesTBest) {
  TIntV NodesS, NodesT;
  WBest = -INFINITY;

  // iterations for rerunning
  for (int i = 0; i < Iters; ++i) {
    double W;
    NodesS.Clr();
    NodesT.Clr();
    // uncomment to start with random nodes in S and T
    for (int i = TInt::GetRnd(Graph->GetNodes()); i > 0; --i) {
      NodesS.AddMerged(Graph->GetRndNId());
    }
    for (int i = TInt::GetRnd(Graph->GetNodes()); i > 0; --i) {
      NodesT.AddMerged(Graph->GetRndNId());
    }
    // uncomment to start with all but one node
    /*int SkipSNId = TInt::GetRnd(Graph->GetNodes());
    int SkipTNId = TInt::GetRnd(Graph->GetNodes());
    for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
      if (NI.GetId() != SkipSNId) {
        NodesS.AddMerged(NI.GetId());
      }
      if (NI.GetId() != SkipTNId) {
        NodesT.AddMerged(NI.GetId());
      }
    }*/
    GroupExtract(Graph, Steps, W, NodesS, NodesT);

    if (W > WBest) {
      WBest = W;
      NodesSBest = NodesS;
      NodesTBest = NodesT;
    }
  }

  return WBest;
}
