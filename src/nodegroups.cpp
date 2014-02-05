/**
 * nodegroups - Main code for node group structures
 *
 * Implementation of the *node group extraction framework* enables the exploration of node group structures of different networks, such as communities, modules, core/periphery, hubs&spokes, or similar structures. Description of the algorithm can be found in:
 *
 * - L. Šubelj, N. Blagus, and M. Bajec, "Group extraction for real-world networks: The case of communities, modules, and hubs and spokes," in Proc. of NetSci '13, 2013, p. 152.
 * - L. Šubelj, S. Žitnik, N. Blagus, and M. Bajec, "Node mixing and group structure of complex software networks," Adv. Complex Syst., 2014. (in review)
 *
 * @author  gw0 [http://gw.tnode.com/] <gw.2014@tnode.com>
 * @version 0.1
 */

#include "Snap.h"
#include "tgroups.h"


/**
 * Estimate average W in Erdos-Renyi random graph
 *
 * @param Iters Number of iterations of different random graphs and reruns
 * @param Steps Number of iterations of different subgraphs S and T
 * @param N Number of nodes in random graph
 * @param M Number of edges in random graph
 * @return Average value of group critetion W
 */
double WRndErdosRenyi(int Iters, int Steps, int N, int M) {
  PUNGraph GraphER;
  TIntV NodesS, NodesT;
  double W;
  double WRnd = 0.0;

  // iterations of different random graphs
  for (int i = 0; i < Iters; ++i) {
    // Generate Erdos-Renyi random graph
    GraphER = TSnap::GenRndGnm<PUNGraph>(N, M, false);

    // Find best group criterion W
    NodesS.Clr();
    NodesT.Clr();
    GroupExtract(GraphER, Steps, W, NodesS, NodesT);

    WRnd += W;
  }

  return WRnd / Iters;
}


/**
 * Console application entry point
 */
int main(int argc, char* argv[]) {
  // Header
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("nodegroups. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
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
  Graph = TSnap::GetMxWcc(Graph);  // largest weakly-connected component
  //Graph = TSnap::GetMxBiCon(Graph);  // largest bi-connected component
  TGroupV GroupV;
  TGroup g, r;
  double AlphaW = 0.01;
  r.W = INFINITY;
  r.LinksST = 0;
  do {
    // Extract group criterion W and subgraphs S and T
    g.NodesS.Clr();
    g.NodesT.Clr();
    GroupExtractRerunner(Graph, 1000, 10000, g.W, g.NodesS, g.NodesT);
    g.N = Graph->GetNodes();
    g.M = Graph->GetEdges();
    g.NodesSLen = g.NodesS.Len();
    g.NodesTLen = g.NodesT.Len();
    g.LinksST = GroupLinks(Graph, g.NodesS, g.NodesT);
    g.Tau = GroupTau(Graph, g.NodesS, g.NodesT);
    g.ModularityS = TSnap::GetModularity(Graph, g.NodesS, g.N);

    // Recompute corresponding Erdos-Renyi random graph
    if (g.W < r.W || r.LinksST < g.LinksST) {
      r = g;  // remember when it was built
      r.W = WRndErdosRenyi(100, 10000, r.N, r.M);
      if (r.W < 0.0)
        r.W = 0.0;
    }

    // Print status
    printf("%d  W=%-12.6f N=%-5d M=%-5d |S|=%-5d |T|=%-5d L(S,T)=%-5.0f Tau=%-8.6f Modu=%-9.6f r.W=%-12.6f\n", GroupV.Len(), g.W, g.N, g.M, g.NodesSLen, g.NodesTLen, g.LinksST, g.Tau, g.ModularityS, r.W);
    GroupV.Add(g);

    // Delete edges between S and T
    GroupLinks(Graph, g.NodesS, g.NodesT, true);
    Graph = TSnap::GetMxWcc(Graph);  // largest weakly-connected component
  } while(g.W * (1.0 - AlphaW) > r.W);

  // Output status and groups
  FILE *F = fopen(OutFNm.CStr(), "wt");
  fprintf(F, "# Input: %s\n", InFNm.CStr());
  fprintf(F, "# Nodes: %d    Edges: %d\n", GroupV[0].N, GroupV[0].M);
  fprintf(F, "# Groups: %d\n", GroupV.Len());
  for (int j = 0; j < GroupV.Len(); ++j) {
    TGroup& g = GroupV[j];
    fprintf(F, "#   %d  W=%-12.6f N=%-5d M=%-5d |S|=%-5d |T|=%-5d L(S,T)=%-5.0f Tau=%-8.6f Modu=%-9.6f\n", j, g.W, g.N, g.M, g.NodesSLen, g.NodesTLen, g.LinksST, g.Tau, g.ModularityS);
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
