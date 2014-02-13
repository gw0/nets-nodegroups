/**
 * nodegroups - Main code for node group structures
 *
 * Implementation of the *node group extraction framework* (into S, T) enables the exploration of node group structures of different networks, such as communities, modules, core/periphery, hubs&spokes, or similar structures. Description of the algorithm can be found in:
 *
 * - L. Šubelj, N. Blagus, and M. Bajec, "Group extraction for real-world networks: The case of communities, modules, and hubs and spokes," in Proc. of NetSci '13, 2013, p. 152.
 * - L. Šubelj, S. Žitnik, N. Blagus, and M. Bajec, "Node mixing and group structure of complex software networks," Adv. Complex Syst., 2014. (in review)
 *
 * The adopted group extraction framework extracts groups from a simple undirected graph sequentially. An optimization method (currently random-restart hill climbing) is used to maximize the group criterion *W(S,T)* and extract group S with the corresponding linking pattern T. After extraction edges between S and T are removed and the whole process repeated on the largest weakly-connected component until the group criterion W is larger than expected on a Erdös-Rényi random graph.
 *
 * @author  gw0 [http://gw.tnode.com/] <gw.2014@tnode.com>
 * @version (see group_h_VERSION)
 */

#include "Snap.h"
#include "group.h"


/**
 * Console application entry point
 */
int main(int argc, char* argv[]) {
  // Header
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("nodegroups. Build: %.2f, %s, %s. Time: %s", group_h_VERSION, __TIME__, __DATE__, TExeTm::GetCurTm()), 1);
  TExeTm ExeTm;
  Try

  // Parameters
  const TStr PrefixFNm = Env.GetIfArgPrefixStr("-o:", "graph", "Input and output file name prefix (can be overriden)");
  const TStr InFNm = Env.GetIfArgPrefixStr("-i:", PrefixFNm + ".edgelist", "Input graph edges (undirected edge per line)");
  const TStr LabelFNm = Env.GetIfArgPrefixStr("-l:", PrefixFNm + ".labels", "Optional input node labels (node ID, node label)");
  const TStr OutFNm = Env.GetIfArgPrefixStr("-og:", PrefixFNm + ".groups", "Output group assignments (for S and T)");
  const TInt OptRestarts = Env.GetIfArgPrefixInt("-n:", DEF_OptRestarts, "Number of restarts of the optimization algorithm");
  const TInt OptMxSteps = Env.GetIfArgPrefixInt("-sm:", DEF_OptMxSteps, "Maximal number of steps in each optimization run");
  const TFlt OptStopSteps = Env.GetIfArgPrefixFlt("-sw:", DEF_OptStopSteps, "Stop optimization if no W improvement in steps");
  const TInt OptInitSample = Env.GetIfArgPrefixInt("-ss:", DEF_OptInitSample, "Initial random-sample size of S ant T (0=random)");
  const TInt RndRestarts = Env.GetIfArgPrefixInt("-rn:", DEF_RndRestarts, "Number of restarts on Erdos-Renyi random graphs");
  const TFlt RndRecompW = Env.GetIfArgPrefixFlt("-rf:", DEF_RndRecompW, "Force recomputation on random graphs if relative W difference smaller");
  const TFlt RndStopW = Env.GetIfArgPrefixFlt("-rw:", DEF_RndStopW, "Stop group extraction if relative W difference smaller");

  // Input
  PUNGraph Graph = TSnap::LoadEdgeList<PUNGraph>(InFNm, false);
  TIntStrH NIdLabelH;
  if (TFile::Exists(LabelFNm)) {  // optional labels
    TSsParser Ss(LabelFNm, ssfTabSep);
    while (Ss.Next()) {
      if (Ss.GetFlds() > 0) {
        NIdLabelH.AddDat(Ss.GetInt(0), Ss.GetFld(1));
      }
    }
  }

  // Run node group extraction framework
  TGroupSTV GroupV;
  GroupExtractFramework(GroupV, Graph, OptRestarts, OptMxSteps, OptStopSteps, OptInitSample, RndRestarts, RndRecompW, RndStopW);

  // Output
  FILE *F = fopen(OutFNm.CStr(), "wt");
  fprintf(F, "# Input: %s\n", InFNm.CStr());
  fprintf(F, "# Nodes: %d    Edges: %d\n", GroupV[0].N, GroupV[0].M);
  fprintf(F, "# Groups: %d\n", GroupV.Len());
  for (int j = 0; j < GroupV.Len(); ++j) {
    TGroupST& G = GroupV[j];
    fprintf(F, "#   %-3d %s\n", j, G.GetStr().CStr());
  }
  fprintf(F, "# NId\tGroupS\tGroupT\tNLabel\n");
  for (int j = 0; j < GroupV.Len(); ++j) {
    for (int i = 0; i < GroupV[j].SubSNIdV.Len(); ++i) {
      int NId = GroupV[j].SubSNIdV[i].Val;
      TStr NIdLabel = TStr("-");
      if (NIdLabelH.IsKey(NId))
        NIdLabel = NIdLabelH.GetDat(NId);

      if (GroupV[j].SubTNIdV.IsInBin(NId)) {  // in S and T
        fprintf(F, "%d\t%d\t%d\t%s\n", NId, j, j, NIdLabel.CStr());
      } else {  // in S, but not T
        fprintf(F, "%d\t%d\t-1\t%s\n", NId, j, NIdLabel.CStr());
      }
    }
    for (int i = 0; i < GroupV[j].SubTNIdV.Len(); ++i) {
      int NId = GroupV[j].SubTNIdV[i].Val;
      TStr NIdLabel = TStr("-");
      if (NIdLabelH.IsKey(NId))
        NIdLabel = NIdLabelH.GetDat(NId);

      if (!GroupV[j].SubSNIdV.IsInBin(NId)) {  // in T, but not S
        fprintf(F, "%d\t-1\t%d\t%s\n", NId, j, NIdLabel.CStr());
      }
    }
  }
  fclose(F);

  // Footer
  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
