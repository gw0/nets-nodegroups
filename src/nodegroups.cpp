/**
 * nets-nodegroups - Research node group structure of networks
 *
 * To research and explore node group structure of various networks, we introduced the *group type parameter Tau* and implemented a *node group extraction framework*. The framework is capable of identifying groups, such as communities, modules, core/periphery, hubs&spokes, or similar structures. Detailed description of the algorithm can be found in:
 *
 * - L. Šubelj, N. Blagus, and M. Bajec, "Group extraction for real-world networks: The case of communities, modules, and hubs and spokes," in Proc. of NetSci '13, 2013, p. 152.
 * - L. Šubelj, S. Žitnik, N. Blagus, and M. Bajec, “Node mixing and group structure of complex software networks,” Adv. Complex Syst., vol. 17, 2014.
 *
 * The adopted network node group extraction framework extracts groups from a simple undirected graph sequentially. An optimization method (currently random-restart hill climbing) is used to maximize the group criterion *W(S,T)* and extract group *S* with the corresponding linking pattern *T*. After extraction edges between *S* and *T* are removed and the whole process repeated on the largest weakly-connected component until the group criterion *W* is larger than expected on a Erdös-Rényi random graph.
 *
 * @author  gw0 [http://gw.tnode.com/] <gw.2014@tnode.com>
 * @version (see group_h_VERSION)
 */

#include "Snap.h"
#include "group.h"


/**
 * Output ST-group assignments (*.groups)
 *
 * @param OutFNm Output file name
 * @param OutMode Output fopen mode (eg. "wt")
 * @param argc Number of program parameters
 * @param argv Program parameters
 * @param InFNm Input file name
 * @param GroupV List of ST-group extraction results
 */
void OutputGroups(const TStr& OutFNm, const TStr& OutMode, const int argc, /*const*/ char* argv[], const TStr& InFNm, const int GraphN, const int GraphM, const TIntStrH& NIdLabelH, /*const*/ TGroupSTV& GroupV) {
  FILE* F = fopen(OutFNm.CStr(), OutMode.CStr());
  fprintf(F, "# nets-nodegroups. Build: %.2f, %s, %s. Time: %s\n", group_h_VERSION, __TIME__, __DATE__, TExeTm::GetCurTm());
  fprintf(F, "#");
  for (int i = 0; i < argc; ++i) {
    fprintf(F, " %s", argv[i]);
  }
  fprintf(F, "\n");
  fprintf(F, "# Input: %s  Nodes: %d  Edges: %d\n", InFNm.CStr(), GraphN, GraphM);
  fprintf(F, "# Group results: %d\n", GroupV.Len());
  for (int j = 0; j < GroupV.Len(); ++j) {
    TGroupST& G = GroupV[j];
    fprintf(F, "#   %-3d %s\n", j, G.GetStr().CStr());
  }
  fprintf(F, "# NId\tGroupS\tGroupT\tNLabel\n");
  for (int j = 0; j < GroupV.Len(); ++j) {
    // in S and T
    for (int i = 0; i < GroupV[j].SubSNIdV.Len(); ++i) {
      int NId = GroupV[j].SubSNIdV[i].Val;
      TStr NIdLabel = TStr("-");
      if (NIdLabelH.IsKey(NId))
        NIdLabel = NIdLabelH.GetDat(NId);

      if (GroupV[j].SubTNIdV.IsInBin(NId)) {
        fprintf(F, "%d\t%d\t%d\t%s\n", NId, j, j, NIdLabel.CStr());
      }
    }
    // in S, but not T
    for (int i = 0; i < GroupV[j].SubSNIdV.Len(); ++i) {
      int NId = GroupV[j].SubSNIdV[i].Val;
      TStr NIdLabel = TStr("-");
      if (NIdLabelH.IsKey(NId))
        NIdLabel = NIdLabelH.GetDat(NId);

      if (!GroupV[j].SubTNIdV.IsInBin(NId)) {
        fprintf(F, "%d\t%d\t-1\t%s\n", NId, j, NIdLabel.CStr());
      }
    }
    // in T, but not S
    for (int i = 0; i < GroupV[j].SubTNIdV.Len(); ++i) {
      int NId = GroupV[j].SubTNIdV[i].Val;
      TStr NIdLabel = TStr("-");
      if (NIdLabelH.IsKey(NId))
        NIdLabel = NIdLabelH.GetDat(NId);

      if (!GroupV[j].SubSNIdV.IsInBin(NId)) {
        fprintf(F, "%d\t-1\t%d\t%s\n", NId, j, NIdLabel.CStr());
      }
    }
  }
  fclose(F);
}


/**
 * Output ST-group extraction results summary (*.groupssum)
 *
 * @param OutSumFNm Output file name
 * @param OutMode Output fopen mode (eg. "wt")
 * @param argc Number of program parameters
 * @param argv Program parameters
 * @param InFNm Input file name
 * @param GroupV List of ST-group extraction results
 */
void OutputGroupsSum(const TStr& OutSumFNm, const TStr& OutMode, const int argc, /*const*/ char* argv[], const TStr& InFNm, const int GraphN, const int GraphM, const TIntStrH& NIdLabelH, /*const*/ TGroupSTV& GroupV) {
  FILE* F = fopen(OutSumFNm.CStr(), OutMode.CStr());
  fprintf(F, "# nets-nodegroups. Build: %.2f, %s, %s. Time: %s\n", group_h_VERSION, __TIME__, __DATE__, TExeTm::GetCurTm());
  fprintf(F, "#");
  for (int i = 0; i < argc; ++i) {
    fprintf(F, " %s", argv[i]);
  }
  fprintf(F, "\n");
  fprintf(F, "# Input: %s  Nodes: %d  Edges: %d\n", InFNm.CStr(), GraphN, GraphM);
  fprintf(F, "# Group results: %d\n", GroupV.Len());
  fprintf(F, "N\tM\tN_S\tM_S\tN_T\tM_T\tN_ST\tM_ST\tL_ST\tL_STc\tW\tTau\tMod_S\tMod_T\tType\n");
  for (int j = 0; j < GroupV.Len(); ++j) {
    TGroupST& G = GroupV[j];
    fprintf(F, "%s\n", G.GetStr(21).CStr());
  }
  fclose(F);
}


/**
 * Console application entry point
 *
 * Return codes:
 *   - `-1`: error, file not found or crash
 *   - `0` : success, groups extracted
 *   - `1` : no groups extracted 
 */
int main(int argc, char* argv[]) {
  // Header
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("nets-nodegroups. Build: %.2f, %s, %s. Time: %s", group_h_VERSION, __TIME__, __DATE__, TExeTm::GetCurTm()), 1);
  TExeTm ExeTm;
  int ret = -1;
  Try

  // Parameters
  const TStr PrefixFNm = Env.GetIfArgPrefixStr("-o:", "graph", "Prefix for all file names (simplified usage)");
  const TStr InFNm = Env.GetIfArgPrefixStr("-i:", PrefixFNm + ".edgelist", "Input file with graph edges (undirected edge per line)");
  const TStr LabelFNm = Env.GetIfArgPrefixStr("-l:", PrefixFNm + ".labels", "Input file (optional) with node labels (node ID, node label)");
  const TStr OutFNm = Env.GetIfArgPrefixStr("-og:", PrefixFNm + ".groups", "Output file with ST-group assignments");
  const TStr OutSumFNm = Env.GetIfArgPrefixStr("-os:", PrefixFNm + ".groupssum", "Output file with only ST-group extraction summary");
  const TInt OptRestarts = Env.GetIfArgPrefixInt("-n:", DEF_OptRestarts, "Number of restarts of the optimization algorithm");
  const TInt OptMxSteps = Env.GetIfArgPrefixInt("-sm:", DEF_OptMxSteps, "Maximal number of steps in each optimization run");
  const TInt OptStopSteps = Env.GetIfArgPrefixInt("-sw:", DEF_OptStopSteps, "Stop optimization if no W improvement in steps");
  const TInt OptInitSample = Env.GetIfArgPrefixInt("-ss:", DEF_OptInitSample, "Initial random-sample size of S ant T (0=random)");
  const TInt FinishCnt = Env.GetIfArgPrefixInt("-fn:", DEF_FinishCnt, "Finish after extracting so many groups (turn off random graphs)");
  const TFlt FinishRndW = Env.GetIfArgPrefixFlt("-fw:", DEF_FinishRndW, "Finish if W smaller than top percentile on random graphs");
  const TInt RndGraphs = Env.GetIfArgPrefixInt("-rg:", DEF_RndGraphs, "Random graphs (Erdos-Renyi) to construct for estimating W");
  const TInt RndRestarts = Env.GetIfArgPrefixInt("-rn:", DEF_RndRestarts, "Random graph restarts of the optimization algorithm");
  const TFlt RndRecompW = Env.GetIfArgPrefixFlt("-rf:", DEF_RndRecompW, "Random graph reestimation of W if relative difference smaller");

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

  // Run ST-group extraction framework
  TGroupSTV GroupV;
  GroupExtractFramework(GroupV, Graph, OptRestarts, OptMxSteps, OptStopSteps, OptInitSample, FinishCnt, FinishRndW, RndGraphs, RndRestarts, RndRecompW);
  ret = !(GroupV.Len() > 0);

  // Output ST-group assignments (*.groups)
  if(OutFNm.Len() > 0) {
    OutputGroups(OutFNm, "wt", argc, argv, InFNm, Graph->GetNodes(), Graph->GetEdges(), NIdLabelH, GroupV);
  }

  // Output ST-group extraction results summary (*.groupssum)
  if(OutSumFNm.Len() > 0) {
    OutputGroupsSum(OutSumFNm, "wt", argc, argv, InFNm, Graph->GetNodes(), Graph->GetEdges(), NIdLabelH, GroupV);
  }

  // Footer
  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return ret;
}
