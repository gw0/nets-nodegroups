/**
 * nodegroups/ex02 - Analyze random Erdos-Renyi graph probability distribution
 *
 * Example:
 *   ./ex02 -dd:ex02_lwcc_n100_m1000.distw -dn:100 -dm:1000 -rg:100000 -rn:1
 *
 * @author  gw0 [http://gw.tnode.com/] <gw.2014@tnode.com>
 * @version (see group_h_VERSION)
 */

#include "Snap.h"
#include "../group.h"


/**
 * Estimate samples of random Erdos-Renyi graphs
 *
 * @param[out] GroupV List of results of node group extraction (into S, T)
 * @param N Number of nodes in random graphs
 * @param M Number of edges in random graphs
 * @param RndGraphs Number of different Erdos-Renyi random graphs
 * @param RndRestarts Number of restarts on Erdos-Renyi random graphs
 * @param OptMxSteps Maximal number of steps in each optimization run
 * @param OptStopSteps Stop optimization if no W improvement in steps
 * @param OptInitSample Initial random-sample size of S ant T (0 for random)
 * @return Number of estimated samples
 */
int GroupExtractRndGnms(TGroupSTV& GroupV, int N, int M, int RndGraphs, int OptMxSteps, int OptStopSteps, int OptInitSample, int RndRestarts, bool RndMxWcc) {
  PUNGraph GraphER;

  for (int i = 0; i < RndGraphs; ++i) {
    // Generate a Erdos-Renyi random graph
    GraphER = TSnap::GenRndGnm<PUNGraph>(N, M, false);
    if (RndMxWcc) {  // largest weakly-connected component
      GraphER = TSnap::GetMxWcc(GraphER);
    }

    // Find node group extraction (into S, T)
    TGroupST R = {};
    GroupExtractRestarter(R, GraphER, RndRestarts, OptMxSteps, OptStopSteps, OptInitSample);

    // Add to list
    GroupV.AddSorted(R);
  }

  return GroupV.Len();
}


/**
 * Console application entry point
 */
int main(int argc, char* argv[]) {
  // Header
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("ex02. Build: %.2f, %s, %s. Time: %s", group_h_VERSION, __TIME__, __DATE__, TExeTm::GetCurTm()), 1);
  TExeTm ExeTm;
  Try

  // Parameters
  const TStr RndDistFNm = Env.GetIfArgPrefixStr("-dd:", "", "Output file with random graphs distribution of W for given N and M");
  const TInt RndN = Env.GetIfArgPrefixInt("-dn:", 100, "Number of nodes in random graphs");
  const TInt RndM = Env.GetIfArgPrefixInt("-dm:", 1000, "Number of edges in random graphs");
  const TInt OptMxSteps = Env.GetIfArgPrefixInt("-sm:", DEF_OptMxSteps, "Maximal number of steps in each optimization run");
  const TInt OptStopSteps = Env.GetIfArgPrefixInt("-sw:", DEF_OptStopSteps, "Stop optimization if no W improvement in steps");
  const TInt OptInitSample = Env.GetIfArgPrefixInt("-ss:", DEF_OptInitSample, "Initial random-sample size of S ant T (0=random)");
  const TInt RndGraphs = Env.GetIfArgPrefixInt("-rg:", DEF_RndGraphs, "Number of different Erdos-Renyi random graphs");
  const TInt RndRestarts = Env.GetIfArgPrefixInt("-rn:", DEF_RndRestarts, "Number of restarts on each Erdos-Renyi random graph");
  const TBool RndMxWcc = Env.GetIfArgPrefixBool("-rc:", true, "Always strip to largest weakly-connected component");
  if (RndDistFNm.Empty())
    exit(-1);

  // Experiment
  TGroupSTV GroupV;
  GroupExtractRndGnms(GroupV, RndN, RndM, RndGraphs, OptMxSteps, OptStopSteps, OptInitSample, RndRestarts, RndMxWcc);

  // Output
  FILE *F = fopen(RndDistFNm.CStr(), "wa");
  fprintf(F, "# ex02. Build: %.2f, %s, %s. Time: %s\n", group_h_VERSION, __TIME__, __DATE__, TExeTm::GetCurTm());
  fprintf(F, "#");
  for (int i = 0; i < argc; ++i) {
    fprintf(F, " %s", argv[i]);
  }
  fprintf(F, "\n");
  fprintf(F, "# Graphs: %d  Nodes: %d  Edges: %d\n", GroupV.Len(), GroupV[0].N, GroupV[0].M);
  fprintf(F, "N\tM\tN_S\tM_S\tN_T\tM_T\tN_ST\tM_ST\tL_ST\tL_STc\tW\tTau\tMod_S\tMod_T\n");
  for (int j = 0; j < GroupV.Len(); ++j) {
    TGroupST& G = GroupV[j];
    fprintf(F, "%s\n", G.GetStr(21).CStr());
  }
  fprintf(F, "#\n");
  fclose(F);

  // Footer
  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
