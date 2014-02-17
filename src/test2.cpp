/**
 * nodegroups - Test2: Analyze random Erdos-Renyi graph probability distribution
 *
 * Example:
 *   ./test2 -dd:../rnd.distw -dn:100 -dm:1000 -rn:100000
 *
 * @author  gw0 [http://gw.tnode.com/] <gw.2014@tnode.com>
 * @version (see group_h_VERSION)
 */

#include "Snap.h"
#include "group.h"


/**
 * Estimate samples of random Erdos-Renyi graphs
 *
 * @param[out] GroupV List of results of node group extraction (into S, T)
 * @param N Number of nodes in random graphs
 * @param M Number of edges in random graphs
 * @param RndRestarts Number of restarts on Erdos-Renyi random graphs
 * @param OptMxSteps Maximal number of steps in each optimization run
 * @param OptStopSteps Stop optimization if no W improvement in steps
 * @param OptInitSample Initial random-sample size of S ant T (0 for random)
 * @return Number of estimated samples
 */
int GroupExtractRndGnms(TGroupSTV& GroupV, int N, int M, int RndRestarts/*=DEF_RndRestarts*/, int OptMxSteps/*=DEF_OptMxSteps*/, int OptStopSteps/*=DEF_OptStopSteps*/, int OptInitSample/*=DEF_OptInitSample*/) {
  PUNGraph GraphER;

  for (int i = 0; i < RndRestarts; ++i) {
    // Generate a Erdos-Renyi random graph
    GraphER = TSnap::GenRndGnm<PUNGraph>(N, M, false);

    // Find node group extraction (into S, T)
    TGroupST R = {};
    GroupExtractSingle(R, GraphER, OptMxSteps, OptStopSteps, OptInitSample);

    // Add to list
    GroupV.AddMerged(R);
  }

  return GroupV.Len();
}


/**
 * Console application entry point
 */
int main(int argc, char* argv[]) {
  // Header
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("test2. Build: %.2f, %s, %s. Time: %s", group_h_VERSION, __TIME__, __DATE__, TExeTm::GetCurTm()), 1);
  TExeTm ExeTm;
  Try

  // Parameters
  const TStr RndDistFNm = Env.GetIfArgPrefixStr("-dd:", "rand.distw", "Output file with random graphs distribution of W for given N and M");
  const TInt RndN = Env.GetIfArgPrefixInt("-dn:", 100, "Number of nodes in random graphs");
  const TInt RndM = Env.GetIfArgPrefixInt("-dm:", 1000, "Number of edges in random graphs");
  const TInt OptMxSteps = Env.GetIfArgPrefixInt("-sm:", DEF_OptMxSteps, "Maximal number of steps in each optimization run");
  const TFlt OptStopSteps = Env.GetIfArgPrefixFlt("-sw:", DEF_OptStopSteps, "Stop optimization if no W improvement in steps");
  const TInt OptInitSample = Env.GetIfArgPrefixInt("-ss:", DEF_OptInitSample, "Initial random-sample size of S ant T (0=random)");
  const TInt RndRestarts = Env.GetIfArgPrefixInt("-rn:", DEF_RndRestarts, "Number of restarts on Erdos-Renyi random graphs");

  // Test2
  TGroupSTV GroupV;
  GroupExtractRndGnms(GroupV, RndN, RndM, RndRestarts, OptMxSteps, OptStopSteps, OptInitSample);

  // Output
  FILE *F = fopen(RndDistFNm.CStr(), "wa");
  fprintf(F, "# RndRestarts: %d\n", GroupV.Len());
  fprintf(F, "# Nodes: %d    Edges: %d\n", GroupV[0].N, GroupV[0].M);
  fprintf(F, "N\tM\tW\tDetails\n");
  for (int j = 0; j < GroupV.Len(); ++j) {
    TGroupST& G = GroupV[j];
    fprintf(F, "%d\t%d\t%.4f\t# %s\n", G.N, G.M, G.W, G.GetStr().CStr());
  }
  fprintf(F, "#\n");
  fclose(F);

  // Footer
  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
