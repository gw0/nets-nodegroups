/**
 * nodegroups/ex01 - Compute W for given labels in groups S and T
 *
 * Example:
 *   ./ex01 -o:../../data_nets/football -is:"Florida State,North Carolina,Virginia,Duke,Georgia Tech,Clemson,Maryland,North Carolina State,Wake Forest" -it:"Florida State,North Carolina,Virginia,Duke,Clemson,Maryland,Georgia Tech,North Carolina State,Wake Forest"
 *
 * @author  gw0 [http://gw.tnode.com/] <gw.2014@tnode.com>
 * @version (see group_h_VERSION)
 */

#include "Snap.h"
#include "../group.h"


void ComputeWForLabelST(TGroupST& G, const PUNGraph& Graph, const TIntStrH& NIdLabelH, const TStr& LabelSStr, const TStr& LabelTStr) {
  TStr CurStr, TmpStr;

  // Prepare inverse mapping (till first space)
  TStrIntH LabelNIdH;
  for (TIntStrH::TIter I = NIdLabelH.BegI(); I != NIdLabelH.EndI(); I++) {
    I.GetDat().SplitOnCh(CurStr, ' ', TmpStr);
    LabelNIdH.AddDat(CurStr, I.GetKey());
  }

  // Populate group S (match without spaces)
  printf("\nGroup S:\n");
  TStrV LabelSStrV;
  LabelSStr.SplitOnAllCh(',', LabelSStrV);
  for (TStrV::TIter I = LabelSStrV.BegI(); I != LabelSStrV.EndI(); I++) {
    CurStr = *I;
    CurStr.DelChAll(' ');

    G.SubSNIdV.AddMerged(LabelNIdH.GetDat(CurStr));
    printf("%d\t%s\t%s\n", LabelNIdH.GetDat(CurStr).Val, CurStr.CStr(), NIdLabelH.GetDat(LabelNIdH.GetDat(CurStr).Val).CStr());
  }

  // Populate group T (match without spaces)
  printf("\nGroup T:\n");
  TStrV LabelTStrV;
  LabelTStr.SplitOnAllCh(',', LabelTStrV);
  for (TStrV::TIter I = LabelTStrV.BegI(); I != LabelTStrV.EndI(); I++) {
    CurStr = *I;
    CurStr.DelChAll(' ');

    G.SubTNIdV.AddMerged(LabelNIdH.GetDat(CurStr));
    printf("%d\t%s\t%s\n", LabelNIdH.GetDat(CurStr).Val, CurStr.CStr(), NIdLabelH.GetDat(LabelNIdH.GetDat(CurStr).Val).CStr());
  }

  // Populate all ST-group extraction results
  G.RecomputeAll(Graph, G.SubSNIdV, G.SubTNIdV);

  // Recompute on corresponding Erdos-Renyi random graphs
  TGroupST R;
  TGroupSTV GroupERV;
  GroupExtractRndGnms(GroupERV, G.N, G.M);
  R = GroupERV[(int)(GroupERV.Len() * (100.0 - DEF_FinishRndW) / 100.0)];
  if (R.W < 0.0)
    R.W = 0.0;

  // Print status
  printf("\n");
  printf("foo %s\n", G.GetStr().CStr());
  printf("  r %s\n", R.GetStr().CStr());
}


/**
 * Console application entry point
 */
int main(int argc, char* argv[]) {
  // Header
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("ex01. Build: %.2f, %s, %s. Time: %s", group_h_VERSION, __TIME__, __DATE__, TExeTm::GetCurTm()), 1);
  TExeTm ExeTm;
  Try

  // Parameters
  const TStr PrefixFNm = Env.GetIfArgPrefixStr("-o:", "graph", "Input and output file name prefix (can be overriden)");
  const TStr InFNm = Env.GetIfArgPrefixStr("-i:", PrefixFNm + ".edgelist", "Input graph edges (undirected edge per line)");
  const TStr LabelFNm = Env.GetIfArgPrefixStr("-l:", PrefixFNm + ".labels", "Optional input node labels (node ID, node label)");
  const TStr OutFNm = Env.GetIfArgPrefixStr("-og:", PrefixFNm + ".groups", "Output group assignments (for S and T)");
  const TStr LabelSStr = Env.GetIfArgPrefixStr("-is:", "", "Comma separated labels for group S");
  const TStr LabelTStr = Env.GetIfArgPrefixStr("-it:", "", "Comma separated labels for group T");

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

  // Experiment
  TGroupSTV GroupV;
  TGroupST G = {};
  ComputeWForLabelST(G, Graph, NIdLabelH, LabelSStr, LabelTStr);
  GroupV.Add(G);

  // Output
  FILE *F = fopen(OutFNm.CStr(), "wt");
  fprintf(F, "# ex01. Build: %.2f, %s, %s. Time: %s\n", group_h_VERSION, __TIME__, __DATE__, TExeTm::GetCurTm());
  fprintf(F, "#");
  for (int i = 0; i < argc; ++i) {
    fprintf(F, " %s", argv[i]);
  }
  fprintf(F, "\n");
  fprintf(F, "# Input: %s  Nodes: %d  Edges: %d\n", InFNm.CStr(), GroupV[0].N, GroupV[0].M);
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
