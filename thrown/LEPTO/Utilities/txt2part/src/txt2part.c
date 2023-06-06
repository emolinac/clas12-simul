/* -=======================================================- *
$Id: txt2part.c,v 1.26 2006/09/17 01:03:01 avto Exp $
$Author: avto $
$Revision: 1.26 $
$Date: 2006/09/17 01:03:01 $
* -=======================================================- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <errno.h>
#include <time.h>
#include "ntypes.h"
#include "bostypes.h"
#include "pid.h"
#include "particleType.h"
#include "pdgutil.h"
#include "scalers.h"

/* declare the bos common */
BOSbank bcs_;
BOSbank wcs_;

void PrintUsage(char *processName);
mchd_t *make_mchd(head_t head, int reactype, float Sigma, float w, float q2, vector4_t photon, vector4_t target);

int main(int argc, char *argv[]) {

  int OutputUnitNo = 7;
  int MaxBanks = 1000;
  int i;
  int group = 0;
  int NEvnt = 0;
  int NPart;
  int DiscardLines = 0;
  int MakeMCTAGR = 0;
  int MakeMCTK = 0;
  int MakePART = 1;
  int use_energy = 0;

  float Charge[NUM_PARTICLES];
  char *argptr;
  char *BOSoutfile = "part.evt";
  char ChJunk[200];
  char out[300];
  clasPART_t *PART = NULL;
  clasMCTK_t *MCTK = NULL;
  clasMCVX_t *MCVX = NULL;
  part_t part[PART_SIZE];
  mctk_t mctk[PART_SIZE];
  mcvx_t mcvx[PART_SIZE];
  head_t head[30];
  clasHEAD_t *HEAD = NULL;
  float Egamma, Eelectron, toffset; /*stuff for making TAGR*/
  int RunNo = 0;
  int Read_Sigma = 0;

  int reactype = 1;
  float Sigma;
  float w = -9999.;
  float q2 = -9999.;
  vector4_t phot = {-9999., -9999., -9999., -9999.};
  vector4_t targ = {-9999., -9999., -9999., -9999};

  for (i = 1; i < argc; ++i) {
    argptr = argv[i];
    if (*(argptr = argv[i]) == '-') {

      switch (*(++argptr)) {
        case 'h':
          PrintUsage(argv[0]);
          exit(1);
          break;
        case 'o':
          if (*(++argptr))
            BOSoutfile = argptr;
          else
            PrintUsage(argv[0]);
          break;
        case 'd':
          DiscardLines = atoi(++argptr);
          break;
        case 'm':
          MakePART = 0;
          MakeMCTK = 1;
          break;
        case 'p':
          MakePART = 1;
          break;
        case 'T':
          MakeMCTAGR = 1;
          break;
        case 't':
          MakeMCTAGR = 1;
          use_energy = 1;
          break;
        case 'P':
          group = atoi(++argptr);
          break;
        case 'R':
          RunNo = atoi(++argptr);
          break;
        case 'x':
          Read_Sigma = 1;
          break;
        default:
          PrintUsage(argv[0]);
          exit(1);
          break;
      }
    }
  }

  if (OutputUnitNo) {
    fprintf(stderr, "Output file: %s\n", BOSoutfile);
    unlink(BOSoutfile);
    sprintf(out, "OPEN AUTOOUTPUT UNIT=7 FILE=\"%s\" READWRITE STATUS=NEW RECL=3600", BOSoutfile);
    if (!fparm_c(out)) {
      fprintf(stderr, "%s: Unable to open file \'%s\': %s\n\n", argv[0], out, strerror(errno));
      exit(1);
    }
  }

  bnames_(&MaxBanks);
  initbos();
  configure_banks(stderr, 0);

  bankList(&bcs_, "E=", "HEADPARTTAGRMCHDMCTKMCVX");

  /* get rid of the first three lines of GenBod spew */
  for (i = 0; i < DiscardLines; i++) {
    fgets(ChJunk, 200, stdin);
  }
  if (MakeMCTAGR) { /*set up TAGR stuff*/
    fill_ebin_map();
    fill_tbin_map();
    fill_E_T_map();
  }

  while (fscanf(stdin, "%i", &NPart) == 1) {
    memset((void *)&(part[0]), 0, PART_SIZE * sizeof(part_t));
    memset((void *)&(mctk[0]), 0, PART_SIZE * sizeof(mctk_t));
    memset((void *)&(mcvx[0]), 0, PART_SIZE * sizeof(mcvx_t));
    memset((void *)&(head[0]), 0, 30 * sizeof(head_t));

    if (Read_Sigma) {
      if (fscanf(stdin, "%e", &Sigma) == 1) {
      }

      else {
        fprintf(stderr, "txt2part: has encountered a format input error in event %d.\n,", NEvnt + 1);
        exit(1);
      }
    }

    if (MakeMCTAGR) {
      if (fscanf(stdin, "%f %f %g", &Eelectron, &Egamma, &toffset) == 3) {
        make_mctagr(&bcs_, Eelectron, Egamma, toffset, use_energy);
      } else {
        fprintf(stderr, "A: txt2part: has encountered a format input error in event %d.\n,", NEvnt + 1);
        exit(1);
      }
    }

    if (MakePART) {
      for (i = 0; i < NPart; i++) {
        if (fscanf(stdin, "%i %f %f %f %f", &part[i].pid, &part[i].p.t, &part[i].p.space.x, &part[i].p.space.y, &part[i].p.space.z) == 5 &&
            fscanf(stdin, "%f %f %f %f", &part[i].vert.x, &part[i].vert.y, &part[i].vert.z, &part[i].qpid) == 4) {
          part[i].q = gchrg(part[i].pid);
        } else {
          fprintf(stderr, "B: txt2part: has encountered a format input error in event %d.\n,", NEvnt + 1);
          exit(1);
        }
      }
    }

    /*------------- additional part  for MCTK bank A.T. --------------*/

    if (MakeMCTK) {
      for (i = 0; i < NPart; i++) {

        if (fscanf(stdin, "%i %f %f %f %f", &mctk[i].id, &mctk[i].cx, &mctk[i].cy, &mctk[i].cz, &mctk[i].pmom) == 5 &&
            fscanf(stdin, "%f %f", &mctk[i].mass, &mctk[i].charge) == 2 &&
            fscanf(stdin, "%f %f %f %f %i", &mcvx[i].x, &mcvx[i].y, &mcvx[i].z, &mcvx[i].tof, &mcvx[i].flag) == 5) {
          mctk[i].beg_vtx = 1;
          mctk[i].end_vtx = 0;
          mctk[i].flag = 5;
          if (mctk[i].id == 22) {
            mctk[i].beg_vtx = 0;
            mctk[i].end_vtx = 1;
            mctk[i].flag = 16;
          }

        } else {
          fprintf(stderr, "txt2part: has encountered a format input error for MCTK Bank in event %d.\n,", NEvnt + 1);
          exit(1);
        }
      }
    }

    if (NPart) {
      if ((HEAD = makeBank(&bcs_, "HEAD", 0, sizeof(head_t) / 4, 1))) {
        time_t secs;
        time(&secs);
        head[0].nevent = ++NEvnt;

        if (RunNo == 0) RunNo = 1;

        head[0].nrun = RunNo;
        head[0].evtclass = 15;

        head[0].type = -4;
        head[0].time = secs;
        head[0].trigbits = 1;
        HEAD->head[0] = head[0];
        if (MakePART && (PART = makeBank(&bcs_, "PART", group, sizeof(part_t) / 4, NPart))) {
          for (i = 0; i < NPart; i++) {
            PART->part[i] = part[i];
          }
        }
        if (MakeMCTK && (MCTK = makeBank(&bcs_, "MCTK", group, sizeof(mctk_t) / 4, NPart))) {
          for (i = 0; i < NPart; i++) {
            MCTK->mctk[i] = mctk[i];
          }
        }
        if (MakeMCTK && (MCVX = makeBank(&bcs_, "MCVX", group, sizeof(mcvx_t) / 4, NPart))) {
          for (i = 0; i < NPart; i++) {
            MCVX->mcvx[i] = mcvx[i];
          }
        }
        if (Read_Sigma) {
          make_mchd(HEAD->head[0], reactype, Sigma, w, q2, phot, targ);
        }
      }
      putBOS(&bcs_, OutputUnitNo, "E");
    }
    /* tidy up */
    dropAllBanks(&bcs_, "E");
    cleanBanks(&bcs_);
  }
  sprintf(out, "CLOSE AUTOOUTPUT UNIT=7");
  fparm_c(out);

  return 0;
}

void PrintUsage(char *processName) {
  fprintf(stderr, "   \n");
  fprintf(stderr, "%s: text to PART BOS bank converter.\n", processName);
  fprintf(stderr, "input:\tstdin text file of particle 4-vectors\n");
  fprintf(stderr, "output:\tBOS file with PART or MCTK bank\n");
  fprintf(stderr, "options:\n");
  fprintf(stderr, "\t-m            \tmake MCTK bank \n");
  fprintf(stderr, "\t-p            \tmake PART bank \n");
  fprintf(stderr, "\t-h          \tprint this Help message.\n");
  fprintf(stderr, "\t-o<outfile> \tOutput BOS file name.\n");
  fprintf(stderr, "\t-d<int>     \tDiscard the first n lines.\n");
  fprintf(stderr, "\t-T          \tmake TAGR bank: requires additional txt\n");
  fprintf(stderr, "\t-t          \tsame as -T except do not simulate\n");
  fprintf(stderr, "\t            \tthe effect of e-bins for egamma\n");
  fprintf(stderr, "\t-P#         \tspecify sector # in PART bank\n");
  fprintf(stderr, "\t-x          \tmake partial MCHD bank with the CrossSection specified below\n");
  fprintf(stderr, "\t-R<int>     \tmake HEAD bank for a given run-number\n");

  fprintf(stderr, "example:\n");
  fprintf(stderr, "\ttxt2part -oEventGenerator.evt < EventGenerator.dat\n");
  fprintf(stderr, "The format of the input text data is as follows:\n");
  fprintf(stderr, "\t --- For PART BANK ---:\n");
  fprintf(stderr, "\tN\n");
  fprintf(stderr, "\tCrossSection      /* only used if -x option chosen*/\n");
  fprintf(stderr, "\tE_electron, E_gamma, T_offset  /*only used if -T or -t option chosen*/\n");
  fprintf(stderr, "\tpid Pt Px Py Pz\n");
  fprintf(stderr, "\tVx Vy Vz\n");
  fprintf(stderr, "\t      ...\n");
  fprintf(stderr,
          "Where N is the number of particles in an event, pid is the geant particle id, \nPt is the time component of the four-momentum, "
          "and Px, Py, Pz are the space \ncomponents.  Vx, Vy, and Vz are the vertex location of the particle.\nE_electon is the incident "
          "beam energy on the radiator, E_photon is the energy\n of the photon, and T_offset is the time offset due to the event vertex "
          "not \nbeing at the origin (ie toffset = -z/c).  E_electron,E_gamma, and T_offset \nare only used if the -T option is used.\n\n");

  fprintf(stderr, "\t --- For MCTK BANK ---:\n");
  fprintf(stderr, "\tN\n");
  fprintf(stderr, "\tCrossSection      /* only used if -x option chosen*/\n");
  fprintf(stderr, "\tPDG_id cx cy cz momentum\n");
  fprintf(stderr, "\tmass charge\n");
  fprintf(stderr, "\tVx Vy Vz T_offset flag\n\n");

  fprintf(stderr,
          "Where N is the number of particles in an event, PDG_id is the PDG particle id, \n cx cy cz cosines of 3-momentum, momentum is "
          "3-momentum.\n  Vx, Vy, and Vz are the vertex location of the particle.\n T_offset is the time offset due to the event vertex "
          "not \nbeing at the origin (ie toffset = -z/c), and falg is vertex flag\n");
}

mchd_t *make_mchd(head_t head, int reactype, float Sigma, float w, float q2, vector4_t photon, vector4_t target) {
  clasMCHD_t *MCHD = NULL;
  if ((MCHD = makeBank(&bcs_, "MCHD", 0, sizeof(mchd_t) / 4, 1))) {
    MCHD->mchd[0].nrun = head.nrun;
    MCHD->mchd[0].nevent = head.nevent;
    MCHD->mchd[0].time = head.time;
    MCHD->mchd[0].type = head.type;
    MCHD->mchd[0].reactype = reactype;
    MCHD->mchd[0].weight = Sigma;
    MCHD->mchd[0].w = w;
    MCHD->mchd[0].q2 = q2;
    MCHD->mchd[0].w = w;
    MCHD->mchd[0].e_phot = photon.t;
    MCHD->mchd[0].px_phot = photon.space.x;
    MCHD->mchd[0].py_phot = photon.space.y;
    MCHD->mchd[0].pz_phot = photon.space.z;
    MCHD->mchd[0].e_targ = target.t;
    MCHD->mchd[0].px_targ = target.space.x;
    MCHD->mchd[0].py_targ = target.space.y;
    MCHD->mchd[0].pz_targ = target.space.z;
    return &(MCHD->mchd[0]);
  }
  return NULL;
}
