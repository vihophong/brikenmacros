#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include <time.h>
#include <unistd.h>
#include <math.h>
#include <sys/stat.h>


#include "TTree.h"
#include "TFile.h"

#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"

#include "TSpectrum.h"
#include "TF1.h"


using namespace std;

int file_error(char *error_type, char *filename)
{
  /* write error message */
  /* cannot perform operation error_type on file filename */

  if (strlen(error_type) + strlen(filename) > 58) {
    printf("ERROR - cannot %s file\n%s\n", error_type, filename);
  } else {
    printf("ERROR - cannot %s file %s\n", error_type, filename);
  }
  return 0;
} /* file_error */


int setext(char *filnam, const char *cext, int filnam_len)
{
  /* set default extension of filename filnam to cext
     leading spaces are first removed from filnam
     if extension is present, it is left unchanged
     if no extension is present, cext is used
     returned value pointer to the dot of the .ext
     cext should include the dot plus a three-letter extension */

  int nc, iext;

  /* remove leading spaces from filnam */
  nc = strlen(filnam);
  if (nc > filnam_len) nc = filnam_len;
  while (nc > 0 && filnam[0] == ' ') {
    memmove(filnam, filnam+1, nc--);
    filnam[nc] = '\0';
  }
  /* remove trailing spaces from filnam */
  while (nc > 0 && filnam[nc-1] == ' ') {
    filnam[--nc] = '\0';
  }
  /* look for file extension in filnam
     if there is none, put it to cext */
  iext = 0;
  if (nc > 0) {
    for (iext = nc-1;
         (iext > 0 &&
          filnam[iext] != ']' &&
          filnam[iext] != '/' &&
          filnam[iext] != ':');
         iext--) {
      if (filnam[iext] == '.') return iext;
    }
    iext = nc;
  }
  strncpy(&filnam[iext], cext, filnam_len - iext - 1);
  return iext;
} /* setext */

int inq_file(char *filename, int *reclen)
{
  /* inquire for file existence and record length in longwords
     returns 0 for file not exists, 1 for file exists */

  int  ext;
  char jfile[80];
  struct stat statbuf;

  *reclen = 0;
  if (stat(filename, &statbuf)) return 0;

  ext = 0;
  strncpy(jfile, filename, 80);
  ext = setext(jfile, "    ", 80);
  if (!strcmp(&jfile[ext], ".mat") ||
      !strcmp(&jfile[ext], ".MAT") ||
      !strcmp(&jfile[ext], ".esc") ||
      !strcmp(&jfile[ext], ".ESC")) {
    *reclen = 2048;
  } else if (!strcmp(&jfile[ext], ".spn") ||
             !strcmp(&jfile[ext], ".SPN") ||
             !strcmp(&jfile[ext], ".m4b") ||
             !strcmp(&jfile[ext], ".M4B") ||
             !strcmp(&jfile[ext], ".e4k") ||
             !strcmp(&jfile[ext], ".E4K")) {
    *reclen = 4096;
  } else if (!strcmp(&jfile[ext], ".cub") ||
             !strcmp(&jfile[ext], ".CUB")) {
    *reclen = 256;
  } else if (!strcmp(&jfile[ext], ".2dp") ||
             !strcmp(&jfile[ext], ".2DP")) {
    if (statbuf.st_size <= 0) {
      *reclen = 0;
    } else {
      *reclen = (int) (0.5 + sqrt((float) (statbuf.st_size/4)));
    }
  }
  return 1;
} /* inq_file */



FILE *open_new_file(char *filename, int force_open)
{
  /* safely open a new file
     filename: name of file to be opened
     force_open = 0 : allow return value NULL for no file opened
     force_open = 1 : require that a file be opened */

  int  j, nc, jext, fn_error = 0;
  char tfn[80], *owf;
  FILE *file = NULL;
  file = fopen(filename, "w+");
  return file;
} /* open_new_file */


int put_file_rec(FILE *fd, void *data, int numbytes)
{
  /* write one fortran-unformatted style binary record into data */
  /* returns 1 for error */

#ifdef VMS  /* vms */
  int   j1;
  short rh[2];
  char  *buf;

  buf = data;
  j1 = numbytes;
  if (numbytes <= 2042) {
    rh[0] = numbytes + 2; rh[1] = 3;
  } else {
    rh[0] = 2044; rh[1] = 1;
    while (j1 > 2042) {
      if (fwrite(rh, 2, 2, fd) != 2 ||
          fwrite(buf, 2042, 1, fd) != 1) return 1;
       rh[1] = 0; j1 -= 2042; buf += 2042;
    }
    rh[0] = j1 + 2; rh[1] = 2;
  }
  if (fwrite(rh, 2, 2, fd) != 2 ||
      fwrite(buf, j1, 1, fd) != 1) return 1;
  /* if numbytes is odd, write an extra (padding) byte */
  if (2*(numbytes>>1) != numbytes) {
    j1 = 0;
    fwrite(&j1, 1, 1, fd);
  }

#else /* unix */

  if (fwrite(&numbytes, 4, 1, fd) != 1 ||
      fwrite(data, numbytes, 1, fd) != 1 ||
      fwrite(&numbytes, 4, 1, fd) != 1) return 1;
#endif
  return 0;
} /*put_file_rec */


int wspec(char *filnam, float *spec, int idim)
{
  /* write spectra in gf3 format
     filnam = name of file to be created and written
     spec = spectrum of length idim */

  char buf[32];
  int  j, c1 = 1, rl = 0;
  char namesp[8];
  FILE *file;

  j = setext(filnam, ".spe", 80);

  if (!(file = open_new_file(filnam, 0))) return 1;

  strncpy(namesp, filnam, 8);
  if (j < 8) memset(&namesp[j], ' ', 8-j);


  /* WRITE(1) NAMESP,IDIM,1,1,1 */
  /* WRITE(1) SPEC */
#define W(a,b) { memcpy(buf + rl, a, b); rl += b; }
  W(namesp,8); W(&idim,4); W(&c1,4); W(&c1,4); W(&c1,4);
#undef W
  if (put_file_rec(file, buf, rl) ||
      put_file_rec(file, spec, 4*idim)) {
    printf("write to %s", filnam);
    fclose(file);
    return 1;
  }
  fclose(file);
  return 0;
} /* wspec */


void tospe(char* infile, char* histname, char* outfile)
{

    TFile *f1 = TFile::Open(infile);
    TH1F* hin=(TH1F*) f1->Get(histname);

//    std::ofstream ofs(outfile);

    Int_t k=0;
    float spec[hin->GetNbinsX()];
    //ofs<<histname<<"\t"<<hin->GetNbinsX()<<"\t1\t1\t1"<<endl;
    for (Int_t i=0;i< hin->GetNbinsX();i++){
        //ofs<<setw(10)<<hin->GetBinContent(i+1);
        spec[i]=hin->GetBinContent(i+1);
        k++;
        if (k>7){
            k=0;
            //ofs<<endl;
        }
    }
    wspec(outfile,spec,(int) hin->GetNbinsX());
    //ofs<<endl;
    //ofs.close();
}
