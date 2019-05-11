#include <fstream>
#include <sstream>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraphErrors.h>

#include <TMarker.h>


void normalize(char* infilename,char* outfilename,Double_t runtime=4601.66 ,Double_t act=7346.5,Double_t act_relerr=1.9,Double_t factor=10)
{
    act=act*runtime;
    Double_t temp=0;
    Int_t k=0;
    std::ifstream ifs(infilename);
    Double_t arr[4];
    std::ofstream ofs(outfilename);



    while (!ifs.eof()){
        ifs>>temp;
        arr[k]=temp;
        k++;
        if (k==4){
            arr[3]=sqrt(arr[3]*arr[3]/arr[2]/arr[2]+act_relerr*act_relerr/10000)*arr[2]*act/100./factor;
            arr[2]=arr[2]*act/100./factor;
            k=0;
            ofs<<std::setw(16)<<arr[0];
            for (Int_t i=1;i<4;i++)
            ofs<<std::setw(12)<<arr[i];
            ofs<<endl;
        }
    }

}
