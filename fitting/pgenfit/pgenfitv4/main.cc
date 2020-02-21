#include <iostream>
#include "unbinfit.hh"

int main(int argc, char *argv[])
{
    /*
    unbinfit* fit=new unbinfit;
    char inputRootFile[1000];
    sprintf(inputRootFile,"testdata.root");
    fit->Init(argv[1],inputRootFile);
    fit->generateRoofitEvaluate();
    */
    if (argc==2) {
        unbinfit* fit=new unbinfit;
        char inpparms[1000];
        char inputRootFile[1000];
        sprintf(inpparms,"parmsex.txt");
        sprintf(inputRootFile,"testdata.root");
        //keep default start time (0.08s)
        fit->Init(inpparms,inputRootFile);
        fit->setOutputFile(argv[1]);
        fit->Run();
    }else if(argc==6){
        unbinfit* fit=new unbinfit;
        fit->setStartTime(atof(argv[4]));
        fit->Init(argv[1],argv[2]);
        fit->setOutputFile(argv[3]);
        fit->setNumberOfMC(atoi(argv[5]));
        fit->Run();
    }else{
        std::cout<<"check inputs!"<<std::endl;
        return 0;
    }
    return 0;
}
