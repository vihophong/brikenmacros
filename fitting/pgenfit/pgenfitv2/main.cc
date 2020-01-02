#include <iostream>
#include "unbinfit.hh"

int main(int argc, char *argv[])
{
    unbinfit* fit;
    if (argc==2) {
        fit=new unbinfit;
        char inpparms[1000];
        char outpparms[1000];
        sprintf(inpparms,"parmsex.txt");
        sprintf(outpparms,"testdata.root");
        fit->Init(inpparms,outpparms);
        fit->Run(argv[1]);
    }else if(argc==4){
        fit=new unbinfit;
        fit->Init(argv[1],argv[2]);
        fit->Run(argv[3]);
    }else{
        std::cout<<"check inputs!"<<std::endl;
        return 0;
    }
    delete fit;
    return 0;
}