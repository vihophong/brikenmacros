#include <bits/stdc++.h>
#include <chrono>

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TThread.h"
#include "TStopwatch.h"

std::multimap < unsigned long long, int > fbigripsMap;
std::multimap < unsigned long long, int >::iterator fbigripsMap_it;
std::multimap < unsigned long long, int >::iterator fbigripsMap_it2;
std::multimap < unsigned long long, int >::iterator fbigripsMap_it3;

using namespace std;
void maptest()
{
    for (Int_t i=0;i<10000;i++){
        fbigripsMap.insert(make_pair(i,i));
    }
    auto start = chrono::high_resolution_clock::now();
    ios_base::sync_with_stdio(false);
    fbigripsMap_it2 = fbigripsMap.lower_bound(7);
    auto end = chrono::high_resolution_clock::now();

    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *=1e-9;
    cout<<fixed<<time_taken<<setprecision(9)<<endl;
}
