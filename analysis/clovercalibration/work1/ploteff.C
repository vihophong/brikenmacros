#include <fstream>
#include <sstream>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLatex.h>

#include <TMarker.h>


Double_t fnc_eff(Double_t *x, Double_t *pars) {

    double f, g, r, f1, f2, f3, x1, x2, y1, y2, y3;
    double dg = 0.0, d1, df1 = 0.0, df2 = 0.0;

      /* this eval is for use with effit */
      /* calculate the fit using present values of the pars */

      x1 = log(x[0] / pars[7]);
      x2 = log(x[0] / pars[8]);
      g  = pars[6];
      f1 = pars[0] + pars[1]*x1 + pars[2]*x1*x1;
      f2 = pars[3] + pars[4]*x2 + pars[5]*x2*x2;
      if (f1 <= f2) {
        f = f1;
        r = f1 / f2;
      } else {
        f = f2;
        r = f2 / f1;
      }
      if (r <= 1e-6) {
        df1 = exp(f);
        return pars[9]*exp(f);

      } else {
        y1 = pow(r, g);
        y2 = y1 + 1.0;
        d1 = -1.0 / g;
        y3 = pow(y2, d1);
        f3 = exp(f * y3);
        return pars[9]*f3;
      }
}

TGraphErrors* gr;
TF1* efff;

void caleff(char* infilename,char* sinname,char* sinname2,Double_t norm=1.)
{
    std::ifstream infile(infilename);
    std::string line;
    Int_t k=0;
    Int_t j=0;
    Double_t arr[10];



    while (std::getline(infile, line))
    {
        std::istringstream iss(line);

        if (k>0){
            for (Int_t i=0;i<5;i++){
                if (j<10){
                    iss >> arr[j];
                    j++;
                }
            }
        }
        k++;
    }


    efff=new TF1("fnc_eff",fnc_eff,1,2000,10);
    for (Int_t i=0;i<9;i++){
        efff->SetParameter(i,arr[i]);
    }
    efff->SetParError(3,arr[9]);
    //! normalized factor
    efff->FixParameter(9,norm);

    efff->SetNpx(2000);


    Double_t eff,err;
    Double_t n,c,dc,a,da,e,de,intens,dint;
    Int_t np=0;
    Double_t en[1000];
    Double_t effi[1000];
    Double_t effierr[1000];

    std::ifstream sinfile(sinname);
    k=0;
    while (std::getline(sinfile, line))
    {
        std::istringstream iss(line);
        if (k>0){
            iss>>n>>c>>dc>>a>>da>>e>>de>>intens>>dint;
            if (a == 0.0f) continue;
            if (intens == 0.0f) break;
            eff = a / intens;
            err = eff * sqrt(da*da/(a*a) + dint*dint/(intens*intens));

            eff=eff*norm;
            err=err*norm;

            cout<<e<<"\t"<<eff<<"\t"<<err<<endl;
            en[np]=e;
            effi[np]=eff;
            effierr[np]=err;
            np++;
        }
        k++;
    }

    std::ifstream sinfile2(sinname2);
    k=0;
    while (std::getline(sinfile2, line))
    {
        std::istringstream iss(line);
        if (k>0){
            iss>>n>>c>>dc>>a>>da>>e>>de>>intens>>dint;
            if (a == 0.0f) continue;
            if (intens == 0.0f) break;
            eff = a / intens;
            err = eff * sqrt(da*da/(a*a) + dint*dint/(intens*intens));

            eff=eff*norm;
            err=err*norm;

            cout<<e<<"\t"<<eff<<"\t"<<err<<endl;
            en[np]=e;
            effi[np]=eff;
            effierr[np]=err;
            np++;
        }
        k++;
    }

    gr=new TGraphErrors(np,en,effi,0,effierr);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(1);
    gr->SetMarkerSize(1);
    gr->SetLineColor(1);
    gr->SetLineWidth(1);
}

void ploteffsingle(Int_t nch,Double_t norm=1.){
    char tempchar1[1000];
    sprintf(tempchar1,"spec%i.aef",nch);
    char tempchar2[1000];
    sprintf(tempchar2,"spec%i.sin",nch);
    char tempchar3[1000];
    sprintf(tempchar3,"spec%iba.sin",nch);

    caleff(tempchar1,tempchar2,tempchar3,norm);
    cout<<tempchar1<<endl;
    cout<<tempchar2<<endl;

    TCanvas* c1=new TCanvas("c1","c1",900,700);

    gr->Draw("AP");
    efff->Draw("same");
}


void ploteffall(Double_t inpute=1000.){

    Double_t norm[8];
    Double_t relerreff[8];
    TGraphErrors* grall[8];
    TF1* efffall[8];

    for (Int_t i=0;i<8;i++){
        norm[i]=0.1;
        char tempchar1[1000];
        sprintf(tempchar1,"spec%i.aef",i+1);
        char tempchar2[1000];
        sprintf(tempchar2,"spec%i.sin",i+1);
        char tempchar3[1000];
        sprintf(tempchar3,"spec%iba.sin",i+1);
        caleff(tempchar1,tempchar2,tempchar3,norm[i]);
        grall[i]=(TGraphErrors*)gr->Clone();
        efffall[i]=(TF1*)efff->Clone();
    }
    for (Int_t i=0;i<8;i++) {
        relerreff[i]=exp(efffall[i]->GetParError(3))-1;
        cout<<"Percentage of Efficiency error on channel "<<i<<" ="<<relerreff[i]*100<<endl;
    }


    TCanvas* d4c=new TCanvas("d4","d4",900,700);
    d4c->Divide(2,2);
    for (Int_t i=0;i<4;i++){
        d4c->cd(i+1);
        efffall[i]->SetLineColor(i+1);
        if (i==0) grall[i]->SetTitle(Form("D4 black"));
        else if (i==1) grall[i]->SetTitle(Form("D4 red"));
        else if (i==2) grall[i]->SetTitle(Form("D4 green"));
        else if (i==3) grall[i]->SetTitle(Form("D4 blue"));
        grall[i]->Draw("AP");
        efffall[i]->Draw("same");
    }

    TCanvas* g7c=new TCanvas("g7","g7",900,700);
    g7c->Divide(2,2);
    for (Int_t i=4;i<8;i++){
        g7c->cd(i-4+1);
        efffall[i]->SetLineColor(i-4+1);
        if (i-4==0) grall[i]->SetTitle(Form("G7 black"));
        else if (i-4==1) grall[i]->SetTitle(Form("G7 red"));
        else if (i-4==2) grall[i]->SetTitle(Form("G7 green"));
        else if (i-4==3) grall[i]->SetTitle(Form("G7 blue"));
        grall[i]->Draw("AP");
        efffall[i]->Draw("same");
    }

    Double_t sume[2000];
    Double_t sumeff[2000];
    for (Int_t i=1;i<2001;i++){
        sume[i]=i;
        sumeff[i]=0;
        for (Int_t j=0;j<8;j++){
            sumeff[i]+=efffall[j]->Eval(i);
        }
    }

    TGraph* fcnsum=new TGraph(2000,sume,sumeff);
    TCanvas* c2=new TCanvas("c2","c2",900,700);
    c2->cd();
    fcnsum->SetLineWidth(3);
    fcnsum->SetLineColor(6);
    fcnsum->Draw("AL");


    Double_t pntsum[1000];
    Double_t pntsumeff[1000];
    Double_t pntsumefferr[1000];


    Double_t peaklib[]={121.783,
                        244.692,
                        344.276,
                        367.789,
                        411.115,
                        443.976,
                        688.678,
                        778.903,
                        867.388,
                        964.131,
                        1212.95,
                        1299.12,
                        1408.01,
                        53.156,
                        80.311,
                        276.404,
                        302.858,
                        356.014,
                        383.859
                       };

    Int_t np=0;
    for (Int_t m=0;m<19;m++){
        pntsumeff[m]=0;
        pntsumefferr[m]=0;
    }
    for (Int_t j=0;j<8;j++){
        Double_t* x=grall[j]->GetX();
        Double_t* y=grall[j]->GetY();
        Double_t* yerr=grall[j]->GetEY();
        np=0;
        for (Int_t i=0;i<grall[j]->GetN();i++){
            for (Int_t m=0;m<19;m++){
                if (x[i]==peaklib[m]){
                    pntsum[np]=x[i];
                    pntsumeff[np]+=y[i];
                    pntsumefferr[np]+=yerr[i]*yerr[i];
                    np++;
                }
            }
        }
    }

    cout<<" Summary"<<endl;
    for (Int_t m=0;m<19;m++){
        pntsumefferr[m]=sqrt(pntsumefferr[m]);
        cout<<pntsum[m]<<"\t"<<pntsumeff[m]<<"\t"<<pntsumefferr[m]<<endl;
    }

    TGraphErrors* grpntsum=new TGraphErrors(19,pntsum,pntsumeff,0,pntsumefferr);
    grpntsum->SetMarkerSize(1.2);
    grpntsum->SetMarkerStyle(20);
    grpntsum->SetLineColor(1);
    grpntsum->SetLineWidth(1);

    grpntsum->Draw("P same");

    Double_t effeval=0;
    Double_t deffeval=0;

    for (Int_t j=0;j<8;j++){
        Double_t abseff=efffall[j]->Eval(inpute);
        deffeval+=abseff*relerreff[j]*abseff*relerreff[j];
        effeval+=abseff;
    }
    deffeval=sqrt(deffeval);


    TMarker* mr=new TMarker(inpute,effeval,1);
    mr->SetMarkerSize(2.);
    mr->SetMarkerStyle(22);
    mr->SetMarkerColor(3);
    mr->Draw("same");

    TLatex* tag=new TLatex(inpute+10.,effeval,Form("%10.2f keV, %.2f %% #pm %.2f",inpute,effeval,deffeval));
    tag->SetTextAngle(90);
    tag->SetTextSize(0.025);
    tag->Draw("same");


    cout<<"efficiency at "<<inpute<<" is "<<effeval<<" +/- "<<deffeval<<endl;

    TFile* f0=new TFile("hallgccurve.root","recreate");
    grpntsum->SetName("hgc_exp");
    grpntsum->Write();
    fcnsum->SetName("fgc_fit");
    fcnsum->Write();
    f0->Close();

}


