#include <TH1F.h>
#include <TGraph.h>
#include <TTree.h>
#include <string>
#include <gsl/gsl_cdf.h>
#include <cmath>
#include <vector>
#include <stdlib.h>

using namespace std;

double Na = 6.022 * pow(10,23);
double T_2nubb_82Se = 9.39 * pow(10,19);
double T_2nubb_150Nd = 9.11 * pow(10,18);
double exposure = 17.5;
double A_208Tl = 2*31.5;
//double A_208Tl = 50*31.5;
double A_214Bi = 10*31.5;
//double A_214Bi = 200*31.5;
double mass_mol_82Se = 0.082; // kg/mol
double mass_mol_150Nd = 0.150;
double _nbr_events_0nubb_ = 0;
double _nbr_events_2nubb_ = 0;
double _nbr_events_2nubb_2MeV_ = 0;
double _nbr_events_214Bi_ = 0;
double _nbr_events_208Tl_ = 0;

//output file
std::string const final_rate("final_rate.txt");
std::ofstream final_flux(final_rate.c_str());


TH1F *henergy_sum(string isotope, bool field, string process, double xmin, double xmax, int nbins);
TH1F *hefficiency(string isotope, bool field, string process, double xmin, double xmax, int nbins, TH1F *histo_energy_sum = 0);
TH1F *hN_background(string isotope, bool field, string process, double xmin, double xmax, int nbins, TH1F *histo_energy_sum = 0,TH1F *histo_efficiency = 0);
double WindowMethodFindExpSigEvts(Double_t B);
Double_t sensitivity_FC(string isotope, int bin_emin, int bin_emax, TH1F *histo_energy_0nubb, TH1F *histo_tot);
TH2F *h2sensitivity(string isotope, bool field, double i_min, double i_max, double j_min, double j_max, double xmin, double xmax, int nbins,  TH1F *histo_energy_sum_0nubb = 0, TH1F *histo_energy_sum_2nubb_2MeV = 0, TH1F *histo_energy_sum_208Tl = 0, TH1F *histo_energy_sum_214Bi = 0,  TCanvas *c_sensitivity_spectrum = 0);
double ErrorStatEfficiency(string process, string isotope, double efficiency);
double ErrorStatbdf(string process, string isotope, double nbr_bdf, double efficiency);
search_ROI get_ROI(TH2F *histo_demie_vie = 0);
double mbb_min(double T12_max, string isotope);
double mbb_max(double T12_max, string isotope);


struct search_ROI{
  double Einf_ROI;
  double Esup_ROI;
  double T12_max;
};


void efficiency(string isotope, bool field, string activities){

  gStyle->SetPaintTextFormat("4.2e");

  TH1F *Nbackground_tot_spectrum = new TH1F("h_nbr_bdf_tot","Total expected number of background events",80,0,4);
  TH1F *energy_spectrum_0nubb = henergy_sum(isotope, field,"0nubb",0,4,80);
  TH1F *efficiency_spectrum_0nubb = hefficiency(isotope, field,"0nubb",0,4,80,energy_spectrum_0nubb);

  ///Drawing

  TCanvas *c_energy_spectrum = new TCanvas("canvas1","canvas1");
  TCanvas *c_efficiency_spectrum = new TCanvas("canvas2","canvas2");
  TCanvas *c_Nbackground_spectrum = new TCanvas("canvas3","canvas3");
  TCanvas *c_sensitivity_spectrum = new TCanvas("canvas4","canvas4");
  gStyle->SetOptStat(kFALSE);


  if (field && isotope == "82Se") {


    TH1F *energy_spectrum_2nubb = henergy_sum(isotope, field,"2nubb",0,4,80);
    TH1F *efficiency_spectrum_2nubb = hefficiency(isotope, field,"2nubb",0,4,80,energy_spectrum_2nubb);
    TH1F *Nbackground_spectrum_2nubb = hN_background(isotope, field,"2nubb",0,4,80,energy_spectrum_2nubb,efficiency_spectrum_2nubb);

    TH1F *energy_spectrum_2nubb_2MeV = henergy_sum(isotope, field,"2nubb_2MeV",0,4,80);
    TH1F *efficiency_spectrum_2nubb_2MeV = hefficiency(isotope, field,"2nubb_2MeV",0,4,80,energy_spectrum_2nubb_2MeV);
    TH1F *Nbackground_spectrum_2nubb_2MeV = hN_background(isotope, field,"2nubb_2MeV",0,4,80,energy_spectrum_2nubb_2MeV,efficiency_spectrum_2nubb_2MeV);

    TH1F *energy_spectrum_214Bi = henergy_sum(isotope, field,"214Bi",0,4,80);
    TH1F *efficiency_spectrum_214Bi = hefficiency(isotope, field,"214Bi",0,4,80,energy_spectrum_214Bi);
    TH1F *Nbackground_spectrum_214Bi = hN_background(isotope, field,"214Bi",0,4,80,energy_spectrum_214Bi,efficiency_spectrum_214Bi);

    TH1F *energy_spectrum_208Tl = henergy_sum(isotope, field,"208Tl",0,4,80);
    TH1F *efficiency_spectrum_208Tl = hefficiency(isotope, field,"208Tl",0,4,80,energy_spectrum_208Tl);
    TH1F *Nbackground_spectrum_208Tl = hN_background(isotope, field,"208Tl",0,4,80,energy_spectrum_208Tl,efficiency_spectrum_208Tl);


    c_energy_spectrum->cd();
    gPad->SetLogy();
    energy_spectrum_0nubb->Draw();
    energy_spectrum_0nubb->SetTitle("Total energy spectrum");
    energy_spectrum_0nubb->GetXaxis()->SetTitle("Etot");
    energy_spectrum_0nubb->GetYaxis()->SetTitle("#Event");
    energy_spectrum_0nubb->SetLineColor(1);

    energy_spectrum_2nubb->Draw("SAME");
    energy_spectrum_2nubb->SetLineColor(2);

    energy_spectrum_2nubb_2MeV->Draw("HIST SAME");
    energy_spectrum_2nubb_2MeV->SetLineColor(3);

    energy_spectrum_214Bi->Draw("SAME");
    energy_spectrum_214Bi->SetLineColor(4);

    energy_spectrum_208Tl->Draw("SAME");
    energy_spectrum_208Tl->SetLineColor(6);

    energy_spectrum_0nubb->GetYaxis()->SetRangeUser(1,pow(10,6));

    auto legend1 = new TLegend(0.75,0.6,0.89,0.89);
    legend1->AddEntry(energy_spectrum_0nubb,"0nubb","l");
    legend1->AddEntry(energy_spectrum_2nubb,"2nubb","l");
    legend1->AddEntry(energy_spectrum_2nubb_2MeV,"2nubb_2MeV","l");
    legend1->AddEntry(energy_spectrum_214Bi,"214Bi","l");
    legend1->AddEntry(energy_spectrum_208Tl,"208Tl","l");
    legend1->Draw();


    c_efficiency_spectrum->cd();
    gPad->SetLogy();
    efficiency_spectrum_0nubb->Draw();
    efficiency_spectrum_0nubb->SetTitle("Efficiency spectrum");
    efficiency_spectrum_0nubb->GetXaxis()->SetTitle("E>E_min");
    efficiency_spectrum_0nubb->GetYaxis()->SetTitle("Efficiency");
    efficiency_spectrum_0nubb->GetYaxis()->SetRangeUser(pow(10,-7),1);
    efficiency_spectrum_0nubb->SetLineColor(1);

    efficiency_spectrum_2nubb->Draw("SAME");
    efficiency_spectrum_2nubb->SetLineColor(2);

    efficiency_spectrum_2nubb_2MeV->Draw("SAME");
    efficiency_spectrum_2nubb_2MeV->SetLineColor(3);

    efficiency_spectrum_208Tl->Draw("SAME");
    efficiency_spectrum_208Tl->SetLineColor(6);


    efficiency_spectrum_214Bi->Draw("SAME");
    efficiency_spectrum_214Bi->SetLineColor(4);

    auto legend2 = new TLegend(0.75,0.6,0.89,0.89);
    legend2->AddEntry(efficiency_spectrum_0nubb,"0nubb","l");
    legend2->AddEntry(efficiency_spectrum_2nubb_2MeV,"2nubb_2MeV","l");
    legend2->AddEntry(efficiency_spectrum_2nubb,"2nubb","l");
    legend2->AddEntry(efficiency_spectrum_214Bi,"214Bi","l");
    legend2->AddEntry(efficiency_spectrum_208Tl,"208Tl","l");
    legend2->Draw();

    c_Nbackground_spectrum->cd();
    gPad->SetLogy();

    Nbackground_spectrum_2nubb->Draw();
    Nbackground_spectrum_2nubb->SetTitle("Expected number of background events");
    Nbackground_spectrum_2nubb->GetXaxis()->SetTitle("E>E_min");
    Nbackground_spectrum_2nubb->GetYaxis()->SetTitle("#Background");
    Nbackground_spectrum_2nubb->GetYaxis()->SetRangeUser(pow(10,-3),pow(10,6));
    Nbackground_spectrum_2nubb->SetLineColor(2);

    Nbackground_spectrum_2nubb_2MeV->Draw("SAME");
    Nbackground_spectrum_2nubb_2MeV->SetTitle("Nbr background events");
    Nbackground_spectrum_2nubb_2MeV->GetXaxis()->SetTitle("E>E_min");
    Nbackground_spectrum_2nubb_2MeV->GetYaxis()->SetTitle("Nbackground");
    Nbackground_spectrum_2nubb_2MeV->SetLineColor(3);

    Nbackground_spectrum_214Bi->Draw("SAME");
    Nbackground_spectrum_214Bi->SetLineColor(7);

    Nbackground_spectrum_208Tl->Draw("SAME");
    Nbackground_spectrum_208Tl->SetLineColor(6);

    Nbackground_tot_spectrum->Add(Nbackground_spectrum_2nubb,Nbackground_spectrum_2nubb_2MeV,1,0.0439913);
    Nbackground_tot_spectrum->Add(Nbackground_spectrum_208Tl);
    Nbackground_tot_spectrum->Add(Nbackground_spectrum_214Bi);
    Nbackground_tot_spectrum->Draw("SAME");

    auto legend3= new TLegend(0.75,0.6,0.89,0.89);

    legend3->AddEntry(Nbackground_spectrum_2nubb,"2nubb","l");
    legend3->AddEntry(Nbackground_spectrum_2nubb_2MeV,"2nubb_2MeV","l");
    legend3->AddEntry(Nbackground_spectrum_214Bi,"214Bi","l");
    legend3->AddEntry(Nbackground_spectrum_208Tl,"208Tl","l");
    legend3->AddEntry(Nbackground_tot_spectrum,"TOTAL","l");
    legend3->Draw();


    TH2F *sensitivity_spectrum = new TH2F("h2_ROI","low bound vs up bound",10,2.45,2.95,20,2.45,3.45);

    Double_t Nbackground = 0;
    Double_t Nbackground_2nubb = 0;
    Double_t Nbackground_2nubb_2MeV = 0;
    Double_t Nbackground_208Tl = 0;
    Double_t Nbackground_214Bi = 0;
    Double_t sensitivity = 0;
    Double_t efficiency_0nubb = 0;
    Double_t efficiency_2nubb = 0;
    Double_t efficiency_2nubb_2MeV = 0;
    Double_t efficiency_208Tl = 0;
    Double_t efficiency_214Bi = 0;
    Double_t expectedSignalEventLimit = 0;
    Double_t Nbackground_222Rn = 0;

    for (int i=Nbackground_tot_spectrum->GetXaxis()->FindBin(2.45);i<=Nbackground_tot_spectrum->GetXaxis()->FindBin(2.95);i++) {
      for (int j=i+1;j<=Nbackground_tot_spectrum->GetXaxis()->FindBin(3.45);j++) {

        efficiency_0nubb += energy_spectrum_0nubb->Integral(i,j);
        efficiency_0nubb *= 1./pow(10,7);

        efficiency_2nubb_2MeV += energy_spectrum_2nubb_2MeV->Integral(i,j);
        efficiency_2nubb_2MeV *= 1./pow(10,7);
        Nbackground_2nubb_2MeV = ((Na*log(2))/mass_mol_82Se)*((efficiency_2nubb_2MeV*exposure)/T_2nubb_82Se);

        efficiency_214Bi += energy_spectrum_214Bi->Integral(i,j);
        efficiency_214Bi *= 1./pow(10,7);
        Nbackground_214Bi = A_214Bi*efficiency_214Bi*exposure;

        efficiency_208Tl += energy_spectrum_208Tl->Integral(i,j);
        efficiency_208Tl *= 1./pow(10,7);
        Nbackground_208Tl = A_208Tl*efficiency_208Tl*exposure;

        Nbackground_222Rn=Nbackground_214Bi*4.39;

        Nbackground = Nbackground_2nubb_2MeV+Nbackground_208Tl+Nbackground_214Bi+Nbackground_222Rn;
        if (Nbackground > 200.) {
          expectedSignalEventLimit = TMath::Sqrt(Nbackground);
        }
        else {
          expectedSignalEventLimit = WindowMethodFindExpSigEvts(Nbackground);
        }

        sensitivity = ((Na*log(2))/mass_mol_82Se)*((efficiency_0nubb*exposure)/expectedSignalEventLimit);

        sensitivity_spectrum->SetBinContent(i+1-Nbackground_tot_spectrum->GetXaxis()->FindBin(2.45),j+2-Nbackground_tot_spectrum->GetXaxis()->FindBin(2.5),sensitivity);

      }
    }

    cout << "T12 = " << get_ROI(sensitivity_spectrum).T12_max << endl;
    cout << "mbb = " << "[" << mbb_min(get_ROI(sensitivity_spectrum).T12_max,isotope) << "," << mbb_max(get_ROI(sensitivity_spectrum).T12_max,isotope) << "]" << endl;

    double efficiency_0nubb_ROI = energy_spectrum_0nubb->Integral(Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Einf_ROI),Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Esup_ROI));
    cout << "0nubb " << "[" << get_ROI(sensitivity_spectrum).Einf_ROI << "," << get_ROI(sensitivity_spectrum).Esup_ROI << "] = " << efficiency_0nubb_ROI/pow(10,7) << endl;


    double efficiency_2nubb_2MeV_ROI = energy_spectrum_2nubb_2MeV->Integral(Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Einf_ROI),Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Esup_ROI));
    double Nbackground_2nubb_2MeV_ROI = ((Na*log(2))/mass_mol_82Se)*((efficiency_2nubb_2MeV_ROI*exposure)/T_2nubb_82Se);
    cout << "2nubb_2MeV " << "[" << get_ROI(sensitivity_spectrum).Einf_ROI << "," << get_ROI(sensitivity_spectrum).Esup_ROI << "] = " << Nbackground_2nubb_2MeV_ROI/pow(10,7) << endl;


    double efficiency_208Tl_ROI = energy_spectrum_208Tl->Integral(Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Einf_ROI),Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Esup_ROI));
    double Nbackground_208Tl_ROI = A_208Tl*efficiency_208Tl_ROI*exposure;
    cout << "208Tl " << "[" << get_ROI(sensitivity_spectrum).Einf_ROI << "," << get_ROI(sensitivity_spectrum).Esup_ROI << "] = " << Nbackground_208Tl_ROI/pow(10,7) << endl;


    double efficiency_214Bi_ROI = energy_spectrum_214Bi->Integral(Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Einf_ROI),Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Esup_ROI));
    double Nbackground_214Bi_ROI = A_214Bi*efficiency_214Bi_ROI*exposure;
    cout << "214Bi " << "[" << get_ROI(sensitivity_spectrum).Einf_ROI << "," << get_ROI(sensitivity_spectrum).Esup_ROI << "] = " << Nbackground_214Bi_ROI/pow(10,7) << endl;



    c_sensitivity_spectrum->cd();
    TGaxis::SetMaxDigits(2);
    sensitivity_spectrum->Draw("colz TEXT");
    sensitivity_spectrum->SetTitle("Expected sensitivity");
    sensitivity_spectrum->GetXaxis()->SetTitle("Inf ROI (MeV)");
    sensitivity_spectrum->GetYaxis()->SetTitle("Sup ROI (MeV)");


    // TH2F *histo2_sensitivity = h2sensitivity(isotope, field, 2.45, 2.95, 2.5, 3.45, 0, 4, 80, energy_spectrum_0nubb, energy_spectrum_2nubb_2MeV, energy_spectrum_208Tl, energy_spectrum_214Bi, c_sensitivity_spectrum);

    energy_spectrum_0nubb->SaveAs("root_outpute/energy_spectrum_with_B_82Se_0nubb.root");
    energy_spectrum_2nubb->SaveAs("root_outpute/energy_spectrum_with_B_82Se_2nubb.root");
    energy_spectrum_2nubb_2MeV->SaveAs("root_outpute/energy_spectrum_with_B_82Se_2nubb_2MeV.root");
    energy_spectrum_208Tl->SaveAs("root_outpute/energy_spectrum_with_B_82Se_208Tl.root");
    energy_spectrum_214Bi->SaveAs("root_outpute/energy_spectrum_with_B_82Se_214Bi.root");

    efficiency_spectrum_0nubb->SaveAs("root_outpute/efficiency_spectrum_with_B_82Se_0nubb.root");
    efficiency_spectrum_2nubb->SaveAs("root_outpute/efficiency_spectrum_with_B_82Se_2nubb.root");
    efficiency_spectrum_2nubb_2MeV->SaveAs("root_outpute/efficiency_spectrum_with_B_82Se_2nubb_2MeV.root");
    efficiency_spectrum_208Tl->SaveAs("root_outpute/efficiency_spectrum_with_B_82Se_208Tl.root");
    efficiency_spectrum_214Bi->SaveAs("root_outpute/efficiency_spectrum_with_B_82Se_214Bi.root");

    Nbackground_spectrum_2nubb->SaveAs("root_outpute/Nbackground_spectrum_with_B_82Se_2nubb.root");
    Nbackground_spectrum_2nubb_2MeV->SaveAs("root_outpute/Nbackground_spectrum_with_B_82Se_2nubb_2MeV.root");
    Nbackground_spectrum_208Tl->SaveAs("root_outpute/Nbackground_spectrum_with_B_82Se_208Tl.root");
    Nbackground_spectrum_214Bi->SaveAs("root_outpute/Nbackground_spectrum_with_B_82Se_214Bi.root");

  }

  else if (field && isotope == "150Nd") {

    TH1F *energy_spectrum_2nubb = henergy_sum(isotope, field,"2nubb",0,4,80);
    TH1F *efficiency_spectrum_2nubb = hefficiency(isotope, field,"2nubb",0,4,80,energy_spectrum_2nubb);
    TH1F *Nbackground_spectrum_2nubb = hN_background(isotope, field,"2nubb",0,4,80,energy_spectrum_2nubb,efficiency_spectrum_2nubb);

    TH1F *energy_spectrum_214Bi = henergy_sum(isotope, field,"214Bi",0,4,80);
    TH1F *efficiency_spectrum_214Bi = hefficiency(isotope, field,"214Bi",0,4,80,energy_spectrum_214Bi);
    TH1F *Nbackground_spectrum_214Bi = hN_background(isotope, field,"214Bi",0,4,80,energy_spectrum_214Bi,efficiency_spectrum_214Bi);

    TH1F *energy_spectrum_208Tl = henergy_sum(isotope, field,"208Tl",0,4,80);
    TH1F *efficiency_spectrum_208Tl = hefficiency(isotope, field,"208Tl",0,4,80,energy_spectrum_208Tl);
    TH1F *Nbackground_spectrum_208Tl = hN_background(isotope, field,"208Tl",0,4,80,energy_spectrum_208Tl,efficiency_spectrum_208Tl);

    // TH2F *sensitivity = hSensitivity(isotope, field,"2nubb_2MeV",0,4,80,energy_spectrum_2nubb_2MeV,efficiency_spectrum_2nubb_2MeV,Nbackground_spectrum_2nubb_2MeV);

    c_energy_spectrum->cd();
    gPad->SetLogy();
    energy_spectrum_0nubb->Draw();
    energy_spectrum_0nubb->SetTitle("Total energy spectrum");
    energy_spectrum_0nubb->GetXaxis()->SetTitle("Etot");
    energy_spectrum_0nubb->GetYaxis()->SetTitle("#Event");
    energy_spectrum_0nubb->GetYaxis()->SetRangeUser(1,pow(10,6));
    energy_spectrum_0nubb->SetLineColor(1);

    energy_spectrum_2nubb->Draw("SAME");
    energy_spectrum_2nubb->SetLineColor(2);

    energy_spectrum_214Bi->Draw("SAME");
    energy_spectrum_214Bi->SetLineColor(4);

    energy_spectrum_208Tl->SetLineColor(6);
    energy_spectrum_208Tl->Draw("SAME");

    auto legend1 = new TLegend(0.75,0.6,0.89,0.89);
    legend1->AddEntry(energy_spectrum_0nubb,"0nubb","l");
    legend1->AddEntry(energy_spectrum_2nubb,"2nubb","l");
    legend1->AddEntry(energy_spectrum_214Bi,"214Bi","l");
    legend1->AddEntry(energy_spectrum_208Tl,"208Tl","l");
    legend1->Draw();


    c_efficiency_spectrum->cd();
    gPad->SetLogy();
    efficiency_spectrum_0nubb->Draw();
    efficiency_spectrum_0nubb->SetTitle("Efficiency spectrum");
    efficiency_spectrum_0nubb->GetXaxis()->SetTitle("E>E_min");
    efficiency_spectrum_0nubb->GetYaxis()->SetTitle("Efficiency");
    efficiency_spectrum_0nubb->GetYaxis()->SetRangeUser(pow(10,-7),1);
    efficiency_spectrum_0nubb->SetLineColor(1);

    efficiency_spectrum_2nubb->Draw("SAME");

    efficiency_spectrum_214Bi->Draw("SAME");
    efficiency_spectrum_214Bi->SetLineColor(4);

    efficiency_spectrum_208Tl->Draw("SAME");
    efficiency_spectrum_208Tl->SetLineColor(6);

    efficiency_spectrum_2nubb->SetLineColor(2);


    auto legend2 = new TLegend(0.75,0.6,0.89,0.89);
    legend2->AddEntry(efficiency_spectrum_0nubb,"0nubb","l");
    legend2->AddEntry(efficiency_spectrum_2nubb,"2nubb","l");
    legend2->AddEntry(efficiency_spectrum_214Bi,"214Bi","l");
    legend2->AddEntry(efficiency_spectrum_208Tl,"208Tl","l");
    legend2->Draw();

    c_Nbackground_spectrum->cd();
    gPad->SetLogy();

    Nbackground_spectrum_2nubb->Draw();
    Nbackground_spectrum_2nubb->GetYaxis()->SetRangeUser(pow(10,-3),pow(10,6));
    Nbackground_spectrum_2nubb->SetTitle("Expected number of background events");
    Nbackground_spectrum_2nubb->GetXaxis()->SetTitle("E>E_min");
    Nbackground_spectrum_2nubb->GetYaxis()->SetTitle("#Background");
    Nbackground_spectrum_2nubb->SetLineColor(2);

    Nbackground_spectrum_214Bi->Draw("SAME");
    Nbackground_spectrum_214Bi->SetLineColor(7);

    Nbackground_spectrum_208Tl->Draw("SAME");
    Nbackground_spectrum_208Tl->SetLineColor(6);

    Nbackground_tot_spectrum->Add(Nbackground_spectrum_2nubb);
    Nbackground_tot_spectrum->Add(Nbackground_spectrum_208Tl);
    Nbackground_tot_spectrum->Add(Nbackground_spectrum_214Bi);
    Nbackground_tot_spectrum->Draw("SAME");


    auto legend3= new TLegend(0.75,0.6,0.89,0.89);
    legend3->AddEntry(Nbackground_spectrum_2nubb,"2nubb","l");
    legend3->AddEntry(Nbackground_spectrum_214Bi,"214Bi","l");
    legend3->AddEntry(Nbackground_spectrum_208Tl,"208Tl","l");
    legend3->AddEntry(Nbackground_tot_spectrum,"TOTAL","l");
    legend3->Draw();

    efficiency_spectrum_0nubb->SaveAs("root_outpute/efficiency_spectrum_with_B_150Nd_0nubb.root");
    efficiency_spectrum_2nubb->SaveAs("root_outpute/efficiency_spectrum_with_B_150Nd_2nubb.root");
    efficiency_spectrum_208Tl->SaveAs("root_outpute/efficiency_spectrum_with_B_150Nd_208Tl.root");
    efficiency_spectrum_214Bi->SaveAs("root_outpute/efficiency_spectrum_with_B_150Nd_214Bi.root");


    TH2F *sensitivity_spectrum = new TH2F("h2_ROI","low bound vs up bound",10,2.85,3.35,20,2.85,3.85);

    Double_t Nbackground = 0;
    Double_t Nbackground_2nubb = 0;
    Double_t Nbackground_208Tl = 0;
    Double_t Nbackground_214Bi = 0;
    Double_t Nbackground_222Rn = 0;
    Double_t sensitivity = 0;
    Double_t efficiency_0nubb = 0;
    Double_t efficiency_2nubb = 0;
    Double_t efficiency_2nubb_2MeV = 0;
    Double_t efficiency_208Tl = 0;
    Double_t efficiency_214Bi = 0;
    Double_t expectedSignalEventLimit = 0;
    int cointer = 0;

    for (int i=Nbackground_tot_spectrum->GetXaxis()->FindBin(2.85);i<=Nbackground_tot_spectrum->GetXaxis()->FindBin(3.35);i++) {
      for (int j=i+1;j<=Nbackground_tot_spectrum->GetXaxis()->FindBin(3.85);j++) {

        efficiency_0nubb += energy_spectrum_0nubb->Integral(i,j);
        efficiency_0nubb *= 1./9574121;

        efficiency_2nubb += energy_spectrum_2nubb->Integral(i,j);
        efficiency_2nubb *= 1./9656794;
        Nbackground_2nubb = ((Na*log(2))/mass_mol_150Nd)*((efficiency_2nubb*exposure)/T_2nubb_150Nd);

        efficiency_214Bi += energy_spectrum_214Bi->Integral(i,j);
        efficiency_214Bi *= 1./9778522;
        Nbackground_214Bi = A_214Bi*efficiency_214Bi*exposure;

        efficiency_208Tl += energy_spectrum_208Tl->Integral(i,j);
        efficiency_208Tl *= 1./9773130;
        Nbackground_208Tl = A_208Tl*efficiency_208Tl*exposure;

        Nbackground_222Rn = Nbackground_214Bi*4.39;

        Nbackground = Nbackground_2nubb+Nbackground_208Tl+Nbackground_214Bi+Nbackground_222Rn;
        if (Nbackground > 200.) {
          expectedSignalEventLimit = TMath::Sqrt(Nbackground);
        }
        else {
          expectedSignalEventLimit = WindowMethodFindExpSigEvts(Nbackground);
        }


        sensitivity = ((Na*log(2))/mass_mol_150Nd)*((efficiency_0nubb*exposure)/expectedSignalEventLimit);

        sensitivity_spectrum->SetBinContent(i+1-Nbackground_tot_spectrum->GetXaxis()->FindBin(2.85),j+2-Nbackground_tot_spectrum->GetXaxis()->FindBin(2.9),sensitivity);

      }
    }

    cout << "T12 = " << get_ROI(sensitivity_spectrum).T12_max << endl;
    cout << "mbb = " << "[" << mbb_min(get_ROI(sensitivity_spectrum).T12_max,isotope) << "," << mbb_max(get_ROI(sensitivity_spectrum).T12_max,isotope) << "]" << endl;

    double efficiency_0nubb_ROI = energy_spectrum_0nubb->Integral(Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Einf_ROI),Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Esup_ROI));
    cout << "0nubb " << "[" << get_ROI(sensitivity_spectrum).Einf_ROI << "," << get_ROI(sensitivity_spectrum).Esup_ROI << "] = " << efficiency_0nubb_ROI/9574121 << endl;


    double efficiency_2nubb_2MeV_ROI = energy_spectrum_2nubb->Integral(Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Einf_ROI),Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Esup_ROI));
    double Nbackground_2nubb_2MeV_ROI = ((Na*log(2))/mass_mol_82Se)*((efficiency_2nubb_2MeV_ROI*exposure)/T_2nubb_82Se);
    cout << "2nubb_2MeV " << "[" << get_ROI(sensitivity_spectrum).Einf_ROI << "," << get_ROI(sensitivity_spectrum).Esup_ROI << "] = " << Nbackground_2nubb_2MeV_ROI/9656794 << endl;


    double efficiency_208Tl_ROI = energy_spectrum_208Tl->Integral(Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Einf_ROI),Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Esup_ROI));
    double Nbackground_208Tl_ROI = A_208Tl*efficiency_208Tl_ROI*exposure;
    cout << "208Tl " << "[" << get_ROI(sensitivity_spectrum).Einf_ROI << "," << get_ROI(sensitivity_spectrum).Esup_ROI << "] = " << Nbackground_208Tl_ROI/9773130 << endl;

    double efficiency_214Bi_ROI = energy_spectrum_214Bi->Integral(Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Einf_ROI),Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Esup_ROI));
    double Nbackground_214Bi_ROI = A_214Bi*efficiency_214Bi_ROI*exposure;
    cout << "214Bi " << "[" << get_ROI(sensitivity_spectrum).Einf_ROI << "," << get_ROI(sensitivity_spectrum).Esup_ROI << "] = " << Nbackground_214Bi_ROI/9778522 << endl;

    c_sensitivity_spectrum->cd();
    TGaxis::SetMaxDigits(2);
    sensitivity_spectrum->Draw("colz TEXT");
    sensitivity_spectrum->SetTitle("Expected sensitivity");
    sensitivity_spectrum->GetXaxis()->SetTitle("Sup ROI (MeV)");
    sensitivity_spectrum->GetYaxis()->SetTitle("Inf ROI (MeV)");

  }

  else if (!field && isotope == "82Se") {

    TH1F *energy_spectrum_2nubb_2MeV = henergy_sum(isotope, field,"2nubb_2MeV",0,4,80);
    TH1F *efficiency_spectrum_2nubb_2MeV = hefficiency(isotope, field,"2nubb_2MeV",0,4,80,energy_spectrum_2nubb_2MeV);
    TH1F *Nbackground_spectrum_2nubb_2MeV = hN_background(isotope, field,"2nubb_2MeV",0,4,80,energy_spectrum_2nubb_2MeV,efficiency_spectrum_2nubb_2MeV);

    TH1F *energy_spectrum_214Bi = henergy_sum(isotope, field,"214Bi",0,4,80);
    TH1F *efficiency_spectrum_214Bi = hefficiency(isotope, field,"214Bi",0,4,80,energy_spectrum_214Bi);
    TH1F *Nbackground_spectrum_214Bi = hN_background(isotope, field,"214Bi",0,4,80,energy_spectrum_214Bi,efficiency_spectrum_214Bi);

    TH1F *energy_spectrum_208Tl = henergy_sum(isotope, field,"208Tl",0,4,80);
    TH1F *efficiency_spectrum_208Tl = hefficiency(isotope, field,"208Tl",0,4,80,energy_spectrum_208Tl);
    TH1F *Nbackground_spectrum_208Tl = hN_background(isotope, field,"208Tl",0,4,80,energy_spectrum_208Tl,efficiency_spectrum_208Tl);

    // TH2F *sensitivity = hSensitivity(isotope, field,"2nubb_2MeV",0,4,80,energy_spectrum_2nubb_2MeV,efficiency_spectrum_2nubb_2MeV,Nbackground_spectrum_2nubb_2MeV);

    c_energy_spectrum->cd();
    gPad->SetLogy();
    energy_spectrum_0nubb->Draw("HIST");
    energy_spectrum_0nubb->SetTitle("Total energy spectrum");
    energy_spectrum_0nubb->GetXaxis()->SetTitle("Etot");
    energy_spectrum_0nubb->GetYaxis()->SetTitle("#Event");
    energy_spectrum_0nubb->GetYaxis()->SetRangeUser(1,pow(10,6));
    energy_spectrum_0nubb->SetLineColor(1);

    energy_spectrum_2nubb_2MeV->Draw("HIST SAME");
    energy_spectrum_2nubb_2MeV->SetLineColor(3);

    energy_spectrum_214Bi->Draw("HIST SAME");
    energy_spectrum_214Bi->SetLineColor(4);

    energy_spectrum_208Tl->Draw("HIST SAME");
    energy_spectrum_208Tl->SetLineColor(6);


    auto legend1 = new TLegend(0.75,0.6,0.89,0.89);
    legend1->AddEntry(energy_spectrum_0nubb,"0nubb","l");
    legend1->AddEntry(energy_spectrum_2nubb_2MeV,"2nubb_2MeV","l");
    legend1->AddEntry(energy_spectrum_214Bi,"214Bi","l");
    legend1->AddEntry(energy_spectrum_208Tl,"208Tl","l");
    legend1->Draw();


    c_efficiency_spectrum->cd();
    gPad->SetLogy();
    efficiency_spectrum_0nubb->Draw();
    efficiency_spectrum_0nubb->SetTitle("Efficiency spectrum");
    efficiency_spectrum_0nubb->GetXaxis()->SetTitle("E>E_min");
    efficiency_spectrum_0nubb->GetYaxis()->SetTitle("Efficiency");
    efficiency_spectrum_0nubb->GetYaxis()->SetRangeUser(pow(10,-7),1);
    efficiency_spectrum_0nubb->SetLineColor(1);

    efficiency_spectrum_2nubb_2MeV->Draw("SAME");
    efficiency_spectrum_2nubb_2MeV->SetLineColor(3);

    efficiency_spectrum_214Bi->Draw("SAME");
    efficiency_spectrum_214Bi->SetLineColor(4);

    efficiency_spectrum_208Tl->Draw("SAME");
    efficiency_spectrum_208Tl->SetLineColor(6);

    auto legend2 = new TLegend(0.75,0.6,0.89,0.89);
    legend2->AddEntry(efficiency_spectrum_0nubb,"0nubb","l");
    legend2->AddEntry(efficiency_spectrum_2nubb_2MeV,"2nubb_2MeV","l");
    legend2->AddEntry(efficiency_spectrum_214Bi,"214Bi","l");
    legend2->AddEntry(efficiency_spectrum_208Tl,"208Tl","l");
    legend2->Draw();

    c_Nbackground_spectrum->cd();
    gPad->SetLogy();

    Nbackground_spectrum_2nubb_2MeV->Draw();
    Nbackground_spectrum_2nubb_2MeV->SetTitle("Expected number of background events");
    Nbackground_spectrum_2nubb_2MeV->GetXaxis()->SetTitle("E>E_min");
    Nbackground_spectrum_2nubb_2MeV->GetYaxis()->SetTitle("#Background");
    Nbackground_spectrum_2nubb_2MeV->GetYaxis()->SetRangeUser(pow(10,-3),pow(10,6));
    Nbackground_spectrum_2nubb_2MeV->SetLineColor(3);

    Nbackground_spectrum_214Bi->Draw("SAME");
    Nbackground_spectrum_214Bi->SetLineColor(7);

    Nbackground_spectrum_208Tl->Draw("SAME");
    Nbackground_spectrum_208Tl->SetLineColor(6);

    Nbackground_tot_spectrum->Add(Nbackground_spectrum_2nubb_2MeV);
    Nbackground_tot_spectrum->Add(Nbackground_spectrum_208Tl);
    Nbackground_tot_spectrum->Add(Nbackground_spectrum_214Bi);
    Nbackground_tot_spectrum->Draw("SAME");


    auto legend3= new TLegend(0.75,0.6,0.89,0.89);
    legend3->AddEntry(Nbackground_spectrum_2nubb_2MeV,"2nubb_2MeV","l");
    legend3->AddEntry(Nbackground_spectrum_214Bi,"214Bi","l");
    legend3->AddEntry(Nbackground_spectrum_208Tl,"208Tl","l");
    legend3->AddEntry(Nbackground_tot_spectrum,"TOTAL","l");
    legend3->Draw();



    TH2F *sensitivity_spectrum = new TH2F("h2_ROI","low bound vs up bound",10,2.45,2.95,20,2.45,3.45);

    Double_t Nbackground = 0;
    Double_t Nbackground_2nubb = 0;
    Double_t Nbackground_2nubb_2MeV = 0;
    Double_t Nbackground_208Tl = 0;
    Double_t Nbackground_214Bi = 0;
    Double_t Nbackground_222Rn;
    if (activities == "zero") {
      Nbackground_222Rn = 0;
    }
    else {
      Nbackground_222Rn = 0.39;
    }
    Double_t sensitivity = 0;
    Double_t efficiency_0nubb = 0;
    Double_t efficiency_2nubb = 0;
    Double_t efficiency_2nubb_2MeV = 0;
    Double_t efficiency_208Tl = 0;
    Double_t efficiency_214Bi = 0;
    Double_t expectedSignalEventLimit = 0;

    for (int i=Nbackground_tot_spectrum->GetXaxis()->FindBin(2.45);i<=Nbackground_tot_spectrum->GetXaxis()->FindBin(2.95);i++) {
      for (int j=i+1;j<=Nbackground_tot_spectrum->GetXaxis()->FindBin(3.45);j++) {

        efficiency_0nubb += energy_spectrum_0nubb->Integral(i,j);
        efficiency_0nubb *= 1./pow(10,7);

        efficiency_2nubb_2MeV += energy_spectrum_2nubb_2MeV->Integral(i,j);
        efficiency_2nubb_2MeV *= 1./pow(10,7);
        Nbackground_2nubb_2MeV = ((Na*log(2))/mass_mol_82Se)*((efficiency_2nubb_2MeV*exposure)/T_2nubb_82Se);

        efficiency_214Bi += energy_spectrum_214Bi->Integral(i,j);
        efficiency_214Bi *= 1./pow(10,7);
        Nbackground_214Bi = A_214Bi*efficiency_214Bi*exposure;
        Nbackground_222Rn = Nbackground_214Bi*4.39;

        efficiency_208Tl += energy_spectrum_208Tl->Integral(i,j);
        efficiency_208Tl *= 1./pow(10,7);
        Nbackground_208Tl = A_208Tl*efficiency_208Tl*exposure;

        Nbackground = Nbackground_2nubb_2MeV+Nbackground_208Tl+Nbackground_214Bi;
        if (Nbackground > 200.) {
          expectedSignalEventLimit = TMath::Sqrt(Nbackground);
        }
        else {
          expectedSignalEventLimit = WindowMethodFindExpSigEvts(Nbackground);
        }


        sensitivity = ((Na*log(2))/mass_mol_82Se)*((efficiency_0nubb*exposure)/expectedSignalEventLimit);

        sensitivity_spectrum->SetBinContent(i+1-Nbackground_tot_spectrum->GetXaxis()->FindBin(2.45),j+2-Nbackground_tot_spectrum->GetXaxis()->FindBin(2.5),sensitivity);

      }
    }

    cout << "T12 = " << get_ROI(sensitivity_spectrum).T12_max << endl;
    cout << "mbb = " << "[" << mbb_min(get_ROI(sensitivity_spectrum).T12_max,isotope) << "," << mbb_max(get_ROI(sensitivity_spectrum).T12_max,isotope) << "]" << endl;

    double efficiency_0nubb_ROI = energy_spectrum_0nubb->Integral(Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Einf_ROI),Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Esup_ROI));
    cout << "0nubb " << "[" << get_ROI(sensitivity_spectrum).Einf_ROI << "," << get_ROI(sensitivity_spectrum).Esup_ROI << "] = " << efficiency_0nubb_ROI/pow(10,7) << endl;


    double efficiency_2nubb_2MeV_ROI = energy_spectrum_2nubb_2MeV->Integral(Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Einf_ROI),Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Esup_ROI));
    double Nbackground_2nubb_2MeV_ROI = ((Na*log(2))/mass_mol_82Se)*((efficiency_2nubb_2MeV_ROI*exposure)/T_2nubb_82Se);
    cout << "2nubb_2MeV " << "[" << get_ROI(sensitivity_spectrum).Einf_ROI << "," << get_ROI(sensitivity_spectrum).Esup_ROI << "] = " << Nbackground_2nubb_2MeV_ROI/pow(10,7) << endl;


    double efficiency_208Tl_ROI = energy_spectrum_208Tl->Integral(Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Einf_ROI),Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Esup_ROI));
    double Nbackground_208Tl_ROI = A_208Tl*efficiency_208Tl_ROI*exposure;
    cout << "208Tl " << "[" << get_ROI(sensitivity_spectrum).Einf_ROI << "," << get_ROI(sensitivity_spectrum).Esup_ROI << "] = " << Nbackground_208Tl_ROI/pow(10,7) << endl;


    double efficiency_214Bi_ROI = energy_spectrum_214Bi->Integral(Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Einf_ROI),Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Esup_ROI));
    double Nbackground_214Bi_ROI = A_214Bi*efficiency_214Bi_ROI*exposure;
    cout << "214Bi " << "[" << get_ROI(sensitivity_spectrum).Einf_ROI << "," << get_ROI(sensitivity_spectrum).Esup_ROI << "] = " << Nbackground_214Bi_ROI/pow(10,7) << endl;


    c_sensitivity_spectrum->cd();
    TGaxis::SetMaxDigits(2);
    sensitivity_spectrum->Draw("colz TEXT");
    sensitivity_spectrum->SetTitle("Expected sensitivity");
    sensitivity_spectrum->GetXaxis()->SetTitle("Sup ROI (MeV)");
    sensitivity_spectrum->GetYaxis()->SetTitle("Inf ROI (MeV)");


    energy_spectrum_0nubb->SaveAs("root_outpute/energy_spectrum_without_B_82Se_0nubb.root");
    energy_spectrum_2nubb_2MeV->SaveAs("root_outpute/energy_spectrum_without_B_82Se_2nubb_2MeV.root");
    energy_spectrum_208Tl->SaveAs("root_outpute/energy_spectrum_without_B_82Se_208Tl.root");
    energy_spectrum_214Bi->SaveAs("root_outpute/energy_spectrum_without_B_82Se_214Bi.root");

    efficiency_spectrum_0nubb->SaveAs("root_outpute/efficiency_spectrum_without_B_82Se_0nubb.root");
    efficiency_spectrum_2nubb_2MeV->SaveAs("root_outpute/efficiency_spectrum_without_B_82Se_2nubb_2MeV.root");
    efficiency_spectrum_208Tl->SaveAs("root_outpute/efficiency_spectrum_without_B_82Se_208Tl.root");
    efficiency_spectrum_214Bi->SaveAs("root_outpute/efficiency_spectrum_without_B_82Se_214Bi.root");

    Nbackground_spectrum_2nubb_2MeV->SaveAs("root_outpute/Nbackground_spectrum_without_B_82Se_2nubb_2MeV.root");
    Nbackground_spectrum_208Tl->SaveAs("root_outpute/Nbackground_spectrum_without_B_82Se_208Tl.root");
    Nbackground_spectrum_214Bi->SaveAs("root_outpute/Nbackground_spectrum_without_B_82Se_214Bi.root");

  }

  else if (!field && isotope == "150Nd") {

    TH1F *energy_spectrum_2nubb = henergy_sum(isotope, field,"2nubb",0,4,80);
    TH1F *efficiency_spectrum_2nubb = hefficiency(isotope, field,"2nubb",0,4,80,energy_spectrum_2nubb);
    TH1F *Nbackground_spectrum_2nubb = hN_background(isotope, field,"2nubb",0,4,80,energy_spectrum_2nubb,efficiency_spectrum_2nubb);

    TH1F *energy_spectrum_214Bi = henergy_sum(isotope, field,"214Bi",0,4,80);
    TH1F *efficiency_spectrum_214Bi = hefficiency(isotope, field,"214Bi",0,4,80,energy_spectrum_214Bi);
    TH1F *Nbackground_spectrum_214Bi = hN_background(isotope, field,"214Bi",0,4,80,energy_spectrum_214Bi,efficiency_spectrum_214Bi);

    TH1F *energy_spectrum_208Tl = henergy_sum(isotope, field,"208Tl",0,4,80);
    TH1F *efficiency_spectrum_208Tl = hefficiency(isotope, field,"208Tl",0,4,80,energy_spectrum_208Tl);
    TH1F *Nbackground_spectrum_208Tl = hN_background(isotope, field,"208Tl",0,4,80,energy_spectrum_208Tl,efficiency_spectrum_208Tl);

    // TH2F *sensitivity = hSensitivity(isotope, field,"2nubb_2MeV",0,4,80,energy_spectrum_2nubb_2MeV,efficiency_spectrum_2nubb_2MeV,Nbackground_spectrum_2nubb_2MeV);

    c_energy_spectrum->cd();
    gPad->SetLogy();
    energy_spectrum_0nubb->Draw("HIST");
    energy_spectrum_0nubb->SetTitle("Total energy spectrum");
    energy_spectrum_0nubb->GetXaxis()->SetTitle("Etot");
    energy_spectrum_0nubb->GetYaxis()->SetTitle("#Event");
    energy_spectrum_0nubb->GetYaxis()->SetRangeUser(1,pow(10,6));
    energy_spectrum_0nubb->SetLineColor(1);

    energy_spectrum_2nubb->Draw("HIST SAME");
    energy_spectrum_2nubb->SetLineColor(2);

    energy_spectrum_214Bi->Draw("HIST SAME");
    energy_spectrum_214Bi->SetLineColor(4);

    energy_spectrum_208Tl->SetLineColor(6);
    energy_spectrum_208Tl->Draw("HIST SAME");

    auto legend1 = new TLegend(0.75,0.6,0.89,0.89);
    legend1->AddEntry(energy_spectrum_0nubb,"0nubb","l");
    legend1->AddEntry(energy_spectrum_2nubb,"2nubb","l");
    legend1->AddEntry(energy_spectrum_214Bi,"214Bi","l");
    legend1->AddEntry(energy_spectrum_208Tl,"208Tl","l");
    legend1->Draw();


    c_efficiency_spectrum->cd();
    gPad->SetLogy();
    efficiency_spectrum_0nubb->Draw();
    efficiency_spectrum_0nubb->SetTitle("Efficiency spectrum");
    efficiency_spectrum_0nubb->GetXaxis()->SetTitle("E>E_min");
    efficiency_spectrum_0nubb->GetYaxis()->SetTitle("Efficiency");
    efficiency_spectrum_0nubb->GetYaxis()->SetRangeUser(pow(10,-7),1);
    efficiency_spectrum_0nubb->SetLineColor(1);

    efficiency_spectrum_2nubb->Draw("SAME");
    efficiency_spectrum_2nubb->SetLineColor(2);

    efficiency_spectrum_214Bi->Draw("SAME");
    efficiency_spectrum_214Bi->SetLineColor(4);

    efficiency_spectrum_208Tl->Draw("SAME");
    efficiency_spectrum_208Tl->SetLineColor(6);



    auto legend2 = new TLegend(0.75,0.6,0.89,0.89);
    legend2->AddEntry(efficiency_spectrum_0nubb,"0nubb","l");
    legend2->AddEntry(efficiency_spectrum_2nubb,"2nubb","l");
    legend2->AddEntry(efficiency_spectrum_214Bi,"214Bi","l");
    legend2->AddEntry(efficiency_spectrum_208Tl,"208Tl","l");
    legend2->Draw();

    c_Nbackground_spectrum->cd();
    gPad->SetLogy();

    Nbackground_spectrum_2nubb->Draw();
    Nbackground_spectrum_2nubb->SetTitle("Expected number of background events");
    Nbackground_spectrum_2nubb->GetXaxis()->SetTitle("E>E_min");
    Nbackground_spectrum_2nubb->GetYaxis()->SetTitle("#Background");
    Nbackground_spectrum_2nubb->GetYaxis()->SetRangeUser(pow(10,-3),pow(10,6));
    Nbackground_spectrum_2nubb->SetLineColor(2);

    Nbackground_spectrum_214Bi->Draw("SAME");
    Nbackground_spectrum_214Bi->SetLineColor(7);

    Nbackground_spectrum_208Tl->Draw("SAME");
    Nbackground_spectrum_208Tl->SetLineColor(6);

    Nbackground_tot_spectrum->Add(Nbackground_spectrum_2nubb);
    Nbackground_tot_spectrum->Add(Nbackground_spectrum_208Tl);
    Nbackground_tot_spectrum->Add(Nbackground_spectrum_214Bi);
    Nbackground_tot_spectrum->Draw("SAME");


    auto legend3= new TLegend(0.75,0.6,0.89,0.89);
    legend3->AddEntry(Nbackground_spectrum_2nubb,"2nubb","l");
    legend3->AddEntry(Nbackground_spectrum_214Bi,"214Bi","l");
    legend3->AddEntry(Nbackground_spectrum_208Tl,"208Tl","l");
    legend3->AddEntry(Nbackground_tot_spectrum,"TOTAL","l");
    legend3->Draw();


    TH2F *sensitivity_spectrum = new TH2F("h2_ROI","low bound vs up bound",10,2.85,3.35,20,2.85,3.85);

    Double_t Nbackground = 0;
    Double_t Nbackground_2nubb = 0;
    Double_t Nbackground_208Tl = 0;
    Double_t Nbackground_214Bi = 0;
    Double_t Nbackground_222Rn = 0;
    Double_t sensitivity = 0;
    Double_t efficiency_0nubb = 0;
    Double_t efficiency_2nubb = 0;
    Double_t efficiency_2nubb_2MeV = 0;
    Double_t efficiency_208Tl = 0;
    Double_t efficiency_214Bi = 0;
    Double_t expectedSignalEventLimit = 0;
    int cointer = 0;

    for (int i=Nbackground_tot_spectrum->GetXaxis()->FindBin(2.85);i<=Nbackground_tot_spectrum->GetXaxis()->FindBin(3.35);i++) {
      for (int j=i+1;j<=Nbackground_tot_spectrum->GetXaxis()->FindBin(3.85);j++) {

        efficiency_0nubb += energy_spectrum_0nubb->Integral(i,j);
        efficiency_0nubb *= 1./9574121;

        efficiency_2nubb += energy_spectrum_2nubb->Integral(i,j);
        efficiency_2nubb *= 1./9574121;
        Nbackground_2nubb = ((Na*log(2))/mass_mol_150Nd)*((efficiency_2nubb*exposure)/T_2nubb_150Nd);

        efficiency_214Bi += energy_spectrum_214Bi->Integral(i,j);
        efficiency_214Bi *= 1./9574121;
        //Nbackground_214Bi = 2.3;
        //Nbackground_222Rn = 10.1;
        Nbackground_214Bi = A_214Bi*efficiency_214Bi*exposure;
        Nbackground_222Rn = Nbackground_214Bi*4.39;


        efficiency_208Tl += energy_spectrum_208Tl->Integral(i,j);
        efficiency_208Tl *= 1./9574121;
        Nbackground_208Tl = A_208Tl*efficiency_208Tl*exposure;

        Nbackground = Nbackground_2nubb+Nbackground_208Tl+Nbackground_214Bi+Nbackground_222Rn;
        if (Nbackground > 200.) {
          expectedSignalEventLimit = TMath::Sqrt(Nbackground);
        }
        else {
          expectedSignalEventLimit = WindowMethodFindExpSigEvts(Nbackground);
        }


        sensitivity = ((Na*log(2))/mass_mol_150Nd)*((efficiency_0nubb*exposure)/expectedSignalEventLimit);

        sensitivity_spectrum->SetBinContent(i+1-Nbackground_tot_spectrum->GetXaxis()->FindBin(2.85),j+2-Nbackground_tot_spectrum->GetXaxis()->FindBin(2.9),sensitivity);

      }
    }

    cout << "T12 = " << get_ROI(sensitivity_spectrum).T12_max << endl;
    cout << "mbb = " << "[" << mbb_min(get_ROI(sensitivity_spectrum).T12_max,isotope) << "," << mbb_max(get_ROI(sensitivity_spectrum).T12_max,isotope) << "]" << endl;

    double efficiency_0nubb_ROI = energy_spectrum_0nubb->Integral(Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Einf_ROI),Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Esup_ROI));
    cout << "0nubb " << "[" << get_ROI(sensitivity_spectrum).Einf_ROI << "," << get_ROI(sensitivity_spectrum).Esup_ROI << "] = " << efficiency_0nubb_ROI/9574121 << endl;


    double efficiency_2nubb_2MeV_ROI = energy_spectrum_2nubb->Integral(Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Einf_ROI),Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Esup_ROI));
    double Nbackground_2nubb_2MeV_ROI = ((Na*log(2))/mass_mol_82Se)*((efficiency_2nubb_2MeV_ROI*exposure)/T_2nubb_82Se);
    cout << "2nubb_2MeV " << "[" << get_ROI(sensitivity_spectrum).Einf_ROI << "," << get_ROI(sensitivity_spectrum).Esup_ROI << "] = " << Nbackground_2nubb_2MeV_ROI/9656794 << endl;


    double efficiency_208Tl_ROI = energy_spectrum_208Tl->Integral(Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Einf_ROI),Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Esup_ROI));
    double Nbackground_208Tl_ROI = A_208Tl*efficiency_208Tl_ROI*exposure;
    cout << "208Tl " << "[" << get_ROI(sensitivity_spectrum).Einf_ROI << "," << get_ROI(sensitivity_spectrum).Esup_ROI << "] = " << Nbackground_208Tl_ROI/9773130 << endl;


    double efficiency_214Bi_ROI = energy_spectrum_214Bi->Integral(Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Einf_ROI),Nbackground_tot_spectrum->GetXaxis()->FindBin(get_ROI(sensitivity_spectrum).Esup_ROI));
    double Nbackground_214Bi_ROI = A_214Bi*efficiency_214Bi_ROI*exposure;
    cout << "214Bi " << "[" << get_ROI(sensitivity_spectrum).Einf_ROI << "," << get_ROI(sensitivity_spectrum).Esup_ROI << "] = " << Nbackground_214Bi_ROI/9778522 << endl;

    c_sensitivity_spectrum->cd();
    TGaxis::SetMaxDigits(2);
    sensitivity_spectrum->Draw("colz TEXT");
    sensitivity_spectrum->SetTitle("Expected sensitivity");
    sensitivity_spectrum->GetXaxis()->SetTitle("Sup ROI (MeV)");
    sensitivity_spectrum->GetYaxis()->SetTitle("Inf ROI (MeV)");

  }




  string file_path_Canvas1;

  if (field && isotope == "82Se") {
    c_energy_spectrum->SaveAs("plots/energy_spectrum_with_B_82Se.pdf");
    c_efficiency_spectrum->SaveAs("plots/efficiency_spectrum_with_B_82Se.pdf");
    c_Nbackground_spectrum->SaveAs("plots/Nbackground_spectrum_with_B_82Se.pdf");
    c_sensitivity_spectrum->SaveAs("plots/sensitivity_spectrum_with_B_82Se.pdf");
  }

  else if (!field &&  isotope == "82Se") {
    c_energy_spectrum->SaveAs("plots/energy_spectrum_without_B_82Se.pdf");
    c_efficiency_spectrum->SaveAs("plots/efficiency_spectrum_without_B_82Se.pdf");
    c_Nbackground_spectrum->SaveAs("plots/Nbackground_spectrum_without_B_82Se.pdf");
    c_sensitivity_spectrum->SaveAs("plots/sensitivity_spectrum_without_B_82Se.pdf");
  }

  else if (field &&  isotope == "150Nd") {
    c_energy_spectrum->SaveAs("plots/energy_spectrum_with_B_150Nd.pdf");
    c_efficiency_spectrum->SaveAs("plots/efficiency_spectrum_with_B_150Nd.pdf");
    c_Nbackground_spectrum->SaveAs("plots/Nbackground_spectrum_with_B_150Nd.pdf");
    c_sensitivity_spectrum->SaveAs("plots/sensitivity_spectrum_with_B_150Nd.pdf");
  }

  else if (!field &&  isotope == "150Nd") {
    c_energy_spectrum->SaveAs("plots/energy_spectrum_without_B_150Nd.pdf");
    c_efficiency_spectrum->SaveAs("plots/efficiency_spectrum_without_B_150Nd.pdf");
    c_Nbackground_spectrum->SaveAs("plots/Nbackground_spectrum_without_B_150Nd.pdf");
    c_sensitivity_spectrum->SaveAs("plots/sensitivity_spectrum_without_B_150Nd.pdf");
  }


}

TH1F *henergy_sum(string isotope, bool field, string process, double xmin, double xmax, int nbins){

  string file_path;
  if (field == 1) {
    file_path = "$WORKDIR/Sensitivity/"+isotope+string("/with_B/")+process+string("/root_file_")+process+string("_");
  }
  else {
    file_path = "$WORKDIR/Sensitivity/"+isotope+string("/without_B/")+process+string("/root_file_")+process+string("_");
  }
  string file_extension = ".root";
  string file_path_new;
  string file_adress;
  string tree_adress;
  TList *list = new TList;


  TFile *file_events = new TFile("events.root","UPDATE");
  TTree *tree_events = new TTree("TreeEvents","Tree storing all selected events");


  ///Loop on all .root output files
  for (int i = 0; i < 10; ++i) {

    stringstream ss;
    ss << i;
    string str = ss.str();
    file_path_new = file_path+str+file_extension;
    file_adress = string("f")+str;
    tree_adress = string("tree")+str;
    const char* c_file_path_new = file_path_new.c_str();
    file_events = TFile::Open(c_file_path_new);
    if (file_events->IsOpen()) {
      tree_events = (TTree*)file_events->Get("calorimeter_hit");
      list->Add(tree_events);
    }
    else cout << "Failed opening root file" << endl;
  }

  TH1F * histo = new TH1F("","",30,-0.5,2);

  ///Merge all trees in "newtree"
  TFile *outputFile = TFile::Open("outputfile.root","RECREATE");
  TTree *newtree = TTree::MergeTrees(list);
  int event_counter;
  int event_counter_ytrue;
  Double_t time_difference_E;
  Double_t energy_sum;
  Double_t time_Emin;
  Double_t time_Emax;
  Double_t probability;
  Double_t length_Emin;
  Double_t length_Emax;
  Double_t minimal_energy;
  Double_t maximal_energy;
  Double_t sigma_time_Emin;
  Double_t sigma_time_Emax;

  newtree->SetBranchAddress("event_counter",&event_counter);
  if (process == "150Nd") {
    newtree->SetBranchAddress("event_counter_ytrue",&event_counter_ytrue);
  }
  newtree->SetBranchAddress("time_difference_E",&time_difference_E);
  newtree->SetBranchAddress("time_Emin",&time_Emin);
  newtree->SetBranchAddress("time_Emax",&time_Emax);
  newtree->SetBranchAddress("probability",&probability);
  newtree->SetBranchAddress("energy_sum",&energy_sum);
  newtree->SetBranchAddress("length_Emin",&length_Emin);
  newtree->SetBranchAddress("length_Emax",&length_Emax);
  newtree->SetBranchAddress("minimal_energy",&minimal_energy);
  newtree->SetBranchAddress("maximal_energy",&maximal_energy);
  newtree->SetBranchAddress("sigma_time_Emin",&sigma_time_Emin);
  newtree->SetBranchAddress("sigma_time_Emax",&sigma_time_Emax);


  ///Create histograms to be filled with reconstructed data outputs
  Long64_t nentries = newtree->GetEntries();

  int nbr_entries = 0;

  TH1F *henergy_sum = new TH1F("energy_sum","Energy sum histogram",nbins,xmin,xmax);

  for (Long64_t i=0;i<newtree->GetEntries();i++) {
    newtree->GetEntry(i);
    histo->Fill(sigma_time_Emin);
    henergy_sum->Fill(energy_sum);
    nbr_entries++;
    if (process == "0nubb") {
      _nbr_events_0nubb_ = i+1;
      final_flux << "Event # " << event_counter << endl;
    }
    if (process == "2nubb") {
      _nbr_events_2nubb_ = i+1;
    }
    if (process == "2nubb_2MeV") {
      _nbr_events_2nubb_2MeV_ = i+1;
    }
    if (process == "214Bi") {
      _nbr_events_214Bi_ = i+1;
    }
    if (process == "208Tl") {
      _nbr_events_208Tl_ = i+1;
    }
  }
  if (process == "0nubb") {
    cout << process << "  " << setprecision(9) << _nbr_events_0nubb_ << endl;
  }
  if (process == "2nubb") {
    cout << process << "  " << setprecision(9) << _nbr_events_2nubb_ << endl;
  }
  if (process == "2nubb_2MeV") {
    cout << process << "  " << setprecision(9) << _nbr_events_2nubb_2MeV_ << endl;
  }
  if (process == "214Bi") {
    cout << process << "  " << setprecision(9) << _nbr_events_214Bi_ << endl;
  }
  if (process == "208Tl") {
    cout << process << "  " << setprecision(9) << _nbr_events_208Tl_ << endl;
  }

  if (process == "2nubb_2MeV") {
    if (isotope == "82Se") {
      henergy_sum->Scale(0.0439913);
    }
    else if (isotope == "150Nd"){
      henergy_sum->Scale(0.0986686);
    }
  }
  return henergy_sum;
}

TH1F *hefficiency(string isotope, bool field, string process, double xmin, double xmax, int nbins, TH1F *histo_energy_sum = 0){

  TH1F *hefficiency = new TH1F("efficiency","Efficiency histogram",nbins,xmin,xmax);

  if (!histo_energy_sum) {
    histo_energy_sum = henergy_sum(isotope,field, process,xmin,xmax,nbins);
  }

  int n_bins = histo_energy_sum->GetNbinsX();

  for (int emin=0; emin<n_bins; emin++){
    double efficiency = 0;
    for (int bin_i=emin+1; bin_i<=n_bins; bin_i++){
      efficiency += histo_energy_sum->GetBinContent(bin_i);
    }
    // if (efficiency < 2.3){
    //   efficiency = 2.3;
    // }
    if (isotope == "82Se") {
      efficiency *= 1./pow(10,7);
    }
    if (isotope == "150Nd") {
      if (process == "0nubb") {
        efficiency *= 1./9574121;
      }
      else if (process == "2nubb"|process == "2nubb_2MeV") {
        efficiency *= 1./9656794;
      }
      else if (process == "208Tl") {
        efficiency *= 1./9773130;
      }
      else if (process == "214Bi") {
        efficiency *= 1./9778522;
      }
      else {
        cout  << process << " Unknown process" << endl;
      }
    }
    hefficiency->SetBinContent(emin+1, efficiency);
  }
  return hefficiency;
}


TH1F *hN_background(string isotope, bool field, string process, double xmin, double xmax, int nbins, TH1F *histo_energy_sum = 0,TH1F *histo_efficiency = 0){

  TH1F *hnbr_bdf = new TH1F("nbr_bdf","Nbackground histogram",nbins,xmin,xmax);

  if (!histo_energy_sum) {
    histo_energy_sum = henergy_sum(isotope, field,process,xmin,xmax,nbins);
  }

  if (!histo_efficiency) {
    histo_efficiency = hefficiency(isotope, field,process,xmin,xmax,nbins,histo_energy_sum);
  }

  int n_bins = histo_efficiency->GetNbinsX();

  double nbr_bdf = 0;
  for (int emin=1; emin<=n_bins; emin++){
    double efficiency = 0;
    efficiency = histo_efficiency->GetBinContent(emin);

    if (process == "2nubb_2MeV") {
      if (isotope == "82Se") {
        nbr_bdf = ((Na*log(2))/mass_mol_82Se)*((efficiency*exposure)/T_2nubb_82Se);
      }
      else if (isotope == "150Nd") {
        nbr_bdf = ((Na*log(2))/mass_mol_150Nd)*((efficiency*exposure)/T_2nubb_150Nd);
      }
    }

    else if (process == "2nubb") {
      if (isotope == "82Se") {
        nbr_bdf = ((Na*log(2))/mass_mol_82Se)*((efficiency*exposure)/T_2nubb_82Se);
      }
      else if (isotope == "150Nd") {
        nbr_bdf = ((Na*log(2))/mass_mol_150Nd)*((efficiency*exposure)/T_2nubb_150Nd);
      }
    }

    else if (process == "214Bi") {
      nbr_bdf = A_214Bi*efficiency*exposure;
    }
    else if (process == "208Tl") {
      nbr_bdf = A_208Tl*efficiency*exposure;
    }
    else {
      cout << "Unknown process" << endl;
    }
    hnbr_bdf->SetBinContent(emin,nbr_bdf);

  }
  return hnbr_bdf;
}

double WindowMethodFindExpSigEvts(Double_t B){

  //Find S using CL(S+B)/CL(B) = 0.1 with N_Obs = B
  Double_t Likelihood = 1.0;
  Double_t S = 0;
  Int_t NEvts = (int)B;

  while(Likelihood > 0.1){
    S += 0.001;
    Float_t CLsb = 0; Float_t CLb = 0;
    for (Int_t i=0; i<=NEvts; i++) {
      CLsb+=TMath::Poisson(i,S+B);
      CLb+=TMath::Poisson(i,B);
    }
    Likelihood =  CLsb/CLb;
  }
  return S;
}

Double_t sensitivity_FC(string isotope, int bin_emin, int bin_emax, TH1F *histo_energy_0nubb, TH1F *histo_tot) {
  Double_t sensitivity = 0;
  Double_t efficiency = 0;
  TFeldmanCousins f;
  Double_t Nobserved = 0;

  Double_t uplim = 0;
  Double_t lowlim = 0;

  Double_t Nbackground = 0;
  Nbackground = histo_tot->GetBinContent(bin_emin);
  Double_t fractpart = modf (Nbackground , &Nobserved);

  efficiency += histo_energy_0nubb->Integral(bin_emin,bin_emax);

  uplim = f.CalculateUpperLimit(Nobserved, Nbackground);

  if (isotope == "82Se") {
    sensitivity = ((Na*log(2))/mass_mol_82Se)*((efficiency*exposure)/uplim);
  }
  else {
    sensitivity = ((Na*log(2))/mass_mol_150Nd)*((efficiency*exposure)/uplim);
  }
  cout << "efficiency = " << efficiency << endl;
  cout << "sensitivity = " << sensitivity << endl;
  cout << "Nbackground = " << Nbackground << endl;
  cout << "Nobserved = " << Nobserved << endl;
  return sensitivity;
}


// Double_t sensitivity_Poisson(string isotope, int bin_emin, int bin_emax, TH1F *histo_energy_0nubb, TH1F *histo_tot) {
//   Double_t Nbackground = histo_tot->Integral(bin_emin,bin_emax);
//   Double_t expectedSignalEventLimit = WindowMethodFindExpSigEvts(Nbackground);


//   // for (int emin=bin_emin; emin<bin_emax; emin++){
//   //   double efficiency = 0;
//   //   for (int bin_i=emin+1; bin_i<=bin_emax; bin_i++){
//   //     efficiency += histo_energy_sum->GetBinContent(bin_i);
//   //   }
//   //   efficiency *= 1./pow(10,7);
//   //   hefficiency->SetBinContent(emin, efficiency);
//   // }

//   efficiency += histo_energy_0nubb->Integral(bin_emin,bin_emax);

//   if (isotope == "82Se") {
//     Double_t sensitivity = ((Na*log(2))/mass_mol_82Se)*((efficiency*exposure)/expectedSignalEventLimit);
//   }
//   else {
//     Double_t sensitivity = ((Na*log(2))/mass_mol_150Nd)*((efficiency*exposure)/expectedSignalEventLimit);
//   }
//   // cout << "efficiency = " << efficiency << endl;
//   // cout << "sensitivity = " << sensitivity << endl;
//   // cout << "Nbackground = " << Nbackground << endl;
//   return sensitivity;
// }


TH2F *h2sensitivity(string isotope, bool field, double i_min, double i_max, double j_min, double j_max, double xmin, double xmax, int nbins,  TH1F *histo_energy_sum_0nubb = 0, TH1F *histo_energy_sum_2nubb_2MeV = 0, TH1F *histo_energy_sum_208Tl = 0, TH1F *histo_energy_sum_214Bi = 0,  TCanvas *c_sensitivity_spectrum = 0){

  if (!histo_energy_sum_0nubb) {
    histo_energy_sum_0nubb = henergy_sum(isotope,field,"0nubb",xmin,xmax,nbins);
  }
  if (!histo_energy_sum_2nubb_2MeV) {
    histo_energy_sum_2nubb_2MeV = henergy_sum(isotope,field,"2nubb_2MeV",xmin,xmax,nbins);
  }
  if (!histo_energy_sum_208Tl) {
    histo_energy_sum_208Tl = henergy_sum(isotope,field,"208Tl",xmin,xmax,nbins);
  }
  if (!histo_energy_sum_214Bi) {
    histo_energy_sum_214Bi = henergy_sum(isotope,field,"214Bi",xmin,xmax,nbins);
  }

  c_sensitivity_spectrum = new TCanvas("canvas4","canvas4");
  TH2F *sensitivity_spectrum = new TH2F("h2_ROI","low bound vs up bound",10,i_min,i_max,20,j_min,j_max);

  Double_t Nbackground = 0;
  Double_t Nbackground_2nubb_2MeV = 0;
  Double_t Nbackground_208Tl = 0;
  Double_t Nbackground_214Bi = 0;
  Double_t sensitivity = 0;
  double efficiency_0nubb = 0;
  double efficiency_2nubb_2MeV = 0;
  double efficiency_208Tl = 0;
  double efficiency_214Bi = 0;
  double expectedSignalEventLimit = 0;

  for (int i=histo_energy_sum_0nubb->GetXaxis()->FindBin(i_min);i<=histo_energy_sum_0nubb->GetXaxis()->FindBin(i_max);i++) {
    for (int j=histo_energy_sum_0nubb->GetXaxis()->FindBin(j_min);j<=histo_energy_sum_0nubb->GetXaxis()->FindBin(j_max);j++) {

      efficiency_0nubb += histo_energy_sum_0nubb->Integral(i,j);
      if (isotope == "82Se") {
        efficiency_0nubb *= 1./pow(10,7);
      }
      else if (isotope == "150Nd") {
        efficiency_0nubb *= 1./9574121;
      }

      efficiency_2nubb_2MeV += histo_energy_sum_2nubb_2MeV->Integral(i,j);
      if (isotope == "82Se") {
        efficiency_2nubb_2MeV *= 1./pow(10,7);
        Nbackground_2nubb_2MeV = ((Na*log(2))/mass_mol_82Se)*((efficiency_2nubb_2MeV*exposure)/T_2nubb_82Se);
      }
      else if (isotope == "150Nd") {
        efficiency_2nubb_2MeV *= 1./9574121;
        Nbackground_2nubb_2MeV = ((Na*log(2))/mass_mol_150Nd)*((efficiency_2nubb_2MeV*exposure)/T_2nubb_150Nd);
      }

      efficiency_214Bi += histo_energy_sum_214Bi->Integral(i,j);
      if (isotope == "82Se") {
        efficiency_2nubb_2MeV *= 1./pow(10,7);
      }
      else if (isotope == "150Nd") {
        efficiency_2nubb_2MeV *= 1./9574121;
      }
      Nbackground_214Bi = A_214Bi*efficiency_214Bi*exposure;

      efficiency_208Tl += histo_energy_sum_208Tl->Integral(i,j);
      if (isotope == "82Se") {
        efficiency_2nubb_2MeV *= 1./pow(10,7);
      }
      else if (isotope == "150Nd") {
        efficiency_2nubb_2MeV *= 1./9574121;
      }
      Nbackground_208Tl = A_208Tl*efficiency_208Tl*exposure;

      Nbackground = Nbackground_2nubb_2MeV+Nbackground_208Tl+Nbackground_214Bi;

      expectedSignalEventLimit = WindowMethodFindExpSigEvts(Nbackground);

      if (isotope == "82Se") {
        sensitivity = ((Na*log(2))/mass_mol_82Se)*((efficiency_0nubb*exposure)/expectedSignalEventLimit);
      }
      else if (isotope == "150Nd") {
        sensitivity = ((Na*log(2))/mass_mol_150Nd)*((efficiency_0nubb*exposure)/expectedSignalEventLimit);
      }

      sensitivity_spectrum->SetBinContent(i+1-histo_energy_sum_0nubb->GetXaxis()->FindBin(i_min),j+2-histo_energy_sum_0nubb->GetXaxis()->FindBin(j_min),sensitivity);
    }
  }

  c_sensitivity_spectrum->cd();
  TGaxis::SetMaxDigits(2);
  sensitivity_spectrum->Draw("colz TEXT");
  sensitivity_spectrum->SetTitle("Expected sensitivity");
  sensitivity_spectrum->GetXaxis()->SetTitle("Inf ROI (MeV)");
  sensitivity_spectrum->GetYaxis()->SetTitle("Sup ROI (MeV)");

  return sensitivity_spectrum;
}


double ErrorStatEfficiency(string process, string isotope, double efficiency){
  double statistic_error = 0;

  if (isotope == "82Se") {
    statistic_error = sqrt(efficiency/pow(10,7));
  }
  if (isotope == "150Nd") {
    if (process == "0nubb") {
      statistic_error = sqrt(efficiency/9574121);
    }
    if (process == "2nubb_2MeV") {
      statistic_error = sqrt(efficiency/9656794);
    }
    if (process == "2nubb") {
      statistic_error = sqrt(efficiency/9656794);
    }
    else if (process == "214Bi") {
      statistic_error = sqrt(efficiency/9778522);
    }
    else if (process == "208Tl") {
      statistic_error = sqrt(efficiency/9773130);
    }
  }
  return statistic_error;
}

double ErrorStatbdf(string process, string isotope, double nbr_bdf, double efficiency){
  double statistic_error = 0;

  statistic_error = (ErrorStatEfficiency(process, isotope, efficiency)*nbr_bdf)/efficiency;
  return statistic_error;
}

search_ROI get_ROI(TH2F *histo_demie_vie = 0){
  search_ROI bornes;
  double demie_vie_max = 0;
  double Einf_ROI;
  double Esup_ROI;

  for (int i=0; i<histo_demie_vie->GetNbinsX(); i++) {
    double Einf = histo_demie_vie->GetXaxis()->GetBinLowEdge(i+1);
    for (int j=i+1; j<histo_demie_vie->GetNbinsY()-1; j++) {
      double Esup = histo_demie_vie->GetYaxis()->GetBinLowEdge(j+1);
      if (demie_vie_max < histo_demie_vie->GetBinContent(i+1, j+1)) {
        demie_vie_max = histo_demie_vie->GetBinContent(i+1, j+1);
        bornes.Einf_ROI = Einf;
        bornes.Esup_ROI = Esup;
        bornes.T12_max = demie_vie_max;
      }
    }
  }
  return bornes;
}

double mbb_min(double T12_max, string isotope){
  double matrix_element;
  double phase_space;
  double mbb;
  if (isotope == "82Se") {
    matrix_element = 5.4;
    phase_space = 10.16*pow(10,-15);
  }
  else if (isotope == "150Nd") {
    matrix_element = 5.6;
    phase_space = 63.03*pow(10,-15);
  }
  mbb = (511*pow(10,3)/(pow(1.27,2)*matrix_element))*sqrt(1/(T12_max*phase_space));
  return mbb;
}

double mbb_max(double T12_max, string isotope){
  double matrix_element;
  double phase_space;
  if (isotope == "82Se") {
    matrix_element = 2.79;
    phase_space = 10.16*pow(10,-15);
  }
  else if (isotope == "150Nd") {
    matrix_element = 1.71;
    phase_space = 63.03*pow(10,-15);
  }
  mbb = (511*pow(10,3)/(pow(1.27,2)*matrix_element))*sqrt(1/(T12_max*phase_space));
  return mbb;
}
