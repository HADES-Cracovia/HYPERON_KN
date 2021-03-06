#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <format_histograms.C>

bool is_in(int n, int signals[])
{
  size_t myArraySize = sizeof(signals) / sizeof(int);
  int *end = signals + myArraySize;
  // find the value 0:
  int *result = std::find(signals, end, n);
  if(signals[myArraySize]==n)
    return 1;

  if (result != end)
    return 1;
  else
    return 0;
}

TH1F* renolmalize(TH1F* hist, double scale)
{
  TH1F* temp =(TH1F*)hist->Clone("hist_renolmalize");
  temp->Rebin(scale);
  temp->Scale(1/scale);
  return temp;
}

int sum_background(TH1F* &hist, TH1F* back1, TH1F* back2)
{
  if(back1->GetNbinsX()!=back2->GetNbinsX())
    {
      cout<<"Wrong amout of bins in background histograms!"<<endl;
      return 0;
    }
  if(hist->GetNbinsX()!=back2->GetNbinsX())
    {
      cout<<"Wrong amout of bins in sum histogram!"<<endl;
      return 0;
    }

  for (Int_t j=1; j<back1->GetNbinsX()+1; ++j)
    {
      //cout<<j;
      hist->SetBinContent(j, 2*TMath::Sqrt(back1->GetBinContent(j)*back2->GetBinContent(j)));
      //hist->SetBinError(j, TMath::Sqrt( back1->GetBinError(j)*back1->GetBinError(j)+back2->GetBinError(j)*back2->GetBinError(j) ));
    }
  return 1;
}

void sethist(TH1F* hist, int color=1, int rebin=1, int style=1, double scale=1)
{
  if(color!=0 && color!=10)//avoid white color
    hist->SetLineColor(color);
  else
    hist->SetLineColor(40);
  
  hist->Rebin(rebin);
  double binW=hist->GetBinWidth(10); //After re-binning
  hist->SetLineStyle(style);
  hist->Scale(scale/binW);
  hist->GetYaxis()->SetTitle("#frac{d #sigma}{d M} [#frac{#mu b}{MeV}]");
  //hist->Smooth();
}

int calcBg(TH1F* &histBG, TH1F** hists, int channels, int signal[])
{
  if(!is_in(0,signal))
    histBG=(TH1F*)hists[0]->Clone("background");
  else
    histBG=(TH1F*)hists[1]->Clone("background");

  histBG->Reset();
  
  for(int i=0;i<channels;i++)
    {
      if(is_in(i,signal))
	{
	  //cout<<"channel "<<i<<" is a signal channel"<<endl;
	  continue;
	}
      else
	histBG->Add(hists[i]);
    }
  return 1;
}

int sumAll(TH1F* &histBG, TH1F** hists, int channels, int signal[])
{
  histBG=(TH1F*)hists[0]->Clone("background");
  histBG->Reset();
  
  for(int i=0;i<channels;i++)
    histBG->Add(hists[i]);
    
  return 1;
}

int sumSignals(TH1F* &histSum, TH1F** hists, int channels, int signal[])
{
  histSum=(TH1F*)hists[0]->Clone("SumOfSignals");
  histSum->Reset();
  
  for(int i=0;i<channels;i++)
    if(is_in(i,signal))
      {
	//cout<<"channel "<<i<<" is a signal channel"<<endl;
	histSum->Add(hists[i],1);
      }
  return 1;
}

int L_CM_pictures_M3()
{
  const int fn=14;//number of files in "infile"
  const int signal_ch[]={1,8,9};//list of signal channels in "infile" file, starting from 0
  const int ur_background[]={10,11,12,13};
  std::ifstream infile("files_list_k.dat");//list of histograms
  std::string file;
  //std::string directory="/lustre/nyx/hades/user/iciepal/Lambda1520_ic/";//directory comon for all files
  std::string directory="./results_newGeometryM3/";//directory comon for all files
  int n=0;
  //write everything to file
  //TFile *MyFile = new TFile("output.root","recreate");
  //if ( MyFile->IsOpen() )
  //printf("File opened successfully\n");

  const float scale_factor=50000*1000;
  float scale[]={
    1840./scale_factor/10,//channel 01 //10 times bgger statistics
    64*0.0085*1.35/137/scale_factor,//channel 48
    64/5.34/scale_factor,//channel 49
    300/scale_factor,//channel 04
    300/scale_factor,//channel 02
    43/scale_factor,//channel 40
    10.0/scale_factor /10,//channel 42 //10 times bgger statistics
    7.0/scale_factor/10,//channel 44 //10 times bgger statistics
    30*0.00054*1.35/137/scale_factor,//channel 50
    53.15*0.012*1.35*1/137/scale_factor,//channel 52
    1.242e-1/scale_factor / 8.065167e-10,//channel 55 weight 8.065167e-10 set in PLUTO
    24e-2/scale_factor / 8.065167e-10,//channel 56 weight 8.065167e-10 set in PLUTO
    29e-4/scale_factor/ 8.065167e-10,//channel 57 weight 8.065167e-10 set in PLUTO
    18.6e-4/scale_factor/ 8.065167e-10//channel 58 weight 8.065167e-10 set in PLUTO
  };//ub
  
  TH1F *hL1520massDistZLpi0[fn];
  TH1F *hL1520massFTDistZLpi0[fn];
  TH1F *hL1520massFinalRLpi0[fn];
  TH1F *hL1520massFTFinalRLpi0[fn];
  TH1F *hL1520massFinalpi0[fn];
  TH1F *hL1520massFTFinalpi0[fn];
  
  TH1F *hL1520massFinalRLpi0_L[fn];
  TH1F *hL1520massDistZLRLpi0_L[fn];
  
  TH1F *hDLmassFinalRL_L[fn];
  TH1F *hDLmassDistZL_L[fn];
  TH1F *hDLmassDistZL[fn];
  TH1F *hDLmassFTDistZL[fn];
  TH1F *hL1520massDistZLpi0_emem[fn];
  TH1F *hL1520massDistZLpi0_epep[fn];
  TH1F *hDLmassDistZL_emem[fn];
  TH1F *hDLmassDistZL_epep[fn];


  
  TH1F *hL1520mass_background;
  TH1F *hL1520mass_background_L;
  TH1F *hDLmass_sum;
  TH1F *hDLmass_background;
  TH1F *hDLmass_background_L;
  TH1F *hL1520mass_background_epep;
  TH1F *hDLmass_background_epep;
  TH1F *hL1520mass_background_emem;
  TH1F *hDLmass_background_emem;
  TH1F *hL1520mass_real_backgroud;
  TH1F *hDLmass_sum_background;
  TH1F *hL1520mass_sum_CBbackground;  

  TH1F *hL1520mass_sum_all_signals;
  TH1F *hDLmass_sum_all_signals;
  TH1F *hDLmass_sum_right_vertex;

  TH1F *hDLmass_delta;
  TH1F *hL1520_delta;
  //read file wth list of histograms
  if(!infile)
    {
      cout<<"Can't open the file"<<endl;
      return 1;
    }

  //loop over files names and clone histograms*********************************************************
  while (std::getline(infile,file))
    {
      cout<<"file number: "<<n<<" ";
      cout<<file<<endl;
      
      const char* file_name=(directory+file).c_str();
      TFile* hist_file=new TFile(file_name);
      if(hist_file->IsZombie())
	{
	  cout<<"can't open a file!"<<endl<<file_name<<endl;
	  n++;
	  continue;
	}
      hist_file->cd();
      
      hL1520massDistZLpi0[n]= (TH1F*)hist_file->Get("hL1520massDistZLpi0")->Clone();
      //hL1520massFTDistZLpi0[n]= (TH1F*)hist_file->Get("hL1520massFTDistZLpi0")->Clone();
      hL1520massFinalRLpi0[n]= (TH1F*)hist_file->Get("hL1520massFinalRLpi0")->Clone();
      hL1520massFTFinalRLpi0[n]= (TH1F*)hist_file->Get("hL1520massFTFinalRLpi0")->Clone();
      hL1520massFinalpi0[n]= (TH1F*)hist_file->Get("hL1520massFinalpi0")->Clone();
      hL1520massFTFinalpi0[n]= (TH1F*)hist_file->Get("hL1520massFTFinalpi0")->Clone();
      
      hL1520massFinalRLpi0_L[n]=(TH1F*)hist_file->Get("hL1520massFinalRLpi0_L")->Clone();
      hL1520massDistZLRLpi0_L[n]=(TH1F*)hist_file->Get("hL1520massDistZLRLpi0_L")->Clone();

      hDLmassFinalRL_L[n]=(TH1F*)hist_file->Get("hDLmassFinalRL_L")->Clone();
      hDLmassDistZL_L[n]=(TH1F*)hist_file->Get("hDLmassDistZL_L")->Clone();
      hDLmassDistZL[n]= (TH1F*)hist_file->Get("hDLmassDistZL")->Clone();
      hDLmassFTDistZL[n]= (TH1F*)hist_file->Get("hDLmassFTDistZL")->Clone();

      hL1520massDistZLpi0_epep[n]= (TH1F*)hist_file->Get("hL1520massDistZLpi0_epep")->Clone();
      hL1520massDistZLpi0_emem[n]= (TH1F*)hist_file->Get("hL1520massDistZLpi0_emem")->Clone();

      hDLmassDistZL_epep[n]= (TH1F*)hist_file->Get("hDLmassDistZL_epep")->Clone();
      hDLmassDistZL_emem[n]= (TH1F*)hist_file->Get("hDLmassDistZL_emem")->Clone();
      
      n++;
    }
  //end of reading histograms*************************************************************************

  //Print resoults************************************************************************************
  for(int k=0; k<n;k++)
    {
      int bins=10;
      //set colors
          
      sethist(hL1520massFinalpi0[k],k,bins,1,scale[k]);
      sethist(hL1520massFTFinalpi0[k],k,bins,1,scale[k]);
      sethist(hL1520massFinalRLpi0_L[k],k,bins,2,scale[k]);
      sethist(hL1520massDistZLRLpi0_L[k],k,bins,2,scale[k]);
      sethist(hL1520massDistZLpi0[k],k,bins,1,scale[k]);
      //sethist(hL1520massFTDistZLpi0[k],k,bins,1,scale[k]);
      
      sethist(hDLmassFinalRL_L[k],k,bins,2,scale[k]);
      sethist(hDLmassDistZL_L[k],k,bins,2,scale[k]);
      sethist(hDLmassDistZL[k],k,bins,1,scale[k]);
      sethist(hDLmassFTDistZL[k],k,bins,1,scale[k]);

      sethist(hL1520massDistZLpi0_emem[k],k,bins,1,scale[k]);
      sethist(hL1520massDistZLpi0_epep[k],k,bins,1,scale[k]);
      sethist(hDLmassDistZL_emem[k],k,bins,1,scale[k]);
      sethist(hDLmassDistZL_epep[k],k,bins,1,scale[k]);

      //sum FW with HADES
      
      //hL1520massDistZLpi0[k]->Add(hL1520massFTDistZLpi0[k]);
      hL1520massFinalRLpi0[k]->Add(hL1520massFTFinalRLpi0[k]);
      hL1520massFinalpi0[k]->Add(hL1520massFTFinalpi0[k]);
      hDLmassDistZL[k]->Add(hDLmassFTDistZL[k]);
      
    }

  cout<<"attempt to create sum of background channels"<<endl;
  //calculate sum of all background channels
  
  calcBg(hL1520mass_background, hL1520massDistZLpi0, fn, signal_ch);
  calcBg(hDLmass_background, hDLmassDistZL, fn, signal_ch);
  //calcBg(hL1520mass_background_emem, hL1520massDistZLpi0_emem,fn,signal_ch);
  //calcBg(hL1520mass_background_epep, hL1520massDistZLpi0_epep,fn,signal_ch);
    
  //sum all channels
  
  //sumAll(hL1520mass_background, hL1520massDistZLpi0, fn, signal_ch);
  //sumAll(hDLmass_background, hDLmassDistZL, fn, signal_ch);
  sumAll(hL1520mass_background_emem, hL1520massDistZLpi0_emem,fn,signal_ch);
  sumAll(hL1520mass_background_epep, hL1520massDistZLpi0_epep,fn,signal_ch);
  sumAll(hDLmass_background_emem, hDLmassDistZL_emem,fn,signal_ch);
  sumAll(hDLmass_background_epep, hDLmassDistZL_epep,fn,signal_ch);
  sumAll(hDLmass_sum,hDLmassDistZL,fn,signal_ch);
  sumAll(hDLmass_sum_right_vertex,hDLmassDistZL_L,fn,signal_ch);
  
  //hL1520mass_real_backgroud=(TH1F*)hL1520mass_background->Clone("hL1520mass_real_backgroud");
  //hL1520mass_real_backgroud->Reset();
  //sum_background(hL1520mass_real_backgroud,hL1520mass_background_emem,hL1520mass_background_epep);
  hDLmass_sum_background=(TH1F*)hDLmass_background_epep->Clone("hDLmass_sum_background");
  hDLmass_sum_background->Reset();
  hL1520mass_sum_CBbackground=(TH1F*)hL1520mass_background_epep->Clone("hDLmass_sum_background");
  hL1520mass_sum_CBbackground->Reset();
  sum_background(hDLmass_sum_background,hDLmass_background_epep,hDLmass_background_epep);
  sum_background(hL1520mass_sum_CBbackground,hL1520mass_background_epep,hL1520mass_background_emem);
  sumSignals(hL1520mass_sum_all_signals,hL1520massDistZLpi0,fn,signal_ch);
  sumSignals(hDLmass_sum_all_signals,hDLmassDistZL,fn,signal_ch);
  sumSignals(hDLmass_delta,hDLmassDistZL,fn,ur_background);
  sumSignals(hL1520_delta,hL1520massDistZLpi0,fn,ur_background);
  
  cout<<"all histograms set"<<endl;
  
  TLegend *legend = new TLegend(0.1,0.2,0.99,0.9);
  legend->AddEntry(hL1520massFinalpi0[1],"p K+ #Lambda(1520)[#Lambda(1115) e+ e-] 130#mub","l");
  legend->AddEntry(hL1520massFinalRLpi0_L[1],"p K+ #Lambda(1520)[#Lambda(1115) e+ e-] 130#mub - true signal","l");
  legend->AddEntry(hL1520massFinalpi0[5],"p K+ #Lambda(1115) #pi^{0}  100#mub","l");
  legend->AddEntry(hL1520massFinalpi0[6],"p K+ #Lambda(1115) 2#pi^{0} 20#mub","l");
  legend->AddEntry(hL1520massFinalpi0[7],"p K+ #Lambda(1115) 3#pi^{0} 7#mub","l");
  legend->AddEntry(hL1520massFinalpi0[0],"p p #pi^{+} #pi^{-} #pi^{0} 1840#mub","l");
  legend->AddEntry(hL1520massFinalpi0[4],"p p #pi^{+} #pi^{-} 2#pi^{0} 300#mub","l");
  legend->AddEntry(hL1520massFinalpi0[3],"p n 2#pi^{+} #pi^{-} #pi^{0} 200#mub","l");
  legend->AddEntry(hL1520massFinalpi0[2],"L1520 decays 130#mub","l");
  legend->AddEntry(hL1520massFinalpi0[8],"p K+ #Lambda(1405)[#Lambda(1115) e+ e-] 173#mub","l");
  legend->AddEntry(hL1520massFinalpi0[9],"p K+ #Sigma(1385)[#Lambda(1115) e+ e-] 64#mub","l");
  legend->AddEntry(hL1520massFinalpi0[10],"p #Delta^{+}[p e+ e-] #pi^{+} #pi^{-} #mub","l");
  legend->AddEntry(hL1520massFinalpi0[11],"p #Delta^{0}[n e+ e-] #pi^{+} #pi^{-} #pi^{-} #mub","l");
  legend->AddEntry(hL1520massFinalpi0[12],"#Delta^{+}[p e+ e-] #Lambda K^{+} #mub","l");
  legend->AddEntry(hL1520massFinalpi0[13],"#Delta^{0}[n e+ e-] #Lambda K^{+} #pi^{+} #mub","l");
  //Draw all things
  TCanvas *cPictures = new TCanvas("cPictures","cPictures");

  cPictures->Divide(3,2);
  double ymin=1e-12; //min value for y axis
  double ymax=1e-4;
    
  cPictures->cd(1);
  gPad->SetLogy();
  hL1520massDistZLpi0[0]->GetYaxis()->SetRangeUser(ymin,10e4);
  hL1520massDistZLRLpi0_L[0]->GetYaxis()->SetRangeUser(ymin,10e4);
  for(int x=0;x<n;x++)
    {
      //if(x==signal_ch-1)
	hL1520massDistZLpi0[x]->Draw("same");
	hL1520massDistZLRLpi0_L[x]->Draw("same");
    }
  hL1520mass_background->Draw("same");
  hL1520mass_background->SetLineWidth(2);
  hL1520mass_background->SetLineColor(kRed);

  cPictures->cd(2);
  gPad->SetLogy();
  hDLmassDistZL[0]->GetYaxis()->SetRangeUser(ymin,ymax);
  hDLmassDistZLRL_L[0]->GetYaxis()->SetRangeUser(ymin,ymax);
  for(int x=0;x<n;x++)
    {
      hDLmassDistZL[x]->Draw("same");
      //hDLmassDistZLRL_L[x]->Draw("same");
    }
  hDLmass_background->Draw("same");
  hDLmass_sum->Draw("same");
  
  cPictures->cd(3);
  gPad->SetLogy();
  hDLmass_background_epep->Draw("same");
  hDLmass_background_epep->SetLineColor(kBlue);
  hDLmass_background_emem->Draw("same");
  hDLmass_background_emem->SetLineColor(kRed);
  hDLmass_sum_background->Draw("same");
  hDLmass_sum->Draw("same");
  hDLmass_sum->SetLineWidth(2);

  cPictures->cd(4);
  gPad->SetLogy();
  hDLmassDistZL_emem[0]->GetYaxis()->SetRangeUser(ymin,ymax);
  for(int x=0; x<n; x++)
    {
      //hL1520massDistZLpi0_epep[x]->Draw("same");
      hDLmassDistZL_emem[x]->Draw("same");
    }
  hDLmass_background_emem->SetLineStyle(2);
  hDLmass_background_emem->Draw("same");
  
  
  cPictures->cd(5);
  gPad->SetLogy();
  hDLmassDistZL_epep[0]->GetYaxis()->SetRangeUser(ymin,ymax);
  for(int x=0; x<n; x++)
    {
      hDLmassDistZL_epep[x]->Draw("same");
      //hL1520massDistZLpi0_emem[x]->Draw("same");
    }
  hDLmass_background_epep->SetLineStyle(2);
  hDLmass_background_epep->Draw("same");
   
  
  cPictures->cd(6);
  legend->Draw();

  //Including formfactor for di-lepton
  //**********************************
  const TF1* fFormFactor=new TF1("fFormFactor","(1/(1-(x/700)**2))**2",0,420); //in MeV
  TH1F* hDLmass_FF=(TH1F*)hDLmass_sum_all_signals->Clone("hDLmass_FF");
  TH1F* hL1520mass_FF=(TH1F*)hL1520mass_sum_all_signals->Clone("hDLmass_FF");
  double ff_const;//ratio between resoult with and without FF
  double mianownik;
  double licznik;
  double m_min=140;
  double m_max=420;

  hDLmass_FF->Multiply(fFormFactor);
  licznik=hDLmass_FF->Integral(hDLmass_FF->FindBin(m_min),hDLmass_FF->FindBin(m_max));
  mianownik=hDLmass_sum_all_signals->Integral(hDLmass_sum_all_signals->FindBin(m_min),hDLmass_sum_all_signals->FindBin(m_max));
  ff_const=licznik/mianownik;
  cout<<"************"<<endl<<"FF scaling factor="<<ff_const<<endl<<"***********"<<endl;
  hL1520mass_FF->Scale(ff_const);
  //**********************************

  gStyle->SetOptStat(0);
  
  //L1520 SPECTRUM
  
  TCanvas* cFinalL1520=new TCanvas("cFinalL1520","cFinalL1520");
  cFinalL1520->cd();

  hL1520mass_background->Add(hL1520_delta,-1);
  hL1520mass_background->GetXaxis()->SetTitle("M_{#Lambda^{0} e^{+} e^{-}} [MeV]");
  
  TH1F* sum_renormalize=renolmalize(hL1520mass_background,4);
  TH1F* CM_background=renolmalize(hL1520mass_sum_CBbackground,4);
  //hL1520mass_FF->Draw("same");
  //hL1520mass_FF->SetLineStyle(2);
  //hL1520mass_FF->SetLineColor(kGreen);
  hL1520mass_sum_all_signals->GetXaxis()->SetTitle("M_{#Lambda^{0} e^{+} e^{-}} [MeV]");
  hL1520mass_sum_all_signals->Draw("same");
  hL1520mass_sum_all_signals->SetLineWidth(2);
  hL1520mass_sum_all_signals->SetLineColor(kGreen);
  hL1520mass_sum_all_signals->SetTitle("Y #rightarrow #Lambda e^{+}e^{-}");
  hL1520mass_sum_all_signals->GetXaxis()->SetRangeUser(1200,1800);
  hL1520mass_sum_all_signals->GetXaxis()->SetLabelFont(42);
  hL1520mass_sum_all_signals->GetXaxis()->SetNdivisions(505);
  hL1520mass_sum_all_signals->GetXaxis()->SetLabelSize(0.05);
  hL1520mass_sum_all_signals->GetXaxis()->SetTitleSize(0.05);
  hL1520mass_sum_all_signals->GetXaxis()->SetTitleOffset(1.1);
  hL1520mass_sum_all_signals->GetXaxis()->SetTitleFont(42);
  hL1520mass_sum_all_signals->GetYaxis()->SetLabelFont(42);
  hL1520mass_sum_all_signals->GetYaxis()->SetLabelSize(0.05);
  hL1520mass_sum_all_signals->GetYaxis()->SetTitleOffset(1.15);
  hL1520mass_sum_all_signals->GetYaxis()->SetTitleSize(0.05);
  hL1520mass_sum_all_signals->GetYaxis()->SetTitleFont(42);
  
  sum_renormalize->Draw("psame");
  sum_renormalize->Smooth();
  //CM_background->Draw("same");
  //CM_background->SetLineColor(kBlue);
  hL1520massDistZLpi0[1]->Draw("same");
  hL1520massDistZLpi0[8]->Draw("same");
  hL1520massDistZLpi0[9]->Draw("same");
  hL1520_delta->Draw("same");
  hL1520_delta->SetLineWidth(2);

  format_l1520(hL1520massDistZLpi0[1]);
  format_s1385(hL1520massDistZLpi0[8]);
  format_l1405(hL1520massDistZLpi0[9]);
  format_cb(sum_renormalize);

  TLegend *leg = new TLegend(0.7,0.6,0.9,0.95,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.027);
  TLegendEntry *entry=leg->AddEntry(hL1520mass_sum_all_signals,"Y #rightarrow #Lambda e^{+}e^{-}","lpf");
  entry->SetFillStyle(1001);
  entry->SetTextFont(42);
  entry=leg->AddEntry(sum_renormalize,"CB","lpf");
  entry->SetTextFont(42);
  entry=leg->AddEntry(hL1520massDistZLpi0[1],"#Lambda(1520) #rightarrow #Lambda e^{+}e^{-}","lpf");
  entry->SetTextFont(42);
  entry=leg->AddEntry(hL1520massDistZLpi0[9],"#Lambda(1405) #rightarrow #Lambda e^{+}e^{-}","lpf");
  entry->SetFillStyle(1001);
  entry->SetTextFont(42);
  entry=leg->AddEntry(hL1520massDistZLpi0[8],"#Sigma(1385) #rightarrow #Lambda e^{+}e^{-}","lpf");
  entry->SetTextFont(42);
  leg->Draw();
  
  
  cFinalL1520->Update();

  //dI-LEPTON SPECTRUM
  
  TCanvas* cFinalDL_cb=new TCanvas("cFinalDL_cb","cFinalDL_cb");
  gPad->SetLogy();
  TH1F* hDLmass_CB=(TH1F*)hDLmass_sum->Clone("CBbackground");
  hDLmass_CB->Add(hDLmass_sum_right_vertex,-1);
  hDLmass_CB->Smooth();
  TH1F* hDLmass_background_renolmalize=renolmalize(hDLmass_CB,2);
  TH1F* hDLmass_sum_background_re=renolmalize(hDLmass_sum_background,2);
  TH1F* hDLmass_sum_renolmalize=renolmalize(hDLmass_sum,2);

  hDLmass_sum_right_vertex->Add(hDLmass_delta,-1);
  hDLmass_sum_right_vertex->Add(hDLmass_sum_all_signals,-1);
  hDLmass_sum_right_vertex->GetXaxis()->SetTitle("M^{inv}_{e^{+} e^{-}} #left[MeV c^{-2} #right]");

  //set graphical properties for first histogram on canvas
  hDLmass_sum_right_vertex->SetTitle("Y #rightarrow #Lambda e^{+}e^{-}");
  hDLmass_sum_right_vertex->GetXaxis()->SetRangeUser(0,500);
  hDLmass_sum_right_vertex->GetXaxis()->SetLabelFont(42);
  hDLmass_sum_right_vertex->GetXaxis()->SetNdivisions(505);
  hDLmass_sum_right_vertex->GetXaxis()->SetLabelSize(0.05);
  hDLmass_sum_right_vertex->GetXaxis()->SetTitleSize(0.05);
  hDLmass_sum_right_vertex->GetXaxis()->SetTitleOffset(1.1);
  hDLmass_sum_right_vertex->GetXaxis()->SetTitleFont(42);
  hDLmass_sum_right_vertex->GetYaxis()->SetTitle("#frac{d#sigma}{dM_{e^{+}e^{-}}} #left[#frac{#mub}{MeV c^{-2}} #right]");
  hDLmass_sum_right_vertex->GetYaxis()->SetLabelFont(42);
  hDLmass_sum_right_vertex->GetYaxis()->SetLabelSize(0.05);
  hDLmass_sum_right_vertex->GetYaxis()->SetTitleOffset(1.35);
  hDLmass_sum_right_vertex->GetYaxis()->SetTitleSize(0.05);
  hDLmass_sum_right_vertex->GetYaxis()->SetTitleFont(42);

  
  hDLmass_sum_right_vertex->Draw("same");
  hDLmass_sum_right_vertex->SetLineColor(kMagenta);
  hDLmass_sum_right_vertex->SetLineStyle(1);
  hDLmass_sum_right_vertex->SetLineWidth(3);

  hDLmass_background_renolmalize->Draw("psame");
  //hDLmass_background_renolmalize->SetLineWidth(1);
  //hDLmass_background_renolmalize->SetLineColor(kRed);
  //hDLmass_background_renolmalize->SetMarkerColor(kRed);
  //hDLmass_background_renolmalize->SetMarkerStyle(8);
  //hDLmass_background_renolmalize->SetMarkerSize(2.4);


  hDLmass_sum_all_signals->Draw("same");
  hDLmass_sum_all_signals->SetLineColor(kGreen);
  hDLmass_sum_all_signals->SetLineWidth(4);
  hDLmass_sum_all_signals->SetFillStyle(3144);
  //hDLmass_sum_background_re->Draw("same");
  //hDLmass_sum_background_re->SetLineColor(kBlue);
  hDLmass_FF->SetLineWidth(0);
  hDLmass_FF->Draw("same");
  hDLmass_FF->SetLineStyle(3);
  hDLmass_FF->SetLineColor(kGreen);
  hDLmass_delta->Draw("same");
  hDLmass_delta->SetLineWidth(3);
  
  //draw signal channels
  int jmax=sizeof(signal_ch)/sizeof(signal_ch[0]);
  for(int j=0;j<jmax;j++)
    {
      hDLmassDistZL[signal_ch[j]]->Draw("same");
      hDLmassDistZL[signal_ch[j]]->SetLineWidth(3);
    }
  format_l1520(hDLmassDistZL[signal_ch[0]]);
  format_s1385(hDLmassDistZL[signal_ch[1]]);
  format_l1405(hDLmassDistZL[signal_ch[2]]);
  format_cb(hDLmass_background_renolmalize);

  cFinalDL_cb->Update();
  //MyFile->Close();
  //if (!( MyFile->IsOpen()) )
  //printf("File closed successfully\n");
  
  double x_min=1500;
  double x_max=1540;

  double signal=hL1520massFinalRLpi0_L[1]->Integral(hL1520massFinalRLpi0_L[1]->FindBin(x_min),hL1520massFinalRLpi0_L[1]->FindBin(x_max));
  double background=hL1520mass_background->Integral(hL1520mass_background->FindBin(x_min),hL1520mass_background->FindBin(x_max));
  
  cout<<"*********** Final raport ***********"<<endl;
  cout<<"Signal window: "<<x_min<<" - "<<x_max<<"MeV"<<endl;
  cout<<"signal: "<<signal<<" background: "<<background<<endl;
  cout<<"S/B ratio :"<< signal/background << endl;
  cout<<"significance (S/Sqrt(S+B))"<<signal/TMath::Sqrt(signal+background)<<endl;
    
  //end of printing results***************************************************************************
    
  return 1;
}
