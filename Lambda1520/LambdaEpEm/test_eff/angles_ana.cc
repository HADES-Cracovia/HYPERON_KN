#include "angles_ana.h"

#include "hgeantfwdet.h"
#include "fwdetdef.h"
#include "hfwdetstrawcalsim.h"
#include "hfwdetcand.h"
#include "hfwdetcandsim.h"
#include "hparticlecandsim.h"
#include "hparticlecand.h"
#include "hparticletool.h"

#include <TCanvas.h>
#include <TStyle.h>
#include <sstream>

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";
#define PI TMath::Pi()

using namespace std;


HGeomVector calcPrimVertex_Track_Mother(const std::vector<HParticleCand *>cands, const HGeomVector & DecayVertex, const HGeomVector & dirMother, int trackA_num, int trackB_num);
HGeomVector calcPrimVertex_Track_Mother(const std::vector<HParticleCand *>cands, const HGeomVector & beamVector, const HGeomVector & DecayVertex, const HGeomVector & dirMother, int trackA_num, int trackB_num);

Int_t fwdet_tests(HLoop * loop, const AnaParameters & anapars)
{
  if (!loop->setInput(""))
    {                                                    // reading file structure
      std::cerr << "READBACK: ERROR : cannot read input !" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  TStopwatch timer;
  timer.Reset();
  timer.Start();

  //////////////////////////////////////////////////////////////////////////////
  //      Fast tree builder for creating of ntuples                            //
  //////////////////////////////////////////////////////////////////////////////

  loop->printCategories();    // print all categories found in input + status

  //create categorys - a containers designed to store data. In TBrowser you see all of them, all important for you are hire
  HCategory * fCatGeantKine = nullptr;
  fCatGeantKine = HCategoryManager::getCategory(catGeantKine, kTRUE, "catGeantKine");
  if (!fCatGeantKine)
    {
      cout << "No catGeantKine!" << endl;
      exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

  HCategory * fFwDetStrawCal = nullptr;
  fFwDetStrawCal = HCategoryManager::getCategory(catFwDetStrawCal, kTRUE, "catFwDetStrawCalSim");
  if (!fFwDetStrawCal)
    {
      cout << "No catFwDetStrawCal!" << endl;
      exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

  HCategory * fFwDetCandSim = nullptr;
  fFwDetCandSim = HCategoryManager::getCategory(catFwDetCand, kTRUE, "catFwDetCand");
  if (!fFwDetCandSim)
    {
      cout << "No catFwDetCand!" << endl;
      //exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }


  HCategory * fCatParticleCandSim= nullptr;
  fCatParticleCandSim = HCategoryManager::getCategory(catParticleCand, kTRUE, "catParticleCand");
  if(!fCatParticleCandSim)
    {
      cout<< "No catParticleCandSim!"<<endl;
    }

  Int_t entries = loop->getEntries();
  //setting numbers of events regarding the input number of events by the user
  if (anapars.events < entries and anapars.events >= 0 ) entries = anapars.events;
  
  // specify output file
  TFile * output_file = TFile::Open(anapars.outfile, "RECREATE");
  output_file->cd();
  //
  cout << "NEW ROOT TREE " << endl;
  //
  //crete histograms
  TH1F* hk_particleid=new TH1F("hk_particleid","ParticleID",100,0,50);
  TH1F* hk_particleid_primary=new TH1F("hk_particleid_primary","ParticleID for particles from primary vertex",100,0,50);
  TH1F* hk_minv_epem_all_pv=new TH1F("hk_minv_epem_all_pv","M^{inv}_{e^{+} e^{-}} for all leptons from primary reaction",1000,0,1000);
  TH1F* hk_minv_epem_all_pv_acc=new TH1F("hk_minv_epem_all_pv_acc","M^{inv}_{e^{+} e^{-}} for leptons from primary reaction in acceptance",1000,0,1000);
  
TH1F* hpk_particleid=new TH1F("hpk_particleid","ParticleID for hParticleCand",100,0,50);
  TH1F* hpk_particleid_primary=new TH1F("hpk_particleid_primary","ParticleID for particles from primary vertex for hParticleCand",100,0,50);
  TH1F* hpk_minv_epem_all_pv=new TH1F("hpk_minv_epem_all_pv","M^{inv}_{e^{+} e^{-}} for all leptons from primary reaction for hParticleCand",1000,0,1000);
  TH1F* hpk_minv_epem_all_pv_oa=new TH1F("hpk_minv_epem_all_pv_oa","M^{inv}_{e^{+} e^{-}} for all leptons from primary reaction for hParticleCand",1000,0,1000);
  TH1F *hpk_ring_match_quolity=new TH1F("hpk_ring_match_quolity","Maching quality",100,0,20);  
  TH1F *hpk_opening_angle=new TH1F("hpk_opening_angle","Openieng angle for e^{+} e^{-} pair",180,0,90);  
  //main loop
  for (Int_t i = 0; i < entries; i++)                   
    {
      loop->nextEvent(i);         // get next event. categories will be cleared before
      if(i%5000==0)
	cout<<"event no. "<<i<<endl;
      //create containers to load data
      HParticleCandSim* particlecand =nullptr;
      HParticleCandSim* particlecand2 =nullptr;
      HFwDetCandSim* fwdetstrawvec = nullptr;
      HGeantKine* kine=nullptr;
      HGeantKine* kine2=nullptr;
      HFwDetStrawCalSim* strawcal=nullptr;
      int vcnt, gknt, hpartn;
      double phisim,phisim0,phisim1, phirec,thetasim, thetasim1,thetasim2,thetasim0, thetarec, xsim, ysim, xrec, yrec, zsts1, rsim,rsim1,rsim2;
      HParticleTool tool;
      //FW Detector iterate over fFwDetCandSim category
      /*if (fFwDetCandSim)
	{
	  vcnt = fFwDetCandSim->getEntries();//number of candidates in CatVectorCandSim
	  gknt = fCatGeantKine->getEntries();//numbers of candidates in fCatGeantKine
	  for (int j = 0; j < vcnt; ++j)
	    {
	      fwdetstrawvec = HCategoryManager::getObject(fwdetstrawvec, fFwDetCandSim, j);
	      	     
	    }
	}
      */  
      //Kine iterate over fCatGeantKine category
      if(fCatGeantKine)
	{
	  gknt=fCatGeantKine->getEntries();
	  for(int i=0;i<gknt;i++)//lepton 1
	    {
	      kine = HCategoryManager::getObject(kine, fCatGeantKine,i);
	      hk_particleid->Fill(kine->getID());

	      if(kine->getMechanism() ==0)
		hk_particleid_primary->Fill(kine->getID());

	      for(int j=i;j<gknt;j++)//lepton 2
		{
		  
		  kine2 = HCategoryManager::getObject(kine2, fCatGeantKine,j);

		  if(((kine->getID()==2 && kine2->getID()==3)||(kine->getID()==3 && kine2->getID()==2))
		     && kine->getMechanism() ==0
		     && kine2->getMechanism()==0
		     )//e+e- from PLUTO
		    {
		      TLorentzVector l1;
		      TLorentzVector l2;
		      TLorentzVector ldi;
		      HGeomVector p1;
		      HGeomVector p2;
		      kine->getMomentum(p1);
		      kine2->getMomentum(p2);
		      l1.SetPxPyPzE(p1.getX(),p1.getY(),p1.getZ(),kine->getE());
		      l2.SetPxPyPzE(p2.getX(),p2.getY(),p2.getZ(),kine2->getE());
		      ldi=l1+l2;

		      hk_minv_epem_all_pv->Fill(ldi.M());

		      if(kine->isInAcceptance(4,4,4,4,1,0) && kine2->isInAcceptance(4,4,4,4,1,0)) //tutaj mam problem bo na razie wsyztko zwraca false
			{
			  hk_minv_epem_all_pv_acc->Fill(ldi.M());
			}
		    }
		}
	    }
	}
      //iterate ofer ParticleCandSim
      if(fCatParticleCandSim)
	{
	  hpartn = fCatParticleCandSim->getEntries();//number of candidates in CatVectorCandSim
	  for (int j = 0; j < hpartn; ++j)//lepton 1
	    {
	      particlecand = HCategoryManager::getObject(particlecand, fCatParticleCandSim, j);
	      if(!particlecand->isFlagBit(kIsUsed))
		continue;
	      
	      hpk_particleid->Fill(particlecand->getGeantPID());
	      if(particlecand->getGeantCreationMechanism() ==0)
		hpk_particleid_primary->Fill(particlecand->getGeantPID());

	      for (int i = j; i < hpartn; ++i)//lepton 2
		{
		  particlecand2 = HCategoryManager::getObject(particlecand2, fCatParticleCandSim, i);
		  if(!particlecand2->isFlagBit(kIsUsed))
		    continue;
	      
		  if(((particlecand->getGeantPID()==2 && particlecand2->getGeantPID()==3)
		      ||(particlecand->getGeantPID()==3 && particlecand2->getGeantPID()==2))
		     && particlecand->getGeantCreationMechanism() ==0
		     && particlecand2->getGeantCreationMechanism()==0
		     )//e+e- from PLUTO
		    {
		      TLorentzVector l1;
		      TLorentzVector l2;
		      TLorentzVector ldi;
		      double oa_lepton;

		      l1.SetPxPyPzE(particlecand->getGeantxMom(),
				    particlecand->getGeantyMom(),
				    particlecand->getGeantzMom(),
				    TMath::Sqrt(TMath::Power(particlecand->getGeantTotalMom(),2)
						+TMath::Power(particlecand->getGeantGenweight(),2))
				    );
		      l2.SetPxPyPzE(particlecand2->getGeantxMom(),
				    particlecand2->getGeantyMom(),
				    particlecand2->getGeantzMom(),
				    TMath::Sqrt(TMath::Power(particlecand2->getGeantTotalMom(),2)
						+TMath::Power(particlecand2->getGeantGenweight(),2))
				    );
		     
		      ldi=l1+l2;
		      oa_lepton=tool.getOpeningAngle(particlecand,particlecand2);

		      	
		      hpk_ring_match_quolity->Fill(particlecand->getRichMatchingQuality());
		      hpk_ring_match_quolity->Fill(particlecand2->getRichMatchingQuality());
		      hpk_minv_epem_all_pv->Fill(ldi.M());
		      hpk_opening_angle->Fill(oa_lepton);

		      if(oa_lepton>4)
			hpk_minv_epem_all_pv_oa->Fill(ldi.M());
		    }
		}
	    }
	}
        
     	
    } // end eventloop
  //***********************************************************************************

  //drawing histograms
  hk_minv_epem_all_pv->Draw();
  hk_minv_epem_all_pv_acc->Draw("same");

  
  //save histograms
  output_file->Write();
  //output_file->Close();
  cout << "writing root tree done" << endl;

  timer.Stop();
  timer.Print();

  return 0;
}
