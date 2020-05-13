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

Int_t getMotherIndex(HGeantKine* particle)
{
  Int_t trackID=particle->getTrack();
  HGeantKine* particleParent=particle->getParent(trackID);
  Int_t parentID=0;
  if(particleParent!=0)
    parentID=particleParent->getID();
  return parentID;
      
  //return particle->getGeneratorInfo1();
}



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
  const int nbin=100;
  const int maxE=1000;
  TH1F* hk_particleid=new TH1F("hk_particleid","ParticleID",100,0,50);
  TH1F* hk_particleid_primary=new TH1F("hk_particleid_primary","ParticleID for particles from primary vertex",100,0,50);
  TH1F* hk_particleid_L1116=new TH1F("hk_particleid_L1116","ParticleID for particles from L1116",100,0,50);
  TH1F* hk_minv_epem_all_pv=new TH1F("hk_minv_epem_all_pv","M^{inv}_{e^{+} e^{-}} for all leptons from primary reaction",nbin,0,maxE);
  TH1F* hk_minv_epem_all_pv_acc=new TH1F("hk_minv_epem_all_pv_acc","M^{inv}_{e^{+} e^{-}} for leptons from primary reaction in acceptance",nbin,0,maxE);
  
  TH1F* hpk_particleid=new TH1F("hpk_particleid","ParticleID for hParticleCand",100,0,50);
  TH1F* hpk_mult=new TH1F("hpk_mult","Multiplicity of tracks in HADES, only creationmechanizm==0",30,0,10);
  TH1F* hpk_particleid_primary=new TH1F("hpk_particleid_primary","ParticleID for particles from primary vertex for hParticleCand",100,0,50);
  TH1F* hpk_minv_epem_all_pv=new TH1F("hpk_minv_epem_all_pv","M^{inv}_{e^{+} e^{-}} for all leptons from primary reaction for hParticleCand",nbin,0,maxE);
  TH1F* hpk_minv_epem_all_pv_oa=new TH1F("hpk_minv_epem_all_pv_oa","M^{inv}_{e^{+} e^{-}} for all leptons from primary reaction for hParticleCand",nbin,0,maxE);
  TH1F* hpk_eff_epem_all_pv_oa=new TH1F("hpk_eff_epem_all_pv_oa","Efficiency for all leptons from primary reaction for hParticleCand",nbin,0,maxE);
  TH1F* hpk_minv_epem_all_pv_oa_pi=new TH1F("hpk_minv_epem_all_pv_oa_pi","M^{inv}_{e^{+} e^{-}} for leptons from primary reaction and #pi in HADES",nbin,0,maxE);
  TH1F* hpk_eff_epem_all_pv_oa_pi=new TH1F("hpk_eff_epem_all_pv_oa_pi","Efficiency for leptons from primary reaction and #pi in HADES",nbin,0,maxE);
  TH1F* hpk_minv_epem_all_pv_oa_pi_p=new TH1F("hpk_minv_epem_all_pv_oa_pi_p","M^{inv}_{e^{+} e^{-}} for leptons from primary reaction and #pi in HADES and p in HADES or FwDet",nbin,0,maxE);
  TH1F* hpk_eff_epem_all_pv_oa_pi_p=new TH1F("hpk_eff_epem_all_pv_oa_pi_p","Efficiency for leptons from primary reaction and #pi in HADES and p in HADES or FwDet",nbin,0,maxE);
  TH1F* hpk_minv_epem_all_pv_oa_pi_p_cut=new TH1F("hpk_minv_epem_all_pv_oa_pi_p_cut","M^{inv}_{e^{+} e^{-}} for leptons from primary reaction and #pi in HADES and p in HADES or FwDet",nbin,0,maxE);
  TH1F *hpk_ring_match_quolity=new TH1F("hpk_ring_match_quolity","Maching quality",100,0,20);  
  TH1F *hpk_opening_angle=new TH1F("hpk_opening_angle","Openieng angle for e^{+} e^{-} pair",180,0,90);

  TH2F *hpkep_inAcceptance=new TH2F("hpkep_inAcceptance","e^{+} in HADES acceptance;p[MeV];#theta",100,0,1500,100,0,100);
  TH2F *hpkem_inAcceptance=new TH2F("hpkem_inAcceptance","e^{-} in HADES acceptance;p[MeV];#theta",100,0,1500,100,0,100);


  
  const int zmin=-60;
  const int zmax=800;
  const int zn=200;
  TH1F *hk_pim_vertex=new TH1F("hk_pim_vertex","Z-coordinate for #pi^{-} vertex, from hGeantKine",zn,zmin,zmax);
  TH1F *hpk_pim_vertex=new TH1F("hpk_pim_vertex","Z-coordinate for #pi^{-} vertex, from hParticleCand",zn,zmin,zmax);
  TH1F *hpk_pim_vertex_eff=new TH1F("hpk_pim_vertex_eff","#pi^{-} efficiency in function of Z vertex coordinate",zn,zmin,zmax);

  TH1F *hk_pim_vertex_scan=new TH1F("hk_pim_vertex_scan","Z-coordinate for #pi^{-} vertex, from hGeantKine",zn,zmin,zmax);
  TH1F *hpk_pim_vertex_scan=new TH1F("hpk_pim_vertex_scan","Z-coordinate for #pi^{-} vertex, from hParticleCand",zn,zmin,zmax);
  TH1F *hpk_pim_vertex_eff_scan=new TH1F("hpk_pim_vertex_eff_scan","#pi^{-} efficiency in function of Z vertex coordinate",zn,zmin,zmax);
  //main loop
  for (Int_t i = 0; i < entries; i++)                   
    {
      loop->nextEvent(i);         // get next event. categories will be cleared before
      if(i%5000==0)
	cout<<"event no. "<<i<<endl;
      //create containers to load data
      HParticleCandSim* particlecand =nullptr;
      HParticleCandSim* particlecand2 =nullptr;
      HParticleCandSim* particlecand_pi =nullptr;
      HParticleCandSim* particlecand_p =nullptr;
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

	      if(getMotherIndex(kine)==18)
		hk_particleid_L1116->Fill(kine->getID());
	      if(kine->getMechanism() ==0)
		hk_particleid_primary->Fill(kine->getID());
	      if(kine->getID()==9 && getMotherIndex(kine)==18)//pi- from L1116
		{
		  Float_t vx;
		  Float_t vy;
		  Float_t vz;
		  kine->getVertex(vx,vy,vz);
		  hk_pim_vertex->Fill(vz);
		  for(int n=0;n<=zn;n++)
		    {
		      if(vz<(double)zmin+(double)n*((double)zmax-(double)zmin)/(double)zn)
			hk_pim_vertex_scan->Fill((double)zmin+(double)n*((double)zmax-(double)zmin)/(double)zn-0.0001);
		    }
		}
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
	  int hpart_mult=0;
	  for (int j = 0; j < hpartn; ++j)//lepton 1
	    {
	      particlecand = HCategoryManager::getObject(particlecand, fCatParticleCandSim, j);
	      if(!particlecand->isFlagBit(kIsUsed))
		continue;
	      
	      hpk_particleid->Fill(particlecand->getGeantPID());
	      if(particlecand->getGeantCreationMechanism() ==0)
		{
		  hpk_particleid_primary->Fill(particlecand->getGeantPID());
		  hpart_mult++;
		}
	      if(particlecand->getGeantCreationMechanism() ==0 && particlecand->getGeantPID()==2)//e+
		hpkep_inAcceptance->Fill(particlecand->getGeantTotalMom(),particlecand->getTheta());
	      if(particlecand->getGeantCreationMechanism() ==0 && particlecand->getGeantPID()==3)//e-
		hpkem_inAcceptance->Fill(particlecand->getGeantTotalMom(),particlecand->getTheta());
	      if(particlecand->getGeantPID()==9 && particlecand->getGeantParentPID()==18)//pi- from Lambda
		{
		  hpk_pim_vertex->Fill(particlecand->getGeantzVertex());

		  for(int n=0;n<=zn;n++)
		    {
		      if(particlecand->getGeantzVertex()<(double)zmin+(double)n*((double)zmax-(double)zmin)/(double)zn)
			hpk_pim_vertex_scan->Fill((double)zmin+(double)n*((double)zmax-(double)zmin)/(double)zn-0.0001);
		    }
		}
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
			{
			  hpk_minv_epem_all_pv_oa->Fill(ldi.M());

			  for(int pi=0;pi < hpartn; ++pi)//pion in hades
			    {
			      particlecand_pi = HCategoryManager::getObject(particlecand_pi, fCatParticleCandSim, pi);

			      if(!(particlecand_pi->isFlagBit(kIsUsed)
				   && particlecand_pi->getGeantPID()==9
				   && particlecand_pi->getGeantParentPID()==18))
				continue;
			      
			      hpk_minv_epem_all_pv_oa_pi->Fill(ldi.M());
			      for(int p=0;p < hpartn; ++p)//proton in hades
				{
				  particlecand_p = HCategoryManager::getObject(particlecand_p, fCatParticleCandSim, p);
				  if(!(particlecand_p->isFlagBit(kIsUsed)
				       && particlecand_p->getGeantPID()==14
				       && particlecand_p->getGeantParentPID()==18))
				    continue;
				  hpk_minv_epem_all_pv_oa_pi_p->Fill(ldi.M());
				}//end of p in hades

			      if (fFwDetCandSim)//FwDet class
				{
				  vcnt = fFwDetCandSim->getEntries();//number of candidates in CatVectorCandSim
				  for (int f = 0; f < vcnt; ++f)//proton in fwdet
				    {
				      fwdetstrawvec = HCategoryManager::getObject(fwdetstrawvec, fFwDetCandSim, f);
				      if(!(//fwdetstrawvec->isFlagBit(kIsUsed)
					    fwdetstrawvec->getGeantPID()==14
					   && fwdetstrawvec->getGeantParentPID()==18))
					continue;
				      hpk_minv_epem_all_pv_oa_pi_p->Fill(ldi.M());
				    }//end of proton in fwDet
				}//end of Fwdet class
			    }//end of pi
			}// end of OA
		    }//end of second lepton
		}//end of first lepton
	      hpk_mult->Fill(hpart_mult);
	    }
	}
        
     	
    } // end eventloop
  //***********************************************************************************

  //drawing histograms
  hk_minv_epem_all_pv->Draw();
  hk_minv_epem_all_pv_acc->Draw("same");

   //calculate efficiency histograms

  double scaleup=400/348;
  hpk_eff_epem_all_pv_oa->Divide(hpk_minv_epem_all_pv_oa,hk_minv_epem_all_pv,1,scaleup);
  hpk_eff_epem_all_pv_oa_pi->Divide(hpk_minv_epem_all_pv_oa_pi,hk_minv_epem_all_pv,1,scaleup);
  hpk_eff_epem_all_pv_oa_pi_p->Divide(hpk_minv_epem_all_pv_oa_pi_p,hk_minv_epem_all_pv,1,scaleup);

  hpk_pim_vertex_eff->Divide(hpk_pim_vertex,hk_pim_vertex,100.0/64.0,scaleup);
  hpk_pim_vertex_eff_scan->Divide(hpk_pim_vertex_scan,hk_pim_vertex_scan,100.0/64.0,scaleup);

  TCanvas* cZeff=new TCanvas("cZeff","cZeff");
  cZeff->Divide(2);
  cZeff->cd(1);
  hk_pim_vertex->Draw();
  hpk_pim_vertex->Draw("same");
  cZeff->cd(2);
  hpk_pim_vertex_eff->Draw();

  cZeff->Write();

  TCanvas* cZeff_scan=new TCanvas("cZeff_scan","cZeff_scan");
  cZeff_scan->Divide(2);
  cZeff_scan->cd(1);
  hk_pim_vertex_scan->Draw();
  hpk_pim_vertex_scan->Draw("same");
  cZeff_scan->cd(2);
  hpk_pim_vertex_eff_scan->Draw();

  cZeff_scan->Write();
  
  TCanvas* cLeptonAcceptance=new TCanvas("cLeptonAcceptance","Leptons' acceptance");
  cLeptonAcceptance->Divide(2);
  cLeptonAcceptance->cd(1);
  hpkep_inAcceptance->Draw("colz"); 
  cLeptonAcceptance->cd(2);
  hpkem_inAcceptance->Draw("colz");

  cLeptonAcceptance->Write();
 
  //save histograms
  output_file->Write();
  //output_file->Close();
  cout << "writing root tree done" << endl;

  timer.Stop();
  timer.Print();

  return 0;
}
