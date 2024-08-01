

#include "ecal_tests.h"
#include "emcdef.h"

#include <hades.h>
#include <hdst.h>
#include <hparasciifileio.h>
#include <hspectrometer.h>
#include <hparticlecand.h>
#include "hemcneutralcand.h"
#include "hemcneutralcandsim.h"
#include "hemcclustersim.h"
#include "hemccluster.h"
#include "hparticlecandsim.h"
#include "hgeantkine.h"

#include <hgeomvolume.h>
#include <hgeomcompositevolume.h>
#include <hgeomvector.h>
#include "hparticleevtinfo.h"
#include "hparticletracksorter.h"
#include "hphysicsconstants.h"
#include "htool.h"
#include "htime.h"

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TH2I.h>
#include <TH3I.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TMath.h>

#include <algorithm>
#include <cstdlib>
#include <vector>

#define PI 3.14159265

#define PR(x)                                                                                      \
    std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

using namespace std;

const double D2R = 1.74532925199432955e-02;
const double R2D = 57.2957795130823229;





//Int_t ecal_tests(HLoop* loop, const AnaParameters& anapars)
Int_t ecal_tests(TString inputlist, TString outfile, TString eventsList, Int_t nev)
{
      HLoop* loop = new HLoop(kTRUE);


      
      //Long64_t read = 500000000;
      //TFile::SetReadaheadSize(read);
      //loop->setTreeCacheSize(read);

      
      if (inputlist.EndsWith(".list")) { loop->addFilesList(inputlist); }
    else
    {
        loop->addMultFiles(inputlist);
    }

   // Check if loop was properly initialized
    if (!loop->setInput("-*,+HParticleCand,+HEmcNeutralCand,+HEmcCluster"))

      { // reading file structure
        std::cerr << "READBACK: ERROR : cannot read input !" << std::endl;
        std::exit(EXIT_FAILURE);
    }


    gHades->setBeamTimeID(HADES::kFeb22);
    HBeamTime::printBeamTimeInfo();

 


      Int_t entries = loop->getEntries(); // Number of entries in loop




    // Timer for checking analysis time
    TStopwatch timer;
    timer.Reset(); // Reset timer
    timer.Start(); // Start timer (T0)


    // Hades Particle Candidates
    HCategory* fParticleCand = HCategoryManager::getCategory(catParticleCand, kTRUE, "catParticleCandSim");
    if (!fParticleCand) { cout << "No catParticleCand!" << endl; }

    HCategory* fEmcNeutralCand = HCategoryManager::getCategory(catEmcNeutralCand, kTRUE, "catEmcNeutralCandSim");
    if(!fEmcNeutralCand){cout << "No catEmcNeutralCand!" << endl;}
    
    HCategory* fEmcCluster = HCategoryManager::getCategory(catEmcCluster, 0, "catEmcCluster");
    if(!fEmcCluster){cout << "No catEmcCluster!" << endl;}

    //HCategory * fStart2Hit = HCategoryManager::getCategory(catStart2Hit, kTRUE, "catStart2Hit");
    //if (!fStart2Hit) { cout << "No catStart2Hit!" << endl; }


    HEnergyLossCorrPar dEdxCorr;
    dEdxCorr.setDefaultPar("feb22");    
    //----------------------------------------------------------------------------------------------
    // Setting parameters for loop over events
    //----------------------------------------------------------------------------------------------


    //----------------------------------------------------------------------------------------------
    // Specifying output file
    //----------------------------------------------------------------------------------------------

    TFile* output_file = new TFile(outfile.Data(), "RECREATE");
    output_file->cd();
    cout << "NEW ROOT TREE " << endl;

    //----------------------------------------------------------------------------------------------
  
    TH2F *hbeta_mom= new TH2F("hbeta_mom","hbeta_mom",1000,-4000.,4000.,700,0.,1.4);  ///poprawki na straty energii
    TH2F *hbeta_mom_p= new TH2F("hbeta_mom_p","hbeta_mom_p",1000,0.,4000.,700,0.,1.4);
    TH2F *hbeta_mom_pip= new TH2F("hbeta_mom_pip","hbeta_mom_pip",1000,0.,4000.,700,0.,1.4);
    TH2F *hbeta_mom_pim= new TH2F("hbeta_mom_pim","hbeta_mom_pim",1000,-4000.,0.,700,0.,1.4);

    TH2F* hmass_mom = new TH2F("hmass_mom","hmass_mom; p*q [MeV/c]; mass",1000,-2000,4000,1000,0,2000);
    TH1F *hmass=new TH1F("hmass","hmass",1000,-2000,2000);
    
    TH1F *hbeta=new TH1F("hbeta","hbeta",100,0,2);  //dodac katy i pedy czastek 
    TH1F *hg_energy=new TH1F("hg_energy","hg_energy",2000,0,2000);
	TH1F *hMM= new TH1F("missing_mass_lambda", "missing mass lambda [MeV]",1000,0.,4000.); //z prawej wzgorze od pi-2
	TH1F *hIM= new TH1F("inv_mass_kaon+prot", "inv mass kaon + proton [MeV]",1000,200.,600.); //ogon z lewej od pi-2:w
    TH1F *hMM_pim2 = new TH1F("missing_mass_pim2", "missing mass kaon + proton daje pim2 [MeV]",1000,0.,600.);
	TH1F *hIM_kaon = new TH1F("inv_mass_kaon", "inv mass kaonu [MeV]",1000,200,600.);
    //*********************************************************
    //****************************************
    //PROT PID - mass vs mom

    TCutG *cutP = new TCutG("cutP",11); 
    cutP->SetPoint(0,188.0792,1187.5);
    cutP->SetPoint(1,1559.183,1378.472);
    cutP->SetPoint(2,2203.665,1517.361);
    cutP->SetPoint(3,2677.548,1579.861);
    cutP->SetPoint(4,2886.057,1392.361);
    cutP->SetPoint(5,2544.861,1177.083);
    cutP->SetPoint(6,2052.022,888.8889);
    cutP->SetPoint(7,188.0792,593.75);
    cutP->SetPoint(8,86.98397,642.3611);
    cutP->SetPoint(9,200.7161,1190.972);
    cutP->SetPoint(10,188.0792,1187.5);

    //****************************************

    //*******************************************************
    const double mp=938.27231;
    const double mpim=140.;
    
    vector<TLorentzVector> lv_neutr, lv_neutr1, lv_prot, lv_pip, lv_pim ;

    
    const double beam_energy = 1170; //sqrt(p^2 +m^2)
    const double beam_momentum = sqrt(beam_energy*beam_energy-mp*mp);

    TLorentzVector proj(0,0,beam_momentum, beam_energy); //porusza sie wzdloz z
    TLorentzVector targ(0,0,0, mp); // x,y,z =0 spoczynek
    TLorentzVector beam(0,0,0,0); //stan poczatkowy
    beam = proj + targ;
    
    string currentFileName;


      for (int i = 0; i < entries; i++) // event loop
	{
        
        if (i % 10000 == 0)
        {
            printf("Event nr.: %d\n", i);
        }

        //Int_t nbytes =
	  loop->nextEvent(i); // get next event. categories will be cleared before

	HTool::printProgress(i, nev, 1, "Analysis :");

	/*
       TString tempFileName;
        Int_t dayOfYear;
        Int_t hour;
        Int_t min;
        Int_t second;
        if (loop->isNewFile(tempFileName))
        {
            TString type;
            Int_t year;
            Int_t eb;
            currentFileName = tempFileName;
            HTime::splitFileName(HTime::stripFileName(tempFileName), type, year, dayOfYear, hour, min, second, eb);
        }

        Float_t valueMinute = dayOfYear + hour / 24. + min / 24. / 60.;
        Float_t valueSecond = valueMinute + second / 24. / 60. / 60.;
	*/
	
	lv_neutr.clear();
	lv_prot.clear();
	lv_pip.clear();
	lv_pim.clear();

	//HStart2Hit * fstart = nullptr;
	//fstart = (HStart2Hit *) fStart2Hit->getObject(0);
	//if (!fstart || fstart->getCorrFlag() == -1) continue;
	//fCorrFlag=-1 iTOF
	//fCorrFlag=0 only LGAD
	//fCorrFlag>=0
	
	
        HEventHeader* event_header = NULL;
        if (!(event_header = gHades->getCurrentEvent()->getHeader())) continue;

        Int_t TBit = (Int_t)event_header->getTBit();
        Double_t VertexX = event_header->getVertexReco().getX();
        Double_t VertexY = event_header->getVertexReco().getY();
        Double_t VertexZ = event_header->getVertexReco().getZ();   //czy vertex pokrywa sie nam z rozkladem tarczy

	if(VertexZ<-200 || VertexZ>0) continue;  //jakies defaultowe ciecie, najlepiej przed continue narysowac sobie rozklad ,,z"
	//if(VertexZ<-148. || VertexZ>-97.) continue;



	//         Int_t DF     = (Int_t) event_header->getDownscalingFlag();
        //         Int_t SeqNum = (Int_t) event_header->getEventSeqNumber();
        //         Int_t TDec   = (Int_t) event_header->getTriggerDecision();

	//****************************************************
       
	Int_t nNeutral_ev = fEmcNeutralCand->getEntries();

	
	
	for (int j = 0; j < nNeutral_ev; ++j)
	  {
	  
	    HEmcNeutralCandSim* neutr_cand = HCategoryManager::getObject(neutr_cand, fEmcNeutralCand, j);

		//Float_t dist  = neutr_cand->getDistanceToEmc();
		Int_t ind=neutr_cand->getEmcClusterIndex();


		HEmcCluster *cl=nullptr;
		cl=HCategoryManager::getObject(cl, fEmcCluster, ind);
		Int_t cl_size = cl->getNCells();

		Int_t sec = cl->getSector();
		Int_t cel = cl->getCell();		
		Int_t pid=neutr_cand->getPID();

		if(cel<33) continue;

		Double_t energy  = cl->getEnergy();
		Double_t tof  =cl->getTime();
		Double_t beta = neutr_cand->getBeta();
		Float_t theta = cl->getTheta();
		Float_t phi = cl->getPhi();

		
		  
		hbeta->Fill(beta);
		hg_energy->Fill(energy);

		TLorentzVector lvg1, lvg;
		//lvg1.SetXYZM(trackVec.getX(),trackVec.getY(),trackVec.getZ(),0);
	       
		lvg1 = *neutr_cand;
		if (energy>100 && pid==1){

		  lv_neutr.push_back(lvg1);	    

		}



	  }
	//*********************************************************

	if (fParticleCand)
	      {
		Int_t nPart_ev  = fParticleCand->getEntries();
		//cout<<"---------->>> "<<endl;
		for (int j = 0; j < nPart_ev; ++j)  //NAJWAZNIEJSZA CZESC
		  {
		    HParticleCandSim* fparticlecand = HCategoryManager::getObject(fparticlecand, fParticleCand, j);
		   
		    Float_t theta = fparticlecand->getTheta(); //wszystko po rekonstrukcji
		    Float_t phi = fparticlecand->getPhi();
		    Float_t mom = fparticlecand->getMomentum();
		    Float_t beta = fparticlecand->getBeta();
		    Float_t tof = fparticlecand->getTof();
		    Int_t sec = fparticlecand->getSector();
		    //Float_t d2meta=fparticlecand->getDistanceToMetaHit();
		    Int_t system = fparticlecand->getSystem();
		    Float_t charge=fparticlecand->getCharge(); 
		    Float_t mass= fparticlecand->getMass();
		    
		    //if(!fparticlecand->isFlagBit(kIsUsed))continue;
		    //if(fparticlecand->isFlagBit(kIsLepton))continue;
		    int geant_id = fparticlecand->getGeantPID();
                    int geant_info1 = fparticlecand->getGeantGeninfo1();
                    int parent_id=fparticlecand->getGeantParentPID();


		    //if(parent_id==-1) cout<< i << " " << geant_id <<" "<< geant_info1   <<" "<<parent_id<<endl;
		    
		    if(fparticlecand->isFlagBit(kIsUsed) && fparticlecand->getGeantParentTrackNum()==0 && fparticlecand->getGeantParentPID()==-1)
		      {         ///trajektorie jakies nw                            ///te 2 warunki wybieraja pierwotne czastki, nie obchodza nas wtorne, odbite np
			//cout<<fparticlecand->getGeantPID()<<endl;
			//cout<< i << " " << geant_id <<" "<< geant_info1   <<" "<<parent_id<<endl;
			
			hbeta_mom->Fill(charge*mom,beta);
			hmass_mom->Fill(charge*mom,mass);
			hmass->Fill(charge*mass);


			//if(charge>0 &&  cutP->IsInside(charge*mom,mass)){//PID protony
			if(fparticlecand->getGeantPID()==14){//PID protony  //spr czy jest to proton

			  Double_t deltaMom = dEdxCorr.getDeltaMom(14, mom, theta); 
			  fparticlecand->setMomentum(mom + deltaMom);                                     
			  fparticlecand->calc4vectorProperties(HPhysicsConstants::mass(14));
			  hbeta_mom_p->Fill(charge*mom,beta);

  
			  //part_pos.push_back(fparticlecand);
			  TLorentzVector lv1;
			  lv1=*fparticlecand;
			  lv_prot.push_back(lv1);
			}			

			
			//if(charge>0 &&  cutPIP->IsInside(charge*mom,mass)){//PID pi+
			if(fparticlecand->getGeantPID()==8){//PID pi+

			  Double_t deltaMom = dEdxCorr.getDeltaMom(8, mom, theta); 
			  fparticlecand->setMomentum(mom + deltaMom);                                     
			  fparticlecand->calc4vectorProperties(HPhysicsConstants::mass(8));
			  hbeta_mom_pip->Fill(charge*mom,beta);

			  
			  //cout<<charge*mom<<" "<<beta<<endl;
			  TLorentzVector lv1a;
			  lv1a=*fparticlecand;
			  lv_pip.push_back(lv1a);
			  
			}			

			//if(charge<0 && mass < 550){//PID pi-
			//if(charge<0 && cutPIM->IsInside(charge*mom,mass)){//PID pi-
			if(fparticlecand->getGeantPID()==9){//PID pi-

			  Double_t deltaMom = dEdxCorr.getDeltaMom(9, mom, theta); 
			  fparticlecand->setMomentum(mom + deltaMom);                                     
			  fparticlecand->calc4vectorProperties(HPhysicsConstants::mass(9));

			  TLorentzVector lv1b;
			  //lv1b=TLV_p(theta, phi, mom, 140.);
			  lv1b=*fparticlecand;

			  lv_pim.push_back(lv1b);
			  hbeta_mom_pim->Fill(charge*mom,beta);

			}			
			//******************************************	
		      }
		  }//fParticleCand
	      }//fParticleCand   //MOJEEEEEEEEE!!!!!!!!!!!!!
		//*****************************************************
		//***** pim+p -> K [pi+pi-] + Lambda [p pi-]
		
		  float mm_L, mm2, invM_kp, invM_pipi, mm_pipip; //mm - missing mass lambdy, invM - masa inwariantna kaonu 
		 //TLorentzVector lv3;

		  if(lv_pip.size() && lv_pim.size() && lv_prot.size()){ //czy sa, czy nie wypadly (protoniki)
		  for (int i=0;i<lv_pim.size();i++){
		  for (int j=0;j<lv_pip.size();j++){
		  
		  mm_L = (beam - lv_pip[j] - lv_pim[i]).M();//MeV 
		  mm2 = mm_L * 1e-6; //GeV missing mass lambdy :)
		  hMM->Fill(mm_L);

		  invM_pipi = (lv_pip[j] + lv_pim[i]).M(); //masa inv pim+pip (kaon)
          hIM_kaon->Fill(invM_pipi);		

		  if(mm_L>1170 && mm_L<1420)
			{ //po masie lambdy robimy ciecie od 1150 do 1450 ogon z lewej to pi2-(?)
			 
			invM_kp = (lv_pip[j] + lv_pim[i] + lv_prot[1]).M(); //kaon z protonikie		  	
			hIM->Fill(invM_kp);

				//if(invM_pipi > 480 && invM_pipi <520)
				//{
				//	mm_pipip = (beam - lv_pip[j] - lv_pim[i] - lv_prot[1]).M();
				//	hMM_pim2->Fill(mm_pipip);
				//}
			}

		  if(invM_pipi > 480 && invM_pipi < 520) 
		   {
				mm_pipip = (beam - lv_pip[j] - lv_pim[i] - lv_prot[0]).M();
				hMM_pim2->Fill(mm_pipip); //missing mass kaonu + protonu co daje mi pik o masie pi-2
		   }
 
		}
	      }
	    }
	    
	    //****************************************
	    
	    //**********************************************	    
	    }  //end event loop

	    
    output_file->cd();

    hbeta_mom->Write();
    hmass_mom->Write();
    hmass->Write();
    hbeta_mom_p->Write();
    hbeta_mom_pip->Write();
    hbeta_mom_pim->Write();
	hIM->Write();
	hMM->Write();
	hMM_pim2->Write();    
    hIM_kaon->Write();
    output_file->Close();
    cout << "writing root tree done" << endl;
    //*********************************
    //*******************************

    timer.Stop();
    timer.Print();

    return 0;
}
