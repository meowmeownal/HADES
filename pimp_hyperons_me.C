#include "TH2.h"
#include "TH2.h"
#include "TH3.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLorentzVector.h" 
#include "TMath.h"
#include "TRandom3.h"

void pimp_hyperons_me(){
  PUtils::SetSeed(0);

    makeDistributionManager()->Enable("strangeness");


  TFile *f1= new TFile("out_pimp_test.root","recreate");

///---------------	KATY	----------------------
  TH1F * hist_2_p = new TH1F ("mom_p","mom p",500,0,1500);
  TH1F * hist_2_pi_min1 = new TH1F ("mom_pi_min1","mom pi min 1",500,0,1500);
  TH1F * hist_2_pi_min2 = new TH1F ("mom_pi_min2","mom pi min 2",500,0,1500);
  TH1F * hist_2_pi_plus = new TH1F ("mom_pi_plus","mom pi plus",500,0,1500);
  TH1F * hist_2_kaon = new TH1F ("mom_kaon","mom kaon",500,0,1500);
  TH1F * hist_2_lambda = new TH1F ("mom_lambda","mom lambda",500,0,1500);

///--------------	PEDY	-----------------------------
  TH1F * hist_1_p = new TH1F ("theta_proton","th p",500,0,80);
  TH1F * hist_1_pi_min1 = new TH1F ("theta_pi_min1","pi minus th 1",500,0,180);
  TH1F * hist_1_pi_min2 = new TH1F ("theta_pi_min2","pi minus th 2",500,0,180);
  TH1F * hist_1_pi_plus = new TH1F ("theta_pi_plus","pi plus th",500,0,180);  
  TH1F * hist_1_kaon = new TH1F ("theta_kaon", "kaon th",500,0,100);
  TH1F * hist_1_lambda = new TH1F ("theta_lambda","lambda th",500,0,30);

///-------------	PEDY OD KATOW               ----------------------
  TH1F * hist_3_p = new TH1F ("mom_prot","mom_prot",500,0,4500);
  TH1F * hist_3_pi_min1 = new TH1F ("mom_th_pi_min1","mom(th) pion minus 1",500,0,1000);
  TH1F * hist_3_pi_min2 = new TH1F ("mom_th_pi_min2","mom(th) pion minus 2",500,0,1000);
  TH1F * hist_3_pi_plus = new TH1F ("mom_th_pi_plus","mom(th) pion plus",500,0,1000);
  TH1F * hist_3_kaon = new TH1F ("mom_th_kaon", "mom(th) kaon",500,0,1000);
  TH1F * hist_3_lambda = new TH1F ("mom_th_lambda","mom(th) lambda",500,0,1000);

//-------------		CM LAMBDA	---------------------------
  TH1F *hist_cms_mom_p = new TH1F("cms_mom_p", "CMS lambda mom p", 500, 0, 80);
  TH1F *hist_cms_mom_pi_min1 = new TH1F("cms_mom_pi_min1", "CMS lambda mom pi min 1", 500, 0, 180);
  TH1F *hist_cms_mom_pi_min2 = new TH1F("cms_mom_pi_min2", "CMS lambda mom pi min 2", 500, 0, 180);
  TH1F *hist_cms_mom_pi_plus = new TH1F("cms_mom_pi_plus", "CMS lambda mom pi plus", 500, 0, 100);
  TH1F *hist_cms_mom_kaon = new TH1F("cms_mom_kaon", "CMS lambda mom kaon", 500, 0, 30);

  TH1F *hist_cms_th_p = new TH1F("cms_mom_p", "CMS lambda mom p", 500, 0, 1500);
  TH1F *hist_cms_th_pi_min1 = new TH1F("cms_th_pi_min1", "CMS lambda th pi min 1", 500, 0, 1500);
  TH1F *hist_cms_th_pi_min2 = new TH1F("cms_th_pi_min2", "CMS lambda th pi min 2", 500, 0, 1500);
  TH1F *hist_cms_th_pi_plus = new TH1F("cms_th_pi_plus", "CMS lambda th pi plus", 500, 0, 1500);
  TH1F *hist_cms_th_kaon = new TH1F("cms_th_kaon", "CMS lambda th kaon", 500, 0, 1500);

 //############################################
  PReaction my_reaction2(1.16,"pi-","p","K0S [pi+ pi-] Lambda [p pi-]","pimp_test",1,0,1,1);

 

 my_reaction2.Do("thp = ([p]->Theta() * 180.)/TMath::Pi();");
 my_reaction2.Do("th_pion_min1 = ([pi-,1]->Theta() * 180.)/TMath::Pi();");
 my_reaction2.Do("th_pion_min2 = ([pi-,2]->Theta() * 180.)/TMath::Pi();");
 my_reaction2.Do("th_pion_plus = ([pi+]->Theta() * 180.)/TMath::Pi();");
 my_reaction2.Do("th_kaon = ([K0S]->Theta() * 180.)/TMath::Pi();");
 my_reaction2.Do("th_lambda = ([Lambda]->Theta() * 180.)/TMath::Pi();");

 my_reaction2.Do("mom_p = [p]->P()*1000.;");
 my_reaction2.Do("mom_pi_min1 = [pi-,1]->P()*1000.;");
 my_reaction2.Do("mom_pi_min2 = [pi-,2]->P()*1000.;");
 my_reaction2.Do("mom_pi_plus = [pi+]->P()*1000.;");
 my_reaction2.Do("mom_kaon = [K0S]->P()*1000.;");
 my_reaction2.Do("mom_lambda = [Lambda]->P()*1000.;");
 
//my_reaction2.Do("beam1=[p+p];");
 //my_reaction2.Do("beam2=P3E(0,0,5.3565,6.376);"); 
 //my_reaction2.Do("mm1=(beam1-([p,1]+[p,2]))->M()*1000;");
 //my_reaction2.Do("mm2=(beam2-([p,1]+[p,2]))->M()*1000;");


 my_reaction2.Do("invMpippim=([pi-]+[pi+])->M();");
 //my_reaction2.Do("echo $invMpippim;"); //wypluwa

///-------------------		CMS		-------------------------------

 my_reaction2.Do("boost_beta = -[Lambda]->BoostVector();");

 my_reaction2.Do("proton_cms = [p]->Boost(boost_beta);");
 my_reaction2.Do("pi_min1_cms = [pi-,1]->Boost(boost_beta);");
 my_reaction2.Do("pi_min2_cms = [pi-,2]->Boost(boost_beta);");
 my_reaction2.Do("pi_plus_cms = [pi+]->Boost(boost_beta);");
 my_reaction2.Do("kaon_cms = [K0S]->Boost(boost_beta);");

 my_reaction2.Do("mom_p_cms = proton_cms.P()*1000.;");
 my_reaction2.Do("mom_pi_min1_cms = pi_min1_cms.P()*1000.;");
 my_reaction2.Do("mom_pi_min2_cms = pi_min2_cms.P()*1000.;");
 my_reaction2.Do("mom_pi_plus_cms = pi_plus_cms.P()*1000.;");
 my_reaction2.Do("mom_kaon_cms = kaon_cms.P()*1000.;");

 my_reaction2.Do("thp_cms = (proton_cms.Theta() * 180.)/TMath::Pi();");
 my_reaction2.Do("th_pion_min1_cms = (pi_min1_cms.Theta() * 180.)/TMath::Pi();");
 my_reaction2.Do("th_pion_min2_cms = (pi_min2_cms.Theta() * 180.)/TMath::Pi();");
 my_reaction2.Do("th_pion_plus_cms = (pi_plus_cms.Theta() * 180.)/TMath::Pi();");
 my_reaction2.Do("th_kaon_cms = (kaon_cms.Theta() * 180.)/TMath::Pi();");

//--------------------------------------------------------------------------------

// ped od katow
 my_reaction2.Do(hist_3_p,"if thp>18 && thp<90;_x=mom_p;");
 my_reaction2.Do(hist_3_pi_min1,"if th_pion_min1>18 && th_pion_min1<90;_x=mom_pi_min1;");
 my_reaction2.Do(hist_3_pi_min2,"if th_pion_min2>18 && th_pion_min2<90;_x=mom_pi_min2;");
 my_reaction2.Do(hist_3_pi_plus,"if th_pion_plus>18 && th_pion_plus<90;_x=mom_pi_plus;");
 my_reaction2.Do(hist_3_kaon,"if th_kaon>18 && th_kaon<90;_x=mom_kaon;");
 my_reaction2.Do(hist_3_lambda,"if th_lambda>18 && th_lambda<90;_x=mom_lambda;");

// katy
 my_reaction2.Do(hist_1_pi_min1,"_x=th_pion_min1;");
 my_reaction2.Do(hist_1_pi_min2,"_x=th_pion_min2;");
 my_reaction2.Do(hist_1_pi_plus,"_x=th_pion_plus;");
 my_reaction2.Do(hist_1_kaon,"_x=th_kaon;");
 my_reaction2.Do(hist_1_lambda,"_x=th_lambda;");
 my_reaction2.Do(hist_1_p,"_x=thp;");

// pedy
 my_reaction2.Do(hist_2_pi_min1,"_x=mom_pi_min1;");
 my_reaction2.Do(hist_2_pi_min2,"_x=mom_pi_min2;");
 my_reaction2.Do(hist_2_pi_plus,"_x=mom_pi_plus;");
 my_reaction2.Do(hist_2_kaon,"_x=mom_kaon;");
 my_reaction2.Do(hist_2_lambda,"_x=mom_lambda;");
 my_reaction2.Do(hist_2_p,"_x=mom_p;");

// CMS 
 my_reaction2.Do(hist_cms_th_p, "_x = th_p_cms;");
 my_reaction2.Do(hist_cms_th_pi_min1, "_x = th_pion_min1_cms;");
 my_reaction2.Do(hist_cms_mom_pi_min2, "_x = mom_pi_min2_cms;");
 my_reaction2.Do(hist_cms_th_pi_plus, "_x = th_pion_plus_cms;");
 my_reaction2.Do(hist_cms_th_kaon, "_x = th_kaon_cms;");

 my_reaction2.Do(hist_cms_mom_p, "_x = mom_p_cms;");
 my_reaction2.Do(hist_cms_mom_pi_min1, "_x = mom_pi_min1_cms;");
 my_reaction2.Do(hist_cms_mom_pi_min2, "_x = mom_pi_min2_cms;");
 my_reaction2.Do(hist_cms_mom_pi_plus, "_x = mom_pi_plus_cms;");
 my_reaction2.Do(hist_cms_mom_kaon, "_x = mom_kaon_cms;"); 

 my_reaction2.Print();
 // my_reaction2.Loop(1000000);
 my_reaction2.Loop(10000);

 //###########################################

 //TCanvas *c1 = new TCanvas("theta","theta");
 //c1->Divide(3,1);

 //**********************************

 f1->cd();

 hist_1_p->Write();
 hist_1_pi_min1->Write();
 hist_1_pi_min2->Write();
 hist_1_pi_plus->Write();
 hist_1_kaon->Write();
 hist_1_lambda->Write();

 hist_2_p->Write();
 hist_2_pi_min1->Write();
 hist_2_pi_min2->Write();
 hist_2_pi_plus->Write();
 hist_2_kaon->Write();
 hist_2_lambda->Write();

 hist_3_p->Write();
 hist_3_pi_min1->Write();
 hist_3_pi_min2->Write();
 hist_3_pi_plus->Write();
 hist_3_kaon->Write();
 hist_3_lambda->Write();

 hist_cms_mom_p->Write();
 hist_cms_mom_pi_min1->Write();
 hist_cms_mom_pi_min2->Write();
 hist_cms_mom_pi_plus->Write();
 hist_cms_mom_kaon->Write();

 hist_cms_th_p->Write();
 hist_cms_th_pi_min1->Write();
 hist_cms_th_pi_min2->Write();
 hist_cms_th_pi_plus->Write();
 hist_cms_th_kaon->Write();  

 //c1->Write();
 


 
 f1->Close();
 
 


  
}
