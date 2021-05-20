// Encapsulate all of Vladimir Bashkirov pedestal and gain on-the-fly
// calibration
// R.P. Johnson   5/22/2016

#include "pedGainCalib.h"
#include "TFile.h"
#include "TROOT.h"
#include "TFitResultPtr.h"
pedGainCalib::~pedGainCalib()
{
}

pedGainCalib::pedGainCalib(TFile* root, int pedMin[5], float oldPed[5], float t1, float t2, float t3, float t4):RootFile(root)
{  
  theConfig = pCTconfig::GetInstance();
  // Two ranges in t occupied by unobstructed (empty) protons
  emtrng1 = t1; // negative t side
  emtrng2 = t2;
  emtrng3 = t3; // positive t side
  emtrng4 = t4;
  for (int stage = 0; stage < 5; stage++) {
    Ped[stage] = oldPed[stage];
    GainFac[stage] = 1.0;

    cout << "pedGainCalib setting the default pedestal for stage " << stage << " to " << Ped[stage] << endl;
  }
  // Define histograms for pedestal calculation
  for (int i =0; i<5; i++) hPed[i] = new TH1D(Form("PedestalStage_%i",i), Form("Pedestal region for stage %i",i), 400, pedMin[i], pedMin[i] +400*5);
  for (int i =0; i<5; i++) hPed2D[i] = new TH2D(Form("PedestalStage_vsTime_%i",i), Form("Pedestal region for stage %i vs Time",i), 400, pedMin[i], pedMin[i] +400*5, 20, 0, 5e9);
  for (int i =0; i<5; i++) hTotFil[i] = new TH1D(Form("FullADCStage_Filtered_%i",i), Form("Full ADC for stage %i",i), 400, pedMin[i], pedMin[i] +400*5);
  for (int i =0; i<5; i++) hTotUnFil[i] = new TH1D(Form("FullADCStage_Unfiltered_%i",i), Form("Full ADC for stage %i",i), 400, pedMin[i], pedMin[i] +400*50);


  // Define histograms for pedestal calculation
  for (int i =0; i<5; i++) hPed_In[i] = new TH1D(Form("PedestalStage_%i_In",i), Form("Pedestal region for stage %i",i), 400, pedMin[i], pedMin[i] +400*5);
  for (int i =0; i<5; i++) hPed_Out[i] = new TH1D(Form("PedestalStage_%i_Out",i), Form("Pedestal region for stage %i",i), 400, pedMin[i], pedMin[i] +400*5);

  // Define histograms for gain calibration
  if (theConfig->item_str["partType"] == "H") {
    for (int i =0; i<5; i++) hEnrg[i] = new TH1D(Form("EnergyDistribution_%i",i), Form("Energy Distribution for stage %i",i), 800, 15, 15 + 800*0.175);
    hEnrgTot = new TH1D("SumStageEnergies", "Sum of stage energies", 800, 0, 0 + 0.3*800);
  }
  else {
    for (int i =0; i<5; i++) hEnrg[i] = new TH1D(Form("EnergyDistribution_%i",i), Form("Energy Distribution for stage %i",i), 800, 15, 15 + 800*0.9);
    hEnrgTot = new TH1D("SumStageEnergies", "Sum of stage energies", 800, 0, 0 + 1.2*800);
  }
  // Profile plot to make sure that the phantom does not extend into the regions used for gain calibration
  hProfT = new TProfile2D("Stage0EnergyProfile", "Stage 0 energy profile", 100, -150, -150 + 100*3.0, 100, -50., -50 + 1.0*100);
  hTedet = new TH1D("T_Ions_GainRecalib", "T of ions used for gain recalibration", 100, -150, -150 + 100*3.0);
  RootFile->mkdir("Pedestals");
}

void pedGainCalib::FillPeds(pCTraw &rawEvt) { // Called for each raw event read in from the input data file
  if(first){ startTime = rawEvt.time_tag; first=false;}
  // Accumulate histograms for pedestal analysis  
  hPed[0]->Fill(rawEvt.enrg_fpga[0].pulse_sum[0]);
  hPed[1]->Fill(rawEvt.enrg_fpga[0].pulse_sum[1]);
  hPed[2]->Fill(rawEvt.enrg_fpga[0].pulse_sum[2]);
  hPed[3]->Fill(rawEvt.enrg_fpga[1].pulse_sum[0]);
  hPed[4]->Fill(rawEvt.enrg_fpga[1].pulse_sum[1]);
  hPed2D[0]->Fill(rawEvt.enrg_fpga[0].pulse_sum[0],rawEvt.time_tag - startTime);
  hPed2D[1]->Fill(rawEvt.enrg_fpga[0].pulse_sum[1],rawEvt.time_tag - startTime);
  hPed2D[2]->Fill(rawEvt.enrg_fpga[0].pulse_sum[2],rawEvt.time_tag - startTime);
  hPed2D[3]->Fill(rawEvt.enrg_fpga[1].pulse_sum[0],rawEvt.time_tag - startTime);
  hPed2D[4]->Fill(rawEvt.enrg_fpga[1].pulse_sum[1],rawEvt.time_tag - startTime);
  hTotUnFil[0]->Fill(rawEvt.enrg_fpga[0].pulse_sum[0]);
  hTotUnFil[1]->Fill(rawEvt.enrg_fpga[0].pulse_sum[1]);
  hTotUnFil[2]->Fill(rawEvt.enrg_fpga[0].pulse_sum[2]);
  hTotUnFil[3]->Fill(rawEvt.enrg_fpga[1].pulse_sum[0]);
  hTotUnFil[4]->Fill(rawEvt.enrg_fpga[1].pulse_sum[1]);
}
  
void pedGainCalib::ResetHist(){
  for(int i =0; i<5; i++){
    hTotFil[i]->Reset();
    hTotUnFil[i]->Reset();
    hPed[i]->Reset();
    hPed_Out[i]->Reset();
    hPed_In[i]->Reset();
    hEnrg[i]->Reset();
  }
  hEnrgTot->Reset();
  hTedet->Reset();
  hProfT->Reset();
}

void pedGainCalib::FillADC(int ADC[5]){ for(int i =0; i<5; i++) hTotFil[i]->Fill(ADC[i]);}

void pedGainCalib::WriteHist(){ // For analysis
  inFileName_s = theConfig->item_str["inputFileName"];
  inFileName_s = inFileName_s.substr(inFileName_s.find_last_of("\\/") +1,  inFileName_s.size());
  RootFile->mkdir(Form("Pedestals/%s",inFileName_s.c_str()));
  RootFile->cd(Form("Pedestals/%s",inFileName_s.c_str()));
  for(int i =0; i<5; i++) {
    hTotFil[i]->Write("", TObject::kOverwrite);
    hTotUnFil[i]->Write("", TObject::kOverwrite);
    hPed[i]->Write("", TObject::kOverwrite);
    hPed_Out[i]->Write("", TObject::kOverwrite);
    hPed_In[i]->Write("", TObject::kOverwrite);
    hPed2D[i]->Write("",TObject::kOverwrite);
    hEnrg[i]->Write("", TObject::kOverwrite);
    hTotFil[i]->Write("", TObject::kOverwrite);
  }
  hEnrgTot->Write("",TObject::kOverwrite);
  hTedet->Write("", TObject::kOverwrite);
  hProfT->Write("", TObject::kOverwrite);  
}

void pedGainCalib::GetPeds() {
  // Calculate the pedestal
  for (int stage = 0; stage < 5; stage++) {
    float xLow, xHigh;
    float max, xmax, xpeak;
    int imax = hPed[stage]->GetMaximumBin();
    max   = hPed[stage]->GetMaximum();
    xmax  = hPed[stage]->GetBinCenter(imax);
    TFitResultPtr r = hPed[stage]->Fit("gaus","QNS","",xmax-10,xmax+10); // Q means quiet, R is using TF1 range
    if(int(r)==0) // Everything passed
    {
      Ped[stage]  = r->Parameter(1);// mean
      stdPed[stage]  = r->Parameter(2);// std
      FivepercentPed[stage]=stdPed[stage]; 
      for(int i=imax; i<=imax+30; i++) {
	if(hPed[stage]->GetBinContent(i)>=max*0.05) FivepercentPed[stage] = hPed[stage]->GetBinCenter(i)-Ped[stage]; 
      }
    }    
    else
    {
	cout <<  "WARNING: Fit of peds did not converge, will manually compute the peak from FWHM" << endl; 
	int xlow,xhigh;
	for(int i = imax-100; i<=imax+100; i++){ // First loop to find FWHM (LastBinAbove sensitive to rise of the ADC counts at low E)
	    if(i<imax && hPed[stage]->GetBinContent(i)<=max/2) xlow = i; // will stay at the last value 
	    if(i>imax && hPed[stage]->GetBinContent(i)>=max/2) xhigh = i; // will stay at teh last value
	}
	stdPed[stage] = (hPed[stage]->GetBinCenter(xhigh)-hPed[stage]->GetBinCenter(xlow))/2.355; // FWHM to sigma for Gauss
	int cnt=0; 
	Ped[stage]=0;
	for(int i = xlow; i<xhigh; i++){
		cnt += hPed[stage]->GetBinContent(i);
		Ped[stage] += hPed[stage]->GetBinContent(i)*hPed[stage]->GetBinCenter(i);
	}
	if(cnt>0) Ped[stage]/=cnt; 
	else Ped[stage] = 0;
     }
     if(Ped[stage]==0) cout << "getPeds ERROR: Could not find a peak in the pedestal for stage" << stage << endl;
     else cout<< inFileName_s <<" "<<hPed[stage]->Integral()<< " getPeds: measured pedestal for stage " << stage << " is " << Ped[stage] << " with a width of "<<stdPed[stage]<<endl;
  }
}

void pedGainCalib::FillGains(float V[4], float T[4], float Ene[5], int phSum[5]) {// Called for each event in the temporary file of proton histories -- gain recalibration
  float Esum = 0.;
  //between t1 and t2 or between t3 and t4
  if ((T[0] > emtrng1 && T[0] < emtrng2) || (T[0] > emtrng3 && T[0] < emtrng4) &&
	(T[2] > emtrng1 && T[2] < emtrng2) || (T[2] > emtrng3 && T[2] < emtrng4) ) { //protons outside of the phantom region
    if (fabs(V[3])< 30.) { // not outside the detector [3]
      for (int stage = 0; stage < 5; stage++) {
        hEnrg[stage]->Fill(Ene[stage]);
        Esum += Ene[stage];
	hPed_Out[stage]->Fill(phSum[stage]);
      }
      hEnrgTot->Fill(Esum);
      hTedet->Fill(T[3]); //[3]
    }
  }
  else if(T[0]<-50){ // Analyze full energy particles inside the phantom region [0]
    if (fabs(V[3]) < 40.) for(int stage =0; stage<5; stage++)hPed_In[stage]->Fill(phSum[stage]);//[3]
  }    
  
  if (Ene[0] > 10.) hProfT->Fill(T[3], V[3], Ene[0]); //[3]
}
void pedGainCalib::GetGains(TVcorrection *TVcorr) {
                            
  // Calculate the gain
 for (int stage = 0; stage < 5; stage++) {
    float xLow, xHigh;
    float max, xmax, xpeak,min;
    int imax;
    max   =  hEnrg[stage]->GetMaximum();
    min   =  hEnrg[stage]->GetMinimum();
    imax  = hEnrg[stage]->GetMaximumBin();
    xmax  =  hEnrg[stage]->GetBinCenter(imax);
    TFitResultPtr r =hEnrg[stage]->Fit("gaus","QNS","",xmax-10,xmax+10); // Q means quiet,
    if(int(r)==0) // Everything passed
      {
	Peak[stage]  = r->Parameter(1);// mean
	stdGain[stage]   = r->Parameter(2);//std
      }
    else
      {
	cout <<  "WARNING: Fit of peds did not converge, will manually compute the peak from FWHM" << endl; 
	int xlow,xhigh;
	for(int i = imax-50; i<=imax+50; i++){ // First loop to find FWHM (LastBinAbove sensitive to rise of the ADC counts at low E)
	    if(i<imax && hEnrg[stage]->GetBinContent(i)<=max/2) xlow = i; // will stay at the last value 
	    if(i>imax && hEnrg[stage]->GetBinContent(i)>=max/2) xhigh = i; // will stay at teh last value
	}
	stdGain[stage] = (hEnrg[stage]->GetBinCenter(xhigh)-hEnrg[stage]->GetBinCenter(xlow))/2.355; // FWHM to sigma for Gauss
	int cnt=0; 
	Peak[stage]=0;
	for(int i = xlow; i<xhigh; i++){
		cnt += hEnrg[stage]->GetBinContent(i);
		Peak[stage] += hEnrg[stage]->GetBinContent(i)*hEnrg[stage]->GetBinCenter(i);
	}
	if(cnt>0) Peak[stage]/=cnt; 
	else Peak[stage]=0;
      }
      //else Peak[stage] = 0.0;
    cout << "getGains: measured peak location for stage " << stage << " is " << Peak[stage] << " and a width of "<<stdGain[stage]<<endl;
    if (Peak[stage] > 0.01) GainFac[stage] = TVcorr->Eempt[stage] / Peak[stage];
    else {
      cout << "getGains ERROR: unable to find a peak position; leaving the gains unchanged. **********" << endl;
      GainFac[stage] = 1.0;
    }
  }

}
