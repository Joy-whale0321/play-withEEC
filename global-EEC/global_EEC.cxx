#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PlotFile;
#endif

#ifndef __CINT__
#include <stdio.h>
#include "stdlib.h"
#include <iostream>
#include "Riostream.h"
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <fstream>
#include <complex>
#include <vector>
#include <algorithm>

#include "math.h"
#include "string.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TMath.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TGraph.h"
#endif

#include "AMPT.h"

using namespace std;

#define PI 3.1415926

void hadron_loop(AMPT* track,int i);
void try2correlate(vector<TVector3> Particles_v,TH1D* h1_global_eec);
void writeHistograms(char *outFile);

vector<TVector3> LightB;
vector<TVector3> LightM;
vector<TVector3> HeavyM;

double eta_cut=1.0;//MLCC
double delta_eta_cut=0.5;//MLCC                              
double pt_low_cut=0.3;//MLCC
double pt_high_cut=5.0;//MLCC

TH1D *h1_global_eec_allparticle;

int main(int argc, char **argv)
{
	if(argc!=3) return 0;

    TChain *chain1 = new TChain("AMPT");

    char *FileInput1=0;
	char *FileOutput=0;

    FileInput1 = argv[1];
	FileOutput = argv[2];

    int fileNumber = 0;
    char FileList[512];

    ifstream* inputStream1 = new ifstream;

    inputStream1->open(FileInput1);
    if (!(inputStream1))
    {
		printf("can not open FileInput1 list file\n");
		return 0;
    }

    for (;inputStream1->good();)
    {
		inputStream1->getline(FileList,512);
		if  ( inputStream1->good() )
		{
	    	TFile *ftmp1 = new TFile(FileList);
	    	if(!ftmp1||!(ftmp1->IsOpen())||!(ftmp1->GetNkeys()))
	    	{
				printf("input1 zpct0 file %s error in opening!!!\n",FileList);
	    	}
	    	else
	    	{
				printf("input1 zpct0 read in file %s\n",FileList);
				chain1->Add(FileList);
				fileNumber++;
	    	}
	    	delete ftmp1;
		}
    }
    printf(" files read in %d\n",fileNumber/3);


	AMPT* ampt1 = new AMPT(chain1);

	int mNEvents1 = (int)chain1->GetEntries();//number of events

	TH1::SetDefaultSumw2(1);

 	// cout<<"No. zpct0="<<mNEvents1<<" No. afterART="<<mNEvents2<<" No. after zpc="<<mNEvents3<<endl;

	h1_global_eec_allparticle = new TH1D("h1_global_eec_allparticle","h1_global_eec_allparticle",11000,-1.0, 10.0);

	for(int mNEventCount = 0; mNEventCount < mNEvents1; mNEventCount++)
	{
		chain1->GetEntry(mNEventCount);
		int nmult1=(int)ampt1->Event_multi;

		for(int i=0; i<nmult1; i++)
		{
		   	hadron_loop(ampt1,i);					  
		}

    	std::vector<TVector3> AllParticles;
    	AllParticles.insert(AllParticles.end(), LightB.begin(), LightB.end());
    	AllParticles.insert(AllParticles.end(), LightM.begin(), LightM.end());
    	AllParticles.insert(AllParticles.end(), HeavyM.begin(), HeavyM.end());

		try2correlate(AllParticles,h1_global_eec_allparticle);

		LightB.clear();
		LightM.clear();
		HeavyM.clear();
	}//end event loop

	writeHistograms(FileOutput);

	return 0;
}

void hadron_loop(AMPT* track,int i)
{
	//in >> id >> px >> py >> pz >> am >> others1 >>others2 >>others3 >>others4;//input track(particle) information
	
	int id = (int)track->ID[i];
	//double index = (double)track->Indx[i];
	double px = (double)track->Px[i];
	double py = (double)track->Py[i];
	double pz = (double)track->Pz[i];
	double x=(double)track->X[i];
	double y=(double)track->Y[i];
	double z=(double)track->Z[i];
	double am = (double)track->Mass[i];
	TVector3 mP(px,py,pz);

	if(mP.Pt()<0.3 || mP.Pt()>99.0) return;      

	double theta=acos(pz/sqrt(px*px+py*py+pz*pz));
	double pseorap=-log(tan(theta/2.));
	double pt=sqrt(px*px+py*py);
	double phi=atan2(py,px); //momentum coordinate MLCC!!

	if(fabs(pseorap)<eta_cut && (abs(id)==2212) && pt<5.0)//选出light重子存到数组
	{
		TVector3 lightbaryon_v(pseorap,phi,pt);//eta/phi/pt
		LightB.push_back(lightbaryon_v);
	}

	if(fabs(pseorap)<eta_cut && (abs(id)==321||abs(id)==211) && pt<5.0)//选出light介子存到数组
	{
		TVector3 lightmeson_v(pseorap,phi,pt);
		LightM.push_back(lightmeson_v);
	}

	if(fabs(pseorap)<eta_cut && (abs(id)==421||abs(id)==411||abs(id)==413) && pt<5.0)//选出D介子存到数组
	{
		TVector3 heavymeson_v(pseorap,phi,pt);
		HeavyM.push_back(heavymeson_v);
	}

	return;
}

void try2correlate(vector<TVector3> Particles_v,TH1D* h1_global_eec)
{
  	if((int)(Particles_v.size())==0) return;

	int number_pair = 0;

	for(int i=0;i<(int)(Particles_v.size());i++)
  	{
		for(int j=0;j<(int)(Particles_v.size());j++)
		{
		    if(i!=j && fabs(Particles_v[i].X()-Particles_v[j].X())>delta_eta_cut && Particles_v[i].Z()>pt_low_cut && Particles_v[i].Z()<pt_high_cut && Particles_v[j].Z()>pt_low_cut && Particles_v[j].Z()<pt_high_cut)
			{
				number_pair = number_pair + 1;
			}
		}
	}

	for(int i=0;i<(int)(Particles_v.size());i++)
  	{
		for(int j=0;j<(int)(Particles_v.size());j++)
		{
		    if(i!=j && fabs(Particles_v[i].X()-Particles_v[j].X())>delta_eta_cut && Particles_v[i].Z()>pt_low_cut && Particles_v[i].Z()<pt_high_cut && Particles_v[j].Z()>pt_low_cut && Particles_v[j].Z()<pt_high_cut)
			{
				if(fabs(Particles_v[i].Y()-Particles_v[j].Y())>PI)
				{
					double delta_eta = fabs(Particles_v[i].X()-Particles_v[j].X());
					double delta_phi = 2.0*PI - fabs(Particles_v[i].Y()-Particles_v[j].Y());
					double RL = sqrt(delta_phi*delta_phi + delta_eta*delta_eta);

					double energy_factor = Particles_v[i].Z() * Particles_v[j].Z();

					h1_global_eec->Fill(RL,energy_factor);
				}
				else
				{
					double delta_eta = fabs(Particles_v[i].X()-Particles_v[j].X());
					double delta_phi = fabs(Particles_v[i].Y()-Particles_v[j].Y());
					double RL = sqrt(delta_phi*delta_phi + delta_eta*delta_eta);

					double energy_factor = Particles_v[i].Z() * Particles_v[j].Z();

					h1_global_eec->Fill(RL,energy_factor);
				}
			}
		}
  	}

  	return;
}

void writeHistograms(char *outFile)//histogram write to file
{
	char file[256];
	sprintf(file,"%s",outFile);
	TFile *f = new TFile(file,"RECREATE");

	h1_global_eec_allparticle->Write();

	f->Write(); ////!!!!!!!!
	f->Close();
	delete f;
	return;

}



