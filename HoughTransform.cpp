/* ===================================================================
* IR3ana: the analysis framework for IR3 test
* Copyright (C) 2022  IR3 test software group
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <https://www.gnu.org/licenses/>.
* ======================================================================
*
* Author: Federico Zangari
* Email: federico.zangari@cern.ch
*/
#include "HoughTransform.h"
#include "IR3Custom.h"
#include "IR3TaskManager.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include <stdlib.h>
#include <stdio.h>
#include <chrono>

using namespace std::chrono;

/// Set the default address of the static instance of your algorithm to be 0
HoughTransform* HoughTransform::houghtransform = 0;

HoughTransform::HoughTransform(){

    /// Set the name of your algorithm
    SetModuleName("HoughTransform");

    /// Set address of the instance to be "this"
    houghtransform = this;

    /// Get an instance of IR3TaskManager and tell the task manager which variable to be delivered. The last 
    /// parameter is the type of your variable, which can be set as I(int), D(double) and S(string). If you 
    /// need to use the variables, you’d better set default values. Once the variables are delivered to
    /// task manager you can change them conveniently by modify the jobOption file.
    IR3TaskManager* taskManager = IR3TaskManager::Instance();
    
    //taskManager->DeclareVar(this, "Max", &max, D);
    //taskManager->DeclareVar(this, "Min", &min, D);
    taskManager->DeclareVar(this, "N_pix_x", &n_pix_x, I);
    taskManager->DeclareVar(this, "N_pix_y", &n_pix_y, I);
    taskManager->DeclareVar(this, "Pix_dim", &pix_dim, D);
    taskManager->DeclareVar(this, "Threshold", &threshold, D);
    taskManager->DeclareVar(this, "D", &D_trackers, D);
    /// Instance of IR3IO
    ir3io = IR3IO::Instance();
}

HoughTransform::~HoughTransform(){
    delete houghtransform;
}

/// The function to get the instance of your algorithm
HoughTransform* HoughTransform::Instance(){
    if(!houghtransform){
        houghtransform = new HoughTransform();
    }
    return houghtransform;
}

bool HoughTransform::Initialize(){

    /// Open an output file and create a histogram. Note: the output file here only for test and the IR3IO will responsible for the output in future.
    outfile = new TFile("hist_Lc_1000_hough.root", "RECREATE");

    mat_dim_half_x = pix_dim * n_pix_x / 2. ;
    mat_dim_half_y = pix_dim * n_pix_y / 2. ;

    

    
    energ_depo = new double*[n_pix_x];

    // initialize the matrices
    for(int i=0; i<n_pix_x; i++){
        
        energ_depo[i] = new double[n_pix_y];

        for(int j=0; j<n_pix_y; j++){

            
            energ_depo[i][j]  = 0.;
        }
    }
    
    

    return true;
}

bool HoughTransform::Execute(){
    auto start = high_resolution_clock::now();
    Event_number++;

    /// Get the hit collection of the tracker trk_id-trk_nmr
    vector<hit*>* v = ir3io->GetHits(1, 1);
    //  event info
    cout << "Event number \t" << Event_number << endl;
    cout << "Number of hits in this event: " << v->size() << endl;
    double x, y, z;
    int i, j;
    double size = 0.002;

    string canvas_title = "c_" + to_string(Event_number);
    TCanvas *c = new TCanvas(canvas_title.c_str(), canvas_title.c_str(), 900, 600);
    
    vector<TF1*>* vec = new vector<TF1*>;
    int kk = 0;

    for(int q = 1; q < 9; q++){
        cout << "Tracker " << q << endl;
        if(q == 2) v = ir3io->GetHits(1, 2);
        if(q == 3) v = ir3io->GetHits(1, 3);
        if(q == 4) v = ir3io->GetHits(1, 4);
        if(q == 5) v = ir3io->GetHits(2, 1);
        if(q == 6) v = ir3io->GetHits(2, 2);
        if(q == 7) v = ir3io->GetHits(2, 3);
        if(q == 8) v = ir3io->GetHits(2, 4);

        if (v->size() != 0){ for(Int_t k = 0; k < v->size(); k++){
            //assign the hit to its position in the pixel grid
            x = v->at(k)->position.X();
            y = v->at(k)->position.Y();
            cout << "x real: " << x << endl;
            i = floor(x/pix_dim + n_pix_x/2.);
            j = floor(-y/pix_dim + n_pix_y/2.);
            
            //check if the hit is inside the tracker matrix
            if(x > -mat_dim_half_x && x < mat_dim_half_x && y > -mat_dim_half_y && y < mat_dim_half_y){
           
                // evaluate total energy deposit in the pixel
                energ_depo[i][j] += v->at(k)->energyDeposit;
            }
        }}
        
        for(int ii = 0 ; ii< n_pix_x; ii++){
            for(int jj = 0 ; jj < n_pix_y ; jj++){
                
                // check if the pixel is triggered
                if(energ_depo[ii][jj] > threshold ) {
                    x = -mat_dim_half_x + ii*pix_dim + pix_dim/2.;
                    
                    if(q < 5) z = 20.+70.+1000.-D_trackers-50. + (q-1)*D_trackers/3. + 0.1;
                    if(q > 4) z = 20.+70.+1000.+ 1700. + 50. + (q-5)*D_trackers/3. + 0.1;
                    cout << "z: " << z <<endl;
                    cout << "x: " << x <<endl;
                    //cout << "-z: "<< -z << endl;
                    //TF1 *f = new TF1("f","[0]*x+[1]",-1.125,1.125);
                    TF1 *f = new TF1("f","[0]*x+[1]",-size, size);
                    f->SetParameter(0,-z);
                    f->SetParameter(1,x);
                    f->SetParError(1,2);

                    vec->push_back(f);
                }
                // reset the matrices                
                energ_depo[ii][jj]  = 0.;
            }
        }
    }
    
    
    
    cout << vec->size() << endl;
    double n = 1000;
    double step = 2.*size/n;

    for(Int_t k = 0; k < vec->size(); k++){
        //if(k%3 == 0) vec->at(k)->SetLineColor(kRed);
        //if(k%3 == 1) vec->at(k)->SetLineColor(kGreen);
        //if(k%3 == 2) vec->at(k)->SetLineColor(kBlue);
        if(k==0){
            vec->at(0)->SetTitle("Hough transform");
            vec->at(0)->GetXaxis()->SetTitle("angular coefficient m_{x}");
            vec->at(0)->GetYaxis()->SetTitle("quota q_{x} [mm]");
            vec->at(0)->GetYaxis()->SetRangeUser(-1,1);
            vec->at(0)->Draw();
        }

        vec->at(k)->Draw("same");
        TGraphErrors *gr = new TGraphErrors(n);
        for (int ii=0; ii<n; ii++){
            double m = -size + step*ii;
            double q = vec->at(k)->GetParameter(0)*m + vec->at(k)->GetParameter(1);
            gr->SetPoint(ii, m, q);
            gr->SetPointError(ii, 0, pix_dim/2.);
        }
        gr->SetFillColor(kRed);
        //gr->SetFillStyle(3010);
        gr->Draw("same3");
    }
     
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Time duration of this loop: " << setprecision(4) << duration.count() << " milliseconds "<< endl << endl;
    canvas_title = canvas_title +  ".png";
    c->SaveAs(canvas_title.c_str());
    c->Write();
    return true;
}

bool HoughTransform::Finalize(){

    
    
    outfile->Close();

    delete []energ_depo;

    return true;
}