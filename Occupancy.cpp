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
#include "Occupancy.h"
#include "IR3Custom.h"
#include "IR3TaskManager.h"
#include "TCanvas.h"
#include <stdlib.h>
#include <stdio.h>
#include <chrono>

using namespace std::chrono;

/// Set the default address of the static instance of your algorithm to be 0
Occupancy* Occupancy::occupancy = 0;

Occupancy::Occupancy(){

    /// Set the name of your algorithm
    SetModuleName("Occupancy");

    /// Set address of the instance to be "this"
    occupancy = this;

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
    taskManager->DeclareVar(this, "Trk_id", &trk_id, I);
    taskManager->DeclareVar(this, "Trk_nmr", &trk_nmr, I);
    /// Instance of IR3IO
    ir3io = IR3IO::Instance();
}

Occupancy::~Occupancy(){
    delete occupancy;
}

/// The function to get the instance of your algorithm
Occupancy* Occupancy::Instance(){
    if(!occupancy){
        occupancy = new Occupancy();
    }
    return occupancy;
}

bool Occupancy::Initialize(){

    /// Open an output file and create a histogram. Note: the output file here only for test and the IR3IO will responsible for the output in future.
    outfile = new TFile("hist_Lc_1000.root", "UPDATE"); // root file with Lambda _c
    //outfile = new TFile("hist_mb_1000.root", "UPDATE");   //form mb events

    mat_dim_half_x = pix_dim * n_pix_x / 2. ;
    mat_dim_half_y = pix_dim * n_pix_y / 2. ;

    

    grid = new vector<hit*>**[n_pix_x];
    energ_depo = new double*[n_pix_x];
    trigged = new int*[n_pix_x];

    // initialize the matrices
    for(int i=0; i<n_pix_x; i++){
        grid[i] = new vector<hit*>*[n_pix_y];
        energ_depo[i] = new double[n_pix_y];
        trigged[i] = new int[n_pix_y];

        for(int j=0; j<n_pix_y; j++){

            grid[i][j] = new vector<hit*>();
            energ_depo[i][j]  = 0.;
            trigged[i][j]  = 0;
        }
    }
    
    

    return true;
}

bool Occupancy::Execute(){
    auto start = high_resolution_clock::now();
    Event_number++;

    /// Get the hit collection of the tracker trk_id-trk_nmr
    vector<hit*>* Hits = ir3io->GetHits(trk_id, trk_nmr);
    
 
    //  event info
    cout << "Event number \t" << Event_number << endl;
    cout << "Number of hits in this event: " << Hits->size() << endl;
    double x, y;
    int i, j;
    
    if (Hits->size() != 0){ for(Int_t k = 0; k < Hits->size(); k++){
        //assign the hit to its position in the pixel grid
        x = Hits->at(k)->position.X();
        y = Hits->at(k)->position.Y();
        i = floor(x/pix_dim + n_pix_x/2.);
        j = floor(-y/pix_dim + n_pix_y/2.);
        
        //check if the hit is inside the tracker matrix
        if(x > -mat_dim_half_x && x < mat_dim_half_x && y > -mat_dim_half_y && y < mat_dim_half_y){
            grid[i][j]->push_back( Hits->at(k));
        
            // evaluate total energy deposit in the pixel
            energ_depo[i][j] += grid[i][j]->back()->energyDeposit;
        }
    }}
    
    for(int ii = 0 ; ii< n_pix_x; ii++){
        for(int jj = 0 ; jj < n_pix_y ; jj++){
            
            // check if the pixel is triggered
            if(energ_depo[ii][jj] > threshold) {
                trigged[ii][jj] ++;
            }

            // reset the matrices
            grid[ii][jj]->clear();
            energ_depo[ii][jj]  = 0.;
        }
    }
/*
    for(int i = 0 ; i< n_pix_x; i++){
    //parallel_for (0, n_pix_x, [] int i){
        for(int j = 0 ; j < n_pix_y ; j++){

            // fill the grid with the hits in the considered pixel
            if (Hits->size() != 0){ for(Int_t k = 0; k < Hits->size(); k++){
                if( Hits->at(k)->position.X() < -mat_dim_half_x + (i+1) * pix_dim  && Hits->at(k)->position.X() > -mat_dim_half_x + i * pix_dim &&
                    Hits->at(k)->position.Y() < mat_dim_half_y - j * pix_dim       && Hits->at(k)->position.Y() > mat_dim_half_y - (j+1) * pix_dim   ){
                        grid[i][j]->push_back( Hits->at(k));
                    }    
            }}

            // evaluate total energy deposit in the pixel
            if (Hits->size() != 0) for(int k = 0; k < grid[i][j]->size(); k++){
                  energ_depo[i][j] += grid[i][j]->at(k)->energyDeposit;
            }

            // check if the pixel is triggered
            if(energ_depo[i][j] > threshold) {
                trigged[i][j] ++;
            }

            // reset the matrices
            grid[i][j]->clear();
            energ_depo[i][j]  = 0.;
        }
    }
*/
     
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Time duration of this loop: " << setprecision(4) << duration.count() << " milliseconds "<< endl << endl;
    
    return true;
}

bool Occupancy::Finalize(){

    //string title = "Occupancy % of each pixel, tracker " + to_string(trk_id) + "_" + to_string(trk_nmr) + ", number events " + to_string(Event_number);
    string title = "Occupancy % of each pixel, tracker " + to_string(trk_id) + "_" + to_string(trk_nmr) + ", " + to_string(Event_number) +" #Lambda_{c}^{+} events";
    string title2 = "Occupancy" + to_string(trk_id) + "_" + to_string(trk_nmr) + "_" + to_string(Event_number);
    occupancy_pix = new TH2F(title2.c_str(), title.c_str() , n_pix_x, -mat_dim_half_x, mat_dim_half_x, n_pix_y, -mat_dim_half_y, mat_dim_half_y);

    // fill the 2D Histo
    for(int i=0; i<n_pix_x; i++){
        for(int j=0; j<n_pix_y; j++){
            occupancy_pix->SetBinContent(i, n_pix_y-1-j, trigged[i][j] * 100 / static_cast<double>(Event_number) );     
        }
    }
    occupancy_pix->Write();
    
    outfile->Close();

    delete []grid;
    delete []energ_depo;
    delete []trigged;

    return true;
}