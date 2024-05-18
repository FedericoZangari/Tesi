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
* Author: Han Miao
* Email: miaohan@ihep.ac.cn
*/
#include "DrawPosition.h"
#include "IR3Custom.h"
#include "IR3TaskManager.h"
#include "TCanvas.h"

/// Set the default address of the static instance of your algorithm to be 0
DrawPosition* DrawPosition::drawPosition = 0;

DrawPosition::DrawPosition(){

    /// Set the name of your algorithm
    SetModuleName("DrawPosition");

    /// Set address of the instance to be "this"
    drawPosition = this;

    /// Get an instance of IR3TaskManager and tell the task manager which variable to be delivered. The last 
    /// parameter is the type of your variable, which can be set as I(int), D(double) and S(string). If you 
    /// need to use the variables, you’d better set default values. Once the variables are delivered to
    /// task manager you can change them conveniently by modify the jobOption file.
    IR3TaskManager* taskManager = IR3TaskManager::Instance();
    taskManager->DeclareVar(this, "Max", &max, D);
    taskManager->DeclareVar(this, "Min", &min, D);

    /// Instance of IR3IO
    ir3io = IR3IO::Instance();
}

DrawPosition::~DrawPosition(){
    delete drawPosition;
}

/// The function to get the instance of your algorithm
DrawPosition* DrawPosition::Instance(){
    if(!drawPosition){
        drawPosition = new DrawPosition();
    }
    return drawPosition;
}

bool DrawPosition::Initialize(){

    /// Open an output file and create a histogram. Note: the output file here only for test and the IR3IO will responsible for the output in future.
    outfile = new TFile("hist.root", "RECREATE");
    positionY = new TH1F("positionY", "Y position of tacker 2_2", 100, min, max);

    return true;
}

bool DrawPosition::Execute(){
    /// Get the hit collection of the tracker 2-2
    vector<hit*>* hits2_2 = ir3io->GetHits(2, 2);
    positionY->Fill(hits2_2->at(0)->position.Y());

    //cout<<"Y position: "<<hits2_2->at(0)->position.Y()<<endl;

    return true;
}

bool DrawPosition::Finalize(){

    positionY->Write();

    outfile->Close();
    return true;
}