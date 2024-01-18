// Copyright (C) 2024 Paolo Massioni paolo.massioni@insa-lyon.fr

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

// If you use this work for research purposes, please cite 
// Hiba Houmsi, Paolo Massioni, Federico Bribiesca Argomedo, Romain Delpoux. "Embedded 
// controller optimization for efficient electric motor drive," 2023 IEEE Vehicle Power
// and Propulsion (VPPC 2023), Oct 2023, Milan, Italy. 

// This code shows an example of embedded solver for a Semi-Definite Programming problem,
// a Linear Matrix Inequality (LMI) feasibility problem from control engineering, 
// namely the so-called regional pole placement for a state-space system with state
// feedback, as introduced in M. Chilali, P. Gahinet and P. Apkarian, "Robust pole 
// placement in LMI regions," IEEE Transactions on Automatic Control, vol. 44, no. 12,
// pp. 2257-2270, Dec. 1999, doi: 10.1109/9.811208. The solver is a simple implementation
// of an interior point method as described in Stephen Boyd, Laurent El Ghaoui, E. Feron, 
// and V. Balakrishnan, "Linear Matrix Inequalities in System and Control Theory",
// Volume 15 of Studies in Applied Mathematics, Society for Industrial and Applied 
// Mathematics (SIAM), 1994. 

// The code provides a solution to a specific control engineering promblem, namely
// the electric motor control problem presented in Hiba Houmsi, Paolo Massioni, 
// Federico Bribiesca Argomedo, Romain Delpoux, "Embedded controller optimization for 
// efficient electric motor drive," 2023 IEEE Vehicle Power and Propulsion (VPPC 2023),
// Oct 2023, Milan, Italy. This code has been developed in the context of this
// last work. To sumarize, the code considers a state-space dynamics of the kind
// dx/dt = Ax + Bu, and it computes a controller gain K such that for a state 
// feedbck law u = Kx the closed-loop systems features poles in the regions described 
// by the parameters ALPHAMAX, ALPHAMIN and BETA (see the paper for details).

// This code requires the use of the eigen c++ library. To my knowledge, this is not
// supported by all Arduino boards. This code has been tested to work on Arduino 
// MKR1000 and Arduino MKR1010 Wifi.

#include <ArduinoLMI.h>

#define MATRIXA -1874.3,-0.0264,0.0,  3960.0,-1.0,0.0,  0.0,-1.0,0.0 // problem-specific constant (A matrix)
#define MATRIXA1 -100.0,0.0,0.0, 0.0,0.0,0.0,  0.0,0.0,0.0 // problem-specific constant (A matrix)
#define MATRIXB 2857.1,0.0,0.0 // problem-specific constant (B matrix)
#define MATRIXN 3  // problem-specific constant (number of state variables)
#define MATRIXM 1  // problem-specific constant (number of input variables)
#define ALPHAMAX 1200.0 // problem-specific constant
#define ALPHAMIN 400.0 // problem-specific constant
#define BETA 1.5 // problem-specific constant

void print_mtxf(const Eigen::MatrixXf& K);  

void setup() {
 Serial.begin(9600);
  while (!Serial) {
    ; // wait for serial port to connect
  }
 Serial.println("Ready."); 
 Serial.println();
}


void loop() {
  long timer;
  Serial.println("Send a key through serial com in order to start the computation.");
  Serial.println();

  while(!Serial.available());
  while(Serial.available()) {
     Serial.read();
   }    

  Serial.println("Computing...");
  Serial.println();

  //MatrixXf A2(N,N); // test matrix, delete when finished

  MatrixXf A(MATRIXN,MATRIXN); // problem-specific definition
  MatrixXf B(MATRIXN,MATRIXM); // problem-specific definition 
  A << MATRIXA;
  B << MATRIXB;
  
  /// nominal regional pole placement example
  timer=millis();
  Eigen::MatrixXf K= RegionalPolePlacement(A, B, ALPHAMAX, ALPHAMIN, BETA);  
  timer=millis()-timer;
  Serial.print("Elapsed time (ms): " );
  Serial.println(timer);   
  Serial.println(); 
  Serial.println("K");
  print_mtxf(K);          
  Serial.println("A+B*K");
  print_mtxf(A+B*K);          
  Serial.println("Copy and paste this matrix into your favourite computational software to verify that for");
  Serial.println("each eigenvalue lambda, -ALPHAMAX < Re(lambda) < -ALPHAMIN, |Im(lambda)| < BETA |Re(lambda)|");
  Serial.println();

  /// robust regional pole placement example
  //MatrixXf A1(MATRIXN,MATRIXN); // problem-specific definition
  //A1 << MATRIXA1;
  //timer=millis();
  //Eigen::MatrixXf Kr= RegionalPolePlacementRobust1Param(A, A1, B, ALPHAMAX, ALPHAMIN, BETA, 1,0);  
  //timer=millis()-timer;
  //Serial.print("Elapsed time (ms): " );
  //Serial.println(timer);   
  //Serial.println();
  //Serial.println("A+B*K");
  //print_mtxf(A+B*Kr);      // p=0
  //print_mtxf(A+A1+B*Kr);   // p=1
  
  /// gain scheduling regional pole placement example
  //MatrixXf A1(MATRIXN,MATRIXN); // problem-specific definition
  //A1 << MATRIXA1;
  //timer=millis();
  //Eigen::MatrixXf** Kp= RegionalPolePlacementGainScheduling1Param(A, A1, B, ALPHAMAX, ALPHAMIN, BETA,1,0);  
  //timer=millis()-timer;
  //Serial.print("Elapsed time (ms): " );
  //Serial.println(timer);   
  //Serial.println();
  //Serial.println("A+B*K");
  //print_mtxf(A+B*(*Kp[0]));            // p=0
  //print_mtxf(A+A1+B*(*Kp[1]+*Kp[0]));  // p=1
  //delete Kp[1];
  //delete Kp[0];
}


// PRINT MATRIX (float type)
// By: randomvibe
//-----------------------------
void print_mtxf(const Eigen::MatrixXf& X)  
{
    int i, j, nrow, ncol;
    
    nrow = X.rows();
    ncol = X.cols();

    Serial.print("nrow: "); Serial.println(nrow);
    Serial.print("ncol: "); Serial.println(ncol);       
    Serial.println();
    
    for (i=0; i<nrow; i++)
    {   
        for (j=0; j<ncol; j++)
        {
            Serial.print(X(i,j), 6);   // print 6 decimal places
            Serial.print(", ");
        }
        Serial.println();
    }
    Serial.println();
}
