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

#include <ArduinoEigenDense.h>
using namespace Eigen;    // Eigen related statement; simplifies syntax for declaration of matrices

#define TOLERANCE 0.001 // solver tolerance (>0)
#define FEASTOL 0.1// solver tolerance for feasibility (>0)
#define THETA 0.05 // solver convergence parameter (>0)


void LMISolver(MatrixXf *F[], MatrixXf *V[], VectorXf& x, int sizeunk, int *sizecon, int NUNK, int NCON) {
  MatrixXf* FF[NCON];
  VectorXf grad(sizeunk);
  MatrixXf hess(sizeunk,sizeunk);


  int it=0;
  int itin=0;
  float l0=(1.0+THETA)*x(0);

  for (int i=0;i<NCON;i++) {
    FF[i]=new MatrixXf(sizecon[i],sizecon[i]);
  }
  
  for(int bailout0=0;!bailout0; it++) {
    VectorXf x0=x;
    
    for(int bailout1=0;!bailout1; itin++) {
      VectorXf x00=x;
      grad.setZero();
      hess.setZero();

      // compute gradient and hessian
      for (int c=0; c<NCON; c++) {
        FF[c]->noalias()=*F[c*(sizeunk+1)+sizeunk];
        for (int u=0; u<sizeunk; u++) {
          *FF[c]+=x(u)*(*F[c*(sizeunk+1)+u]);
        }
              
      }
      for (int c=0; c<NCON; c++) {      
        MatrixXf FFINVc=FF[c]->inverse();
        for (int u=0; u<sizeunk; u++) {
                V[c*(sizeunk+1)+u]->noalias()=FFINVc*(*F[c*(sizeunk+1)+u]);
                grad(u)=grad(u)-(V[c*(sizeunk+1)+u]->trace());
        } 
      }
      grad(0)+= 1.0/(l0-x(0));

      for (int u1=0; u1<sizeunk; u1++) {
        for (int u2=u1; u2<sizeunk; u2++) {
          for (int c=0; c<NCON; c++) {
            MatrixXf temp=(*V[c*(sizeunk+1)+u1])*(*V[c*(sizeunk+1)+u2]);
            hess(u1,u2)+=temp.trace();
            if ((u1!=u2) && (c==NCON-1)){
              hess(u2,u1)=hess(u1,u2);
            }            
          }
        }
      }
      hess(0,0)+= 1.0/sq(l0-x(0));
 //    Serial.print("hessian: ");  
 //       print_mtxf(hess);
      x+= -hess.llt().solve(grad);
      //Serial.print("x0: ");
      //Serial.print(x(0));
      VectorXf delta1=x-x00;
      //Serial.print(" delta1: ");      
      //Serial.println(      delta1.squaredNorm());
      if ((delta1.squaredNorm()<TOLERANCE*TOLERANCE )|| (x(0)<-FEASTOL)|| (itin>=200)){
        bailout1=1;
      }
    }
    VectorXf delta0=x-x0;
    //Serial.print(" delta2: ");      
    //Serial.println(      delta0.squaredNorm());
    if ((delta0.squaredNorm()<TOLERANCE*TOLERANCE )|| (x(0)<-FEASTOL)|| (itin>=200)){
      bailout0=1;      
    } else {
      l0=THETA*(x0(0))+(1-THETA)*x(0);
    }   
  }
  for (int c=0;c<NCON; c++){     
    delete FF[c];    
  }
  return;
}

Eigen::MatrixXf RegionalPolePlacement( Eigen::MatrixXf A,  Eigen::MatrixXf B, float amax, float amin, float beta)  {
  int NUNK=2;
  int NCON=4;
  MatrixXf* unk[NUNK];  // vector of pointers to matrix unknowns 
  float **unkpoint1, **unkpoint2;  // pointers to values of matrix unknowns (2 for symmetric matrices)
  char sym[NUNK]={1,0};            // symmetric matrix indicator (1 yes, 0 no)
  int start[NUNK];                 // start position of each matrix unknown (MAYBE NOT USED ANYMORE)
  int sizeunk; // size of unknown vector (first unknown is lambda)
  int counter;
  int MATRIXN;
  int MATRIXM;
  float lambda=0.0;

  MATRIXN=A.cols();
  MATRIXM=B.cols();

  int sizecon[NCON]={MATRIXN,MATRIXN,2*MATRIXN,MATRIXN}; // size of linear matrix inequality constraints

  unk[0]=new MatrixXf(MATRIXN,MATRIXN);
  unk[1]=new MatrixXf(MATRIXM,MATRIXN);
  
  start[0]=1; // first unknown is lambda
  for(int i=1;i<NUNK;i++) {
    if (sym[i-1]){
      start[i]= start[i-1]+unk[i-1]->rows()*(unk[i-1]->rows()+1)/2;
    }
    else {
      start[i]=start[i-1]+unk[i-1]->cols()*unk[i-1]->rows();
    }
  }
  if (sym[NUNK-1]) {
    sizeunk=start[NUNK-1]+unk[NUNK-1]->rows()*(unk[NUNK-1]->rows()+1)/2;
  }
  else {
    sizeunk=start[NUNK-1]+unk[NUNK-1]->cols()*unk[NUNK-1]->rows();
  }   
  
  unkpoint1 = new float*[sizeunk];  
  unkpoint2 = new float*[sizeunk];  
  counter=1;

  for (int i=0; i<NUNK;i++) { 
    for (int lin=0; lin<unk[i]->rows(); lin++ ){
      if (sym[i]) {
        for (int col=0; col<=lin; col++ ){
           unkpoint1[counter]=&(unk[i]->coeffRef(lin,col));          
           unkpoint2[counter]=&(unk[i]->coeffRef(col,lin));  
           *(unkpoint1[counter])=0.0;
           *(unkpoint2[counter])=0.0;
           counter++;
        }
      } 
      else {
        for (int col=0; col<unk[i]->cols(); col++ ){
          unkpoint1[counter]=&(unk[i]->coeffRef(lin,col));  
          unkpoint2[counter]=&(unk[i]->coeffRef(lin,col));
          *(unkpoint1[counter])=0.0; 
          counter++;             
        }
      }
    }
  }

  unkpoint1[0]=&lambda;
  unkpoint2[0]=&lambda;

  VectorXf x(sizeunk); // vector of all unknowns 
  MatrixXf* F[NCON*(sizeunk+1)];
  MatrixXf* V[NCON*(sizeunk+1)];
 
 
  x.setZero();
  x(0)=1.0; //initialisation of lambda
  
  for(int i=sizeunk; i>=0;i--){

    //i = sizeunk is the constant term F0
    //i = 0 is lambda term
    for (int j=0;j<NCON;j++) {

      V[j*(sizeunk+1)+i]=new MatrixXf(sizecon[j],sizecon[j]);
      F[j*(sizeunk+1)+i]=new MatrixXf(sizecon[j],sizecon[j]);
    }
    
    if (i<sizeunk){
      *(unkpoint1[i])=1.0;
      *(unkpoint2[i])=1.0;  
    }
    
    ///////////////////////////////// WRITE LMIS HERE /////////////////////////////////    
    MatrixXf X=*(unk[0]);   // this is suboptimal, can be skipped! just replace directly the unknowns
    MatrixXf Y=*(unk[1]);    
    MatrixXf I(MATRIXN,MATRIXN);    
    I.setIdentity();
    V[0*(sizeunk+1)+i]->noalias()=X-0.1*I;
    V[1*(sizeunk+1)+i]->noalias()=-B*Y-Y.transpose()*B.transpose()-2*amin*X-X*A.transpose()-A*X; 
    V[2*(sizeunk+1)+i]->block(0,0,MATRIXN,MATRIXN)=-beta*(X*A.transpose()+A*X+Y.transpose()*B.transpose()+B*Y); //
    V[2*(sizeunk+1)+i]->block(0,MATRIXN,MATRIXN,MATRIXN)=-A*X+X*A.transpose()-B*Y+Y.transpose()*B.transpose(); //
    V[2*(sizeunk+1)+i]->block(MATRIXN,MATRIXN,MATRIXN,MATRIXN)= -beta*(X*A.transpose()+A*X+Y.transpose()*B.transpose()+B*Y); //
    V[2*(sizeunk+1)+i]->block(MATRIXN,0,MATRIXN,MATRIXN)=A*X-X*A.transpose()+B*Y-Y.transpose()*B.transpose(); //
    V[3*(sizeunk+1)+i]->noalias()=B*Y+Y.transpose()*B.transpose()+2*amax*X+X*A.transpose()+A*X; //
    
    /////////////////////////////// END WRITE LMIS HERE ///////////////////////////////    
    
    

    if (i<sizeunk){
      *(unkpoint1[i])=0.0;
      *(unkpoint2[i])=0.0;  
      for (int j=0;j<NCON;j++) {
        F[j*(sizeunk+1)+i]->noalias()=*V[j*(sizeunk+1)+i]-*V[j*(sizeunk+1)+sizeunk];
      }
    } else {
      for (int j=0;j<NCON;j++) {  // this is executed first and it initialises F
        F[j*(sizeunk+1)+i]->noalias()=*V[j*(sizeunk+1)+i];
 //       es.compute(*F[j][i], false);              /// These three lines of code force the compiler to include
 //       if (-es.eigenvalues().minCoeff()>x(0)) {  /// eigenvalue computations in the sketch.  Remove to save
  //        x(0)=-es.eigenvalues().minCoeff();      /// memory (but make sure initialisation of [0] is ok!)
 //       }
      }
    }
  }

  for (int j=0;j<NCON;j++) {
    F[j*(sizeunk+1)+0]->setIdentity();
  }

  ////// solver starts here!!!!! //////
  LMISolver(F, V,  x,  sizeunk, sizecon,  NUNK,  NCON);

  for (int i=1;i<sizeunk;i++){
    *(unkpoint1[i])=x(i);
    *(unkpoint2[i])=x(i);  
  }
  
  MatrixXf K=(*unk[1])*(unk[0]->inverse());
  
  if (x(0)>=0) {// infeasible problem
    K.setZero();
  }

  delete unkpoint1;
  delete unkpoint2;
  for (int i=0;i<NUNK; i++){
    delete unk[i];
  }

  for (int i=0;i<sizeunk+1; i++){
    for (int c=0;c<NCON; c++){
      delete F[c*(sizeunk+1)+i];    
      delete V[c*(sizeunk+1)+i];  
    }
  }
  return K;
}


Eigen::MatrixXf RegionalPolePlacementRobust1Param( Eigen::MatrixXf A0,  Eigen::MatrixXf A1, Eigen::MatrixXf B, float amax, float amin, float beta, float pmax, float pmin)  {
  int NUNK=2;
  int NCON=7;
  MatrixXf* unk[NUNK];  // vector of pointers to matrix unknowns 
  float **unkpoint1, **unkpoint2;  // pointers to values of matrix unknowns (2 for symmetric matrices)
  char sym[NUNK]={1,0};            // symmetric matrix indicator (1 yes, 0 no)
  int start[NUNK];                 // start position of each matrix unknown (MAYBE NOT USED ANYMORE)
  int sizeunk; // size of unknown vector (first unknown is lambda)
  int counter;
  int MATRIXN;
  int MATRIXM;
  float lambda=0.0;

  MATRIXN=A0.cols();
  MATRIXM=B.cols();

  int sizecon[NCON]={MATRIXN,MATRIXN,2*MATRIXN,MATRIXN,MATRIXN,2*MATRIXN,MATRIXN}; // size of linear matrix inequality constraints

  unk[0]=new MatrixXf(MATRIXN,MATRIXN);
  unk[1]=new MatrixXf(MATRIXM,MATRIXN);
  
  start[0]=1; // first unknown is lambda
  for(int i=1;i<NUNK;i++) {
    if (sym[i-1]){
      start[i]= start[i-1]+unk[i-1]->rows()*(unk[i-1]->rows()+1)/2;
    }
    else {
      start[i]=start[i-1]+unk[i-1]->cols()*unk[i-1]->rows();
    }
  }
  if (sym[NUNK-1]) {
    sizeunk=start[NUNK-1]+unk[NUNK-1]->rows()*(unk[NUNK-1]->rows()+1)/2;
  }
  else {
    sizeunk=start[NUNK-1]+unk[NUNK-1]->cols()*unk[NUNK-1]->rows();
  }   
  
  unkpoint1 = new float*[sizeunk];  
  unkpoint2 = new float*[sizeunk];  
  counter=1;

  for (int i=0; i<NUNK;i++) { 
    for (int lin=0; lin<unk[i]->rows(); lin++ ){
      if (sym[i]) {
        for (int col=0; col<=lin; col++ ){
           unkpoint1[counter]=&(unk[i]->coeffRef(lin,col));          
           unkpoint2[counter]=&(unk[i]->coeffRef(col,lin));  
           *(unkpoint1[counter])=0.0;
           *(unkpoint2[counter])=0.0;
           counter++;
        }
      } 
      else {
        for (int col=0; col<unk[i]->cols(); col++ ){
          unkpoint1[counter]=&(unk[i]->coeffRef(lin,col));  
          unkpoint2[counter]=&(unk[i]->coeffRef(lin,col));
          *(unkpoint1[counter])=0.0; 
          counter++;             
        }
      }
    }
  }

  unkpoint1[0]=&lambda;
  unkpoint2[0]=&lambda;

  VectorXf x(sizeunk); // vector of all unknowns 
  MatrixXf* F[NCON*(sizeunk+1)];
  MatrixXf* V[NCON*(sizeunk+1)];
 
 
  x.setZero();
  x(0)=1.0; //initialisation of lambda
  
  for(int i=sizeunk; i>=0;i--){

    //i = sizeunk is the constant term F0
    //i = 0 is lambda term
    for (int j=0;j<NCON;j++) {

      V[j*(sizeunk+1)+i]=new MatrixXf(sizecon[j],sizecon[j]);
      F[j*(sizeunk+1)+i]=new MatrixXf(sizecon[j],sizecon[j]);
    }
    
    if (i<sizeunk){
      *(unkpoint1[i])=1.0;
      *(unkpoint2[i])=1.0;  
    }
    
    ///////////////////////////////// WRITE LMIS HERE /////////////////////////////////    
    MatrixXf X=*(unk[0]);   // this is suboptimal, can be skipped! just replace directly the unknowns
    MatrixXf Y=*(unk[1]);    
    MatrixXf I(MATRIXN,MATRIXN);    
    I.setIdentity();
    V[0*(sizeunk+1)+i]->noalias()=X-0.1*I;
    MatrixXf A = A0+pmax*A1;
    V[1*(sizeunk+1)+i]->noalias()=-B*Y-Y.transpose()*B.transpose()-2*amin*X-X*A.transpose()-A*X; 
    V[2*(sizeunk+1)+i]->block(0,0,MATRIXN,MATRIXN)=-beta*(X*A.transpose()+A*X+Y.transpose()*B.transpose()+B*Y); //
    V[2*(sizeunk+1)+i]->block(0,MATRIXN,MATRIXN,MATRIXN)=-A*X+X*A.transpose()-B*Y+Y.transpose()*B.transpose(); //
    V[2*(sizeunk+1)+i]->block(MATRIXN,MATRIXN,MATRIXN,MATRIXN)= -beta*(X*A.transpose()+A*X+Y.transpose()*B.transpose()+B*Y); //
    V[2*(sizeunk+1)+i]->block(MATRIXN,0,MATRIXN,MATRIXN)=A*X-X*A.transpose()+B*Y-Y.transpose()*B.transpose(); //
    V[3*(sizeunk+1)+i]->noalias()=B*Y+Y.transpose()*B.transpose()+2*amax*X+X*A.transpose()+A*X; //
    A = A0+pmin*A1;
    V[4*(sizeunk+1)+i]->noalias()=-B*Y-Y.transpose()*B.transpose()-2*amin*X-X*A.transpose()-A*X; 
    V[5*(sizeunk+1)+i]->block(0,0,MATRIXN,MATRIXN)=-beta*(X*A.transpose()+A*X+Y.transpose()*B.transpose()+B*Y); //
    V[5*(sizeunk+1)+i]->block(0,MATRIXN,MATRIXN,MATRIXN)=-A*X+X*A.transpose()-B*Y+Y.transpose()*B.transpose(); //
    V[5*(sizeunk+1)+i]->block(MATRIXN,MATRIXN,MATRIXN,MATRIXN)= -beta*(X*A.transpose()+A*X+Y.transpose()*B.transpose()+B*Y); //
    V[5*(sizeunk+1)+i]->block(MATRIXN,0,MATRIXN,MATRIXN)=A*X-X*A.transpose()+B*Y-Y.transpose()*B.transpose(); //
    V[6*(sizeunk+1)+i]->noalias()=B*Y+Y.transpose()*B.transpose()+2*amax*X+X*A.transpose()+A*X; // 

    /////////////////////////////// END WRITE LMIS HERE ///////////////////////////////    
    
    

    if (i<sizeunk){
      *(unkpoint1[i])=0.0;
      *(unkpoint2[i])=0.0;  
      for (int j=0;j<NCON;j++) {
        F[j*(sizeunk+1)+i]->noalias()=*V[j*(sizeunk+1)+i]-*V[j*(sizeunk+1)+sizeunk];
      }
    } else {
      for (int j=0;j<NCON;j++) {  // this is executed first and it initialises F
        F[j*(sizeunk+1)+i]->noalias()=*V[j*(sizeunk+1)+i];
 //       es.compute(*F[j][i], false);              /// These three lines of code force the compiler to include
 //       if (-es.eigenvalues().minCoeff()>x(0)) {  /// eigenvalue computations in the sketch.  Remove to save
  //        x(0)=-es.eigenvalues().minCoeff();      /// memory (but make sure initialisation of [0] is ok!)
 //       }
      }
    }
  }

  for (int j=0;j<NCON;j++) {
    F[j*(sizeunk+1)+0]->setIdentity();
  }

  ////// solver starts here!!!!! //////
  LMISolver(F, V,  x,  sizeunk, sizecon,  NUNK,  NCON);

  for (int i=1;i<sizeunk;i++){
    *(unkpoint1[i])=x(i);
    *(unkpoint2[i])=x(i);  
  }
  
  MatrixXf K=(*unk[1])*(unk[0]->inverse());
  
  if (x(0)>=0) {// infeasible problem
    K.setZero();
  }

  delete unkpoint1;
  delete unkpoint2;
  for (int i=0;i<NUNK; i++){
    delete unk[i];
  }

  for (int i=0;i<sizeunk+1; i++){
    for (int c=0;c<NCON; c++){
      delete F[c*(sizeunk+1)+i];    
      delete V[c*(sizeunk+1)+i];  
    }
  }
  return K;
}




Eigen::MatrixXf** RegionalPolePlacementGainScheduling1Param( Eigen::MatrixXf A0,  Eigen::MatrixXf A1, Eigen::MatrixXf B, float amax, float amin, float beta, float pmax, float pmin)  {
  int NUNK=3;
  int NCON=7;
  MatrixXf* unk[NUNK];  // vector of pointers to matrix unknowns 
  float **unkpoint1, **unkpoint2;  // pointers to values of matrix unknowns (2 for symmetric matrices)
  char sym[NUNK]={1,0,0};            // symmetric matrix indicator (1 yes, 0 no)
  int start[NUNK];                 // start position of each matrix unknown (MAYBE NOT USED ANYMORE)
  int sizeunk; // size of unknown vector (first unknown is lambda)
  int counter;
  int MATRIXN;
  int MATRIXM;
  float lambda=0.0;

  MATRIXN=A0.cols();
  MATRIXM=B.cols();

  int sizecon[NCON]={MATRIXN,MATRIXN,2*MATRIXN,MATRIXN,MATRIXN,2*MATRIXN,MATRIXN}; // size of linear matrix inequality constraints

  unk[0]=new MatrixXf(MATRIXN,MATRIXN);
  unk[1]=new MatrixXf(MATRIXM,MATRIXN);
  unk[2]=new MatrixXf(MATRIXM,MATRIXN);
  
  start[0]=1; // first unknown is lambda
  for(int i=1;i<NUNK;i++) {
    if (sym[i-1]){
      start[i]= start[i-1]+unk[i-1]->rows()*(unk[i-1]->rows()+1)/2;
    }
    else {
      start[i]=start[i-1]+unk[i-1]->cols()*unk[i-1]->rows();
    }
  }
  if (sym[NUNK-1]) {
    sizeunk=start[NUNK-1]+unk[NUNK-1]->rows()*(unk[NUNK-1]->rows()+1)/2;
  }
  else {
    sizeunk=start[NUNK-1]+unk[NUNK-1]->cols()*unk[NUNK-1]->rows();
  }   
  
  unkpoint1 = new float*[sizeunk];  
  unkpoint2 = new float*[sizeunk];  
  counter=1;

  for (int i=0; i<NUNK;i++) { 
    for (int lin=0; lin<unk[i]->rows(); lin++ ){
      if (sym[i]) {
        for (int col=0; col<=lin; col++ ){
           unkpoint1[counter]=&(unk[i]->coeffRef(lin,col));          
           unkpoint2[counter]=&(unk[i]->coeffRef(col,lin));  
           *(unkpoint1[counter])=0.0;
           *(unkpoint2[counter])=0.0;
           counter++;
        }
      } 
      else {
        for (int col=0; col<unk[i]->cols(); col++ ){
          unkpoint1[counter]=&(unk[i]->coeffRef(lin,col));  
          unkpoint2[counter]=&(unk[i]->coeffRef(lin,col));
          *(unkpoint1[counter])=0.0; 
          counter++;             
        }
      }
    }
  }

  unkpoint1[0]=&lambda;
  unkpoint2[0]=&lambda;

  VectorXf x(sizeunk); // vector of all unknowns 
  MatrixXf* F[NCON*(sizeunk+1)];
  MatrixXf* V[NCON*(sizeunk+1)];
 
 
  x.setZero();
  x(0)=1.0; //initialisation of lambda
  
  for(int i=sizeunk; i>=0;i--){

    //i = sizeunk is the constant term F0
    //i = 0 is lambda term
    for (int j=0;j<NCON;j++) {

      V[j*(sizeunk+1)+i]=new MatrixXf(sizecon[j],sizecon[j]);
      F[j*(sizeunk+1)+i]=new MatrixXf(sizecon[j],sizecon[j]);
    }
    
    if (i<sizeunk){
      *(unkpoint1[i])=1.0;
      *(unkpoint2[i])=1.0;  
    }
    
    ///////////////////////////////// WRITE LMIS HERE /////////////////////////////////    
    MatrixXf X=*(unk[0]);   // this is suboptimal, can be skipped! just replace directly the unknowns
    MatrixXf Y0=*(unk[1]);    
    MatrixXf Y1=*(unk[2]);   
    MatrixXf I(MATRIXN,MATRIXN);    
    I.setIdentity();
    V[0*(sizeunk+1)+i]->noalias()=X-0.1*I;
    MatrixXf A = A0+pmax*A1;
    MatrixXf Y = Y0+pmax*Y1;
    V[1*(sizeunk+1)+i]->noalias()=-B*Y-Y.transpose()*B.transpose()-2*amin*X-X*A.transpose()-A*X; 
    V[2*(sizeunk+1)+i]->block(0,0,MATRIXN,MATRIXN)=-beta*(X*A.transpose()+A*X+Y.transpose()*B.transpose()+B*Y); //
    V[2*(sizeunk+1)+i]->block(0,MATRIXN,MATRIXN,MATRIXN)=-A*X+X*A.transpose()-B*Y+Y.transpose()*B.transpose(); //
    V[2*(sizeunk+1)+i]->block(MATRIXN,MATRIXN,MATRIXN,MATRIXN)= -beta*(X*A.transpose()+A*X+Y.transpose()*B.transpose()+B*Y); //
    V[2*(sizeunk+1)+i]->block(MATRIXN,0,MATRIXN,MATRIXN)=A*X-X*A.transpose()+B*Y-Y.transpose()*B.transpose(); //
    V[3*(sizeunk+1)+i]->noalias()=B*Y+Y.transpose()*B.transpose()+2*amax*X+X*A.transpose()+A*X; //
    A = A0+pmin*A1;
    Y = Y0+pmin*Y1;
    V[4*(sizeunk+1)+i]->noalias()=-B*Y-Y.transpose()*B.transpose()-2*amin*X-X*A.transpose()-A*X; 
    V[5*(sizeunk+1)+i]->block(0,0,MATRIXN,MATRIXN)=-beta*(X*A.transpose()+A*X+Y.transpose()*B.transpose()+B*Y); //
    V[5*(sizeunk+1)+i]->block(0,MATRIXN,MATRIXN,MATRIXN)=-A*X+X*A.transpose()-B*Y+Y.transpose()*B.transpose(); //
    V[5*(sizeunk+1)+i]->block(MATRIXN,MATRIXN,MATRIXN,MATRIXN)= -beta*(X*A.transpose()+A*X+Y.transpose()*B.transpose()+B*Y); //
    V[5*(sizeunk+1)+i]->block(MATRIXN,0,MATRIXN,MATRIXN)=A*X-X*A.transpose()+B*Y-Y.transpose()*B.transpose(); //
    V[6*(sizeunk+1)+i]->noalias()=B*Y+Y.transpose()*B.transpose()+2*amax*X+X*A.transpose()+A*X; // 

    /////////////////////////////// END WRITE LMIS HERE ///////////////////////////////    
    
    

    if (i<sizeunk){
      *(unkpoint1[i])=0.0;
      *(unkpoint2[i])=0.0;  
      for (int j=0;j<NCON;j++) {
        F[j*(sizeunk+1)+i]->noalias()=*V[j*(sizeunk+1)+i]-*V[j*(sizeunk+1)+sizeunk];
      }
    } else {
      for (int j=0;j<NCON;j++) {  // this is executed first and it initialises F
        F[j*(sizeunk+1)+i]->noalias()=*V[j*(sizeunk+1)+i];
 //       es.compute(*F[j][i], false);              /// These three lines of code force the compiler to include
 //       if (-es.eigenvalues().minCoeff()>x(0)) {  /// eigenvalue computations in the sketch.  Remove to save
  //        x(0)=-es.eigenvalues().minCoeff();      /// memory (but make sure initialisation of [0] is ok!)
 //       }
      }
    }
  }

  for (int j=0;j<NCON;j++) {
    F[j*(sizeunk+1)+0]->setIdentity();
  }

  ////// solver starts here!!!!! //////
  LMISolver(F, V,  x,  sizeunk, sizecon,  NUNK,  NCON);

  for (int i=1;i<sizeunk;i++){
    *(unkpoint1[i])=x(i);
    *(unkpoint2[i])=x(i);  
  }
  
  MatrixXf** K;
  K = new MatrixXf*[2];
  K[0] = new MatrixXf(MATRIXM,MATRIXN);
  K[1] = new MatrixXf(MATRIXM,MATRIXN); 
  
  *K[0] = (*unk[1])*(unk[0]->inverse());
  *K[1] = (*unk[2])*(unk[0]->inverse());

  if (x(0)>=0) {// infeasible problem
    K[0]->setZero();
    K[1]->setZero();
  }

  delete unkpoint1;
  delete unkpoint2;
  for (int i=0;i<NUNK; i++){
    delete unk[i];
  }

  for (int i=0;i<sizeunk+1; i++){
    for (int c=0;c<NCON; c++){
      delete F[c*(sizeunk+1)+i];    
      delete V[c*(sizeunk+1)+i];  
    }
  }
  return K;
}




