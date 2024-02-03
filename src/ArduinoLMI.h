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

void LMISolver(Eigen::MatrixXf *F[], Eigen::MatrixXf *V[], Eigen::VectorXf& x, int sizeunk, int *sizecon, int NUNK, int NCON);
Eigen::MatrixXf RegionalPolePlacement(Eigen::MatrixXf& A, Eigen::MatrixXf& B, float amax, float amin, float beta);  
Eigen::MatrixXf RegionalPolePlacementRobust1Param(Eigen::MatrixXf& A0, Eigen::MatrixXf& A1, Eigen::MatrixXf& B, float amax, float amin, float beta, float p1max, float p1min);  
Eigen::MatrixXf** RegionalPolePlacementGainScheduling1Param(Eigen::MatrixXf& A0,  Eigen::MatrixXf& A1, Eigen::MatrixXf& B, float amax, float amin, float beta, float pmax, float pmin);
Eigen::MatrixXf DecayRate(Eigen::MatrixXf& A, Eigen::MatrixXf& B, float amin);
Eigen::MatrixXf DecayRateRobust1Param(Eigen::MatrixXf& A0, Eigen::MatrixXf& A1, Eigen::MatrixXf& B, float amin,float p1max, float p1min);  
Eigen::MatrixXf** DecayRateGainScheduling1Param(Eigen::MatrixXf& A0,  Eigen::MatrixXf& A1, Eigen::MatrixXf& B, float amin, float pmax, float pmin);
Eigen::MatrixXf DecayRateRobustNormBound(Eigen::MatrixXf& A0, Eigen::MatrixXf& B, Eigen::MatrixXf& E, Eigen::MatrixXf& F, float amin);  




