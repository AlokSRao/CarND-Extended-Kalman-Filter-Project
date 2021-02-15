#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  if(estimations.size()==0)
  {
  	cout << "Invalid estimation data" << endl;
    return rmse;
  }
  else if(ground_truth.size()==0)
  {
    cout << "Invalid ground_truth data" << endl;
    return rmse;
  }
  else if(estimations.size() != ground_truth.size())
  {
    cout << "size mismatch between estimation and ground truth" << endl;
    return rmse;
  }
  
  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {

  VectorXd residual = estimations[i] - ground_truth[i];

  // coefficient-wise multiplication
  residual = residual.array()*residual.array();
  rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;
  
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Jacobian_H(3,4);

  if ( x_state.size() != 4 ) {
    cout << "ERROR - CalculateJacobian () - State vector wron dimensions." << endl;
    return Jacobian_H;
  }
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px*px+py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if(fabs(c1) < 0.0001){
		cout << "ERROR - CalculateJacobian () - Division by Zero" << endl;
		return Jacobian_H;
	}

	//compute the Jacobian matrix
	Jacobian_H << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Jacobian_H;
}
