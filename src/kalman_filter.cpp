#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_ ; 
  P_ = F_ * P_ * F_.transpose() + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {

  VectorXd z_predicted = H_ * x_;
  VectorXd y = z - z_predicted;
  UpdateStep(y);
  
    
}
void KalmanFilter::UpdateEKF(const VectorXd &z) {

  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  float rho = sqrt(px*px + py*py);
  float theta = atan2(py, px);
  if(rho==0)
  {
  	return;
  }
  float rho_dot = (px*vx + py*vy) / rho;
  VectorXd H_of_x = VectorXd(3);
  H_of_x << rho, theta, rho_dot;
  VectorXd y = z - H_of_x;
  //Making sure that the heading is within pi to -pi
  while ( y(1) > M_PI || y(1) < -M_PI ) {
    if ( y(1) > M_PI ) {
      y(1) -= M_PI;
    } else {
      y(1) += M_PI;
    }
  }
  UpdateStep(y);
}

void KalmanFilter::UpdateStep(const VectorXd &y){

  //Common update step for kalman filter and extended Kalman Filter
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}
