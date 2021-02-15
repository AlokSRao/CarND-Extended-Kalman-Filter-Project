#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
	//Initialising P matrix
  	ekf_.P_ = MatrixXd(4, 4);
  	ekf_.P_ <<  1,0,0,0,
  				0,1,0,0,
  				0,0,1000,
  				0,0,0,0,1000;
  
  H_laser_  << 1,0,0,0,
  				0,1,0,0;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {


    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      
      cout << "EKF : RADAR measurement " << endl;
      float rho = measurement_pack.raw_measurements_[0]; // radius
  	  float phi = measurement_pack.raw_measurements_[1]; // angle with x axis
  	  float rho_dot = measurement_pack.raw_measurements_[2]; // velocity in polar
      
      float px = rho * cos(phi);
      float py = rho * sin(phi);
      float vx = rho_dot * cos(phi);
      float vy = rho_dot * sin(phi);
      
      //avoid any divide by zero in future
      if ( px < 0.0001 ) {
        px = 0.0001;
      }
  	  
      if ( py < 0.0001 ) {
        py = 0.0001;
      }
      
      ekf_.x_ << px,py,vx,vy;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      cout << "EKF : LASER measurement " << endl;
      float px = measurement_pack.raw_measurements_[0]; 
  	  float py = measurement_pack.raw_measurements_[1];
      
      ekf_.x_ << px,py,0,0;
    }
	// Saving first timestamp in seconds
    previous_timestamp_ = measurement_pack.timestamp_ ;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

   float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
   previous_timestamp_ = measurement_pack.timestamp_;

    // State transition matrix update
   ekf_.F_ = MatrixXd(4, 4);
   ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;
  
  // Noise covariance matrix computation
  // Noise values from the task
  float noise_ax = 15.0;
  float noise_ay = 15.0;

  float dt_2 = dt * dt; //dt^2
  float dt_3 = dt_2 * dt; //dt^3
  float dt_4 = dt_3 * dt; //dt^4
  float dt_4_4 = dt_4 / 4; //dt^4/4
  float dt_3_2 = dt_3 / 2; //dt^3/2
  
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4_4 * noise_ax, 0, dt_3_2 * noise_ax, 0,
	         0, dt_4_4 * noise_ay, 0, dt_3_2 * noise_ay,
	         dt_3_2 * noise_ax, 0, dt_2 * noise_ax, 0,
 	         0, dt_3_2 * noise_ay, 0, dt_2 * noise_ay;

  
  ekf_.Predict();

  /**
   * Update
   */

 
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
  	ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    ekf_.H_ = H_laser_;
  	ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
  

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
