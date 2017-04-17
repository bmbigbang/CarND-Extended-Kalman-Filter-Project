#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
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
  R_laser_ << 0.001, 0,
              0, 0.001;

  //measurement covariance matrix - radar
  R_radar_ << 100, 0, 0,
              0, 1000, 0,
              0, 0, 70;

  // process noise covariance matrix
  ekf_.Q_ = MatrixXd(4, 4);

  // initialize state transition matrix
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  // initialize values for the object covariance matrix
  ekf_.P_ = MatrixXd(4, 4);

  ekf_.I_ = MatrixXd::Identity(4, 4);

  // initialize place holders for kalman filter calculation
  MatrixXd S_;
  S_ = MatrixXd(4, 4);

  MatrixXd K_;
  K_ = MatrixXd(4, 4);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.m_ = VectorXd(3);
    ekf_.P_ << 1, 0, 2, 0,
               0, 1, 0, 2,
               2, 0, 100000, 0,
               0, 2, 0, 100000;
    // set the linear laser measurement matrix
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      //set the state with the initial location and velocity

      double ro = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double ro_dot = measurement_pack.raw_measurements_[2];
      ekf_.x_ << ro * cos(phi), ro * sin(phi), ro_dot * cos(phi), ro_dot * sin(phi);

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //set the state with the initial location and zero velocity
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
   */
  //dt - expressed in seconds
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  float noise_ax = 9;
  float noise_ay = 9;
  ekf_.Q_ << (noise_ax * pow(dt, 4)) / 4 , 0, (noise_ax * pow(dt, 3)) / 2, 0,
              0, (noise_ay * pow(dt, 4)) / 4, 0, (noise_ay * pow(dt, 3)) / 2,
              (noise_ax * pow(dt, 3)) / 2, 0, (noise_ax * pow(dt, 2)), 0,
              0, (noise_ay * pow(dt, 3)) / 2, 0, (noise_ay * pow(dt, 2));
  // avoid predicting twice if the measurements coincide
  if ( dt > 0.001 )
  {
    ekf_.Predict();
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    float px = ekf_.x_[0];
    float py = ekf_.x_[1];
    float vx = ekf_.x_[2];
    float vy = ekf_.x_[3];

    //check division by zero
    if (px <= 0.00001) {
      px = 0.00001;
    }
    if (py <= 0.00001) {
      py = 0.00001;
    }

    float rho = sqrt(px * px + py * py);
    float phi = atan(py / px);
    float rho_dot = (px * vx + py * vy) / rho;

    ekf_.m_ << rho, phi, rho_dot;

    ekf_.UpdateEKF(measurement_pack.raw_measurements_, ekf_.m_ , Hj_, R_radar_);
  } else {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_, H_laser_, R_laser_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
