#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  /**
    * predict the state
  */
  // assume no external motion u for the position prediciton x_

  x_ = (F_ * x_);
  P_ = (F_ * P_ * F_.transpose()) + Q_;
}

void KalmanFilter::Update(const VectorXd &z, MatrixXd &H_laser_, MatrixXd &R_laser_) {
  /**
    * update the state by using Kalman Filter equations
  */
  MatrixXd Ht_ = H_laser_.transpose();
  // KF Measurement update step
  S_ = H_laser_ * P_ * Ht_ + R_laser_;
  K_ = P_ * Ht_ * S_.inverse();

  x_ = x_ + (K_ * (z - (H_laser_ * x_)));
  P_ = (I_ - (K_ * H_laser_)) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z, VectorXd &m_, MatrixXd &Hj_, MatrixXd &R_radar_) {
  /**
    * update the state by using Extended Kalman Filter equations
  */
  MatrixXd Ht_ = Hj_.transpose();
  S_ = Hj_ * P_ * Ht_ + R_radar_;
  K_ = P_ * Ht_ * S_.inverse();

  x_ = x_ + (K_ * (z - m_));
  P_ = (I_ - (K_ * Hj_)) * P_;
}
