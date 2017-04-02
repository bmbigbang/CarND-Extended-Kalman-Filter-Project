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
  P_ = F_ * P_ * F_.transpose();
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */

  I_ = MatrixXd::Identity(4, 4);
  S_ = MatrixXd(4, 4);
  K_ = MatrixXd(4, 4);

  // KF Measurement update step
  S_ = H_ * P_ * H_.transpose() + R_;
  K_ = P_ * H_.transpose() * S_.inverse();

  x_ = x_ + (K_ * (z - (H_ * x_)));
  P_ = (I_ - (K_ * H_)) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
}
