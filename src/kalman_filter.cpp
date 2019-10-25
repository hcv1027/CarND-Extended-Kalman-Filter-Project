#include "kalman_filter.h"
#include <cmath>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_laser_in,
                        MatrixXd &R_radar_in, MatrixXd &Q_in) {
  x_ = x_in;                // 4x1
  P_ = P_in;                // 4x4
  F_ = F_in;                // 4x4
  H_ = H_in;                // 2x4
  R_laser_ = R_laser_in;    // 2x2
  R_radar_in = R_radar_in;  // 3x3
  Q_ = Q_in;                // 4x4
}

void KalmanFilter::Predict() {
  /**
   * Predict the state
   */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * Update the state by using Kalman Filter equations
   */
  VectorXd y = z - H_ * x_;
  MatrixXd H_T = H_.transpose();
  MatrixXd S = H_ * P_ * H_T + R_laser_;
  MatrixXd S_I = S.inverse();
  MatrixXd K = P_ * H_T * S_I;

  x_ = x_ + K * y;
  P_ = (MatrixXd::Identity(4, 4) - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * Update the state by using Extended Kalman Filter equations
   */
  auto h = [&]() {
    Eigen::VectorXd temp = Eigen::VectorXd(3);
    temp(0) = std::sqrt(std::pow(x_(0), 2) + std::pow(x_(1), 2));
    temp(1) = std::atan2(x_(1), x_(0));
    if (std::fabs(temp(0)) < 0.0001) {
      std::cout << "UpdateEKF() - Error - Division by Zero" << std::endl;
      temp(2) = 0.0;
    } else {
      temp(2) = (x_(0) * x_(2) + x_(1) * x_(3)) / temp(0);
    }
    return temp;
  };

  VectorXd y = z - h();
  // Check if y(1) is out of the range [-pi, pi]
  y(1) = tools.AdjustRadian(y(1));
  MatrixXd Hj = tools.CalculateJacobian(x_);  // 3x4
  MatrixXd H_T = Hj.transpose();              // 4x3
  MatrixXd S = Hj * P_ * H_T + R_radar_;      // 3x4 4x4 4x3 + 3x3 = 3x3
  MatrixXd S_I = S.inverse();                 // 3x3
  MatrixXd K = P_ * H_T * S_I;                // 4x4 4x4 3x3 = 4x3

  x_ = x_ + K * y;
  P_ = (MatrixXd::Identity(4, 4) - K * Hj) * P_;
}
