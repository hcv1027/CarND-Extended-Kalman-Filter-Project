#include "tools.h"
#include <cmath>
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

const float Tools::DOUBLE_PI = 2 * 3.14159265358979323846;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i = 0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];

    // coefficient-wise multiplication
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse / estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) {
  /**
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3, 4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // check division by zero
  float d1 = px * px + py * py;
  if (std::fabs(d1) < 0.0001) {
    std::cout << "CalculateJacobian() - Error - Division by Zero" << std::endl;
    return Hj;
  }
  float d2 = std::pow(d1, 0.5);
  float d3 = std::pow(d1, 1.5);

  // compute the Jacobian matrix
  Hj(0, 0) = px / d2;
  Hj(0, 1) = py / d2;
  Hj(0, 2) = 0;
  Hj(0, 3) = 0;
  Hj(1, 0) = -py / d1;
  Hj(1, 1) = px / d1;
  Hj(1, 2) = 0;
  Hj(1, 3) = 0;
  Hj(2, 0) = py * (vx * py - vy * px) / d3;
  Hj(2, 1) = px * (vy * px - vx * py) / d3;
  Hj(2, 2) = px / d2;
  Hj(2, 3) = py / d2;

  return Hj;
}

float Tools::AdjustRadian(float radian) {
  if (std::fabs(radian) > DOUBLE_PI && radian > 0.0) {
    do {
      radian -= DOUBLE_PI;
    } while (std::fabs(radian) > DOUBLE_PI);
  } else if (std::fabs(radian) > DOUBLE_PI && radian < 0.0) {
    do {
      radian += DOUBLE_PI;
    } while (std::fabs(radian) > DOUBLE_PI);
  }
  return radian;
}
