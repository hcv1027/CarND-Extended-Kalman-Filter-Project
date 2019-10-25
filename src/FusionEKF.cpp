#include "FusionEKF.h"
#include <cmath>
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

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0, 0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0, 0, 0.0009, 0, 0, 0, 0.09;

  /**
   * Finish initializing the FusionEKF.
   * Set the process and measurement noises
   */
  noise_ax_ = 9.0;
  noise_ay_ = 9.0;
  ekf_.F_ = MatrixXd::Identity(4, 4);
  // ekf_.F_ << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
  ekf_.Q_ = MatrixXd::Zero(4, 4);
  // ekf_.Q_ << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  ekf_.H_ = MatrixXd(2, 4);
  ekf_.H_ << 1, 0, 0, 0, 0, 1, 0, 0;
  ekf_.R_laser_ = R_laser_;
  ekf_.R_radar_ = R_radar_;
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
    /**
     * Initialize the state ekf_.x_ with the first measurement.
     * Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state.
      float rho = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float rho_dot = measurement_pack.raw_measurements_(2);
      float sin_phi = std::sin(phi);
      float cos_phi = std::cos(phi);
      float px = cos_phi * rho;
      float py = sin_phi * rho;
      float vx = cos_phi * rho_dot;
      float vy = sin_phi * rho_dot;
      ekf_.x_ << px, py, vx, vy;
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_(0),
          measurement_pack.raw_measurements_(1), 0, 0;
    }
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1000, 0, 0, 0, 0, 1000;

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  float dt_2 = std::pow(dt, 2);
  float dt_3 = std::pow(dt, 3) / 2;
  float dt_4 = std::pow(dt, 4) / 4;
  ekf_.Q_(0, 0) = dt_4 * noise_ax_;
  ekf_.Q_(0, 2) = dt_3 * noise_ax_;
  ekf_.Q_(1, 1) = dt_4 * noise_ay_;
  ekf_.Q_(1, 3) = dt_3 * noise_ay_;
  ekf_.Q_(2, 0) = dt_3 * noise_ax_;
  ekf_.Q_(2, 2) = dt_2 * noise_ax_;
  ekf_.Q_(3, 1) = dt_3 * noise_ay_;
  ekf_.Q_(3, 3) = dt_2 * noise_ay_;
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
