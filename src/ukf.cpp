#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

	
  if (is_initialized_){
  	/**
  	TODO:

  	Complete this function! Make sure you switch between lidar and radar
  	measurements.
  	*/
	  
	if( meas_package.sensor_type_ == MeasurementPackage::RADAR)
  	{
  		printf("RADAR \n");
  		use_laser_ = false;
  		use_radar_ = true;
	}	
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
	{
  		printf("LASER\n");
  		use_laser_ = true;
  		use_radar_ = false;
  	}
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3- n_x_;

  // Create sigma points
  MatrixXd Xsig_ = MatrixXd(n_x_, 2 * n_x_ + 1);
  
  // Calculate square root of P
  MatrixXd A_ = P_.llt().matrixL(); 
  
  VectorXd x_k = x_;
  Xsig_.col(0) = x_k;
  int p1, p2 = 0;
  for (int i=0; i< n_x_; i++)
  {
  	p1 = i + 1;
  	p2 = i + 5;
  	Xsig_.col(p1) = x_ + A_.col(i) * sqrt(lambda_ - n_x_);
  	Xsig_.col(p2) = x_ - A_.col(i) * sqrt(lambda_ - n_x_);
  	//printf("X_sig(%d): %f",p1,Xsig_.col(p1));
  	//printf("X_sig(%d): %f",p2,Xsig_.col(p2));
  }
  
  // Create augumented mean vector
  VectorXd x_aug_ = VectorXd(n_aug_);
  
  // Create augmented covariance matrix
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
  
  // Create sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  // Create augmented mean state
  x_aug_.head(5) = x_;
  x_aug_(6) = 0;
  Xsig_aug_.col(0) = x_aug_;
  
  // Create augumeted covariance matirx
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_,n_x_) = P_;
  P_aug_(n_x_+1,n_x_+1) = std_a_ * std_a_;
  P_aug_(n_x_+2,n_x_+2) = std_yawdd_ * std_yawdd_ ;
  
  // Create square root of augmented covarance matrix
  MatrixXd A_aug_ = P_aug_.llt().matrixL();
  
  // Create augmented sigma points
  Xsig_aug_.col(0) = x_aug_;
  for (int i=0 ; i< 2 * n_aug_ + 1 ; i++)
  {
  	p1 = i + 1;
  	p2 = i + 5;
  	Xsig_aug_.col(p1) = x_aug_ + A_aug_.col(i) * sqrt(lambda_ - n_x_);
  	Xsig_aug_.col(p2) = x_aug_ - A_aug_.col(i) * sqrt(lambda_ - n_x_);
  	
  }
  
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  
  // Predict sigma points
  // Avoid division by zero
  // write predicted sigma points into the right column
  for (int i=0; i < 2 * n_x_+1 ; i++)
  {
  	double px_ = Xsig_aug_(0,i);
  	double py_ = Xsig_aug_(1,i);
  	double vk_ = Xsig_aug_(2,i);
  	double yaw_ = Xsig_aug_(3,i);
  	double yawd_ = Xsig_aug_(4,i);
  	double n_a_ = Xsig_aug_(5,i);
  	double n_yawdd_ = Xsig_aug_(6,i);
  	
  	double px_p, py_p=0;
  	
  	if (fabs(yawd_) > 0.001)
  	{
  		px_p = px_ + vk_ / yawd_ * (sin((yaw_ + (yawd_ * delta_t)) - sin(yaw_)));
  		py_p = py_ + vk_ / yawd_ * (cos(yaw_) - cos(yaw_ + (yawd_ * delta_t)));
  	}
  	else
  	{
  		px_p = px_ + vk_ * cos(yaw_) * delta_t;
  		py_p = py_ + vk_ * sin(yaw_) * delta_t;  	
  	}
  	
  	// Add noise
  	px_p = px_p + 0.5 * delta_t * delta_t * cos(yaw_) * n_a_;
  	py_p = py_p + 0.5 * delta_t * delta_t * sin(yaw_) * n_yawdd_;
  	
  	double v_p_ = vk_;
  	v_p_ = v_p_ + delta_t * n_a_;
  	
  	double yaw_p_ = yaw_p_;
  	yaw_p_ = yaw_p_ + yawd_ * delta_t + 0.5 * delta_t * delta_t * n_yawdd_;
  	
  	double yawd_p_ = yawd_;
  	yawd_p_ = yawd_p_ + delta_t * n_yawdd_; 
  	
  	//Write the predicted sigma points into right column
  	Xsig_pred_(0,i) = px_p;
  	Xsig_pred_(1,i) = py_p;
  	Xsig_pred_(2,i) = v_p_;
  	Xsig_pred_(3,i) = yaw_p_;
  	Xsig_pred_(4,i) = yawd_p_; 
  }
  
  VectorXd weights_ = VectorXd(2*n_aug_+1);
  // Set weights
  
  weights_(0) =  lambda_ / (lambda_ + n_aug_);
  for( int i=0; i < 2 * n_aug_ + 1; i++ )
  {
  	//Iterate over sigma points
  	weights_(i) = 0.5 / (lambda_ + n_aug_);   
  } 
  
  //predict state mean
  x_.fill(0.0);
  for( int i=0; i < 2 * n_aug_ + 1; i++ )
  {
  	//Iterate over sigma points
  	x_ = x_ + weights_(i) * Xsig_pred_.col(i) ;   
  }
  
  //predict covariance matrix
  P_.fill(0.0);
  for( int i=0; i < 2 * n_aug_ + 1; i++ )
  {
  	//Iterate over sigma points
  	
  	VectorXd x_diff_ =  Xsig_pred_.col(i) - x_; 	
  	
  	// Angel normalization
  	while (x_diff_(3) > M_PI) x_diff_(3) -= 2.*M_PI;
  	while (x_diff_(3) < M_PI) x_diff_(3) += 2.*M_PI;
  	P_ = P_ + weights_(i) * x_diff_ * x_diff_.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
