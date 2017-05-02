#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
int g_count=0;

inline double constrainAngle(double x)
{
     while (x > M_PI)  x -= 2.*M_PI;
     while (x < -M_PI) x += 2.*M_PI;
     
     /* Alternate normalization techniques */    
     /*x = fmod(x,2.*M_PI);
      if (x<0)
        x+=2.*M_PI;

     /* if (fabs(x) > M_PI)
      {
	x -= round(x/(2.0*M_PI)) * (2.0*M_PI);
      }*/
      
      //if(x > M_PI) x -= round(x/(2.*M_PI))*2.*M_PI;
      //if(x < -M_PI) x += round(x/(2.*M_PI))*2.*M_PI;
 
    return x;
}
/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
 
  n_x_ = 5;		// State dimension
  n_aug_ = 7;		// Augmented state dimension
  
  x_ = VectorXd(n_x_);	// initial state vector
  P_ = MatrixXd(n_x_, n_x_);	// initial covariance matrix

  std_a_ = 0.2811;  	// Process noise standard deviation longitudinal acceleration in m/s^2
  std_yawdd_ = 0.331;     // Process noise standard deviation yaw acceleration in rad/s^2
  std_laspx_ = 0.0599;	// Laser measurement noise standard deviation position1 in m
  std_laspy_ = 0.0537;   // Laser measurement noise standard deviation position2 in m
  std_radr_ = 1.28632;	// Radar measurement noise standard deviation radius in m
  std_radphi_ = 0.03;	// Radar measurement noise standard deviation angle in rad
  std_radrd_ = 0.3;	// Radar measurement noise standard deviation radius change in m/s
  use_laser_ = true;	// If this is false, laser measurements will ne ignored (except during init)
  use_radar_ = true;	// If this is false , radar measurements will be ignored (excet during init)
  
  is_initialized_ = false;	// initially set to false, set to true in first call of ProcessMeasurement
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  
  
  lambda_ = 3 - n_aug_;	// Sigma point spreading parameter
   
  sigma_points = 2 * n_aug_ + 1;

  Xsig_pred_ = MatrixXd(n_x_, sigma_points);	// Create predicted sigma points matrix
  
  n_z_radar_ = 3;	// Set measurement dimension, radar can measure r, phi and r_dot
  n_z_laser_ = 2;  	// Set measurment dimension, laser can measure x, y
  NIS_radar_ = 0.0;	// Initialize the current NIS for radar
  NIS_laser_ = 0.0;	// Initialize the current NIS for laser
  
  // Set weights
  weights_ = VectorXd(sigma_points);
  weights_(0) =  lambda_ / (lambda_ + n_aug_);
  for( int i=1; i < sigma_points; i++ )
  {
  	//Iterate over sigma points
  	weights_(i) = 0.5 / (lambda_ + n_aug_);   
  }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

	
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  	
  
  if( meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == false)
  {
   	return;
  }	
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == false)
  {
   	return;
  }
  
  
	  
  //first measurement
  	
  if (!is_initialized_){	

  	
  	//done initializing , no need to predict or update
	is_initialized_ = true;

	P_ << 1000.0, -1.0, 0.0, 1.0, -0.0020,
	      1.0, 5.9, 0.05, 0.5, 0.06,
	      -0.955, 0.22, 0.95, 0.0, 0.0,
	      -10.0, 0.05, 0.0, 0.1, 0.0,
	      -5.0, 0.0, 0.0, 0.0, 0.5;

	double rho = meas_package.raw_measurements_[0];
	double phi = meas_package.raw_measurements_[1];
	double yaw = 0;
	double px, py=0.0;
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
	{
		px = rho * cos(phi);
		py = rho * sin(phi);
		yaw = meas_package.raw_measurements_[2];
		
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
	{
		px = meas_package.raw_measurements_[0];
		py = meas_package.raw_measurements_[1];
		yaw = 0;
		phi =0;		
	}
	if (fabs(px) < 0.001)
	{
		px = 0.001;
	}
	if (fabs(py) < 0.001)
	{
		py = 0.001;
	} 
	x_ << px,py,yaw,0,0;

	// Calculating R
	 R_radar = MatrixXd(n_z_radar_,n_z_radar_);
	 R_radar << (std_radr_*std_radr_), 0, 0,
 		    0, (std_radphi_ * std_radphi_), 0,
 		    0, 0, (std_radrd_ * std_radrd_);

        // Calculating R
  	R_laser = MatrixXd(n_z_laser_ , n_z_laser_);
	R_laser << std_laspx_*std_laspx_,0,
       	     0,std_laspy_*std_laspy_; 		
	
	cout << "Initialized :" << endl;
	
	time_us_ = meas_package.timestamp_;
	
 }
 else
 {	     

  	g_count++;  	
	//Calculate time
  	double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;  //dt express in seconds
  	time_us_ = meas_package.timestamp_;
  		
	/*****************************************************************************
	   *  Prediction
	   ****************************************************************************/
	  /*if(g_count >= 100)
	  {
	  	return;
	  	//g_count++;
	  }*/
	  Prediction(dt);  
	  	  

	/*****************************************************************************
	   *  Update
	   ****************************************************************************/
	  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true)
	  {
	  	//cout << "Update RADAR: Debug:" << endl;
	  	
	  	if (meas_package.raw_measurements_[0] > 0.001)
	  	 	UpdateRadar(meas_package);
	  }
	  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true )
	  {
	  	double px = meas_package.raw_measurements_[0];
		double py = meas_package.raw_measurements_[1];
	  	//cout << "Update LASER: Debug:" << endl;
	  	if (((sqrt(px*px+py*py))) > 0.001)
	  	 	UpdateLidar(meas_package);
	  }
	  
  }
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  

  cout << "******************** Before prediction: *********************" << endl;
  cout << "x_ : " << endl;
  cout << x_ << endl;
  cout << "P_ : " << endl;
  cout << P_ << endl;

 
  // Create augumented mean vector
  VectorXd x_aug_ = VectorXd(n_aug_);
  
  // Create augmented covariance matrix
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
  
  // Create sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //cout << " Creating augumented mean state" << endl;  
  // Create augmented mean state
  x_aug_.head(n_x_) = x_;
  x_aug_(n_x_) = 0;
  x_aug_(n_x_ + 1) = 0;
       
  // Create augumeted covariance matirx
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_,n_x_) = P_;
  P_aug_(n_x_,n_x_) = std_a_ * std_a_;
  P_aug_(n_x_+1,n_x_+1) = std_yawdd_ * std_yawdd_ ;
    
  cout << "g_count :" <<g_count << endl;
  
  // Create square root of augmented covarance matrix
  MatrixXd A_aug_ = P_aug_.llt().matrixL();
  //cout<< "A_aug_: " << A_aug_ << endl;

  // Create augmented sigma points
  Xsig_aug_.col(0) = x_aug_;
  for (int i=0 ; i< n_aug_ ; i++)
  {
  	Xsig_aug_.col(i+1) = x_aug_ + sqrt(lambda_ + n_aug_) * A_aug_.col(i);
  	Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * A_aug_.col(i);
  }

  // Predict sigma points
  // Avoid division by zero
  // write predicted sigma points into the right column
  for (int i=0; i < sigma_points ; i++)
  {
  	double p_x_ = Xsig_aug_(0,i);
  	double p_y_ = Xsig_aug_(1,i);
  	double vk_ = Xsig_aug_(2,i);
  	double yaw_ = Xsig_aug_(3,i);
  	double yawd_ = Xsig_aug_(4,i);
  	double n_a_ = Xsig_aug_(5,i);
  	double n_yawdd_ = Xsig_aug_(6,i);
  	
  	double px_p, py_p=0;
  	
  	if (fabs(yawd_) > 0.001)
  	{
  		px_p = p_x_ + ((vk_ / yawd_) * (sin(yaw_ + (yawd_ * delta_t)) - sin(yaw_)));
  		py_p = p_y_ + ((vk_ / yawd_) * (cos(yaw_) - cos(yaw_ + (yawd_ * delta_t))));
  	}
  	else
  	{
  		px_p = p_x_ + vk_ * cos(yaw_) * delta_t;
  		py_p = p_y_ + vk_ * sin(yaw_) * delta_t;  
  	}

  	double v_p_ = vk_;
  	double yaw_p_ = yaw_ + yawd_ * delta_t;
  	double yawd_p_ = yawd_;
  	
  	// Add noise
  	px_p +=  0.5 * delta_t * delta_t * cos(yaw_) * n_a_;
  	py_p +=  0.5 * delta_t * delta_t * sin(yaw_) * n_a_;
  	v_p_ +=  delta_t * n_a_;
  	
  	yaw_p_ += 0.5 * delta_t * delta_t * n_yawdd_;
  	yawd_p_ += delta_t * n_yawdd_; 
  	
  	//Write the predicted sigma points into right column
  	Xsig_pred_(0,i) = px_p;
  	Xsig_pred_(1,i) = py_p;
  	Xsig_pred_(2,i) = v_p_;
  	Xsig_pred_(3,i) = yaw_p_;
  	Xsig_pred_(4,i) = yawd_p_;

  }

   
  //predict state mean
  x_.fill(0.0);
  for( int i=0; i < sigma_points; i++ )
  {
  	//Iterate over sigma points
  	x_ += weights_(i) * Xsig_pred_.col(i) ;   
  }
  
  //predict covariance matrix
  P_.fill(0.0);
  for( int i=0; i < sigma_points; i++ )
  {
  	//Iterate over sigma points
  	
  	VectorXd x_diff_ =  Xsig_pred_.col(i) - x_; 	
  	x_diff_(3)  = constrainAngle(x_diff_(3));		// Angel normalization
	P_ += weights_(i) * x_diff_ * x_diff_.transpose();
  	
  }
  cout << "******************* After prediction: ****************" << endl;
  cout << "x_ : " << endl;
  cout <<  x_ << endl;
  cout << "P_ : " << endl;
  cout <<  P_ << endl;
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
    
  //Create Matrix for sigma points to measurement space
  MatrixXd Zsig_l_ = MatrixXd(n_z_laser_, sigma_points);
  
  // Create vector for mean predicted measurement 
  VectorXd z_pred_l_ = VectorXd(n_z_laser_);
  
  // Create mean predicted covariance martix
  MatrixXd S = MatrixXd(n_z_laser_, n_z_laser_);
  Zsig_l_.fill(0.0);
  for (int i=0; i< sigma_points; i++)
  {
	Zsig_l_(0,i) = Xsig_pred_(0,i);
	Zsig_l_(1,i) = Xsig_pred_(1,i);
  } 
  
  z_pred_l_.fill(0.0);
  for (int i=0; i < sigma_points; i++)
  {
  	z_pred_l_ += weights_(i)*Zsig_l_.col(i);
  }
    
   // Calculating measurement covariance matrix S
   S.fill(0.0);
   for(int i =0; i< sigma_points; i++)
   {
   	VectorXd z_diff = Zsig_l_.col(i) - z_pred_l_;
   	S += weights_(i) * z_diff * z_diff.transpose();
   }
   S += R_laser;
   
   VectorXd z_laser_ = VectorXd(n_z_laser_);
   z_laser_ << meas_package.raw_measurements_;
   
   // Create Matrix for corrlation matrix
   MatrixXd Tc = MatrixXd(n_x_, n_z_laser_);
   Tc.fill(0.0);
      for (int i=0; i < sigma_points; i++)
   {
  	VectorXd z_diff = Zsig_l_.col(i) - z_pred_l_;	   	// Residual
   	VectorXd x_diff = Xsig_pred_.col(i) - x_;	   	// State difference
   	Tc += weights_(i)*x_diff*z_diff.transpose();
   }
   
   // Calculate Kalman gain
   MatrixXd k = MatrixXd(n_x_, n_z_laser_); 
   k = Tc * S.inverse();
  
   VectorXd z_diff_l_ = z_laser_ - z_pred_l_;  				   // Measurement residual
   
   // update mean state and covariance matrix
   x_ += k * z_diff_l_;
   P_ -= k * S * k.transpose();
   
   NIS_laser_ = z_diff_l_.transpose() * S.inverse() * z_diff_l_;          // Compute NIS
   
 cout << "******************Laser Measurement:*************" << endl;
 cout<< "x_: "<< x_ << endl;
 cout<< "P_: "<< P_ << endl;
 cout<< "NIS_laser_: "<< NIS_laser_ << endl;}

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
  
  
  // create matrix for sigma points to measurement space
  MatrixXd Zsig_r_ = MatrixXd(n_z_radar_, sigma_points);
  
 // mean predicted measurement
 VectorXd z_pred_r_ = VectorXd(n_z_radar_);
 
 // Mean covariance matrix S
 MatrixXd S  = MatrixXd(n_z_radar_,n_z_radar_);
 //cout << "Update :" <<endl;
 Zsig_r_.fill(0.0);
 
 for (int i=0; i < sigma_points ; i++)
 {
 	double px_ = Xsig_pred_(0,i);
 	double py_ = Xsig_pred_(1,i);
 	double v_ = Xsig_pred_(2,i);
 	double yaw_ = Xsig_pred_(3,i);

 	double v1 = cos(yaw_) * v_;
 	double v2 = sin(yaw_) * v_;
 	
 	double row_ = sqrt(px_ *px_ + py_ * py_);
 	double psi_ = atan2(py_,px_);
 	double yawd_= ((px_ * v1) + (py_ * v2 ))/row_;

 	if ( fabs(row_) > 0.001)
 	{
 		Zsig_r_(0,i) = row_;
 		Zsig_r_(1,i) = psi_;
 		Zsig_r_(2,i) = yawd_;
 		
 	}
 	else
 	{
 		Zsig_r_(0,i) = sqrt(0.001);
 		Zsig_r_(1,i) = atan2(0.001,0.001);
 		Zsig_r_(2,i) = yawd_;
 	} 	
 }
  
 // Calculate mean predicted measurement
 z_pred_r_.fill(0.0);
 for (int i=0; i < sigma_points; i++)
 {
 	z_pred_r_ += weights_(i)*Zsig_r_.col(i);
 }
 
 // Calculate meaurement covariance matrix S
 S.fill(0.0);
 for (int i=0 ; i < sigma_points; i++)
 {
 	VectorXd z_diff = Zsig_r_.col(i) - z_pred_r_;
	z_diff(1)  = constrainAngle(z_diff(1));			// Angle Normalization
 	S +=  weights_(i) * z_diff *z_diff.transpose();
 }
 
 S +=  R_radar;

 VectorXd z_radar = VectorXd(n_z_radar_);			// Create a vector for incoming radar measurement
 z_radar << meas_package.raw_measurements_;

 MatrixXd Tc = MatrixXd(n_x_,n_z_radar_);		// Create a matrix for correlation matrix
 Tc.fill(0.0);

 //Calculate cross correlataion matrix
 for (int i=0; i< sigma_points; i++)
 {
 	
 	VectorXd z_diff_ = Zsig_r_.col(i) - z_pred_r_;		//residual 
	z_diff_(1)  = constrainAngle(z_diff_(1));		// Angle normalization

 	VectorXd x_diff_ = Xsig_pred_.col(i) - x_;	 	//state difference
  	x_diff_(3)  = constrainAngle(x_diff_(3));		// Angle normalization

 	Tc += weights_(i) * x_diff_ * z_diff_.transpose();
 	
 }
 
 MatrixXd k = Tc * S.inverse();			// Calculate the kalman gain

 VectorXd z_diff_r = z_radar - z_pred_r_;		// Measurement residual
 z_diff_r(1)  = constrainAngle(z_diff_r(1));		// Angle Normalization
 
 
 // Update state mean and covariance matrix
 x_ +=  k * z_diff_r;
 P_ -=  k * S * k.transpose();
 
 // Compute NIS
 NIS_radar_ = z_diff_r.transpose() * S.inverse() * z_diff_r;

 cout << "******************Radar Measurement:*************" << endl;
 cout<< "x_: "<< x_ << endl;
 cout<< "P_: "<< P_ << endl;
 cout<< "NIS_radar_: "<< NIS_radar_ << endl;
  	   	 
}


