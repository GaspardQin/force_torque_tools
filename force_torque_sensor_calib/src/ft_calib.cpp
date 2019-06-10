/*
 *  ft_calib.cpp
 *
 *  Modified By Chuan QIN, July 2019
 * 	Logs: 
 * 		* Remove the require for IMU sensor, the TF of manipulator is used instead.
 *    * The base frame's inclination to ground is not required, only the local gravity acceleration is required (usually 9.81)
 *    * In order to accelerate the calibration SVD process, the force and torque are calibrated separately
 * 
 * 
 * 	Original file infomation:
 *  Created on: Sep 26, 2012
 *  Authors:   Francisco Viña
 *            fevb <at> kth.se
 */

/* Copyright (c) 2012, Francisco Viña, CVAP, KTH
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
      * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
      * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
      * Neither the name of KTH nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL KTH BE LIABLE FOR ANY
   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <ros/ros.h>
#include <force_torque_sensor_calib/ft_calib.h>
#include <Eigen/Dense>
#include "tf2_eigen/tf2_eigen.h"

namespace Calibration{

FTCalib::FTCalib(double m_local_gravitational_acceleration):m_local_gravitational_acceleration(m_local_gravitational_acceleration)
{
	m_num_meas = 0;
}

FTCalib::~FTCalib(){

}

void FTCalib::addMeasurement(const geometry_msgs::TransformStamped &tf_sensor_to_base,
		const geometry_msgs::WrenchStamped &ft_raw)
{
	if(tf_sensor_to_base.header.frame_id != ft_raw.header.frame_id)
	{
		ROS_ERROR("TF's source frame (%s) is not the same as the ft raw expressed frame (%s)!",
				tf_sensor_to_base.header.frame_id.c_str(), ft_raw.header.frame_id.c_str());
		return;
	}

	m_num_meas++;

	Eigen::Matrix<double, 3,6> A_force = getForceMeasurementMatrix(tf_sensor_to_base);
	Eigen::Matrix<double, 3,6> A_biased_torque = getBiasedTorqueMeasurementMatrix(ft_raw);
	
	Eigen::Vector3d measured_force, measured_torque;
	measured_force(0) = ft_raw.wrench.force.x;
	measured_force(1) = ft_raw.wrench.force.y;
	measured_force(2) = ft_raw.wrench.force.z;

	measured_torque(0) = ft_raw.wrench.torque.x;
	measured_torque(1) = ft_raw.wrench.torque.y;
	measured_torque(2) = ft_raw.wrench.torque.z;

	if(m_num_meas==1)
	{
		stacked_A_force = A_force;
		stacked_A_biased_torque = A_biased_torque;
		stacked_measured_force = measured_force;
		stacked_measured_torque = measured_torque;
	}

	else
	{
		Eigen::MatrixXd stacked_A_force_tmp = stacked_A_force;
		Eigen::MatrixXd stacked_A_biased_torque_tmp = stacked_A_biased_torque;
		
		Eigen::VectorXd stacked_measured_force_tmp = stacked_measured_force;
		Eigen::VectorXd stacked_measured_torque_tmp = stacked_measured_torque;

		stacked_A_force.resize(m_num_meas*3, 6);
		stacked_A_biased_torque.resize(m_num_meas*3, 6);

		stacked_measured_force.resize(m_num_meas*3);
		stacked_measured_torque.resize(m_num_meas*3);

		stacked_A_force.topRows((m_num_meas-1)*3) = stacked_A_force_tmp;
		stacked_A_biased_torque.topRows((m_num_meas-1)*3) = stacked_A_biased_torque_tmp;
		stacked_measured_force.topRows((m_num_meas-1)*3) = stacked_measured_force_tmp;
		stacked_measured_torque.topRows((m_num_meas-1)*3) = stacked_measured_torque_tmp;
		

		stacked_A_force.bottomRows(3) = A_force;
		stacked_A_biased_torque.bottomRows(3) = A_biased_torque;
		stacked_measured_force.bottomRows(3) = measured_force;
		stacked_measured_torque.bottomRows(3) = measured_torque;
	}


}
Eigen::Matrix<double, 3,6> FTCalib::getForceMeasurementMatrix(const geometry_msgs::TransformStamped &tf_sensor_to_base){
	Eigen::Matrix<double, 3, 6> A_force; // the force measurement matrix
	// [fx, fy, fz]^t = A_force * [g_x_in_base * m, g_y_in_base * m, g_z_in_base * m, f_bias_x, f_bias_y, f_bias_z]^t

	Eigen::Matrix3d R_sensor_to_base; // rotation from sensor frame to base frame

	// translate from geometry_msgs::TransformStamped to eigen matrix
	Eigen::Quaternion<double> q;
	Eigen::fromMsg(tf_sensor_to_base.transform.rotation, q);
	q.normalize();
	R_sensor_to_base = q.toRotationMatrix();

	A_force.block<3,3>(0,0) = R_sensor_to_base;
	A_force.block<3,3>(0,3).setIdentity();
	return A_force;
}

Eigen::Matrix<double, 3, 6> FTCalib::getBiasedTorqueMeasurementMatrix(const geometry_msgs::WrenchStamped &ft_raw){
	Eigen::Matrix<double, 3, 6> A_torque_biased; // the torque measurement matrix, which is biased because the force measurement is not unbiased
	// [tx, ty, tz]^t = A_torque * [x_center, y_center, z_center, t_bias_x, t_bias_y, t_bias_z]
	geometry_msgs::Vector3 force = ft_raw.wrench.force;
	A_torque_biased << 			0,	 force.z,	 force.y,	 1,	 0,	 0,
						 force.z, 					0,   force.x,	 0,	 1,	 0,
						 force.y, 	 force.x,					0,	 0,	 0,	 1;
	return A_torque_biased;
}

// Least squares to estimate the FT sensor parameters
Eigen::VectorXd FTCalib::getCalib()
{
	Eigen::VectorXd force_calib_params(6); //[g_x_in_base * mass, g_y_in_base * mass, g_z_in_base * mass, f_bias_x, f_bias_y, f_bias_z]

	force_calib_params = stacked_A_force.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(stacked_measured_force);

	Eigen::MatrixXd unbiased_stacked_A_torque = stacked_A_biased_torque;
	Eigen::Matrix<double,3,6> bias_matrix; 
	Eigen::Vector3d force_bias = force_calib_params.segment<3>(3);
	bias_matrix << 						0,	force_bias(2),	force_bias(1),	0,	0,	0,
								force_bias(2),							0,	force_bias(0),	0,	0,	0,
								force_bias(1), 	force_bias(0), 							0,	0,	0,	0;

	for(int i=0; i<m_num_meas; i++){
		unbiased_stacked_A_torque.block<3,6>(i*3,0) -= bias_matrix;
	}

	Eigen::VectorXd torque_calib_params(6); //[x_center_in_sensor_frame, y_center_in_sensor_frame, z_center_in_sensor_frame, t_bias_x, t_bias_y, t_bias_z] 
	torque_calib_params = unbiased_stacked_A_torque.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(stacked_measured_torque);

	// ft_calib_params: [mass, 
	//									 g_x_in_base, g_y_in_base, g_z_in_base,
	//									 x_center_in_sensor_frame, y_center_in_sensor_frame, z_center_in_sensor_frame, 
	//									 f_bias_x, f_bias_y, f_bias_z, 
	//									 t_bias_x, t_bias_y, t_bias_z]
	Eigen::VectorXd ft_calib_params;
	ft_calib_params.resize(13);
	ft_calib_params(0) = force_calib_params.head<3>().norm() / m_local_gravitational_acceleration; //mass
	ft_calib_params(1) = force_calib_params(0) / ft_calib_params(0); //g_x_in_base
	ft_calib_params(2) = force_calib_params(1) / ft_calib_params(0); //g_y_in_base
	ft_calib_params(3) = force_calib_params(2) / ft_calib_params(0); //g_z_in_base

	ft_calib_params(4) = torque_calib_params(0); //x_center_in_sensor_frame
	ft_calib_params(5) = torque_calib_params(1); //y_center_in_sensor_frame
	ft_calib_params(6) = torque_calib_params(2); //z_center_in_sensor_frame

	ft_calib_params(7) = force_calib_params(3); // f_bias_x
	ft_calib_params(8) = force_calib_params(4); // f_bias_y
	ft_calib_params(9) = force_calib_params(5); // f_bias_z

	ft_calib_params(10) = torque_calib_params(3); //t_bias_x
	ft_calib_params(11) = torque_calib_params(4); //t_bias_y
	ft_calib_params(12) = torque_calib_params(5); //t_bias_z
	 
	return ft_calib_params;
}



}
