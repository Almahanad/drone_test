=================================================================================================
File: ./all_cpp_files.cpp
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_gazebo/hector_gazebo_thermal_camera/src/gazebo_ros_thermal_depth_camera_plugin.cpp
//=================================================================================================
// Copyright (c) 2013, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include "gazebo_ros_thermal_camera.cpp"

// gazebo stuff
#include <gazebo/plugins/DepthCameraPlugin.hh>

namespace gazebo
{

template <>
void GazeboRosThermalCamera_<DepthCameraPlugin>::LoadImpl(sensors::SensorPtr _parent, sdf::ElementPtr _sdf)
{
  this->camera_ = this->DepthCameraPlugin::depthCamera;
}

template class GazeboRosThermalCamera_<DepthCameraPlugin>;
typedef GazeboRosThermalCamera_<DepthCameraPlugin> GazeboRosThermalCamera;

// Register this plugin with the simulator
GZ_REGISTER_SENSOR_PLUGIN(GazeboRosThermalCamera)

}
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_gazebo/hector_gazebo_thermal_camera/src/gazebo_ros_thermal_camera_plugin.cpp
//=================================================================================================
// Copyright (c) 2013, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include "gazebo_ros_thermal_camera.cpp"

// gazebo stuff
#include <gazebo/plugins/CameraPlugin.hh>

namespace gazebo
{

template <>
void GazeboRosThermalCamera_<CameraPlugin>::LoadImpl(sensors::SensorPtr _parent, sdf::ElementPtr _sdf)
{
  this->camera_ = this->CameraPlugin::camera;
}

template class GazeboRosThermalCamera_<CameraPlugin>;
typedef GazeboRosThermalCamera_<CameraPlugin> GazeboRosThermalCamera;

// Register this plugin with the simulator
GZ_REGISTER_SENSOR_PLUGIN(GazeboRosThermalCamera)

}
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_gazebo/hector_gazebo_thermal_camera/src/gazebo_ros_thermal_camera.cpp
//=================================================================================================
// Copyright (c) 2012, Stefan Kohlbrecher, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Simulation, Systems Optimization and Robotics
//       group, TU Darmstadt nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================
/**
 * Copy of the CameraSensor/DepthCameraSensor plugin with minor changes
 */

/*
 *  Gazebo - Outdoor Multi-Robot Simulator
 *  Copyright (C) 2003
 *     Nate Koenig & Andrew Howard
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include <hector_gazebo_thermal_camera/gazebo_ros_thermal_camera.h>

#include <gazebo/sensors/Sensor.hh>
#include <gazebo/sensors/SensorTypes.hh>

#include <sensor_msgs/image_encodings.h>

#include <gazebo/gazebo_config.h>

namespace gazebo
{

////////////////////////////////////////////////////////////////////////////////
// Constructor
template <class Base>
GazeboRosThermalCamera_<Base>::GazeboRosThermalCamera_()
{
}

////////////////////////////////////////////////////////////////////////////////
// Destructor
template <class Base>
GazeboRosThermalCamera_<Base>::~GazeboRosThermalCamera_()
{
}

template <class Base>
void GazeboRosThermalCamera_<Base>::Load(sensors::SensorPtr _parent, sdf::ElementPtr _sdf)
{
  Base::Load(_parent, _sdf);
  // copying from CameraPlugin into GazeboRosCameraUtils
  this->parentSensor_ = this->parentSensor;
  this->width_ = this->width;
  this->height_ = this->height;
  this->depth_ = this->depth;
  this->format_ = this->format;

  this->image_connect_count_ = boost::shared_ptr<int>(new int);
  *this->image_connect_count_ = 0;
  this->image_connect_count_lock_ = boost::shared_ptr<boost::mutex>(new boost::mutex);
  this->was_active_ = boost::shared_ptr<bool>(new bool);
  *this->was_active_ = false;

  LoadImpl(_parent, _sdf);
  GazeboRosCameraUtils::Load(_parent, _sdf);
}

////////////////////////////////////////////////////////////////////////////////
// Update the controller
template <class Base>
void GazeboRosThermalCamera_<Base>::OnNewFrame(const unsigned char *_image,
    unsigned int _width, unsigned int _height, unsigned int _depth,
    const std::string &_format)
{
  if (!this->initialized_ || this->height_ <=0 || this->width_ <=0)
    return;

#if (GAZEBO_MAJOR_VERSION > 6)
  this->sensor_update_time_ = this->parentSensor_->LastUpdateTime();
#else
  this->sensor_update_time_ = this->parentSensor_->GetLastUpdateTime();
#endif

  if (!this->parentSensor->IsActive())
  {
    if ((*this->image_connect_count_) > 0)
      // do this first so there's chance for sensor to run 1 frame after activate
      this->parentSensor->SetActive(true);
  }
  else
  {
    if ((*this->image_connect_count_) > 0)
    {
#if (GAZEBO_MAJOR_VERSION >= 8)
      common::Time cur_time = this->world_->SimTime();
#else
      common::Time cur_time = this->world_->GetSimTime();
#endif
      if (cur_time - this->last_update_time_ >= this->update_period_)
      {
        this->PutCameraData(_image);
        this->PublishCameraInfo();
        this->last_update_time_ = cur_time;
      }
    }
  }
}

template <class Base>
void GazeboRosThermalCamera_<Base>::OnNewImageFrame(const unsigned char *_image,
    unsigned int _width, unsigned int _height, unsigned int _depth,
    const std::string &_format)
{
  OnNewFrame(_image, _width, _height, _depth, _format);
}

////////////////////////////////////////////////////////////////////////////////
// Put camera_ data to the interface
template <class Base>
void GazeboRosThermalCamera_<Base>::PutCameraData(const unsigned char *_src, common::Time &last_update_time)
{
  this->sensor_update_time_ = last_update_time;
  this->PutCameraData(_src);
}

template <class Base>
void GazeboRosThermalCamera_<Base>::PutCameraData(const unsigned char *_src)
{
  if (!this->initialized_ || this->height_ <=0 || this->width_ <=0)
    return;

  this->lock_.lock();

  // copy data into image
  this->image_msg_.header.frame_id = this->frame_name_;
  this->image_msg_.header.stamp.sec = this->sensor_update_time_.sec;
  this->image_msg_.header.stamp.nsec = this->sensor_update_time_.nsec;

  /// don't bother if there are no subscribers
  if ((*this->image_connect_count_) > 0)
  {
    this->image_msg_.width = this->width_;
    this->image_msg_.height = this->height_;
    this->image_msg_.encoding = sensor_msgs::image_encodings::MONO8;
    this->image_msg_.step = this->image_msg_.width;

    size_t size = this->width_ * this->height_;

    std::vector<uint8_t>& data (this->image_msg_.data);
    data.resize(size);

    size_t img_index = 0;

    for (size_t i = 0; i < size; ++i){
      if ((_src[img_index] >254) && (_src[img_index+1] < 1) && (_src[img_index+2] < 1)){
        //RGB [255,0,0] translates to white (white hot)
        data[i]= 255;
      }else{
        //Everything else is written to the MONO8 output image much darker
        data[i]= (_src[img_index] + _src[img_index+1] + _src[img_index+2]) /8 ;
      }
      img_index += 3;
    }

    // publish to ros
    this->image_pub_.publish(this->image_msg_);
  }

  this->lock_.unlock();
}

}
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_gazebo/hector_gazebo_plugins/src/gazebo_ros_gps.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_gazebo_plugins/gazebo_ros_gps.h>
#include <gazebo/physics/physics.hh>

// WGS84 constants
static const double equatorial_radius = 6378137.0;
static const double flattening = 1.0/298.257223563;
static const double excentrity2 = 2*flattening - flattening*flattening;

// default reference position
static const double DEFAULT_REFERENCE_LATITUDE  = 49.9;
static const double DEFAULT_REFERENCE_LONGITUDE = 8.9;
static const double DEFAULT_REFERENCE_HEADING   = 0.0;
static const double DEFAULT_REFERENCE_ALTITUDE  = 0.0;

namespace gazebo {

GazeboRosGps::GazeboRosGps()
{
}

////////////////////////////////////////////////////////////////////////////////
// Destructor
GazeboRosGps::~GazeboRosGps()
{
  updateTimer.Disconnect(updateConnection);

  dynamic_reconfigure_server_position_.reset();
  dynamic_reconfigure_server_velocity_.reset();
  dynamic_reconfigure_server_status_.reset();

  node_handle_->shutdown();
  delete node_handle_;
}

////////////////////////////////////////////////////////////////////////////////
// Load the controller
void GazeboRosGps::Load(physics::ModelPtr _model, sdf::ElementPtr _sdf)
{
  world = _model->GetWorld();

  // load parameters
  if (!_sdf->HasElement("robotNamespace"))
    namespace_.clear();
  else
    namespace_ = _sdf->GetElement("robotNamespace")->GetValue()->GetAsString();

  if (!_sdf->HasElement("bodyName"))
  {
    link = _model->GetLink();
    link_name_ = link->GetName();
  }
  else {
    link_name_ = _sdf->GetElement("bodyName")->GetValue()->GetAsString();
    link = _model->GetLink(link_name_);
  }

  if (!link)
  {
    ROS_FATAL("GazeboRosGps plugin error: bodyName: %s does not exist\n", link_name_.c_str());
    return;
  }

  // default parameters
  frame_id_ = "/world";
  fix_topic_ = "fix";
  velocity_topic_ = "fix_velocity";

  reference_latitude_  = DEFAULT_REFERENCE_LATITUDE;
  reference_longitude_ = DEFAULT_REFERENCE_LONGITUDE;
  reference_heading_   = DEFAULT_REFERENCE_HEADING * M_PI/180.0;
  reference_altitude_  = DEFAULT_REFERENCE_ALTITUDE;

  fix_.status.status  = sensor_msgs::NavSatStatus::STATUS_FIX;
  fix_.status.service = 0;

  if (_sdf->HasElement("frameId"))
    frame_id_ = _sdf->GetElement("frameId")->GetValue()->GetAsString();

  if (_sdf->HasElement("topicName"))
    fix_topic_ = _sdf->GetElement("topicName")->GetValue()->GetAsString();

  if (_sdf->HasElement("velocityTopicName"))
    velocity_topic_ = _sdf->GetElement("velocityTopicName")->GetValue()->GetAsString();

  if (_sdf->HasElement("referenceLatitude"))
    _sdf->GetElement("referenceLatitude")->GetValue()->Get(reference_latitude_);

  if (_sdf->HasElement("referenceLongitude"))
    _sdf->GetElement("referenceLongitude")->GetValue()->Get(reference_longitude_);

  if (_sdf->HasElement("referenceHeading"))
    if (_sdf->GetElement("referenceHeading")->GetValue()->Get(reference_heading_))
      reference_heading_ *= M_PI/180.0;

  if (_sdf->HasElement("referenceAltitude"))
    _sdf->GetElement("referenceAltitude")->GetValue()->Get(reference_altitude_);

  if (_sdf->HasElement("status")) {
    int status = fix_.status.status;
    if (_sdf->GetElement("status")->GetValue()->Get(status))
      fix_.status.status = static_cast<sensor_msgs::NavSatStatus::_status_type>(status);
  }

  if (_sdf->HasElement("service")) {
    unsigned int service = fix_.status.service;
    if (_sdf->GetElement("service")->GetValue()->Get(service))
      fix_.status.service = static_cast<sensor_msgs::NavSatStatus::_service_type>(service);
  }

  fix_.header.frame_id = frame_id_;
  velocity_.header.frame_id = frame_id_;

  position_error_model_.Load(_sdf);
  velocity_error_model_.Load(_sdf, "velocity");

  // calculate earth radii
  double temp = 1.0 / (1.0 - excentrity2 * sin(reference_latitude_ * M_PI/180.0) * sin(reference_latitude_ * M_PI/180.0));
  double prime_vertical_radius = equatorial_radius * sqrt(temp);
  radius_north_ = prime_vertical_radius * (1 - excentrity2) * temp;
  radius_east_  = prime_vertical_radius * cos(reference_latitude_ * M_PI/180.0);

  // Make sure the ROS node for Gazebo has already been initialized
  if (!ros::isInitialized())
  {
    ROS_FATAL_STREAM("A ROS node for Gazebo has not been initialized, unable to load plugin. "
      << "Load the Gazebo system plugin 'libgazebo_ros_api_plugin.so' in the gazebo_ros package)");
    return;
  }

  node_handle_ = new ros::NodeHandle(namespace_);
  fix_publisher_ = node_handle_->advertise<sensor_msgs::NavSatFix>(fix_topic_, 10);
  velocity_publisher_ = node_handle_->advertise<geometry_msgs::Vector3Stamped>(velocity_topic_, 10);

  // setup dynamic_reconfigure servers
  dynamic_reconfigure_server_position_.reset(new dynamic_reconfigure::Server<SensorModelConfig>(ros::NodeHandle(*node_handle_, fix_topic_ + "/position")));
  dynamic_reconfigure_server_velocity_.reset(new dynamic_reconfigure::Server<SensorModelConfig>(ros::NodeHandle(*node_handle_, fix_topic_ + "/velocity")));
  dynamic_reconfigure_server_status_.reset(new dynamic_reconfigure::Server<GNSSConfig>(ros::NodeHandle(*node_handle_, fix_topic_ + "/status")));
  dynamic_reconfigure_server_position_->setCallback(boost::bind(&SensorModel3::dynamicReconfigureCallback, &position_error_model_, _1, _2));
  dynamic_reconfigure_server_velocity_->setCallback(boost::bind(&SensorModel3::dynamicReconfigureCallback, &velocity_error_model_, _1, _2));
  dynamic_reconfigure_server_status_->setCallback(boost::bind(&GazeboRosGps::dynamicReconfigureCallback, this, _1, _2));

  Reset();

  // connect Update function
  updateTimer.setUpdateRate(4.0);
  updateTimer.Load(world, _sdf);
  updateConnection = updateTimer.Connect(boost::bind(&GazeboRosGps::Update, this));
}

void GazeboRosGps::Reset()
{
  updateTimer.Reset();
  position_error_model_.reset();
  velocity_error_model_.reset();
}

void GazeboRosGps::dynamicReconfigureCallback(GazeboRosGps::GNSSConfig &config, uint32_t level)
{
  using sensor_msgs::NavSatStatus;
  if (level == 1) {
    if (!config.STATUS_FIX) {
      fix_.status.status = NavSatStatus::STATUS_NO_FIX;
    } else {
      fix_.status.status = (config.STATUS_SBAS_FIX ? NavSatStatus::STATUS_SBAS_FIX : 0) |
                           (config.STATUS_GBAS_FIX ? NavSatStatus::STATUS_GBAS_FIX : 0);
    }
    fix_.status.service = (config.SERVICE_GPS     ? NavSatStatus::SERVICE_GPS : 0) |
                          (config.SERVICE_GLONASS ? NavSatStatus::SERVICE_GLONASS : 0) |
                          (config.SERVICE_COMPASS ? NavSatStatus::SERVICE_COMPASS : 0) |
                          (config.SERVICE_GALILEO ? NavSatStatus::SERVICE_GALILEO : 0);
  } else {
    config.STATUS_FIX      = (fix_.status.status != NavSatStatus::STATUS_NO_FIX);
    config.STATUS_SBAS_FIX = (fix_.status.status & NavSatStatus::STATUS_SBAS_FIX);
    config.STATUS_GBAS_FIX = (fix_.status.status & NavSatStatus::STATUS_GBAS_FIX);
    config.SERVICE_GPS     = (fix_.status.service & NavSatStatus::SERVICE_GPS);
    config.SERVICE_GLONASS = (fix_.status.service & NavSatStatus::SERVICE_GLONASS);
    config.SERVICE_COMPASS = (fix_.status.service & NavSatStatus::SERVICE_COMPASS);
    config.SERVICE_GALILEO = (fix_.status.service & NavSatStatus::SERVICE_GALILEO);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Update the controller
void GazeboRosGps::Update()
{
#if (GAZEBO_MAJOR_VERSION >= 8)
  common::Time sim_time = world->SimTime();
  double dt = updateTimer.getTimeSinceLastUpdate().Double();

  ignition::math::Pose3d pose = link->WorldPose();

  ignition::math::Vector3d velocity = velocity_error_model_(link->WorldLinearVel(), dt);
  ignition::math::Vector3d position = position_error_model_(pose.Pos(), dt);
#else
  common::Time sim_time = world->GetSimTime();
  double dt = updateTimer.getTimeSinceLastUpdate().Double();

  math::Pose pose = link->GetWorldPose();

  gazebo::math::Vector3 velocity = velocity_error_model_(link->GetWorldLinearVel(), dt);
  gazebo::math::Vector3 position = position_error_model_(pose.pos, dt);
#endif

  // An offset error in the velocity is integrated into the position error for the next timestep.
  // Note: Usually GNSS receivers have almost no drift in the velocity signal.
  position_error_model_.setCurrentDrift(position_error_model_.getCurrentDrift() + dt * velocity_error_model_.getCurrentDrift());

  fix_.header.stamp = ros::Time(sim_time.sec, sim_time.nsec);
  velocity_.header.stamp = fix_.header.stamp;

#if (GAZEBO_MAJOR_VERSION >= 8)
  fix_.latitude  = reference_latitude_  + ( cos(reference_heading_) * position.X() + sin(reference_heading_) * position.Y()) / radius_north_ * 180.0/M_PI;
  fix_.longitude = reference_longitude_ - (-sin(reference_heading_) * position.X() + cos(reference_heading_) * position.Y()) / radius_east_  * 180.0/M_PI;
  fix_.altitude  = reference_altitude_  + position.Z();
  velocity_.vector.x =  cos(reference_heading_) * velocity.X() + sin(reference_heading_) * velocity.Y();
  velocity_.vector.y = -sin(reference_heading_) * velocity.X() + cos(reference_heading_) * velocity.Y();
  velocity_.vector.z = velocity.Z();

  fix_.position_covariance_type = sensor_msgs::NavSatFix::COVARIANCE_TYPE_DIAGONAL_KNOWN;
  fix_.position_covariance[0] = position_error_model_.drift.X()*position_error_model_.drift.X() + position_error_model_.gaussian_noise.X()*position_error_model_.gaussian_noise.X();
  fix_.position_covariance[4] = position_error_model_.drift.Y()*position_error_model_.drift.Y() + position_error_model_.gaussian_noise.Y()*position_error_model_.gaussian_noise.Y();
  fix_.position_covariance[8] = position_error_model_.drift.Z()*position_error_model_.drift.Z() + position_error_model_.gaussian_noise.Z()*position_error_model_.gaussian_noise.Z();
#else
  fix_.latitude  = reference_latitude_  + ( cos(reference_heading_) * position.x + sin(reference_heading_) * position.y) / radius_north_ * 180.0/M_PI;
  fix_.longitude = reference_longitude_ - (-sin(reference_heading_) * position.x + cos(reference_heading_) * position.y) / radius_east_  * 180.0/M_PI;
  fix_.altitude  = reference_altitude_  + position.z;
  velocity_.vector.x =  cos(reference_heading_) * velocity.x + sin(reference_heading_) * velocity.y;
  velocity_.vector.y = -sin(reference_heading_) * velocity.x + cos(reference_heading_) * velocity.y;
  velocity_.vector.z = velocity.z;

  fix_.position_covariance_type = sensor_msgs::NavSatFix::COVARIANCE_TYPE_DIAGONAL_KNOWN;
  fix_.position_covariance[0] = position_error_model_.drift.x*position_error_model_.drift.x + position_error_model_.gaussian_noise.x*position_error_model_.gaussian_noise.x;
  fix_.position_covariance[4] = position_error_model_.drift.y*position_error_model_.drift.y + position_error_model_.gaussian_noise.y*position_error_model_.gaussian_noise.y;
  fix_.position_covariance[8] = position_error_model_.drift.z*position_error_model_.drift.z + position_error_model_.gaussian_noise.z*position_error_model_.gaussian_noise.z;
#endif

  fix_publisher_.publish(fix_);
  velocity_publisher_.publish(velocity_);
}

// Register this plugin with the simulator
GZ_REGISTER_MODEL_PLUGIN(GazeboRosGps)

} // namespace gazebo
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_gazebo/hector_gazebo_plugins/src/gazebo_ros_imu.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_gazebo_plugins/gazebo_ros_imu.h>
#include <gazebo/common/Events.hh>
#include <gazebo/physics/physics.hh>

#include <ros/console.h>

namespace gazebo
{

// #define DEBUG_OUTPUT
#ifdef DEBUG_OUTPUT
  #include <geometry_msgs/PoseStamped.h>
  static ros::Publisher debugPublisher;
#endif // DEBUG_OUTPUT

////////////////////////////////////////////////////////////////////////////////
// Constructor
GazeboRosIMU::GazeboRosIMU()
{
}

////////////////////////////////////////////////////////////////////////////////
// Destructor
GazeboRosIMU::~GazeboRosIMU()
{
  updateTimer.Disconnect(updateConnection);

  dynamic_reconfigure_server_accel_.reset();
  dynamic_reconfigure_server_rate_.reset();
  dynamic_reconfigure_server_yaw_.reset();

  node_handle_->shutdown();
#ifdef USE_CBQ
  callback_queue_thread_.join();
#endif
  delete node_handle_;
}

////////////////////////////////////////////////////////////////////////////////
// Load the controller
void GazeboRosIMU::Load(physics::ModelPtr _model, sdf::ElementPtr _sdf)
{
  // Get the world name.
  world = _model->GetWorld();

  // load parameters
  if (_sdf->HasElement("robotNamespace"))
    namespace_ = _sdf->GetElement("robotNamespace")->GetValue()->GetAsString();
  else
    namespace_.clear();

  if (_sdf->HasElement("bodyName"))
  {
    link_name_ = _sdf->GetElement("bodyName")->GetValue()->GetAsString();
    link = _model->GetLink(link_name_);
  }
  else
  {
    link = _model->GetLink();
    link_name_ = link->GetName();
  }

  // assert that the body by link_name_ exists
  if (!link)
  {
    ROS_FATAL("GazeboRosIMU plugin error: bodyName: %s does not exist\n", link_name_.c_str());
    return;
  }

  // default parameters
  frame_id_ = link_name_;
  topic_ = "imu";

  if (_sdf->HasElement("frameId"))
    frame_id_ = _sdf->GetElement("frameId")->GetValue()->GetAsString();

  if (_sdf->HasElement("topicName"))
    topic_ = _sdf->GetElement("topicName")->GetValue()->GetAsString();

  if (_sdf->HasElement("biasTopicName"))
    bias_topic_ = _sdf->GetElement("biasTopicName")->GetValue()->GetAsString();
  else
    bias_topic_ = (!topic_.empty() ? topic_ + "/bias" : "");

  if (_sdf->HasElement("serviceName"))
    serviceName = _sdf->GetElement("serviceName")->GetValue()->GetAsString();
  else
    serviceName = topic_ + "/calibrate";

  accelModel.Load(_sdf, "accel");
  rateModel.Load(_sdf, "rate");
  yawModel.Load(_sdf, "yaw");

  // also use old configuration variables from gazebo_ros_imu
  if (_sdf->HasElement("gaussianNoise")) {
    double gaussianNoise;
    if (_sdf->GetElement("gaussianNoise")->GetValue()->Get(gaussianNoise) && gaussianNoise != 0.0) {
      accelModel.gaussian_noise = gaussianNoise;
      rateModel.gaussian_noise  = gaussianNoise;
    }
  }

  if (_sdf->HasElement("xyzOffset")) {
#if (GAZEBO_MAJOR_VERSION >= 8)
    this->offset_.Pos() = _sdf->Get<ignition::math::Vector3d>("xyzOffset");
#else
    this->offset_.pos = _sdf->Get<math::Vector3>("xyzOffset");
#endif
  } else {
    ROS_INFO("imu plugin missing <xyzOffset>, defaults to 0s");
#if (GAZEBO_MAJOR_VERSION >= 8)
    this->offset_.Pos() = ignition::math::Vector3d(0, 0, 0);
#else
    this->offset_.pos = math::Vector3(0, 0, 0);
#endif
  }

  if (_sdf->HasElement("rpyOffset")) {
#if (GAZEBO_MAJOR_VERSION >= 8)
    this->offset_.Rot() = _sdf->Get<ignition::math::Quaterniond>("rpyOffset");
#else
    this->offset_.rot = _sdf->Get<math::Vector3>("rpyOffset");
#endif
  } else {
    ROS_INFO("imu plugin missing <rpyOffset>, defaults to 0s");
#if (GAZEBO_MAJOR_VERSION >= 8)
    this->offset_.Rot() = ignition::math::Quaterniond(0, 0, 0);
#else
    this->offset_.rot = math::Vector3(0, 0, 0);
#endif
  }

  // fill in constant covariance matrix
#if (GAZEBO_MAJOR_VERSION >= 8)
  imuMsg.angular_velocity_covariance[0] = rateModel.gaussian_noise.X()*rateModel.gaussian_noise.X();
  imuMsg.angular_velocity_covariance[4] = rateModel.gaussian_noise.Y()*rateModel.gaussian_noise.Y();
  imuMsg.angular_velocity_covariance[8] = rateModel.gaussian_noise.Z()*rateModel.gaussian_noise.Z();
  imuMsg.linear_acceleration_covariance[0] = accelModel.gaussian_noise.X()*accelModel.gaussian_noise.X();
  imuMsg.linear_acceleration_covariance[4] = accelModel.gaussian_noise.Y()*accelModel.gaussian_noise.Y();
  imuMsg.linear_acceleration_covariance[8] = accelModel.gaussian_noise.Z()*accelModel.gaussian_noise.Z();
#else
  imuMsg.angular_velocity_covariance[0] = rateModel.gaussian_noise.x*rateModel.gaussian_noise.x;
  imuMsg.angular_velocity_covariance[4] = rateModel.gaussian_noise.y*rateModel.gaussian_noise.y;
  imuMsg.angular_velocity_covariance[8] = rateModel.gaussian_noise.z*rateModel.gaussian_noise.z;
  imuMsg.linear_acceleration_covariance[0] = accelModel.gaussian_noise.x*accelModel.gaussian_noise.x;
  imuMsg.linear_acceleration_covariance[4] = accelModel.gaussian_noise.y*accelModel.gaussian_noise.y;
  imuMsg.linear_acceleration_covariance[8] = accelModel.gaussian_noise.z*accelModel.gaussian_noise.z;
#endif

  // Make sure the ROS node for Gazebo has already been initialized
  if (!ros::isInitialized())
  {
    ROS_FATAL_STREAM("A ROS node for Gazebo has not been initialized, unable to load plugin. "
      << "Load the Gazebo system plugin 'libgazebo_ros_api_plugin.so' in the gazebo_ros package.");
    return;
  }

  node_handle_ = new ros::NodeHandle(namespace_);

  // if topic name specified as empty, do not publish (then what is this plugin good for?)
  if (!topic_.empty())
    pub_ = node_handle_->advertise<sensor_msgs::Imu>(topic_, 10);
  if (!bias_topic_.empty())
    bias_pub_ = node_handle_->advertise<sensor_msgs::Imu>(bias_topic_, 10);

#ifdef DEBUG_OUTPUT
  debugPublisher = rosnode_->advertise<geometry_msgs::PoseStamped>(topic_ + "/pose", 10);
#endif // DEBUG_OUTPUT

  // advertise services for calibration and bias setting
  if (!serviceName.empty())
    srv_ = node_handle_->advertiseService(serviceName, &GazeboRosIMU::ServiceCallback, this);

  accelBiasService = node_handle_->advertiseService(topic_ + "/set_accel_bias", &GazeboRosIMU::SetAccelBiasCallback, this);
  rateBiasService  = node_handle_->advertiseService(topic_ + "/set_rate_bias", &GazeboRosIMU::SetRateBiasCallback, this);

  // setup dynamic_reconfigure servers
  if (!topic_.empty()) {
    dynamic_reconfigure_server_accel_.reset(new dynamic_reconfigure::Server<SensorModelConfig>(ros::NodeHandle(*node_handle_, topic_ + "/accel")));
    dynamic_reconfigure_server_rate_.reset(new dynamic_reconfigure::Server<SensorModelConfig>(ros::NodeHandle(*node_handle_, topic_ + "/rate")));
    dynamic_reconfigure_server_yaw_.reset(new dynamic_reconfigure::Server<SensorModelConfig>(ros::NodeHandle(*node_handle_, topic_ + "/yaw")));
    dynamic_reconfigure_server_accel_->setCallback(boost::bind(&SensorModel3::dynamicReconfigureCallback, &accelModel, _1, _2));
    dynamic_reconfigure_server_rate_->setCallback(boost::bind(&SensorModel3::dynamicReconfigureCallback, &rateModel, _1, _2));
    dynamic_reconfigure_server_yaw_->setCallback(boost::bind(&SensorModel::dynamicReconfigureCallback, &yawModel, _1, _2));
  }

#ifdef USE_CBQ
  // start custom queue for imu
  callback_queue_thread_ = boost::thread( boost::bind( &GazeboRosIMU::CallbackQueueThread,this ) );
#endif

  Reset();

  // connect Update function
  updateTimer.Load(world, _sdf);
  updateConnection = updateTimer.Connect(boost::bind(&GazeboRosIMU::Update, this));
}

void GazeboRosIMU::Reset()
{
  updateTimer.Reset();

#if (GAZEBO_MAJOR_VERSION >= 8)
  orientation = ignition::math::Quaterniond();
#else
  orientation = math::Quaternion();
#endif
  velocity = 0.0;
  accel = 0.0;

  accelModel.reset();
  rateModel.reset();
  yawModel.reset();
}

////////////////////////////////////////////////////////////////////////////////
// returns true always, imu is always calibrated in sim
bool GazeboRosIMU::ServiceCallback(std_srvs::Empty::Request &req,
                                        std_srvs::Empty::Response &res)
{
  boost::mutex::scoped_lock scoped_lock(lock);
#if (GAZEBO_MAJOR_VERSION >= 8)
  rateModel.reset(ignition::math::Vector3d(0.0, 0.0, 0.0));
#else
  rateModel.reset(math::Vector3(0.0, 0.0, 0.0));
#endif
  return true;
}

bool GazeboRosIMU::SetAccelBiasCallback(hector_gazebo_plugins::SetBias::Request &req, hector_gazebo_plugins::SetBias::Response &res)
{
  boost::mutex::scoped_lock scoped_lock(lock);
#if (GAZEBO_MAJOR_VERSION >= 8)
  accelModel.reset(ignition::math::Vector3d(req.bias.x, req.bias.y, req.bias.z));
#else
  accelModel.reset(math::Vector3(req.bias.x, req.bias.y, req.bias.z));
#endif
  return true;
}

bool GazeboRosIMU::SetRateBiasCallback(hector_gazebo_plugins::SetBias::Request &req, hector_gazebo_plugins::SetBias::Response &res)
{
  boost::mutex::scoped_lock scoped_lock(lock);
#if (GAZEBO_MAJOR_VERSION >= 8)
  rateModel.reset(ignition::math::Vector3d(req.bias.x, req.bias.y, req.bias.z));
#else
  rateModel.reset(math::Vector3(req.bias.x, req.bias.y, req.bias.z));
#endif
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// Update the controller
void GazeboRosIMU::Update()
{
  // Get Time Difference dt
#if (GAZEBO_MAJOR_VERSION >= 8)
  common::Time cur_time = world->SimTime();
#else
  common::Time cur_time = world->GetSimTime();
#endif
  double dt = updateTimer.getTimeSinceLastUpdate().Double();
  boost::mutex::scoped_lock scoped_lock(lock);

  // Get Pose/Orientation
#if (GAZEBO_MAJOR_VERSION >= 8)
  ignition::math::Pose3d pose = link->WorldPose();
  // ignition::math::Vector3d pos = pose.pos + this->offset_.pos;
  ignition::math::Quaterniond rot = this->offset_.Rot() * pose.Rot();
#else
  math::Pose pose = link->GetWorldPose();
  // math::Vector3 pos = pose.pos + this->offset_.pos;
  math::Quaternion rot = this->offset_.rot * pose.rot;
#endif
  rot.Normalize();

  // get Gravity
#if (GAZEBO_MAJOR_VERSION >= 8)
  gravity = world->Gravity();
  double gravity_length = gravity.Length();
  ROS_DEBUG_NAMED("gazebo_ros_imu", "gravity_world = [%g %g %g]", gravity.X(), gravity.Y(), gravity.Z());
#else
  gravity = world->GetPhysicsEngine()->GetGravity();
  double gravity_length = gravity.GetLength();
  ROS_DEBUG_NAMED("gazebo_ros_imu", "gravity_world = [%g %g %g]", gravity.x, gravity.y, gravity.z);
#endif

  // get Acceleration and Angular Rates
  // the result of GetRelativeLinearAccel() seems to be unreliable (sum of forces added during the current simulation step)?
  //accel = myBody->GetRelativeLinearAccel(); // get acceleration in body frame
#if (GAZEBO_MAJOR_VERSION >= 8)
  ignition::math::Vector3d temp = link->WorldLinearVel(); // get velocity in world frame
#else
  math::Vector3 temp = link->GetWorldLinearVel(); // get velocity in world frame
#endif
  if (dt > 0.0) accel = rot.RotateVectorReverse((temp - velocity) / dt - gravity);
  velocity = temp;

  // calculate angular velocity from delta quaternion
  // note: link->GetRelativeAngularVel() sometimes return nan?
  // rate  = link->GetRelativeAngularVel(); // get angular rate in body frame
#if (GAZEBO_MAJOR_VERSION >= 8)
  ignition::math::Quaterniond delta = this->orientation.Inverse() * rot;
#else
  math::Quaternion delta = this->orientation.GetInverse() * rot;
#endif
  this->orientation = rot;
  if (dt > 0.0) {
#if (GAZEBO_MAJOR_VERSION >= 8)
    rate = this->offset_.Rot().Inverse()
           * (2.0 * acos(std::max(std::min(delta.W(), 1.0), -1.0)) * ignition::math::Vector3d(delta.X(), delta.Y(), delta.Z()).Normalize() / dt);
#else
    rate = this->offset_.rot.GetInverse()
           * (2.0 * acos(std::max(std::min(delta.w, 1.0), -1.0)) * math::Vector3(delta.x, delta.y, delta.z).Normalize() / dt);
#endif
  }

  // update sensor models
  accel = accelModel(accel, dt);
  rate  = rateModel(rate, dt);
  yawModel.update(dt);
#if (GAZEBO_MAJOR_VERSION >= 8)
  ROS_DEBUG_NAMED("gazebo_ros_imu", "Current bias errors: accel = [%g %g %g], rate = [%g %g %g], yaw = %g",
                 accelModel.getCurrentBias().X(), accelModel.getCurrentBias().Y(), accelModel.getCurrentBias().Z(),
                 rateModel.getCurrentBias().X(), rateModel.getCurrentBias().Y(), rateModel.getCurrentBias().Z(),
                 yawModel.getCurrentBias());
  ROS_DEBUG_NAMED("gazebo_ros_imu", "Scale errors: accel = [%g %g %g], rate = [%g %g %g], yaw = %g",
                 accelModel.getScaleError().X(), accelModel.getScaleError().Y(), accelModel.getScaleError().Z(),
                 rateModel.getScaleError().X(), rateModel.getScaleError().Y(), rateModel.getScaleError().Z(),
                 yawModel.getScaleError());
#else
  ROS_DEBUG_NAMED("gazebo_ros_imu", "Current bias errors: accel = [%g %g %g], rate = [%g %g %g], yaw = %g",
                 accelModel.getCurrentBias().x, accelModel.getCurrentBias().y, accelModel.getCurrentBias().z,
                 rateModel.getCurrentBias().x, rateModel.getCurrentBias().y, rateModel.getCurrentBias().z,
                 yawModel.getCurrentBias());
  ROS_DEBUG_NAMED("gazebo_ros_imu", "Scale errors: accel = [%g %g %g], rate = [%g %g %g], yaw = %g",
                 accelModel.getScaleError().x, accelModel.getScaleError().y, accelModel.getScaleError().z,
                 rateModel.getScaleError().x, rateModel.getScaleError().y, rateModel.getScaleError().z,
                 yawModel.getScaleError());
#endif

  // apply accelerometer and yaw drift error to orientation (pseudo AHRS)
#if (GAZEBO_MAJOR_VERSION >= 8)
  ignition::math::Vector3d accelDrift = pose.Rot().RotateVector(accelModel.getCurrentBias());
  double yawError = yawModel.getCurrentBias();
  ignition::math::Quaterniond orientationError(
    ignition::math::Quaterniond(cos(yawError/2), 0.0, 0.0, sin(yawError/2)) *                                         // yaw error
    ignition::math::Quaterniond(1.0, 0.5 * accelDrift.Y() / gravity_length, 0.5 * -accelDrift.X() / gravity_length, 0.0)  // roll and pitch error
  );
#else
  math::Vector3 accelDrift = pose.rot.RotateVector(accelModel.getCurrentBias());
  double yawError = yawModel.getCurrentBias();
  math::Quaternion orientationError(
    math::Quaternion(cos(yawError/2), 0.0, 0.0, sin(yawError/2)) *                                         // yaw error
    math::Quaternion(1.0, 0.5 * accelDrift.y / gravity_length, 0.5 * -accelDrift.x / gravity_length, 0.0)  // roll and pitch error
  );
#endif
  orientationError.Normalize();
  rot = orientationError * rot;

  // copy data into pose message
  imuMsg.header.frame_id = frame_id_;
  imuMsg.header.stamp.sec = cur_time.sec;
  imuMsg.header.stamp.nsec = cur_time.nsec;

  // orientation quaternion
#if (GAZEBO_MAJOR_VERSION >= 8)
  imuMsg.orientation.x = rot.X();
  imuMsg.orientation.y = rot.Y();
  imuMsg.orientation.z = rot.Z();
  imuMsg.orientation.w = rot.W();
#else
  imuMsg.orientation.x = rot.x;
  imuMsg.orientation.y = rot.y;
  imuMsg.orientation.z = rot.z;
  imuMsg.orientation.w = rot.w;
#endif

  // pass angular rates
#if (GAZEBO_MAJOR_VERSION >= 8)
  imuMsg.angular_velocity.x    = rate.X();
  imuMsg.angular_velocity.y    = rate.Y();
  imuMsg.angular_velocity.z    = rate.Z();
#else
  imuMsg.angular_velocity.x    = rate.x;
  imuMsg.angular_velocity.y    = rate.y;
  imuMsg.angular_velocity.z    = rate.z;
#endif

  // pass accelerations
#if (GAZEBO_MAJOR_VERSION >= 8)
  imuMsg.linear_acceleration.x    = accel.X();
  imuMsg.linear_acceleration.y    = accel.Y();
  imuMsg.linear_acceleration.z    = accel.Z();
#else
  imuMsg.linear_acceleration.x    = accel.x;
  imuMsg.linear_acceleration.y    = accel.y;
  imuMsg.linear_acceleration.z    = accel.z;
#endif

  // fill in covariance matrix
  imuMsg.orientation_covariance[8] = yawModel.gaussian_noise*yawModel.gaussian_noise;
  if (gravity_length > 0.0) {
#if (GAZEBO_MAJOR_VERSION >= 8)
    imuMsg.orientation_covariance[0] = accelModel.gaussian_noise.X()*accelModel.gaussian_noise.X()/(gravity_length*gravity_length);
    imuMsg.orientation_covariance[4] = accelModel.gaussian_noise.Y()*accelModel.gaussian_noise.Y()/(gravity_length*gravity_length);
#else
    imuMsg.orientation_covariance[0] = accelModel.gaussian_noise.x*accelModel.gaussian_noise.x/(gravity_length*gravity_length);
    imuMsg.orientation_covariance[4] = accelModel.gaussian_noise.y*accelModel.gaussian_noise.y/(gravity_length*gravity_length);
#endif
  } else {
    imuMsg.orientation_covariance[0] = -1;
    imuMsg.orientation_covariance[4] = -1;
  }

  // publish to ros
  pub_.publish(imuMsg);
  ROS_DEBUG_NAMED("gazebo_ros_imu", "Publishing IMU data at t = %f", cur_time.Double());

  // publish bias
  if (bias_pub_) {
    biasMsg.header = imuMsg.header;
#if (GAZEBO_MAJOR_VERSION >= 8)
    biasMsg.orientation.x = orientationError.X();
    biasMsg.orientation.y = orientationError.Y();
    biasMsg.orientation.z = orientationError.Z();
    biasMsg.orientation.w = orientationError.W();
    biasMsg.angular_velocity.x = rateModel.getCurrentBias().X();
    biasMsg.angular_velocity.y = rateModel.getCurrentBias().Y();
    biasMsg.angular_velocity.z = rateModel.getCurrentBias().Z();
    biasMsg.linear_acceleration.x = accelModel.getCurrentBias().X();
    biasMsg.linear_acceleration.y = accelModel.getCurrentBias().Y();
    biasMsg.linear_acceleration.z = accelModel.getCurrentBias().Z();
#else
    biasMsg.orientation.x = orientationError.x;
    biasMsg.orientation.y = orientationError.y;
    biasMsg.orientation.z = orientationError.z;
    biasMsg.orientation.w = orientationError.w;
    biasMsg.angular_velocity.x = rateModel.getCurrentBias().x;
    biasMsg.angular_velocity.y = rateModel.getCurrentBias().y;
    biasMsg.angular_velocity.z = rateModel.getCurrentBias().z;
    biasMsg.linear_acceleration.x = accelModel.getCurrentBias().x;
    biasMsg.linear_acceleration.y = accelModel.getCurrentBias().y;
    biasMsg.linear_acceleration.z = accelModel.getCurrentBias().z;
#endif
    bias_pub_.publish(biasMsg);
  }

  // debug output
#ifdef DEBUG_OUTPUT
  if (debugPublisher) {
    geometry_msgs::PoseStamped debugPose;
    debugPose.header = imuMsg.header;
    debugPose.header.frame_id = "/map";
    debugPose.pose.orientation.w = imuMsg.orientation.w;
    debugPose.pose.orientation.x = imuMsg.orientation.x;
    debugPose.pose.orientation.y = imuMsg.orientation.y;
    debugPose.pose.orientation.z = imuMsg.orientation.z;
#if (GAZEBO_MAJOR_VERSION >= 8)
    ignition::math::Pose3d pose = link->WorldPose();
    debugPose.pose.position.x = pose.Pos().X();
    debugPose.pose.position.y = pose.Pos().Y();
    debugPose.pose.position.z = pose.Pos().Z();
#else
    math::Pose pose = link->GetWorldPose();
    debugPose.pose.position.x = pose.pos.x;
    debugPose.pose.position.y = pose.pos.y;
    debugPose.pose.position.z = pose.pos.z;
#endif
    debugPublisher.publish(debugPose);
  }
#endif // DEBUG_OUTPUT
}

#ifdef USE_CBQ
void GazeboRosIMU::CallbackQueueThread()
{
  static const double timeout = 0.01;

  while (rosnode_->ok())
  {
    callback_queue_.callAvailable(ros::WallDuration(timeout));
  }
}
#endif

// Register this plugin with the simulator
GZ_REGISTER_MODEL_PLUGIN(GazeboRosIMU)

} // namespace gazebo
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_gazebo/hector_gazebo_plugins/src/servo_plugin.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_gazebo_plugins/servo_plugin.h>
#include <gazebo/common/Events.hh>
#include <gazebo/physics/physics.hh>
#include <gazebo/gazebo_config.h>

#include <sensor_msgs/JointState.h>

#if (GAZEBO_MAJOR_VERSION > 1) || (GAZEBO_MINOR_VERSION >= 2)
  #define RADIAN Radian
#else
  #define RADIAN GetAsRadian
#endif

namespace gazebo {

GZ_REGISTER_MODEL_PLUGIN(ServoPlugin)

enum
{
  FIRST = 0, SECOND = 1, THIRD = 2
};

enum
{
  xyz, zyx
};

// Constructor
ServoPlugin::ServoPlugin()
{
  rosnode_ = 0;
  transform_listener_ = 0;
}

// Destructor
ServoPlugin::~ServoPlugin()
{
#if (GAZEBO_MAJOR_VERSION >= 8)
  updateConnection.reset();
#else
  event::Events::DisconnectWorldUpdateBegin(updateConnection);
#endif
  delete transform_listener_;
  rosnode_->shutdown();
  delete rosnode_;
}

// Load the controller
void ServoPlugin::Load(physics::ModelPtr _model, sdf::ElementPtr _sdf)
{
  // Get the world name.
  world = _model->GetWorld();

  // default parameters
  topicName = "drive";
  jointStateName = "joint_states";
  robotNamespace.clear();
  controlPeriod = 0;
  proportionalControllerGain = 8.0;
  derivativeControllerGain = 0.0;
  maximumVelocity = 0.0;
  maximumTorque = 0.0;

  // load parameters
  if (_sdf->HasElement("robotNamespace")) robotNamespace = _sdf->Get<std::string>("robotNamespace");
  if (_sdf->HasElement("topicName")) topicName = _sdf->Get<std::string>("topicName");
  if (_sdf->HasElement("jointStateName")) jointStateName = _sdf->Get<std::string>("jointStateName");
  if (_sdf->HasElement("firstServoName")) servo[FIRST].name = _sdf->Get<std::string>("firstServoName");
#if (GAZEBO_MAJOR_VERSION >= 8)
  if (_sdf->HasElement("firstServoAxis")) servo[FIRST].axis = _sdf->Get<ignition::math::Vector3d>("firstServoAxis");
#else
  if (_sdf->HasElement("firstServoAxis")) servo[FIRST].axis = _sdf->Get<math::Vector3>("firstServoAxis");
#endif
  if (_sdf->HasElement("secondServoName")) servo[SECOND].name = _sdf->Get<std::string>("secondServoName");
#if (GAZEBO_MAJOR_VERSION >= 8)
  if (_sdf->HasElement("secondServoAxis")) servo[SECOND].axis = _sdf->Get<ignition::math::Vector3d>("secondServoAxis");
#else
  if (_sdf->HasElement("secondServoAxis")) servo[SECOND].axis = _sdf->Get<math::Vector3>("secondServoAxis");
#endif
  if (_sdf->HasElement("thirdServoName")) servo[THIRD].name = _sdf->Get<std::string>("thirdServoName");
#if (GAZEBO_MAJOR_VERSION >= 8)
  if (_sdf->HasElement("thirdServoAxis")) servo[THIRD].axis = _sdf->Get<ignition::math::Vector3d>("thirdServoAxis");
#else
  if (_sdf->HasElement("thirdServoAxis")) servo[THIRD].axis = _sdf->Get<math::Vector3>("thirdServoAxis");
#endif
  if (_sdf->HasElement("proportionalControllerGain")) proportionalControllerGain = _sdf->Get<double>("proportionalControllerGain");
  if (_sdf->HasElement("derivativeControllerGain")) derivativeControllerGain = _sdf->Get<double>("derivativeControllerGain");
  if (_sdf->HasElement("maxVelocity")) maximumVelocity = _sdf->Get<double>("maxVelocity");
  if (_sdf->HasElement("torque")) maximumTorque = _sdf->Get<double>("torque");

  double controlRate = 0.0;
  if (_sdf->HasElement("controlRate")) controlRate = _sdf->Get<double>("controlRate");
  controlPeriod = controlRate > 0.0 ? 1.0/controlRate : 0.0;

  servo[FIRST].joint  = _model->GetJoint(servo[FIRST].name);
  servo[SECOND].joint = _model->GetJoint(servo[SECOND].name);
  servo[THIRD].joint  = _model->GetJoint(servo[THIRD].name);

  if (!servo[FIRST].joint)
    gzthrow("The controller couldn't get first joint");

  countOfServos = 1;
  if (servo[SECOND].joint) {
    countOfServos = 2;
    if (servo[THIRD].joint) {
      countOfServos = 3;
    }
  }
  else {
    if (servo[THIRD].joint) {
      gzthrow("The controller couldn't get second joint, but third joint was loaded");
    }
  }

  if (!ros::isInitialized()) {
    int argc = 0;
    char** argv = NULL;
    ros::init(argc, argv, "gazebo", ros::init_options::NoSigintHandler | ros::init_options::AnonymousName);
  }

  rosnode_ = new ros::NodeHandle(robotNamespace);

  transform_listener_ = new tf::TransformListener();
  transform_listener_->setExtrapolationLimit(ros::Duration(1.0));

  if (!topicName.empty()) {
    ros::SubscribeOptions so = ros::SubscribeOptions::create<geometry_msgs::QuaternionStamped>(topicName, 1,
                                                                                               boost::bind(&ServoPlugin::cmdCallback, this, _1),
                                                                                               ros::VoidPtr(), &queue_);
    sub_ = rosnode_->subscribe(so);
  }

  if (!jointStateName.empty()) {
    jointStatePub_ = rosnode_->advertise<sensor_msgs::JointState>(jointStateName, 10);
  }

  joint_state.header.frame_id = transform_listener_->resolve(_model->GetLink()->GetName());

  // New Mechanism for Updating every World Cycle
  // Listen to the update event. This event is broadcast every
  // simulation iteration.
  updateConnection = event::Events::ConnectWorldUpdateBegin(
      boost::bind(&ServoPlugin::Update, this));
}

// Initialize the controller
void ServoPlugin::Init()
{
  Reset();
}

// Reset
void ServoPlugin::Reset()
{
  // Reset orientation
  current_cmd.reset();

  enableMotors = true;

  servo[FIRST].velocity = 0;
  servo[SECOND].velocity = 0;
  servo[THIRD].velocity = 0;

#if (GAZEBO_MAJOR_VERSION >= 8)
  prevUpdateTime = world->SimTime();
#else
  prevUpdateTime = world->GetSimTime();
#endif
}

// Update the controller
void ServoPlugin::Update()
{
  // handle callbacks
  queue_.callAvailable();

  common::Time stepTime;
#if (GAZEBO_MAJOR_VERSION >= 8)
  stepTime = world->SimTime() - prevUpdateTime;
#else
  stepTime = world->GetSimTime() - prevUpdateTime;
#endif

  if (controlPeriod == 0.0 || stepTime > controlPeriod) {
    CalculateVelocities();
    publish_joint_states();
#if (GAZEBO_MAJOR_VERSION >= 8)
    prevUpdateTime = world->SimTime();
#else
    prevUpdateTime = world->GetSimTime();
#endif
  }

  if (enableMotors)
  {
    servo[FIRST].joint->SetVelocity(0, servo[FIRST].velocity);
    if (countOfServos > 1) {
      servo[SECOND].joint->SetVelocity(0, servo[SECOND].velocity);
      if (countOfServos > 2) {
        servo[THIRD].joint->SetVelocity(0, servo[THIRD].velocity);
      }
    }

#if (GAZEBO_MAJOR_VERSION > 4)
    servo[FIRST].joint->SetEffortLimit(0, maximumTorque);
    if (countOfServos > 1) {
      servo[SECOND].joint->SetEffortLimit(0, maximumTorque);
      if (countOfServos > 2) {
        servo[THIRD].joint->SetEffortLimit(0, maximumTorque);
      }
    }
#else
    servo[FIRST].joint->SetMaxForce(0, maximumTorque);
    if (countOfServos > 1) {
      servo[SECOND].joint->SetMaxForce(0, maximumTorque);
      if (countOfServos > 2) {
        servo[THIRD].joint->SetMaxForce(0, maximumTorque);
      }
    }
#endif
  } else {
#if (GAZEBO_MAJOR_VERSION > 4)
    servo[FIRST].joint->SetEffortLimit(0, 0.0);
    if (countOfServos > 1) {
      servo[SECOND].joint->SetEffortLimit(0, 0.0);
      if (countOfServos > 2) {
        servo[THIRD].joint->SetEffortLimit(0, 0.0);
      }
    }
#else
    servo[FIRST].joint->SetMaxForce(0, 0.0);
    if (countOfServos > 1) {
      servo[SECOND].joint->SetMaxForce(0, 0.0);
      if (countOfServos > 2) {
        servo[THIRD].joint->SetMaxForce(0, 0.0);
      }
    }
#endif
  }
}

void ServoPlugin::CalculateVelocities()
{
  tf::StampedTransform transform;
  boost::mutex::scoped_lock lock(mutex);

  if(!current_cmd){
    geometry_msgs::QuaternionStamped *default_cmd = new geometry_msgs::QuaternionStamped;
    default_cmd->header.frame_id = "base_stabilized";
    default_cmd->quaternion.w = 1;
    current_cmd.reset(default_cmd);
  }

  try{
    // ros::Time simTime(world->GetSimTime().sec, world->GetSimTime().nsec);
    transform_listener_->lookupTransform("base_link", current_cmd->header.frame_id, ros::Time(0), transform);
  }
  catch (tf::TransformException ex){
    ROS_DEBUG("%s",ex.what());
    servo[FIRST].velocity = 0.0;
    servo[SECOND].velocity = 0.0;
    servo[THIRD].velocity = 0.0;
    return;
  }

  rotation_.Set(current_cmd->quaternion.w, current_cmd->quaternion.x, current_cmd->quaternion.y, current_cmd->quaternion.z);

#if (GAZEBO_MAJOR_VERSION >= 8)
  ignition::math::Quaterniond quat(transform.getRotation().getW(),transform.getRotation().getX(),transform.getRotation().getY(),transform.getRotation().getZ());
#else
  math::Quaternion quat(transform.getRotation().getW(),transform.getRotation().getX(),transform.getRotation().getY(),transform.getRotation().getZ());
#endif

  rotation_ = quat * rotation_;

  double temp[5];
  double desAngle[3];
  double actualAngle[3] = {0.0, 0.0, 0.0};
  double actualVel[3] = {0.0, 0.0, 0.0};

  //TODO use countOfServos for calculation
  rotationConv = 99;
  orderOfAxes[0] = 99;
  orderOfAxes[1] = 99;
  orderOfAxes[2] = 99;

  switch(countOfServos) {
    case 2:
#if (GAZEBO_MAJOR_VERSION >= 8)
      if ((servo[FIRST].axis.Z() == 1) && (servo[SECOND].axis.Y() == 1)) {
#else
      if ((servo[FIRST].axis.z == 1) && (servo[SECOND].axis.y == 1)) {
#endif
        rotationConv = zyx;
        orderOfAxes[0] = 0;
        orderOfAxes[1] = 1;
      }
      else {
#if (GAZEBO_MAJOR_VERSION >= 8)
        if ((servo[FIRST].axis.X() == 1) && (servo[SECOND].axis.Y() == 1)) {
#else
        if ((servo[FIRST].axis.x == 1) && (servo[SECOND].axis.y == 1)) {
#endif
          rotationConv = xyz;
          orderOfAxes[0] = 0;
          orderOfAxes[1] = 1;
        }
      }
      break;

    case 3:
#if (GAZEBO_MAJOR_VERSION >= 8)
      if ((servo[FIRST].axis.Z() == 1) && (servo[SECOND].axis.Y() == 1) && (servo[THIRD].axis.X() == 1)) {
#else
      if ((servo[FIRST].axis.z == 1) && (servo[SECOND].axis.y == 1) && (servo[THIRD].axis.x == 1)) {
#endif
        rotationConv = zyx;
        orderOfAxes[0] = 0;
        orderOfAxes[1] = 1;
        orderOfAxes[2] = 2;
      }
#if (GAZEBO_MAJOR_VERSION >= 8)
      else if ((servo[FIRST].axis.X() == 1) && (servo[SECOND].axis.Y() == 1) && (servo[THIRD].axis.Z() == 1)) {
#else
      else if ((servo[FIRST].axis.x == 1) && (servo[SECOND].axis.y == 1) && (servo[THIRD].axis.z == 1)) {
#endif
        rotationConv = xyz;
        orderOfAxes[0] = 0;
        orderOfAxes[1] = 1;
        orderOfAxes[2] = 2;
      }
      break;

    case 1:
#if (GAZEBO_MAJOR_VERSION >= 8)
      if (servo[FIRST].axis.Y() == 1) {
#else
      if (servo[FIRST].axis.y == 1) {
#endif
         rotationConv = xyz;
         orderOfAxes[0] = 1;
      }
      break;

    default:
      gzthrow("Something went wrong. The count of servos is greater than 3");
      break;
  }

  switch(rotationConv)  {
    case zyx:
#if (GAZEBO_MAJOR_VERSION >= 8)
      temp[0] =  2*(rotation_.X()*rotation_.Y() + rotation_.W()*rotation_.Z());
      temp[1] =     rotation_.W()*rotation_.W() + rotation_.X()*rotation_.X() - rotation_.Y()*rotation_.Y() - rotation_.Z()*rotation_.Z();
      temp[2] = -2*(rotation_.X()*rotation_.Z() - rotation_.W()*rotation_.Y());
      temp[3] =  2*(rotation_.Y()*rotation_.Z() + rotation_.W()*rotation_.X());
      temp[4] =     rotation_.W()*rotation_.W() - rotation_.X()*rotation_.X() - rotation_.Y()*rotation_.Y() + rotation_.Z()*rotation_.Z();
#else
      temp[0] =  2*(rotation_.x*rotation_.y + rotation_.w*rotation_.z);
      temp[1] =     rotation_.w*rotation_.w + rotation_.x*rotation_.x - rotation_.y*rotation_.y - rotation_.z*rotation_.z;
      temp[2] = -2*(rotation_.x*rotation_.z - rotation_.w*rotation_.y);
      temp[3] =  2*(rotation_.y*rotation_.z + rotation_.w*rotation_.x);
      temp[4] =     rotation_.w*rotation_.w - rotation_.x*rotation_.x - rotation_.y*rotation_.y + rotation_.z*rotation_.z;
#endif
      break;

    case xyz:
#if (GAZEBO_MAJOR_VERSION >= 8)
      temp[0] =  -2*(rotation_.Y()*rotation_.Z() - rotation_.W()*rotation_.X());
      temp[1] =      rotation_.W()*rotation_.W() - rotation_.X()*rotation_.X() - rotation_.Y()*rotation_.Y() + rotation_.Z()*rotation_.Z();
      temp[2] =   2*(rotation_.X()*rotation_.Z() + rotation_.W()*rotation_.Y());
      temp[3] =  -2*(rotation_.X()*rotation_.Y() - rotation_.W()*rotation_.Z());
      temp[4] =      rotation_.W()*rotation_.W() + rotation_.X()*rotation_.X() - rotation_.Y()*rotation_.Y() - rotation_.Z()*rotation_.Z();
#else
      temp[0] =  -2*(rotation_.y*rotation_.z - rotation_.w*rotation_.x);
      temp[1] =      rotation_.w*rotation_.w - rotation_.x*rotation_.x - rotation_.y*rotation_.y + rotation_.z*rotation_.z;
      temp[2] =   2*(rotation_.x*rotation_.z + rotation_.w*rotation_.y);
      temp[3] =  -2*(rotation_.x*rotation_.y - rotation_.w*rotation_.z);
      temp[4] =      rotation_.w*rotation_.w + rotation_.x*rotation_.x - rotation_.y*rotation_.y - rotation_.z*rotation_.z;
#endif
      break;

    default:
      gzthrow("joint order " << rotationConv << " not supported");
      break;
  }

  desAngle[0] = atan2(temp[0], temp[1]);
  desAngle[1] = asin(temp[2]);
  desAngle[2] = atan2(temp[3], temp[4]);

#if (GAZEBO_MAJOR_VERSION >= 8)
  actualAngle[FIRST] = servo[FIRST].joint->Position(0);
#else
  actualAngle[FIRST] = servo[FIRST].joint->GetAngle(0).RADIAN();
#endif
  actualVel[FIRST] = servo[FIRST].joint->GetVelocity(0);
  ROS_DEBUG_NAMED("servo_plugin", "%s servo angle: %f - %f", servo[FIRST].name.c_str(), desAngle[orderOfAxes[FIRST]], actualAngle[FIRST]);
  servo[FIRST].velocity = ( proportionalControllerGain*(desAngle[orderOfAxes[FIRST]] - actualAngle[FIRST]) - derivativeControllerGain*actualVel[FIRST]);
  if (maximumVelocity > 0.0 && fabs(servo[FIRST].velocity) > maximumVelocity) servo[FIRST].velocity = (servo[FIRST].velocity > 0 ? maximumVelocity : -maximumVelocity);

  if (countOfServos > 1) {
#if (GAZEBO_MAJOR_VERSION >= 8)
    actualAngle[SECOND] = servo[SECOND].joint->Position(0);
#else
    actualAngle[SECOND] = servo[SECOND].joint->GetAngle(0).RADIAN();
#endif
    actualVel[SECOND] = servo[SECOND].joint->GetVelocity(0);
    ROS_DEBUG_NAMED("servo_plugin", "%s servo angle: %f - %f", servo[SECOND].name.c_str(), desAngle[orderOfAxes[SECOND]], actualAngle[SECOND]);
    servo[SECOND].velocity = ( proportionalControllerGain*(desAngle[orderOfAxes[SECOND]] - actualAngle[SECOND]) - derivativeControllerGain*actualVel[SECOND]);
    if (maximumVelocity > 0.0 && fabs(servo[SECOND].velocity) > maximumVelocity) servo[SECOND].velocity = (servo[SECOND].velocity > 0 ? maximumVelocity : -maximumVelocity);

    if (countOfServos == 3) {
#if (GAZEBO_MAJOR_VERSION >= 8)
      actualAngle[THIRD] = servo[THIRD].joint->Position(0);
#else
      actualAngle[THIRD] = servo[THIRD].joint->GetAngle(0).RADIAN();
#endif
      actualVel[THIRD] = servo[THIRD].joint->GetVelocity(0);
      ROS_DEBUG_NAMED("servo_plugin", "%s servo angle: %f - %f", servo[THIRD].name.c_str(), desAngle[orderOfAxes[THIRD]], actualAngle[THIRD]);
      servo[THIRD].velocity = ( proportionalControllerGain*(desAngle[orderOfAxes[THIRD]] - actualAngle[THIRD]) - derivativeControllerGain*actualVel[THIRD]);
      if (maximumVelocity > 0.0 && fabs(servo[THIRD].velocity) > maximumVelocity) servo[THIRD].velocity = (servo[THIRD].velocity > 0 ? maximumVelocity : -maximumVelocity);
    }
  }

  // Changed motors to be always on, which is probably what we want anyway
  enableMotors = true; //myIface->data->cmdEnableMotors > 0;
}

// NEW: Store the velocities from the ROS message
void ServoPlugin::cmdCallback(const geometry_msgs::QuaternionStamped::ConstPtr& cmd_msg)
{
  boost::mutex::scoped_lock lock(mutex);
  current_cmd = cmd_msg;
}

void ServoPlugin::publish_joint_states()
{
  if (!jointStatePub_) return;

#if (GAZEBO_MAJOR_VERSION >= 8)
  joint_state.header.stamp.sec = (world->SimTime()).sec;
  joint_state.header.stamp.nsec = (world->SimTime()).nsec;
#else
  joint_state.header.stamp.sec = (world->GetSimTime()).sec;
  joint_state.header.stamp.nsec = (world->GetSimTime()).nsec;
#endif

  joint_state.name.resize(countOfServos);
  joint_state.position.resize(countOfServos);
  joint_state.velocity.resize(countOfServos);
  joint_state.effort.resize(countOfServos);

  for (unsigned int i = 0; i < countOfServos; i++) {
    joint_state.name[i] = servo[i].joint->GetName();
#if (GAZEBO_MAJOR_VERSION >= 8)
    joint_state.position[i] = servo[i].joint->Position(0);
#else
    joint_state.position[i] = servo[i].joint->GetAngle(0).RADIAN();
#endif
    joint_state.velocity[i] = servo[i].joint->GetVelocity(0);
    joint_state.effort[i] = servo[i].joint->GetForce(0u);
  }

  jointStatePub_.publish(joint_state);
}

} // namespace gazebo
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_gazebo/hector_gazebo_plugins/src/gazebo_ros_force_based_move.cpp
/*
 * Copyright 2015 Stefan Kohlbrecher, TU Darmstadt
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
*/

/*
 * Desc: Simple model controller that uses a twist message to exert
 *       forces on a robot, resulting in motion. Based on the
 *       planar_move plugin by Piyush Khandelwal.
 * Author: Stefan Kohlbrecher
 * Date: 06 August 2015
 */

#include <hector_gazebo_plugins/gazebo_ros_force_based_move.h>

namespace gazebo 
{

  GazeboRosForceBasedMove::GazeboRosForceBasedMove() {}

  GazeboRosForceBasedMove::~GazeboRosForceBasedMove() {}

  // Load the controller
  void GazeboRosForceBasedMove::Load(physics::ModelPtr parent,
      sdf::ElementPtr sdf) 
  {

    parent_ = parent;

    /* Parse parameters */

    robot_namespace_ = "";
    if (!sdf->HasElement("robotNamespace")) 
    {
      ROS_INFO("PlanarMovePlugin missing <robotNamespace>, "
          "defaults to \"%s\"", robot_namespace_.c_str());
    }
    else 
    {
      robot_namespace_ = 
        sdf->GetElement("robotNamespace")->Get<std::string>();
    }

    command_topic_ = "cmd_vel";
    if (!sdf->HasElement("commandTopic")) 
    {
      ROS_WARN("PlanarMovePlugin (ns = %s) missing <commandTopic>, "
          "defaults to \"%s\"", 
          robot_namespace_.c_str(), command_topic_.c_str());
    } 
    else 
    {
      command_topic_ = sdf->GetElement("commandTopic")->Get<std::string>();
    }

    odometry_topic_ = "odom";
    if (!sdf->HasElement("odometryTopic")) 
    {
      ROS_WARN("PlanarMovePlugin (ns = %s) missing <odometryTopic>, "
          "defaults to \"%s\"", 
          robot_namespace_.c_str(), odometry_topic_.c_str());
    } 
    else 
    {
      odometry_topic_ = sdf->GetElement("odometryTopic")->Get<std::string>();
    }

    odometry_frame_ = "odom";
    if (!sdf->HasElement("odometryFrame")) 
    {
      ROS_WARN("PlanarMovePlugin (ns = %s) missing <odometryFrame>, "
          "defaults to \"%s\"",
          robot_namespace_.c_str(), odometry_frame_.c_str());
    }
    else 
    {
      odometry_frame_ = sdf->GetElement("odometryFrame")->Get<std::string>();
    }
    
    
    torque_yaw_velocity_p_gain_ = 100.0;
    force_x_velocity_p_gain_ = 10000.0;
    force_y_velocity_p_gain_ = 10000.0;
    
    if (sdf->HasElement("yaw_velocity_p_gain"))
      (sdf->GetElement("yaw_velocity_p_gain")->GetValue()->Get(torque_yaw_velocity_p_gain_));

    if (sdf->HasElement("x_velocity_p_gain"))
      (sdf->GetElement("x_velocity_p_gain")->GetValue()->Get(force_x_velocity_p_gain_));

    if (sdf->HasElement("y_velocity_p_gain"))
      (sdf->GetElement("y_velocity_p_gain")->GetValue()->Get(force_y_velocity_p_gain_));
      
    ROS_INFO_STREAM("ForceBasedMove using gains: yaw: " << torque_yaw_velocity_p_gain_ <<
                                                 " x: " << force_x_velocity_p_gain_ <<
                                                 " y: " << force_y_velocity_p_gain_ << "\n");

    robot_base_frame_ = "base_footprint";
    if (!sdf->HasElement("robotBaseFrame")) 
    {
      ROS_WARN("PlanarMovePlugin (ns = %s) missing <robotBaseFrame>, "
          "defaults to \"%s\"",
          robot_namespace_.c_str(), robot_base_frame_.c_str());
    } 
    else 
    {
      robot_base_frame_ = sdf->GetElement("robotBaseFrame")->Get<std::string>();
    }

    ROS_INFO_STREAM("robotBaseFrame for force based move plugin: " << robot_base_frame_  << "\n");

    this->link_ = parent->GetLink(robot_base_frame_);

    odometry_rate_ = 20.0;
    if (!sdf->HasElement("odometryRate")) 
    {
      ROS_WARN("PlanarMovePlugin (ns = %s) missing <odometryRate>, "
          "defaults to %f",
          robot_namespace_.c_str(), odometry_rate_);
    } 
    else 
    {
      odometry_rate_ = sdf->GetElement("odometryRate")->Get<double>();
    } 

    this->publish_odometry_tf_ = true;
    if (!sdf->HasElement("publishOdometryTf")) {
      ROS_WARN("PlanarMovePlugin Plugin (ns = %s) missing <publishOdometryTf>, defaults to %s",
               this->robot_namespace_.c_str(), this->publish_odometry_tf_ ? "true" : "false");
    } else {
      this->publish_odometry_tf_ = sdf->GetElement("publishOdometryTf")->Get<bool>();
    }
 
#if (GAZEBO_MAJOR_VERSION >= 8)
    last_odom_publish_time_ = parent_->GetWorld()->SimTime();
    last_odom_pose_ = parent_->WorldPose();
#else
    last_odom_publish_time_ = parent_->GetWorld()->GetSimTime();
    last_odom_pose_ = parent_->GetWorldPose();
#endif
    x_ = 0;
    y_ = 0;
    rot_ = 0;
    alive_ = true;

    odom_transform_.setIdentity();

    // Ensure that ROS has been initialized and subscribe to cmd_vel
    if (!ros::isInitialized()) 
    {
      ROS_FATAL_STREAM("PlanarMovePlugin (ns = " << robot_namespace_
        << "). A ROS node for Gazebo has not been initialized, "
        << "unable to load plugin. Load the Gazebo system plugin "
        << "'libgazebo_ros_api_plugin.so' in the gazebo_ros package)");
      return;
    }
    rosnode_.reset(new ros::NodeHandle(robot_namespace_));

    ROS_DEBUG("OCPlugin (%s) has started!", 
        robot_namespace_.c_str());

    tf_prefix_ = tf::getPrefixParam(*rosnode_);

    if (publish_odometry_tf_)
      transform_broadcaster_.reset(new tf::TransformBroadcaster());

    // subscribe to the odometry topic
    ros::SubscribeOptions so =
      ros::SubscribeOptions::create<geometry_msgs::Twist>(command_topic_, 1,
          boost::bind(&GazeboRosForceBasedMove::cmdVelCallback, this, _1),
          ros::VoidPtr(), &queue_);

    vel_sub_ = rosnode_->subscribe(so);
    odometry_pub_ = rosnode_->advertise<nav_msgs::Odometry>(odometry_topic_, 1);

    // start custom queue for diff drive
    callback_queue_thread_ = 
      boost::thread(boost::bind(&GazeboRosForceBasedMove::QueueThread, this));

    // listen to the update event (broadcast every simulation iteration)
    update_connection_ = 
      event::Events::ConnectWorldUpdateBegin(
          boost::bind(&GazeboRosForceBasedMove::UpdateChild, this));

  }

  // Update the controller
  void GazeboRosForceBasedMove::UpdateChild()
  {
    boost::mutex::scoped_lock scoped_lock(lock);
#if (GAZEBO_MAJOR_VERSION >= 8)
    ignition::math::Pose3d pose = parent_->WorldPose();

    ignition::math::Vector3d angular_vel = parent_->WorldAngularVel();

    double error = angular_vel.Z() - rot_;

    link_->AddTorque(ignition::math::Vector3d(0.0, 0.0, -error * torque_yaw_velocity_p_gain_));

    float yaw = pose.Rot().Yaw();

    ignition::math::Vector3d linear_vel = parent_->RelativeLinearVel();

    link_->AddRelativeForce(ignition::math::Vector3d((x_ - linear_vel.X())* force_x_velocity_p_gain_,
                                                     (y_ - linear_vel.Y())* force_y_velocity_p_gain_,
                                                     0.0));
#else
    math::Pose pose = parent_->GetWorldPose();

    math::Vector3 angular_vel = parent_->GetWorldAngularVel();

    double error = angular_vel.z - rot_;

    link_->AddTorque(math::Vector3(0.0, 0.0, -error * torque_yaw_velocity_p_gain_));

    float yaw = pose.rot.GetYaw();

    math::Vector3 linear_vel = parent_->GetRelativeLinearVel();

    link_->AddRelativeForce(math::Vector3((x_ - linear_vel.x)* force_x_velocity_p_gain_,
                                          (y_ - linear_vel.y)* force_y_velocity_p_gain_,
                                          0.0));
#endif
    //parent_->PlaceOnNearestEntityBelow();
    //parent_->SetLinearVel(math::Vector3(
    //      x_ * cosf(yaw) - y_ * sinf(yaw),
    //      y_ * cosf(yaw) + x_ * sinf(yaw),
    //      0));
    //parent_->SetAngularVel(math::Vector3(0, 0, rot_));

    if (odometry_rate_ > 0.0) {
#if (GAZEBO_MAJOR_VERSION >= 8)
      common::Time current_time = parent_->GetWorld()->SimTime();
#else
      common::Time current_time = parent_->GetWorld()->GetSimTime();
#endif
      double seconds_since_last_update = 
        (current_time - last_odom_publish_time_).Double();
      if (seconds_since_last_update > (1.0 / odometry_rate_)) {
        publishOdometry(seconds_since_last_update);
        last_odom_publish_time_ = current_time;
      }
    }
  }

  // Finalize the controller
  void GazeboRosForceBasedMove::FiniChild() {
    alive_ = false;
    queue_.clear();
    queue_.disable();
    rosnode_->shutdown();
    callback_queue_thread_.join();
  }

  void GazeboRosForceBasedMove::cmdVelCallback(
      const geometry_msgs::Twist::ConstPtr& cmd_msg) 
  {
    boost::mutex::scoped_lock scoped_lock(lock);
    x_ = cmd_msg->linear.x;
    y_ = cmd_msg->linear.y;
    rot_ = cmd_msg->angular.z;
  }

  void GazeboRosForceBasedMove::QueueThread()
  {
    static const double timeout = 0.01;
    while (alive_ && rosnode_->ok()) 
    {
      queue_.callAvailable(ros::WallDuration(timeout));
    }
  }

  void GazeboRosForceBasedMove::publishOdometry(double step_time)
  {

    ros::Time current_time = ros::Time::now();
    std::string odom_frame = tf::resolve(tf_prefix_, odometry_frame_);
    std::string base_footprint_frame = 
      tf::resolve(tf_prefix_, robot_base_frame_);

#if (GAZEBO_MAJOR_VERSION >= 8)
    ignition::math::Vector3d angular_vel = parent_->RelativeAngularVel();
    ignition::math::Vector3d linear_vel = parent_->RelativeLinearVel();

    odom_transform_= odom_transform_ * this->getTransformForMotion(linear_vel.X(), angular_vel.Z(), step_time);

    tf::poseTFToMsg(odom_transform_, odom_.pose.pose);
    odom_.twist.twist.angular.z = angular_vel.Z();
    odom_.twist.twist.linear.x  = linear_vel.X();
#else
    math::Vector3 angular_vel = parent_->GetRelativeAngularVel();
    math::Vector3 linear_vel = parent_->GetRelativeLinearVel();

    odom_transform_= odom_transform_ * this->getTransformForMotion(linear_vel.x, angular_vel.z, step_time);

    tf::poseTFToMsg(odom_transform_, odom_.pose.pose);
    odom_.twist.twist.angular.z = angular_vel.z;
    odom_.twist.twist.linear.x  = linear_vel.x;
#endif

    odom_.header.stamp = current_time;
    odom_.header.frame_id = odom_frame;
    odom_.child_frame_id = base_footprint_frame;

    if (transform_broadcaster_.get()){
      transform_broadcaster_->sendTransform(
          tf::StampedTransform(odom_transform_, current_time, odom_frame,
              base_footprint_frame));
    }
    
    odom_.pose.covariance[0] = 0.001;
    odom_.pose.covariance[7] = 0.001;
    odom_.pose.covariance[14] = 1000000000000.0;
    odom_.pose.covariance[21] = 1000000000000.0;
    odom_.pose.covariance[28] = 1000000000000.0;
    
#if (GAZEBO_MAJOR_VERSION >= 8)
    if (std::abs(angular_vel.Z()) < 0.0001) {
#else
    if (std::abs(angular_vel.z) < 0.0001) {
#endif
      odom_.pose.covariance[35] = 0.01;
    }else{
      odom_.pose.covariance[35] = 100.0;
    }

    odom_.twist.covariance[0] = 0.001;
    odom_.twist.covariance[7] = 0.001;
    odom_.twist.covariance[14] = 0.001;
    odom_.twist.covariance[21] = 1000000000000.0;
    odom_.twist.covariance[28] = 1000000000000.0;

#if (GAZEBO_MAJOR_VERSION >= 8)
    if (std::abs(angular_vel.Z()) < 0.0001) {
#else
    if (std::abs(angular_vel.z) < 0.0001) {
#endif
      odom_.twist.covariance[35] = 0.01;
    }else{
      odom_.twist.covariance[35] = 100.0;
    }



    odometry_pub_.publish(odom_);
  }


  tf::Transform GazeboRosForceBasedMove::getTransformForMotion(double linear_vel_x, double angular_vel, double timeSeconds) const
  {
    tf::Transform tmp;
    tmp.setIdentity();


    if (std::abs(angular_vel) < 0.0001) {
      //Drive straight
      tmp.setOrigin(tf::Vector3(static_cast<double>(linear_vel_x*timeSeconds), 0.0, 0.0));
    } else {
      //Follow circular arc
      double distChange = linear_vel_x * timeSeconds;
      double angleChange = angular_vel * timeSeconds;

      double arcRadius = distChange / angleChange;

      tmp.setOrigin(tf::Vector3(std::sin(angleChange) * arcRadius,
                                arcRadius - std::cos(angleChange) * arcRadius,
                                0.0));
      tmp.setRotation(tf::createQuaternionFromYaw(angleChange));
    }

    return tmp;
  }

  GZ_REGISTER_MODEL_PLUGIN(GazeboRosForceBasedMove)
}

=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_gazebo/hector_gazebo_plugins/src/gazebo_ros_sonar.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_gazebo_plugins/gazebo_ros_sonar.h>
#include <gazebo/common/Events.hh>
#include <gazebo/physics/physics.hh>
#include <gazebo/sensors/RaySensor.hh>

#include <limits>

#include <gazebo/gazebo_config.h>

namespace gazebo {

GazeboRosSonar::GazeboRosSonar()
{
}

////////////////////////////////////////////////////////////////////////////////
// Destructor
GazeboRosSonar::~GazeboRosSonar()
{
  updateTimer.Disconnect(updateConnection);
  sensor_->SetActive(false);

  dynamic_reconfigure_server_.reset();

  node_handle_->shutdown();
  delete node_handle_;
}

////////////////////////////////////////////////////////////////////////////////
// Load the controller
void GazeboRosSonar::Load(sensors::SensorPtr _sensor, sdf::ElementPtr _sdf)
{
  // Get then name of the parent sensor
#if (GAZEBO_MAJOR_VERSION > 6)
  sensor_ = std::dynamic_pointer_cast<sensors::RaySensor>(_sensor);
#else
  sensor_ = boost::dynamic_pointer_cast<sensors::RaySensor>(_sensor);
#endif
  if (!sensor_)
  {
    gzthrow("GazeboRosSonar requires a Ray Sensor as its parent");
    return;
  }

  // Get the world name.
#if (GAZEBO_MAJOR_VERSION > 6)
  std::string worldName = sensor_->WorldName();
#else
  std::string worldName = sensor_->GetWorldName();
#endif
  world = physics::get_world(worldName);

  // default parameters
  namespace_.clear();
  topic_ = "sonar";
  frame_id_ = "/sonar_link";

  // load parameters
  if (_sdf->HasElement("robotNamespace"))
    namespace_ = _sdf->GetElement("robotNamespace")->GetValue()->GetAsString();

  if (_sdf->HasElement("frameId"))
    frame_id_ = _sdf->GetElement("frameId")->GetValue()->GetAsString();

  if (_sdf->HasElement("topicName"))
    topic_ = _sdf->GetElement("topicName")->GetValue()->GetAsString();

  sensor_model_.Load(_sdf);

  range_.header.frame_id = frame_id_;
  range_.radiation_type = sensor_msgs::Range::ULTRASOUND;
#if (GAZEBO_MAJOR_VERSION > 6)
  range_.field_of_view = std::min(fabs((sensor_->AngleMax() - sensor_->AngleMin()).Radian()), fabs((sensor_->VerticalAngleMax() - sensor_->VerticalAngleMin()).Radian()));
  range_.max_range = sensor_->RangeMax();
  range_.min_range = sensor_->RangeMin();
#else
  range_.field_of_view = std::min(fabs((sensor_->GetAngleMax() - sensor_->GetAngleMin()).Radian()), fabs((sensor_->GetVerticalAngleMax() - sensor_->GetVerticalAngleMin()).Radian()));
  range_.max_range = sensor_->GetRangeMax();
  range_.min_range = sensor_->GetRangeMin();
#endif

  // Make sure the ROS node for Gazebo has already been initialized
  if (!ros::isInitialized())
  {
    ROS_FATAL_STREAM("A ROS node for Gazebo has not been initialized, unable to load plugin. "
      << "Load the Gazebo system plugin 'libgazebo_ros_api_plugin.so' in the gazebo_ros package)");
    return;
  }

  node_handle_ = new ros::NodeHandle(namespace_);
  publisher_ = node_handle_->advertise<sensor_msgs::Range>(topic_, 1);

  // setup dynamic_reconfigure server
  dynamic_reconfigure_server_.reset(new dynamic_reconfigure::Server<SensorModelConfig>(ros::NodeHandle(*node_handle_, topic_)));
  dynamic_reconfigure_server_->setCallback(boost::bind(&SensorModel::dynamicReconfigureCallback, &sensor_model_, _1, _2));

  Reset();

  // connect Update function
  updateTimer.setUpdateRate(10.0);
  updateTimer.Load(world, _sdf);
  updateConnection = updateTimer.Connect(boost::bind(&GazeboRosSonar::Update, this));

  // activate RaySensor
  sensor_->SetActive(true);
}

void GazeboRosSonar::Reset()
{
  updateTimer.Reset();
  sensor_model_.reset();
}

////////////////////////////////////////////////////////////////////////////////
// Update the controller
void GazeboRosSonar::Update()
{
#if (GAZEBO_MAJOR_VERSION >= 8)
  common::Time sim_time = world->SimTime();
#else
  common::Time sim_time = world->GetSimTime();
#endif
  double dt = updateTimer.getTimeSinceLastUpdate().Double();

  // activate RaySensor if it is not yet active
  if (!sensor_->IsActive()) sensor_->SetActive(true);

#if (GAZEBO_MAJOR_VERSION >= 8)
  range_.header.stamp.sec  = (world->SimTime()).sec;
  range_.header.stamp.nsec = (world->SimTime()).nsec;
#else
  range_.header.stamp.sec  = (world->GetSimTime()).sec;
  range_.header.stamp.nsec = (world->GetSimTime()).nsec;
#endif

  // find ray with minimal range
  range_.range = std::numeric_limits<sensor_msgs::Range::_range_type>::max();
#if (GAZEBO_MAJOR_VERSION > 6)
  int num_ranges = sensor_->LaserShape()->GetSampleCount() * sensor_->LaserShape()->GetVerticalSampleCount();
#else
  int num_ranges = sensor_->GetLaserShape()->GetSampleCount() * sensor_->GetLaserShape()->GetVerticalSampleCount();
#endif
  for(int i = 0; i < num_ranges; ++i) {
#if (GAZEBO_MAJOR_VERSION > 6)
    double ray = sensor_->LaserShape()->GetRange(i);
#else
    double ray = sensor_->GetLaserShape()->GetRange(i);
#endif
    if (ray < range_.range) range_.range = ray;
  }

  // add Gaussian noise (and limit to min/max range)
  if (range_.range < range_.max_range) {
    range_.range = sensor_model_(range_.range, dt);
    if (range_.range < range_.min_range) range_.range = range_.min_range;
    if (range_.range > range_.max_range) range_.range = range_.max_range;
  }

  publisher_.publish(range_);
}

// Register this plugin with the simulator
GZ_REGISTER_SENSOR_PLUGIN(GazeboRosSonar)

} // namespace gazebo
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_gazebo/hector_gazebo_plugins/src/diffdrive_plugin_multi_wheel.cpp
/*
    Copyright (c) 2010, Daniel Hewlett, Antons Rebguns
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the <organization> nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY Antons Rebguns <email> ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL Antons Rebguns <email> BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/*
 * \file  gazebo_ros_diff_drive.cpp
 *
 * \brief A differential drive plugin for gazebo. Based on the diffdrive plugin 
 * developed for the erratic robot (see copyright notice above). The original
 * plugin can be found in the ROS package gazebo_erratic_plugins.
 *
 * \author  Piyush Khandelwal (piyushk@gmail.com)
 *
 * $ Id: 06/21/2013 11:23:40 AM piyushk $
 */


/**
 * A diff drive plugin supporting multiple wheels per vehicle side. Based on
 * existing plugins as stated above this notice.
 */

/*
    Copyright (c) 2014, Stefan Kohlbrecher
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the <organization> nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY Antons Rebguns <email> ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL Antons Rebguns <email> BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include <algorithm>
#include <assert.h>

#include <hector_gazebo_plugins/diffdrive_plugin_multi_wheel.h>

#if (GAZEBO_MAJOR_VERSION < 8)
#include <gazebo/math/gzmath.hh>
#endif
#include <sdf/sdf.hh>

#include <ros/ros.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_listener.h>
#include <geometry_msgs/Twist.h>
#include <nav_msgs/GetMap.h>
#include <nav_msgs/Odometry.h>
#include <boost/bind.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <gazebo/gazebo_config.h>

namespace gazebo {

  enum {
    RIGHT,
    LEFT,
  };

  GazeboRosDiffDriveMultiWheel::GazeboRosDiffDriveMultiWheel() {}

  // Destructor
  GazeboRosDiffDriveMultiWheel::~GazeboRosDiffDriveMultiWheel() {
    delete rosnode_;
    delete transform_broadcaster_;
  }

  // Load the controller
  void GazeboRosDiffDriveMultiWheel::Load(physics::ModelPtr _parent, sdf::ElementPtr _sdf) {

    this->parent = _parent;
    this->world = _parent->GetWorld();

    this->robot_namespace_ = "";
    if (!_sdf->HasElement("robotNamespace")) {
      ROS_INFO("GazeboRosDiffDriveMultiWheel Plugin missing <robotNamespace>, defaults to \"%s\"", 
          this->robot_namespace_.c_str());
    } else {
      this->robot_namespace_ = 
        _sdf->GetElement("robotNamespace")->Get<std::string>() + "/";
    }

    //this->left_joint_names_ = "left_joint";
    if (!_sdf->HasElement("leftJoints")) {
      gzthrow("Have to specify space separated left side joint names via <leftJoints> tag!");
    } else {
      std::string joint_string = _sdf->GetElement("leftJoints")->Get<std::string>();
      boost::split( joint_names_[LEFT], joint_string, boost::is_any_of(" ") );
    }

    //this->right_joint_names_ = "right_joint";
    if (!_sdf->HasElement("rightJoints")) {
      gzthrow("Have to specify space separated right side joint names via <rightJoints> tag!");
    } else {
      std::string joint_string = _sdf->GetElement("rightJoints")->Get<std::string>();
      boost::split( joint_names_[RIGHT], joint_string, boost::is_any_of(" ") );
    }

    this->wheel_separation_ = 0.34;
    if (!_sdf->HasElement("wheelSeparation")) {
      ROS_WARN("GazeboRosDiffDriveMultiWheel Plugin (ns = %s) missing <wheelSeparation>, defaults to %f",
          this->robot_namespace_.c_str(), this->wheel_separation_);
    } else {
      this->wheel_separation_ = 
        _sdf->GetElement("wheelSeparation")->Get<double>();
    }

    this->wheel_diameter_ = 0.15;
    if (!_sdf->HasElement("wheelDiameter")) {
      ROS_WARN("GazeboRosDiffDriveMultiWheel Plugin (ns = %s) missing <wheelDiameter>, defaults to %f",
          this->robot_namespace_.c_str(), this->wheel_diameter_);
    } else {
      this->wheel_diameter_ = _sdf->GetElement("wheelDiameter")->Get<double>();
    }

    this->torque = 5.0;
    if (!_sdf->HasElement("torque")) {
      ROS_WARN("GazeboRosDiffDriveMultiWheel Plugin (ns = %s) missing <torque>, defaults to %f",
          this->robot_namespace_.c_str(), this->torque);
    } else {
      this->torque = _sdf->GetElement("torque")->Get<double>();
    }

    this->command_topic_ = "cmd_vel";
    if (!_sdf->HasElement("commandTopic")) {
      ROS_WARN("GazeboRosDiffDriveMultiWheel Plugin (ns = %s) missing <commandTopic>, defaults to \"%s\"",
          this->robot_namespace_.c_str(), this->command_topic_.c_str());
    } else {
      this->command_topic_ = _sdf->GetElement("commandTopic")->Get<std::string>();
    }

    this->odometry_topic_ = "odom";
    if (!_sdf->HasElement("odometryTopic")) {
      ROS_WARN("GazeboRosDiffDriveMultiWheel Plugin (ns = %s) missing <odometryTopic>, defaults to \"%s\"",
          this->robot_namespace_.c_str(), this->odometry_topic_.c_str());
    } else {
      this->odometry_topic_ = _sdf->GetElement("odometryTopic")->Get<std::string>();
    }

    this->odometry_frame_ = "odom";
    if (!_sdf->HasElement("odometryFrame")) {
      ROS_WARN("GazeboRosDiffDriveMultiWheel Plugin (ns = %s) missing <odometryFrame>, defaults to \"%s\"",
          this->robot_namespace_.c_str(), this->odometry_frame_.c_str());
    } else {
      this->odometry_frame_ = _sdf->GetElement("odometryFrame")->Get<std::string>();
    }

    this->robot_base_frame_ = "base_footprint";
    if (!_sdf->HasElement("robotBaseFrame")) {
      ROS_WARN("GazeboRosDiffDriveMultiWheel Plugin (ns = %s) missing <robotBaseFrame>, defaults to \"%s\"",
          this->robot_namespace_.c_str(), this->robot_base_frame_.c_str());
    } else {
      this->robot_base_frame_ = _sdf->GetElement("robotBaseFrame")->Get<std::string>();
    }

    this->update_rate_ = 100.0;
    if (!_sdf->HasElement("updateRate")) {
      ROS_WARN("GazeboRosDiffDriveMultiWheel Plugin (ns = %s) missing <updateRate>, defaults to %f",
          this->robot_namespace_.c_str(), this->update_rate_);
    } else {
      this->update_rate_ = _sdf->GetElement("updateRate")->Get<double>();
    }


    this->publish_odometry_tf_ = true;
    if (!_sdf->HasElement("publishOdometryTf")) {
      ROS_WARN("GazeboRosDiffDriveMultiWheel Plugin (ns = %s) missing <publishOdometryTf>, defaults to %s",
               this->robot_namespace_.c_str(), this->publish_odometry_tf_ ? "true" : "false");
    } else {
      this->publish_odometry_tf_ = _sdf->GetElement("publishOdometryTf")->Get<bool>();
    }

    this->publish_odometry_msg_ = true;
    if (!_sdf->HasElement("publishOdometryMsg")) {
      ROS_WARN("GazeboRosDiffDriveMultiWheel Plugin (ns = %s) missing <publishOdometryMsg>, defaults to %s",
               this->robot_namespace_.c_str(), this->publish_odometry_msg_ ? "true" : "false");
    } else {
      this->publish_odometry_msg_ = _sdf->GetElement("publishOdometryMsg")->Get<bool>();
    }



    // Initialize update rate stuff
    if (this->update_rate_ > 0.0) {
      this->update_period_ = 1.0 / this->update_rate_;
    } else {
      this->update_period_ = 0.0;
    }
#if (GAZEBO_MAJOR_VERSION >= 8)
    last_update_time_ = this->world->SimTime();
#else
    last_update_time_ = this->world->GetSimTime();
#endif

    // Initialize velocity stuff
    wheel_speed_[RIGHT] = 0;
    wheel_speed_[LEFT] = 0;

    x_ = 0;
    rot_ = 0;
    alive_ = true;

    for (size_t side = 0; side < 2; ++side){
      for (size_t i = 0; i < joint_names_[side].size(); ++i){
        joints_[side].push_back(this->parent->GetJoint(joint_names_[side][i]));
        if (!joints_[side][i]){
          char error[200];
          snprintf(error, 200,
                   "GazeboRosDiffDriveMultiWheel Plugin (ns = %s) couldn't get hinge joint named \"%s\"",
                   this->robot_namespace_.c_str(), joint_names_[side][i].c_str());
          gzthrow(error);
        }
#if (GAZEBO_MAJOR_VERSION > 4)
        joints_[side][i]->SetEffortLimit(0, torque);
#else
        joints_[side][i]->SetMaxForce(0, torque);
#endif
      }
    }

    // Make sure the ROS node for Gazebo has already been initialized
    if (!ros::isInitialized())
    {
      ROS_FATAL_STREAM("A ROS node for Gazebo has not been initialized, unable to load plugin. "
        << "Load the Gazebo system plugin 'libgazebo_ros_api_plugin.so' in the gazebo_ros package)");
      return;
    }

    rosnode_ = new ros::NodeHandle(this->robot_namespace_);

    ROS_INFO("Starting GazeboRosDiffDriveMultiWheel Plugin (ns = %s)!", this->robot_namespace_.c_str());

    tf_prefix_ = tf::getPrefixParam(*rosnode_);
    transform_broadcaster_ = new tf::TransformBroadcaster();

    // ROS: Subscribe to the velocity command topic (usually "cmd_vel")
    ros::SubscribeOptions so =
      ros::SubscribeOptions::create<geometry_msgs::Twist>(command_topic_, 1,
          boost::bind(&GazeboRosDiffDriveMultiWheel::cmdVelCallback, this, _1),
          ros::VoidPtr(), &queue_);

    cmd_vel_subscriber_ = rosnode_->subscribe(so);

    odometry_publisher_ = rosnode_->advertise<nav_msgs::Odometry>(odometry_topic_, 1);

    // start custom queue for diff drive
    this->callback_queue_thread_ = 
      boost::thread(boost::bind(&GazeboRosDiffDriveMultiWheel::QueueThread, this));

    // listen to the update event (broadcast every simulation iteration)
    this->update_connection_ = 
      event::Events::ConnectWorldUpdateBegin(
          boost::bind(&GazeboRosDiffDriveMultiWheel::UpdateChild, this));

  }

  // Update the controller
  void GazeboRosDiffDriveMultiWheel::UpdateChild() {
#if (GAZEBO_MAJOR_VERSION >= 8)
    common::Time current_time = this->world->SimTime();
#else
    common::Time current_time = this->world->GetSimTime();
#endif
    double seconds_since_last_update = 
      (current_time - last_update_time_).Double();
    if (seconds_since_last_update > update_period_) {

      if (this->publish_odometry_tf_ || this->publish_odometry_msg_){
        publishOdometry(seconds_since_last_update);
      }

      // Update robot in case new velocities have been requested
      getWheelVelocities();
      //joints[LEFT]->SetVelocity(0, wheel_speed_[LEFT] / wheel_diameter_);
      //joints[RIGHT]->SetVelocity(0, wheel_speed_[RIGHT] / wheel_diameter_);

      for (size_t side = 0; side < 2; ++side){
        for (size_t i = 0; i < joints_[side].size(); ++i){
          joints_[side][i]->SetVelocity(0, wheel_speed_[side] / (0.5 * wheel_diameter_));
        }
      }

      last_update_time_+= common::Time(update_period_);

    }
  }

  // Finalize the controller
  void GazeboRosDiffDriveMultiWheel::FiniChild() {
    alive_ = false;
    queue_.clear();
    queue_.disable();
    rosnode_->shutdown();
    callback_queue_thread_.join();
  }

  void GazeboRosDiffDriveMultiWheel::getWheelVelocities() {
    boost::mutex::scoped_lock scoped_lock(lock);

    double vr = x_;
    double va = rot_;

    wheel_speed_[LEFT] = vr - va * wheel_separation_ / 2.0;
    wheel_speed_[RIGHT] = vr + va * wheel_separation_ / 2.0;
  }

  void GazeboRosDiffDriveMultiWheel::cmdVelCallback(
      const geometry_msgs::Twist::ConstPtr& cmd_msg) {

    boost::mutex::scoped_lock scoped_lock(lock);
    x_ = cmd_msg->linear.x;
    rot_ = cmd_msg->angular.z;
  }

  void GazeboRosDiffDriveMultiWheel::QueueThread() {
    static const double timeout = 0.01;

    while (alive_ && rosnode_->ok()) {
      queue_.callAvailable(ros::WallDuration(timeout));
    }
  }

  void GazeboRosDiffDriveMultiWheel::publishOdometry(double step_time) {
    ros::Time current_time = ros::Time::now();
    std::string odom_frame = 
      tf::resolve(tf_prefix_, odometry_frame_);
    std::string base_footprint_frame = 
      tf::resolve(tf_prefix_, robot_base_frame_);

    // getting data for base_footprint to odom transform
#if (GAZEBO_MAJOR_VERSION >= 8)
    ignition::math::Pose3d pose = this->parent->WorldPose();

    tf::Quaternion qt(pose.Rot().X(), pose.Rot().Y(), pose.Rot().Z(), pose.Rot().W());
    tf::Vector3 vt(pose.Pos().X(), pose.Pos().Y(), pose.Pos().Z());
#else
    math::Pose pose = this->parent->GetWorldPose();

    tf::Quaternion qt(pose.rot.x, pose.rot.y, pose.rot.z, pose.rot.w);
    tf::Vector3 vt(pose.pos.x, pose.pos.y, pose.pos.z);
#endif

    tf::Transform base_footprint_to_odom(qt, vt);

    if (this->publish_odometry_tf_){
      transform_broadcaster_->sendTransform(
            tf::StampedTransform(base_footprint_to_odom, current_time,
                                 odom_frame, base_footprint_frame));
    }

    // publish odom topic
#if (GAZEBO_MAJOR_VERSION >= 8)
    odom_.pose.pose.position.x = pose.Pos().X();
    odom_.pose.pose.position.y = pose.Pos().Y();

    odom_.pose.pose.orientation.x = pose.Rot().X();
    odom_.pose.pose.orientation.y = pose.Rot().Y();
    odom_.pose.pose.orientation.z = pose.Rot().Z();
    odom_.pose.pose.orientation.w = pose.Rot().W();
#else
    odom_.pose.pose.position.x = pose.pos.x;
    odom_.pose.pose.position.y = pose.pos.y;

    odom_.pose.pose.orientation.x = pose.rot.x;
    odom_.pose.pose.orientation.y = pose.rot.y;
    odom_.pose.pose.orientation.z = pose.rot.z;
    odom_.pose.pose.orientation.w = pose.rot.w;
#endif
    odom_.pose.covariance[0] = 0.00001;
    odom_.pose.covariance[7] = 0.00001;
    odom_.pose.covariance[14] = 1000000000000.0;
    odom_.pose.covariance[21] = 1000000000000.0;
    odom_.pose.covariance[28] = 1000000000000.0;
    odom_.pose.covariance[35] = 0.001;

    // get velocity in /odom frame
#if (GAZEBO_MAJOR_VERSION >= 8)
    ignition::math::Vector3d linear;
    linear = this->parent->WorldLinearVel();
    odom_.twist.twist.angular.z = this->parent->WorldAngularVel().Z();
#else
    math::Vector3 linear;
    linear = this->parent->GetWorldLinearVel();
    odom_.twist.twist.angular.z = this->parent->GetWorldAngularVel().z;
#endif

    // convert velocity to child_frame_id (aka base_footprint)
#if (GAZEBO_MAJOR_VERSION >= 8)
    float yaw = pose.Rot().Yaw();
    odom_.twist.twist.linear.x = cosf(yaw) * linear.X() + sinf(yaw) * linear.Y();
    odom_.twist.twist.linear.y = cosf(yaw) * linear.Y() - sinf(yaw) * linear.X();
#else
    float yaw = pose.rot.GetYaw();
    odom_.twist.twist.linear.x = cosf(yaw) * linear.x + sinf(yaw) * linear.y;
    odom_.twist.twist.linear.y = cosf(yaw) * linear.y - sinf(yaw) * linear.x;
#endif

    odom_.header.stamp = current_time;
    odom_.header.frame_id = odom_frame;
    odom_.child_frame_id = base_footprint_frame;

    if (this->publish_odometry_msg_){
      odometry_publisher_.publish(odom_);
    }
  }

  GZ_REGISTER_MODEL_PLUGIN(GazeboRosDiffDriveMultiWheel)
}
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_gazebo/hector_gazebo_plugins/src/diffdrive_plugin_6w.cpp
//=================================================================================================
// Copyright (c) 2011, Stefan Kohlbrecher, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Simulation, Systems Optimization and Robotics
//       group, TU Darmstadt nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

/**
 * Based on diffdrive_plugin by Nathan Koenig, Andrew Howard and Daniel Hewlett
 */

#include <algorithm>
#include <assert.h>

#include <hector_gazebo_plugins/diffdrive_plugin_6w.h>
#include <gazebo/common/Events.hh>
#include <gazebo/physics/physics.hh>

#include <ros/ros.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_listener.h>
#include <geometry_msgs/Twist.h>
#include <nav_msgs/Odometry.h>
#include <boost/bind.hpp>

#include <gazebo/gazebo_config.h>

namespace gazebo {

enum
{
  FRONT_LEFT,
  FRONT_RIGHT,
  MID_LEFT,
  MID_RIGHT,
  REAR_LEFT,
  REAR_RIGHT,
  NUM_WHEELS
};

// Constructor
DiffDrivePlugin6W::DiffDrivePlugin6W()
{
}

// Destructor
DiffDrivePlugin6W::~DiffDrivePlugin6W()
{
#if (GAZEBO_MAJOR_VERSION >= 8)
  updateConnection.reset();
#else
  event::Events::DisconnectWorldUpdateBegin(updateConnection);
#endif
  delete transform_broadcaster_;
  rosnode_->shutdown();
  callback_queue_thread_.join();
  delete rosnode_;
}

// Load the controller
void DiffDrivePlugin6W::Load(physics::ModelPtr _model, sdf::ElementPtr _sdf)
{
  world = _model->GetWorld();

  // default parameters
  namespace_.clear();
  topic_ = "cmd_vel";
  wheelSep = 0.34;
  wheelDiam = 0.15;
  torque = 10.0;

  // load parameters
  if (_sdf->HasElement("robotNamespace"))
    namespace_ = _sdf->GetElement("robotNamespace")->GetValue()->GetAsString();

  if (_sdf->HasElement("topicName"))
    topic_ = _sdf->GetElement("topicName")->GetValue()->GetAsString();

  if (_sdf->HasElement("bodyName"))
  {
    link_name_ = _sdf->GetElement("bodyName")->GetValue()->GetAsString();
    link = _model->GetLink(link_name_);
  } else {
    link = _model->GetLink();
    link_name_ = link->GetName();
  }

  // assert that the body by link_name_ exists
  if (!link)
  {
    ROS_FATAL("DiffDrivePlugin6W plugin error: bodyName: %s does not exist\n", link_name_.c_str());
    return;
  }

  if (_sdf->HasElement("frontLeftJoint"))  joints[FRONT_LEFT]  = _model->GetJoint(_sdf->GetElement("frontLeftJoint")->GetValue()->GetAsString());
  if (_sdf->HasElement("frontRightJoint")) joints[FRONT_RIGHT] = _model->GetJoint(_sdf->GetElement("frontRightJoint")->GetValue()->GetAsString());
  if (_sdf->HasElement("midLeftJoint"))    joints[MID_LEFT]    = _model->GetJoint(_sdf->GetElement("midLeftJoint")->GetValue()->GetAsString());
  if (_sdf->HasElement("midRightJoint"))   joints[MID_RIGHT]   = _model->GetJoint(_sdf->GetElement("midRightJoint")->GetValue()->GetAsString());
  if (_sdf->HasElement("rearLeftJoint"))   joints[REAR_LEFT]   = _model->GetJoint(_sdf->GetElement("rearLeftJoint")->GetValue()->GetAsString());
  if (_sdf->HasElement("rearRightJoint"))  joints[REAR_RIGHT]  = _model->GetJoint(_sdf->GetElement("rearRightJoint")->GetValue()->GetAsString());

  if (!joints[FRONT_LEFT])  ROS_FATAL("diffdrive_plugin_6w error: The controller couldn't get front left joint");
  if (!joints[FRONT_RIGHT]) ROS_FATAL("diffdrive_plugin_6w error: The controller couldn't get front right joint");
  if (!joints[MID_LEFT])    ROS_FATAL("diffdrive_plugin_6w error: The controller couldn't get mid left joint");
  if (!joints[MID_RIGHT])   ROS_FATAL("diffdrive_plugin_6w error: The controller couldn't get mid right joint");
  if (!joints[REAR_LEFT])   ROS_FATAL("diffdrive_plugin_6w error: The controller couldn't get rear left joint");
  if (!joints[REAR_RIGHT])  ROS_FATAL("diffdrive_plugin_6w error: The controller couldn't get rear right joint");

  if (_sdf->HasElement("wheelSeparation"))
    _sdf->GetElement("wheelSeparation")->GetValue()->Get(wheelSep);

  if (_sdf->HasElement("wheelDiameter"))
    _sdf->GetElement("wheelDiameter")->GetValue()->Get(wheelDiam);

  if (_sdf->HasElement("torque"))
    _sdf->GetElement("torque")->GetValue()->Get(torque);

  // Make sure the ROS node for Gazebo has already been initialized
  if (!ros::isInitialized())
  {
    ROS_FATAL_STREAM("A ROS node for Gazebo has not been initialized, unable to load plugin. "
      << "Load the Gazebo system plugin 'libgazebo_ros_api_plugin.so' in the gazebo_ros package)");
    return;
  }

  rosnode_ = new ros::NodeHandle(namespace_);

  tf_prefix_ = tf::getPrefixParam(*rosnode_);
  transform_broadcaster_ = new tf::TransformBroadcaster();

  // ROS: Subscribe to the velocity command topic (usually "cmd_vel")
  ros::SubscribeOptions so =
      ros::SubscribeOptions::create<geometry_msgs::Twist>(topic_, 1,
                                                          boost::bind(&DiffDrivePlugin6W::cmdVelCallback, this, _1),
                                                          ros::VoidPtr(), &queue_);
  sub_ = rosnode_->subscribe(so);
  pub_ = rosnode_->advertise<nav_msgs::Odometry>("odom", 1);

  callback_queue_thread_ = boost::thread(boost::bind(&DiffDrivePlugin6W::QueueThread, this));

  Reset();

  // New Mechanism for Updating every World Cycle
  // Listen to the update event. This event is broadcast every
  // simulation iteration.
  updateConnection = event::Events::ConnectWorldUpdateBegin(
      boost::bind(&DiffDrivePlugin6W::Update, this));
}

// Initialize the controller
void DiffDrivePlugin6W::Reset()
{
  enableMotors = true;

  for (size_t i = 0; i < 2; ++i){
    wheelSpeed[i] = 0;
  }

#if (GAZEBO_MAJOR_VERSION >= 8)
  prevUpdateTime = world->SimTime();
#else
  prevUpdateTime = world->GetSimTime();
#endif

  x_ = 0;
  rot_ = 0;
  alive_ = true;

  // Reset odometric pose
  odomPose[0] = 0.0;
  odomPose[1] = 0.0;
  odomPose[2] = 0.0;

  odomVel[0] = 0.0;
  odomVel[1] = 0.0;
  odomVel[2] = 0.0;
}

// Update the controller
void DiffDrivePlugin6W::Update()
{
  // TODO: Step should be in a parameter of this function
  double d1, d2;
  double dr, da;
  common::Time stepTime;

  GetPositionCmd();

  //stepTime = World::Instance()->GetPhysicsEngine()->GetStepTime();
#if (GAZEBO_MAJOR_VERSION >= 8)
  stepTime = world->SimTime() - prevUpdateTime;
  prevUpdateTime = world->SimTime();
#else
  stepTime = world->GetSimTime() - prevUpdateTime;
  prevUpdateTime = world->GetSimTime();
#endif

  // Distance travelled by front wheels
  d1 = stepTime.Double() * wheelDiam / 2 * joints[MID_LEFT]->GetVelocity(0);
  d2 = stepTime.Double() * wheelDiam / 2 * joints[MID_RIGHT]->GetVelocity(0);

  dr = (d1 + d2) / 2;
  da = (d1 - d2) / wheelSep;

  // Compute odometric pose
  odomPose[0] += dr * cos(odomPose[2]);
  odomPose[1] += dr * sin(odomPose[2]);
  odomPose[2] += da;

  // Compute odometric instantaneous velocity
  odomVel[0] = dr / stepTime.Double();
  odomVel[1] = 0.0;
  odomVel[2] = da / stepTime.Double();

  if (enableMotors)
  {
    joints[FRONT_LEFT]->SetVelocity(0, wheelSpeed[0] / (wheelDiam / 2.0));
    joints[MID_LEFT]->SetVelocity(0, wheelSpeed[0] / (wheelDiam / 2.0));
    joints[REAR_LEFT]->SetVelocity(0, wheelSpeed[0] / (wheelDiam / 2.0));

    joints[FRONT_RIGHT]->SetVelocity(0, wheelSpeed[1] / (wheelDiam / 2.0));
    joints[MID_RIGHT]->SetVelocity(0, wheelSpeed[1] / (wheelDiam / 2.0));
    joints[REAR_RIGHT]->SetVelocity(0, wheelSpeed[1] / (wheelDiam / 2.0));

#if (GAZEBO_MAJOR_VERSION > 4)
    joints[FRONT_LEFT]->SetEffortLimit(0, torque);
    joints[MID_LEFT]->SetEffortLimit(0, torque);
    joints[REAR_LEFT]->SetEffortLimit(0, torque);

    joints[FRONT_RIGHT]->SetEffortLimit(0, torque);
    joints[MID_RIGHT]->SetEffortLimit(0, torque);
    joints[REAR_RIGHT]->SetEffortLimit(0, torque);
#else
    joints[FRONT_LEFT]->SetMaxForce(0, torque);
    joints[MID_LEFT]->SetMaxForce(0, torque);
    joints[REAR_LEFT]->SetMaxForce(0, torque);

    joints[FRONT_RIGHT]->SetMaxForce(0, torque);
    joints[MID_RIGHT]->SetMaxForce(0, torque);
    joints[REAR_RIGHT]->SetMaxForce(0, torque);
#endif
  }

  //publish_odometry();
}

// NEW: Now uses the target velocities from the ROS message, not the Iface 
void DiffDrivePlugin6W::GetPositionCmd()
{
  lock.lock();

  double vr, va;

  vr = x_; //myIface->data->cmdVelocity.pos.x;
  va = -rot_; //myIface->data->cmdVelocity.yaw;

  //std::cout << "X: [" << x_ << "] ROT: [" << rot_ << "]" << std::endl;

  // Changed motors to be always on, which is probably what we want anyway
  enableMotors = true; //myIface->data->cmdEnableMotors > 0;

  //std::cout << enableMotors << std::endl;

  wheelSpeed[0] = vr + va * wheelSep / 2;
  wheelSpeed[1] = vr - va * wheelSep / 2;

  lock.unlock();
}

// NEW: Store the velocities from the ROS message
void DiffDrivePlugin6W::cmdVelCallback(const geometry_msgs::Twist::ConstPtr& cmd_msg)
{
  //std::cout << "BEGIN CALLBACK\n";

  lock.lock();

  x_ = cmd_msg->linear.x;
  rot_ = cmd_msg->angular.z;

  lock.unlock();

  //std::cout << "END CALLBACK\n";
}

// NEW: custom callback queue thread
void DiffDrivePlugin6W::QueueThread()
{
  static const double timeout = 0.01;

  while (alive_ && rosnode_->ok())
  {
    //    std::cout << "CALLING STUFF\n";
    queue_.callAvailable(ros::WallDuration(timeout));
  }
}

// NEW: Update this to publish odometry topic
void DiffDrivePlugin6W::publish_odometry()
{
  // get current time
#if (GAZEBO_MAJOR_VERSION >= 8)
  ros::Time current_time_((world->SimTime()).sec, (world->SimTime()).nsec); 
#else
  ros::Time current_time_((world->GetSimTime()).sec, (world->GetSimTime()).nsec); 
#endif

  // getting data for base_footprint to odom transform
#if (GAZEBO_MAJOR_VERSION >= 8)
  ignition::math::Pose3d pose = link->WorldPose();
  ignition::math::Vector3d velocity = link->WorldLinearVel();
  ignition::math::Vector3d angular_velocity = link->WorldAngularVel();

  tf::Quaternion qt(pose.Rot().X(), pose.Rot().Y(), pose.Rot().Z(), pose.Rot().W());
  tf::Vector3 vt(pose.Pos().X(), pose.Pos().Y(), pose.Pos().Z());
#else
  math::Pose pose = link->GetWorldPose();
  math::Vector3 velocity = link->GetWorldLinearVel();
  math::Vector3 angular_velocity = link->GetWorldAngularVel();

  tf::Quaternion qt(pose.rot.x, pose.rot.y, pose.rot.z, pose.rot.w);
  tf::Vector3 vt(pose.pos.x, pose.pos.y, pose.pos.z);
#endif
  tf::Transform base_footprint_to_odom(qt, vt);

  transform_broadcaster_->sendTransform(tf::StampedTransform(base_footprint_to_odom,
                                                            current_time_,
                                                            "odom",
                                                            "base_footprint"));

  // publish odom topic
#if (GAZEBO_MAJOR_VERSION >= 8)
  odom_.pose.pose.position.x = pose.Pos().X();
  odom_.pose.pose.position.y = pose.Pos().Y();

  odom_.pose.pose.orientation.x = pose.Rot().X();
  odom_.pose.pose.orientation.y = pose.Rot().Y();
  odom_.pose.pose.orientation.z = pose.Rot().Z();
  odom_.pose.pose.orientation.w = pose.Rot().W();

  odom_.twist.twist.linear.x = velocity.X();
  odom_.twist.twist.linear.y = velocity.Y();
  odom_.twist.twist.angular.z = angular_velocity.Z();
#else
  odom_.pose.pose.position.x = pose.pos.x;
  odom_.pose.pose.position.y = pose.pos.y;

  odom_.pose.pose.orientation.x = pose.rot.x;
  odom_.pose.pose.orientation.y = pose.rot.y;
  odom_.pose.pose.orientation.z = pose.rot.z;
  odom_.pose.pose.orientation.w = pose.rot.w;

  odom_.twist.twist.linear.x = velocity.x;
  odom_.twist.twist.linear.y = velocity.y;
  odom_.twist.twist.angular.z = angular_velocity.z;
#endif

  odom_.header.frame_id = tf::resolve(tf_prefix_, "odom");
  odom_.child_frame_id = "base_footprint";
  odom_.header.stamp = current_time_;

  pub_.publish(odom_);
}

GZ_REGISTER_MODEL_PLUGIN(DiffDrivePlugin6W)

} // namespace gazebo
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_gazebo/hector_gazebo_plugins/src/gazebo_ros_magnetic.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_gazebo_plugins/gazebo_ros_magnetic.h>
#include <gazebo/common/Events.hh>
#include <gazebo/physics/physics.hh>

static const double DEFAULT_MAGNITUDE           = 1.0;
static const double DEFAULT_REFERENCE_HEADING   = 0.0;
static const double DEFAULT_DECLINATION         = 0.0;
static const double DEFAULT_INCLINATION         = 60.0;

namespace gazebo {

GazeboRosMagnetic::GazeboRosMagnetic()
{
}

////////////////////////////////////////////////////////////////////////////////
// Destructor
GazeboRosMagnetic::~GazeboRosMagnetic()
{
  updateTimer.Disconnect(updateConnection);

  dynamic_reconfigure_server_.reset();

  node_handle_->shutdown();
  delete node_handle_;
}

////////////////////////////////////////////////////////////////////////////////
// Load the controller
void GazeboRosMagnetic::Load(physics::ModelPtr _model, sdf::ElementPtr _sdf)
{
  world = _model->GetWorld();

  // load parameters
  if (_sdf->HasElement("robotNamespace"))
    namespace_ = _sdf->GetElement("robotNamespace")->GetValue()->GetAsString();
  else
    namespace_.clear();

  if (!_sdf->HasElement("topicName"))
    topic_ = "magnetic";
  else
    topic_ = _sdf->GetElement("topicName")->Get<std::string>();


  if (_sdf->HasElement("bodyName"))
  {
    link_name_ = _sdf->GetElement("bodyName")->GetValue()->GetAsString();
    link = _model->GetLink(link_name_);
  }
  else
  {
    link = _model->GetLink();
    link_name_ = link->GetName();
  }

  if (!link)
  {
    ROS_FATAL("GazeboRosMagnetic plugin error: bodyName: %s does not exist\n", link_name_.c_str());
    return;
  }

  // default parameters
  frame_id_ = link_name_;
  magnitude_ = DEFAULT_MAGNITUDE;
  reference_heading_ = DEFAULT_REFERENCE_HEADING * M_PI/180.0;
  declination_ = DEFAULT_DECLINATION * M_PI/180.0;
  inclination_ = DEFAULT_INCLINATION * M_PI/180.0;

  if (_sdf->HasElement("frameId"))
    frame_id_ = _sdf->GetElement("frameId")->GetValue()->GetAsString();

  if (_sdf->HasElement("magnitude"))
    _sdf->GetElement("magnitude")->GetValue()->Get(magnitude_);

  if (_sdf->HasElement("referenceHeading"))
    if (_sdf->GetElement("referenceHeading")->GetValue()->Get(reference_heading_))
      reference_heading_ *= M_PI/180.0;

  if (_sdf->HasElement("declination"))
    if (_sdf->GetElement("declination")->GetValue()->Get(declination_))
      declination_ *= M_PI/180.0;

  if (_sdf->HasElement("inclination"))
    if (_sdf->GetElement("inclination")->GetValue()->Get(inclination_))
      inclination_ *= M_PI/180.0;

  // Note: Gazebo uses NorthWestUp coordinate system, heading and declination are compass headings
  magnetic_field_.header.frame_id = frame_id_;
#if (GAZEBO_MAJOR_VERSION >= 8)
  magnetic_field_world_.X() = magnitude_ *  cos(inclination_) * cos(reference_heading_ - declination_);
  magnetic_field_world_.Y() = magnitude_ *  cos(inclination_) * sin(reference_heading_ - declination_);
  magnetic_field_world_.Z() = magnitude_ * -sin(inclination_);
#else
  magnetic_field_world_.x = magnitude_ *  cos(inclination_) * cos(reference_heading_ - declination_);
  magnetic_field_world_.y = magnitude_ *  cos(inclination_) * sin(reference_heading_ - declination_);
  magnetic_field_world_.z = magnitude_ * -sin(inclination_);
#endif

  sensor_model_.Load(_sdf);

  // start ros node
  if (!ros::isInitialized())
  {
    int argc = 0;
    char** argv = NULL;
    ros::init(argc,argv,"gazebo",ros::init_options::NoSigintHandler|ros::init_options::AnonymousName);
  }

  node_handle_ = new ros::NodeHandle(namespace_);
  publisher_ = node_handle_->advertise<geometry_msgs::Vector3Stamped>(topic_, 1);

  // setup dynamic_reconfigure server
  dynamic_reconfigure_server_.reset(new dynamic_reconfigure::Server<SensorModelConfig>(ros::NodeHandle(*node_handle_, topic_)));
  dynamic_reconfigure_server_->setCallback(boost::bind(&SensorModel3::dynamicReconfigureCallback, &sensor_model_, _1, _2));

  Reset();

  // connect Update function
  updateTimer.Load(world, _sdf);
  updateConnection = updateTimer.Connect(boost::bind(&GazeboRosMagnetic::Update, this));
}

void GazeboRosMagnetic::Reset()
{
  updateTimer.Reset();
  sensor_model_.reset();
}

////////////////////////////////////////////////////////////////////////////////
// Update the controller
void GazeboRosMagnetic::Update()
{
#if (GAZEBO_MAJOR_VERSION >= 8)
  common::Time sim_time = world->SimTime();
  double dt = updateTimer.getTimeSinceLastUpdate().Double();

  ignition::math::Pose3d pose = link->WorldPose();
  ignition::math::Vector3d field = sensor_model_(pose.Rot().RotateVectorReverse(magnetic_field_world_), dt);

  magnetic_field_.header.stamp = ros::Time(sim_time.sec, sim_time.nsec);
  magnetic_field_.vector.x = field.X();
  magnetic_field_.vector.y = field.Y();
  magnetic_field_.vector.z = field.Z();
#else
  common::Time sim_time = world->GetSimTime();
  double dt = updateTimer.getTimeSinceLastUpdate().Double();

  math::Pose pose = link->GetWorldPose();
  math::Vector3 field = sensor_model_(pose.rot.RotateVectorReverse(magnetic_field_world_), dt);

  magnetic_field_.header.stamp = ros::Time(sim_time.sec, sim_time.nsec);
  magnetic_field_.vector.x = field.x;
  magnetic_field_.vector.y = field.y;
  magnetic_field_.vector.z = field.z;
#endif

  publisher_.publish(magnetic_field_);
}

// Register this plugin with the simulator
GZ_REGISTER_MODEL_PLUGIN(GazeboRosMagnetic)

} // namespace gazebo
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_gazebo/hector_gazebo_plugins/src/reset_plugin.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_gazebo_plugins/reset_plugin.h>
#include <gazebo/common/Events.hh>

#include <std_msgs/String.h>

namespace gazebo {

GazeboResetPlugin::GazeboResetPlugin()
{
}

////////////////////////////////////////////////////////////////////////////////
// Destructor
GazeboResetPlugin::~GazeboResetPlugin()
{
  node_handle_->shutdown();
  delete node_handle_;
}

////////////////////////////////////////////////////////////////////////////////
// Load the controller
void GazeboResetPlugin::Load(physics::ModelPtr _model, sdf::ElementPtr _sdf)
{
  // Make sure the ROS node for Gazebo has already been initialized
  if (!ros::isInitialized())
  {
    ROS_FATAL_STREAM("A ROS node for Gazebo has not been initialized, unable to load plugin. "
      << "Load the Gazebo system plugin 'libgazebo_ros_api_plugin.so' in the gazebo_ros package)");
    return;
  }

  node_handle_ = new ros::NodeHandle;
  publisher_ = node_handle_->advertise<std_msgs::String>("/syscommand", 1);
}

////////////////////////////////////////////////////////////////////////////////
// Reset the controller
void GazeboResetPlugin::Reset()
{
  std_msgs::String command;
  command.data = "reset";
  publisher_.publish(command);
}

// Register this plugin with the simulator
GZ_REGISTER_MODEL_PLUGIN(GazeboResetPlugin)

} // namespace gazebo
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/global_reference.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/global_reference.h>
#include <hector_pose_estimation/state.h>
#include <cmath>

#include <ros/console.h>

using namespace std;

namespace hector_pose_estimation {

GlobalReference::GlobalReference()
{
  parameters().add("reference_latitude",  reference_latitude_  = std::numeric_limits<double>::quiet_NaN());
  parameters().add("reference_longitude", reference_longitude_ = std::numeric_limits<double>::quiet_NaN());
  parameters().add("reference_altitude",  reference_altitude_  = std::numeric_limits<double>::quiet_NaN());
  parameters().add("reference_heading",   heading_.value       = std::numeric_limits<double>::quiet_NaN());

  reset();
}

const GlobalReferencePtr &GlobalReference::Instance()
{
  static GlobalReferencePtr instance;
  if (!instance) { instance.reset(new GlobalReference); }
  return instance;
}

void GlobalReference::reset()
{
  position_ = Position();
  heading_ = Heading();
  radius_ = Radius();

  // parameters are in degrees
  position_.latitude  = reference_latitude_ * M_PI/180.0;
  position_.longitude = reference_longitude_ * M_PI/180.0;
  position_.altitude  = reference_altitude_;
  heading_.value      = reference_heading_;

  updated();
}

ParameterList& GlobalReference::parameters() {
  return parameters_;
}

void GlobalReference::updated(bool intermediate) {
  // calculate earth radii
  if (hasPosition()) {
    radius_ = Radius(position_.latitude);
  }

  // calculate sin and cos of the heading reference
  if (hasHeading()) {
    sincos(heading_.value, &heading_.sin, &heading_.cos);
  }

  // execute update callbacks
  if (!intermediate) {
    for(std::list<UpdateCallback>::iterator cb = update_callbacks_.begin(); cb != update_callbacks_.end(); ++cb)
      (*cb)();
  }
}

GlobalReference::Heading::Heading(double heading) : value(heading)
{
  sincos(heading, &sin, &cos);
}

Quaternion GlobalReference::Heading::quaternion() const
{
  double sin_2, cos_2;
  sincos(0.5 * value, &sin_2, &cos_2);
  return Quaternion(cos_2, 0.0, 0.0, -sin_2);
}

GlobalReference::Radius::Radius(double latitude) {
  // WGS84 constants
  static const double equatorial_radius = 6378137.0;
  static const double flattening = 1.0/298.257223563;
  static const double excentrity2 = 2*flattening - flattening*flattening;

  double temp = 1.0 / (1.0 - excentrity2 * sin(latitude) * sin(latitude));
  double prime_vertical_radius = equatorial_radius * sqrt(temp);
  north = prime_vertical_radius * (1 - excentrity2) * temp;
  east  = prime_vertical_radius * cos(latitude);
}

void GlobalReference::fromWGS84(double latitude, double longitude, double &x, double &y) {
  if (!hasPosition()) {
    x = 0.0;
    y = 0.0;
    return;
  }
  double north = radius_.north * (latitude  - position_.latitude);
  double east  = radius_.east  * (longitude - position_.longitude);
  fromNorthEast(north, east, x, y);
}

void GlobalReference::toWGS84(double x, double y, double &latitude, double &longitude) {
  if (!hasPosition()) {
    latitude  = 0.0;
    longitude = 0.0;
    return;
  }

  double north, east;
  toNorthEast(x, y, north, east);
  latitude  = position_.latitude  + north / radius_.north;
  longitude = position_.longitude + east  / radius_.east;
}

void GlobalReference::fromNorthEast(double north, double east, double &x, double &y) {
  if (!hasHeading()) {
    x = 0.0;
    y = 0.0;
    return;
  }
  x = north * heading_.cos + east * heading_.sin;
  y = north * heading_.sin - east * heading_.cos;
}

void GlobalReference::toNorthEast(double x, double y, double &north, double &east) {
  if (!hasHeading()) {
    north = 0.0;
    east = 0.0;
    return;
  }
  north = x * heading_.cos + y * heading_.sin;
  east  = x * heading_.sin - y * heading_.cos;
}

GlobalReference& GlobalReference::setPosition(double latitude, double longitude, bool intermediate /* = false */) {
  position_.latitude = latitude;
  position_.longitude = longitude;
  if (!intermediate) ROS_INFO("Set new reference position to %f deg N / %f deg E", this->position().latitude * 180.0/M_PI, this->position().longitude * 180.0/M_PI);
  updated(intermediate);
  return *this;
}

GlobalReference& GlobalReference::setHeading(double heading, bool intermediate /* = false */) {
  heading_.value = heading;
  if (!intermediate) ROS_INFO("Set new reference heading to %.1f degress", this->heading() * 180.0 / M_PI);
  updated(intermediate);
  return *this;
}

GlobalReference& GlobalReference::setAltitude(double altitude, bool intermediate /* = false */) {
  position_.altitude = altitude;
  if (!intermediate) ROS_INFO("Set new reference altitude to %.2f m", this->position().altitude);
  updated(intermediate);
  return *this;
}

GlobalReference& GlobalReference::setCurrentPosition(const State& state, double new_latitude, double new_longitude) {
  State::ConstPositionType position = state.getPosition();

  // set reference to new latitude/longitude first (intermediate reference)
  setPosition(new_latitude, new_longitude, true);

  // convert current position back to WGS84 using the new reference position
  // and reset the reference position so that current position in x/y coordinates remains the same
  // This will work unless the radii at the origin and the x/y position of the robot differ too much
  toWGS84(-position.x(), -position.y(), new_latitude, new_longitude);
  setPosition(new_latitude, new_longitude);

  return *this;
}

GlobalReference& GlobalReference::setCurrentHeading(const State& state, double new_heading) {
  // get current yaw angle
  double current_yaw = state.getYaw();
  State::ConstPositionType position = state.getPosition();

  // get current position in WGS84
  double current_latitude, current_longitude;
  if (hasPosition()) {
    toWGS84(position.x(), position.y(), current_latitude, current_longitude);
  }

  // set the new reference heading
  setHeading(new_heading - (-current_yaw));

  // set the new reference position so that current position in WGS84 coordinates remains the same as before
  if (hasPosition()) {
    setCurrentPosition(state, current_latitude, current_longitude);
  }

  return *this;
}

GlobalReference& GlobalReference::setCurrentAltitude(const State& state, double new_altitude) {
  State::ConstPositionType position = state.getPosition();
  setAltitude(new_altitude - position.z());
  return *this;
}

void GlobalReference::getGeoPose(geographic_msgs::GeoPose& geopose) const
{
  Quaternion orientation(heading().quaternion());
  geopose.orientation.w = orientation.w();
  geopose.orientation.x = orientation.x();
  geopose.orientation.y = orientation.y();
  geopose.orientation.z = orientation.z();
  geopose.position.latitude  = position().latitude  * 180.0/M_PI;
  geopose.position.longitude = position().longitude * 180.0/M_PI;
  geopose.position.altitude  = position().altitude;
}

bool GlobalReference::getWorldToNavTransform(geometry_msgs::TransformStamped& transform, const std::string &world_frame, const std::string &nav_frame, const ros::Time& stamp) const
{
  // Return transform from world_frame (defined by parameters reference_latitude_, reference_longitude_, reference_altitude_ and reference_heading_)
  // to the nav_frame (this reference)
  if (isnan(reference_latitude_) ||
      isnan(reference_longitude_) ||
      isnan(reference_altitude_) ||
      isnan(reference_heading_)) {
    return false;
  }

  transform.header.stamp = stamp;
  transform.header.frame_id = world_frame;
  transform.child_frame_id = nav_frame;

  Radius reference_radius(reference_latitude_ * M_PI/180.0);
  double north = reference_radius.north * (position_.latitude  - reference_latitude_ * M_PI/180.0);
  double east  = reference_radius.east  * (position_.longitude - reference_longitude_ * M_PI/180.0);
  Heading reference_heading(reference_heading_ * M_PI/180.0);
  transform.transform.translation.x = north * reference_heading.cos + east * reference_heading.sin;
  transform.transform.translation.y = north * reference_heading.sin - east * reference_heading.cos;
  transform.transform.translation.z = position_.altitude - reference_altitude_;
  double heading_diff = heading().value - reference_heading;
  transform.transform.rotation.w =  cos(heading_diff / 2.);
  transform.transform.rotation.x =  0.0;
  transform.transform.rotation.y =  0.0;
  transform.transform.rotation.z = -sin(heading_diff / 2.);

  return true;
}

void GlobalReference::addUpdateCallback(const UpdateCallback &cb)
{
  update_callbacks_.push_back(cb);
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/parameters.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/parameters.h>
#include <hector_pose_estimation/ros/parameters.h>

#include <boost/algorithm/string.hpp>

#include <hector_pose_estimation/matrix.h>

#include <iostream>

namespace hector_pose_estimation {

namespace {
  struct null_deleter {
    void operator()(void const *) const {}
  };
}

ParameterList& ParameterList::add(ParameterPtr const& parameter) {
  erase(parameter->key);
  push_back(parameter);
  return *this;
}

ParameterList& ParameterList::add(ParameterList const& other) {
  for(ParameterList::const_iterator it = other.begin(); it != other.end(); ++it) push_back(*it);
  return *this;
}

ParameterList& ParameterList::add(Parameter& alias, const std::string& key) {
  if (!key.empty()) alias.key = key;
  return add(ParameterPtr(&alias, null_deleter()));
}

ParameterList& ParameterList::copy(const std::string& prefix, ParameterList const& parameters) {
  for(ParameterList::const_iterator it = parameters.begin(); it != parameters.end(); ++it) {
    ParameterPtr copy((*it)->clone());
    if (!copy) continue;
    if (!prefix.empty()) copy->key = prefix + copy->key;
    push_back(copy);
  }
  return *this;
}

ParameterList& ParameterList::copy(ParameterList const& parameters) {
  copy(std::string(), parameters);
  return *this;
}

ParameterPtr const& ParameterList::get(const std::string& key) const {
  for(const_iterator it = begin(); it != end(); ++it) {
    if ((*it)->key == key) {
      return *it;
    }
  }
  throw std::runtime_error("parameter not found");
}

ParameterList::iterator ParameterList::erase(const std::string& key) {
  iterator it = begin();
  for(; it != end(); ++it) {
    if ((*it)->key == key) return erase(it);
  }
  return it;
}

template <typename T>
struct ParameterRegistryROS::Handler
{
  bool operator()(const ParameterPtr& parameter, ros::NodeHandle& nh, bool set_all = false) {
    try {
      ParameterT<T> p(*parameter);
      std::string param_key(boost::algorithm::to_lower_copy(parameter->key));
      if (!nh.getParam(param_key, p.value())) {
        if (set_all) {
          nh.setParam(param_key, p.value());
          ROS_DEBUG_STREAM("Registered parameter " << param_key << " with new value " << p.value());
        }
      } else {
        ROS_DEBUG_STREAM("Found parameter " << param_key << " with value " << p.value());
      }
      return true;
    } catch(std::bad_cast&) {
      return false;
    }
  }
};

template <>
struct ParameterRegistryROS::Handler<ColumnVector>
{
  bool operator()(const ParameterPtr& parameter, ros::NodeHandle& nh, bool set_all = false) {
    try {
      ParameterT<ColumnVector> p(*parameter);
      std::string param_key(boost::algorithm::to_lower_copy(parameter->key));
      XmlRpc::XmlRpcValue vector;
      if (!nh.getParam(param_key, vector)) {
        if (set_all) {
          /// nh.setParam(param_key, p.value);
          ROS_DEBUG_STREAM("Not registered vector parameter " << param_key << ". Using defaults.");
        }
      } else {
        if (vector.getType() != XmlRpc::XmlRpcValue::TypeArray) {
          ROS_WARN_STREAM("Found parameter " << param_key << ", but it's not an array!");
          return false;
        }
        p.value().resize(vector.size());
        for(int i = 0; i < vector.size(); ++i) p.value()[i] = vector[i];
        ROS_DEBUG_STREAM("Found parameter " << param_key << " with value " << p.value());
      }
      return true;
    } catch(std::bad_cast&) {
      return false;
    }
  }
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vector) {
  os << "[";
  for(typename std::vector<T>::const_iterator it = vector.begin(); it != vector.end(); ++it) {
    if (it != vector.begin()) os << ", ";
    os << *it;
  }
  os << "]";
  return os;
}

template <typename T>
struct ParameterRegistryROS::Handler< std::vector<T> >
{
  bool operator()(const ParameterPtr& parameter, ros::NodeHandle& nh, bool set_all = false) {
    try {
      ParameterT< std::vector<T> > p(*parameter);
      std::string param_key(boost::algorithm::to_lower_copy(parameter->key));
      XmlRpc::XmlRpcValue vector;
      if (!nh.getParam(param_key, vector)) {
        if (set_all) {
          /// nh.setParam(param_key, p.value);
          ROS_DEBUG_STREAM("Not registered vector parameter " << param_key << ". Using defaults.");
        }
      } else {
        if (vector.getType() != XmlRpc::XmlRpcValue::TypeArray) {
          ROS_WARN_STREAM("Found parameter " << param_key << ", but it's not an array!");
          return false;
        }
        p.value().resize(vector.size());
        for(int i = 0; i < vector.size(); ++i) p.value()[i] = vector[i];
        ROS_DEBUG_STREAM("Found parameter " << param_key << " with value " << p.value());
      }
      return true;
    } catch(std::bad_cast&) {
      return false;
    }
  }
};

ParameterRegistryROS::ParameterRegistryROS(ros::NodeHandle nh)
  : nh_(nh)
  , set_all_(false)
{
  nh_.getParam("set_all_parameters", set_all_);
}

void ParameterRegistryROS::operator ()(ParameterPtr parameter) {
  // call initialize recursively for ParameterList parameters
  if (parameter->hasType<ParameterList>()) {
    ParameterList with_prefix;
    with_prefix.copy(parameter->key + "/", parameter->as<ParameterList>());
    with_prefix.initialize(*this);
    return;
  }

  ROS_DEBUG_STREAM("Registering ROS parameter " << parameter->key);

  if (Handler<std::string>()(parameter, nh_, set_all_) ||
      Handler<double>()(parameter, nh_, set_all_) ||
      Handler<std::vector<double> >()(parameter, nh_, set_all_) ||
      Handler<int>()(parameter, nh_, set_all_) ||
      Handler<bool>()(parameter, nh_, set_all_) ||
      Handler<ColumnVector>()(parameter, nh_, set_all_)
     ) {
    return;
  }

  ROS_ERROR("Parameter %s has unknown type %s!", parameter->key.c_str(), parameter->type());
}

void ParameterList::initialize(ParameterRegisterFunc func) const {
  for(const_iterator it = begin(); it != end(); ++it) {
    const ParameterPtr& parameter = *it;
    if (parameter->empty()) continue;
    if (parameter->isAlias()) continue;

    func(*it);
  }
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/measurements/baro.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/measurements/baro.h>
#include <hector_pose_estimation/filter/set_filter.h>

#include <boost/bind.hpp>

namespace hector_pose_estimation {

template class Measurement_<BaroModel>;

BaroModel::BaroModel()
{
  stddev_ = 1.0;
  qnh_ = 1013.25;
  parameters().add("qnh", qnh_);
}

BaroModel::~BaroModel() {}

void BaroModel::getExpectedValue(MeasurementVector& y_pred, const State& state)
{
  y_pred(0) = qnh_ * pow(1.0 - (0.0065 * (state.getPosition().z() + getElevation())) / 288.15, 5.255);
}

void BaroModel::getStateJacobian(MeasurementMatrix& C, const State& state, bool)
{
  if (state.position()) {
    state.position()->cols(C)(0,Z) = qnh_ * 5.255 * pow(1.0 - (0.0065 * (state.getPosition().z() + getElevation())) / 288.15, 4.255) * (-0.0065 / 288.15);
  }
}

double BaroModel::getAltitude(const BaroUpdate& update)
{
  return 288.15 / 0.0065 * (1.0 - pow(update.getVector()(0) / qnh_, 1.0/5.255));
}

BaroUpdate::BaroUpdate() : qnh_(0) {}
BaroUpdate::BaroUpdate(double pressure) : qnh_(0) { *this = pressure; }
BaroUpdate::BaroUpdate(double pressure, double qnh) : qnh_(qnh) { *this = pressure; }

Baro::Baro(const std::string &name)
  : Measurement_<BaroModel>(name)
  , HeightBaroCommon(this)
{
  parameters().add("auto_elevation", auto_elevation_);
}

void Baro::onReset()
{
  HeightBaroCommon::onReset();
}

bool Baro::prepareUpdate(State &state, const Update &update) {
  if (update.qnh() != 0) setQnh(update.qnh());
  // Note: boost::bind is not real-time safe!
  setElevation(resetElevation(state, boost::bind(&BaroModel::getAltitude, getModel(), update)));
  return true;
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/measurements/poseupdate.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/measurements/poseupdate.h>
#include <hector_pose_estimation/pose_estimation.h>

#include <Eigen/Core>

#include <boost/bind.hpp>

namespace hector_pose_estimation {

PoseUpdate::PoseUpdate(const std::string& name)
  : Measurement(name)
{
  fixed_alpha_ = 0.0;
  fixed_beta_  = 0.0;
  interpret_covariance_as_information_matrix_ = true;

  max_time_difference_ = 1.0;
  predict_pose_ = true;

  jump_on_max_error_ = true;

  fixed_position_xy_stddev_ = 0.0;
  fixed_position_z_stddev_ = 0.0;
  fixed_yaw_stddev_ = 0.0;

  fixed_velocity_xy_stddev_ = 0.0;
  fixed_velocity_z_stddev_ = 0.0;
  fixed_angular_rate_xy_stddev_ = 0.0;
  fixed_angular_rate_z_stddev_ = 0.0;

  max_position_xy_error_ = 3.0; // 3 sigma
  max_position_z_error_ = 3.0; // 3 sigma
  max_yaw_error_ = 3.0; // 3 sigma

  max_velocity_xy_error_ = 3.0; // 3 sigma
  max_velocity_z_error_ = 3.0; // 3 sigma
  max_angular_rate_xy_error_ = 3.0; // 3 sigma
  max_angular_rate_z_error_ = 3.0; // 3 sigma

  parameters().add("fixed_alpha", fixed_alpha_);
  parameters().add("fixed_beta", fixed_beta_);
  parameters().add("interpret_covariance_as_information_matrix", interpret_covariance_as_information_matrix_);
  parameters().add("max_time_difference", max_time_difference_);
  parameters().add("predict_pose", predict_pose_);
  parameters().add("jump_on_max_error", jump_on_max_error_);

  parameters().add("fixed_position_xy_stddev", fixed_position_xy_stddev_);
  parameters().add("fixed_position_z_stddev", fixed_position_z_stddev_);
  parameters().add("fixed_yaw_stddev", fixed_yaw_stddev_);
  parameters().add("fixed_velocity_xy_stddev", fixed_velocity_xy_stddev_);
  parameters().add("fixed_velocity_z_stddev", fixed_velocity_z_stddev_);
  parameters().add("fixed_angular_rate_xy_stddev", fixed_angular_rate_xy_stddev_);
  parameters().add("fixed_angular_rate_z_stddev", fixed_angular_rate_z_stddev_);
  parameters().add("max_position_xy_error", max_position_xy_error_ );
  parameters().add("max_position_z_error", max_position_z_error_);
  parameters().add("max_yaw_error", max_yaw_error_);
  parameters().add("max_velocity_xy_error", max_velocity_xy_error_ );
  parameters().add("max_velocity_z_error", max_velocity_z_error_);
  parameters().add("max_angular_rate_xy_error", max_angular_rate_xy_error_ );
  parameters().add("max_angular_rate_z_error", max_angular_rate_z_error_);
}

PoseUpdate::~PoseUpdate()
{
}

bool PoseUpdate::updateImpl(const MeasurementUpdate &update_)
{
  Update const &update = static_cast<Update const &>(update_);

  while (update.pose) {
    // convert incoming update information to Eigen
    Eigen::Vector3d update_pose(update.pose->pose.pose.position.x, update.pose->pose.pose.position.y, update.pose->pose.pose.position.z);
    Eigen::Quaterniond update_orientation(update.pose->pose.pose.orientation.w, update.pose->pose.pose.orientation.x, update.pose->pose.pose.orientation.y, update.pose->pose.pose.orientation.z);
    Eigen::Vector3d update_euler;

    // information is the information matrix if interpret_covariance_as_information_matrix_ is true and a covariance matrix otherwise
    // zero elements are counted as zero information in any case
    SymmetricMatrix6 information(Eigen::Map<const SymmetricMatrix6>(update.pose->pose.covariance.data()));

    ROS_DEBUG_STREAM_NAMED("poseupdate", "PoseUpdate: x = [ " << filter()->state().getVector().transpose() << " ], P = [ " << filter()->state().getCovariance() << " ]" << std::endl
                                      << "update: pose = [ " << update_pose.transpose() << " ], rpy = [ " << update_euler.transpose() << " ], information = [ " << information << " ]");
    ROS_DEBUG_STREAM_NAMED("poseupdate", "dt = " << (filter()->state().getTimestamp() - update.pose->header.stamp).toSec() << " s");

    // predict update pose using the estimated velocity and degrade information
    if (!update.pose->header.stamp.isZero()) {
      double dt = (filter()->state().getTimestamp() - update.pose->header.stamp).toSec();
      if (dt < 0.0) {
        ROS_DEBUG_STREAM_NAMED("poseupdate", "Ignoring pose update as it has a negative time difference: dt = " << dt << "s");
        break;

      } else if (max_time_difference_ > 0.0 && dt >= max_time_difference_) {
        ROS_DEBUG_STREAM_NAMED("poseupdate", "Ignoring pose update as the time difference is too large: dt = " << dt << "s");
        break;

      } else if (max_time_difference_ > 0.0){
        if (interpret_covariance_as_information_matrix_)
          information = information * (1.0 - dt/max_time_difference_);
        else
          information = information / (1.0 - dt/max_time_difference_);
      }

      if (predict_pose_) {
        State::ConstVelocityType state_velocity(filter()->state().getVelocity());
        update_pose = update_pose + state_velocity * dt;

        State::ConstRateType state_rate(filter()->state().getRate());
        update_orientation = update_orientation * Eigen::Quaterniond().fromRotationVector(state_rate * dt);
      }
    }

    // Calculate euler angles
    {
        const Eigen::Quaterniond &q = update_orientation;
        /* roll  = */ update_euler(0) = atan2(2*(q.y()*q.z() + q.w()*q.x()), q.w()*q.w() - q.x()*q.x() - q.y()*q.y() + q.z()*q.z());
        /* pitch = */ update_euler(1) = -asin(2*(q.x()*q.z() - q.w()*q.y()));
        /* yaw   = */ update_euler(2) = atan2(2*(q.x()*q.y() + q.w()*q.z()), q.w()*q.w() + q.x()*q.x() - q.y()*q.y() - q.z()*q.z());
    }

    // update PositionXY
    if (information(0,0) > 0.0 || information(1,1) > 0.0) {
      // fetch observation matrix H and current state x
      PositionXYModel::MeasurementMatrix H(position_xy_model_.getDimension(), filter()->state().getCovarianceDimension());
      PositionXYModel::MeasurementVector x(position_xy_model_.getDimension());
      position_xy_model_.getStateJacobian(H, filter()->state(), true);
      position_xy_model_.getExpectedValue(x, filter()->state());

      PositionXYModel::MeasurementVector y(update_pose.segment<2>(0));
      PositionXYModel::NoiseVariance Iy(information.block<2,2>(0,0));

      // invert Iy if information is a covariance matrix
      if (!interpret_covariance_as_information_matrix_) Iy = Iy.inverse().eval();

      // fixed_position_xy_stddev_ = 1.0;
      if (fixed_position_xy_stddev_ != 0.0) {
        Iy.setZero();
        Iy(0,0) = Iy(1,1) = 1.0 / (fixed_position_xy_stddev_*fixed_position_xy_stddev_);
      }

      ROS_DEBUG_STREAM_NAMED("poseupdate", "Position Update: ");
      ROS_DEBUG_STREAM_NAMED("poseupdate", "      x = [" << x.transpose() << "], H = [ " << H << " ], Px = [" <<  (H * filter()->state().P() * H.transpose()) << "], Ix = [ " << (H * filter()->state().P() * H.transpose()).inverse() << "]");
      ROS_DEBUG_STREAM_NAMED("poseupdate", "      y = [" << y.transpose() << "], Iy = [ " << Iy << " ]");
      double innovation = updateInternal(filter()->state(), Iy, y - x, H, "position_xy", max_position_xy_error_, boost::bind(&PositionXYModel::updateState, &position_xy_model_, _1, _2));
      position_xy_model_.getExpectedValue(x, filter()->state());
      ROS_DEBUG_STREAM_NAMED("poseupdate", " ==> xy = [" << x << "], Pxy = [ " << (H * filter()->state().P() * H.transpose()) << " ], innovation = " << innovation);

      status_flags_ |= STATE_POSITION_XY;
    }

    // update PositionZ
    if (information(2,2) > 0.0) {
      // fetch observation matrix H and current state x
      PositionZModel::MeasurementMatrix H(position_z_model_.getDimension(), filter()->state().getCovarianceDimension());
      PositionZModel::MeasurementVector x(position_z_model_.getDimension());
      position_z_model_.getStateJacobian(H, filter()->state(), true);
      position_z_model_.getExpectedValue(x, filter()->state());

      PositionZModel::MeasurementVector y(update_pose.segment<1>(2));
      PositionZModel::NoiseVariance Iy(information.block<1,1>(2,2));

      // invert Iy if information is a covariance matrix
      if (!interpret_covariance_as_information_matrix_) Iy = Iy.inverse().eval();

      // fixed_position_z_stddev_ = 1.0;
      if (fixed_position_z_stddev_ != 0.0) {
        Iy.setZero();
        Iy(0,0) = 1.0 / (fixed_position_z_stddev_*fixed_position_z_stddev_);
      }

      ROS_DEBUG_STREAM_NAMED("poseupdate", "Height Update: ");
      ROS_DEBUG_STREAM_NAMED("poseupdate", "      x = " << x(0) << ", H = [ " << H << " ], Px = [" <<  (H * filter()->state().P() * H.transpose()) << "], Ix = [ " << (H * filter()->state().P() * H.transpose()).inverse() << "]");
      ROS_DEBUG_STREAM_NAMED("poseupdate", "      y = " << y(0) << ", Iy = [ " << Iy << " ]");
      double innovation = updateInternal(filter()->state(), Iy, y - x, H, "position_z", max_position_z_error_, boost::bind(&PositionZModel::updateState, &position_z_model_, _1, _2));
      position_z_model_.getExpectedValue(x, filter()->state());
      ROS_DEBUG_STREAM_NAMED("poseupdate", " ==> xy = " << x(0) << ", Pxy = [ " << (H * filter()->state().P() * H.transpose()) << " ], innovation = " << innovation);

      status_flags_ |= STATE_POSITION_Z;
    }

    // update Yaw
    if (information(5,5) > 0.0) {
      YawModel::MeasurementMatrix H(yaw_model_.getDimension(), filter()->state().getCovarianceDimension());
      YawModel::MeasurementVector x(yaw_model_.getDimension());
      yaw_model_.getStateJacobian(H, filter()->state(), true);
      yaw_model_.getExpectedValue(x, filter()->state());

      YawModel::MeasurementVector y(update_euler.tail<1>());
      YawModel::NoiseVariance Iy(information.block<1,1>(5,5));

      // invert Iy if information is a covariance matrix
      if (!interpret_covariance_as_information_matrix_) Iy = Iy.inverse().eval();

      // fixed_yaw_stddev_ = 5.0 * M_PI/180.0;
      if (fixed_yaw_stddev_ != 0.0) {
        Iy.setZero();
        Iy(0,0) = 1.0 / (fixed_yaw_stddev_*fixed_yaw_stddev_);
      }

      ROS_DEBUG_STREAM_NAMED("poseupdate", "Yaw Update: ");
      ROS_DEBUG_STREAM_NAMED("poseupdate", "      x = " << x(0) * 180.0/M_PI << "°, H = [ " << H << " ], Px = [" <<  (H * filter()->state().P() * H.transpose()) << "], Ix = [ " << (H * filter()->state().P() * H.transpose()).inverse() << "]");
      ROS_DEBUG_STREAM_NAMED("poseupdate", "      y = " << y(0) * 180.0/M_PI << "°, Iy = [ " << Iy << " ]");

      YawModel::MeasurementVector error(y - x);
      error(0) = error(0) - 2.0*M_PI * round(error(0) / (2.0*M_PI));

      double innovation = updateInternal(filter()->state(), Iy, error, H, "yaw", max_yaw_error_, boost::bind(&YawModel::updateState, &yaw_model_, _1, _2));
      yaw_model_.getExpectedValue(x, filter()->state());
      ROS_DEBUG_STREAM_NAMED("poseupdate", " ==> xy = " << x(0) * 180.0/M_PI << "°, Pxy = [ " << (H * filter()->state().P() * H.transpose()) << " ], innovation = " << innovation);

      status_flags_ |= STATE_YAW;
    }

    break;
  }

  while (update.twist) {
    // convert incoming update information to Eigen
    Eigen::Vector3d update_linear(update.twist->twist.twist.linear.x, update.twist->twist.twist.linear.y, update.twist->twist.twist.linear.z);
    Eigen::Vector3d update_angular(update.twist->twist.twist.angular.x, update.twist->twist.twist.angular.y, update.twist->twist.twist.angular.z);

    // information is the information matrix if interpret_covariance_as_information_matrix_ is true and a covariance matrix otherwise
    // zero elements are counted as zero information in any case
    SymmetricMatrix6 information(Eigen::Map<const SymmetricMatrix6>(update.twist->twist.covariance.data()));

    ROS_DEBUG_STREAM_NAMED("poseupdate", "TwistUpdate:  state = [ " << filter()->state().getVector().transpose() << " ], P = [ " << filter()->state().getCovariance() << " ]" << std::endl
                                      << "     update: linear = [ " << update_linear.transpose() << " ], angular = [ " << update_angular.transpose() << " ], information = [ " << information << " ]");
    ROS_DEBUG_STREAM_NAMED("poseupdate", "                 dt = " << (filter()->state().getTimestamp() - update.twist->header.stamp).toSec() << " s");

    // degrade information if the time difference is too large
    if (!update.twist->header.stamp.isZero()) {
      double dt = (filter()->state().getTimestamp() - update.twist->header.stamp).toSec();
      if (dt < 0.0) {
        ROS_DEBUG_STREAM_NAMED("poseupdate", "Ignoring twist update as it has a negative time difference: dt = " << dt << "s");
        break;

      } else if (max_time_difference_ > 0.0 && dt >= max_time_difference_) {
        ROS_DEBUG_STREAM_NAMED("poseupdate", "Ignoring twist update as the time difference is too large: dt = " << dt << "s");
        break;

      } else if (max_time_difference_ > 0.0){
        if (interpret_covariance_as_information_matrix_)
          information = information * (1.0 - dt/max_time_difference_);
        else
          information = information / (1.0 - dt/max_time_difference_);
      }
    }

    // fetch observation matrix H and current state x
    TwistModel::MeasurementMatrix H(twist_model_.getDimension(), filter()->state().getCovarianceDimension());
    TwistModel::MeasurementVector x(twist_model_.getDimension());
    twist_model_.getStateJacobian(H, filter()->state(), true);
    twist_model_.getExpectedValue(x, filter()->state());

    TwistModel::MeasurementVector y(twist_model_.getDimension());
    TwistModel::NoiseVariance Iy(information);
    y.segment<3>(0) = update_linear;
    y.segment<3>(3) = update_angular;

    // invert Iy if information is a covariance matrix
    if (!interpret_covariance_as_information_matrix_) {
      ROS_DEBUG_NAMED("poseupdate", "Twist updates with covariance matrices are currently not supported");
      break;
    }

    // update VelocityXY
    if (information(0,0) > 0.0 || information(0,0) > 0.0) {
      status_flags_ |= STATE_VELOCITY_XY;

      // fixed_velocity_xy_stddev_ = 1.0;
      if (fixed_velocity_xy_stddev_ != 0.0) {
        for(int i = 0; i < 6; ++i) Iy(0,i) = Iy(1,i) = Iy(i,0) = Iy(i,1) = 0.0;
        Iy(0,0) = Iy(1,1) = 1.0 / (fixed_velocity_xy_stddev_*fixed_velocity_xy_stddev_);
      }
    }

    // update VelocityZ
    if (information(2,2) > 0.0) {
      status_flags_ |= STATE_VELOCITY_Z;

      // fixed_velocity_z_stddev_ = 1.0;
      if (fixed_velocity_z_stddev_ != 0.0) {
          for(int i = 0; i < 6; ++i) Iy(2,i) = Iy(i,2) = 0.0;
        Iy(2,2) = 1.0 / (fixed_velocity_z_stddev_*fixed_velocity_z_stddev_);
      }
    }

    // update RateXY
    if (information(3,3) > 0.0 || information(4,4) > 0.0) {
      status_flags_ |= STATE_RATE_XY;

      // fixed_angular_rate_xy_stddev_ = 1.0;
      if (fixed_angular_rate_xy_stddev_ != 0.0) {
        for(int i = 0; i < 6; ++i) Iy(3,i) = Iy(3,i) = Iy(i,4) = Iy(i,4) = 0.0;
        Iy(4,4) = Iy(5,5) = 1.0 / (fixed_angular_rate_xy_stddev_*fixed_angular_rate_xy_stddev_);
      }
    }

    // update RateZ
    if (information(5,5) > 0.0) {
      status_flags_ |= STATE_RATE_Z;

      // fixed_angular_rate_z_stddev_ = 1.0;
      if (fixed_angular_rate_z_stddev_ != 0.0) {
        for(int i = 0; i < 6; ++i) Iy(5,i) = Iy(i,5) = 0.0;
        Iy(5,5) = 1.0 / (fixed_angular_rate_z_stddev_*fixed_angular_rate_z_stddev_);
      }
    }

    ROS_DEBUG_STREAM_NAMED("poseupdate", "Twist Update: ");
    ROS_DEBUG_STREAM_NAMED("poseupdate", "      x = [" << x.transpose() << "], H = [ " << H << " ], Px = [" <<  (H * filter()->state().P() * H.transpose()) << "], Ix = [ " << (H * filter()->state().P() * H.transpose()).inverse() << "]");
    ROS_DEBUG_STREAM_NAMED("poseupdate", "      y = [" << y.transpose() << "], Iy = [ " << Iy << " ]");
    double innovation = updateInternal(filter()->state(), Iy, y - x, H, "twist", 0.0);
    twist_model_.getExpectedValue(x, filter()->state());
    ROS_DEBUG_STREAM_NAMED("poseupdate", " ==> xy = [" << x.transpose() << "], Pxy = [ " << (H * filter()->state().P() * H.transpose()) << " ], innovation = " << innovation);

    break;
  }

  // already done in Measurement::update()
  // filter()->state().updated();
  return true;
}

double PoseUpdate::calculateOmega(const SymmetricMatrix &Ix, const SymmetricMatrix &Iy) {
  double tr_x = Ix.trace();
  double tr_y = Iy.trace();
  return tr_y / (tr_x + tr_y);
}

template <typename MeasurementVector, typename MeasurementMatrix, typename NoiseVariance>
double PoseUpdate::updateInternal(State &state, const NoiseVariance &Iy, const MeasurementVector &error, const MeasurementMatrix &H, const std::string& text, const double max_error, JumpFunction jump_function) {
  NoiseVariance H_Px_HT(H * state.P() * H.transpose());

  if (H_Px_HT.determinant() <= 0) {
    ROS_DEBUG_STREAM("Ignoring poseupdate for " << text << " as the a-priori state covariance is zero!");
    return 0.0;
  }
  NoiseVariance Ix(H_Px_HT.inverse().eval());

  ROS_DEBUG_STREAM_NAMED("poseupdate", "H = [" << H << "]");
  ROS_DEBUG_STREAM_NAMED("poseupdate", "Ix = [" << Ix << "]");

  double alpha = fixed_alpha_, beta = fixed_beta_;
  if (alpha == 0.0 && beta == 0.0) {
    beta = calculateOmega(Ix, Iy);
    alpha = 1.0 - beta;

//    if (beta > 0.8) {
//      ROS_DEBUG_STREAM("Reducing update variance for " << text << " due to high information difference between Ix = [" << Ix << "] and Iy = [" << Iy << "]");
//      beta = 0.8;
//      alpha = 1.0 - beta;
//    }
  }
  ROS_DEBUG_STREAM_NAMED("poseupdate", "alpha = " << alpha << ", beta = " << beta);

  if (max_error > 0.0) {
    double error2 = (error.transpose() * Ix * (Ix + Iy).inverse() * Iy * error)(0);
    if (error2 > max_error * max_error) {
      if (!jump_on_max_error_ || !jump_function) {
        ROS_WARN_STREAM_NAMED("poseupdate", "Ignoring poseupdate for " << text << " as the error [ " << error.transpose() << " ], |error| = " << sqrt(error2) << " sigma exceeds max_error!");
        return 0.0;
      } else {
        ROS_WARN_STREAM_NAMED("poseupdate", "Update for " << text << " with error [ " << error.transpose() << " ], |error| = " << sqrt(error2) << " sigma exceeds max_error!");
        jump_function(state, error);
        return 0.0;
      }
    }
  }

//  SymmetricMatrix Ii(Ix * (alpha - 1) + Iy * beta);
//  double innovation = Ii.determinant();
//  ROS_DEBUG_STREAM_NAMED("poseupdate", "Ii = [" << Ii << "], innovation = " << innovation);

  // S_1 is equivalent to S^(-1) = (H*P*H^T + R)^(-1) in the standard Kalman gain
  NoiseVariance S_1(Ix - Ix * (Ix * alpha + Iy * beta).inverse() * Ix);
  typename Matrix_<State::Covariance::ColsAtCompileTime, MeasurementMatrix::RowsAtCompileTime>::type P_HT((H * state.P()).transpose());
  ROS_DEBUG_STREAM_NAMED("poseupdate", "P*HT = [" << (P_HT) << "]");

  double innovation = S_1.determinant();
  state.P() = state.P() - P_HT * S_1 * P_HT.transpose(); // may invalidate Px if &Pxy == &Px
  state.P().assertSymmetric();
  state.update(P_HT * Iy * beta * error);
  // state.x() = state.x() + P_HT * Iy * beta * error;

  ROS_DEBUG_STREAM_NAMED("poseupdate", "K = [" << (P_HT * Iy * beta) << "]");
  ROS_DEBUG_STREAM_NAMED("poseupdate", "dx = [" << (P_HT * Iy * beta * error).transpose() << "]");

  return innovation;
}

void PositionXYModel::getExpectedValue(MeasurementVector &y_pred, const State &state) {
  y_pred = state.getPosition().head<2>();
}

void PositionXYModel::getStateJacobian(MeasurementMatrix &C, const State &state, bool init) {
  if (init) {
    if (state.position()) {
      state.position()->cols(C)(0,X)   = 1.0;
      state.position()->cols(C)(1,Y)   = 1.0;
    }
  }
}

void PositionXYModel::updateState(State &state, const ColumnVector &diff) const {
  if (state.position()) {
    state.position()->vector().head<2>() += diff;
  }
}

void PositionZModel::getExpectedValue(MeasurementVector &y_pred, const State &state) {
  y_pred(0) = state.getPosition().z();
}

void PositionZModel::getStateJacobian(MeasurementMatrix &C, const State &state, bool init) {
  if (init && state.position()) {
    state.position()->cols(C)(0,Z)   = 1.0;
  }
}

void PositionZModel::updateState(State &state, const ColumnVector &diff) const {
  if (state.position()) {
    state.position()->vector().segment<1>(Z) += diff;
  }
}

void YawModel::getExpectedValue(MeasurementVector &y_pred, const State &state) {
  y_pred(0) = state.getYaw();
}

void YawModel::getStateJacobian(MeasurementMatrix &C, const State &state, bool init) {
  if (init && state.orientation()) {
    state.orientation()->cols(C)(0,Z) = 1.0;
  }
}

void YawModel::updateState(State &state, const ColumnVector &diff) const {
  Eigen::Quaterniond::Matrix3 rotation(Eigen::AngleAxisd(diff(0), Eigen::Vector3d::UnitZ()));

  Eigen::MatrixXd S(Eigen::MatrixXd::Identity(state.getCovarianceDimension(), state.getCovarianceDimension()));

  if (state.orientation()) {
//    S.block(state.orientation()->getCovarianceIndex(), state.orientation()->getCovarianceIndex(), 4, 4) <<
//      /* QUATERNION_X: */  rotation.w(), -rotation.z(),  rotation.y(), rotation.x(),
//      /* QUATERNION_Y: */  rotation.z(),  rotation.w(), -rotation.x(), rotation.y(),
//      /* QUATERNION_Z: */ -rotation.y(),  rotation.x(),  rotation.w(), rotation.z(),
//      /* QUATERNION_W: */ -rotation.x(), -rotation.y(), -rotation.z(), rotation.w();
    S.block(state.orientation()->getCovarianceIndex(), state.orientation()->getCovarianceIndex(), 3, 3) = rotation.transpose();
    state.updateOrientation(ColumnVector3(0.0, 0.0, -diff(0)));
  }

  if (state.velocity()) {
    S.block(state.velocity()->getCovarianceIndex(), state.velocity()->getCovarianceIndex(), 3, 3) = rotation.transpose();
    state.velocity()->vector() = rotation.transpose() * state.velocity()->vector();
  }

// Rate vector is in body frame. No need to rotate it.
//  if (state.rate()) {
//    S.block(state.rate()->getCovarianceIndex(), state.rate()->getCovarianceIndex(), 3, 3) = rotation.transpose();
//    state.rate()->vector() = rotation.transpose() * state.rate()->vector();
//  }

  ROS_DEBUG_STREAM_NAMED("poseupdate", "Jump yaw by " << (diff(0) * 180.0/M_PI) << " degrees. rotation = [" << rotation << "], S = [" << S << "].");

  // update covariance matrix P
  state.P() = S * state.P() * S.transpose();
}

void TwistModel::getExpectedValue(MeasurementVector &y_pred, const State &state) {
  y_pred.segment<3>(0) = state.getVelocity();
  y_pred.segment<3>(3) = state.getRate();
}

void TwistModel::getStateJacobian(MeasurementMatrix &C, const State &state, bool init) {
  if (init && state.velocity()) {
    state.velocity()->cols(C)(0,X) = 1.0;
    state.velocity()->cols(C)(1,Y) = 1.0;
    state.velocity()->cols(C)(2,Z) = 1.0;
  }

  if (init && state.rate()) {
    state.rate()->cols(C)(3,X) = 1.0;
    state.rate()->cols(C)(4,Y) = 1.0;
    state.rate()->cols(C)(5,Z) = 1.0;
  }
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/measurements/gravity.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/measurements/gravity.h>
#include <hector_pose_estimation/pose_estimation.h>
#include <hector_pose_estimation/filter/set_filter.h>

namespace hector_pose_estimation {

template class Measurement_<GravityModel>;

GravityModel::GravityModel()
  : gravity_(MeasurementVector::Zero())
{
  parameters().add("stddev", stddev_, 1.0);
  parameters().add("use_bias", use_bias_, std::string("accelerometer_bias"));
}

GravityModel::~GravityModel() {}

bool GravityModel::init(PoseEstimation &estimator, Measurement &measurement, State &state) {
  if (!use_bias_.empty()) {
    bias_ = state.getSubState<3,3>(use_bias_);
    if (!bias_) {
      ROS_ERROR("Could not find bias substate '%s' during initialization of gravity measurement '%s'.", use_bias_.c_str(), measurement.getName().c_str());
      return false;
    }
  } else {
    bias_.reset();
  }

  setGravity(estimator.parameters().getAs<double>("gravity_magnitude"));
  return true;
}


void GravityModel::getMeasurementNoise(NoiseVariance& R, const State&, bool init)
{
  if (init) {
    R(0,0) = R(1,1) = R(2,2) = pow(stddev_, 2);
  }
}

void GravityModel::getExpectedValue(MeasurementVector& y_pred, const State& state)
{
  const State::RotationMatrix &R = state.R();
  y_pred = -R.row(2).transpose() * gravity_.z();
  if (bias_) {
    y_pred += bias_->getVector();
  }
}

void GravityModel::getStateJacobian(MeasurementMatrix& C, const State& state, bool)
{
  const State::RotationMatrix &R = state.R();
  if (state.orientation()) {
//    C(0,state.getOrientationCovarianceIndex() + W) =  gravity_.z()*2*q.y();
//    C(0,state.getOrientationCovarianceIndex() + X) = -gravity_.z()*2*q.z();
//    C(0,state.getOrientationCovarianceIndex() + Y) =  gravity_.z()*2*q.w();
//    C(0,state.getOrientationCovarianceIndex() + Z) = -gravity_.z()*2*q.x();
//    C(1,state.getOrientationCovarianceIndex() + W) = -gravity_.z()*2*q.x();
//    C(1,state.getOrientationCovarianceIndex() + X) = -gravity_.z()*2*q.w();
//    C(1,state.getOrientationCovarianceIndex() + Y) = -gravity_.z()*2*q.z();
//    C(1,state.getOrientationCovarianceIndex() + Z) = -gravity_.z()*2*q.y();
//    C(2,state.getOrientationCovarianceIndex() + W) = -gravity_.z()*2*q.w();
//    C(2,state.getOrientationCovarianceIndex() + X) =  gravity_.z()*2*q.x();
//    C(2,state.getOrientationCovarianceIndex() + Y) =  gravity_.z()*2*q.y();
//    C(2,state.getOrientationCovarianceIndex() + Z) = -gravity_.z()*2*q.z();

     state.orientation()->cols(C)(X,X) = -gravity_.z() * R(1,0);
     state.orientation()->cols(C)(X,Y) =  gravity_.z() * R(0,0);
     state.orientation()->cols(C)(Y,X) = -gravity_.z() * R(1,1);
     state.orientation()->cols(C)(Y,Y) =  gravity_.z() * R(0,1);
     state.orientation()->cols(C)(Z,X) = -gravity_.z() * R(1,2);
     state.orientation()->cols(C)(Z,Y) =  gravity_.z() * R(0,2);
  }

//  Only the bias component in direction of the gravity is observable, under the assumption that we do not accelerate vertically.
  if (bias_) {
    bias_->cols(C) = R.row(2).transpose() * R.row(2);
  }
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/measurements/gps.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/measurements/gps.h>
#include <hector_pose_estimation/global_reference.h>
#include <hector_pose_estimation/filter/set_filter.h>

namespace hector_pose_estimation {

template class Measurement_<GPSModel>;

GPSModel::GPSModel()
{
  position_stddev_ = 10.0;
  velocity_stddev_ = 1.0;
  parameters().add("position_stddev", position_stddev_);
  parameters().add("velocity_stddev", velocity_stddev_);
}

GPSModel::~GPSModel() {}

void GPSModel::getMeasurementNoise(NoiseVariance& R, const State&, bool init)
{
  if (init) {
    R(0,0) = R(1,1) = pow(position_stddev_, 2);
    R(2,2) = R(3,3) = pow(velocity_stddev_, 2);
  }
}

bool GPSModel::prepareUpdate(State &state, const MeasurementUpdate &update)
{
  state.getRotationMatrix(R);
  return true;
}

void GPSModel::getExpectedValue(MeasurementVector& y_pred, const State& state)
{
  y_pred(0) = state.getPosition().x();
  y_pred(1) = state.getPosition().y();
  y_pred(2) = state.getVelocity().x();
  y_pred(3) = state.getVelocity().y();
}

void GPSModel::getStateJacobian(MeasurementMatrix& C, const State& state, bool init)
{
  if (!init) return; // C is time-constant

  if (state.position()) {
    state.position()->cols(C)(0,X) = 1.0;
    state.position()->cols(C)(1,Y) = 1.0;
  }

  if (state.velocity()) {
    state.velocity()->cols(C)(2,X) = 1.0;
    state.velocity()->cols(C)(3,Y) = 1.0;
  }
}

GPS::GPS(const std::string &name)
  : Measurement_<GPSModel>(name)
  , auto_reference_(true)
  , y_(4)
{
  parameters().add("auto_reference", auto_reference_);
}

GPS::~GPS()
{}

void GPS::onReset() {
  reference_.reset();
}

GPSModel::MeasurementVector const& GPS::getVector(const GPSUpdate &update, const State&) {
  if (!reference_) {
    y_.setConstant(0.0/0.0);
    return y_;
  }

  reference_->fromWGS84(update.latitude, update.longitude, y_(0), y_(1));
  reference_->fromNorthEast(update.velocity_north, update.velocity_east, y_(2), y_(3));

  return y_;
}

bool GPS::prepareUpdate(State &state, const Update &update) {
  // reset reference position if GPS has not been updated for a while
  if (timedout()) reference_.reset();

  // find new reference position
  if (reference_ != GlobalReference::Instance()) {
    reference_ = GlobalReference::Instance();
    if (!auto_reference_ && !reference_->hasPosition()) {
      ROS_ERROR("Cannot use GPS measurements if no reference latitude/longitude is set and %s/auto_reference parameter is false.", name_.c_str());
      return false;
    }
    if (auto_reference_) reference_->setCurrentPosition(state, update.latitude, update.longitude);
  }

  return true;
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/measurements/heading.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/measurements/heading.h>
#include <hector_pose_estimation/filter/set_filter.h>

namespace hector_pose_estimation {

template class Measurement_<HeadingModel>;

HeadingModel::HeadingModel()
{
  parameters().add("stddev", stddev_, 10.0*M_PI/180.0);
}

HeadingModel::~HeadingModel() {}

void HeadingModel::getMeasurementNoise(NoiseVariance& R, const State&, bool init)
{
  if (init) {
    R(0,0) = pow(stddev_, 2);
  }
}

void HeadingModel::getExpectedValue(MeasurementVector& y_pred, const State& state)
{
  y_pred(0) = state.getYaw();
}

void HeadingModel::getStateJacobian(MeasurementMatrix& C, const State& state, bool)
{
  if (state.orientation()) {
    state.orientation()->cols(C)(0,Z) = 1.0;
  }
}

void HeadingModel::limitError(MeasurementVector &error) {
  error(0) = remainder(error(0), 2 * M_PI);
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/measurements/height.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/measurements/height.h>
#include <hector_pose_estimation/pose_estimation.h>
#include <hector_pose_estimation/global_reference.h>
#include <hector_pose_estimation/filter/set_filter.h>

#include <ros/console.h>

namespace hector_pose_estimation {

template class Measurement_<HeightModel>;

HeightModel::HeightModel()
{
  stddev_ = 10.0;
  elevation_ = 0.0;
  parameters().add("stddev", stddev_);
}

HeightModel::~HeightModel() {}

void HeightModel::getMeasurementNoise(NoiseVariance& R, const State&, bool init)
{
  if (init) {
    R(0,0) = pow(stddev_, 2);
  }
}

void HeightModel::getExpectedValue(MeasurementVector& y_pred, const State& state)
{
  y_pred(0) = state.getPosition().z() + getElevation();
}

void HeightModel::getStateJacobian(MeasurementMatrix& C, const State& state, bool init)
{
  if (!init) return; // C is time-constant

  if (state.position()) {
    state.position()->cols(C)(0,Z) = 1.0;
  }
}

HeightBaroCommon::HeightBaroCommon(Measurement* parent)
  : auto_elevation_(true)
  , elevation_initialized_(false)
{
}

HeightBaroCommon::~HeightBaroCommon() {}

void HeightBaroCommon::onReset() {
  elevation_initialized_ = false;
}

double HeightBaroCommon::resetElevation(const State &state, boost::function<double()> altitude_func) {
  if (!elevation_initialized_) {
    if (auto_elevation_) GlobalReference::Instance()->setCurrentAltitude(state, altitude_func());
    elevation_initialized_ = true;
  }

  return GlobalReference::Instance()->position().altitude;
}

Height::Height(const std::string &name)
  : Measurement_<HeightModel>(name)
  , HeightBaroCommon(this)
{
  parameters().add("auto_elevation", auto_elevation_);
}

void Height::onReset() {
  HeightBaroCommon::onReset();
}

template <typename T> struct functor_wrapper
{
  functor_wrapper(const T& value) : value(value) {}
  T& operator()() { return value; }
  const T& operator()() const { return value; }
private:
  T value;
};

bool Height::prepareUpdate(State &state, const Update &update) {
  setElevation(resetElevation(state, functor_wrapper<double>(update.getVector()(0))));
  return true;
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/measurements/magnetic.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/measurements/magnetic.h>
#include <hector_pose_estimation/filter/set_filter.h>

#include <Eigen/Geometry>

namespace hector_pose_estimation {

template class Measurement_<MagneticModel>;

MagneticModel::MagneticModel()
  : declination_(0.0), inclination_(60.0 * M_PI/180.0), magnitude_(0.0)
{
  parameters().add("stddev", stddev_, 1.0);
  parameters().add("declination", declination_);
  parameters().add("inclination", inclination_);
  parameters().add("magnitude", magnitude_);
}

MagneticModel::~MagneticModel() {}

bool MagneticModel::init(PoseEstimation &estimator, Measurement &measurement, State &state)
{
  updateMagneticField();
  return true;
}

void MagneticModel::setReference(const GlobalReference::Heading &reference_heading) {
  magnetic_field_reference_.x() = reference_heading.cos * magnetic_field_north_.x() - reference_heading.sin * magnetic_field_north_.y();
  magnetic_field_reference_.y() = reference_heading.sin * magnetic_field_north_.x() + reference_heading.cos * magnetic_field_north_.y();
  magnetic_field_reference_.z() = magnetic_field_north_.z();
}

void MagneticModel::getMeasurementNoise(NoiseVariance& R, const State&, bool init)
{
  if (init) {
    R(0,0) = R(1,1) = R(2,2) = pow(stddev_, 2);
  }
}

void MagneticModel::getExpectedValue(MeasurementVector& y_pred, const State& state)
{
  const State::RotationMatrix &R = state.R();
  y_pred = R.transpose() * magnetic_field_reference_;
}

void MagneticModel::getStateJacobian(MeasurementMatrix& C, const State& state, bool)
{
  if (state.orientation()) {
    const State::RotationMatrix &R = state.R();
    state.orientation()->cols(C)(X,Z) = R(0,0) * magnetic_field_reference_.y() - R(1,0) * magnetic_field_reference_.x();
    state.orientation()->cols(C)(Y,Z) = R(0,1) * magnetic_field_reference_.y() - R(1,1) * magnetic_field_reference_.x();
    state.orientation()->cols(C)(Z,Z) = R(0,2) * magnetic_field_reference_.y() - R(1,2) * magnetic_field_reference_.x();
  }
}

double MagneticModel::getMagneticHeading(const State& state, const MeasurementVector &y) const {
  MeasurementVector y_nav;
  y_nav = state.R() * y;
  return atan2(y_nav.y(), y_nav.x()) - state.getYaw();
}

double MagneticModel::getTrueHeading(const State& state, const MeasurementVector &y) const {
  return getMagneticHeading(state, y) + declination_;
}

void MagneticModel::updateMagneticField()
{
  double cos_inclination, sin_inclination;
  sincos(inclination_, &sin_inclination, &cos_inclination);

  double cos_declination, sin_declination;
  sincos(declination_, &sin_declination, &cos_declination);

  // return normalized magnetic field if magnitude is zero
  double magnitude = magnitude_;
  if (magnitude == 0.0) magnitude = 1.0;

  magnetic_field_north_.x() = magnitude * (cos_inclination *   cos_declination);
  magnetic_field_north_.y() = magnitude * (cos_inclination * (-sin_declination));
  magnetic_field_north_.z() = magnitude * (-sin_inclination);
}

Magnetic::Magnetic(const std::string &name)
  : Measurement_<MagneticModel>(name)
  , auto_heading_(true)
  , deviation_(3)
{
  deviation_.setZero();
  parameters().add("auto_heading", auto_heading_);
  parameters().add("deviation", deviation_);
}

void Magnetic::onReset() {
  reference_.reset();
}

const MagneticModel::MeasurementVector& Magnetic::getVector(const Magnetic::Update& update, const State& state) {
  y_ = Measurement_<MagneticModel>::getVector(update, state) + deviation_;
  if (getModel()->hasMagnitude()) return y_;

  double norm = y_.norm();
  if (norm < 1e-5) {
    y_.setZero();
  } else {
    y_ = y_ / norm;
  }
  return y_;
}

//const MagneticModel::NoiseVariance& Magnetic::getVariance(const Magnetic::Update& update, const State& state) {
//  if (getModel()->hasMagnitude()) return Measurement_<MagneticModel>::getVariance(update, state);

//  R_ = Measurement_<MagneticModel>::getVariance(update, state);
//  double norm = Measurement_<MagneticModel>::getVector(update, state).norm();
//  if (norm < 1e-5) {
//    R_ = NoiseVariance(1.0, 1.0, 1.0);
//  } else {
//    R_ =  R_ / (norm*norm);
//  }
//  return R_;
//}

bool Magnetic::prepareUpdate(State &state, const Update &update) {
  // reset reference position if Magnetic has not been updated for a while
  if (timedout()) reference_.reset();

  if (reference_ != GlobalReference::Instance()) {
    reference_ = GlobalReference::Instance();
    if (auto_heading_) reference_->setCurrentHeading(state, getModel()->getTrueHeading(state, update.getVector()));
  }

  getModel()->setReference(reference_->heading());
  return true;
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/measurements/rate.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/measurements/rate.h>
#include <hector_pose_estimation/filter/set_filter.h>

namespace hector_pose_estimation {

template class Measurement_<RateModel>;

RateModel::RateModel()
{
  parameters().add("stddev", stddev_, 10.0 * M_PI/180.0);
  parameters().add("use_bias", use_bias_, std::string("gyro_bias"));
}

RateModel::~RateModel() {}

bool RateModel::init(PoseEstimation &estimator, Measurement &measurement, State &state)
{
  if (!use_bias_.empty()) {
    bias_ = state.getSubState<3,3>(use_bias_);
    if (!bias_) {
      ROS_ERROR("Could not find bias substate '%s' during initialization of rate measurement '%s'.", use_bias_.c_str(), measurement.getName().c_str());
      return false;
    }
  } else {
    bias_.reset();
  }

  return true;
}

void RateModel::getMeasurementNoise(NoiseVariance& R, const State&, bool init)
{
  if (init) {
    R(0,0) = R(1,1) = R(2,2) = pow(stddev_, 2);
  }
}

void RateModel::getExpectedValue(MeasurementVector& y_pred, const State& state)
{
  y_pred = state.getRate();

  if (bias_) {
    y_pred += bias_->getVector();
  }
}

void RateModel::getStateJacobian(MeasurementMatrix &C, const State &state, bool init)
{
  if (!init) return;

  if (state.rate()) {
   state.rate()->cols(C).setIdentity();
  }

  if (bias_) {
    bias_->cols(C).setIdentity();
  }
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/measurements/zerorate.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/measurements/zerorate.h>
#include <hector_pose_estimation/system/imu_model.h>
#include <hector_pose_estimation/filter/set_filter.h>

namespace hector_pose_estimation {

template class Measurement_<ZeroRateModel>;

ZeroRateModel::ZeroRateModel()
{
  parameters().add("stddev", stddev_, 90.0*M_PI/180.0);
  parameters().add("use_bias", use_bias_, std::string("gyro_bias"));
}

ZeroRateModel::~ZeroRateModel() {}

bool ZeroRateModel::init(PoseEstimation &estimator, Measurement &measurement, State &state)
{
  if (!use_bias_.empty()) {
    bias_ = state.getSubState<3,3>(use_bias_);
    if (!bias_) {
      ROS_ERROR("Could not find bias substate '%s' during initialization of zero rate pseudo measurement '%s'.", use_bias_.c_str(), measurement.getName().c_str());
      return false;
    }
  } else {
    bias_.reset();
  }

  if (!bias_ && !state.rate()) {
    ROS_WARN("Pseudo updating with zero rate is a no-op, as the state does not contain rates nor biases.");
    // return false;
  }

  return true;
}

void ZeroRateModel::getMeasurementNoise(NoiseVariance& R, const State&, bool init)
{
  if (init) {
    R(0,0) = pow(stddev_, 2);
  }
}

void ZeroRateModel::getExpectedValue(MeasurementVector& y_pred, const State& state)
{
  y_pred(0) = state.getRate().z();

  if (!state.rate() && bias_) {
    y_pred(0) += bias_->getVector().z();
  }
}

void ZeroRateModel::getStateJacobian(MeasurementMatrix& C, const State& state, bool)
{
  if (state.rate()) {
    state.rate()->cols(C)(0,Z) = 1.0;
  } else if (bias_) {
    bias_->cols(C)(0,Z) = 1.0;
  }
}

const ZeroRateModel::MeasurementVector* ZeroRateModel::getFixedMeasurementVector() const
{
  static MeasurementVector zero(MeasurementVector::Zero());
  return &zero;
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/state.cpp
//=================================================================================================
// Copyright (c) 2013, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/state.h>
#include <hector_pose_estimation/substate.h>

#include <ros/console.h>

namespace hector_pose_estimation {

FullState::FullState()
{
  orientation_ = addSubState<OrientationStateType::VectorDimension,OrientationStateType::CovarianceDimension>("orientation");
  rate_ = addSubState<RateStateType::VectorDimension,RateStateType::CovarianceDimension>("rate");
  position_ = addSubState<PositionStateType::VectorDimension,PositionStateType::CovarianceDimension>("position");
  velocity_ = addSubState<VelocityStateType::VectorDimension,VelocityStateType::CovarianceDimension>("velocity");
  construct();
}

FullState::~FullState() {}

OrientationPositionVelocityState::OrientationPositionVelocityState()
{
  orientation_ = addSubState<OrientationStateType::VectorDimension,OrientationStateType::CovarianceDimension>("orientation");
  position_ = addSubState<PositionStateType::VectorDimension,PositionStateType::CovarianceDimension>("position");
  velocity_ = addSubState<VelocityStateType::VectorDimension,VelocityStateType::CovarianceDimension>("velocity");
  construct();
}

OrientationPositionVelocityState::~OrientationPositionVelocityState() {}

OrientationOnlyState::OrientationOnlyState()
{
  orientation_ = addSubState<OrientationStateType::VectorDimension,OrientationStateType::CovarianceDimension>("orientation");
  construct();
}

OrientationOnlyState::~OrientationOnlyState() {}

PositionVelocityState::PositionVelocityState()
{
  position_ = addSubState<PositionStateType::VectorDimension,PositionStateType::CovarianceDimension>("position");
  velocity_ = addSubState<VelocityStateType::VectorDimension,VelocityStateType::CovarianceDimension>("velocity");
  construct();
}

PositionVelocityState::~PositionVelocityState() {}

State::State() {}

State::~State() {}

void State::construct()
{
  base_.reset(new BaseState(*this, getVectorDimension(), getCovarianceDimension()));
  reset();
}

void State::reset()
{
  // reset status flags
  system_status_ = 0;
  measurement_status_ = 0;

  // reset pseudo-states
  fake_rate_ = Vector::Zero(3,1);
  fake_orientation_ = Vector::Zero(4,1);
  fake_position_ = Vector::Zero(3,1);
  fake_velocity_ = Vector::Zero(3,1);
  fake_acceleration_ = Vector::Zero(3,1);

  // reset state
  vector_.setZero();
  covariance_.setZero();
  fake_orientation_.w() = 1.0;
  if (orientation()) orientation()->vector().w() = 1.0;

  R_valid_ = false;
}

bool State::valid() const {
  return (vector_ == vector_);
}

void State::updated()
{
  normalize();
  R_valid_ = false;
}

State::ConstOrientationType State::getOrientation() const   { return (orientation()  ? orientation()->getVector()  : ConstOrientationType(fake_orientation_, 0)); }
State::ConstRateType State::getRate() const                 { return (rate()         ? rate()->getVector()         : ConstRateType(fake_rate_, 0)); }
State::ConstPositionType State::getPosition() const         { return (position()     ? position()->getVector()     : ConstPositionType(fake_position_, 0)); }
State::ConstVelocityType State::getVelocity() const         { return (velocity()     ? velocity()->getVector()     : ConstVelocityType(fake_velocity_, 0)); }
State::ConstAccelerationType State::getAcceleration() const { return (acceleration() ? acceleration()->getVector() : ConstAccelerationType(fake_acceleration_, 0)); }

void State::getRotationMatrix(RotationMatrix &R) const
{
  Quaternion q(getOrientation());
  R = q.toRotationMatrix();
//  R << (q.w()*q.w()+q.x()*q.x()-q.y()*q.y()-q.z()*q.z()), (2.0*q.x()*q.y()-2.0*q.w()*q.z()),                 (2.0*q.x()*q.z()+2.0*q.w()*q.y()),
//       (2.0*q.x()*q.y()+2.0*q.w()*q.z()),                 (q.w()*q.w()-q.x()*q.x()+q.y()*q.y()-q.z()*q.z()), (2.0*q.y()*q.z()-2.0*q.w()*q.x()),
//       (2.0*q.x()*q.z()-2.0*q.w()*q.y()),                 (2.0*q.y()*q.z()+2.0*q.w()*q.x()),                 (q.w()*q.w()-q.x()*q.x()-q.y()*q.y()+q.z()*q.z());
}

const State::RotationMatrix &State::R() const {
  if (!R_valid_) {
    getRotationMatrix(R_);
    R_valid_ = true;
  }
  return R_;
}

double State::getYaw() const
{
  ConstOrientationType q(getOrientation());
  return atan2(2*q.x()*q.y() + 2*q.w()*q.z(), q.x()*q.x() + q.w()*q.w() - q.z()*q.z() - q.y()*q.y());
}

void State::getEuler(double &roll, double &pitch, double &yaw) const
{
  ConstOrientationType q(getOrientation());
  roll  = atan2(2*(q.y()*q.z() + q.w()*q.x()), q.w()*q.w() - q.x()*q.x() - q.y()*q.y() + q.z()*q.z());
  pitch = -asin(2*(q.x()*q.z() - q.w()*q.y()));
  yaw   = atan2(2*(q.x()*q.y() + q.w()*q.z()), q.w()*q.w() + q.x()*q.x() - q.y()*q.y() - q.z()*q.z());
}

ColumnVector3 State::getEuler() const
{
  ColumnVector3 euler;
  getEuler(euler(0), euler(1), euler(2));
  return euler;
}

void State::update(const Vector &vector_update)
{
  if (orientation()) {
    // if vector_update has only n - 1 elements, we have to use covariance indices
    int orientation_index, orientation_size;
    if (vector_update.size() == getVectorDimension() - 1) {
      orientation_index = orientation()->getCovarianceIndex();
      orientation_size  = orientation()->getCovarianceDimension();
    } else {
      assert(vector_update.size() == getVectorDimension());
      orientation_index = orientation()->getVectorIndex();
      orientation_size  = orientation()->getVectorDimension();
    }

    // add everything before orientation part
    if (orientation_index > 0) {
      int length = orientation_index;
      x().head(length) += vector_update.head(length);
    }

    // add everything after orientation part
    if (orientation_index + orientation_size < vector_update.size()) {
      int length = vector_update.size() - orientation_index - orientation_size;
      x().tail(length) += vector_update.tail(length);
    }

    // update orientation
    updateOrientation(vector_update.segment<3>(orientation_index));

  } else {
    // simply add vectors
    x() += vector_update;
  }
}

void State::updateOrientation(const ColumnVector3 &rotation_vector)
{
  if (!orientation()) return;

//  Eigen::Quaterniond q(orientation()->vector().data());
//  q = Eigen::Quaterniond(1.0, 0.5 * rotation_vector.x(), 0.5 * rotation_vector.y(), 0.5 * rotation_vector.z()) * q;
//  orientation()->vector() = q.normalized().coeffs();

  // Eigen::Map<Eigen::Quaterniond> q(orientation()->vector().data());
  Eigen::Quaterniond q(orientation()->vector().data());
  q = Eigen::Quaterniond().fromRotationVector(rotation_vector) * q;
  orientation()->vector() = q.coeffs();

  R_valid_ = false;
}

bool State::inSystemStatus(SystemStatus test_status) const {
  return (getSystemStatus() & test_status) == test_status;
}

bool State::setSystemStatus(SystemStatus new_status) {
  if (new_status == system_status_) return true;

  // iterate through StatusCallbacks
  for(std::vector<SystemStatusCallback>::const_iterator it = status_callbacks_.begin(); it != status_callbacks_.end(); ++it)
    if (!(*it)(new_status)) return false;

  SystemStatus set = new_status & ~system_status_;
  SystemStatus cleared = system_status_ & ~new_status;
  if (set)     ROS_INFO_STREAM("Set system status " << getSystemStatusString(new_status, set));
  if (cleared) ROS_INFO_STREAM("Cleared system status " << getSystemStatusString(cleared, cleared));

  system_status_ = new_status;
  return true;
}

bool State::setMeasurementStatus(SystemStatus new_measurement_status) {
  SystemStatus set = new_measurement_status & ~measurement_status_;
  SystemStatus cleared = measurement_status_ & ~new_measurement_status;
  if (set)     ROS_INFO_STREAM("Set measurement status " << getSystemStatusString(new_measurement_status, set));
  if (cleared) ROS_INFO_STREAM("Cleared measurement status " << getSystemStatusString(cleared, cleared));

  measurement_status_ = new_measurement_status;
  return true;
}

bool State::updateSystemStatus(SystemStatus set, SystemStatus clear) {
  return setSystemStatus((system_status_ & ~clear) | set);
}

bool State::updateMeasurementStatus(SystemStatus set, SystemStatus clear) {
  return setMeasurementStatus((measurement_status_ & ~clear) | set);
}

void State::addSystemStatusCallback(const SystemStatusCallback& callback) {
//  for(std::vector<SystemStatusCallback>::const_iterator it = status_callbacks_.begin(); it != status_callbacks_.end(); ++it)
//    if (*it == callback) return;
  status_callbacks_.push_back(callback);
}

void State::normalize() {
  if (orientation()) {
    double s = 1.0 / orientation()->vector().norm();
    orientation()->vector() = orientation()->vector() * s;
  }
}

void State::setOrientation(const Quaternion& orientation)
{
  setOrientation(orientation.coeffs());
}

void State::setRollPitch(const Quaternion& q)
{
  ScalarType roll, pitch;
  roll  = atan2(2*(q.y()*q.z() + q.w()*q.x()), q.w()*q.w() - q.x()*q.x() - q.y()*q.y() + q.z()*q.z());
  pitch = -asin(2*(q.x()*q.z() - q.w()*q.y()));
  setRollPitch(roll, pitch);
}

void State::setRollPitch(ScalarType roll, ScalarType pitch)
{
  ScalarType yaw = getYaw();
  fake_orientation_ = Quaternion(Eigen::AngleAxis<ScalarType>(yaw, ColumnVector3::UnitZ()) * Eigen::AngleAxis<ScalarType>(pitch, ColumnVector3::UnitY()) * Eigen::AngleAxis<ScalarType>(roll, ColumnVector3::UnitX())).coeffs();
}

void State::setYaw(const Quaternion& orientation)
{
  const Quaternion &q = orientation;
  double yaw = atan2(2*(q.x()*q.y() + q.w()*q.z()), q.w()*q.w() + q.x()*q.x() - q.y()*q.y() - q.z()*q.z());
  setYaw(yaw);
}

void State::setYaw(ScalarType yaw)
{
  ColumnVector3 euler = getEuler();
  fake_orientation_ = Quaternion(Eigen::AngleAxis<ScalarType>(yaw, ColumnVector3::UnitZ()) * Eigen::AngleAxis<ScalarType>(euler(1), ColumnVector3::UnitY()) * Eigen::AngleAxis<ScalarType>(euler(2), ColumnVector3::UnitX())).coeffs();
}

template class SubState::initializer<Dynamic,Dynamic>;
template class SubState_<Dynamic,Dynamic>;

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/types.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/types.h>

namespace hector_pose_estimation {

std::string getSystemStatusString(const SystemStatus& status, const SystemStatus& asterisk_status) {
  std::string result;
  static const char* const str[] = {
    "ALIGNMENT", "DEGRADED", "READY", 0,
    "ROLLPITCH", "YAW", "PSEUDO_ROLLPITCH", "PSEUDO_YAW",
    "RATE_XY", "RATE_Z", "PSEUDO_RATE_XY", "PSEUDO_RATE_Z",
    "VELOCITY_XY", "VELOCITY_Z", "PSEUDO_VELOCITY_XY", "PSEUDO_VELOCITY_Z",
    "POSITION_XY", "POSITION_Z", "PSEUDO_POSITION_XY", "PSEUDO_POSITION_Z",
  };

  if (asterisk_status) {
      for(unsigned int i = 0; i < sizeof(str)/sizeof(*str); ++i) {
      if (asterisk_status & (1 << i)) {
        result += "*" + std::string(str[i]) + " ";
      }
    }

    if (asterisk_status != status) result += "(";
  }

  for(unsigned int i = 0; i < sizeof(str)/sizeof(*str); ++i) {
    if (status & (1 << i)) {
      if (asterisk_status & (1 << i)) continue;
      result += std::string(str[i]) + " ";
    }
  }
  if (result.size() > 0) result.resize(result.size() - 1);

  if (asterisk_status && (asterisk_status != status)) {
    result += ")";
  }

  return result;
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/filter/ekf.cpp
//=================================================================================================
// Copyright (c) 2013, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/filter/ekf.h>
#include <hector_pose_estimation/system.h>

#include <boost/pointer_cast.hpp>

#ifdef USE_HECTOR_TIMING
  #include <hector_diagnostics/timing.h>
#endif

namespace hector_pose_estimation {
namespace filter {

EKF::EKF(State &state)
  : Filter(state)
{}

EKF::~EKF()
{}

bool EKF::init(PoseEstimation &estimator)
{
  x_diff = State::Vector(state_.getVectorDimension());
  A = State::SystemMatrix(state_.getCovarianceDimension(), state_.getCovarianceDimension());
  Q = State::Covariance(state_.getCovarianceDimension(), state_.getCovarianceDimension());
  return true;
}

bool EKF::preparePredict(double dt)
{
  x_diff.setZero();
  A.setIdentity();
  Q.setZero();
  return Filter::preparePredict(dt);
}

bool EKF::predict(const SystemPtr& system, double dt)
{
  if (!Filter::predict(system, dt)) return false;
  EKF::Predictor *predictor = boost::dynamic_pointer_cast<EKF::Predictor>(system->predictor());
  x_diff += predictor->x_diff;
  A += predictor->A;
  Q += predictor->Q;
  return true;
}

bool EKF::doPredict(double dt) {
  ROS_DEBUG_NAMED("ekf.prediction", "EKF prediction (dt = %f):", dt);

  ROS_DEBUG_STREAM_NAMED("ekf.prediction", "A      = [" << std::endl << A << "]");
  ROS_DEBUG_STREAM_NAMED("ekf.prediction", "Q      = [" << std::endl << Q << "]");

#ifdef USE_HECTOR_TIMING
  { hector_diagnostics::TimingSection section("predict.ekf.covariance");
#endif
  state().P() = A * state().P() * A.transpose() + Q;
  state().P().assertSymmetric();

#ifdef USE_HECTOR_TIMING
  }
  { hector_diagnostics::TimingSection section("predict.ekf.state");
#endif
  state().update(x_diff);

#ifdef USE_HECTOR_TIMING
  }
#endif

  ROS_DEBUG_STREAM_NAMED("ekf.prediction", "x_pred = [" << state().getVector().transpose() << "]");
  ROS_DEBUG_STREAM_NAMED("ekf.prediction", "P_pred = [" << std::endl << state().getCovariance() << "]");

  Filter::doPredict(dt);
  return true;
}

} // namespace filter
} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/filter.cpp
//=================================================================================================
// Copyright (c) 2013, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/filter.h>
#include <hector_pose_estimation/pose_estimation.h>

#ifdef USE_HECTOR_TIMING
  #include <hector_diagnostics/timing.h>
#endif

namespace hector_pose_estimation {

Filter::Filter(State &state)
  : state_(state)
{
}

Filter::~Filter()
{
}

bool Filter::init(PoseEstimation& estimator)
{
  return true;
}

void Filter::cleanup()
{
}

void Filter::reset()
{
  state_.reset();
}

bool Filter::preparePredict(double)
{
  return true;
}

bool Filter::predict(const Systems& systems, double dt) {
  bool result = true;

#ifdef USE_HECTOR_TIMING
  hector_diagnostics::TimingSection section("predict");
#endif

  if (!preparePredict(dt)) return false;

  // Iterate through system models. For an EKF, this will populate the x_diff vector, A and Q matrices.
  for(Systems::iterator it = systems.begin(); it != systems.end(); it++) {
    const SystemPtr& system = *it;
    result &= predict(system, dt);
  }

  // Call the filter's global predict method. This will actually calculate the updated state vector and variance.
  result &= doPredict(dt);

  return result;
}

bool Filter::predict(const SystemPtr& system, double dt) {
#ifdef USE_HECTOR_TIMING
  hector_diagnostics::TimingSection section("predict." + system->getName());
#endif
  return system->update(dt);
}

bool Filter::doPredict(double dt) {
  // already done in System::update()
  // state_.updated();
  return true;
}

bool Filter::prepareCorrect()
{
  return true;
}

bool Filter::correct(const Measurements& measurements) {
  bool result = true;

#ifdef USE_HECTOR_TIMING
  hector_diagnostics::TimingSection section("correct");
#endif

  if (!prepareCorrect()) return false;

  // Iterate through measurement models. This will process the correction step directly.
  for(Measurements::iterator it = measurements.begin(); it != measurements.end(); it++) {
    const MeasurementPtr& measurement = *it;
    result &= correct(measurement);
  }

  // Call the filter's global correct method. No-op for EKF.
  result &= doCorrect();

  return result;
}

bool Filter::correct(const MeasurementPtr& measurement) {
#ifdef USE_HECTOR_TIMING
  hector_diagnostics::TimingSection section("correct." + measurement->getName());
#endif
  return measurement->process();
}

bool Filter::doCorrect() {
  // already done in Measurement::update()
  // state_.updated();
  return true;
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/system.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/system.h>
#include <hector_pose_estimation/filter/ekf.h>

namespace hector_pose_estimation {

System::System(const std::string &name)
  : name_(name)
  , status_flags_(0)
{
}

System::~System() {
}

void System::getPrior(State &state) const
{
  getModel()->getPrior(state);
}

bool System::init(PoseEstimation& estimator, State& state)
{
  if (!getModel() || !getModel()->init(estimator, *this, state)) return false;
  return true;
}

void System::cleanup()
{
  if (getModel()) getModel()->cleanup();
}

void System::reset(State& state)
{
  if (getModel()) getModel()->reset(state);
  status_flags_ = 0;
}

bool System::active(const State& state) {
  bool active = (!getModel() || getModel()->active(state));
  if (!active) status_flags_ = 0;
  return active;
}

bool System::update(double dt) {
  if (!filter() || !active(filter()->state())) return false;

  if (getModel()) status_flags_ = getModel()->getStatusFlags(filter()->state());
  if (!this->updateImpl(dt)) return false;
  filter()->state().updated();

  updated();
  return true;
}

void System::updated() {
}

bool System::limitState(State &state) {
  return getModel()->limitState(state);
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/system/ground_vehicle_model.cpp
//=================================================================================================
// Copyright (c) 2013, Johannes Meyer and contributors, Technische Universitat Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/system/ground_vehicle_model.h>
#include <hector_pose_estimation/pose_estimation.h>
#include <hector_pose_estimation/filter/set_filter.h>

#include <limits>

namespace hector_pose_estimation {

template class System_<GroundVehicleModel>;

GroundVehicleModel::GroundVehicleModel()
{
  gain_ = 1.0;
  base_height_ = 0.0;
  min_height_ = -std::numeric_limits<double>::quiet_NaN();
  max_height_ =  std::numeric_limits<double>::quiet_NaN();

  parameters().add("gain", gain_);
  parameters().add("base_height", base_height_);
  parameters().add("min_height", min_height_);
  parameters().add("max_height", max_height_);

  // derivative of the 3rd column of the rotation matrix
  dR3 <<  0.0, 1.0, 0.0,
         -1.0, 0.0, 0.0,
          0.0, 0.0, 0.0;
}

GroundVehicleModel::~GroundVehicleModel()
{
}

void GroundVehicleModel::getPrior(State &state)
{
  GenericQuaternionSystemModel::getPrior(state);
  if (state.position()) state.position()->vector().z() = base_height_;
}

void GroundVehicleModel::getDerivative(StateVector& x_dot, const State& state)
{
  // forward to GenericQuaternionSystemModel
  GenericQuaternionSystemModel::getDerivative(x_dot, state);

  const State::RotationMatrix &R = state.R();
  State::ConstVelocityType v(state.getVelocity());

  // Update the body z velocity towards 0
  if (state.velocity()) {
    // v_z_body = R.col(2).dot(v)
    state.velocity()->segment(x_dot) += -gain_ * R.col(2) * (R.col(2).dot(v));
  }
}

void GroundVehicleModel::getStateJacobian(SystemMatrix& A, const State& state, bool init)
{
  GenericQuaternionSystemModel::getStateJacobian(A, state, init);

  const State::RotationMatrix &R = state.R();
  State::ConstVelocityType v(state.getVelocity());

  if (state.velocity()) {
    state.velocity()->block(A) += -gain_ * R.col(2) * R.col(2).transpose();

    if (state.orientation()) {
      state.velocity()->block(A, *state.orientation()) += -gain_ * (dR3 * (R.col(2).dot(v)) + R.col(2) * (v.transpose() * dR3));
    }
  }
}

SystemStatus GroundVehicleModel::getStatusFlags(const State& state)
{
  SystemStatus flags = GenericQuaternionSystemModel::getStatusFlags(state);
  if (flags & STATE_VELOCITY_XY) {
    flags |= STATE_VELOCITY_Z;
    flags |= STATE_POSITION_Z;
  }
  return flags;
}

bool GroundVehicleModel::limitState(State& state)
{
  bool result = GenericQuaternionSystemModel::limitState(state);
  if (state.position()) {
    if (state.position()->vector().z() < min_height_) {
      state.position()->vector().z() = min_height_;
      result = false;
    }
    if (state.position()->vector().z() > max_height_) {
      state.position()->vector().z() = max_height_;
      result = false;
    }
  }
  return result;
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/system/generic_quaternion_system_model.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer and Martin Nowara, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/system/generic_quaternion_system_model.h>
#include <hector_pose_estimation/filter/set_filter.h>

#include <hector_pose_estimation/pose_estimation.h>

namespace hector_pose_estimation {

template class System_<GenericQuaternionSystemModel>;

GenericQuaternionSystemModel::GenericQuaternionSystemModel()
{
  angular_acceleration_stddev_ = 360.0 * M_PI/180.0;
  rate_stddev_ = 0.0;
  acceleration_stddev_ = 0.0;
  velocity_stddev_ = 0.0;

//  parameters().addAlias("gravity", gravity_);
  parameters().add("angular_acceleration_stddev", angular_acceleration_stddev_);
  parameters().add("rate_stddev", rate_stddev_);
  parameters().add("acceleration_stddev", acceleration_stddev_);
  parameters().add("velocity_stddev", velocity_stddev_);
}

GenericQuaternionSystemModel::~GenericQuaternionSystemModel()
{
}

bool GenericQuaternionSystemModel::init(PoseEstimation& estimator, System &system, State& state)
{
  gravity_ = estimator.parameters().get("gravity_magnitude");

  imu_ = estimator.getInputType<ImuInput>("imu");
  if (imu_ && state.orientation()) {
    gyro_ = estimator.getSystem_<Gyro>("gyro");
    if (!gyro_) {
      gyro_.reset(new Gyro("gyro"));
      estimator.addSystem(gyro_);
    }
  }
  if (imu_ && state.velocity()) {
    accelerometer_ = estimator.getSystem_<Accelerometer>("accelerometer");
    if (!accelerometer_) {
      accelerometer_.reset(new Accelerometer("accelerometer"));
      estimator.addSystem(accelerometer_);
    }

  }

//  // precalculate Q in order to just copy the parts that are active
//  // Note: Q might be too small as we add other substates later.
//  Q_ = State::Covariance::Zero(state.getCovarianceDimension(), state.getCovarianceDimension());
//  if (state.orientation()) {
//    if (!state.rate() && imu_ && gyro_) {
//      gyro_->getModel()->getRateNoise(state.orientation()->block(Q_), state, true);
//    }
//    state.orientation()->block(Q_) += pow(rate_stddev_, 2) * SymmetricMatrix3::Identity();
//  }
//  if (state.rate()) {
//    state.rate()->block(Q_) = pow(angular_acceleration_stddev_, 2) * SymmetricMatrix3::Identity();
//  }
//  if (state.position()) {
//    state.position()->block(Q_) = pow(velocity_stddev_, 2) * SymmetricMatrix3::Identity();
//  }
//  if (state.velocity()) {
//    if (!state.acceleration() && imu_ && accelerometer_) {
//      accelerometer_->getModel()->getAccelerationNoise(state.velocity()->block(Q_), state, true);
//    }
//    state.velocity()->block(Q_) += pow(acceleration_stddev_, 2) * SymmetricMatrix3::Identity();
//  }

  return true;
}

void GenericQuaternionSystemModel::getPrior(State &state) {
  if (state.orientation()) {
    state.orientation()->P()(X,X) = 1.0;
    state.orientation()->P()(Y,Y) = 1.0;
    state.orientation()->P()(Z,Z) = 0.0;
  }

  if (state.rate()) {
    state.rate()->P()(X,X) = pow(0.0 * M_PI/180.0, 2);
    state.rate()->P()(Y,Y) = pow(0.0 * M_PI/180.0, 2);
    state.rate()->P()(Z,Z) = pow(0.0 * M_PI/180.0, 2);
  }

  if (state.position()) {
    state.position()->P()(X,X) = 0.0;
    state.position()->P()(Y,Y) = 0.0;
    state.position()->P()(Z,Z) = 0.0;
  }

  if (state.velocity()) {
    state.velocity()->P()(X,X) = 0.0;
    state.velocity()->P()(Y,Y) = 0.0;
    state.velocity()->P()(Z,Z) = 0.0;
  }
}

bool GenericQuaternionSystemModel::prepareUpdate(State& state, double dt)
{
  if (state.rate()) {
    rate_nav_ = state.R() * state.getRate();
  }
  else if (rate_input_) {
    rate_nav_ = state.R() * rate_input_->getVector();
  }
  else if (imu_) {
    if (gyro_) {
      rate_nav_ = state.R() * gyro_->getModel()->getRate(imu_->getRate(), state);
    } else {
      rate_nav_ = state.R() * imu_->getRate();
    }
  } else {
    rate_nav_.setZero();
  }

  if (state.acceleration()) {
    acceleration_nav_ = state.R() * state.getAcceleration();
  }
  else if (force_input_) {
    acceleration_nav_ = state.R() * force_input_->getVector();
  }
  else if (imu_) {
    if (accelerometer_) {
      acceleration_nav_ = state.R() * accelerometer_->getModel()->getAcceleration(imu_->getAcceleration(), state);
    } else {
      acceleration_nav_ = state.R() * imu_->getAcceleration();
    }
  } else {
    acceleration_nav_.setZero();
  }

  ROS_DEBUG_STREAM_NAMED("system", "rate_nav = [" << rate_nav_.transpose() << "]");
  ROS_DEBUG_STREAM_NAMED("system", "acceleration_nav = [" << acceleration_nav_.transpose() << "]");
  return true;
}

void GenericQuaternionSystemModel::getDerivative(StateVector& x_dot, const State& state)
{
  x_dot.setZero();

  if (state.rate()) {
    if (torque_input_) {
      state.rate()->segment(x_dot) = torque_input_->getVector();
    }
  }

  if (state.orientation()) {
    state.orientation()->segment(x_dot).head(3) = rate_nav_;
    if (!(state.getSystemStatus() & STATE_YAW) || (state.getSystemStatus() & STATUS_ALIGNMENT)) {
      state.orientation()->segment(x_dot).z() = 0.0;
    }
  }

  if (state.velocity()) {
    if ((state.getSystemStatus() & STATE_VELOCITY_XY) && !(state.getSystemStatus() & STATUS_ALIGNMENT)) {
      state.velocity()->segment(x_dot)(X) = acceleration_nav_.x();
      state.velocity()->segment(x_dot)(Y) = acceleration_nav_.y();
    }
    if ((state.getSystemStatus() & STATE_VELOCITY_Z) && !(state.getSystemStatus() & STATUS_ALIGNMENT)) {
      state.velocity()->segment(x_dot)(Z) = acceleration_nav_.z();
      if (imu_) {
        state.velocity()->segment(x_dot)(Z) += gravity_;
      }
    }
  }

  if (state.position()) {
    State::ConstVelocityType v(state.getVelocity());
    if ((state.getSystemStatus() & STATE_POSITION_XY) && !(state.getSystemStatus() & STATUS_ALIGNMENT)) {
      state.position()->segment(x_dot)(X) = v.x();
      state.position()->segment(x_dot)(Y) = v.y();
    }
    if ((state.getSystemStatus() & STATE_POSITION_Z) && !(state.getSystemStatus() & STATUS_ALIGNMENT)) {
      state.position()->segment(x_dot)(Z) = v.z();
    }
  }
}

void GenericQuaternionSystemModel::getSystemNoise(NoiseVariance& Q, const State& state, bool init)
{
//  if (init) Q.setZero();
//  Q.topLeftCorner(Q_.rows(), Q_.cols()) = Q_;

//  if (state.orientation()) {
//    if (!(state.getSystemStatus() & STATE_YAW) || (state.getSystemStatus() & STATUS_ALIGNMENT)) {
//      state.orientation()->block(Q)(Z,Z) = 0.0;
//    }
//  }

//  if (state.velocity()) {
//    if (!(state.getSystemStatus() & STATE_VELOCITY_XY) || (state.getSystemStatus() & STATUS_ALIGNMENT)) {
//      state.velocity()->block(Q)(X,X) = 0.0;
//      state.velocity()->block(Q)(Y,Y) = 0.0;
//    }
//    if (!(state.getSystemStatus() & STATE_VELOCITY_Z) || (state.getSystemStatus() & STATUS_ALIGNMENT)) {
//      state.velocity()->block(Q)(Z,Z) = 0.0;
//    }
//  }

//  if (state.position()) {
//    if (!(state.getSystemStatus() & STATE_POSITION_XY) || (state.getSystemStatus() & STATUS_ALIGNMENT)) {
//      state.position()->block(Q)(X,X) = 0.0;
//      state.position()->block(Q)(Y,Y) = 0.0;
//    }
//    if (!(state.getSystemStatus() & STATE_POSITION_Z) || (state.getSystemStatus() & STATUS_ALIGNMENT)) {
//      state.position()->block(Q)(Z,Z) = 0.0;
//    }
//  }

  if (!init) return;

  Q.setZero();
  if (state.orientation()) {
    if (!state.rate() && imu_ && gyro_) {
      gyro_->getModel()->getRateNoise(state.orientation()->block(Q), state, init);
    }
    state.orientation()->block(Q) += pow(rate_stddev_, 2) * SymmetricMatrix3::Identity();
  }
  if (state.rate()) {
    state.rate()->block(Q) = pow(angular_acceleration_stddev_, 2) * SymmetricMatrix3::Identity();
  }
  if (state.position()) {
    state.position()->block(Q) = pow(velocity_stddev_, 2) * SymmetricMatrix3::Identity();
  }
  if (state.velocity()) {
    if (!state.acceleration() && imu_ && accelerometer_) {
      accelerometer_->getModel()->getAccelerationNoise(state.velocity()->block(Q), state, init);
    }
    state.velocity()->block(Q) += pow(acceleration_stddev_, 2) * SymmetricMatrix3::Identity();
  }
}

void GenericQuaternionSystemModel::getStateJacobian(SystemMatrix& A, const State& state, bool)
{
  const State::RotationMatrix &R = state.R();
  A.setZero();

  if (state.orientation()) {
    if (state.rate()) {
      state.orientation()->block(A, *state.rate()) = R;
    } else if (imu_ && gyro_) {
      GyroModel::SystemMatrixBlock A_orientation = state.orientation()->rows(A);
      gyro_->getModel()->getRateJacobian(A_orientation, state);
      state.orientation()->rows(A) = R * A_orientation;
    }

    state.orientation()->block(A) += SkewSymmetricMatrix(-rate_nav_);

    if (!(state.getSystemStatus() & STATE_YAW) || (state.getSystemStatus() & STATUS_ALIGNMENT)) {
      state.orientation()->rows(A).row(2).setZero();
      state.orientation()->block(A).col(2).setZero();
    }
  }

  if (state.velocity()) {
    if (state.acceleration()) {
      state.velocity()->block(A, *state.acceleration()) = R;
    } else if (imu_ && accelerometer_) {
      AccelerometerModel::SystemMatrixBlock A_velocity = state.velocity()->rows(A);
      accelerometer_->getModel()->getAccelerationJacobian(A_velocity, state);
      state.velocity()->rows(A) = R * A_velocity;
    }

    if (state.orientation()) {
      state.velocity()->block(A, *state.orientation()) += SkewSymmetricMatrix(-acceleration_nav_);
    }

    if (!(state.getSystemStatus() & STATE_VELOCITY_XY) || (state.getSystemStatus() & STATUS_ALIGNMENT)) {
      state.velocity()->rows(A).topRows(2).setZero();
    }

    if (!(state.getSystemStatus() & STATE_VELOCITY_Z) || (state.getSystemStatus() & STATUS_ALIGNMENT)) {
      state.velocity()->rows(A).row(2).setZero();
    }
  }

  if (state.position() && state.velocity()) {
    state.position()->block(A, *state.velocity()).setIdentity();

    if ((state.getSystemStatus() & STATE_POSITION_XY) && !(state.getSystemStatus() & STATUS_ALIGNMENT)) {
      state.position()->block(A, *state.velocity())(X,X) = 1.0;
      state.position()->block(A, *state.velocity())(Y,Y) = 1.0;
    }
    if ((state.getSystemStatus() & STATE_POSITION_Z) && !(state.getSystemStatus() & STATUS_ALIGNMENT)) {
      state.position()->block(A, *state.velocity())(Z,Z) = 1.0;
    }
  }

}

SystemStatus GenericQuaternionSystemModel::getStatusFlags(const State& state)
{
  SystemStatus flags = state.getMeasurementStatus();
  if (flags & STATE_POSITION_XY) flags |= STATE_VELOCITY_XY;
  if (flags & STATE_POSITION_Z)  flags |= STATE_VELOCITY_Z;
  if (imu_) {
    if (flags & STATE_VELOCITY_XY)      flags |= STATE_ROLLPITCH;
    if (flags & STATE_ROLLPITCH)        flags |= STATE_RATE_XY;
    if (flags & STATE_PSEUDO_ROLLPITCH) flags |= STATE_PSEUDO_RATE_XY;
    if (flags & STATE_YAW)              flags |= STATE_RATE_Z;
    if (flags & STATE_PSEUDO_YAW)       flags |= STATE_PSEUDO_RATE_Z;
  }
  return flags & STATE_MASK;
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/system/imu_model.cpp
//=================================================================================================
// Copyright (c) 2013, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/system/imu_model.h>
#include <hector_pose_estimation/pose_estimation.h>
#include <hector_pose_estimation/filter/set_filter.h>

namespace hector_pose_estimation {

template class System_<GyroModel>;
template class System_<AccelerometerModel>;

static const Matrix3 MinusIdentity = -Matrix3::Identity();

GyroModel::GyroModel()
{
  rate_stddev_ = 1.0 * M_PI/180.0;
  rate_drift_ = 1.0e-1 * M_PI/180.0;
  parameters().add("stddev", rate_stddev_);
  parameters().add("drift", rate_drift_);
}

GyroModel::~GyroModel()
{}

bool GyroModel::init(PoseEstimation& estimator, System &system, State& state)
{
  bias_ = state.addSubState<3,3>(this, system.getName() + "_bias");
  return (bias_.get() != 0);
}

void GyroModel::getPrior(State &state)
{
  bias_->block(state.P()) = 3600./2. * pow(rate_drift_, 2) * SymmetricMatrix3::Identity();
}

void GyroModel::getSystemNoise(NoiseVariance& Q, const State& state, bool init)
{
  if (!init) return;
  bias_->block(Q)(X,X) = bias_->block(Q)(Y,Y) = pow(rate_drift_, 2);
  bias_->block(Q)(Z,Z) = pow(rate_drift_, 2);
}

ColumnVector3 GyroModel::getRate(const ImuInput::RateType& imu_rate, const State& state) const
{
  return imu_rate - bias_->getVector();
}

void GyroModel::getRateJacobian(SystemMatrixBlock& C, const State& state, bool init)
{
  if (!init) return;
  bias_->cols(C) = MinusIdentity;
}

void GyroModel::getRateNoise(CovarianceBlock Q, const State &, bool init)
{
  if (!init) return;
  Q(X,X) = Q(Y,Y) = Q(Z,Z) = pow(rate_stddev_, 2);
}

//void GyroModel::getDerivative(StateVector &x_dot, const State &state)
//{
//  x_dot.setZero();
//  if (state.orientation() && !state.rate()) {
//    state.orientation()->segment(x_dot).head(3) = state.R() * bias_->vector();
//  }
//}

//void GyroModel::getStateJacobian(SystemMatrix& A, const State& state)
//{
//  A.setZero();
//  if (state.orientation() && !state.rate()) {
//    state.orientation()->block(A, *bias_) = state.R();
//  }
//}

AccelerometerModel::AccelerometerModel()
{
  acceleration_stddev_ = 1.0e-2;
  acceleration_drift_ = 1.0e-2;
  parameters().add("stddev", acceleration_stddev_);
  parameters().add("drift", acceleration_drift_);
}

AccelerometerModel::~AccelerometerModel()
{}

bool AccelerometerModel::init(PoseEstimation& estimator, System &system, State& state)
{
  bias_ = state.addSubState<3,3>(this, system.getName() + "_bias");
  return (bias_.get() != 0);
}

void AccelerometerModel::getPrior(State &state)
{
  bias_->block(state.P()) = 3600./2. * pow(acceleration_drift_, 2) * SymmetricMatrix3::Identity();
}

void AccelerometerModel::getSystemNoise(NoiseVariance& Q, const State&, bool init)
{
  if (!init) return;
  bias_->block(Q)(X,X) = bias_->block(Q)(Y,Y) = pow(acceleration_drift_, 2);
  bias_->block(Q)(Z,Z) = pow(acceleration_drift_, 2);
}

ColumnVector3 AccelerometerModel::getAcceleration(const ImuInput::AccelerationType& imu_acceleration, const State& state) const
{
  return imu_acceleration - bias_->getVector();
}

void AccelerometerModel::getAccelerationJacobian(SystemMatrixBlock& C, const State&, bool init)
{
  if (!init) return;
  bias_->cols(C) = MinusIdentity;
}

void AccelerometerModel::getAccelerationNoise(CovarianceBlock Q, const State &, bool init)
{
  if (!init) return;
  Q(X,X) = Q(Y,Y) = Q(Z,Z) = pow(acceleration_stddev_, 2);
}

//void AccelerometerModel::getDerivative(StateVector &x_dot, const State &state)
//{
//  x_dot.setZero();
//  if (state.velocity() && !state.acceleration()) {
//    if (state.getSystemStatus() & STATE_VELOCITY_XY) {
//      state.velocity()->segment(x_dot)(X) = bias_nav_.x();
//      state.velocity()->segment(x_dot)(Y) = bias_nav_.y();
//    }
//    if (state.getSystemStatus() & STATE_VELOCITY_Z) {
//      state.velocity()->segment(x_dot)(Z) = bias_nav_.z();
//    }
//  }
//}

//void AccelerometerModel::getStateJacobian(SystemMatrixBlock& A, const State& state)
//{
//  A.setZero();
//  if (state.velocity() && !state.acceleration()) {
//    const State::RotationMatrix &R = state.R();

//    if (state.getSystemStatus() & STATE_VELOCITY_XY) {
//      state.velocity()->block(A, *bias_).row(X) = R.row(X);
//      state.velocity()->block(A, *bias_).row(Y) = R.row(Y);
//    }
//    if (state.getSystemStatus() & STATE_VELOCITY_Z) {
//      state.velocity()->block(A, *bias_).row(Z) = R.row(Z);
//    }

//    if (state.getSystemStatus() & STATE_VELOCITY_XY) {
//      state.velocity()->block(A, *state.orientation())(X,X) = 0.0;
//      state.velocity()->block(A, *state.orientation())(X,Y) =  bias_nav_.z();
//      state.velocity()->block(A, *state.orientation())(X,Z) = -bias_nav_.y();

//      state.velocity()->block(A, *state.orientation())(Y,X) = -bias_nav_.z();
//      state.velocity()->block(A, *state.orientation())(Y,Y) = 0.0;
//      state.velocity()->block(A, *state.orientation())(Y,Z) =  bias_nav_.x();
//    }

//    if (state.getSystemStatus() & STATE_VELOCITY_Z) {
//      state.velocity()->block(A, *state.orientation())(Z,X) =  bias_nav_.y();
//      state.velocity()->block(A, *state.orientation())(Z,Y) = -bias_nav_.x();
//      state.velocity()->block(A, *state.orientation())(Z,Z) = 0.0;
//    }
//  }
//}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/pose_estimation.cpp
//=================================================================================================

// Copyright (c) 2011, Johannes Meyer and Martin Nowara, TU Darmstadt

// All rights reserved.



// Redistribution and use in source and binary forms, with or without

// modification, are permitted provided that the following conditions are met:

//     * Redistributions of source code must retain the above copyright

//       notice, this list of conditions and the following disclaimer.

//     * Redistributions in binary form must reproduce the above copyright

//       notice, this list of conditions and the following disclaimer in the

//       documentation and/or other materials provided with the distribution.

//     * Neither the name of the Flight Systems and Automatic Control group,

//       TU Darmstadt, nor the names of its contributors may be used to

//       endorse or promote products derived from this software without

//       specific prior written permission.



// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND

// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED

// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE

// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY

// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES

// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;

// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND

// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT

// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS

// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//=================================================================================================



#include <hector_pose_estimation/pose_estimation.h>

#include <hector_pose_estimation/filter/ekf.h>

#include <hector_pose_estimation/global_reference.h>



#include <hector_pose_estimation/system/imu_input.h>

#include <hector_pose_estimation/system/imu_model.h>



#include <boost/weak_ptr.hpp>



namespace hector_pose_estimation {



namespace {

  static PoseEstimation *the_instance = 0;

}



PoseEstimation::PoseEstimation(const SystemPtr& system, const StatePtr& state)

  : state_(state ? state : StatePtr(new OrientationPositionVelocityState))

  , rate_update_(new Rate("rate"))

  , gravity_update_(new Gravity ("gravity"))

  , zerorate_update_(new ZeroRate("zerorate"))

{

  if (!the_instance) the_instance = this;

  if (system) addSystem(system);



  world_frame_ = "/world";

  nav_frame_ = "nav";

  base_frame_ = "base_link";

  stabilized_frame_ = "base_stabilized";

  footprint_frame_ = "base_footprint";

  // position_frame_ = "base_position";

  alignment_time_ = 0.0;

  gravity_ = -9.8065;



  parameters().add("world_frame", world_frame_);

  parameters().add("nav_frame", nav_frame_);

  parameters().add("base_frame", base_frame_);

  parameters().add("stabilized_frame", stabilized_frame_);

  parameters().add("footprint_frame", footprint_frame_);

  parameters().add("position_frame", position_frame_);

  parameters().add(globalReference()->parameters());

  parameters().add("alignment_time", alignment_time_);

  parameters().add("gravity_magnitude", gravity_);



  // add default measurements

  addMeasurement(rate_update_);

  addMeasurement(gravity_update_);

  addMeasurement(zerorate_update_);

}



PoseEstimation::~PoseEstimation()

{

  cleanup();

}



PoseEstimation *PoseEstimation::Instance() {

  if (!the_instance) the_instance = new PoseEstimation();

  return the_instance;

}



bool PoseEstimation::init()

{

#ifdef EIGEN_RUNTIME_NO_MALLOC

  Eigen::internal::set_is_malloc_allowed(true);

#endif



  // initialize global reference

  globalReference()->reset();



  // check if system is initialized

  if (systems_.empty()) return false;



  // create new filter

  filter_.reset(new filter::EKF(*state_));



  // initialize systems (new systems could be added during initialization!)

  for(Systems::iterator it = systems_.begin(); it != systems_.end(); ++it)

    if (!(*it)->init(*this, state())) return false;



  // initialize measurements (new systems could be added during initialization!)

  for(Measurements::iterator it = measurements_.begin(); it != measurements_.end(); ++it)

    if (!(*it)->init(*this, state())) return false;



  // initialize filter

  filter_->init(*this);



  // call setFilter for each system and each measurement

  for(Systems::iterator it = systems_.begin(); it != systems_.end(); ++it)

    (*it)->setFilter(filter_.get());

  for(Measurements::iterator it = measurements_.begin(); it != measurements_.end(); ++it)

    (*it)->setFilter(filter_.get());



  // reset (or initialize) filter and measurements

  reset();



  return true;

}



void PoseEstimation::cleanup()

{

  // cleanup system

  for(Systems::iterator it = systems_.begin(); it != systems_.end(); ++it) (*it)->cleanup();



  // cleanup measurements

  for(Measurements::iterator it = measurements_.begin(); it != measurements_.end(); ++it) (*it)->cleanup();



  // delete filter instance

  if (filter_) filter_.reset();

}



void PoseEstimation::reset()

{

  // check if system is initialized

  if (systems_.empty()) return;



  // set initial status

  if (filter_) filter_->reset();



  // restart alignment

  alignment_start_ = ros::Time();

  if (alignment_time_ > 0) {

    state().setSystemStatus(STATUS_ALIGNMENT);

  }



  // reset systems and measurements

  for(Systems::iterator it = systems_.begin(); it != systems_.end(); ++it) {

    (*it)->reset(state());

    (*it)->getPrior(state());

  }



  for(Measurements::iterator it = measurements_.begin(); it != measurements_.end(); ++it) {

    (*it)->reset(state());

  }



  updated();

}



void PoseEstimation::update(ros::Time new_timestamp)

{

  // check if system is initialized

  if (systems_.empty()) return;



  ros::Duration dt;

  if (!getTimestamp().isZero()) {

    if (new_timestamp.isZero()) new_timestamp = ros::Time::now();

    dt = new_timestamp - getTimestamp();

  }

  setTimestamp(new_timestamp);



  // do the update step

  update(dt.toSec());

}



void PoseEstimation::update(double dt)

{

  // check dt

  if (dt < -1.0)

    reset();

  else if (dt < 0.0)

    return;

  else if (dt > 1.0)

    dt = 1.0;



  // check if system and filter is initialized

  if (systems_.empty() || !filter_) return;



  // filter rate measurement first

  boost::shared_ptr<ImuInput> imu = getInputType<ImuInput>("imu");

  if (imu) {

    // Should the biases already be integrated here?

    // Note: The state set here only has an effect if the state vector does not have a rate/acceleration component.

    state().setRate(imu->getRate());

    state().setAcceleration(imu->getAcceleration() + state().R().row(2).transpose() * gravity_);



    if (state().rate() && rate_update_) {

      rate_update_->update(Rate::Update(imu->getRate()));

    }

  }



  // time update step

  filter_->predict(systems_, dt);



  // pseudo measurement updates (if required)

  if (imu && !(getSystemStatus() & STATE_ROLLPITCH)) {

    gravity_update_->enable();

    gravity_update_->update(Gravity::Update(imu->getAcceleration()));

  } else {

    gravity_update_->disable();

  }

  if (!(getSystemStatus() & STATE_RATE_Z)) {

    zerorate_update_->enable();

    zerorate_update_->update(ZeroRate::Update());

  } else {

    zerorate_update_->disable();

  }



  // measurement updates

  filter_->correct(measurements_);



  // updated hook

  updated();



  // set measurement status and increase timers

  SystemStatus measurement_status = 0;

  for(Measurements::iterator it = measurements_.begin(); it != measurements_.end(); it++) {

    const MeasurementPtr& measurement = *it;

    measurement_status |= measurement->getStatusFlags();

    measurement->increase_timer(dt);

  }

  setMeasurementStatus(measurement_status);



  // set system status

  SystemStatus system_status = 0;

  for(Systems::iterator it = systems_.begin(); it != systems_.end(); it++) {

    const SystemPtr& system = *it;

    system_status |= system->getStatusFlags();

  }

  updateSystemStatus(system_status, STATE_MASK | STATE_PSEUDO_MASK);



  // check for invalid state

  if (!state().valid()) {

    ROS_FATAL("Invalid state, resetting...");

    reset();

    return;

  }



  // switch overall system status

  if (inSystemStatus(STATUS_ALIGNMENT)) {

    if (alignment_start_.isZero()) alignment_start_ = getTimestamp();

    if ((getTimestamp() - alignment_start_).toSec() >= alignment_time_) {

      updateSystemStatus(STATUS_DEGRADED, STATUS_ALIGNMENT);

    }

  } else if (inSystemStatus(STATE_ROLLPITCH | STATE_YAW | STATE_POSITION_XY | STATE_POSITION_Z)) {

    updateSystemStatus(STATUS_READY, STATUS_DEGRADED);

  } else {

    updateSystemStatus(STATUS_DEGRADED, STATUS_READY);

  }





#ifdef EIGEN_RUNTIME_NO_MALLOC

  // No memory allocations allowed after the first update!

  Eigen::internal::set_is_malloc_allowed(false);

#endif

}



void PoseEstimation::updated() {

  for(Systems::iterator it = systems_.begin(); it != systems_.end(); ++it) {

    (*it)->limitState(state());

  }

}



const SystemPtr& PoseEstimation::addSystem(const SystemPtr& system, const std::string& name) {

  if (!name.empty() && system->getName().empty()) system->setName(name);

  parameters().add(system->getName(), system->parameters());

  return systems_.add(system, system->getName());

}



InputPtr PoseEstimation::addInput(const InputPtr& input, const std::string& name)

{

  if (!name.empty()) input->setName(name);

  return inputs_.add(input, input->getName());

}



InputPtr PoseEstimation::setInput(const Input& value, std::string name)

{

  if (name.empty()) name = value.getName();

  InputPtr input = inputs_.get(name);

  if (!input) {

    ROS_WARN("Set input \"%s\", but this input is not registered by any system model.", name.c_str());

    return InputPtr();

  }



  *input = value;

  return input;

}



const MeasurementPtr& PoseEstimation::addMeasurement(const MeasurementPtr& measurement, const std::string& name) {

  if (!name.empty()) measurement->setName(name);

  parameters().add(measurement->getName(), measurement->parameters());

  return measurements_.add(measurement, measurement->getName());

}



const State::Vector& PoseEstimation::getStateVector() {

//  if (state_is_dirty_) {

//    state_ = filter_->PostGet()->ExpectedValueGet();

//    state_is_dirty_ = false;

//  }

  return state().getVector();

}



const State::Covariance& PoseEstimation::getCovariance() {

//  if (covariance_is_dirty_) {

//    covariance_ = filter_->PostGet()->CovarianceGet();

//    covariance_is_dirty_ = false;

//  }

  return state().getCovariance();

}



SystemStatus PoseEstimation::getSystemStatus() const {

  return state().getSystemStatus();

}



SystemStatus PoseEstimation::getMeasurementStatus() const {

  return state().getMeasurementStatus();

}



bool PoseEstimation::inSystemStatus(SystemStatus test_status) const {

  return state().inSystemStatus(test_status);

}



bool PoseEstimation::setSystemStatus(SystemStatus new_status) {

  return state().setSystemStatus(new_status);

}



bool PoseEstimation::setMeasurementStatus(SystemStatus new_measurement_status) {

  return state().setMeasurementStatus(new_measurement_status);

}



bool PoseEstimation::updateSystemStatus(SystemStatus set, SystemStatus clear) {

  return state().updateSystemStatus(set, clear);

}



bool PoseEstimation::updateMeasurementStatus(SystemStatus set, SystemStatus clear) {

  return state().updateMeasurementStatus(set, clear);

}



const ros::Time& PoseEstimation::getTimestamp() const {

  return state().getTimestamp();

}



void PoseEstimation::setTimestamp(const ros::Time& timestamp) {

  state().setTimestamp(timestamp);

}



void PoseEstimation::getHeader(std_msgs::Header& header) {

  header.stamp = getTimestamp();

  header.frame_id = nav_frame_;

}



void PoseEstimation::getState(nav_msgs::Odometry& msg, bool with_covariances) {

  getHeader(msg.header);

  getPose(msg.pose.pose);

  getVelocity(msg.twist.twist.linear);

  getRate(msg.twist.twist.angular);

  msg.child_frame_id = base_frame_;



  // rotate body vectors to nav frame

  ColumnVector3 rate_nav = state().R() * ColumnVector3(msg.twist.twist.angular.x, msg.twist.twist.angular.y, msg.twist.twist.angular.z);

  msg.twist.twist.angular.x = rate_nav.x();

  msg.twist.twist.angular.y = rate_nav.y();

  msg.twist.twist.angular.z = rate_nav.z();



  // fill covariances

  if (with_covariances) {

    Eigen::Map< Eigen::Matrix<geometry_msgs::PoseWithCovariance::_covariance_type::value_type,6,6> >  pose_covariance_msg(msg.pose.covariance.data());

    Eigen::Map< Eigen::Matrix<geometry_msgs::TwistWithCovariance::_covariance_type::value_type,6,6> > twist_covariance_msg(msg.twist.covariance.data());



    // position covariance

    if (state().position()) {

      pose_covariance_msg.block<3,3>(0,0) = state().position()->getCovariance();

    }



    // rotation covariance (world-fixed)

    if (state().orientation()) {

      pose_covariance_msg.block<3,3>(3,3) = state().orientation()->getCovariance();

    }



    // position/orientation cross variance

    if (state().position() && state().orientation()) {

      pose_covariance_msg.block<3,3>(0,3) = state().orientation()->getCrossVariance(*state().position());

      pose_covariance_msg.block<3,3>(3,0) = pose_covariance_msg.block<3,3>(0,3).transpose();

    }



    // velocity covariance

    if (state().velocity()) {

      twist_covariance_msg.block<3,3>(0,0) = state().velocity()->getCovariance();

    }



    // angular rate covariance

    if (state().rate()) {

      twist_covariance_msg.block<3,3>(3,3) = state().rate()->getCovariance();

    }



    // cross velocity/angular_rate variance

    if (state().velocity() && state().rate()) {

      pose_covariance_msg.block<3,3>(0,3) = state().velocity()->getCrossVariance(*state().rate());

      pose_covariance_msg.block<3,3>(3,0) = pose_covariance_msg.block<3,3>(0,3).transpose();

    }

  }

}



void PoseEstimation::getPose(tf::Pose& pose) {

  tf::Quaternion quaternion;

  getPosition(pose.getOrigin());

  getOrientation(quaternion);

  pose.setRotation(quaternion);

}



void PoseEstimation::getPose(tf::Stamped<tf::Pose>& pose) {

  getPose(static_cast<tf::Pose &>(pose));

  pose.stamp_ = getTimestamp();

  pose.frame_id_ = nav_frame_;

}



void PoseEstimation::getPose(geometry_msgs::Pose& pose) {

  getPosition(pose.position);

  getOrientation(pose.orientation);

}



void PoseEstimation::getPose(geometry_msgs::PoseStamped& pose) {

  getHeader(pose.header);

  getPose(pose.pose);

}



void PoseEstimation::getPosition(tf::Point& point) {

  State::ConstPositionType position(state().getPosition());

  point = tf::Point(position.x(), position.y(), position.z());

}



void PoseEstimation::getPosition(tf::Stamped<tf::Point>& point) {

  getPosition(static_cast<tf::Point &>(point));

  point.stamp_ = getTimestamp();

  point.frame_id_ = nav_frame_;

}



void PoseEstimation::getPosition(geometry_msgs::Point& point) {

  State::ConstPositionType position(state().getPosition());

  point.x = position.x();

  point.y = position.y();

  point.z = position.z();

}



void PoseEstimation::getPosition(geometry_msgs::PointStamped& point) {

  getHeader(point.header);

  getPosition(point.point);

}



void PoseEstimation::getGlobal(double &latitude, double &longitude, double &altitude) {

  State::ConstPositionType position(state().getPosition());

  double north =  position.x() * globalReference()->heading().cos - position.y() * globalReference()->heading().sin;

  double east  = -position.x() * globalReference()->heading().sin - position.y() * globalReference()->heading().cos;

  latitude  = globalReference()->position().latitude  + north / globalReference()->radius().north;

  longitude = globalReference()->position().longitude + east  / globalReference()->radius().east;

  altitude  = globalReference()->position().altitude  + position.z();

}



void PoseEstimation::getGlobalPosition(double &latitude, double &longitude, double &altitude) {

  getGlobal(latitude, longitude, altitude);

}



void PoseEstimation::getGlobal(geographic_msgs::GeoPoint& global)

{

  getGlobalPosition(global.latitude, global.longitude, global.altitude);

  global.latitude  *= 180.0/M_PI;

  global.longitude *= 180.0/M_PI;

}



void PoseEstimation::getGlobal(sensor_msgs::NavSatFix& global)

{

  getHeader(global.header);

  global.header.frame_id = world_frame_;



  if ((getSystemStatus() & STATE_POSITION_XY) && globalReference()->hasPosition()) {

    global.status.status = sensor_msgs::NavSatStatus::STATUS_FIX;

  } else {

    global.status.status = sensor_msgs::NavSatStatus::STATUS_NO_FIX;

  }



  getGlobalPosition(global.latitude, global.longitude, global.altitude);

  global.latitude  *= 180.0/M_PI;

  global.longitude *= 180.0/M_PI;



  if (getSystemStatus() & STATE_POSITION_XY) {

    global.status.status = sensor_msgs::NavSatStatus::STATUS_FIX;

  } else {

    global.status.status = sensor_msgs::NavSatStatus::STATUS_NO_FIX;

  }

}



void PoseEstimation::getGlobalPosition(sensor_msgs::NavSatFix& global)

{

  getGlobal(global);

}



void PoseEstimation::getGlobal(geographic_msgs::GeoPoint& position, geometry_msgs::Quaternion& quaternion)

{

  getGlobal(position);

  Quaternion global_orientation = globalReference()->heading().quaternion() * Quaternion(state().getOrientation());

  quaternion.w = global_orientation.w();

  quaternion.x = global_orientation.x();

  quaternion.y = global_orientation.y();

  quaternion.z = global_orientation.z();

}



void PoseEstimation::getGlobal(geographic_msgs::GeoPose& global)

{

  getGlobal(global.position, global.orientation);

}



void PoseEstimation::getOrientation(tf::Quaternion& quaternion) {

  Quaternion orientation(state().getOrientation());

  quaternion = tf::Quaternion(orientation.x(), orientation.y(), orientation.z(), orientation.w());

}



void PoseEstimation::getOrientation(tf::Stamped<tf::Quaternion>& quaternion) {

  getOrientation(static_cast<tf::Quaternion &>(quaternion));

  quaternion.stamp_ = getTimestamp();

  quaternion.frame_id_ = nav_frame_;

}



void PoseEstimation::getOrientation(geometry_msgs::Quaternion& quaternion) {

  Quaternion orientation(state().getOrientation());

  quaternion.w = orientation.w();

  quaternion.x = orientation.x();

  quaternion.y = orientation.y();

  quaternion.z = orientation.z();

}



void PoseEstimation::getOrientation(geometry_msgs::QuaternionStamped& quaternion) {

  getHeader(quaternion.header);

  getOrientation(quaternion.quaternion);

}



void PoseEstimation::getOrientation(double &yaw, double &pitch, double &roll) {

  tf::Quaternion quaternion;

  getOrientation(quaternion);

#ifdef TF_MATRIX3x3_H

  tf::Matrix3x3(quaternion).getRPY(roll, pitch, yaw);

#else

  btMatrix3x3(quaternion).getRPY(roll, pitch, yaw);

#endif

}



void PoseEstimation::getImuWithBiases(geometry_msgs::Vector3& linear_acceleration, geometry_msgs::Vector3& angular_velocity) {

  boost::shared_ptr<const ImuInput>  input     = boost::dynamic_pointer_cast<const ImuInput>(getInput("imu"));

  boost::shared_ptr<const Accelerometer> accel = boost::dynamic_pointer_cast<const Accelerometer>(getSystem("accelerometer"));



  if (input) {

    linear_acceleration.x = input->getAcceleration().x();

    linear_acceleration.y = input->getAcceleration().y();

    linear_acceleration.z = input->getAcceleration().z();

  } else {

    linear_acceleration.x = 0.0;

    linear_acceleration.y = 0.0;

    linear_acceleration.z = 0.0;

  }



  if (accel) {

    linear_acceleration.x -= accel->getModel()->getError().x();

    linear_acceleration.y -= accel->getModel()->getError().y();

    linear_acceleration.z -= accel->getModel()->getError().z();

  }



  getRate(angular_velocity);

}



void PoseEstimation::getVelocity(tf::Vector3& vector) {

  State::ConstVelocityType velocity(state().getVelocity());

  vector = tf::Vector3(velocity.x(), velocity.y(), velocity.z());

}



void PoseEstimation::getVelocity(tf::Stamped<tf::Vector3>& vector) {

  getVelocity(static_cast<tf::Vector3 &>(vector));

  vector.stamp_ = getTimestamp();

  vector.frame_id_ = nav_frame_;

}



void PoseEstimation::getVelocity(geometry_msgs::Vector3& vector) {

  State::ConstVelocityType velocity(state().getVelocity());

  vector.x = velocity.x();

  vector.y = velocity.y();

  vector.z = velocity.z();

}



void PoseEstimation::getVelocity(geometry_msgs::Vector3Stamped& vector) {

  getHeader(vector.header);

  getVelocity(vector.vector);

}



void PoseEstimation::getRate(tf::Vector3& vector) {

  geometry_msgs::Vector3 rate;

  getRate(rate);

  vector = tf::Vector3(rate.x, rate.y, rate.z);

}



void PoseEstimation::getRate(tf::Stamped<tf::Vector3>& vector) {

  getRate(static_cast<tf::Vector3 &>(vector));

  vector.stamp_ = getTimestamp();

  vector.frame_id_ = base_frame_;

}



void PoseEstimation::getRate(geometry_msgs::Vector3& vector) {

  if (state().rate()) {

    State::ConstRateType rate(state().getRate());

    vector.x    = rate.x();

    vector.y    = rate.y();

    vector.z    = rate.z();



  } else {

    boost::shared_ptr<const ImuInput> input = boost::dynamic_pointer_cast<const ImuInput>(getInput("imu"));

    boost::shared_ptr<const Gyro> gyro      = boost::dynamic_pointer_cast<const Gyro>(getSystem("gyro"));



    if (input) {

      vector.x = input->getRate().x();

      vector.y = input->getRate().y();

      vector.z = input->getRate().z();

    } else {

      vector.x = 0.0;

      vector.y = 0.0;

      vector.z = 0.0;

    }



    if (gyro) {

      vector.x -= gyro->getModel()->getError().x();

      vector.y -= gyro->getModel()->getError().y();

      vector.z -= gyro->getModel()->getError().z();

    }

  }

}



void PoseEstimation::getRate(geometry_msgs::Vector3Stamped& vector) {

  getHeader(vector.header);

  getRate(vector.vector);

  vector.header.frame_id = base_frame_;

}



void PoseEstimation::getBias(geometry_msgs::Vector3& angular_velocity, geometry_msgs::Vector3& linear_acceleration) {

  boost::shared_ptr<const Accelerometer> accel = boost::dynamic_pointer_cast<const Accelerometer>(getSystem("accelerometer"));

  boost::shared_ptr<const Gyro> gyro           = boost::dynamic_pointer_cast<const Gyro>(getSystem("gyro"));



  if (gyro) {

    angular_velocity.x = gyro->getModel()->getError().x();

    angular_velocity.y = gyro->getModel()->getError().y();

    angular_velocity.z = gyro->getModel()->getError().z();

  } else {

    angular_velocity.x = 0.0;

    angular_velocity.y = 0.0;

    angular_velocity.z = 0.0;

  }



  if (accel) {

    linear_acceleration.x = accel->getModel()->getError().x();

    linear_acceleration.y = accel->getModel()->getError().y();

    linear_acceleration.z = accel->getModel()->getError().z();

  } else {

    linear_acceleration.x = 0.0;

    linear_acceleration.y = 0.0;

    linear_acceleration.z = 0.0;

  }

}



void PoseEstimation::getBias(geometry_msgs::Vector3Stamped& angular_velocity, geometry_msgs::Vector3Stamped& linear_acceleration) {

  getBias(angular_velocity.vector, linear_acceleration.vector);

  angular_velocity.header.stamp = getTimestamp();

  angular_velocity.header.frame_id = base_frame_;

  linear_acceleration.header.stamp = getTimestamp();

  linear_acceleration.header.frame_id = base_frame_;

}



void PoseEstimation::getTransforms(std::vector<tf::StampedTransform>& transforms) {

  tf::Quaternion orientation;

  tf::Point position;

  getOrientation(orientation);

  getPosition(position);



  tf::Transform transform(orientation, position);

  double y,p,r;

  transform.getBasis().getEulerYPR(y,p,r);



  std::string parent_frame = nav_frame_;



  if(!position_frame_.empty()) {

    tf::Transform position_transform;

    position_transform.getBasis().setIdentity();

    position_transform.setOrigin(tf::Point(position.x(), position.y(), position.z()));

    transforms.push_back(tf::StampedTransform(position_transform, getTimestamp(), parent_frame, position_frame_ ));

  }



  if (!footprint_frame_.empty()) {

    tf::Transform footprint_transform;

    footprint_transform.getBasis().setEulerYPR(y, 0.0, 0.0);

    footprint_transform.setOrigin(tf::Point(position.x(), position.y(), 0.0));

    transforms.push_back(tf::StampedTransform(footprint_transform, getTimestamp(), parent_frame, footprint_frame_));



    parent_frame = footprint_frame_;

    transform = footprint_transform.inverseTimes(transform);

  }



  if (!stabilized_frame_.empty()) {

    tf::Transform stabilized_transform(transform);

#ifdef TF_MATRIX3x3_H

    tf::Matrix3x3 rollpitch_rotation; rollpitch_rotation.setEulerYPR(0.0, p, r);

#else

    btMatrix3x3 rollpitch_rotation; rollpitch_rotation.setEulerYPR(0.0, p, r);

#endif

    stabilized_transform = stabilized_transform * tf::Transform(rollpitch_rotation.inverse());

    transforms.push_back(tf::StampedTransform(stabilized_transform, getTimestamp(), parent_frame, stabilized_frame_));



    parent_frame = stabilized_frame_;

    transform = stabilized_transform.inverseTimes(transform);

  }



  transforms.push_back(tf::StampedTransform(transform, getTimestamp(), parent_frame, base_frame_));



//  transforms.resize(3);



//  transforms[0].stamp_ = getTimestamp();

//  transforms[0].frame_id_ = nav_frame_;

//  transforms[0].child_frame_id_ = footprint_frame_;

//  transforms[0].setOrigin(tf::Point(position.x(), position.y(), 0.0));

//  rotation.setEulerYPR(y,0.0,0.0);

//  transforms[0].setBasis(rotation);



//  transforms[1].stamp_ = getTimestamp();

//  transforms[1].frame_id_ = footprint_frame_;

//  transforms[1].child_frame_id_ = stabilized_frame_;

//  transforms[1].setIdentity();

//  transforms[1].setOrigin(tf::Point(0.0, 0.0, position.z()));



//  transforms[2].stamp_ = getTimestamp();

//  transforms[2].frame_id_ = stabilized_frame_;

//  transforms[2].child_frame_id_ = base_frame_;

//  transforms[2].setIdentity();

//  rotation.setEulerYPR(0.0,p,r);

//  transforms[2].setBasis(rotation);

}



void PoseEstimation::updateWorldToOtherTransform(tf::StampedTransform& world_to_other_transform) {

  world_to_other_transform.frame_id_ = world_frame_;



  double y,p,r;

  world_to_other_transform.getBasis().getEulerYPR(y,p,r);

  if (!(getSystemStatus() & (STATE_ROLLPITCH   | STATE_PSEUDO_ROLLPITCH)))   { r = p = 0.0; }

  if (!(getSystemStatus() & (STATE_YAW         | STATE_PSEUDO_YAW)))         { y = 0.0; }

  if (!(getSystemStatus() & (STATE_POSITION_XY | STATE_PSEUDO_POSITION_XY))) { world_to_other_transform.getOrigin().setX(0.0); world_to_other_transform.getOrigin().setY(0.0); }

  if (!(getSystemStatus() & (STATE_POSITION_Z  | STATE_PSEUDO_POSITION_Z)))  { world_to_other_transform.getOrigin().setZ(0.0); }

  world_to_other_transform.getBasis().setEulerYPR(y, p, r);

}



bool PoseEstimation::getWorldToNavTransform(geometry_msgs::TransformStamped& transform) {

  return globalReference()->getWorldToNavTransform(transform, world_frame_, nav_frame_, getTimestamp());

}



const GlobalReferencePtr &PoseEstimation::globalReference() {

  return GlobalReference::Instance();

}



} // namespace hector_pose_estimation

=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation_core/src/measurement.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/measurement.h>

#include <ros/console.h>

namespace hector_pose_estimation {

Measurement::Measurement(const std::string& name)
  : name_(name)
  , status_flags_(0)
  , enabled_(true)
  , min_interval_(0.0)
  , timeout_(0.0)
  , timer_(0.0)
{
  parameters().add("enabled", enabled_);
  parameters().add("timeout", timeout_);
  parameters().add("min_interval", min_interval_);
  parameters().add("timeout", timeout_);
}

Measurement::~Measurement()
{
}

bool Measurement::init(PoseEstimation& estimator, State& state)
{
  if (getModel() && !getModel()->init(estimator, *this, state)) return false;
  if (!onInit(estimator)) return false;
  return true;
}

void Measurement::cleanup()
{
  if (getModel()) getModel()->cleanup();
  onCleanup();
}

void Measurement::reset(State& state)
{
  queue().clear();
  timer_ = 0;
  status_flags_ = 0;

  if (getModel()) getModel()->reset(state);
  onReset();
}

bool Measurement::active(const State& state) {
  bool active = enabled() && (getModel() ? getModel()->active(state) : !(state.getSystemStatus() & STATUS_ALIGNMENT));
  if (!active) status_flags_ = 0;
  if (min_interval_ > 0.0 && timer_ < min_interval_) return false;
  return active;
}

void Measurement::increase_timer(double dt) {
  timer_ += dt;
}

bool Measurement::timedout() const {
  return timeout_ > 0.0 && timer_ > timeout_;
}

void Measurement::add(const MeasurementUpdate& update) {
  queue().push(update);
}

bool Measurement::process() {
  bool result = true;

  while(!(queue().empty())) {
    result &= update(queue().pop());
  }

  // check for timeout
  if (timedout()) {
    if (status_flags_ > 0) ROS_WARN("Measurement %s timed out.", getName().c_str());
    status_flags_ = 0;
  }
  return result;
}

bool Measurement::update(const MeasurementUpdate &update)
{
  if (!filter() || !active(filter()->state())) return false;

  if (!updateImpl(update)) return false;
  filter()->state().updated();

  timer_ = 0.0;
  if (getModel()) status_flags_ = getModel()->getStatusFlags();
  return true;
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation/src/pose_estimation_nodelet.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/pose_estimation_node.h>
#include <nodelet/nodelet.h>

namespace hector_pose_estimation {

class PoseEstimationNodelet : public PoseEstimationNode, public nodelet::Nodelet
{
public:
  PoseEstimationNodelet(const SystemPtr& system = SystemPtr())
    : PoseEstimationNode(system)
  {}

private:
  void onInit() {
    PoseEstimationNode::init();
  }

  void onReset() {
    PoseEstimationNode::reset();
  }

  void onCleanup() {
    PoseEstimationNode::cleanup();
  }
};

} // namespace hector_pose_estimation

#include <pluginlib/class_list_macros.h>
PLUGINLIB_EXPORT_CLASS(hector_pose_estimation::PoseEstimationNodelet, nodelet::Nodelet)
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation/src/main.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/pose_estimation_node.h>

int main(int argc, char **argv) {
  ros::init(argc, argv, "pose_estimation");
  hector_pose_estimation::PoseEstimationNode node;
  if (!node.init()) return 1;

  ros::spin();

  node.cleanup();
  return 0;
}
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/hector_pose_estimation/src/pose_estimation_node.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/pose_estimation_node.h>
#include <hector_pose_estimation/ros/parameters.h>

#include <hector_pose_estimation/system/generic_quaternion_system_model.h>
#include <hector_pose_estimation/system/imu_input.h>
#include <hector_pose_estimation/system/imu_model.h>
#include <hector_pose_estimation/measurements/poseupdate.h>
#include <hector_pose_estimation/measurements/baro.h>
#include <hector_pose_estimation/measurements/height.h>
#include <hector_pose_estimation/measurements/magnetic.h>
#include <hector_pose_estimation/measurements/gps.h>

#ifdef USE_HECTOR_TIMING
  #include <hector_diagnostics/timing.h>
#endif

namespace hector_pose_estimation {

PoseEstimationNode::PoseEstimationNode(const SystemPtr& system, const StatePtr& state)
  : pose_estimation_(new PoseEstimation(system, state))
  , private_nh_("~")
  , transform_listener_(0)
  , world_nav_transform_updated_(true)
  , world_nav_transform_valid_(false)
  , sensor_pose_roll_(0), sensor_pose_pitch_(0), sensor_pose_yaw_(0)
{
  if (!system) pose_estimation_->addSystem(new GenericQuaternionSystemModel);

  pose_estimation_->addInput(new ImuInput, "imu");
  pose_estimation_->addMeasurement(new PoseUpdate("poseupdate"));
#if defined(USE_HECTOR_UAV_MSGS)
  pose_estimation_->addMeasurement(new Baro("baro"));
#endif
  pose_estimation_->addMeasurement(new Height("height"));
  pose_estimation_->addMeasurement(new Magnetic("magnetic"));
  pose_estimation_->addMeasurement(new GPS("gps"));
}

PoseEstimationNode::~PoseEstimationNode()
{
  cleanup();
  delete pose_estimation_;
  delete transform_listener_;
}

bool PoseEstimationNode::init() {
  // get parameters
  pose_estimation_->parameters().initialize(ParameterRegistryROS(getPrivateNodeHandle()));
  getPrivateNodeHandle().getParam("publish_covariances", publish_covariances_ = false);
  if (getPrivateNodeHandle().getParam("publish_world_map_transform", publish_world_other_transform_ = false)) {
    ROS_WARN("Parameter 'publish_world_map_transform' is deprecated. Use 'publish_world_other_transform' and 'other_frame' instead.");
  }
  getPrivateNodeHandle().getParam("publish_world_other_transform", publish_world_other_transform_);
  if (getPrivateNodeHandle().getParam("map_frame", other_frame_ = std::string())) {
    ROS_WARN("Parameter 'map_frame' is deprecated. Use 'other_frame' instead.");
  }
  getPrivateNodeHandle().getParam("other_frame", other_frame_);
  getPrivateNodeHandle().getParam("publish_world_nav_transform", publish_world_nav_transform_ = false);

  // search tf_prefix parameter
  tf_prefix_ = tf::getPrefixParam(getPrivateNodeHandle());
  if (!tf_prefix_.empty()) ROS_INFO("Using tf_prefix '%s'", tf_prefix_.c_str());

  // initialize pose estimation
  if (!pose_estimation_->init()) {
    ROS_ERROR("Intitialization of pose estimation failed!");
    return false;
  }

  imu_subscriber_        = getNodeHandle().subscribe("raw_imu", 10, &PoseEstimationNode::imuCallback, this);
  ahrs_subscriber_       = getNodeHandle().subscribe("ahrs", 10, &PoseEstimationNode::ahrsCallback, this);
  rollpitch_subscriber_  = getNodeHandle().subscribe("rollpitch", 10, &PoseEstimationNode::rollpitchCallback, this);
#if defined(USE_HECTOR_UAV_MSGS)
  baro_subscriber_       = getNodeHandle().subscribe("altimeter", 10, &PoseEstimationNode::baroCallback, this);
#else
  height_subscriber_     = getNodeHandle().subscribe("pressure_height", 10, &PoseEstimationNode::heightCallback, this);
#endif
  magnetic_subscriber_   = getNodeHandle().subscribe("magnetic", 10, &PoseEstimationNode::magneticCallback, this);

  gps_subscriber_.subscribe(getNodeHandle(), "fix", 10);
  gps_velocity_subscriber_.subscribe(getNodeHandle(), "fix_velocity", 10);
  gps_synchronizer_ = new message_filters::TimeSynchronizer<sensor_msgs::NavSatFix,geometry_msgs::Vector3Stamped>(gps_subscriber_, gps_velocity_subscriber_, 10);
  gps_synchronizer_->registerCallback(&PoseEstimationNode::gpsCallback, this);

  state_publisher_       = getNodeHandle().advertise<nav_msgs::Odometry>("state", 10, false);
  pose_publisher_        = getNodeHandle().advertise<geometry_msgs::PoseStamped>("pose", 10, false);
  velocity_publisher_    = getNodeHandle().advertise<geometry_msgs::Vector3Stamped>("velocity", 10, false);
  imu_publisher_         = getNodeHandle().advertise<sensor_msgs::Imu>("imu", 10, false);
  geopose_publisher_     = getNodeHandle().advertise<geographic_msgs::GeoPose>("geopose", 10, false);
  global_fix_publisher_  = getNodeHandle().advertise<sensor_msgs::NavSatFix>("global", 10, false);
  euler_publisher_       = getNodeHandle().advertise<geometry_msgs::Vector3Stamped>("euler", 10, false);

  angular_velocity_bias_publisher_    = getNodeHandle().advertise<geometry_msgs::Vector3Stamped>("angular_velocity_bias", 10, false);
  linear_acceleration_bias_publisher_ = getNodeHandle().advertise<geometry_msgs::Vector3Stamped>("linear_acceleration_bias", 10, false);
  gps_pose_publisher_                 = getNodeHandle().advertise<geometry_msgs::PoseStamped>("fix/pose", 10, false);
  sensor_pose_publisher_              = getNodeHandle().advertise<geometry_msgs::PoseStamped>("sensor_pose", 10, false);

  poseupdate_subscriber_  = getNodeHandle().subscribe("poseupdate", 10, &PoseEstimationNode::poseupdateCallback, this);
  twistupdate_subscriber_ = getNodeHandle().subscribe("twistupdate", 10, &PoseEstimationNode::twistupdateCallback, this);
  syscommand_subscriber_  = getNodeHandle().subscribe("syscommand", 10, &PoseEstimationNode::syscommandCallback, this);

#ifdef USE_HECTOR_TIMING
  timing_publisher_ = getPrivateNodeHandle().advertise<hector_diagnostics::TimingAggregator>("timing", 10, false);
#endif

  global_reference_publisher_  = getNodeHandle().advertise<geographic_msgs::GeoPose>("global/reference", 1, true);
  pose_estimation_->globalReference()->addUpdateCallback(boost::bind(&PoseEstimationNode::globalReferenceUpdated, this));

  // setup publish_world_nav_transform timer
  if (publish_world_nav_transform_) {
    double period = 0.1;
    getPrivateNodeHandle().getParam("publish_world_nav_transform_period", period);
    publish_world_nav_transform_period_ = ros::Duration(period);
    publish_world_nav_transform_timer_ = getNodeHandle().createTimer(publish_world_nav_transform_period_, &PoseEstimationNode::publishWorldNavTransform, this,
                                                                     /* oneshot = */ false, /* autostart = */ true);
  }

  // publish initial state
  publish();

  return true;
}

void PoseEstimationNode::reset() {
  pose_estimation_->reset();

  sensor_pose_ = geometry_msgs::PoseStamped();
  sensor_pose_roll_  = 0.0;
  sensor_pose_pitch_ = 0.0;
  sensor_pose_yaw_   = 0.0;
}

void PoseEstimationNode::cleanup() {
  if (gps_synchronizer_) {
    delete gps_synchronizer_;
    gps_synchronizer_ = 0;
  }
  publish_world_nav_transform_timer_.stop();

  pose_estimation_->cleanup();
}

void PoseEstimationNode::imuCallback(const sensor_msgs::ImuConstPtr& imu) {
  pose_estimation_->setInput(ImuInput(*imu));
  pose_estimation_->update(imu->header.stamp);

  // calculate roll and pitch purely from acceleration
  if (sensor_pose_publisher_) {
    tf::Vector3 linear_acceleration_body(imu->linear_acceleration.x, imu->linear_acceleration.y, imu->linear_acceleration.z);
    linear_acceleration_body.normalize();
    sensor_pose_roll_  =  atan2(linear_acceleration_body.y(), linear_acceleration_body.z());
    sensor_pose_pitch_ = -asin(linear_acceleration_body.x());
  }

  publish();
}

void PoseEstimationNode::ahrsCallback(const sensor_msgs::ImuConstPtr& ahrs) {
  pose_estimation_->state().setOrientation(Quaternion(ahrs->orientation.w, ahrs->orientation.x, ahrs->orientation.y, ahrs->orientation.z));
  pose_estimation_->setInput(ImuInput(*ahrs));
  pose_estimation_->update(ahrs->header.stamp);
  publish();
}

void PoseEstimationNode::rollpitchCallback(const sensor_msgs::ImuConstPtr& attitude) {
  pose_estimation_->state().setRollPitch(Quaternion(attitude->orientation.w, attitude->orientation.x, attitude->orientation.y, attitude->orientation.z));
  pose_estimation_->setInput(ImuInput(*attitude));
  pose_estimation_->update(attitude->header.stamp);
  publish();
}

#if defined(USE_HECTOR_UAV_MSGS)
void PoseEstimationNode::baroCallback(const hector_uav_msgs::AltimeterConstPtr& altimeter) {
  boost::shared_ptr<Baro> m = boost::static_pointer_cast<Baro>(pose_estimation_->getMeasurement("baro"));
  m->add(Baro::Update(altimeter->pressure, altimeter->qnh));
}

#else
void PoseEstimationNode::heightCallback(const geometry_msgs::PointStampedConstPtr& height) {
  boost::shared_ptr<Height> m = boost::static_pointer_cast<Height>(pose_estimation_->getMeasurement("height"));

  Height::MeasurementVector update;
  update(0) = height->point.z;
  m->add(Height::Update(update));

  if (sensor_pose_publisher_) {
    sensor_pose_.pose.position.z = height->point.z - m->getElevation();
  }
}
#endif

void PoseEstimationNode::magneticCallback(const geometry_msgs::Vector3StampedConstPtr& magnetic) {
  boost::shared_ptr<Magnetic> m = boost::static_pointer_cast<Magnetic>(pose_estimation_->getMeasurement("magnetic"));

  Magnetic::MeasurementVector update;
  update.x() = magnetic->vector.x;
  update.y() = magnetic->vector.y;
  update.z() = magnetic->vector.z;
  m->add(Magnetic::Update(update));

  if (sensor_pose_publisher_) {
    sensor_pose_yaw_ = -(m->getModel()->getTrueHeading(pose_estimation_->state(), update) - pose_estimation_->globalReference()->heading());
  }
}

void PoseEstimationNode::gpsCallback(const sensor_msgs::NavSatFixConstPtr& gps, const geometry_msgs::Vector3StampedConstPtr& gps_velocity) {
  boost::shared_ptr<GPS> m = boost::static_pointer_cast<GPS>(pose_estimation_->getMeasurement("gps"));

  if (gps->status.status == sensor_msgs::NavSatStatus::STATUS_NO_FIX) {
    if (m->getStatusFlags() > 0) m->reset(pose_estimation_->state());
    return;
  }

  GPS::Update update;
  update.latitude = gps->latitude * M_PI/180.0;
  update.longitude = gps->longitude * M_PI/180.0;
  update.velocity_north =  gps_velocity->vector.x;
  update.velocity_east  = -gps_velocity->vector.y;
  m->add(update);

  if (gps_pose_publisher_ || sensor_pose_publisher_) {
    geometry_msgs::PoseStamped gps_pose;
    pose_estimation_->getHeader(gps_pose.header);
    gps_pose.header.stamp = gps->header.stamp;
    GPSModel::MeasurementVector y = m->getVector(update, pose_estimation_->state());

    if (gps_pose_publisher_) {
      gps_pose.pose.position.x = y(0);
      gps_pose.pose.position.y = y(1);
      gps_pose.pose.position.z = gps->altitude - pose_estimation_->globalReference()->position().altitude;
      double track = atan2(gps_velocity->vector.y, gps_velocity->vector.x);
      gps_pose.pose.orientation.w = cos(track/2);
      gps_pose.pose.orientation.z = sin(track/2);
      gps_pose_publisher_.publish(gps_pose);
    }

    sensor_pose_.pose.position.x = y(0);
    sensor_pose_.pose.position.y = y(1);
  }
}

void PoseEstimationNode::poseupdateCallback(const geometry_msgs::PoseWithCovarianceStampedConstPtr& pose) {
  pose_estimation_->getMeasurement("poseupdate")->add(PoseUpdate::Update(pose));

  if (sensor_pose_publisher_) {
    if (pose->pose.covariance[0] > 0)  sensor_pose_.pose.position.x = pose->pose.pose.position.x;
    if (pose->pose.covariance[7] > 0)  sensor_pose_.pose.position.y = pose->pose.pose.position.y;
    if (pose->pose.covariance[14] > 0) sensor_pose_.pose.position.z = pose->pose.pose.position.z;
    tf::Quaternion q;
    double yaw, pitch, roll;
    tf::quaternionMsgToTF(pose->pose.pose.orientation, q);
    tf::Matrix3x3(q).getEulerYPR(yaw, pitch, roll);
    if (pose->pose.covariance[21] > 0) sensor_pose_roll_  = roll;
    if (pose->pose.covariance[28] > 0) sensor_pose_pitch_ = pitch;
    if (pose->pose.covariance[35] > 0) sensor_pose_yaw_   = yaw;
  }
}

void PoseEstimationNode::twistupdateCallback(const geometry_msgs::TwistWithCovarianceStampedConstPtr& twist) {
  pose_estimation_->getMeasurement("poseupdate")->add(PoseUpdate::Update(twist));
}

void PoseEstimationNode::syscommandCallback(const std_msgs::StringConstPtr& syscommand) {
  if (syscommand->data == "reset") {
    ROS_INFO("Resetting pose_estimation");
    pose_estimation_->reset();
    publish();
  }
}

void PoseEstimationNode::globalReferenceUpdated() {
  geographic_msgs::GeoPose geopose;
  pose_estimation_->globalReference()->getGeoPose(geopose);
  global_reference_publisher_.publish(geopose);

  // update world nav transform
  world_nav_transform_updated_ = true;
}

void PoseEstimationNode::publishWorldNavTransform(const ros::TimerEvent &) {
  if (world_nav_transform_updated_) {
    world_nav_transform_valid_ = pose_estimation_->getWorldToNavTransform(world_nav_transform_);
    world_nav_transform_updated_ = false;
  }

  if (world_nav_transform_valid_) {
    world_nav_transform_.header.stamp = ros::Time::now() + publish_world_nav_transform_period_;
    getTransformBroadcaster()->sendTransform(world_nav_transform_);
  }
}

void PoseEstimationNode::publish() {
  if (state_publisher_ && state_publisher_.getNumSubscribers() > 0) {
    nav_msgs::Odometry state;
    pose_estimation_->getState(state, publish_covariances_);
    state_publisher_.publish(state);
  }

  if (pose_publisher_ && pose_publisher_.getNumSubscribers() > 0) {
    geometry_msgs::PoseStamped pose_msg;
    pose_estimation_->getPose(pose_msg);
    pose_publisher_.publish(pose_msg);
  }

  if (imu_publisher_ && imu_publisher_.getNumSubscribers() > 0) {
    sensor_msgs::Imu imu_msg;
    pose_estimation_->getHeader(imu_msg.header);
    pose_estimation_->getOrientation(imu_msg.orientation);
    pose_estimation_->getImuWithBiases(imu_msg.linear_acceleration, imu_msg.angular_velocity);
    imu_publisher_.publish(imu_msg);
  }

  if (velocity_publisher_ && velocity_publisher_.getNumSubscribers() > 0) {
    geometry_msgs::Vector3Stamped velocity_msg;
    pose_estimation_->getVelocity(velocity_msg);
    velocity_publisher_.publish(velocity_msg);
  }

  if (geopose_publisher_ && geopose_publisher_.getNumSubscribers() > 0) {
    geographic_msgs::GeoPose geopose_msg;
    pose_estimation_->getGlobal(geopose_msg);
    geopose_publisher_.publish(geopose_msg);
  }

  if (global_fix_publisher_ && global_fix_publisher_.getNumSubscribers() > 0) {
    sensor_msgs::NavSatFix global_msg;
    pose_estimation_->getGlobal(global_msg);
    global_fix_publisher_.publish(global_msg);
  }

  if (euler_publisher_ && euler_publisher_.getNumSubscribers() > 0) {
    geometry_msgs::Vector3Stamped euler_msg;
    pose_estimation_->getHeader(euler_msg.header);
    pose_estimation_->getOrientation(euler_msg.vector.z, euler_msg.vector.y, euler_msg.vector.x);
    euler_publisher_.publish(euler_msg);
  }

  if ((angular_velocity_bias_publisher_ && angular_velocity_bias_publisher_.getNumSubscribers() > 0) ||
      (linear_acceleration_bias_publisher_ && linear_acceleration_bias_publisher_.getNumSubscribers() > 0)) {
    geometry_msgs::Vector3Stamped angular_velocity_msg, linear_acceleration_msg;
    pose_estimation_->getBias(angular_velocity_msg, linear_acceleration_msg);
    if (angular_velocity_bias_publisher_) angular_velocity_bias_publisher_.publish(angular_velocity_msg);
    if (linear_acceleration_bias_publisher_) linear_acceleration_bias_publisher_.publish(linear_acceleration_msg);
  }

  if (sensor_pose_publisher_ && sensor_pose_publisher_.getNumSubscribers() > 0) {
    pose_estimation_->getHeader(sensor_pose_.header);
    tf::Quaternion orientation;
    orientation.setRPY(sensor_pose_roll_, sensor_pose_pitch_, sensor_pose_yaw_);
    tf::quaternionTFToMsg(orientation, sensor_pose_.pose.orientation);
    sensor_pose_publisher_.publish(sensor_pose_);
  }

  if (getTransformBroadcaster())
  {
    transforms_.clear();

    pose_estimation_->getTransforms(transforms_);

    if (publish_world_other_transform_) {
      tf::StampedTransform world_to_other_transform;
      std::string nav_frame = pose_estimation_->parameters().getAs<std::string>("nav_frame");
      try {
        getTransformListener()->lookupTransform(nav_frame, other_frame_, ros::Time(), world_to_other_transform);
        pose_estimation_->updateWorldToOtherTransform(world_to_other_transform);
        transforms_.push_back(world_to_other_transform);

      } catch (tf::TransformException& e) {
        ROS_WARN("Could not find a transformation from %s to %s to publish the world transformation", nav_frame.c_str(), other_frame_.c_str());
      }
    }

    // resolve tf frames
    for(std::vector<tf::StampedTransform>::iterator t = transforms_.begin(); t != transforms_.end(); t++) {
      t->frame_id_       = tf::resolve(tf_prefix_, t->frame_id_);
      t->child_frame_id_ = tf::resolve(tf_prefix_, t->child_frame_id_);
    }

    getTransformBroadcaster()->sendTransform(transforms_);
  }

#ifdef USE_HECTOR_TIMING
  timing_publisher_.publish(*hector_diagnostics::TimingAggregator::Instance());
#endif
}

tf::TransformBroadcaster *PoseEstimationNode::getTransformBroadcaster() {
  return &transform_broadcaster_;
}

tf::TransformListener *PoseEstimationNode::getTransformListener() {
  if (!transform_listener_) transform_listener_ = new tf::TransformListener();
  return transform_listener_;
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/rtt_hector_pose_estimation/src/parameters.cpp
//=================================================================================================
// Copyright (c) 2013, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include "parameters.h"
#include <hector_pose_estimation/matrix.h>

#include <rtt/Logger.hpp>
#include <rtt/PropertyBag.hpp>
#include <rtt/Property.hpp>

namespace hector_pose_estimation {

ParameterRegistryProperties::ParameterRegistryProperties(RTT::PropertyBag *properties, bool set_all)
  : properties_(properties)
  , set_all_(set_all)
{
}

void ParameterRegistryProperties::operator ()(ParameterPtr parameter) {
  properties_->removeProperty(properties_->getProperty(parameter->key));
  if (parameter->hasType<std::string>()) { properties_->addProperty(parameter->key, parameter->as<std::string>()); return; }
  if (parameter->hasType<double>())      { properties_->addProperty(parameter->key, parameter->as<double>()); return; }
  if (parameter->hasType<int>())         { properties_->addProperty(parameter->key, parameter->as<int>()); return; }
  if (parameter->hasType<bool>())        { properties_->addProperty(parameter->key, parameter->as<bool>()); return; }

  if (parameter->hasType< std::vector<double> >()) {
    // TODO
    // return;
  }

  if (parameter->hasType< std::vector<ColumnVector> >()) {
    // TODO
    // return;
  }

  RTT::log(RTT::Error) << "Could not register parameter " << parameter->key << " due to unknown type " << parameter->type() << "!" << RTT::endlog();
}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/rtt_hector_pose_estimation/src/services.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include "services.h"
#include "parameters.h"

namespace hector_pose_estimation {

SystemService::SystemService(RTT::TaskContext *owner, const SystemPtr& system, const std::string& name)
  : RTT::Service(name.empty() ? system->getName() : name, owner)
{
  system->parameters().initialize(ParameterRegistryProperties(this->properties()));
}

SystemService::~SystemService()
{}

MeasurementService::MeasurementService(RTT::TaskContext *owner, const MeasurementPtr& measurement, const std::string& name)
  : RTT::Service(name.empty() ? measurement->getName() : name, owner)
{
  measurement->parameters().initialize(ParameterRegistryProperties(this->properties()));
}

MeasurementService::~MeasurementService()
{}

} // namespace hector_pose_estimation
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/rtt_hector_pose_estimation/src/taskcontext.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer and Martin Nowara, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_pose_estimation/rtt/taskcontext.h>
#include "services.h"
#include "parameters.h"

#include <hector_pose_estimation/system/generic_quaternion_system_model.h>
#include <hector_pose_estimation/measurements/poseupdate.h>
#include <hector_pose_estimation/measurements/height.h>
#include <hector_pose_estimation/measurements/magnetic.h>
#include <hector_pose_estimation/measurements/gps.h>

#include <hector_pose_estimation/ros/parameters.h>

#include <rtt/Component.hpp>

namespace hector_pose_estimation {

PoseEstimationTaskContext::PoseEstimationTaskContext(const std::string& name, const SystemPtr& system)
  : RTT::TaskContext(name, RTT::TaskContext::PreOperational)
  , PoseEstimation(system)
{
  if (!system) addSystem(System::create(new GenericQuaternionSystemModel));

  this->addEventPort("raw_imu", imu_input_);
  this->addPort("poseupdate", pose_update_input_);
  this->addPort("magnetic", magnetic_input_);
  this->addPort("pressure_height", height_input_);
  this->addPort("fix", gps_input_);
  this->addPort("fix_velocity", gps_velocity_input_);
  this->addPort("state", state_output_);
  this->addPort("imu", imu_output_);
  this->addPort("pose", pose_output_);
  this->addPort("velocity", velocity_output_);
  this->addPort("global", global_position_output_);
  this->addPort("angular_velocity_bias", angular_velocity_bias_output_);
  this->addPort("linear_acceleration_bias", linear_acceleration_bias_output_);
  this->addEventPort("syscommand", command_input_, boost::bind(&PoseEstimationTaskContext::commandCallback, this, _1));

  param_namespace_ = "~";
  this->addProperty("param_namespace", param_namespace_);

  this->addOperation("reset", &PoseEstimationTaskContext::reset, this, RTT::OwnThread);
  this->addOperation("getSystemStatus", &PoseEstimation::getSystemStatus, static_cast<PoseEstimation *>(this), RTT::OwnThread);
  this->addOperation("getMeasurementStatus", &PoseEstimation::getMeasurementStatus, static_cast<PoseEstimation *>(this), RTT::OwnThread);

  PoseEstimation::addMeasurement(new PoseUpdate("PoseUpdate"));
  PoseEstimation::addMeasurement(new Height("Height"));
  PoseEstimation::addMeasurement(new Magnetic("Magnetic"));
  PoseEstimation::addMeasurement(new GPS("GPS"));

  for(Systems::iterator it = systems_.begin(); it != systems_.end(); ++it) {
    this->provides()->addService(RTT::Service::shared_ptr(new SystemService(this, *it)));
  }

  for(Measurements::iterator it = measurements_.begin(); it != measurements_.end(); ++it) {
    this->provides()->addService(RTT::Service::shared_ptr(new MeasurementService(this, *it)));
  }

  PoseEstimation::parameters().initialize(ParameterRegistryROS(ros::NodeHandle(param_namespace_)));
  PoseEstimation::parameters().initialize(ParameterRegistryProperties(this->properties()));
}

PoseEstimationTaskContext::~PoseEstimationTaskContext()
{
  stop();
  cleanup();
}

bool PoseEstimationTaskContext::configureHook()
{
  if (!PoseEstimation::init()) {
    RTT::log(RTT::Error) << "Intitialization of pose estimation failed!" << RTT::endlog();
    return false;
  }

  // publish initial state
  updateOutputs();

  return true;
}

void PoseEstimationTaskContext::cleanupHook()
{
  PoseEstimation::cleanup();
}

bool PoseEstimationTaskContext::startHook()
{
  return true;
}

void PoseEstimationTaskContext::stopHook()
{
}

void PoseEstimationTaskContext::updateHook()
{
  while(magnetic_input_.read(magnetic_) == RTT::NewData && PoseEstimation::getMeasurement("Magnetic"))
  {
    Magnetic::MeasurementVector update(3);
    update.x() = magnetic_.vector.x;
    update.y() = magnetic_.vector.y;
    update.z() = magnetic_.vector.z;
    PoseEstimation::getMeasurement("Magnetic")->add(Magnetic::Update(update));
  }

  while(height_input_.read(height_) == RTT::NewData && PoseEstimation::getMeasurement("Height"))
  {
    Height::Update update;
#ifdef USE_MAV_MSGS
    update = height_.height;
#else
    update = height_.point.z;
#endif
    PoseEstimation::getMeasurement("Height")->add(update);
  }

  while(gps_input_.read(gps_) == RTT::NewData && gps_velocity_input_.read(gps_velocity_) != RTT::NoData && PoseEstimation::getMeasurement("GPS"))
  {
    if (gps_.status.status != sensor_msgs::NavSatStatus::STATUS_NO_FIX)
    {
      GPS::Update update;
      update.latitude = gps_.latitude * M_PI/180.0;
      update.longitude = gps_.longitude * M_PI/180.0;
      update.velocity_north =  gps_velocity_.vector.x;
      update.velocity_east  = -gps_velocity_.vector.y;
      PoseEstimation::getMeasurement("GPS")->add(update);
    }
  }

  while(pose_update_input_.read(pose_update_) == RTT::NewData && PoseEstimation::getMeasurement("PoseUpdate"))
  {
    PoseEstimation::getMeasurement("PoseUpdate")->add(PoseUpdate::Update(pose_update_));
  }

  while(imu_input_.read(imu_in_) == RTT::NewData) {
    setInput(ImuInput(imu_in_));
    PoseEstimation::update(imu_in_.header.stamp);
    updateOutputs();
  }
}

void PoseEstimationTaskContext::updateOutputs()
{
  getState(state_);
  state_output_.write(state_);

  if (imu_output_.connected()) {
    imu_out_.header = state_.header;
    imu_out_.orientation = state_.pose.pose.orientation;
    getImuWithBiases(imu_out_.linear_acceleration, imu_out_.angular_velocity);
    imu_output_.write(imu_out_);
  }

  if (pose_output_.connected()) {
    pose_.header = state_.header;
    pose_.pose = state_.pose.pose;
    pose_output_.write(pose_);
  }

  if (velocity_output_.connected()) {
    velocity_.header = state_.header;
    velocity_.vector = state_.twist.twist.linear;
    velocity_output_.write(velocity_);
  }

  if (global_position_output_.connected()) {
    getGlobalPosition(global_position_);
    global_position_output_.write(global_position_);
  }

  if (angular_velocity_bias_output_.connected() || linear_acceleration_bias_output_.connected()) {
    getHeader(angular_velocity_bias_.header);
    getHeader(linear_acceleration_bias_.header);
    getBias(angular_velocity_bias_, linear_acceleration_bias_);
    angular_velocity_bias_output_.write(angular_velocity_bias_);
    linear_acceleration_bias_output_.write(linear_acceleration_bias_);
  }
}


void PoseEstimationTaskContext::reset() {
  RTT::log(RTT::Info) << "Resetting pose_estimation" << RTT::endlog();
  PoseEstimation::reset();
  updateOutputs();
}

void PoseEstimationTaskContext::commandCallback(RTT::base::PortInterface *) {
  std_msgs::String command;
  if (command_input_.read(command) == RTT::NewData) {
    if (command.data == "reset") reset();
  }
}

} // namespace

ORO_CREATE_COMPONENT( hector_pose_estimation::PoseEstimationTaskContext )
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_localization/message_to_tf/src/message_to_tf.cpp
#include <ros/ros.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/Vector3Stamped.h>
#include <geometry_msgs/TransformStamped.h>
#include <sensor_msgs/Imu.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_listener.h> // for tf::getPrefixParam()
#include <tf/transform_datatypes.h>

#include <topic_tools/shape_shifter.h>

std::string g_odometry_topic;
std::string g_pose_topic;
std::string g_imu_topic;
std::string g_topic;
std::string g_frame_id;
std::string g_footprint_frame_id;
std::string g_position_frame_id;
std::string g_stabilized_frame_id;
std::string g_child_frame_id;

bool g_publish_roll_pitch;

std::string g_tf_prefix;

tf::TransformBroadcaster *g_transform_broadcaster;
ros::Publisher g_pose_publisher;
ros::Publisher g_euler_publisher;

#ifndef TF_MATRIX3x3_H
  typedef btScalar tfScalar;
  namespace tf { typedef btMatrix3x3 Matrix3x3; }
#endif

void addTransform(std::vector<geometry_msgs::TransformStamped>& transforms, const tf::StampedTransform& tf)
{
  transforms.push_back(geometry_msgs::TransformStamped());
  tf::transformStampedTFToMsg(tf, transforms.back());
}

void sendTransform(geometry_msgs::Pose const &pose, const std_msgs::Header& header, std::string child_frame_id = "")
{
  std::vector<geometry_msgs::TransformStamped> transforms;

  tf::StampedTransform tf;
  tf.stamp_ = header.stamp;

  tf.frame_id_ = header.frame_id;
  if (!g_frame_id.empty()) tf.frame_id_ = g_frame_id;
  tf.frame_id_ = tf::resolve(g_tf_prefix, tf.frame_id_);

  if (!g_child_frame_id.empty()) child_frame_id = g_child_frame_id;
  if (child_frame_id.empty()) child_frame_id = "base_link";

  tf::Quaternion orientation;
  tf::quaternionMsgToTF(pose.orientation, orientation);
  tfScalar yaw, pitch, roll;
  tf::Matrix3x3(orientation).getEulerYPR(yaw, pitch, roll);
  tf::Point position;
  tf::pointMsgToTF(pose.position, position);

  // position intermediate transform (x,y,z)
  if( !g_position_frame_id.empty() && child_frame_id != g_position_frame_id) {
    tf.child_frame_id_ = tf::resolve(g_tf_prefix, g_position_frame_id);
    tf.setOrigin(tf::Vector3(position.x(), position.y(), position.z() ));
    tf.setRotation(tf::Quaternion(0.0, 0.0, 0.0, 1.0));
    addTransform(transforms, tf);
  }

  // footprint intermediate transform (x,y,yaw)
  if (!g_footprint_frame_id.empty() && child_frame_id != g_footprint_frame_id) {
    tf.child_frame_id_ = tf::resolve(g_tf_prefix, g_footprint_frame_id);
    tf.setOrigin(tf::Vector3(position.x(), position.y(), 0.0));
    tf.setRotation(tf::createQuaternionFromRPY(0.0, 0.0, yaw));
    addTransform(transforms, tf);

    yaw = 0.0;
    position.setX(0.0);
    position.setY(0.0);
    tf.frame_id_ = tf::resolve(g_tf_prefix, g_footprint_frame_id);
  }

  // stabilized intermediate transform (z)
  if (!g_footprint_frame_id.empty() && child_frame_id != g_stabilized_frame_id) {
    tf.child_frame_id_ = tf::resolve(g_tf_prefix, g_stabilized_frame_id);
    tf.setOrigin(tf::Vector3(0.0, 0.0, position.z()));
    tf.setBasis(tf::Matrix3x3::getIdentity());
    addTransform(transforms, tf);

    position.setZ(0.0);
    tf.frame_id_ = tf::resolve(g_tf_prefix, g_stabilized_frame_id);
  }

  // base_link transform (roll, pitch)
  if (g_publish_roll_pitch) {
    tf.child_frame_id_ = tf::resolve(g_tf_prefix, child_frame_id);
    tf.setOrigin(position);
    tf.setRotation(tf::createQuaternionFromRPY(roll, pitch, yaw));
    addTransform(transforms, tf);
  }

  g_transform_broadcaster->sendTransform(transforms);

  // publish pose message
  if (g_pose_publisher) {
    geometry_msgs::PoseStamped pose_stamped;
    pose_stamped.pose = pose;
    pose_stamped.header = header;
    g_pose_publisher.publish(pose_stamped);
  }

  // publish pose message
  if (g_euler_publisher) {
    geometry_msgs::Vector3Stamped euler_stamped;
    euler_stamped.vector.x = roll;
    euler_stamped.vector.y = pitch;
    euler_stamped.vector.z = yaw;
    euler_stamped.header = header;
    g_euler_publisher.publish(euler_stamped);
  }
}

void odomCallback(nav_msgs::Odometry const &odometry) {
  sendTransform(odometry.pose.pose, odometry.header, odometry.child_frame_id);
}

void poseCallback(geometry_msgs::PoseStamped const &pose) {
  sendTransform(pose.pose, pose.header);
}

void tfCallback(geometry_msgs::TransformStamped const &tf) {
  geometry_msgs::Pose pose;
  pose.position.x = tf.transform.translation.x;
  pose.position.y = tf.transform.translation.y;
  pose.position.z = tf.transform.translation.z;
  pose.orientation = tf.transform.rotation;

  sendTransform(pose, tf.header);
}

void imuCallback(sensor_msgs::Imu const &imu) {
  std::vector<geometry_msgs::TransformStamped> transforms;
  std::string child_frame_id;

  tf::StampedTransform tf;
  tf.stamp_ = imu.header.stamp;

  tf.frame_id_ = tf::resolve(g_tf_prefix, g_stabilized_frame_id);
  if (!g_child_frame_id.empty()) child_frame_id = g_child_frame_id;
  if (child_frame_id.empty()) child_frame_id = "base_link";

  tf::Quaternion orientation;
  tf::quaternionMsgToTF(imu.orientation, orientation);
  tfScalar yaw, pitch, roll;
  tf::Matrix3x3(orientation).getEulerYPR(yaw, pitch, roll);
  tf::Quaternion rollpitch = tf::createQuaternionFromRPY(roll, pitch, 0.0);

  // base_link transform (roll, pitch)
  if (g_publish_roll_pitch) {
    tf.child_frame_id_ = tf::resolve(g_tf_prefix, child_frame_id);
    tf.setOrigin(tf::Vector3(0.0, 0.0, 0.0));
    tf.setRotation(rollpitch);
    addTransform(transforms, tf);
  }

  if (!transforms.empty()) g_transform_broadcaster->sendTransform(transforms);

  // publish pose message
  if (g_pose_publisher) {
    geometry_msgs::PoseStamped pose_stamped;
    pose_stamped.header.stamp = imu.header.stamp;
    pose_stamped.header.frame_id = g_stabilized_frame_id;
    tf::quaternionTFToMsg(rollpitch, pose_stamped.pose.orientation);
    g_pose_publisher.publish(pose_stamped);
  }
}

void multiCallback(topic_tools::ShapeShifter const &input) {
  if (input.getDataType() == "nav_msgs/Odometry") {
    nav_msgs::Odometry::ConstPtr odom = input.instantiate<nav_msgs::Odometry>();
    odomCallback(*odom);
    return;
  }

  if (input.getDataType() == "geometry_msgs/PoseStamped") {
    geometry_msgs::PoseStamped::ConstPtr pose = input.instantiate<geometry_msgs::PoseStamped>();
    poseCallback(*pose);
    return;
  }

  if (input.getDataType() == "sensor_msgs/Imu") {
    sensor_msgs::Imu::ConstPtr imu = input.instantiate<sensor_msgs::Imu>();
    imuCallback(*imu);
    return;
  }

  if (input.getDataType() == "geometry_msgs/TransformStamped") {
    geometry_msgs::TransformStamped::ConstPtr tf = input.instantiate<geometry_msgs::TransformStamped>();
    tfCallback(*tf);
    return;
  }

  ROS_ERROR_THROTTLE(1.0, "message_to_tf received a %s message. Supported message types: nav_msgs/Odometry geometry_msgs/PoseStamped geometry_msgs/TransformStamped sensor_msgs/Imu", input.getDataType().c_str());
}

int main(int argc, char** argv) {
  ros::init(argc, argv, "message_to_tf");

  g_footprint_frame_id = "base_footprint";
  g_stabilized_frame_id = "base_stabilized";
  // g_position_frame_id = "base_position";
  // g_child_frame_id = "base_link";

  ros::NodeHandle priv_nh("~");
  priv_nh.getParam("odometry_topic", g_odometry_topic);
  priv_nh.getParam("pose_topic", g_pose_topic);
  priv_nh.getParam("imu_topic", g_imu_topic);
  priv_nh.getParam("topic", g_topic);
  priv_nh.getParam("frame_id", g_frame_id);
  priv_nh.getParam("footprint_frame_id", g_footprint_frame_id);
  priv_nh.getParam("position_frame_id", g_position_frame_id);
  priv_nh.getParam("stabilized_frame_id", g_stabilized_frame_id);
  priv_nh.getParam("child_frame_id", g_child_frame_id);

  // get topic from the commandline
  if (argc > 1) {
      g_topic = argv[1];
      g_odometry_topic.clear();
      g_pose_topic.clear();
      g_imu_topic.clear();
  }

  g_publish_roll_pitch = true;
  priv_nh.getParam("publish_roll_pitch", g_publish_roll_pitch);

  g_tf_prefix = tf::getPrefixParam(priv_nh);
  g_transform_broadcaster = new tf::TransformBroadcaster;

  ros::NodeHandle node;
  ros::Subscriber sub1, sub2, sub3, sub4;
  int subscribers = 0;
  if (!g_odometry_topic.empty()) {
      sub1 = node.subscribe(g_odometry_topic, 10, &odomCallback);
      subscribers++;
  }
  if (!g_pose_topic.empty()) {
      sub2 = node.subscribe(g_pose_topic, 10, &poseCallback);
      subscribers++;
  }
  if (!g_imu_topic.empty()) {
      sub3 = node.subscribe(g_imu_topic, 10, &imuCallback);
      subscribers++;
  }
  if (!g_topic.empty()) {
      sub4 = node.subscribe(g_topic, 10, &multiCallback);
      subscribers++;
  }

  if (subscribers == 0) {
    ROS_FATAL("Usage: rosrun message_to_tf message_to_tf <topic>");
    return 1;
  } else if (subscribers > 1) {
    ROS_FATAL("More than one of the parameters odometry_topic, pose_topic, imu_topic and topic are set.\n"
              "Please specify exactly one of them or simply add the topic name to the command line.");
    return 1;
  }

  bool publish_pose = true;
  priv_nh.getParam("publish_pose", publish_pose);
  if (publish_pose) {
    std::string publish_pose_topic;
    priv_nh.getParam("publish_pose_topic", publish_pose_topic);

    if (!publish_pose_topic.empty())
      g_pose_publisher = node.advertise<geometry_msgs::PoseStamped>(publish_pose_topic, 10);
    else
      g_pose_publisher = priv_nh.advertise<geometry_msgs::PoseStamped>("pose", 10);
  }

  bool publish_euler = true;
  priv_nh.getParam("publish_euler", publish_euler);
  if (publish_euler) {
    std::string publish_euler_topic;
    priv_nh.getParam("publish_euler_topic", publish_euler_topic);

    if (!publish_euler_topic.empty())
      g_euler_publisher = node.advertise<geometry_msgs::Vector3Stamped>(publish_euler_topic, 10);
    else
      g_euler_publisher = priv_nh.advertise<geometry_msgs::Vector3Stamped>("euler", 10);
  }

  ros::spin();
  delete g_transform_broadcaster;
  return 0;
}
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_quadrotor/hector_quadrotor_gazebo_plugins/src/gazebo_ros_baro.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_quadrotor_gazebo_plugins/gazebo_ros_baro.h>
#include <gazebo/common/Events.hh>
#include <gazebo/physics/physics.hh>
#include <ignition/math.hh>

static const double DEFAULT_ELEVATION = 0.0;
static const double DEFAULT_QNH       = 1013.25;

namespace gazebo {

GazeboRosBaro::GazeboRosBaro()
{
}

////////////////////////////////////////////////////////////////////////////////
// Destructor
GazeboRosBaro::~GazeboRosBaro()
{
  updateTimer.Disconnect(updateConnection);

  dynamic_reconfigure_server_.reset();

  node_handle_->shutdown();
  delete node_handle_;
}

////////////////////////////////////////////////////////////////////////////////
// Load the controller
void GazeboRosBaro::Load(physics::ModelPtr _model, sdf::ElementPtr _sdf)
{
  world = _model->GetWorld();
  link = _model->GetLink();
  link_name_ = link->GetName();

  if (_sdf->HasElement("bodyName"))
  {
    link_name_ = _sdf->GetElement("bodyName")->Get<std::string>();
    link = _model->GetLink(link_name_);
  }

  if (!link)
  {
    ROS_FATAL("gazebo_ros_baro plugin error: bodyName: %s does not exist\n", link_name_.c_str());
    return;
  }

  // default parameters
  namespace_.clear();
  frame_id_ = link->GetName();
  height_topic_ = "pressure_height";
  altimeter_topic_ = "altimeter";
  elevation_ = DEFAULT_ELEVATION;
  qnh_ = DEFAULT_QNH;

  // load parameters
  if (_sdf->HasElement("robotNamespace"))     namespace_ = _sdf->GetElement("robotNamespace")->Get<std::string>();
  if (_sdf->HasElement("frameId"))            frame_id_ = _sdf->GetElement("frameId")->Get<std::string>();
  if (_sdf->HasElement("topicName"))          height_topic_ = _sdf->GetElement("topicName")->Get<std::string>();
  if (_sdf->HasElement("altimeterTopicName")) altimeter_topic_ = _sdf->GetElement("altimeterTopicName")->Get<std::string>();
  if (_sdf->HasElement("elevation"))          elevation_ = _sdf->GetElement("elevation")->Get<double>();
  if (_sdf->HasElement("qnh"))                qnh_ = _sdf->GetElement("qnh")->Get<double>();

  // load sensor model
  sensor_model_.Load(_sdf);

  height_.header.frame_id = frame_id_;

  // Make sure the ROS node for Gazebo has already been initialized
  if (!ros::isInitialized())
  {
    ROS_FATAL_STREAM("A ROS node for Gazebo has not been initialized, unable to load plugin. "
      << "Load the Gazebo system plugin 'libgazebo_ros_api_plugin.so' in the gazebo_ros package)");
    return;
  }

  node_handle_ = new ros::NodeHandle(namespace_);

  // advertise height
  if (!height_topic_.empty()) {
    height_publisher_ = node_handle_->advertise<geometry_msgs::PointStamped>(height_topic_, 10);
  }

  // advertise altimeter
  if (!altimeter_topic_.empty()) {
    altimeter_publisher_ = node_handle_->advertise<hector_uav_msgs::Altimeter>(altimeter_topic_, 10);
  }

  // setup dynamic_reconfigure server
  dynamic_reconfigure_server_.reset(new dynamic_reconfigure::Server<SensorModelConfig>(ros::NodeHandle(*node_handle_, altimeter_topic_)));
  dynamic_reconfigure_server_->setCallback(boost::bind(&SensorModel::dynamicReconfigureCallback, &sensor_model_, _1, _2));

  Reset();

  // connect Update function
  updateTimer.Load(world, _sdf);
  updateConnection = updateTimer.Connect(boost::bind(&GazeboRosBaro::Update, this));
}

void GazeboRosBaro::Reset()
{
  updateTimer.Reset();
  sensor_model_.reset();
}

////////////////////////////////////////////////////////////////////////////////
// Update the controller
void GazeboRosBaro::Update()
{
  common::Time sim_time = world->SimTime();
  double dt = updateTimer.getTimeSinceLastUpdate().Double();

  ignition::math::Pose3 pose = link->WorldPose();
  double height = sensor_model_(pose.Pos()[2], dt);

  if (height_publisher_) {
    height_.header.stamp = ros::Time(sim_time.sec, sim_time.nsec);
    height_.point.z = height;
    height_publisher_.publish(height_);
  }

  if (altimeter_publisher_) {
    altimeter_.header = height_.header;
    altimeter_.altitude = height + elevation_;
    altimeter_.pressure = pow((1.0 - altimeter_.altitude / 44330.0), 5.263157) * qnh_;
    altimeter_.qnh = qnh_;
    altimeter_publisher_.publish(altimeter_);
  }
}

// Register this plugin with the simulator
GZ_REGISTER_MODEL_PLUGIN(GazeboRosBaro)

} // namespace gazebo
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_quadrotor/hector_quadrotor_gazebo_plugins/src/gazebo_quadrotor_simple_controller.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_quadrotor_gazebo_plugins/gazebo_quadrotor_simple_controller.h>
#include <gazebo/common/Events.hh>
#include <gazebo/physics/physics.hh>
#include <ignition/math.hh>
#include <cmath>

#include <geometry_msgs/Wrench.h>

namespace gazebo {

GazeboQuadrotorSimpleController::GazeboQuadrotorSimpleController()
{
}

////////////////////////////////////////////////////////////////////////////////
// Destructor
GazeboQuadrotorSimpleController::~GazeboQuadrotorSimpleController()
{
  //event::Events::DisconnectWorldUpdateBegin(updateConnection);
  this->updateConnection.reset();

  node_handle_->shutdown();
  delete node_handle_;
}

////////////////////////////////////////////////////////////////////////////////
// Load the controller
void GazeboQuadrotorSimpleController::Load(physics::ModelPtr _model, sdf::ElementPtr _sdf)
{
  gzwarn << "The GazeboQuadrotorSimpleController plugin is DEPRECATED in ROS hydro. Please use the twist_controller in package hector_quadrotor_controller instead." << std::endl;

  world = _model->GetWorld();
  link = _model->GetLink();
  link_name_ = link->GetName();

  // default parameters
  namespace_.clear();
  velocity_topic_ = "cmd_vel";
  imu_topic_.clear();
  state_topic_.clear();
  wrench_topic_ = "wrench_out";
  max_force_ = -1;
  auto_engage_ = true;

  // load parameters from sdf
  if (_sdf->HasElement("robotNamespace")) namespace_ = _sdf->GetElement("robotNamespace")->Get<std::string>();
  if (_sdf->HasElement("topicName"))      velocity_topic_ = _sdf->GetElement("topicName")->Get<std::string>();
  if (_sdf->HasElement("imuTopic"))       imu_topic_ = _sdf->GetElement("imuTopic")->Get<std::string>();
  if (_sdf->HasElement("stateTopic"))     state_topic_ = _sdf->GetElement("stateTopic")->Get<std::string>();
  if (_sdf->HasElement("wrenchTopic"))    wrench_topic_ = _sdf->GetElement("wrenchTopic")->Get<std::string>();
  if (_sdf->HasElement("maxForce"))       max_force_ = _sdf->GetElement("maxForce")->Get<double>();
  if (_sdf->HasElement("autoEngage"))     auto_engage_ = _sdf->GetElement("autoEngage")->Get<bool>();

  if (_sdf->HasElement("bodyName") && _sdf->GetElement("bodyName")->GetValue()) {
    link_name_ = _sdf->GetElement("bodyName")->Get<std::string>();
    link = _model->GetLink(link_name_);
  }

  if (!link)
  {
    ROS_FATAL("gazebo_ros_baro plugin error: bodyName: %s does not exist\n", link_name_.c_str());
    return;
  }

  // configure controllers
  controllers_.roll.Load(_sdf, "rollpitch");
  controllers_.pitch.Load(_sdf, "rollpitch");
  controllers_.yaw.Load(_sdf, "yaw");
  controllers_.velocity_x.Load(_sdf, "velocityXY");
  controllers_.velocity_y.Load(_sdf, "velocityXY");
  controllers_.velocity_z.Load(_sdf, "velocityZ");

  // Get inertia and mass of quadrotor body
  inertia = link->GetInertial()->PrincipalMoments();
  mass = link->GetInertial()->Mass();

  // Make sure the ROS node for Gazebo has already been initialized
  if (!ros::isInitialized())
  {
    ROS_FATAL_STREAM("A ROS node for Gazebo has not been initialized, unable to load plugin. "
      << "Load the Gazebo system plugin 'libgazebo_ros_api_plugin.so' in the gazebo_ros package)");
    return;
  }

  node_handle_ = new ros::NodeHandle(namespace_);
  ros::NodeHandle param_handle(*node_handle_, "controller");

  // subscribe command
  param_handle.getParam("velocity_topic", velocity_topic_);
  if (!velocity_topic_.empty())
  {
    ros::SubscribeOptions ops = ros::SubscribeOptions::create<geometry_msgs::Twist>(
      velocity_topic_, 1,
      boost::bind(&GazeboQuadrotorSimpleController::VelocityCallback, this, _1),
      ros::VoidPtr(), &callback_queue_);
    velocity_subscriber_ = node_handle_->subscribe(ops);
  }

  // subscribe imu
  param_handle.getParam("imu_topic", imu_topic_);
  if (!imu_topic_.empty())
  {
    ros::SubscribeOptions ops = ros::SubscribeOptions::create<sensor_msgs::Imu>(
      imu_topic_, 1,
      boost::bind(&GazeboQuadrotorSimpleController::ImuCallback, this, _1),
      ros::VoidPtr(), &callback_queue_);
    imu_subscriber_ = node_handle_->subscribe(ops);

    ROS_INFO_NAMED("quadrotor_simple_controller", "Using imu information on topic %s as source of orientation and angular velocity.", imu_topic_.c_str());
  }

  // subscribe state
  param_handle.getParam("state_topic", state_topic_);
  if (!state_topic_.empty())
  {
    ros::SubscribeOptions ops = ros::SubscribeOptions::create<nav_msgs::Odometry>(
      state_topic_, 1,
      boost::bind(&GazeboQuadrotorSimpleController::StateCallback, this, _1),
      ros::VoidPtr(), &callback_queue_);
    state_subscriber_ = node_handle_->subscribe(ops);

    ROS_INFO_NAMED("quadrotor_simple_controller", "Using state information on topic %s as source of state information.", state_topic_.c_str());
  }

  // advertise wrench
  param_handle.getParam("wrench_topic", wrench_topic_);
  if (!wrench_topic_.empty())
  {
    ros::AdvertiseOptions ops = ros::AdvertiseOptions::create<geometry_msgs::Wrench>(
      wrench_topic_, 10,
      ros::SubscriberStatusCallback(), ros::SubscriberStatusCallback(), ros::VoidConstPtr(), &callback_queue_
    );
    wrench_publisher_ = node_handle_->advertise(ops);
  }

  // engage/shutdown service servers
  {
    ros::AdvertiseServiceOptions ops = ros::AdvertiseServiceOptions::create<std_srvs::Empty>(
      "engage", boost::bind(&GazeboQuadrotorSimpleController::EngageCallback, this, _1, _2),
      ros::VoidConstPtr(), &callback_queue_
    );
    engage_service_server_ = node_handle_->advertiseService(ops);

    ops = ros::AdvertiseServiceOptions::create<std_srvs::Empty>(
      "shutdown", boost::bind(&GazeboQuadrotorSimpleController::ShutdownCallback, this, _1, _2),
      ros::VoidConstPtr(), &callback_queue_
    );
    shutdown_service_server_ = node_handle_->advertiseService(ops);
  }

  // callback_queue_thread_ = boost::thread( boost::bind( &GazeboQuadrotorSimpleController::CallbackQueueThread,this ) );


  Reset();

  // New Mechanism for Updating every World Cycle
  // Listen to the update event. This event is broadcast every
  // simulation iteration.
  controlTimer.Load(world, _sdf);
  updateConnection = event::Events::ConnectWorldUpdateBegin(
      boost::bind(&GazeboQuadrotorSimpleController::Update, this));
}

////////////////////////////////////////////////////////////////////////////////
// Callbacks
void GazeboQuadrotorSimpleController::VelocityCallback(const geometry_msgs::TwistConstPtr& velocity)
{
  velocity_command_ = *velocity;
}

void GazeboQuadrotorSimpleController::ImuCallback(const sensor_msgs::ImuConstPtr& imu)
{
  pose.Rot().Set(imu->orientation.w, imu->orientation.x, imu->orientation.y, imu->orientation.z);
  euler = pose.Rot().Euler();
  angular_velocity = pose.Rot().RotateVector(ignition::math::Vector3(imu->angular_velocity.x, imu->angular_velocity.y, imu->angular_velocity.z));
}

void GazeboQuadrotorSimpleController::StateCallback(const nav_msgs::OdometryConstPtr& state)
{
	ignition::math::Vector3 velocity1(velocity);

  if (imu_topic_.empty()) {
    pose.Pos().Set(state->pose.pose.position.x, state->pose.pose.position.y, state->pose.pose.position.z);
    pose.Rot().Set(state->pose.pose.orientation.w, state->pose.pose.orientation.x, state->pose.pose.orientation.y, state->pose.pose.orientation.z);
    euler = pose.Rot().Euler();
    angular_velocity.Set(state->twist.twist.angular.x, state->twist.twist.angular.y, state->twist.twist.angular.z);
  }

  velocity.Set(state->twist.twist.linear.x, state->twist.twist.linear.y, state->twist.twist.linear.z);

  // calculate acceleration
  double dt = !state_stamp.isZero() ? (state->header.stamp - state_stamp).toSec() : 0.0;
  state_stamp = state->header.stamp;
  if (dt > 0.0) {
    acceleration = (velocity - velocity1) / dt;
  } else {
    acceleration.Set();
  }
}

bool GazeboQuadrotorSimpleController::EngageCallback(std_srvs::Empty::Request &, std_srvs::Empty::Response &)
{
  ROS_INFO_NAMED("quadrotor_simple_controller", "Engaging motors!");
  running_ = true;
  return true;
}

bool GazeboQuadrotorSimpleController::ShutdownCallback(std_srvs::Empty::Request &, std_srvs::Empty::Response &)
{
  ROS_INFO_NAMED("quadrotor_simple_controller", "Shutting down motors!");
  running_ = false;
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// Update the controller
void GazeboQuadrotorSimpleController::Update()
{
  // Get new commands/state
  callback_queue_.callAvailable();

  double dt;
  if (controlTimer.update(dt) && dt > 0.0) {
    // Get Pose/Orientation from Gazebo (if no state subscriber is active)
    if (imu_topic_.empty()) {
      pose = link->WorldPose();
      angular_velocity = link->WorldAngularVel();
      euler = pose.Rot().Euler();
    }
    if (state_topic_.empty()) {
      acceleration = (link->WorldLinearVel() - velocity) / dt;
      velocity = link->WorldLinearVel();
    }

    // Auto engage/shutdown
    if (auto_engage_) {
      if (!running_ && velocity_command_.linear.z > 0.1) {
        running_ = true;
        ROS_INFO_NAMED("quadrotor_simple_controller", "Engaging motors!");
      } else if (running_ && controllers_.velocity_z.i < -1.0 && velocity_command_.linear.z < -0.1 && (velocity[2] > -0.1 && velocity[2] < 0.1)) {
        running_ = false;
        ROS_INFO_NAMED("quadrotor_simple_controller", "Shutting down motors!");
      }
    }

  //  static Time lastDebug;
  //  if ((world->GetSimTime() - lastDebug).Double() > 0.5) {
  //    ROS_DEBUG_STREAM_NAMED("quadrotor_simple_controller", "Velocity:         gazebo = [" << link->GetWorldLinearVel()   << "], state = [" << velocity << "]");
  //    ROS_DEBUG_STREAM_NAMED("quadrotor_simple_controller", "Acceleration:     gazebo = [" << link->GetWorldLinearAccel() << "], state = [" << acceleration << "]");
  //    ROS_DEBUG_STREAM_NAMED("quadrotor_simple_controller", "Angular Velocity: gazebo = [" << link->GetWorldAngularVel() << "], state = [" << angular_velocity << "]");
  //    lastDebug = world->GetSimTime();
  //  }

    // Get gravity
    ignition::math::Vector3 gravity_body = pose.Rot().RotateVector(world->Gravity());
    double gravity = gravity_body.Length();
    double load_factor = gravity * gravity / world->Gravity().Dot(gravity_body);  // Get gravity

    // Rotate vectors to coordinate frames relevant for control
    ignition::math::Quaterniond heading_quaternion(cos(euler.Z()/2),0,0,sin(euler.Z()/2));
    ignition::math::Vector3d velocity_xy = heading_quaternion.RotateVectorReverse(velocity);
    ignition::math::Vector3d acceleration_xy = heading_quaternion.RotateVectorReverse(acceleration);
    ignition::math::Vector3d angular_velocity_body = pose.Rot().RotateVectorReverse(angular_velocity);

    // update controllers
    force.Set(0.0, 0.0, 0.0);
    torque.Set(0.0, 0.0, 0.0);
    if (running_) {
      double pitch_command =  controllers_.velocity_x.update(velocity_command_.linear.x, velocity_xy.X(), acceleration_xy.X(), dt) / gravity;
      double roll_command  = -controllers_.velocity_y.update(velocity_command_.linear.y, velocity_xy.Y(), acceleration_xy.Y(), dt) / gravity;
      torque.X() = inertia.X() *  controllers_.roll.update(roll_command, euler.X(), angular_velocity_body.X(), dt);
      torque.Y() = inertia.Y() *  controllers_.pitch.update(pitch_command, euler.Y(), angular_velocity_body.Y(), dt);
      // torque.x = inertia.x *  controllers_.roll.update(-velocity_command_.linear.y/gravity, euler.x, angular_velocity_body.x, dt);
      // torque.y = inertia.y *  controllers_.pitch.update(velocity_command_.linear.x/gravity, euler.y, angular_velocity_body.y, dt);
      torque.Z() = inertia.Z() *  controllers_.yaw.update(velocity_command_.angular.z, angular_velocity.Z(), 0, dt);
      force.Z()  = mass      * (controllers_.velocity_z.update(velocity_command_.linear.z,  velocity.Z(), acceleration.Z(), dt) + load_factor * gravity);
      if (max_force_ > 0.0 && force.Z() > max_force_) force.Z() = max_force_;
      if (force.Z() < 0.0) force.Z() = 0.0;

    } else {
      controllers_.roll.reset();
      controllers_.pitch.reset();
      controllers_.yaw.reset();
      controllers_.velocity_x.reset();
      controllers_.velocity_y.reset();
      controllers_.velocity_z.reset();
    }

  //  static double lastDebugOutput = 0.0;
  //  if (last_time.Double() - lastDebugOutput > 0.1) {
  //    ROS_DEBUG_NAMED("quadrotor_simple_controller", "Velocity = [%g %g %g], Acceleration = [%g %g %g]", velocity.x, velocity.y, velocity.z, acceleration.x, acceleration.y, acceleration.z);
  //    ROS_DEBUG_NAMED("quadrotor_simple_controller", "Command: linear = [%g %g %g], angular = [%g %g %g], roll/pitch = [%g %g]", velocity_command_.linear.x, velocity_command_.linear.y, velocity_command_.linear.z, velocity_command_.angular.x*180/M_PI, velocity_command_.angular.y*180/M_PI, velocity_command_.angular.z*180/M_PI, roll_command*180/M_PI, pitch_command*180/M_PI);
  //    ROS_DEBUG_NAMED("quadrotor_simple_controller", "Mass: %g kg, Inertia: [%g %g %g], Load: %g g", mass, inertia.x, inertia.y, inertia.z, load_factor);
  //    ROS_DEBUG_NAMED("quadrotor_simple_controller", "Force: [%g %g %g], Torque: [%g %g %g]", force.x, force.y, force.z, torque.x, torque.y, torque.z);
  //    lastDebugOutput = last_time.Double();
  //  }

    // Publish wrench
    if (wrench_publisher_) {
      geometry_msgs::Wrench wrench;
      wrench.force.x = force.X();
      wrench.force.y = force.Y();
      wrench.force.z = force.Z();
      wrench.torque.x = torque.X();
      wrench.torque.y = torque.Y();
      wrench.torque.z = torque.Z();
      wrench_publisher_.publish(wrench);
    }
  }

  // set force and torque in gazebo
  link->AddRelativeForce(force);
  link->AddRelativeTorque(torque - link->GetInertial()->CoG().Cross(force));
}

////////////////////////////////////////////////////////////////////////////////
// Reset the controller
void GazeboQuadrotorSimpleController::Reset()
{
  controllers_.roll.reset();
  controllers_.pitch.reset();
  controllers_.yaw.reset();
  controllers_.velocity_x.reset();
  controllers_.velocity_y.reset();
  controllers_.velocity_z.reset();

  force.Set();
  torque.Set();

  // reset state
  pose.Reset();
  velocity.Set();
  angular_velocity.Set();
  acceleration.Set();
  euler.Set();
  state_stamp = ros::Time();

  running_ = false;
}

////////////////////////////////////////////////////////////////////////////////
// PID controller implementation
GazeboQuadrotorSimpleController::PIDController::PIDController()
{
}

GazeboQuadrotorSimpleController::PIDController::~PIDController()
{
}

void GazeboQuadrotorSimpleController::PIDController::Load(sdf::ElementPtr _sdf, const std::string& prefix)
{
  gain_p = 0.0;
  gain_d = 0.0;
  gain_i = 0.0;
  time_constant = 0.0;
  limit = -1.0;

  if (!_sdf) return;
  // _sdf->PrintDescription(_sdf->GetName());
  if (_sdf->HasElement(prefix + "ProportionalGain")) gain_p = _sdf->GetElement(prefix + "ProportionalGain")->Get<double>();
  if (_sdf->HasElement(prefix + "DifferentialGain")) gain_d = _sdf->GetElement(prefix + "DifferentialGain")->Get<double>();
  if (_sdf->HasElement(prefix + "IntegralGain"))     gain_i = _sdf->GetElement(prefix + "IntegralGain")->Get<double>();
  if (_sdf->HasElement(prefix + "TimeConstant"))     time_constant = _sdf->GetElement(prefix + "TimeConstant")->Get<double>();
  if (_sdf->HasElement(prefix + "Limit"))            limit = _sdf->GetElement(prefix + "Limit")->Get<double>();
}

double GazeboQuadrotorSimpleController::PIDController::update(double new_input, double x, double dx, double dt)
{
  // limit command
  if (limit > 0.0 && fabs(new_input) > limit) new_input = (new_input < 0 ? -1.0 : 1.0) * limit;

  // filter command
  if (dt + time_constant > 0.0) {
    dinput = (new_input - input) / (dt + time_constant);
    input  = (dt * new_input + time_constant * input) / (dt + time_constant);
  }

  // update proportional, differential and integral errors
  p = input - x;
  d = dinput - dx;
  i = i + dt * p;

  // update control output
  output = gain_p * p + gain_d * d + gain_i * i;
  return output;
}

void GazeboQuadrotorSimpleController::PIDController::reset()
{
  input = dinput = 0;
  p = i = d = output = 0;
}

// Register this plugin with the simulator
GZ_REGISTER_MODEL_PLUGIN(GazeboQuadrotorSimpleController)

} // namespace gazebo
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_quadrotor/hector_quadrotor_gazebo_plugins/src/gazebo_quadrotor_aerodynamics.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_quadrotor_gazebo_plugins/gazebo_quadrotor_aerodynamics.h>
#include <hector_quadrotor_model/helpers.h>

#include <gazebo/common/Events.hh>
#include <gazebo/physics/physics.hh>

#include <geometry_msgs/WrenchStamped.h>
#include <ignition/math.hh>
namespace gazebo {

using namespace common;
using namespace ignition::math;
using namespace hector_quadrotor_model;

GazeboQuadrotorAerodynamics::GazeboQuadrotorAerodynamics()
  : node_handle_(0)
{
}

GazeboQuadrotorAerodynamics::~GazeboQuadrotorAerodynamics()
{
  //event::Events::DisconnectWorldUpdateBegin(updateConnection);
  updateConnection.reset();
  if (node_handle_) {
    node_handle_->shutdown();
    if (callback_queue_thread_.joinable())
      callback_queue_thread_.join();
    delete node_handle_;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Load the controller
void GazeboQuadrotorAerodynamics::Load(physics::ModelPtr _model, sdf::ElementPtr _sdf)
{
  world = _model->GetWorld();
  link = _model->GetLink();

  // default parameters
  namespace_.clear();
  param_namespace_ = "quadrotor_aerodynamics";
  wind_topic_ = "/wind";
  wrench_topic_ = "aerodynamics/wrench";
  frame_id_ = link->GetName();

  // load parameters
  if (_sdf->HasElement("robotNamespace")) namespace_ = _sdf->GetElement("robotNamespace")->Get<std::string>();
  if (_sdf->HasElement("paramNamespace")) param_namespace_ = _sdf->GetElement("paramNamespace")->Get<std::string>();
  if (_sdf->HasElement("windTopicName"))  wind_topic_ = _sdf->GetElement("windTopicName")->Get<std::string>();
  if (_sdf->HasElement("wrenchTopic"))    wrench_topic_ = _sdf->GetElement("wrenchTopic")->Get<std::string>();
  if (_sdf->HasElement("frameId"))    frame_id_= _sdf->GetElement("frameId")->Get<std::string>();

  // Make sure the ROS node for Gazebo has already been initialized
  if (!ros::isInitialized())
  {
    ROS_FATAL_STREAM("A ROS node for Gazebo has not been initialized, unable to load plugin. "
      << "Load the Gazebo system plugin 'libgazebo_ros_api_plugin.so' in the gazebo_ros package)");
    return;
  }

  node_handle_ = new ros::NodeHandle(namespace_);

  // get model parameters
  if (!model_.configure(ros::NodeHandle(*node_handle_, param_namespace_))) {
    gzwarn << "[quadrotor_propulsion] Could not properly configure the aerodynamics plugin. Make sure you loaded the parameter file." << std::endl;
    return;
  }

  // subscribe command
  if (!wind_topic_.empty())
  {
    ros::SubscribeOptions ops;
    ops.callback_queue = &callback_queue_;
    ops.initByFullCallbackType<const geometry_msgs::Vector3 &>(
      wind_topic_, 1,
      boost::bind(&QuadrotorAerodynamics::setWind, &model_, _1)
    );
    wind_subscriber_ = node_handle_->subscribe(ops);
  }

  // advertise wrench
  if (!wrench_topic_.empty())
  {
    ros::AdvertiseOptions ops;
    ops.callback_queue = &callback_queue_;
    ops.init<geometry_msgs::WrenchStamped>(wrench_topic_, 10);
    wrench_publisher_ = node_handle_->advertise(ops);
  }

  // callback_queue_thread_ = boost::thread( boost::bind( &GazeboQuadrotorAerodynamics::QueueThread,this ) );

  // New Mechanism for Updating every World Cycle
  // Listen to the update event. This event is broadcast every
  // simulation iteration.
  updateConnection = event::Events::ConnectWorldUpdateBegin(
      boost::bind(&GazeboQuadrotorAerodynamics::Update, this));
}

////////////////////////////////////////////////////////////////////////////////
// Update the controller
void GazeboQuadrotorAerodynamics::Update()
{
  // Get simulator time
  Time current_time = world->SimTime();
  Time dt = current_time - last_time_;
  last_time_ = current_time;
  if (dt <= 0.0) return;

  // Get new commands/state
  callback_queue_.callAvailable();

  // fill input vector u for drag model
  geometry_msgs::Quaternion orientation;
  fromQuaternion(link->WorldPose().Rot(), orientation);
  model_.setOrientation(orientation);

  geometry_msgs::Twist twist;
  fromVector(link->WorldLinearVel(), twist.linear);
  fromVector(link->WorldAngularVel(), twist.angular);
  model_.setTwist(twist);

  // update the model
  model_.update(dt.Double());

  // get wrench from model
  ignition::math::Vector3d force, torque;
  toVector(model_.getWrench().force, force);
  toVector(model_.getWrench().torque, torque);

  // publish wrench
  if (wrench_publisher_) {
    geometry_msgs::WrenchStamped wrench_msg;
    wrench_msg.header.stamp = ros::Time(current_time.sec, current_time.nsec);
    wrench_msg.header.frame_id = frame_id_;
    wrench_msg.wrench = model_.getWrench();
    wrench_publisher_.publish(wrench_msg);
  }

  // set force and torque in gazebo
  link->AddRelativeForce(force);
  link->AddRelativeTorque(torque - link->GetInertial()->CoG().Cross(force));
}

////////////////////////////////////////////////////////////////////////////////
// Reset the controller
void GazeboQuadrotorAerodynamics::Reset()
{
  model_.reset();
}

////////////////////////////////////////////////////////////////////////////////
// custom callback queue thread
void GazeboQuadrotorAerodynamics::QueueThread()
{
  static const double timeout = 0.01;

  while (node_handle_->ok())
  {
    callback_queue_.callAvailable(ros::WallDuration(timeout));
  }
}

// Register this plugin with the simulator
GZ_REGISTER_MODEL_PLUGIN(GazeboQuadrotorAerodynamics)

} // namespace gazebo
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_quadrotor/hector_quadrotor_gazebo_plugins/src/gazebo_quadrotor_propulsion.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_quadrotor_gazebo_plugins/gazebo_quadrotor_propulsion.h>
#include <hector_quadrotor_model/helpers.h>

#include <gazebo/common/Events.hh>
#include <gazebo/physics/physics.hh>

#include <rosgraph_msgs/Clock.h>
#include <geometry_msgs/WrenchStamped.h>
#include <ignition/math.hh>

namespace gazebo {

using namespace common;
using namespace ignition::math;
using namespace hector_quadrotor_model;

GazeboQuadrotorPropulsion::GazeboQuadrotorPropulsion()
  : node_handle_(0)
{
}

GazeboQuadrotorPropulsion::~GazeboQuadrotorPropulsion()
{
  this->updateConnection.reset();
  if (node_handle_) {
    node_handle_->shutdown();
    if (callback_queue_thread_.joinable())
      callback_queue_thread_.join();
    delete node_handle_;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Load the controller
void GazeboQuadrotorPropulsion::Load(physics::ModelPtr _model, sdf::ElementPtr _sdf)
{
  world = _model->GetWorld();
  link = _model->GetLink();

  // default parameters
  namespace_.clear();
  param_namespace_ = "quadrotor_propulsion";
  trigger_topic_ = "quadro/trigger";
  command_topic_ = "command/motor";
  pwm_topic_ = "motor_pwm";
  wrench_topic_ = "propulsion/wrench";
  supply_topic_ = "supply";
  status_topic_ = "motor_status";
  control_tolerance_ = ros::Duration();
  control_delay_ = ros::Duration();
  frame_id_ = link->GetName();

  // load parameters
  if (_sdf->HasElement("robotNamespace"))   namespace_ = _sdf->GetElement("robotNamespace")->Get<std::string>();
  if (_sdf->HasElement("paramNamespace"))   param_namespace_ = _sdf->GetElement("paramNamespace")->Get<std::string>();
  if (_sdf->HasElement("triggerTopic"))     trigger_topic_ = _sdf->GetElement("triggerTopic")->Get<std::string>();
  if (_sdf->HasElement("topicName"))        command_topic_ = _sdf->GetElement("topicName")->Get<std::string>();
  if (_sdf->HasElement("pwmTopicName"))     pwm_topic_ = _sdf->GetElement("pwmTopicName")->Get<std::string>();
  if (_sdf->HasElement("wrenchTopic"))      wrench_topic_ = _sdf->GetElement("wrenchTopic")->Get<std::string>();
  if (_sdf->HasElement("supplyTopic"))      supply_topic_ = _sdf->GetElement("supplyTopic")->Get<std::string>();
  if (_sdf->HasElement("statusTopic"))      status_topic_ = _sdf->GetElement("statusTopic")->Get<std::string>();
  if (_sdf->HasElement("frameId"))          frame_id_= _sdf->GetElement("frameId")->Get<std::string>();

  if (_sdf->HasElement("voltageTopicName")) {
    gzwarn << "[quadrotor_propulsion] Plugin parameter 'voltageTopicName' is deprecated! Plese change your config to use "
           << "'topicName' for MotorCommand messages or 'pwmTopicName' for MotorPWM messages." << std::endl;
    pwm_topic_ = _sdf->GetElement("voltageTopicName")->Get<std::string>();
  }

  // set control timing parameters
  controlTimer.Load(world, _sdf, "control");
  motorStatusTimer.Load(world, _sdf, "motorStatus");
  if (_sdf->HasElement("controlTolerance")) control_tolerance_.fromSec(_sdf->GetElement("controlTolerance")->Get<double>());
  if (_sdf->HasElement("controlDelay"))     control_delay_.fromSec(_sdf->GetElement("controlDelay")->Get<double>());

  // set initial supply voltage
  if (_sdf->HasElement("supplyVoltage"))    model_.setInitialSupplyVoltage(_sdf->GetElement("supplyVoltage")->Get<double>());

  // Make sure the ROS node for Gazebo has already been initialized
  if (!ros::isInitialized())
  {
    ROS_FATAL_STREAM("A ROS node for Gazebo has not been initialized, unable to load plugin. "
      << "Load the Gazebo system plugin 'libgazebo_ros_api_plugin.so' in the gazebo_ros package)");
    return;
  }

  node_handle_ = new ros::NodeHandle(namespace_);

  // get model parameters
  if (!model_.configure(ros::NodeHandle(*node_handle_, param_namespace_))) {
    gzwarn << "[quadrotor_propulsion] Could not properly configure the propulsion plugin. Make sure you loaded the parameter file." << std::endl;
    return;
  }

  // publish trigger
  if (!trigger_topic_.empty())
  {
    ros::AdvertiseOptions ops;
    ops.callback_queue = &callback_queue_;
    ops.init<rosgraph_msgs::Clock>(trigger_topic_, 10);
    trigger_publisher_ = node_handle_->advertise(ops);
  }

  // subscribe voltage command
  if (!command_topic_.empty())
  {
    ros::SubscribeOptions ops;
    ops.callback_queue = &callback_queue_;
    ops.init<hector_uav_msgs::MotorCommand>(command_topic_, 1, boost::bind(&QuadrotorPropulsion::addCommandToQueue, &model_, _1));
    command_subscriber_ = node_handle_->subscribe(ops);
  }

  // subscribe pwm command
  if (!pwm_topic_.empty())
  {
    ros::SubscribeOptions ops;
    ops.callback_queue = &callback_queue_;
    ops.init<hector_uav_msgs::MotorPWM>(pwm_topic_, 1, boost::bind(&QuadrotorPropulsion::addPWMToQueue, &model_, _1));
    pwm_subscriber_ = node_handle_->subscribe(ops);
  }

  // advertise wrench
  if (!wrench_topic_.empty())
  {
    ros::AdvertiseOptions ops;
    ops.callback_queue = &callback_queue_;
    ops.init<geometry_msgs::WrenchStamped>(wrench_topic_, 10);
    wrench_publisher_ = node_handle_->advertise(ops);
  }

  // advertise and latch supply
  if (!supply_topic_.empty())
  {
    ros::AdvertiseOptions ops;
    ops.callback_queue = &callback_queue_;
    ops.latch = true;
    ops.init<hector_uav_msgs::Supply>(supply_topic_, 10);
    supply_publisher_ = node_handle_->advertise(ops);
    supply_publisher_.publish(model_.getSupply());
  }

  // advertise motor_status
  if (!status_topic_.empty())
  {
    ros::AdvertiseOptions ops;
    ops.callback_queue = &callback_queue_;
    ops.init<hector_uav_msgs::MotorStatus>(status_topic_, 10);
    motor_status_publisher_ = node_handle_->advertise(ops);
  }

  // callback_queue_thread_ = boost::thread( boost::bind( &GazeboQuadrotorPropulsion::QueueThread,this ) );

  Reset();

  // New Mechanism for Updating every World Cycle
  // Listen to the update event. This event is broadcast every
  // simulation iteration.
  updateConnection = event::Events::ConnectWorldUpdateBegin(
      boost::bind(&GazeboQuadrotorPropulsion::Update, this));
}

////////////////////////////////////////////////////////////////////////////////
// Update the controller
void GazeboQuadrotorPropulsion::Update()
{
  // Get simulator time
  Time current_time = world->SimTime();
  Time dt = current_time - last_time_;
  last_time_ = current_time;
  if (dt <= 0.0) return;

  // Send trigger
  bool trigger = controlTimer.getUpdatePeriod() > 0.0 ? controlTimer.update() : false;
  if (trigger && trigger_publisher_) {
    rosgraph_msgs::Clock clock;
    clock.clock = ros::Time(current_time.sec, current_time.nsec);
    trigger_publisher_.publish(clock);

    ROS_DEBUG_STREAM_NAMED("quadrotor_propulsion", "Sent a trigger message at t = " << current_time.Double() << " (dt = " << (current_time - last_trigger_time_).Double() << ")");
    last_trigger_time_ = current_time;
  }

  // Get new commands/state
  callback_queue_.callAvailable();

  // Process input queue
  model_.processQueue(ros::Time(current_time.sec, current_time.nsec), control_tolerance_, control_delay_, (model_.getMotorStatus().on && trigger) ? ros::WallDuration(1.0) : ros::WallDuration(), &callback_queue_);

  // fill input vector u for propulsion model
  geometry_msgs::Twist twist;
  fromVector(link->RelativeLinearVel(), twist.linear);
  fromVector(link->RelativeAngularVel(), twist.angular);
  model_.setTwist(twist);

  // update the model
  model_.update(dt.Double());

  // get wrench from model
  ignition::math::Vector3d force, torque;
  toVector(model_.getWrench().force, force);
  toVector(model_.getWrench().torque, torque);

  // publish wrench
  if (wrench_publisher_) {
    geometry_msgs::WrenchStamped wrench_msg;
    wrench_msg.header.stamp = ros::Time(current_time.sec, current_time.nsec);
    wrench_msg.header.frame_id = frame_id_;
    wrench_msg.wrench = model_.getWrench();
    wrench_publisher_.publish(wrench_msg);
  }

  // publish motor status
  if (motor_status_publisher_ && motorStatusTimer.update() /* && current_time >= last_motor_status_time_ + control_period_ */) {
    hector_uav_msgs::MotorStatus motor_status = model_.getMotorStatus();
    motor_status.header.stamp = ros::Time(current_time.sec, current_time.nsec);
    motor_status_publisher_.publish(motor_status);
    last_motor_status_time_ = current_time;
  }

  // publish supply
  if (supply_publisher_ && current_time >= last_supply_time_ + 1.0) {
    supply_publisher_.publish(model_.getSupply());
    last_supply_time_ = current_time;
  }

  // set force and torque in gazebo
  link->AddRelativeForce(force);
  link->AddRelativeTorque(torque - link->GetInertial()->CoG().Cross(force));
}

////////////////////////////////////////////////////////////////////////////////
// Reset the controller
void GazeboQuadrotorPropulsion::Reset()
{
  model_.reset();
  last_time_ = Time();
  last_trigger_time_ = Time();
  last_motor_status_time_ = Time();
  last_supply_time_ = Time();
}

////////////////////////////////////////////////////////////////////////////////
// custom callback queue thread
void GazeboQuadrotorPropulsion::QueueThread()
{
  static const double timeout = 0.01;

  while (node_handle_->ok())
  {
    callback_queue_.callAvailable(ros::WallDuration(timeout));
  }
}

// Register this plugin with the simulator
GZ_REGISTER_MODEL_PLUGIN(GazeboQuadrotorPropulsion)

} // namespace gazebo
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_quadrotor/hector_quadrotor_controller_gazebo/src/quadrotor_hardware_gazebo.cpp
//=================================================================================================
// Copyright (c) 2013, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_quadrotor_controller/quadrotor_hardware_gazebo.h>

#include <geometry_msgs/WrenchStamped.h>
#include <ignition/math.hh>
namespace hector_quadrotor_controller_gazebo {

QuadrotorHardwareSim::QuadrotorHardwareSim()
{
  this->registerInterface(static_cast<QuadrotorInterface *>(this));

  wrench_output_ = addInput<WrenchCommandHandle>("wrench");
  motor_output_ = addInput<MotorCommandHandle>("motor");
}

QuadrotorHardwareSim::~QuadrotorHardwareSim()
{

}

bool QuadrotorHardwareSim::initSim(
    const std::string& robot_namespace,
    ros::NodeHandle model_nh,
    gazebo::physics::ModelPtr parent_model,
    const urdf::Model *const urdf_model,
    std::vector<transmission_interface::TransmissionInfo> transmissions)
{
  ros::NodeHandle param_nh(model_nh, "controller");

  // store parent model pointer
  model_ = parent_model;
  link_ = model_->GetLink();
  ignition::math::Vector3d physics_test = model_->GetWorld()->Gravity();
  //physics_ = model_->GetWorld()->Gravity();
  

  model_nh.param<std::string>("world_frame", world_frame_, "world");
  model_nh.param<std::string>("base_link_frame", base_link_frame_, "base_link");

  // subscribe state
  std::string state_topic;
  param_nh.getParam("state_topic", state_topic);
  if (!state_topic.empty()) {
    ros::SubscribeOptions ops = ros::SubscribeOptions::create<nav_msgs::Odometry>(state_topic, 1, boost::bind(&QuadrotorHardwareSim::stateCallback, this, _1), ros::VoidConstPtr(), &callback_queue_);
    subscriber_state_ = model_nh.subscribe(ops);

    gzlog << "[hector_quadrotor_controller_gazebo] Using topic '" << subscriber_state_.getTopic() << "' as state input for control" << std::endl;
  } else {
    gzlog << "[hector_quadrotor_controller_gazebo] Using ground truth from Gazebo as state input for control" << std::endl;
  }

  // subscribe imu
  std::string imu_topic;
  param_nh.getParam("imu_topic", imu_topic);
  if (!imu_topic.empty()) {
    ros::SubscribeOptions ops = ros::SubscribeOptions::create<sensor_msgs::Imu>(imu_topic, 1, boost::bind(&QuadrotorHardwareSim::imuCallback, this, _1), ros::VoidConstPtr(), &callback_queue_);
    subscriber_imu_ = model_nh.subscribe(ops);
    gzlog << "[hector_quadrotor_controller_gazebo] Using topic '" << subscriber_imu_.getTopic() << "' as imu input for control" << std::endl;
  } else {
    gzlog << "[hector_quadrotor_controller_gazebo] Using ground truth from Gazebo as imu input for control" << std::endl;
  }

  // subscribe motor_status
  {
    ros::SubscribeOptions ops = ros::SubscribeOptions::create<hector_uav_msgs::MotorStatus>("motor_status", 1, boost::bind(&QuadrotorHardwareSim::motorStatusCallback, this, _1), ros::VoidConstPtr(), &callback_queue_);
    subscriber_motor_status_ = model_nh.subscribe(ops);
  }

  // publish wrench
  {
    ros::AdvertiseOptions ops = ros::AdvertiseOptions::create<geometry_msgs::WrenchStamped>("command/wrench", 1, ros::SubscriberStatusCallback(), ros::SubscriberStatusCallback(), ros::VoidConstPtr(), &callback_queue_);
    publisher_wrench_command_ = model_nh.advertise(ops);
  }

  // publish motor command
  {
    ros::AdvertiseOptions ops = ros::AdvertiseOptions::create<hector_uav_msgs::MotorCommand>("command/motor", 1, ros::SubscriberStatusCallback(), ros::SubscriberStatusCallback(), ros::VoidConstPtr(), &callback_queue_);
    publisher_motor_command_ = model_nh.advertise(ops);
  }

  return true;
}

bool QuadrotorHardwareSim::getMassAndInertia(double &mass, double inertia[3]) {
  if (!link_) return false;
  mass = link_->GetInertial()->Mass();
  ignition::math::Vector3d Inertia = link_->GetInertial()->PrincipalMoments();
  inertia[0] = Inertia.X();
  inertia[1] = Inertia.Y();
  inertia[2] = Inertia.Z();
  return true;
}

void QuadrotorHardwareSim::stateCallback(const nav_msgs::OdometryConstPtr &state) {
  // calculate acceleration
  if (!header_.stamp.isZero() && !state->header.stamp.isZero()) {
    const double acceleration_time_constant = 0.1;
    double dt((state->header.stamp - header_.stamp).toSec());
    if (dt > 0.0) {
      acceleration_.x = ((state->twist.twist.linear.x - twist_.linear.x) + acceleration_time_constant * acceleration_.x) / (dt + acceleration_time_constant);
      acceleration_.y = ((state->twist.twist.linear.y - twist_.linear.y) + acceleration_time_constant * acceleration_.y) / (dt + acceleration_time_constant);
      acceleration_.z = ((state->twist.twist.linear.z - twist_.linear.z) + acceleration_time_constant * acceleration_.z) / (dt + acceleration_time_constant);
    }
  }

  header_ = state->header;
  pose_ = state->pose.pose;
  twist_ = state->twist.twist;
}

void QuadrotorHardwareSim::imuCallback(const sensor_msgs::ImuConstPtr &imu) {
  imu_ = *imu;
}

void QuadrotorHardwareSim::motorStatusCallback(const hector_uav_msgs::MotorStatusConstPtr &motor_status) {
  motor_status_ = *motor_status;
}

void QuadrotorHardwareSim::readSim(ros::Time time, ros::Duration period)
{
  // call all available subscriber callbacks now
  callback_queue_.callAvailable();

  // read state from Gazebo
  const double acceleration_time_constant = 0.1;
  gz_pose_             =  link_->WorldPose();
  gz_acceleration_     = ((link_->WorldLinearVel() - gz_velocity_) + acceleration_time_constant * gz_acceleration_) / (period.toSec() + acceleration_time_constant);
  gz_velocity_         =  link_->WorldLinearVel();
  gz_angular_velocity_ =  link_->WorldAngularVel();

  if (!subscriber_state_) {
    header_.frame_id = world_frame_;
    header_.stamp = time;
    pose_.position.x = gz_pose_.Pos().X();
    pose_.position.y = gz_pose_.Pos().Y();
    pose_.position.z = gz_pose_.Pos().Z();
    pose_.orientation.w = gz_pose_.Rot().W();
    pose_.orientation.x = gz_pose_.Rot().X();
    pose_.orientation.y = gz_pose_.Rot().Y();
    pose_.orientation.z = gz_pose_.Rot().Z();
    twist_.linear.x = gz_velocity_.X();
    twist_.linear.y = gz_velocity_.Y();
    twist_.linear.z = gz_velocity_.Z();
    twist_.angular.x = gz_angular_velocity_.X();
    twist_.angular.y = gz_angular_velocity_.Y();
    twist_.angular.z = gz_angular_velocity_.Z();
    acceleration_.x = gz_acceleration_.X();
    acceleration_.y = gz_acceleration_.Y();
    acceleration_.z = gz_acceleration_.Z();
  }

  if (!subscriber_imu_) {
    imu_.orientation.w = gz_pose_.Rot().W();
    imu_.orientation.x = gz_pose_.Rot().X();
    imu_.orientation.y = gz_pose_.Rot().Y();
    imu_.orientation.z = gz_pose_.Rot().Z();

    ignition::math::Vector3d gz_angular_velocity_body = gz_pose_.Rot().RotateVectorReverse(gz_angular_velocity_);
    imu_.angular_velocity.x = gz_angular_velocity_body.X();
    imu_.angular_velocity.y = gz_angular_velocity_body.Y();
    imu_.angular_velocity.z = gz_angular_velocity_body.Z();

        
    ignition::math::Vector3d gz_linear_acceleration_body = gz_pose_.Rot().RotateVectorReverse(gz_acceleration_ - physics_test);
    imu_.linear_acceleration.x = gz_linear_acceleration_body.X();
    imu_.linear_acceleration.y = gz_linear_acceleration_body.Y();
    imu_.linear_acceleration.z = gz_linear_acceleration_body.Z();
  }
}

void QuadrotorHardwareSim::writeSim(ros::Time time, ros::Duration period)
{
  bool result_written = false;

  if (motor_output_->connected() && motor_output_->enabled()) {
    publisher_motor_command_.publish(motor_output_->getCommand());
    result_written = true;
  }

  if (wrench_output_->connected() && wrench_output_->enabled()) {
    geometry_msgs::WrenchStamped wrench;
    wrench.header.stamp = time;
    wrench.header.frame_id = base_link_frame_;
    wrench.wrench = wrench_output_->getCommand();
    publisher_wrench_command_.publish(wrench);

    if (!result_written) {
      ignition::math::Vector3d force(wrench.wrench.force.x, wrench.wrench.force.y, wrench.wrench.force.z);
      ignition::math::Vector3d torque(wrench.wrench.torque.x, wrench.wrench.torque.y, wrench.wrench.torque.z);
      link_->AddRelativeForce(force);
      link_->AddRelativeTorque(torque - link_->GetInertial()->CoG().Cross(force));
    }
  }
}

} // namespace hector_quadrotor_controller_gazebo

#include <pluginlib/class_list_macros.h>
PLUGINLIB_EXPORT_CLASS(hector_quadrotor_controller_gazebo::QuadrotorHardwareSim, gazebo_ros_control::RobotHWSim)
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_quadrotor/hector_quadrotor_model/src/quadrotor_propulsion.cpp
//=================================================================================================
// Copyright (c) 2013, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_quadrotor_model/quadrotor_propulsion.h>
#include <hector_quadrotor_model/helpers.h>

#include <ros/node_handle.h>
#include <ros/callback_queue.h>

#include <boost/array.hpp>

#include "matlab_helpers.h"

//extern "C" {
//  #include "quadrotorPropulsion/quadrotorPropulsion.h"
//  #include "quadrotorPropulsion/quadrotorPropulsion_initialize.h"
//}

namespace hector_quadrotor_model {

struct PropulsionParameters
{
    real_T k_m;
    real_T k_t;
    real_T CT2s;
    real_T CT1s;
    real_T CT0s;
    real_T Psi;
    real_T J_M;
    real_T R_A;
    real_T alpha_m;
    real_T beta_m;
    real_T l_m;

    PropulsionParameters()
      : k_m(0.0)
      , k_t(0.0)
      , CT2s(0.0)
      , CT1s(0.0)
      , CT0s(0.0)
      , Psi(0.0)
      , J_M(std::numeric_limits<real_T>::infinity())
      , R_A(std::numeric_limits<real_T>::infinity())
      , alpha_m(0.0)
      , beta_m(0.0)
      , l_m(0.0)
    {}
};

struct QuadrotorPropulsion::PropulsionModel {
  PropulsionParameters parameters_;
  boost::array<real_T,4>  x;
  boost::array<real_T,4>  x_pred;
  boost::array<real_T,10> u;
  boost::array<real_T,14> y;
};

QuadrotorPropulsion::QuadrotorPropulsion()
{
  // initialize propulsion model
  // quadrotorPropulsion_initialize();
  propulsion_model_ = new PropulsionModel;
}

QuadrotorPropulsion::~QuadrotorPropulsion()
{
  delete propulsion_model_;
}

/*
 * quadrotorPropulsion.c
 *
 * Code generation for function 'quadrotorPropulsion'
 *
 * C source code generated on: Sun Nov  3 13:34:35 2013
 *
 */

/* Include files */
//#include "rt_nonfinite.h"
//#include "motorspeed.h"
//#include "quadrotorPropulsion.h"
//#include "quadrotorPropulsion_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void quadrotorPropulsion(const real_T xin[4], const real_T uin[10], const
  PropulsionParameters parameter, real_T dt, real_T y[14], real_T xpred[4])
{
  real_T v_1[4];
  int32_T i;
  real_T F_m[4];
  real_T U[4];
  real_T M_e[4];
  real_T I[4];
  real_T F[3];
  real_T b_F_m;
  real_T temp;
  real_T d0;

  /*  initialize vectors */
  for (i = 0; i < 4; i++) {
    xpred[i] = xin[i];

    /*  motorspeed */
    v_1[i] = 0.0;
  }

  memset(&y[0], 0, 14U * sizeof(real_T));
  for (i = 0; i < 4; i++) {
    F_m[i] = 0.0;
    U[i] = 0.0;
    M_e[i] = 0.0;
    I[i] = 0.0;
  }

  for (i = 0; i < 3; i++) {
    F[i] = 0.0;
  }

  /*  Input variables */
  U[0] = uin[6];
  U[1] = uin[7];
  U[2] = uin[8];
  U[3] = uin[9];

  /*  Constants */
  v_1[0] = -uin[2] + parameter.l_m * uin[4];
  v_1[1] = -uin[2] - parameter.l_m * uin[3];
  v_1[2] = -uin[2] - parameter.l_m * uin[4];
  v_1[3] = -uin[2] + parameter.l_m * uin[3];

  /*  calculate thrust for all 4 rotors */
  for (i = 0; i < 4; i++) {
    /*  if the flow speed at infinity is negative */
    if (v_1[i] < 0.0) {
      b_F_m = (parameter.CT2s * rt_powd_snf(v_1[i], 2.0) + parameter.CT1s *
               v_1[i] * xin[i]) + parameter.CT0s * rt_powd_snf(xin[i], 2.0);

      /*  if the flow speed at infinity is positive */
    } else {
      b_F_m = (-parameter.CT2s * rt_powd_snf(v_1[i], 2.0) + parameter.CT1s *
               v_1[i] * xin[i]) + parameter.CT0s * rt_powd_snf(xin[i], 2.0);
    }

    /*  sum up all rotor forces */
    /*  Identification of Roxxy2827-34 motor with 10x4.5 propeller */
    /*  temporarily used Expressions */
    temp = (U[i] * parameter.beta_m - parameter.Psi * xin[i]) / (2.0 *
      parameter.R_A);
    temp += sqrt(rt_powd_snf(temp, 2.0) + U[i] * parameter.alpha_m /
                 parameter.R_A);
    d0 = parameter.Psi * temp;

    /*  electrical torque motor 1-4 */
    /*  new version */
    /*  old version */
    /*  fx   = (Psi/R_A*(U-Psi*omega_m) - M_m)/J_M; */
    /*  A    = -(Psi^2/R_A)/J_M; */
    /*  B(1) =  Psi/(J_M*R_A); */
    /*  B(2) = -1/J_M; */
    /*  system outputs. Use euler solver to predict next time step */
    /*  predicted motor speed */
    /*  electric torque */
    /* y = [M_e I]; */
    /*  system jacobian */
    /*  A       = 1 + dt*A; */
    /*  input jacobian */
    /*  B       = A*B*dt; */
    M_e[i] = d0;
    I[i] = temp;
    xpred[i] = xin[i] + dt * (1.0 / parameter.J_M * (d0 - (parameter.k_t * b_F_m
      + parameter.k_m * xin[i])));
    F_m[i] = b_F_m;
    F[2] += b_F_m;
  }

  /*  System output, i.e. force and torque of quadrotor */
  y[0] = 0.0;
  y[1] = 0.0;
  y[2] = F[2];

  /*  torque for rotating quadrocopter around x-axis is the mechanical torque */
  y[3] = (F_m[3] - F_m[1]) * parameter.l_m;

  /*  torque for rotating quadrocopter around y-axis is the mechanical torque */
  y[4] = (F_m[0] - F_m[2]) * parameter.l_m;

  /*  torque for rotating quadrocopter around z-axis is the electrical torque */
  y[5] = ((-M_e[0] - M_e[2]) + M_e[1]) + M_e[3];

  /*  motor speeds (rad/s) */
  for (i = 0; i < 4; i++) {
    y[i + 6] = xpred[i];
  }

  /*  motor current (A) */
  for (i = 0; i < 4; i++) {
    y[i + 10] = I[i];
  }

  /*  M_e(1:4) / Psi; */
}

/* End of code generation (quadrotorPropulsion.c) */

inline void QuadrotorPropulsion::f(const double xin[4], const double uin[10], double dt, double y[14], double xpred[4]) const
{
  quadrotorPropulsion(xin, uin, propulsion_model_->parameters_, dt, y, xpred);
}

bool QuadrotorPropulsion::configure(const ros::NodeHandle &param)
{
  // get model parameters
  if (!param.getParam("k_m",     propulsion_model_->parameters_.k_m)) return false;
  if (!param.getParam("k_t",     propulsion_model_->parameters_.k_t)) return false;
  if (!param.getParam("CT0s",    propulsion_model_->parameters_.CT0s)) return false;
  if (!param.getParam("CT1s",    propulsion_model_->parameters_.CT1s)) return false;
  if (!param.getParam("CT2s",    propulsion_model_->parameters_.CT2s)) return false;
  if (!param.getParam("J_M",     propulsion_model_->parameters_.J_M)) return false;
  if (!param.getParam("l_m",     propulsion_model_->parameters_.l_m)) return false;
  if (!param.getParam("Psi",     propulsion_model_->parameters_.Psi)) return false;
  if (!param.getParam("R_A",     propulsion_model_->parameters_.R_A)) return false;
  if (!param.getParam("alpha_m", propulsion_model_->parameters_.alpha_m)) return false;
  if (!param.getParam("beta_m",  propulsion_model_->parameters_.beta_m)) return false;

#ifndef NDEBUG
  std::cout << "Loaded the following quadrotor propulsion model parameters from namespace " << param.getNamespace() << ":" << std::endl;
  std::cout << "k_m     = " << propulsion_model_->parameters_.k_m << std::endl;
  std::cout << "k_t     = " << propulsion_model_->parameters_.k_t << std::endl;
  std::cout << "CT2s    = " << propulsion_model_->parameters_.CT2s << std::endl;
  std::cout << "CT1s    = " << propulsion_model_->parameters_.CT1s << std::endl;
  std::cout << "CT0s    = " << propulsion_model_->parameters_.CT0s << std::endl;
  std::cout << "Psi     = " << propulsion_model_->parameters_.Psi << std::endl;
  std::cout << "J_M     = " << propulsion_model_->parameters_.J_M << std::endl;
  std::cout << "R_A     = " << propulsion_model_->parameters_.R_A << std::endl;
  std::cout << "l_m     = " << propulsion_model_->parameters_.l_m << std::endl;
  std::cout << "alpha_m = " << propulsion_model_->parameters_.alpha_m << std::endl;
  std::cout << "beta_m  = " << propulsion_model_->parameters_.beta_m << std::endl;
#endif

  initial_voltage_ = 14.8;
  param.getParam("supply_voltage", initial_voltage_);
  reset();

  return true;
}

void QuadrotorPropulsion::reset()
{
  boost::mutex::scoped_lock lock(mutex_);

  propulsion_model_->x.assign(0.0);
  propulsion_model_->x_pred.assign(0.0);
  propulsion_model_->u.assign(0.0);
  propulsion_model_->y.assign(0.0);

  wrench_ = geometry_msgs::Wrench();

  motor_status_ = hector_uav_msgs::MotorStatus();
  motor_status_.voltage.resize(4);
  motor_status_.frequency.resize(4);
  motor_status_.current.resize(4);

  supply_ = hector_uav_msgs::Supply();
  supply_.voltage.resize(1);
  supply_.voltage[0] = initial_voltage_;

  last_command_time_ = ros::Time();

  command_queue_ = std::queue<hector_uav_msgs::MotorPWMConstPtr>(); // .clear();
}

void QuadrotorPropulsion::setVoltage(const hector_uav_msgs::MotorPWM& voltage)
{
  boost::mutex::scoped_lock lock(mutex_);
  last_command_time_ = voltage.header.stamp;

  if (motor_status_.on && voltage.pwm.size() >= 4) {
    propulsion_model_->u[6] = voltage.pwm[0] * supply_.voltage[0] / 255.0;
    propulsion_model_->u[7] = voltage.pwm[1] * supply_.voltage[0] / 255.0;
    propulsion_model_->u[8] = voltage.pwm[2] * supply_.voltage[0] / 255.0;
    propulsion_model_->u[9] = voltage.pwm[3] * supply_.voltage[0] / 255.0;
  } else {
    propulsion_model_->u[6] = 0.0;
    propulsion_model_->u[7] = 0.0;
    propulsion_model_->u[8] = 0.0;
    propulsion_model_->u[9] = 0.0;
  }
}

void QuadrotorPropulsion::setTwist(const geometry_msgs::Twist &twist)
{
  boost::mutex::scoped_lock lock(mutex_);
  propulsion_model_->u[0] = twist.linear.x;
  propulsion_model_->u[1] = -twist.linear.y;
  propulsion_model_->u[2] = -twist.linear.z;
  propulsion_model_->u[3] = twist.angular.x;
  propulsion_model_->u[4] = -twist.angular.y;
  propulsion_model_->u[5] = -twist.angular.z;

  // We limit the input velocities to +-100. Required for numeric stability in case of collisions,
  // where velocities returned by Gazebo can be very high.
  limit(boost::iterator_range<boost::array<real_T,10>::iterator>(&(propulsion_model_->u[0]), &(propulsion_model_->u[6])), -100.0, 100.0);
}

void QuadrotorPropulsion::addCommandToQueue(const hector_uav_msgs::MotorCommandConstPtr& command)
{
  hector_uav_msgs::MotorPWMPtr pwm(new hector_uav_msgs::MotorPWM);
  pwm->header = command->header;
  pwm->pwm.resize(command->voltage.size());
  for(std::size_t i = 0; i < command->voltage.size(); ++i) {
    int temp = round(command->voltage[i] / supply_.voltage[0] * 255.0);
    if (temp < 0)
      pwm->pwm[i] = 0;
    else if (temp > 255)
      pwm->pwm[i] = 255;
    else
      pwm->pwm[i] = temp;
  }
  addPWMToQueue(pwm);
}

void QuadrotorPropulsion::addPWMToQueue(const hector_uav_msgs::MotorPWMConstPtr& pwm)
{
  boost::mutex::scoped_lock lock(command_queue_mutex_);

  if (!motor_status_.on) {
    ROS_WARN_NAMED("quadrotor_propulsion", "Received new motor command. Enabled motors.");
    engage();
  }

  ROS_DEBUG_STREAM_NAMED("quadrotor_propulsion", "Received motor command valid at " << pwm->header.stamp);
  command_queue_.push(pwm);
  command_condition_.notify_all();
}

bool QuadrotorPropulsion::processQueue(const ros::Time &timestamp, const ros::Duration &tolerance, const ros::Duration &delay, const ros::WallDuration &wait, ros::CallbackQueue *callback_queue) {
  boost::mutex::scoped_lock lock(command_queue_mutex_);
  bool new_command = false;

  ros::Time min(timestamp), max(timestamp);
  try {
    min = timestamp - delay - tolerance /* - ros::Duration(dt) */;
  } catch (std::runtime_error &e) {}

  try {
    max = timestamp - delay + tolerance;
  } catch (std::runtime_error &e) {}

  do {

    while(!command_queue_.empty()) {
      hector_uav_msgs::MotorPWMConstPtr new_motor_voltage = command_queue_.front();
      ros::Time new_time = new_motor_voltage->header.stamp;

      if (new_time.isZero() || (new_time >= min && new_time <= max)) {
        setVoltage(*new_motor_voltage);
        command_queue_.pop();
        new_command = true;

      // new motor command is too old:
      } else if (new_time < min) {
        ROS_DEBUG_NAMED("quadrotor_propulsion", "Command received was %fs too old. Discarding.", (new_time - timestamp).toSec());
        command_queue_.pop();

      // new motor command is too new:
      } else {
        break;
      }
    }

    if (command_queue_.empty() && !new_command) {
      if (!motor_status_.on || wait.isZero()) break;

      ROS_DEBUG_NAMED("quadrotor_propulsion", "Waiting for command at simulation step t = %fs... last update was %fs ago", timestamp.toSec(), (timestamp - last_command_time_).toSec());
      if (!callback_queue) {
        if (command_condition_.timed_wait(lock, wait.toBoost())) continue;
      } else {
        lock.unlock();
        callback_queue->callAvailable(wait);
        lock.lock();
        if (!command_queue_.empty()) continue;
      }

      ROS_ERROR_NAMED("quadrotor_propulsion", "Command timed out. Disabled motors.");
      shutdown();
    }

  } while(false);

  if (new_command) {
      ROS_DEBUG_STREAM_NAMED("quadrotor_propulsion", "Using motor command valid at t = " << last_command_time_.toSec() << "s for simulation step at t = " << timestamp.toSec() << "s (dt = " << (timestamp - last_command_time_).toSec() << "s)");
  }

  return new_command;
}

void QuadrotorPropulsion::update(double dt)
{
  if (dt <= 0.0) return;

  ROS_DEBUG_STREAM_NAMED("quadrotor_propulsion", "propulsion.twist:   " << PrintVector<double>(propulsion_model_->u.begin(), propulsion_model_->u.begin() + 6));
  ROS_DEBUG_STREAM_NAMED("quadrotor_propulsion", "propulsion.voltage: " << PrintVector<double>(propulsion_model_->u.begin() + 6, propulsion_model_->u.begin() + 10));
  ROS_DEBUG_STREAM_NAMED("quadrotor_propulsion", "propulsion.x_prior: " << PrintVector<double>(propulsion_model_->x.begin(), propulsion_model_->x.end()));

  checknan(propulsion_model_->u, "propulsion model input");
  checknan(propulsion_model_->x, "propulsion model state");

  // update propulsion model
  f(propulsion_model_->x.data(), propulsion_model_->u.data(), dt, propulsion_model_->y.data(), propulsion_model_->x_pred.data());
  propulsion_model_->x = propulsion_model_->x_pred;

  ROS_DEBUG_STREAM_NAMED("quadrotor_propulsion", "propulsion.x_post:  " << PrintVector<double>(propulsion_model_->x.begin(), propulsion_model_->x.end()));
  ROS_DEBUG_STREAM_NAMED("quadrotor_propulsion", "propulsion.force:   " << PrintVector<double>(propulsion_model_->y.begin() + 0, propulsion_model_->y.begin() + 3) << " " <<
                                                 "propulsion.torque:  " << PrintVector<double>(propulsion_model_->y.begin() + 3, propulsion_model_->y.begin() + 6));

  checknan(propulsion_model_->y, "propulsion model output");

  wrench_.force.x  =  propulsion_model_->y[0];
  wrench_.force.y  = -propulsion_model_->y[1];
  wrench_.force.z  =  propulsion_model_->y[2];
  wrench_.torque.x =  propulsion_model_->y[3];
  wrench_.torque.y = -propulsion_model_->y[4];
  wrench_.torque.z = -propulsion_model_->y[5];

  motor_status_.voltage[0] = propulsion_model_->u[6];
  motor_status_.voltage[1] = propulsion_model_->u[7];
  motor_status_.voltage[2] = propulsion_model_->u[8];
  motor_status_.voltage[3] = propulsion_model_->u[9];

  motor_status_.frequency[0] = propulsion_model_->y[6];
  motor_status_.frequency[1] = propulsion_model_->y[7];
  motor_status_.frequency[2] = propulsion_model_->y[8];
  motor_status_.frequency[3] = propulsion_model_->y[9];
  motor_status_.running = motor_status_.frequency[0] > 1.0 && motor_status_.frequency[1] > 1.0 && motor_status_.frequency[2] > 1.0 && motor_status_.frequency[3] > 1.0;

  motor_status_.current[0] = propulsion_model_->y[10];
  motor_status_.current[1] = propulsion_model_->y[11];
  motor_status_.current[2] = propulsion_model_->y[12];
  motor_status_.current[3] = propulsion_model_->y[13];
}

void QuadrotorPropulsion::engage()
{
  motor_status_.on = true;
}

void QuadrotorPropulsion::shutdown()
{
  motor_status_.on = false;
}

} // namespace hector_quadrotor_model
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_quadrotor/hector_quadrotor_model/src/quadrotor_aerodynamics.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_quadrotor_model/quadrotor_aerodynamics.h>
#include <hector_quadrotor_model/helpers.h>

#include <ros/node_handle.h>
#include <boost/array.hpp>

//extern "C" {
//  #include "quadrotorDrag/quadrotorDrag.h"
//  #include "quadrotorDrag/quadrotorDrag_initialize.h"
//}

#include <boost/array.hpp>
#include <Eigen/Geometry>

#include "matlab_helpers.h"

namespace hector_quadrotor_model {

struct DragParameters
{
    real_T C_wxy;
    real_T C_wz;
    real_T C_mxy;
    real_T C_mz;

    DragParameters()
      : C_wxy(0.0)
      , C_wz(0.0)
      , C_mxy(0.0)
      , C_mz(0.0)
    {}
};

// extern void quadrotorDrag(const real_T uin[6], const DragParameters parameter, real_T dt, real_T y[6]);
struct QuadrotorAerodynamics::DragModel {
  DragParameters parameters_;
  boost::array<real_T,6> u;
  boost::array<real_T,6> y;
};

QuadrotorAerodynamics::QuadrotorAerodynamics()
{
  // initialize drag model
  // quadrotorDrag_initialize();
  drag_model_ = new DragModel;
}

QuadrotorAerodynamics::~QuadrotorAerodynamics()
{
  delete drag_model_;
}

/*
 * quadrotorDrag.c
 *
 * Code generation for function 'quadrotorDrag'
 *
 * C source code generated on: Sun Nov  3 13:34:38 2013
 *
 */

/* Include files */
//#include "rt_nonfinite.h"
//#include "quadrotorDrag.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Definitions */
void quadrotorDrag(const real_T uin[6], const DragParameters parameter, real_T
                   dt, real_T y[6])
{
  int32_T i;
  real_T absoluteVelocity;
  real_T absoluteAngularVelocity;

  /*  initialize vectors */
  for (i = 0; i < 6; i++) {
    y[i] = 0.0;
  }

  /*  Input variables */
  /*  Constants */
  /*  temporarily used vector */
  absoluteVelocity = sqrt((rt_powd_snf(uin[0], 2.0) + rt_powd_snf(uin[1], 2.0))
    + rt_powd_snf(uin[2], 2.0));
  absoluteAngularVelocity = sqrt((rt_powd_snf(uin[3], 2.0) + rt_powd_snf(uin[4],
    2.0)) + rt_powd_snf(uin[5], 2.0));

  /*  system outputs */
  /*  calculate drag force */
  y[0] = parameter.C_wxy * absoluteVelocity * uin[0];
  y[1] = parameter.C_wxy * absoluteVelocity * uin[1];
  y[2] = parameter.C_wz * absoluteVelocity * uin[2];

  /*  calculate draq torque */
  y[3] = parameter.C_mxy * absoluteAngularVelocity * uin[3];
  y[4] = parameter.C_mxy * absoluteAngularVelocity * uin[4];
  y[5] = parameter.C_mz * absoluteAngularVelocity * uin[5];
}

/* End of code generation (quadrotorDrag.c) */

inline void QuadrotorAerodynamics::f(const double uin[10], double dt, double y[14]) const
{
  quadrotorDrag(uin, drag_model_->parameters_, dt, y);
}

bool QuadrotorAerodynamics::configure(const ros::NodeHandle &param)
{
  // get model parameters
  if (!param.getParam("C_wxy", drag_model_->parameters_.C_wxy)) return false;
  if (!param.getParam("C_wz",  drag_model_->parameters_.C_wz)) return false;
  if (!param.getParam("C_mxy", drag_model_->parameters_.C_mxy)) return false;
  if (!param.getParam("C_mz",  drag_model_->parameters_.C_mz)) return false;

#ifndef NDEBUG
  std::cout << "Loaded the following quadrotor drag model parameters from namespace " << param.getNamespace() << ":" << std::endl;
  std::cout << "C_wxy = " << drag_model_->parameters_.C_wxy << std::endl;
  std::cout << "C_wz = "  << drag_model_->parameters_.C_wz << std::endl;
  std::cout << "C_mxy = " << drag_model_->parameters_.C_mxy << std::endl;
  std::cout << "C_mz = "  << drag_model_->parameters_.C_mz << std::endl;
#endif

  reset();
  return true;
}

void QuadrotorAerodynamics::reset()
{
  boost::mutex::scoped_lock lock(mutex_);
  drag_model_->u.assign(0.0);
  drag_model_->y.assign(0.0);

  twist_ = geometry_msgs::Twist();
  wind_ = geometry_msgs::Vector3();
  wrench_ = geometry_msgs::Wrench();
}

void QuadrotorAerodynamics::setOrientation(const geometry_msgs::Quaternion& orientation)
{
  boost::mutex::scoped_lock lock(mutex_);
  orientation_ = orientation;
}

void QuadrotorAerodynamics::setTwist(const geometry_msgs::Twist& twist)
{
  boost::mutex::scoped_lock lock(mutex_);
  twist_ = twist;
}

void QuadrotorAerodynamics::setBodyTwist(const geometry_msgs::Twist& body_twist)
{
  boost::mutex::scoped_lock lock(mutex_);
  Eigen::Quaterniond orientation(orientation_.w, orientation_.x, orientation_.y, orientation_.z);
  Eigen::Matrix<double,3,3> inverse_rotation_matrix(orientation.inverse().toRotationMatrix());

  Eigen::Vector3d body_linear(body_twist.linear.x, body_twist.linear.y, body_twist.linear.z);
  Eigen::Vector3d world_linear(inverse_rotation_matrix * body_linear);
  twist_.linear.x = world_linear.x();
  twist_.linear.y = world_linear.y();
  twist_.linear.z = world_linear.z();

  Eigen::Vector3d body_angular(body_twist.angular.x, body_twist.angular.y, body_twist.angular.z);
  Eigen::Vector3d world_angular(inverse_rotation_matrix * body_angular);
  twist_.angular.x = world_angular.x();
  twist_.angular.y = world_angular.y();
  twist_.angular.z = world_angular.z();
}

void QuadrotorAerodynamics::setWind(const geometry_msgs::Vector3& wind)
{
  boost::mutex::scoped_lock lock(mutex_);
  wind_ = wind;
}

void QuadrotorAerodynamics::update(double dt)
{
  if (dt <= 0.0) return;
  boost::mutex::scoped_lock lock(mutex_);

  // fill input vector u for drag model
  drag_model_->u[0] =  (twist_.linear.x - wind_.x);
  drag_model_->u[1] = -(twist_.linear.y - wind_.y);
  drag_model_->u[2] = -(twist_.linear.z - wind_.z);
  drag_model_->u[3] =  twist_.angular.x;
  drag_model_->u[4] = -twist_.angular.y;
  drag_model_->u[5] = -twist_.angular.z;

  // We limit the input velocities to +-100. Required for numeric stability in case of collisions,
  // where velocities returned by Gazebo can be very high.
  limit(drag_model_->u, -100.0, 100.0);

  // convert input to body coordinates
  Eigen::Quaterniond orientation(orientation_.w, orientation_.x, orientation_.y, orientation_.z);
  Eigen::Matrix<double,3,3> rotation_matrix(orientation.toRotationMatrix());
  Eigen::Map<Eigen::Vector3d> linear(&(drag_model_->u[0]));
  Eigen::Map<Eigen::Vector3d> angular(&(drag_model_->u[3]));
  linear  = rotation_matrix * linear;
  angular = rotation_matrix * angular;

  ROS_DEBUG_STREAM_NAMED("quadrotor_aerodynamics", "aerodynamics.twist:  " << PrintVector<double>(drag_model_->u.begin(), drag_model_->u.begin() + 6));
  checknan(drag_model_->u, "drag model input");

  // update drag model
  f(drag_model_->u.data(), dt, drag_model_->y.data());

  ROS_DEBUG_STREAM_NAMED("quadrotor_aerodynamics", "aerodynamics.force:  " << PrintVector<double>(drag_model_->y.begin() + 0, drag_model_->y.begin() + 3));
  ROS_DEBUG_STREAM_NAMED("quadrotor_aerodynamics", "aerodynamics.torque: " << PrintVector<double>(drag_model_->y.begin() + 3, drag_model_->y.begin() + 6));
  checknan(drag_model_->y, "drag model output");

  // drag_model_ gives us inverted vectors!
  wrench_.force.x  = -( drag_model_->y[0]);
  wrench_.force.y  = -(-drag_model_->y[1]);
  wrench_.force.z  = -(-drag_model_->y[2]);
  wrench_.torque.x = -( drag_model_->y[3]);
  wrench_.torque.y = -(-drag_model_->y[4]);
  wrench_.torque.z = -(-drag_model_->y[5]);
}

} // namespace hector_quadrotor_model
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_quadrotor/hector_quadrotor_controller/src/motor_controller.cpp
//=================================================================================================
// Copyright (c) 2013, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_quadrotor_controller/quadrotor_interface.h>
#include <controller_interface/controller.h>

#include <geometry_msgs/WrenchStamped.h>
#include <std_srvs/Empty.h>

#include <ros/subscriber.h>
#include <ros/callback_queue.h>

namespace hector_quadrotor_controller {

using namespace controller_interface;

class MotorController : public controller_interface::Controller<QuadrotorInterface>
{
public:
  MotorController()
    : node_handle_(0)
  {}

  ~MotorController()
  {
    if (node_handle_) {
      node_handle_->shutdown();
      delete node_handle_;
      node_handle_ = 0;
    }
  }

  bool init(QuadrotorInterface *interface, ros::NodeHandle &root_nh, ros::NodeHandle &controller_nh)
  {
    // get interface handles
    wrench_input_  = interface->addInput<WrenchCommandHandle>("wrench");
    motor_output_  = interface->addOutput<MotorCommandHandle>("motor");

    // initialize NodeHandle
    delete node_handle_;
    node_handle_ = new ros::NodeHandle(root_nh);

    // load parameters
    controller_nh.getParam("force_per_voltage", parameters_.force_per_voltage = 0.559966216);
    controller_nh.getParam("torque_per_voltage", parameters_.torque_per_voltage = 7.98598e-3);
    controller_nh.getParam("lever", parameters_.lever = 0.275);
    root_nh.param<std::string>("base_link_frame", base_link_frame_, "base_link");

    // TODO: calculate these parameters from the quadrotor_propulsion parameters
//    quadrotor_propulsion:
//            k_t: 0.015336864714397
//            k_m: -7.011631909766668e-5

//            CT2s: -1.3077e-2
//            CT1s: -2.5224e-4
//            CT0s:  1.538190483976698e-5

    return true;
  }

  void reset()
  {
    wrench_.wrench.force.x  = 0.0;
    wrench_.wrench.force.y  = 0.0;
    wrench_.wrench.force.z  = 0.0;
    wrench_.wrench.torque.x = 0.0;
    wrench_.wrench.torque.y = 0.0;
    wrench_.wrench.torque.z = 0.0;

    motor_.force.assign(4, 0.0);
    motor_.torque.assign(4, 0.0);
    motor_.frequency.clear();
    motor_.voltage.assign(4, 0.0);
  }

  void wrenchCommandCallback(const geometry_msgs::WrenchStampedConstPtr& command)
  {
    wrench_ = *command;
    if (wrench_.header.stamp.isZero()) wrench_.header.stamp = ros::Time::now();

    // start controller if it not running
    if (!isRunning()) this->startRequest(wrench_.header.stamp);
  }

  void starting(const ros::Time &time)
  {
    reset();
    motor_output_->start();
  }

  void stopping(const ros::Time &time)
  {
    motor_output_->stop();
  }

  void update(const ros::Time& time, const ros::Duration& period)
  {
    // Get wrench command input
    if (wrench_input_->connected() && wrench_input_->enabled()) {
      wrench_.wrench = wrench_input_->getCommand();
    }

    // Update output
    if (wrench_.wrench.force.z > 0.0) {

      double nominal_thrust_per_motor = wrench_.wrench.force.z / 4.0;
      motor_.force[0] =  nominal_thrust_per_motor - wrench_.wrench.torque.y / 2.0 / parameters_.lever;
      motor_.force[1] =  nominal_thrust_per_motor - wrench_.wrench.torque.x / 2.0 / parameters_.lever;
      motor_.force[2] =  nominal_thrust_per_motor + wrench_.wrench.torque.y / 2.0 / parameters_.lever;
      motor_.force[3] =  nominal_thrust_per_motor + wrench_.wrench.torque.x / 2.0 / parameters_.lever;

      double nominal_torque_per_motor = wrench_.wrench.torque.z / 4.0;
      motor_.voltage[0] = motor_.force[0] / parameters_.force_per_voltage + nominal_torque_per_motor / parameters_.torque_per_voltage;
      motor_.voltage[1] = motor_.force[1] / parameters_.force_per_voltage - nominal_torque_per_motor / parameters_.torque_per_voltage;
      motor_.voltage[2] = motor_.force[2] / parameters_.force_per_voltage + nominal_torque_per_motor / parameters_.torque_per_voltage;
      motor_.voltage[3] = motor_.force[3] / parameters_.force_per_voltage - nominal_torque_per_motor / parameters_.torque_per_voltage;

      motor_.torque[0] = motor_.voltage[0] * parameters_.torque_per_voltage;
      motor_.torque[1] = motor_.voltage[1] * parameters_.torque_per_voltage;
      motor_.torque[2] = motor_.voltage[2] * parameters_.torque_per_voltage;
      motor_.torque[3] = motor_.voltage[3] * parameters_.torque_per_voltage;

      if (motor_.voltage[0] < 0.0) motor_.voltage[0] = 0.0;
      if (motor_.voltage[1] < 0.0) motor_.voltage[1] = 0.0;
      if (motor_.voltage[2] < 0.0) motor_.voltage[2] = 0.0;
      if (motor_.voltage[3] < 0.0) motor_.voltage[3] = 0.0;

    } else {
      reset();
    }

    // set wrench output
    motor_.header.stamp = time;
    motor_.header.frame_id = "base_link";
    motor_output_->setCommand(motor_);
  }

private:
  WrenchCommandHandlePtr wrench_input_;
  MotorCommandHandlePtr motor_output_;

  ros::NodeHandle *node_handle_;
  ros::Subscriber wrench_subscriber_;
  ros::ServiceServer engage_service_server_;
  ros::ServiceServer shutdown_service_server_;

  geometry_msgs::WrenchStamped wrench_;
  hector_uav_msgs::MotorCommand motor_;
  std::string base_link_frame_;

  struct {
    double force_per_voltage;     // coefficient for linearized volts to force conversion for a single motor [N / V]
    double torque_per_voltage;    // coefficient for linearized volts to force conversion for a single motor [Nm / V]
    double lever;                 // the lever arm from origin to the motor axes (symmetry assumption) [m]
  } parameters_;
};

} // namespace hector_quadrotor_controller

#include <pluginlib/class_list_macros.h>
PLUGINLIB_EXPORT_CLASS(hector_quadrotor_controller::MotorController, controller_interface::ControllerBase)
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_quadrotor/hector_quadrotor_controller/src/pid.cpp
//=================================================================================================
// Copyright (c) 2013, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_quadrotor_controller/pid.h>
#include <limits>

namespace hector_quadrotor_controller {

PID::state::state()
  : p(std::numeric_limits<double>::quiet_NaN())
  , i(0.0)
  , d(std::numeric_limits<double>::quiet_NaN())
  , input(std::numeric_limits<double>::quiet_NaN())
  , dinput(0.0)
  , dx(std::numeric_limits<double>::quiet_NaN())
{
}

PID::parameters::parameters()
  : enabled(true)
  , time_constant(0.0)
  , k_p(0.0)
  , k_i(0.0)
  , k_d(0.0)
  , limit_i(std::numeric_limits<double>::quiet_NaN())
  , limit_output(std::numeric_limits<double>::quiet_NaN())
{
}

PID::PID()
{}

PID::PID(const parameters& params)
  : parameters_(params)
{}

PID::~PID()
{}

void PID::init(const ros::NodeHandle &param_nh)
{
  param_nh.getParam("enabled", parameters_.enabled);
  param_nh.getParam("k_p", parameters_.k_p);
  param_nh.getParam("k_i", parameters_.k_i);
  param_nh.getParam("k_d", parameters_.k_d);
  param_nh.getParam("limit_i", parameters_.limit_i);
  param_nh.getParam("limit_output", parameters_.limit_output);
  param_nh.getParam("time_constant", parameters_.time_constant);
}

void PID::reset()
{
  state_ = state();
}

template <typename T> static inline T& checknan(T& value)
{
  if (std::isnan(value)) value = T();
  return value;
}

double PID::update(double input, double x, double dx, const ros::Duration& dt)
{
  if (!parameters_.enabled) return 0.0;
  double dt_sec = dt.toSec();

  // low-pass filter input
  if (std::isnan(state_.input)) state_.input = input;
  if (dt_sec + parameters_.time_constant > 0.0) {
    state_.dinput = (input - state_.input) / (dt_sec + parameters_.time_constant);
    state_.input  = (dt_sec * input + parameters_.time_constant * state_.input) / (dt_sec + parameters_.time_constant);
  }

  return update(state_.input - x, dx, dt);
}

double PID::update(double error, double dx, const ros::Duration& dt)
{
  if (!parameters_.enabled) return 0.0;
  if (std::isnan(error)) return 0.0;
  double dt_sec = dt.toSec();

  // integral error
  state_.i += error * dt_sec;
  if (parameters_.limit_i > 0.0)
  {
    if (state_.i >  parameters_.limit_i) state_.i =  parameters_.limit_i;
    if (state_.i < -parameters_.limit_i) state_.i = -parameters_.limit_i;
  }

  // differential error
  if (dt_sec > 0.0 && !std::isnan(state_.p) && !std::isnan(state_.dx)) {
    state_.d = (error - state_.p) / dt_sec + state_.dx - dx;
  } else {
    state_.d = -dx;
  }
  state_.dx = dx;

  // proportional error
  state_.p = error;

  // calculate output...
  double output = parameters_.k_p * state_.p + parameters_.k_i * state_.i + parameters_.k_d * state_.d;
  int antiwindup = 0;
  if (parameters_.limit_output > 0.0)
  {
    if (output >  parameters_.limit_output) { output =  parameters_.limit_output; antiwindup =  1; }
    if (output < -parameters_.limit_output) { output = -parameters_.limit_output; antiwindup = -1; }
  }
  if (antiwindup && (error * dt_sec * antiwindup > 0.0)) state_.i -= error * dt_sec;

  checknan(output);
  return output;
}

double PID::getFilteredControlError(double& filtered_error, double time_constant, const ros::Duration& dt)
{
  double dt_sec = dt.toSec();
  filtered_error = checknan(filtered_error);
  if (dt_sec + time_constant > 0.0) {
    filtered_error = (time_constant * filtered_error + dt_sec * state_.p) / (dt_sec + time_constant);
  }
  return filtered_error;
}

} // namespace hector_quadrotor_controller

=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_quadrotor/hector_quadrotor_controller/src/pose_controller.cpp
//=================================================================================================
// Copyright (c) 2013, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_quadrotor_controller/quadrotor_interface.h>
#include <hector_quadrotor_controller/pid.h>

#include <controller_interface/controller.h>

#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>

#include <ros/subscriber.h>
#include <ros/callback_queue.h>

namespace hector_quadrotor_controller {

using namespace controller_interface;

class PoseController : public controller_interface::Controller<QuadrotorInterface>
{
public:
  PoseController() {}

  ~PoseController() {}

  bool init(QuadrotorInterface *interface, ros::NodeHandle &root_nh, ros::NodeHandle &controller_nh)
  {
    // get interface handles
    pose_         = interface->getPose();
    twist_        = interface->getTwist();

    pose_input_   = interface->addInput<PoseCommandHandle>("pose");
    twist_input_  = interface->addInput<TwistCommandHandle>("pose/twist");
    twist_limit_  = interface->addInput<TwistCommandHandle>("pose/twist_limit");
    twist_output_ = interface->addOutput<TwistCommandHandle>("twist");

    node_handle_ = root_nh;

    // subscribe to commanded pose and velocity
    pose_subscriber_ = node_handle_.subscribe<geometry_msgs::PoseStamped>("command/pose", 1, boost::bind(&PoseController::poseCommandCallback, this, _1));
    twist_subscriber_ = node_handle_.subscribe<geometry_msgs::TwistStamped>("command/twist", 1, boost::bind(&PoseController::twistCommandCallback, this, _1));

    // initialize PID controllers
    pid_.x.init(ros::NodeHandle(controller_nh, "xy"));
    pid_.y.init(ros::NodeHandle(controller_nh, "xy"));
    pid_.z.init(ros::NodeHandle(controller_nh, "z"));
    pid_.yaw.init(ros::NodeHandle(controller_nh, "yaw"));

    return true;
  }

  void reset()
  {
    pid_.x.reset();
    pid_.y.reset();
    pid_.z.reset();
    pid_.yaw.reset();
  }

  void poseCommandCallback(const geometry_msgs::PoseStampedConstPtr& command)
  {
    pose_command_ = *command;
    if (!(pose_input_->connected())) *pose_input_ = &(pose_command_.pose);
    pose_input_->start();

    ros::Time start_time = command->header.stamp;
    if (start_time.isZero()) start_time = ros::Time::now();
    if (!isRunning()) this->startRequest(start_time);
  }

  void twistCommandCallback(const geometry_msgs::TwistStampedConstPtr& command)
  {
    twist_command_ = *command;
    if (!(twist_input_->connected())) *twist_input_ = &(twist_command_.twist);
    twist_input_->start();

    ros::Time start_time = command->header.stamp;
    if (start_time.isZero()) start_time = ros::Time::now();
    if (!isRunning()) this->startRequest(start_time);
  }

  void starting(const ros::Time &time)
  {
    reset();
    twist_output_->start();
  }

  void stopping(const ros::Time &time)
  {
    twist_output_->stop();
  }

  void update(const ros::Time& time, const ros::Duration& period)
  {
    Twist output;

    // check command timeout
    // TODO

    // return if no pose command is available
    if (pose_input_->enabled()) {
      // control horizontal position
      double error_n, error_w;
      HorizontalPositionCommandHandle(*pose_input_).getError(*pose_, error_n, error_w);
      output.linear.x = pid_.x.update(error_n, twist_->twist().linear.x, period);
      output.linear.y = pid_.y.update(error_w, twist_->twist().linear.y, period);

      // control height
      output.linear.z = pid_.z.update(HeightCommandHandle(*pose_input_).getError(*pose_), twist_->twist().linear.z, period);

      // control yaw angle
      output.angular.z = pid_.yaw.update(HeadingCommandHandle(*pose_input_).getError(*pose_), twist_->twist().angular.z, period);
    }

    // add twist command if available
    if (twist_input_->enabled())
    {
      output.linear.x  += twist_input_->getCommand().linear.x;
      output.linear.y  += twist_input_->getCommand().linear.y;
      output.linear.z  += twist_input_->getCommand().linear.z;
      output.angular.x += twist_input_->getCommand().angular.x;
      output.angular.y += twist_input_->getCommand().angular.y;
      output.angular.z += twist_input_->getCommand().angular.z;
    }

    // limit twist
    if (twist_limit_->enabled())
    {
      double linear_xy = sqrt(output.linear.x*output.linear.x + output.linear.y*output.linear.y);
      double limit_linear_xy  = std::max(twist_limit_->get()->linear.x, twist_limit_->get()->linear.y);
      if (limit_linear_xy > 0.0 && linear_xy > limit_linear_xy) {
        output.linear.x *= limit_linear_xy / linear_xy;
        output.linear.y *= limit_linear_xy / linear_xy;
      }
      if (twist_limit_->get()->linear.z > 0.0 && fabs(output.linear.z) > twist_limit_->get()->linear.z) {
        output.linear.z *= twist_limit_->get()->linear.z / fabs(output.linear.z);
      }
      double angular_xy = sqrt(output.angular.x*output.angular.x + output.angular.y*output.angular.y);
      double limit_angular_xy  = std::max(twist_limit_->get()->angular.x, twist_limit_->get()->angular.y);
      if (limit_angular_xy > 0.0 && angular_xy > limit_angular_xy) {
        output.angular.x *= limit_angular_xy / angular_xy;
        output.angular.y *= limit_angular_xy / angular_xy;
      }
      if (twist_limit_->get()->angular.z > 0.0 && fabs(output.angular.z) > twist_limit_->get()->angular.z) {
        output.angular.z *= twist_limit_->get()->angular.z / fabs(output.angular.z);
      }
    }

    // set twist output
    twist_output_->setCommand(output);
  }

private:
  PoseHandlePtr pose_;
  PoseCommandHandlePtr pose_input_;
  TwistHandlePtr twist_;
  TwistCommandHandlePtr twist_input_;
  TwistCommandHandlePtr twist_limit_;
  TwistCommandHandlePtr twist_output_;

  geometry_msgs::PoseStamped pose_command_;
  geometry_msgs::TwistStamped twist_command_;

  ros::NodeHandle node_handle_;
  ros::Subscriber pose_subscriber_;
  ros::Subscriber twist_subscriber_;

  struct {
    PID x;
    PID y;
    PID z;
    PID yaw;
  } pid_;
};

} // namespace hector_quadrotor_controller

#include <pluginlib/class_list_macros.h>
PLUGINLIB_EXPORT_CLASS(hector_quadrotor_controller::PoseController, controller_interface::ControllerBase)
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_quadrotor/hector_quadrotor_controller/src/twist_controller.cpp
//=================================================================================================
// Copyright (c) 2013, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_quadrotor_controller/quadrotor_interface.h>
#include <hector_quadrotor_controller/pid.h>

#include <controller_interface/controller.h>

#include <geometry_msgs/TwistStamped.h>
#include <geometry_msgs/WrenchStamped.h>
#include <std_srvs/Empty.h>

#include <ros/subscriber.h>
#include <ros/callback_queue.h>

#include <boost/thread.hpp>

#include <limits>

namespace hector_quadrotor_controller {

using namespace controller_interface;

class TwistController : public controller_interface::Controller<QuadrotorInterface>
{
public:
  TwistController()
  {}

  ~TwistController()
  {}

  bool init(QuadrotorInterface *interface, ros::NodeHandle &root_nh, ros::NodeHandle &controller_nh)
  {
    // get interface handles
    pose_          = interface->getPose();
    twist_         = interface->getTwist();
    acceleration_  = interface->getAcceleration();
    twist_input_   = interface->addInput<TwistCommandHandle>("twist");
    wrench_output_ = interface->addOutput<WrenchCommandHandle>("wrench");
    node_handle_ = root_nh;

    // subscribe to commanded twist (geometry_msgs/TwistStamped) and cmd_vel (geometry_msgs/Twist)
    twist_subscriber_ = node_handle_.subscribe<geometry_msgs::TwistStamped>("command/twist", 1, boost::bind(&TwistController::twistCommandCallback, this, _1));
    cmd_vel_subscriber_ = node_handle_.subscribe<geometry_msgs::Twist>("cmd_vel", 1, boost::bind(&TwistController::cmd_velCommandCallback, this, _1));

    // engage/shutdown service servers
    engage_service_server_ = node_handle_.advertiseService<std_srvs::Empty::Request, std_srvs::Empty::Response>("engage", boost::bind(&TwistController::engageCallback, this, _1, _2));
    shutdown_service_server_ = node_handle_.advertiseService<std_srvs::Empty::Request, std_srvs::Empty::Response>("shutdown", boost::bind(&TwistController::shutdownCallback, this, _1, _2));

    // initialize PID controllers
    pid_.linear.x.init(ros::NodeHandle(controller_nh, "linear/xy"));
    pid_.linear.y.init(ros::NodeHandle(controller_nh, "linear/xy"));
    pid_.linear.z.init(ros::NodeHandle(controller_nh, "linear/z"));
    pid_.angular.x.init(ros::NodeHandle(controller_nh, "angular/xy"));
    pid_.angular.y.init(ros::NodeHandle(controller_nh, "angular/xy"));
    pid_.angular.z.init(ros::NodeHandle(controller_nh, "angular/z"));

    // load other parameters
    controller_nh.getParam("auto_engage", auto_engage_ = true);
    controller_nh.getParam("limits/load_factor", load_factor_limit = 1.5);
    controller_nh.getParam("limits/force/z", limits_.force.z);
    controller_nh.getParam("limits/torque/xy", limits_.torque.x);
    controller_nh.getParam("limits/torque/xy", limits_.torque.y);
    controller_nh.getParam("limits/torque/z", limits_.torque.z);
    root_nh.param<std::string>("base_link_frame", base_link_frame_, "base_link");

    // get mass and inertia from QuadrotorInterface
    interface->getMassAndInertia(mass_, inertia_);

    command_given_in_stabilized_frame_ = false;

    return true;
  }

  void reset()
  {
    pid_.linear.x.reset();
    pid_.linear.y.reset();
    pid_.linear.z.reset();
    pid_.angular.x.reset();
    pid_.angular.y.reset();
    pid_.angular.z.reset();

    wrench_.wrench.force.x  = 0.0;
    wrench_.wrench.force.y  = 0.0;
    wrench_.wrench.force.z  = 0.0;
    wrench_.wrench.torque.x = 0.0;
    wrench_.wrench.torque.y = 0.0;
    wrench_.wrench.torque.z = 0.0;

    linear_z_control_error_ = 0.0;
    motors_running_ = false;
  }

  void twistCommandCallback(const geometry_msgs::TwistStampedConstPtr& command)
  {
    boost::mutex::scoped_lock lock(command_mutex_);

    command_ = *command;
    if (command_.header.stamp.isZero()) command_.header.stamp = ros::Time::now();
    command_given_in_stabilized_frame_ = false;

    // start controller if it not running
    if (!isRunning()) this->startRequest(command_.header.stamp);
  }

  void cmd_velCommandCallback(const geometry_msgs::TwistConstPtr& command)
  {
    boost::mutex::scoped_lock lock(command_mutex_);

    command_.twist = *command;
    command_.header.stamp = ros::Time::now();
    command_given_in_stabilized_frame_ = true;

    // start controller if it not running
    if (!isRunning()) this->startRequest(command_.header.stamp);
  }

  bool engageCallback(std_srvs::Empty::Request&, std_srvs::Empty::Response&)
  {
    boost::mutex::scoped_lock lock(command_mutex_);

    ROS_INFO_NAMED("twist_controller", "Engaging motors!");
    motors_running_ = true;
    return true;
  }

  bool shutdownCallback(std_srvs::Empty::Request&, std_srvs::Empty::Response&)
  {
    boost::mutex::scoped_lock lock(command_mutex_);

    ROS_INFO_NAMED("twist_controller", "Shutting down motors!");
    motors_running_ = false;
    return true;
  }

  void starting(const ros::Time &time)
  {
    reset();
    wrench_output_->start();
  }

  void stopping(const ros::Time &time)
  {
    wrench_output_->stop();
  }

  void update(const ros::Time& time, const ros::Duration& period)
  {
    boost::mutex::scoped_lock lock(command_mutex_);

    // Get twist command input
    if (twist_input_->connected() && twist_input_->enabled()) {
      command_.twist = twist_input_->getCommand();
      command_given_in_stabilized_frame_ = false;
    }

    // Get current state and command
    Twist command = command_.twist;
    Twist twist = twist_->twist();
    Twist twist_body;
    twist_body.linear =  pose_->toBody(twist.linear);
    twist_body.angular = pose_->toBody(twist.angular);

    // Transform to world coordinates if necessary (yaw only)
    if (command_given_in_stabilized_frame_) {
      double yaw = pose_->getYaw();
      Twist transformed = command;
      transformed.linear.x  = cos(yaw) * command.linear.x  - sin(yaw) * command.linear.y;
      transformed.linear.y  = sin(yaw) * command.linear.x  + cos(yaw) * command.linear.y;
      transformed.angular.x = cos(yaw) * command.angular.x - sin(yaw) * command.angular.y;
      transformed.angular.y = sin(yaw) * command.angular.x + cos(yaw) * command.angular.y;
      command = transformed;
    }

    // Get gravity and load factor
    const double gravity = 9.8065;
    double load_factor = 1. / (  pose_->pose().orientation.w * pose_->pose().orientation.w
                                 - pose_->pose().orientation.x * pose_->pose().orientation.x
                                 - pose_->pose().orientation.y * pose_->pose().orientation.y
                                 + pose_->pose().orientation.z * pose_->pose().orientation.z );
    // Note: load_factor could be NaN or Inf...?
    if (load_factor_limit > 0.0 && !(load_factor < load_factor_limit)) load_factor = load_factor_limit;

    // Auto engage/shutdown
    if (auto_engage_) {
      if (!motors_running_ && command.linear.z > 0.1 && load_factor > 0.0) {
        motors_running_ = true;
        ROS_INFO_NAMED("twist_controller", "Engaging motors!");
      } else if (motors_running_ && command.linear.z < -0.1 /* && (twist.linear.z > -0.1 && twist.linear.z < 0.1) */) {
        double shutdown_limit = 0.25 * std::min(command.linear.z, -0.5);
        if (linear_z_control_error_ > 0.0) linear_z_control_error_ = 0.0; // positive control errors should not affect shutdown
        if (pid_.linear.z.getFilteredControlError(linear_z_control_error_, 5.0, period) < shutdown_limit) {
          motors_running_ = false;
          ROS_INFO_NAMED("twist_controller", "Shutting down motors!");
        } else {
          ROS_DEBUG_STREAM_NAMED("twist_controller", "z control error = " << linear_z_control_error_ << " >= " << shutdown_limit);
        }
      } else {
        linear_z_control_error_ = 0.0;
      }

      // flip over?
      if (motors_running_ && load_factor < 0.0) {
        motors_running_ = false;
        ROS_WARN_NAMED("twist_controller", "Shutting down motors due to flip over!");
      }
    }

    // Update output
    if (motors_running_) {
      Vector3 acceleration_command;
      acceleration_command.x = pid_.linear.x.update(command.linear.x, twist.linear.x, acceleration_->acceleration().x, period);
      acceleration_command.y = pid_.linear.y.update(command.linear.y, twist.linear.y, acceleration_->acceleration().y, period);
      acceleration_command.z = pid_.linear.z.update(command.linear.z, twist.linear.z, acceleration_->acceleration().z, period) + gravity;
      Vector3 acceleration_command_body = pose_->toBody(acceleration_command);

      ROS_DEBUG_STREAM_NAMED("twist_controller", "twist.linear:               [" << twist.linear.x << " " << twist.linear.y << " " << twist.linear.z << "]");
      ROS_DEBUG_STREAM_NAMED("twist_controller", "twist_body.angular:         [" << twist_body.angular.x << " " << twist_body.angular.y << " " << twist_body.angular.z << "]");
      ROS_DEBUG_STREAM_NAMED("twist_controller", "twist_command.linear:       [" << command.linear.x << " " << command.linear.y << " " << command.linear.z << "]");
      ROS_DEBUG_STREAM_NAMED("twist_controller", "twist_command.angular:      [" << command.angular.x << " " << command.angular.y << " " << command.angular.z << "]");
      ROS_DEBUG_STREAM_NAMED("twist_controller", "acceleration:               [" << acceleration_->acceleration().x << " " << acceleration_->acceleration().y << " " << acceleration_->acceleration().z << "]");
      ROS_DEBUG_STREAM_NAMED("twist_controller", "acceleration_command_world: [" << acceleration_command.x << " " << acceleration_command.y << " " << acceleration_command.z << "]");
      ROS_DEBUG_STREAM_NAMED("twist_controller", "acceleration_command_body:  [" << acceleration_command_body.x << " " << acceleration_command_body.y << " " << acceleration_command_body.z << "]");

      wrench_.wrench.torque.x = inertia_[0] * pid_.angular.x.update(-acceleration_command_body.y / gravity, 0.0, twist_body.angular.x, period);
      wrench_.wrench.torque.y = inertia_[1] * pid_.angular.y.update( acceleration_command_body.x / gravity, 0.0, twist_body.angular.y, period);
      wrench_.wrench.torque.z = inertia_[2] * pid_.angular.z.update( command.angular.z, twist.angular.z, 0.0, period);
      wrench_.wrench.force.x  = 0.0;
      wrench_.wrench.force.y  = 0.0;
      wrench_.wrench.force.z  = mass_ * ((acceleration_command.z - gravity) * load_factor + gravity);

      if (limits_.force.z > 0.0 && wrench_.wrench.force.z > limits_.force.z) wrench_.wrench.force.z = limits_.force.z;
      if (wrench_.wrench.force.z <= std::numeric_limits<double>::min()) wrench_.wrench.force.z = std::numeric_limits<double>::min();
      if (limits_.torque.x > 0.0) {
        if (wrench_.wrench.torque.x >  limits_.torque.x) wrench_.wrench.torque.x =  limits_.torque.x;
        if (wrench_.wrench.torque.x < -limits_.torque.x) wrench_.wrench.torque.x = -limits_.torque.x;
      }
      if (limits_.torque.y > 0.0) {
        if (wrench_.wrench.torque.y >  limits_.torque.y) wrench_.wrench.torque.y =  limits_.torque.y;
        if (wrench_.wrench.torque.y < -limits_.torque.y) wrench_.wrench.torque.y = -limits_.torque.y;
      }
      if (limits_.torque.z > 0.0) {
        if (wrench_.wrench.torque.z >  limits_.torque.z) wrench_.wrench.torque.z =  limits_.torque.z;
        if (wrench_.wrench.torque.z < -limits_.torque.z) wrench_.wrench.torque.z = -limits_.torque.z;
      }

      ROS_DEBUG_STREAM_NAMED("twist_controller", "wrench_command.force:       [" << wrench_.wrench.force.x << " " << wrench_.wrench.force.y << " " << wrench_.wrench.force.z << "]");
      ROS_DEBUG_STREAM_NAMED("twist_controller", "wrench_command.torque:      [" << wrench_.wrench.torque.x << " " << wrench_.wrench.torque.y << " " << wrench_.wrench.torque.z << "]");

    } else {
      reset();
    }

    // set wrench output
    wrench_.header.stamp = time;
    wrench_.header.frame_id = base_link_frame_;
    wrench_output_->setCommand(wrench_.wrench);
  }

private:
  PoseHandlePtr pose_;
  TwistHandlePtr twist_;
  AccelerationHandlePtr acceleration_;
  TwistCommandHandlePtr twist_input_;
  WrenchCommandHandlePtr wrench_output_;

  ros::NodeHandle node_handle_;
  ros::Subscriber twist_subscriber_;
  ros::Subscriber cmd_vel_subscriber_;
  ros::ServiceServer engage_service_server_;
  ros::ServiceServer shutdown_service_server_;

  geometry_msgs::TwistStamped command_;
  geometry_msgs::WrenchStamped wrench_;
  bool command_given_in_stabilized_frame_;
  std::string base_link_frame_;

  struct {
    struct {
      PID x;
      PID y;
      PID z;
    } linear, angular;
  } pid_;

  geometry_msgs::Wrench limits_;
  bool auto_engage_;
  double load_factor_limit;
  double mass_;
  double inertia_[3];

  bool motors_running_;
  double linear_z_control_error_;
  boost::mutex command_mutex_;

};

} // namespace hector_quadrotor_controller

#include <pluginlib/class_list_macros.h>
PLUGINLIB_EXPORT_CLASS(hector_quadrotor_controller::TwistController, controller_interface::ControllerBase)
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_quadrotor/hector_quadrotor_controller/src/quadrotor_interface.cpp
//=================================================================================================
// Copyright (c) 2013, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_quadrotor_controller/quadrotor_interface.h>

#include <cmath>

namespace hector_quadrotor_controller {

QuadrotorInterface::QuadrotorInterface()
{}

QuadrotorInterface::~QuadrotorInterface()
{}

bool QuadrotorInterface::enabled(const CommandHandle *handle) const
{
  if (!handle || !handle->connected()) return false;
  std::string resource = handle->getName();
  return enabled_.count(resource) > 0;
}

bool QuadrotorInterface::start(const CommandHandle *handle)
{
  if (!handle || !handle->connected()) return false;
  if (enabled(handle)) return true;
  std::string resource = handle->getName();
  enabled_[resource] = handle;
  ROS_DEBUG_NAMED("quadrotor_interface", "Enabled %s control", resource.c_str());
  return true;
}

void QuadrotorInterface::stop(const CommandHandle *handle)
{
  if (!handle) return;
  if (!enabled(handle)) return;
  std::string resource = handle->getName();
  std::map<std::string, const CommandHandle *>::iterator it = enabled_.lower_bound(resource);
  if (it != enabled_.end() && it->second == handle) enabled_.erase(it);
  ROS_DEBUG_NAMED("quadrotor_interface", "Disabled %s control", resource.c_str());
}

void QuadrotorInterface::disconnect(const CommandHandle *handle)
{
  if (!handle) return;
  std::string resource = handle->getName();
  if (inputs_.count(resource)) {
    const CommandHandlePtr& input = inputs_.at(resource);
    if (input.get() != handle) input->reset();
  }
  if (outputs_.count(resource)) {
    const CommandHandlePtr& output = outputs_.at(resource);
    if (output.get() != handle) output->reset();
  }
}

const Pose *QuadrotorInterface::getPoseCommand()          const { return getCommand<PoseCommandHandle>("pose"); }
const Twist *QuadrotorInterface::getTwistCommand()        const { return getCommand<TwistCommandHandle>("twist"); }
const Wrench *QuadrotorInterface::getWrenchCommand()      const { return getCommand<WrenchCommandHandle>("wrench"); }
const MotorCommand *QuadrotorInterface::getMotorCommand() const { return getCommand<MotorCommandHandle>("motor"); }

void PoseHandle::getEulerRPY(double &roll, double &pitch, double &yaw) const
{
  const Quaternion::_w_type& w = pose().orientation.w;
  const Quaternion::_x_type& x = pose().orientation.x;
  const Quaternion::_y_type& y = pose().orientation.y;
  const Quaternion::_z_type& z = pose().orientation.z;
  roll  =  atan2(2.*y*z + 2.*w*x, z*z - y*y - x*x + w*w);
  pitch = -asin(2.*x*z - 2.*w*y);
  yaw   =  atan2(2.*x*y + 2.*w*z, x*x + w*w - z*z - y*y);
}

double PoseHandle::getYaw() const
{
  const Quaternion::_w_type& w = pose().orientation.w;
  const Quaternion::_x_type& x = pose().orientation.x;
  const Quaternion::_y_type& y = pose().orientation.y;
  const Quaternion::_z_type& z = pose().orientation.z;
  return atan2(2.*x*y + 2.*w*z, x*x + w*w - z*z - y*y);
}

Vector3 PoseHandle::toBody(const Vector3& nav) const
{
  const Quaternion::_w_type& w = pose().orientation.w;
  const Quaternion::_x_type& x = pose().orientation.x;
  const Quaternion::_y_type& y = pose().orientation.y;
  const Quaternion::_z_type& z = pose().orientation.z;
  Vector3 body;
  body.x = (w*w+x*x-y*y-z*z) * nav.x + (2.*x*y + 2.*w*z) * nav.y + (2.*x*z - 2.*w*y) * nav.z;
  body.y = (2.*x*y - 2.*w*z) * nav.x + (w*w-x*x+y*y-z*z) * nav.y + (2.*y*z + 2.*w*x) * nav.z;
  body.z = (2.*x*z + 2.*w*y) * nav.x + (2.*y*z - 2.*w*x) * nav.y + (w*w-x*x-y*y+z*z) * nav.z;
  return body;
}

Vector3 PoseHandle::fromBody(const Vector3& body) const
{
  const Quaternion::_w_type& w = pose().orientation.w;
  const Quaternion::_x_type& x = pose().orientation.x;
  const Quaternion::_y_type& y = pose().orientation.y;
  const Quaternion::_z_type& z = pose().orientation.z;
  Vector3 nav;
  nav.x = (w*w+x*x-y*y-z*z) * body.x + (2.*x*y - 2.*w*z) * body.y + (2.*x*z + 2.*w*y) * body.z;
  nav.y = (2.*x*y + 2.*w*z) * body.x + (w*w-x*x+y*y-z*z) * body.y + (2.*y*z - 2.*w*x) * body.z;
  nav.z = (2.*x*z - 2.*w*y) * body.x + (2.*y*z + 2.*w*x) * body.y + (w*w-x*x-y*y+z*z) * body.z;
  return nav;
}

void HorizontalPositionCommandHandle::getError(const PoseHandle &pose, double &x, double &y) const
{
  getCommand(x, y);
  x -= pose.get()->position.x;
  y -= pose.get()->position.y;
}

double HeightCommandHandle::getError(const PoseHandle &pose) const
{
  return getCommand() - pose.get()->position.z;
}

void HeadingCommandHandle::setCommand(double command)
{
  if (get()) {
    get()->x = 0.0;
    get()->y = 0.0;
    get()->z = sin(command / 2.);
    get()->w = cos(command / 2.);
  }

  if (scalar_) {
    *scalar_ = command;
  }
}

double HeadingCommandHandle::getCommand() const {
  if (scalar_) return *scalar_;
  const Quaternion::_w_type& w = get()->w;
  const Quaternion::_x_type& x = get()->x;
  const Quaternion::_y_type& y = get()->y;
  const Quaternion::_z_type& z = get()->z;
  return atan2(2.*x*y + 2.*w*z, x*x + w*w - z*z - y*y);
}

bool HeadingCommandHandle::update(Pose& command) const {
  if (get()) {
    command.orientation = *get();
    return true;
  }
  if (scalar_) {
    command.orientation.x = 0.0;
    command.orientation.y = 0.0;
    command.orientation.z = sin(*scalar_ / 2.);
    command.orientation.x = cos(*scalar_ / 2.);
    return true;
  }
  return false;
}

double HeadingCommandHandle::getError(const PoseHandle &pose) const {
  static const double M_2PI = 2.0 * M_PI;
  double error = getCommand() - pose.getYaw() + M_PI;
  error -= floor(error / M_2PI) * M_2PI;
  return error - M_PI;
}

bool CommandHandle::enabled()    { return interface_->enabled(this); }
bool CommandHandle::start()      { return interface_->start(this); }
void CommandHandle::stop()       { interface_->stop(this); }
void CommandHandle::disconnect() { interface_->disconnect(this); }

} // namespace hector_quadrotor_controller
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_quadrotor/hector_quadrotor_teleop/src/quadrotor_teleop.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================


#include <ros/ros.h>
#include <sensor_msgs/Joy.h>
#include <geometry_msgs/Twist.h>
#include <hector_uav_msgs/YawrateCommand.h>
#include <hector_uav_msgs/ThrustCommand.h>
#include <hector_uav_msgs/AttitudeCommand.h>

namespace hector_quadrotor
{

class Teleop
{
private:
  ros::NodeHandle node_handle_;
  ros::Subscriber joy_subscriber_;

  ros::Publisher velocity_publisher_, attitude_publisher_, yawrate_publisher_, thrust_publisher_;
  geometry_msgs::Twist velocity_;
  hector_uav_msgs::AttitudeCommand attitude_;
  hector_uav_msgs::ThrustCommand thrust_;
  hector_uav_msgs::YawrateCommand yawrate_;

  struct Axis
  {
    int axis;
    double max;
  };

  struct Button
  {
    int button;
  };

  struct
  {
    Axis x;
    Axis y;
    Axis z;
    Axis yaw;
  } axes_;

  struct
  {
    Button slow;
  } buttons_;

  double slow_factor_;

public:
  Teleop()
  {
    ros::NodeHandle params("~");

    params.param<int>("x_axis", axes_.x.axis, 4);
    params.param<int>("y_axis", axes_.y.axis, 3);
    params.param<int>("z_axis", axes_.z.axis, 2);
    params.param<int>("yaw_axis", axes_.yaw.axis, 1);

    params.param<double>("yaw_velocity_max", axes_.yaw.max, 90.0 * M_PI / 180.0);
    params.param<int>("slow_button", buttons_.slow.button, 1);
    params.param<double>("slow_factor", slow_factor_, 0.2);

    std::string control_mode_str;
    params.param<std::string>("control_mode", control_mode_str, "twist");

    if (control_mode_str == "twist")
    {
      params.param<double>("x_velocity_max", axes_.x.max, 2.0);
      params.param<double>("y_velocity_max", axes_.y.max, 2.0);
      params.param<double>("z_velocity_max", axes_.z.max, 2.0);

      joy_subscriber_ = node_handle_.subscribe<sensor_msgs::Joy>("joy", 1, boost::bind(&Teleop::joyTwistCallback, this, _1));
      velocity_publisher_ = node_handle_.advertise<geometry_msgs::Twist>("cmd_vel", 10);
    }
    else if (control_mode_str == "attitude")
    {
      params.param<double>("x_roll_max", axes_.x.max, 0.35);
      params.param<double>("y_pitch_max", axes_.y.max, 0.35);
      params.param<double>("z_thrust_max", axes_.z.max, 25.0);
      joy_subscriber_ = node_handle_.subscribe<sensor_msgs::Joy>("joy", 1, boost::bind(&Teleop::joyAttitudeCallback, this, _1));
      attitude_publisher_ = node_handle_.advertise<hector_uav_msgs::AttitudeCommand>("command/attitude", 10);
      yawrate_publisher_ = node_handle_.advertise<hector_uav_msgs::YawrateCommand>("command/yawrate", 10);
      thrust_publisher_ = node_handle_.advertise<hector_uav_msgs::ThrustCommand>("command/thrust", 10);
    }

  }

  ~Teleop()
  {
    stop();
  }

  void joyTwistCallback(const sensor_msgs::JoyConstPtr &joy)
  {
    velocity_.linear.x = getAxis(joy, axes_.x);
    velocity_.linear.y = getAxis(joy, axes_.y);
    velocity_.linear.z = getAxis(joy, axes_.z);
    velocity_.angular.z = getAxis(joy, axes_.yaw);
    if (getButton(joy, buttons_.slow.button))
    {
      velocity_.linear.x *= slow_factor_;
      velocity_.linear.y *= slow_factor_;
      velocity_.linear.z *= slow_factor_;
      velocity_.angular.z *= slow_factor_;
    }
    velocity_publisher_.publish(velocity_);
  }

  void joyAttitudeCallback(const sensor_msgs::JoyConstPtr &joy)
  {
    attitude_.roll = getAxis(joy, axes_.x);
    attitude_.pitch = getAxis(joy, axes_.y);
    attitude_publisher_.publish(attitude_);

    thrust_.thrust = getAxis(joy, axes_.z);
    thrust_publisher_.publish(thrust_);

    yawrate_.turnrate = getAxis(joy, axes_.yaw);
    if (getButton(joy, buttons_.slow.button))
    {
      yawrate_.turnrate *= slow_factor_;
    }
    yawrate_publisher_.publish(yawrate_);
  }

  sensor_msgs::Joy::_axes_type::value_type getAxis(const sensor_msgs::JoyConstPtr &joy, Axis axis)
  {
    if (axis.axis == 0)
    {return 0;}
    sensor_msgs::Joy::_axes_type::value_type sign = 1.0;
    if (axis.axis < 0)
    {
      sign = -1.0;
      axis.axis = -axis.axis;
    }
    if ((size_t) axis.axis > joy->axes.size())
    {return 0;}
    return sign * joy->axes[axis.axis - 1] * axis.max;
  }

  sensor_msgs::Joy::_buttons_type::value_type getButton(const sensor_msgs::JoyConstPtr &joy, int button)
  {
    if (button <= 0)
    {return 0;}
    if ((size_t) button > joy->buttons.size())
    {return 0;}
    return joy->buttons[button - 1];
  }

  void stop()
  {
    if(velocity_publisher_.getNumSubscribers() > 0)
    {
      velocity_ = geometry_msgs::Twist();
      velocity_publisher_.publish(velocity_);
    }
    if(attitude_publisher_.getNumSubscribers() > 0)
    {
      attitude_ = hector_uav_msgs::AttitudeCommand();
      attitude_publisher_.publish(attitude_);
    }
    if(thrust_publisher_.getNumSubscribers() > 0)
    {
      thrust_ = hector_uav_msgs::ThrustCommand();
      thrust_publisher_.publish(thrust_);
    }
    if(yawrate_publisher_.getNumSubscribers() > 0)
    {
      yawrate_ = hector_uav_msgs::YawrateCommand();
      yawrate_publisher_.publish(yawrate_);
    }

  }
};

} // namespace hector_quadrotor

int main(int argc, char **argv)
{
  ros::init(argc, argv, "quadrotor_teleop");

  hector_quadrotor::Teleop teleop;
  ros::spin();

  return 0;
}
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_quadrotor/hector_quadrotor_pose_estimation/src/pose_estimation_nodelet.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_quadrotor_pose_estimation/pose_estimation_node.h>
#include <nodelet/nodelet.h>

namespace hector_quadrotor_pose_estimation {

class QuadrotorPoseEstimationNodelet : public QuadrotorPoseEstimationNode, public nodelet::Nodelet
{
public:
  QuadrotorPoseEstimationNodelet(const SystemPtr& system = SystemPtr())
    : QuadrotorPoseEstimationNode(system)
  {}

private:
  void onInit() {
    QuadrotorPoseEstimationNode::init();
  }

  void onReset() {
    QuadrotorPoseEstimationNode::reset();
  }

  void onCleanup() {
    QuadrotorPoseEstimationNode::cleanup();
  }
};

} // namespace hector_quadrotor_pose_estimation

#include <pluginlib/class_list_macros.h>
PLUGINLIB_EXPORT_CLASS(hector_quadrotor_pose_estimation::QuadrotorPoseEstimationNodelet, nodelet::Nodelet)
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_quadrotor/hector_quadrotor_pose_estimation/src/main.cpp
//=================================================================================================
// Copyright (c) 2012, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_quadrotor_pose_estimation/pose_estimation_node.h>

int main(int argc, char **argv) {
  ros::init(argc, argv, "pose_estimation");
  hector_quadrotor_pose_estimation::QuadrotorPoseEstimationNode node;
  if (!node.init()) return 1;

  ros::spin();

  node.cleanup();
  return 0;
}
=================================================================================================
File: ./medical-drone-simulation/hector_quadrotor_sim/hector_quadrotor/hector_quadrotor_pose_estimation/src/pose_estimation_node.cpp
//=================================================================================================
// Copyright (c) 2011, Johannes Meyer, TU Darmstadt
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Flight Systems and Automatic Control group,
//       TU Darmstadt, nor the names of its contributors may be used to
//       endorse or promote products derived from this software without
//       specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================================

#include <hector_quadrotor_pose_estimation/pose_estimation_node.h>
#include <hector_pose_estimation/measurements/baro.h>

namespace hector_quadrotor_pose_estimation {

QuadrotorPoseEstimationNode::QuadrotorPoseEstimationNode(const SystemPtr& system, const StatePtr& state)
  : PoseEstimationNode(system, state)
{
  pose_estimation_->addMeasurement(new Baro("baro"));
}

QuadrotorPoseEstimationNode::~QuadrotorPoseEstimationNode()
{
}

bool QuadrotorPoseEstimationNode::init() {
  if (!PoseEstimationNode::init()) return false;
  baro_subscriber_ = getNodeHandle().subscribe("altimeter", 10, &QuadrotorPoseEstimationNode::baroCallback, this);
  height_subscriber_.shutdown();
  return true;
}

void QuadrotorPoseEstimationNode::baroCallback(const hector_uav_msgs::AltimeterConstPtr& altimeter) {
  pose_estimation_->getMeasurement("baro")->add(Baro::Update(altimeter->pressure, altimeter->qnh));

  if (sensor_pose_publisher_) {
    boost::shared_ptr<Baro> baro = boost::static_pointer_cast<Baro>(pose_estimation_->getMeasurement("baro"));
    sensor_pose_.pose.position.z = baro->getModel()->getAltitude(Baro::Update(altimeter->pressure, altimeter->qnh)) - baro->getElevation();
  }
}

} // namespace hector_quadrotor_pose_estimation
