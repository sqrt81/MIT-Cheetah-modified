#include "DogPhysics.h"

#include <cmath>

using std::cos;
using std::sin;
using std::acos;
using std::atan2;

namespace DogPhysics
{

// length
constexpr double hip_pos_x_ = 0.19; // hip_pos_x_ = hip_pos_x + hip_len_x
constexpr double hip_pos_y_ = 0.049;
constexpr double hip_len_y_ = 0.062;
constexpr double thigh_offset_z_ = - 0.209;
constexpr double shin_offset_z_ = - 0.2;

Eigen::Vector3d ForwardKinematics(const Eigen::Vector3d &joint_pos)
{
    const double c1 = cos(joint_pos(0));
    const double s1 = sin(joint_pos(0));
    const double c2 = cos(joint_pos(1));
    const double s2 = sin(joint_pos(1));
    const double c23 = cos(joint_pos(1) + joint_pos(2));
    const double s23 = sin(joint_pos(1) + joint_pos(2));

    const double shin_z = thigh_offset_z_ * c2 + shin_offset_z_ * c23;
    const double shin_x = thigh_offset_z_ * s2 + shin_offset_z_ * s23;

    return Eigen::Vector3d(
                hip_pos_x_ + shin_x,
                hip_len_y_ * c1 - shin_z * s1 + hip_pos_y_,
                hip_len_y_ * s1 + shin_z * c1);
}

Eigen::Vector3d InverseKinematics(const Eigen::Vector3d &foot_pos,
                                  bool knee_out, bool hip_out)
{
    const double x_h = foot_pos.x() - hip_pos_x_;
    const double y_h = foot_pos.y() - hip_pos_y_;
    const double z_t = foot_pos.z();
    const double cos_11 = hip_len_y_ / sqrt(z_t * z_t + y_h * y_h);

    const double angle_11 = acos(cos_11 < 1. ? cos_11 : 1.);
    const double angle_12 = atan2(z_t, y_h);

    const double q1 = hip_out ? (angle_11 + angle_12) : (angle_12 - angle_11);

    const double w_h = y_h * sin(q1) - z_t * cos(q1);
    const double leg_len_2 = x_h * x_h + w_h * w_h;
    const double cos_3
            = (leg_len_2 - thigh_offset_z_ * thigh_offset_z_
               - shin_offset_z_ * shin_offset_z_)
              / (2 * thigh_offset_z_ * shin_offset_z_);
    const double angle_3 = acos(std::max(std::min(cos_3, 1.), - 1.));

    const double q3 = knee_out ? angle_3 : (- angle_3);
    const double sin_22 = shin_offset_z_ * sin(q3) / sqrt(leg_len_2);

    const double q21 = - atan2(x_h, w_h);
    const double q2 = (q21 < M_PI_2 ? q21 : q21 - 2 * M_PI)
            + asin(std::max(std::min(sin_22, 1.), - 1.));

    return Eigen::Vector3d(q1, q2, q3);
}

void GetLegConfig(const Eigen::Vector3d &joint_pos,
                  bool &hip_out, bool &knee_out)
{
    knee_out = joint_pos(2) > 0;

    const double c1 = cos(joint_pos(0));
    const double s1 = sin(joint_pos(0));
    const double c2 = cos(joint_pos(1));
    const double c23 = cos(joint_pos(1) + joint_pos(2));

    const double shin_z = thigh_offset_z_ * c2 + shin_offset_z_ * c23;
    const double y_h = hip_len_y_ * c1 - shin_z * s1;
    const double z_h = hip_len_y_ * s1 + shin_z * c1;
    const double angle_12 = atan2(z_h, y_h);

    hip_out = joint_pos(0) > angle_12;
}

Eigen::Matrix3d ComputeJacobian(const Eigen::Vector3d &joint_pos)
{
    const double c1 = cos(joint_pos[0]);
    const double s1 = sin(joint_pos[0]);
    const double c2 = cos(joint_pos[1]);
    const double s2 = sin(joint_pos[1]);
    const double c23 = cos(joint_pos[1] + joint_pos[2]);
    const double s23 = sin(joint_pos[1] + joint_pos[2]);

    const double shin_z = thigh_offset_z_ * c2 + shin_offset_z_ * c23;
    const double shin_x = thigh_offset_z_ * s2 + shin_offset_z_ * s23;

    Eigen::Matrix3d jacob;

    jacob(0, 0) = 0.;
    jacob(0, 1) = shin_z;
    jacob(0, 2) = shin_offset_z_ * c23;
    jacob(1, 0) = - shin_z * c1 - hip_len_y_ * s1;
    jacob(1, 1) = shin_x * s1;
    jacob(1, 2) = shin_offset_z_ * s1 * s23;
    jacob(2, 0) = - shin_z * s1 + hip_len_y_ * c1;
    jacob(2, 1) = - shin_x * c1;
    jacob(2, 2) = - shin_offset_z_ * c1 * s23;

    return jacob;
}

Eigen::Vector3d ComputeFootVel(const Eigen::Vector3d &joint_pos,
                               const Eigen::Vector3d &joint_vel)
{
    const Eigen::Matrix3d jacobian = ComputeJacobian(joint_pos);
    return jacobian * joint_vel;
}

Eigen::Vector3d ComputeJointVel(const Eigen::Vector3d& foot_vel,
                                const Eigen::Vector3d& joint_pos)
{
    const Eigen::Matrix3d jacobian = ComputeJacobian(joint_pos);
    Eigen::Vector3d joint_vel;

    if (abs(jacobian.determinant()) < 1e-3)
    {
        joint_vel =
                (jacobian + Eigen::Matrix3d::Identity() * 1e-3).inverse()
                * foot_vel;
    }
    else
    {
        joint_vel = jacobian.inverse() * foot_vel;
    }

    return joint_vel;
}

} /* DogPhysics */
