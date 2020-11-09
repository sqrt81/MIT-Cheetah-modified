#include "OneLegMotion.h"
#include "Math.h"
#include "CubicSplineImpl.hpp"
#include "WBC/WBIC/DogPhysics.h"

void OneLegMotion::SetMotion(int leg, double beg, double end,
                             bool beg_contact, bool end_contact,
                             const Eigen::Vector3d &beg_pos,
                             bool beg_knee_out, bool beg_hip_out,
                             const Eigen::Vector3d &end_pos,
                             bool end_knee_out, bool end_hip_out)
{
    leg_id_ = leg;
    beg_time_ = beg;
    end_time_ = end;
    beg_contact_ = beg_contact;
    end_contact_ = end_contact;

    Eigen::Vector3d beg_pos_ = beg_pos;
    Eigen::Vector3d end_pos_ = end_pos;

    Eigen::Vector3d beg_vel(0., 0., beg_contact_ ? 0.5 : 0.);
    Eigen::Vector3d end_vel(0., 0., end_contact_ ? - 0.5 : 0.);

    if (leg_id_ >= 2) {
        beg_pos_.x() = - beg_pos_.x();
        beg_pos_.x() = - beg_pos_.x();
        end_pos_.x() = - end_pos_.x();
        end_pos_.x() = - end_pos_.x();
    }
    if (leg_id_ % 2 == 0) {
        beg_pos_.y() = - beg_pos_.y();
        beg_pos_.y() = - beg_pos_.y();
        end_pos_.y() = - end_pos_.y();
        end_pos_.y() = - end_pos_.y();
    }

    Eigen::Vector3d beg_jpos = DogPhysics::InverseKinematics(
                beg_pos_, beg_knee_out, beg_hip_out);
    Eigen::Vector3d beg_jvel = DogPhysics::ComputeJointVel(
                beg_vel, beg_jpos);
    Eigen::Vector3d end_jpos = DogPhysics::InverseKinematics(
                end_pos_, end_knee_out, end_hip_out);
    Eigen::Vector3d end_jvel = DogPhysics::ComputeJointVel(
                end_vel, end_jpos);

    joint_spline_ = CubicSpline<Eigen::Vector3d>(
                beg_time_, beg_jpos, beg_jvel, end_time_, end_jpos, end_jvel);
}



void OneLegMotion::Sample(double t,
                          Eigen::Vector3d &pos,
                          Eigen::Vector3d &vel,
                          bool &contact,
                          bool &hip_out, bool &knee_out) const
{
    Eigen::Vector3d jpos;
    Eigen::Vector3d jvel;

    if (t < beg_time_)
    {
        joint_spline_.Sample(beg_time_, jpos, jvel);
        pos = DogPhysics::ForwardKinematics(jpos);
        vel.setZero();
        contact = beg_contact_;
    }
    else if (t > end_time_)
    {
        joint_spline_.Sample(end_time_, jpos, jvel);
        pos = DogPhysics::ForwardKinematics(jpos);
        vel.setZero();
        contact = end_contact_;
    }
    else
    {
        joint_spline_.Sample(t, jpos, jvel);
        pos = DogPhysics::ForwardKinematics(jpos);
        vel = DogPhysics::ComputeFootVel(jpos, jvel);
        contact = false;
    }

    if (leg_id_ >= 2) {
        pos.x() = - pos.x();
        vel.x() = - vel.x();
    }
    if (leg_id_ % 2 == 0) {
        pos.y() = - pos.y();
        vel.y() = - vel.y();
    }

    DogPhysics::GetLegConfig(jpos, hip_out, knee_out);
}

void OneLegMotion::SampleJoint(double t,
                               Eigen::Vector3d &jpos,
                               Eigen::Vector3d &jvel) const
{
    if (t < beg_time_)
    {
        joint_spline_.Sample(beg_time_, jpos, jvel);
        jvel.setZero();
    }
    else if (t > end_time_)
    {
        joint_spline_.Sample(end_time_, jpos, jvel);
        jvel.setZero();
    }
    else
    {
        joint_spline_.Sample(t, jpos, jvel);
    }

    if (leg_id_ % 2 == 0) {
        jpos[0] = - jpos[0];
        jvel[0] = - jvel[0];
    }

    if (leg_id_ < 2) {
        jpos[1] = - jpos[1];
        jpos[2] = - jpos[2];
        jvel[1] = - jvel[1];
        jvel[2] = - jvel[2];
    }
}

Eigen::Vector3d OneLegMotion::GetTorsoOffset(double t) const
{
    t = clamp(t, beg_time_, end_time_);

    const double coef1 = (t - beg_time_) / (end_time_ - beg_time_);
//    const double coef1 = 1.;
    const double coef2 = (end_time_ - t) / (end_time_ - beg_time_);

    double factor = 0;

    if (!beg_contact_)
        factor += coef2;

    if (!end_contact_)
        factor += coef1;

    switch (leg_id_)
    {
    case 0: // FR
        return factor * Eigen::Vector3d(- 0.03, 0.03, 0);
    case 1: // FL
        return factor * Eigen::Vector3d(- 0.03, - 0.03, 0);
    case 2: // BR
        return factor * Eigen::Vector3d(0.03, 0.03, 0);
    case 3: // BL
        return factor * Eigen::Vector3d(0.03, - 0.03, 0);
    default:
        return {0, 0, 0};
    }
}
