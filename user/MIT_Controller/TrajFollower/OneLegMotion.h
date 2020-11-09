#ifndef ONELEGMOTION_H
#define ONELEGMOTION_H

#include <Eigen/Eigen>

#include "CubicSpline.h"

class OneLegMotion
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    OneLegMotion() = default;

    void SetMotion(int leg, double beg, double end,
                   bool beg_contact, bool end_contact,
                   const Eigen::Vector3d &beg_pos,
                   bool beg_knee_out, bool beg_hip_out,
                   const Eigen::Vector3d &end_pos,
                   bool end_knee_out, bool end_hip_out);

    int GetID() const
    {
        return leg_id_;
    }

    void Sample(double t, Eigen::Vector3d &pos, Eigen::Vector3d &vel,
                bool &contact, bool &hip_out, bool &knee_out) const;

    void SampleJoint(double t,
                     Eigen::Vector3d &jpos, Eigen::Vector3d &jvel) const;

    Eigen::Vector3d GetTorsoOffset(double t) const;

private:
    int leg_id_;
    double beg_time_;
    double end_time_;
    bool beg_contact_;
    bool end_contact_;

    CubicSpline<Eigen::Vector3d> joint_spline_;
};

#endif /* ONELEGMOTION_H */
