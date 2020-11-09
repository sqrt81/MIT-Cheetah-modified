#ifndef DOGPHYSICS_H
#define DOGPHYSICS_H

#include <Eigen/Eigen>

namespace DogPhysics
{
    Eigen::Vector3d ForwardKinematics(const Eigen::Vector3d &joint_pos);

    /**
     * @brief InverseKinematics
     * Computes desired joint positions of FL leg
     * with the target foot position.
     * @param foot_pos      target foot position.
     * @param joint_pos     desired joint position.
     * @param knee_out      if the knee curves outwards
     * @param hip_out       if the thigh is at side of the torso
     */
    Eigen::Vector3d InverseKinematics(const Eigen::Vector3d &foot_pos,
                                      bool knee_out, bool hip_out);

    void GetLegConfig(const Eigen::Vector3d &joint_pos,
                      bool &hip_out, bool &knee_out);

    /**
     * @brief ComputeJacobian
     * Computes joint jacobian.
     * @param joint_pos     joint position of front-left leg.
     * @return              matrix J that satisfies dpos = J * dq.
     */
    Eigen::Matrix3d ComputeJacobian(const Eigen::Vector3d &joint_pos);

    /**
     * @brief ComputeFootPos
     * Computes foot vel according to jacobian
     * @param joint_stat    joint pos and vel
     * @return              foot linear velocity
     */
    Eigen::Vector3d ComputeFootVel(const Eigen::Vector3d &joint_pos,
                                   const Eigen::Vector3d &joint_vel);

    /**
     * @brief ComputeJointVel
     * Given joint pos and foot vel, computes joint vel.
     * @param foot_vel      foot velocity
     * @param joint_stat    provides joint pos and stores computed joint vel
     */
    Eigen::Vector3d ComputeJointVel(const Eigen::Vector3d& foot_vel,
                                    const Eigen::Vector3d& joint_pos);

} /* DogPhysics */

#endif /* DOGPHYSICS_H */
