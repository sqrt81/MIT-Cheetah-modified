#include "KinWBC.hpp"
#include <Utilities/Utilities_print.h>
#include <Utilities/pseudoInverse.h>

#include "DogPhysics.h"

template <typename T>
KinWBC<T>::KinWBC(size_t num_qdot)
    : threshold_(0.001), num_qdot_(num_qdot), num_act_joint_(num_qdot - 6) {
  I_mtx = DMat<T>::Identity(num_qdot_, num_qdot_);
}

//template <typename T>
//bool KinWBC<T>::FindConfiguration(
//    const DVec<T>& curr_config, const std::vector<Task<T>*>& task_list,
//    const std::vector<ContactSpec<T>*>& contact_list, DVec<T>& jpos_cmd,
//    DVec<T>& jvel_cmd) {

//  // Contact Jacobian Setup
//  DMat<T> Nc(num_qdot_, num_qdot_); Nc.setIdentity();
//  if(contact_list.size() > 0){
//    DMat<T> Jc, Jc_i;
//    contact_list[0]->getContactJacobian(Jc);
//    size_t num_rows = Jc.rows();

//    for (size_t i(1); i < contact_list.size(); ++i) {
//      contact_list[i]->getContactJacobian(Jc_i);
//      size_t num_new_rows = Jc_i.rows();
//      Jc.conservativeResize(num_rows + num_new_rows, num_qdot_);
//      Jc.block(num_rows, 0, num_new_rows, num_qdot_) = Jc_i;
//      num_rows += num_new_rows;
//    }

//    // Projection Matrix
//    _BuildProjectionMatrix(Jc, Nc);
//  }

//  // First Task
//  DVec<T> delta_q, qdot;
//  DMat<T> Jt, JtPre, JtPre_pinv, N_nx, N_pre;

//  Task<T>* task = task_list[0];
//  task->getTaskJacobian(Jt);
//  JtPre = Jt * Nc;
//  _PseudoInverse(JtPre, JtPre_pinv);

//  delta_q = JtPre_pinv * (task->getPosError());
//  qdot = JtPre_pinv * (task->getDesVel());

//  DVec<T> prev_delta_q = delta_q;
//  DVec<T> prev_qdot = qdot;

//  _BuildProjectionMatrix(JtPre, N_nx);
//  N_pre = Nc * N_nx;

//  for (size_t i(1); i < task_list.size(); ++i) {
//    task = task_list[i];

//    task->getTaskJacobian(Jt);
//    JtPre = Jt * N_pre;

//    _PseudoInverse(JtPre, JtPre_pinv);
//    delta_q =
//        prev_delta_q + JtPre_pinv * (task->getPosError() - Jt * prev_delta_q);
//    qdot = prev_qdot + JtPre_pinv * (task->getDesVel() - Jt * prev_qdot);

//    // For the next task
//    _BuildProjectionMatrix(JtPre, N_nx);
//    N_pre *= N_nx;
//    prev_delta_q = delta_q;
//    prev_qdot = qdot;
//  }
//  for (size_t i(0); i < num_act_joint_; ++i) {
//    jpos_cmd[i] = curr_config[i + 6] + delta_q[i + 6];
//    jvel_cmd[i] = qdot[i + 6];
//  }

//  std::cout << "jpos: " << jpos_cmd.transpose() << std::endl;
//  std::cout << "jvel: " << jvel_cmd.transpose() << std::endl;
//  return true;
//}


// override the original method with a less time consuming one.
template <typename T>
bool KinWBC<T>::FindConfiguration(
    const DVec<T>& curr_config, const std::vector<Task<T>*>& task_list,
    const std::vector<ContactSpec<T>*>& contact_list, DVec<T>& jpos_cmd,
    DVec<T>& jvel_cmd) {

    (void) curr_config;
    (void) contact_list;

    if(task_list.size() < 2)
        return true;
    // 6 total tasks for locomotion ctrl:
    // first two are torso orientation and position,
    // last four are foot position.
    Vec4<T> torso_r = task_list[0]->getDesPos();
    Eigen::Quaternion<T> rot(torso_r[0], torso_r[1],
            torso_r[2], torso_r[3]);
    rot = rot.inverse();
    Vec3<T> rot_vel = task_list[0]->getDesVel();
    Vec3<T> trans = task_list[1]->getDesPos();
    Vec3<T> vel = task_list[1]->getDesVel();

    if(task_list.size() < 6)
        return true;

    for (size_t i(0); i < 4; i++)
    {
        Task<T>& foot_task = *task_list[i + 2];
        Vec3<T> rel_pos = foot_task.getDesPos() - trans;
        Vec3<T> rel_vel = rot * (foot_task.getDesVel() - vel
                                 - rot_vel.cross(rel_pos));
        rel_pos = rot * rel_pos; // transfer to local

        if (i >= 2) {
            rel_pos.x() = - rel_pos.x();
            rel_vel.x() = - rel_vel.x();
        }
        if (i % 2 == 0) {
            rel_pos.y() = - rel_pos.y();
            rel_vel.y() = - rel_vel.y();
        }

//        std::cout << "foot " << i << " : "
//                  << rel_pos.transpose() << std::endl;

        Vec3<T> joint_pos = DogPhysics::InverseKinematics(
                    rel_pos.template cast<double>(), i >= 2, true)
                .template cast<T>();
        Vec3<T> joint_vel = DogPhysics::ComputeJointVel(
                    rel_vel.template cast<double>(),
                    joint_pos.template cast<double>()).template cast<T>();
        if (i % 2 == 0) {
            joint_pos[0] = - joint_pos[0];
            joint_vel[0] = - joint_vel[0];
        }

        if (i < 2) {
            joint_pos[1] = - joint_pos[1];
            joint_pos[2] = - joint_pos[2];
            joint_vel[1] = - joint_vel[1];
            joint_vel[2] = - joint_vel[2];
        }

        jpos_cmd.template segment<3>(i * 3) = joint_pos;
        jvel_cmd.template segment<3>(i * 3) = joint_vel;
    }

//    std::cout << "jpos: " << jpos_cmd.transpose() << std::endl;
//    std::cout << "jvel: " << jvel_cmd.transpose() << std::endl;

    return true;
}

template <typename T>
void KinWBC<T>::_BuildProjectionMatrix(const DMat<T>& J, DMat<T>& N) {
  DMat<T> J_pinv;
  _PseudoInverse(J, J_pinv);
  N = I_mtx - J_pinv * J;
}

template <typename T>
void KinWBC<T>::_PseudoInverse(const DMat<T> J, DMat<T>& Jinv) {
  pseudoInverse(J, threshold_, Jinv);
}

template class KinWBC<float>;
template class KinWBC<double>;
