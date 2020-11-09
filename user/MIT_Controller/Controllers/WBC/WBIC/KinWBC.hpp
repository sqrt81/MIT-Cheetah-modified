#ifndef KINEMATICS_WHOLE_BODY_CONTROL
#define KINEMATICS_WHOLE_BODY_CONTROL

#include <WBC/ContactSpec.hpp>
#include <WBC/Task.hpp>
#include <vector>

template <typename T>
class KinWBC {
 public:
  KinWBC(size_t num_qdot);
  ~KinWBC() {}

  void SetLegConfig(int leg_id, bool knee_out, bool hip_out)
  {
      leg_config[leg_id][0] = knee_out;
      leg_config[leg_id][1] = hip_out;
  }

  bool FindConfiguration(const DVec<T>& curr_config,
                         const std::vector<Task<T>*>& task_list,
                         const std::vector<ContactSpec<T>*>& contact_list,
                         DVec<T>& jpos_cmd, DVec<T>& jvel_cmd);

  DMat<T> Ainv_;

 private:
  void _PseudoInverse(const DMat<T> J, DMat<T>& Jinv);
  void _BuildProjectionMatrix(const DMat<T>& J, DMat<T>& N);

  bool leg_config[4][2];

  double threshold_;
  size_t num_qdot_;
  size_t num_act_joint_;
  DMat<T> I_mtx;
};
#endif
