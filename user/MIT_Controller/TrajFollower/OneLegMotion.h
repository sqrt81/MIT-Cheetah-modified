#ifndef ONELEGMOTION_H
#define ONELEGMOTION_H

#include <Eigen/Eigen>

class OneLegMotion
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    OneLegMotion() = default;

    void SetMotion(int leg, double beg, double end,
                   bool beg_contact, bool end_contact,
                   const Eigen::Vector3d &beg_pos,
                   const Eigen::Vector3d &end_pos);

    int GetID() const
    {
        return leg_id_;
    }

    void Sample(double t, Eigen::Vector3d &pos, Eigen::Vector3d &vel,
                bool &contact) const;

    Eigen::Vector3d GetTorsoOffset(double t) const;

private:
    int leg_id_;
    double beg_time_;
    double critic_time_;
    double end_time_;
    bool beg_contact_;
    bool end_contact_;
    Eigen::Vector3d beg_pos_;
    Eigen::Vector3d critic_offset_;
    Eigen::Vector3d target_local_pos_;
};

#endif /* ONELEGMOTION_H */
