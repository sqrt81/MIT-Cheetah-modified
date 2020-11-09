#ifndef TRAJHOLDER_H
#define TRAJHOLDER_H

#include <vector>
#include <string>

#include <Eigen/Eigen>

#include "CubicSpline.h"

/**
 * @brief The TrajHolder class
 * This class holds the trajectory of torso, and
 * also the foot positions and foot forces.
 * It can load trajectory from a csv file.
 */
class TrajHolder
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    TrajHolder() = default;

    void SampleTorso(double t,
                     Eigen::Vector3d &pos, Eigen::Vector3d &vel,
                     Eigen::Quaterniond &rot, Eigen::Vector3d &rot_vel) const;

    void SampleFoot(double t, int foot_id,
                    Eigen::Vector3d &pos, Eigen::Vector3d &vel,
                    Eigen::Vector3d &acc, Eigen::Vector3d &force,
                    bool& contact) const;

    void LoadFromFile(const std::string &file_name);

    void SetActive(bool active)
    {
        is_active_ = active;
    }

    bool Active() const
    {
        return is_active_;
    }

    void SetCurTime(double time)
    {
        cur_time_ = time;
    }

    double GetCurTime() const
    {
        return cur_time_;
    }

private:
    double traj_beg_time_;
    double sample_time_interval_;

    bool is_active_ = false;
    double cur_time_;

    std::vector<CubicSpline<Eigen::Vector3d>,
                Eigen::aligned_allocator<CubicSpline<Eigen::Vector3d>>> torso_trans_traj_;
    std::vector<CubicSpline<Eigen::Quaterniond, Eigen::Vector3d>,
                Eigen::aligned_allocator<CubicSpline<Eigen::Quaterniond, Eigen::Vector3d>>> torso_rot_traj_;

    std::array<std::vector<CubicSpline<Eigen::Vector3d>,
                           Eigen::aligned_allocator<CubicSpline<Eigen::Vector3d>>>, 4> foot_traj_;
    std::array<std::vector<Eigen::Vector3d,
                           Eigen::aligned_allocator<Eigen::Vector3d>>, 4> foot_force_;
    std::array<std::vector<bool>, 4> foot_contact_;
};

#endif /* TRAJHOLDER_H */
