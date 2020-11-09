#include "TrajHolder.h"

#include "CubicSplineImpl.hpp"
#include "QuatSplineImpl.hpp"
#include "Math.h"

#include <fstream>
#include <sstream>

void TrajHolder::SampleTorso(double t,
                             Eigen::Vector3d &pos,
                             Eigen::Vector3d &vel,
                             Eigen::Quaterniond &rot,
                             Eigen::Vector3d &rot_vel) const
{
    int index = static_cast<int>((t - traj_beg_time_) / sample_time_interval_);

    if (index < 0)
    {
        index = 0;
        t = traj_beg_time_;
    }
    else if (index >= (int)torso_trans_traj_.size())
    {
        index = torso_trans_traj_.size() - 1;
        t = traj_beg_time_ + sample_time_interval_ * torso_trans_traj_.size();
    }

    torso_trans_traj_[index].Sample(t, pos, vel);
    torso_rot_traj_[index].Sample(t, rot, rot_vel);
}

void TrajHolder::SampleFoot(double t, int foot_id,
                            Eigen::Vector3d &pos,
                            Eigen::Vector3d &vel,
                            Eigen::Vector3d &acc,
                            Eigen::Vector3d &force,
                            bool &contact) const
{
    int index = static_cast<int>((t - traj_beg_time_) / sample_time_interval_);

    if (index < 0)
    {
        index = 0;
        t = traj_beg_time_;
    }
    else if (index >= (int)torso_trans_traj_.size())
    {
        index = torso_trans_traj_.size() - 1;
        t = traj_beg_time_ + sample_time_interval_ * torso_trans_traj_.size();
    }

    foot_id = foot_id ^ 1;

    foot_traj_[foot_id][index].Sample(t, pos, vel, acc);
    force = foot_force_[foot_id][index];
    contact = foot_contact_[foot_id][index];
}

void TrajHolder::LoadFromFile(const std::string &file_name)
{
    // read from csv file
    std::ifstream file(file_name, std::ios::in);

    if (!file.is_open())
        return;

    std::string line_str;

    std::vector<double> timestamps;
    std::vector<Eigen::Vector3d> torso_pos;
    std::vector<Eigen::Vector3d> torso_vel;
    std::vector<Eigen::Quaterniond> torso_rot;
    std::vector<Eigen::Vector3d> torso_rot_vel;
    std::vector<Eigen::Vector3d> foot_pos[4];
    std::vector<Eigen::Vector3d> foot_vel[4];

    torso_trans_traj_.clear();
    torso_rot_traj_.clear();

    for (int i = 0; i < 4; i++)
    {
        foot_traj_[i].clear();
        foot_force_[i].clear();
        foot_contact_[i].clear();
    }

    // omit the header
    std::getline(file, line_str);

    while(std::getline(file, line_str))
    {
        std::string str;
        double x, y, z, w;
        std::stringstream ss(line_str);

        getline(ss, str, ','); // time
        getline(ss, str, ','); // time from start
        timestamps.push_back(std::stoll(str) * 1e-9);

        // torso pos
        getline(ss, str, ',');
        x = std::stod(str);
        getline(ss, str, ',');
        y = std::stod(str);
        getline(ss, str, ',');
        z = std::stod(str);
        torso_pos.push_back({x, y, z});

        // torso rot
        getline(ss, str, ',');
        x = std::stod(str);
        getline(ss, str, ',');
        y = std::stod(str);
        getline(ss, str, ',');
        z = std::stod(str);
        getline(ss, str, ',');
        w = std::stod(str);
        torso_rot.push_back({w, x, y, z});
        (void)w;
//        torso_rot.push_back({1, 0, 0, 0});

        // torso vel
        getline(ss, str, ',');
        x = std::stod(str);
        getline(ss, str, ',');
        y = std::stod(str);
        getline(ss, str, ',');
        z = std::stod(str);
        torso_vel.push_back(torso_rot.back() * Eigen::Vector3d(x, y, z));

        // torso rot vel (change into local)
        getline(ss, str, ',');
        x = std::stod(str);
        getline(ss, str, ',');
        y = std::stod(str);
        getline(ss, str, ',');
        z = std::stod(str);
        torso_rot_vel.push_back(Eigen::Vector3d(x, y, z));
//        torso_rot_vel.push_back({0, 0, 0});

        // torso linear acc (not used)
        getline(ss, str, ',');
        getline(ss, str, ',');
        getline(ss, str, ',');

        // torso rot acc (not used)
        getline(ss, str, ',');
        getline(ss, str, ',');
        getline(ss, str, ',');

        for (int i = 0; i < 4; i++)
        {
            // foot i pos
            getline(ss, str, ',');
            x = std::stod(str);
            getline(ss, str, ',');
            y = std::stod(str);
            getline(ss, str, ',');
            z = std::stod(str);
            foot_pos[i].push_back({x, y, z});

            // foot i vel
            getline(ss, str, ',');
            x = std::stod(str);
            getline(ss, str, ',');
            y = std::stod(str);
            getline(ss, str, ',');
            z = std::stod(str);
            foot_vel[i].push_back({x, y, z});

            // foot i acc (not used)
            getline(ss, str, ',');
            getline(ss, str, ',');
            getline(ss, str, ',');
        }

        for (int i = 0; i < 4; i++)
        {
            // foot i force
            getline(ss, str, ',');
            x = std::stod(str);
            getline(ss, str, ',');
            y = std::stod(str);
            getline(ss, str, ',');
            z = std::stod(str);
            foot_force_[i].push_back({x, y, z});
        }

        for (int i = 0; i < 4; i++)
        {
            // foot i contact
            getline(ss, str, ',');
            foot_contact_[i].push_back(std::stoi(str));
        }
    }

    traj_beg_time_ = 0;
    sample_time_interval_ = timestamps[1] - timestamps[0];
    const int len = timestamps.size() - 1;

    for (int i = 0; i < len; i++)
    {
        const double beg_time = timestamps[i];
        const double end_time = timestamps[i + 1];

        torso_trans_traj_.emplace_back(
                    CubicSpline<Eigen::Vector3d>(
                        beg_time, torso_pos[i], torso_vel[i],
                        end_time, torso_pos[i + 1], torso_vel[i + 1]));

        torso_rot_traj_.emplace_back(
                    CubicSpline<Eigen::Quaterniond, Eigen::Vector3d>(
                        beg_time, torso_rot[i], torso_rot_vel[i],
                        end_time, torso_rot[i + 1], torso_rot_vel[i + 1]));

        for (int j = 0; j < 4; j++)
        {
            foot_traj_[j].emplace_back(
                        CubicSpline<Eigen::Vector3d>(
                            beg_time, foot_pos[j][i], foot_vel[j][i],
                            end_time, foot_pos[j][i + 1], foot_vel[j][i + 1]));
        }
    }
}
