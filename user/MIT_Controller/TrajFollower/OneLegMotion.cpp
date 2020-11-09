#include "OneLegMotion.h"
#include "Math.h"

void OneLegMotion::SetMotion(int leg, double beg, double end,
                             bool beg_contact, bool end_contact,
                             const Eigen::Vector3d &beg_pos,
                             const Eigen::Vector3d &end_pos)
{
    leg_id_ = leg;
    beg_time_ = beg;
    end_time_ = end;
    beg_contact_ = beg_contact;
    end_contact_ = end_contact;
    beg_pos_ = beg_pos;
    target_local_pos_ = end_pos;

    if ((beg_pos.z() > 0) ^ (end_pos.z() > 0))
    {
        critic_time_ = beg_time_
                + beg_pos.z() / (end_pos.z() - beg_pos.z())
                * (end_time_ - beg_time_);

        switch (leg_id_)
        {
        case 0:
            critic_offset_ = {0.1, - 0.1, 0.};
            break;
        case 1:
            critic_offset_ = {0.1, 0.1, 0.};
            break;
        case 2:
            critic_offset_ = {- 0.1, - 0.1, 0.};
            break;
        case 3:
            critic_offset_ = {- 0.1, 0.1, 0.};
            break;
        }
    }
    else
    {
        critic_time_ = end_time_ + 1.;
    }
}

void OneLegMotion::Sample(double t,
                          Eigen::Vector3d &pos,
                          Eigen::Vector3d &vel,
                          bool &contact) const
{
    if (t < beg_time_)
    {
        pos = beg_pos_;
        vel.setZero();
        contact = beg_contact_;
        return;
    }

    if (t > end_time_)
    {
        pos = target_local_pos_;
        vel.setZero();
        contact = end_contact_;
        return;
    }

    const double coef1 = (t - beg_time_) / (end_time_ - beg_time_);
    const double coef2 = (end_time_ - t) / (end_time_ - beg_time_);

    pos = target_local_pos_ * coef1 + beg_pos_ * coef2;
    vel = (target_local_pos_ - beg_pos_) / (end_time_ - beg_time_);

    if (t > critic_time_)
    {
        pos += critic_offset_
                * ((end_time_ - t) / (end_time_ - critic_time_));
        vel -= critic_offset_ / (end_time_ - critic_time_);
    }
    else
    {
        pos += critic_offset_
                * ((t - beg_time_) / (critic_time_ - beg_time_));
        vel += critic_offset_ / (critic_time_ - beg_time_);
    }

    contact = false;
}

Eigen::Vector3d OneLegMotion::GetTorsoOffset(double t) const
{
    t = clamp(t, beg_time_, end_time_);

    const double coef1 = (t - beg_time_) / (end_time_ - beg_time_);
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
