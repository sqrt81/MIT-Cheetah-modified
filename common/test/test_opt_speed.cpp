/**
 * Test optimization algorithm speed on simulated problems.
 */

#include "../third-party/Goldfarb_Optimizer/QuadSolver.h"
#include "../third-party/Goldfarb_Optimizer/QuadProg++.hh"
#include "../third-party/qpOASES/include/qpOASES.hpp"
#include "Math/MathUtilities.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <sys/time.h>

TEST(Goldfarb_Optimizer, CompareSpeedWithQPOASES)
{
    constexpr int n_v = 30;
    constexpr int n_ce = 9;
    constexpr int n_ce_sel = 10;
    constexpr int n_ci = 60;
    // set a hessian matrix
    Eigen::MatrixXd H = Eigen::MatrixXd::Random(n_v, n_v);
    H = H.transpose() * H;
    Eigen::VectorXd g = Eigen::VectorXd::Random(n_v);

    // span an available area by several corners
    std::array<Eigen::VectorXd, n_ci> x0;

    for(int i = 0; i < n_ci; i++)
        x0[i] = Eigen::VectorXd::Random(n_v) * 10;

    // set constraints. keep all x0 available to ensure
    // there exists a solution.

    Eigen::MatrixXd ci;
    Eigen::VectorXd li;
    ci = Eigen::MatrixXd::Random(n_ci, n_v);
    li.resize(n_ci);
    // ci * x >= li

    for(int i = 0; i < n_ci; i++)
    {
        double min = 100000;

        for(int j = 0; j < n_ci; j++)
        {
            double res = ci.row(i).dot(x0[j]);

            if(min > res)
                min = res;
        }

        li(i) = min;
    }

    // use pca to compute ce, le.
    Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(n_v, n_v);
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(n_v);

    // only use the first n_ce_sel elements
    for(int i = 0; i < n_ce_sel; i++)
    {
        sum += x0[i];
        cov += x0[i] * x0[i].transpose();
    }

    cov -= sum * sum.transpose() / n_ce_sel;
    Eigen::EigenSolver<Eigen::MatrixXd> es(cov);
    Eigen::MatrixXd vec = es.pseudoEigenvectors();
    Eigen::VectorXd val = es.pseudoEigenvalueMatrix().diagonal();

    Eigen::MatrixXd ce;
    ce.resize(n_ce, n_v);

    for(int i = 0; i < n_ce; i++)
    {
        int max_pos = 0;

        val.maxCoeff(&max_pos);

        val(max_pos) = -1;
        ce.row(i) = vec.col(max_pos).transpose();
    }

    Eigen::VectorXd le = ce * x0[0];

    Eigen::VectorXd res_g = Eigen::VectorXd::Zero(n_v);

    struct timeval t1, t2;
    int timeuse;
//    (void) timeuse;
//    (void) t2;

//    solve_quadprog(H, g, ce.transpose(), - le,
//                   ci.transpose(), - li, res_g);

//    Eigen::MatrixXf Hf = H.cast<float>();
//    Eigen::VectorXf gf = g.cast<float>();
//    Eigen::MatrixXf cef = ce.cast<float>();
//    Eigen::VectorXf lef = le.cast<float>();
//    Eigen::MatrixXf cif = ci.cast<float>();
//    Eigen::VectorXf lif = li.cast<float>();
//    Eigen::VectorXf res_gf(n_v);

    gettimeofday(&t1, NULL);
//    (void) timeuse;
//    (void) t2;

    solve_quadprog(H, g,
                   ce.transpose(), - le,
                   ci.transpose(), - li,
                   res_g);
//    res_g = res_gf.cast<double>();

    gettimeofday(&t2, NULL);
    timeuse = (t2.tv_sec - t1.tv_sec) * 1000000
            + (t2.tv_usec - t1.tv_usec);

    std::cout << "origin goldfarb elapsed: " << timeuse << "us" << std::endl;

    gettimeofday(&t1, NULL);
//    (void) timeuse;
//    (void) t2;

    QuadSolver::SolveQuadProg(H, g,
                   ce.transpose(), - le,
                   ci.transpose(), - li,
                   res_g);
//    res_g = res_gf.cast<double>();

    gettimeofday(&t2, NULL);
    timeuse = (t2.tv_sec - t1.tv_sec) * 1000000
            + (t2.tv_usec - t1.tv_usec);

    std::cout << "goldfarb elapsed: " << timeuse << "us" << std::endl;

    qpOASES::QProblem problem(n_v, n_ci + 2 * n_ce);
    qpOASES::Options opt;
    opt.setToMPC();
    opt.printLevel = qpOASES::PrintLevel::PL_NONE;
    problem.setOptions(opt);

    Eigen::MatrixXd expand_A;
    expand_A.resize(n_ci + 2 * n_ce, n_v);
    expand_A.topRows<n_ci>() = ci;
    expand_A.middleRows<n_ce>(n_ci) = ce;
    expand_A.bottomRows<n_ce>() = - ce;
    expand_A.transposeInPlace();
    Eigen::VectorXd expand_l;
    expand_l.resize(n_ci + 2 * n_ce);
    expand_l.segment<n_ci>(0) = li;
    expand_l.segment<n_ce>(n_ci) = le;
    expand_l.segment<n_ce>(n_ci + n_ce) = - le;

    int nWSR = 10000;
    Eigen::VectorXd res_q = Eigen::VectorXd::Zero(n_v);

    gettimeofday(&t1, NULL);

    int stat = problem.init(H.data(), g.data(), expand_A.data(),
                 nullptr, nullptr, expand_l.data(), nullptr, nWSR);
    problem.getPrimalSolution(res_q.data());

    gettimeofday(&t2, NULL);
    timeuse = (t2.tv_sec - t1.tv_sec) * 1000000
            + (t2.tv_usec - t1.tv_usec);

    std::cout << "qp status: "
              << ((stat == qpOASES::SUCCESSFUL_RETURN) ?
                     "successful" : "failed") << std::endl;
    std::cout << "qpOASES elapsed: " << timeuse << "us" << std::endl;

    EXPECT_TRUE(almostEqual(res_g, res_q, 1e-4))
            << "g: " << res_g.segment(0, 5).transpose() << std::endl
            << "q: " << res_q.segment(0, 5).transpose() << std::endl
            << "error: " << (res_g - res_q).norm();
}
