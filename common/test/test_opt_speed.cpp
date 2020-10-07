/**
 * Test optimization algorithm speed on simulated problems.
 */

#include "../third-party/Goldfarb_Optimizer/QuadSolver.h"
#include "../third-party/qpOASES/include/qpOASES.hpp"
#include "Math/MathUtilities.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <sys/time.h>

TEST(Goldfarb_Optimizer, CompareSpeedWithQPOASES)
{
    // set a hessian matrix
    Eigen::MatrixXd H = Eigen::MatrixXd::Random(150, 150);
    H = H.transpose() * H;
    Eigen::VectorXd g = Eigen::VectorXd::Random(150);

    // span an available area by several corners
    std::array<Eigen::VectorXd, 300> x0;

    for(int i = 0; i < 300; i++)
        x0[i] = Eigen::VectorXd::Random(150) * 10;

    // set constraints. keep all x0 available to ensure
    // there exists a solution.

    Eigen::MatrixXd ci;
    Eigen::VectorXd li;
    ci = Eigen::MatrixXd::Random(300, 150);
    li.resize(300);
    // ci * x >= li

    for(int i = 0; i < 300; i++)
    {
        double min = 100000;

        for(int j = 0; j < 300; j++)
        {
            double res = ci.row(i).dot(x0[j]);

            if(min > res)
                min = res;
        }

        li(i) = min;
    }

    // use pca to compute ce, le.
    Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(150, 150);
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(150);

    // only use the first 30 elements
    for(int i = 0; i < 30; i++)
    {
        sum += x0[i];
        cov += x0[i] * x0[i].transpose();
    }

    cov -= sum * sum.transpose() / 30;
    Eigen::EigenSolver<Eigen::MatrixXd> es(cov);
    Eigen::MatrixXd vec = es.pseudoEigenvectors();
    Eigen::VectorXd val = es.pseudoEigenvalueMatrix().diagonal();

    Eigen::MatrixXd ce;
    ce.resize(25, 150);

    for(int i = 0; i < 25; i++)
    {
        int max_pos = 0;

        val.maxCoeff(&max_pos);

        val(max_pos) = -1;
        ce.row(i) = vec.col(max_pos).transpose();
    }

    Eigen::VectorXd le = ce * x0[0];

    Eigen::VectorXd res_g = Eigen::VectorXd::Zero(150);

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
//    Eigen::VectorXf res_gf(30);

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

    gettimeofday(&t2, NULL);
    timeuse = (t2.tv_sec - t1.tv_sec) * 1000000
            + (t2.tv_usec - t1.tv_usec);

    std::cout << "goldfarb elapsed: " << timeuse << "us" << std::endl;

    qpOASES::QProblem problem(150, 350);
    qpOASES::Options opt;
    opt.setToMPC();
    opt.printLevel = qpOASES::PrintLevel::PL_NONE;
    problem.setOptions(opt);

    Eigen::MatrixXd expand_A;
    expand_A.resize(350, 150);
    expand_A.topRows<300>() = ci;
    expand_A.middleRows<25>(300) = ce;
    expand_A.bottomRows<25>() = - ce;
    expand_A.transposeInPlace();
    Eigen::VectorXd expand_l;
    expand_l.resize(350);
    expand_l.segment<300>(0) = li;
    expand_l.segment<25>(300) = le;
    expand_l.segment<25>(325) = - le;

    int nWSR = 10000;
    Eigen::VectorXd res_q = Eigen::VectorXd::Zero(150);

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
            << "error: " << (res_g - res_q).segment(0, 5).transpose();
}
