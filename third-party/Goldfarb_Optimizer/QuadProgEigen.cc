// Eigen version of quad programming

#include "QuadProg++.hh"
void compute_d(Eigen::VectorXd &d, const Eigen::MatrixXd &J, const Eigen::VectorXd &np);
void update_z(Eigen::VectorXd& z, const Eigen::MatrixXd& J, const Eigen::VectorXd& d, int iq);
void update_r(const Eigen::MatrixXd& R, Eigen::VectorXd& r, const Eigen::VectorXd& d, int iq);
bool add_constraint(Eigen::MatrixXd& R, Eigen::MatrixXd& J, Eigen::VectorXd& d, int& iq, double& R_norm);
void delete_constraint(Eigen::MatrixXd& R, Eigen::MatrixXd& J, Eigen::VectorXi& A, Eigen::VectorXd& u, int n, int p, int& iq, int l);

double solve_quadprog(Eigen::MatrixXd& G, Eigen::VectorXd& g0,
                      const Eigen::MatrixXd& CE, const Eigen::VectorXd& ce0,
                      const Eigen::MatrixXd& CI, const Eigen::VectorXd& ci0,
                      Eigen::VectorXd& x)
{
    std::ostringstream msg;
    int n = G.cols(), p = CE.cols(), m = CI.cols();
    if (G.rows() != n)
    {
        msg << "The matrix G is not a squared matrix (" << G.rows() << " x " << G.cols() << ")";
        throw std::logic_error(msg.str());
    }
    if (CE.rows() != n)
    {
        msg << "The matrix CE is incompatible (incorrect number of rows " << CE.rows() << " , expecting " << n << ")";
        throw std::logic_error(msg.str());
    }
    if (ce0.size() != p)
    {
        msg << "The vector ce0 is incompatible (incorrect dimension " << ce0.size() << ", expecting " << p << ")";
        throw std::logic_error(msg.str());
    }
    if (CI.rows() != n)
    {
        msg << "The matrix CI is incompatible (incorrect number of rows " << CI.rows() << " , expecting " << n << ")";
        throw std::logic_error(msg.str());
    }
    if (ci0.size() != m)
    {
        msg << "The vector ci0 is incompatible (incorrect dimension " << ci0.size() << ", expecting " << m << ")";
        throw std::logic_error(msg.str());
    }

    x.resize(n);

    int i, j, k, l; /* indices */
    int ip; // this is the index of the constraint to be added to the active set

    Eigen::MatrixXd R, J;
    R.resize(n, n);

    Eigen::VectorXd s(m + p), z(n), r(m + p), d(n), np(n), u(m + p), x_old(n), u_old(m + p);
    double f_value, psi, c1, c2, ss, R_norm;

    const double inf = std::numeric_limits<double>::has_infinity ?
                std::numeric_limits<double>::infinity() : 1E300;

    double t, t1, t2; /* t is the step lenght, which is the minimum of the partial step length t1
      * and the full step length t2 */

    Eigen::VectorXi A(m + p), A_old(m + p), iai(m + p);
    int q, iq, iter = 0;
    Eigen::VectorXi iaexcl(m + p);

    /* p is the number of equality constraints */
    /* m is the number of inequality constraints */
    q = 0;  /* size of the active set A (containing the indices of the active constraints) */

    /*
     * Preprocessing phase
     */

    /* compute the trace of the original matrix G */
    c1 = G.trace();

    /* decompose the matrix G in the form L L^T */
    J = G.llt().matrixU().solve(Eigen::MatrixXd::Identity(n, n));
//    J = Eigen::MatrixXd::Identity(n, n);
    d.setZero();
    R.setZero();

    R_norm = 1.0; /* this variable will hold the norm of the matrix R */
    c2 = J.trace();

//    std::cout << "J" << std::endl << J << std::endl;

    /*
      * Find the unconstrained minimizer of the quadratic form 0.5 * x G x + g0 x
     * this is a feasible point in the dual space
     * x = - G^-1 * g0
     */
    x = - J * J.transpose() * g0;
    f_value = g0.dot(x) / 2;

    /* Add equality constraints to the working set A */
    iq = 0;
    for (i = 0; i < p; i++)
    {
      np = CE.col(i);
      compute_d(d, J, np);
      update_z(z, J, d, iq);
      update_r(R, r, d, iq);

//      std::cout << "iq " << iq << std::endl;
//      std::cout << "R" << std::endl << R << std::endl;
//      std::cout << "z " << z.transpose() << std::endl;
//      std::cout << "r " << r.transpose() << std::endl;
//      std::cout << "d " << d.transpose() << std::endl;

      /* compute full step length t2: i.e., the minimum step in primal space s.t. the contraint
        becomes feasible */
      t2 = 0.0;
      if (z.norm() > std::numeric_limits<double>::epsilon()) // i.e. z != 0
        t2 = (- np.dot(x) - ce0[i]) / np.dot(z);

      /* set x = x + t2 * z */
      x += t2 * z;

      /* set u = u+ */
      u[iq] = t2;
      u.topRows(iq) -= t2 * r.topRows(iq);

      /* compute the new solution value */
      f_value += 0.5 * (t2 * t2) * np.dot(z);
      A[i] = - i - 1;

      if (!add_constraint(R, J, d, iq, R_norm))
      {
        // Equality constraints are linearly dependent
        throw std::runtime_error("Constraints are linearly dependent");
        return f_value;
      }
    }

    /* set iai = K \ A */
    for (i = 0; i < m; i++)
      iai[i] = i;

l1:	iter++;

//    std::cout << "x: " << x.transpose() << std::endl << std::endl;

    /* step 1: choose a violated constraint */
    for (i = p; i < iq; i++)
    {
      ip = A[i];
      iai[ip] = -1;
    }

    /* compute s[x] = ci^T * x + ci0 for all elements of K \ A */
    ss = 0.0;
    psi = 0.0; /* this value will contain the sum of all infeasibilities */
    ip = 0; /* ip will be the index of the chosen violated constraint */

    iaexcl.topRows(m).setOnes();
    s.segment(0, m) = CI.transpose() * x + ci0;
    psi = s.cwiseMin(0).sum();

    if (fabs(psi) <= m * std::numeric_limits<double>::epsilon() * c1 * c2 * 100.)
    {
      /* numerically there are not infeasibilities anymore */
      q = iq;

      return f_value;
    }

    /* save old values for u and A */
    u_old.segment(0, iq) = u.segment(0, iq);
    A_old.segment(0, iq) = A.segment(0, iq);
    /* and for x */
    x_old = x;

l2: /* Step 2: check for feasibility and determine a new S-pair */
    for (i = 0; i < m; i++)
    {
      if (s[i] < ss && iai[i] != -1 && iaexcl[i])
      {
        ss = s[i];
        ip = i;
      }
    }

    if (ss >= 0.0)
    {
        q = iq;

        return f_value;
    }

    /* set np = n[ip] */
    np = CI.col(ip);
    /* set u = [u 0]^T */
    u[iq] = 0.0;
    /* add ip to the active set A */
    A[iq] = ip;

//    std::cout << "Trying with constraint " << ip << std::endl;
//    std::cout << "np " << np.transpose() << std::endl;

l2a:/* Step 2a: determine step direction */
    /* compute z = H np: the step direction in the primal space (through J, see the paper) */
    compute_d(d, J, np);
    update_z(z, J, d, iq);
    /* compute N* np (if q > 0): the negative of the step direction in the dual space */
    update_r(R, r, d, iq);

    /* Step 2b: compute step length */
    l = 0;
    /* Compute t1: partial step length (maximum step in dual space without violating dual feasibility */
    t1 = inf; /* +inf */
    /* find the index l s.t. it reaches the minimum of u+[x] / r */
    for (k = p; k < iq; k++)
    {
      if (r[k] > 0.0)
      {
        if (u[k] / r[k] < t1)
          {
            t1 = u[k] / r[k];
            l = A[k];
          }
      }
    }

    /* Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible */
    if (fabs(z.norm()) > std::numeric_limits<double>::epsilon()) // i.e. z != 0
    {
      t2 = -s[ip] / np.dot(z);
      if (t2 < 0) // patch suggested by Takano Akio for handling numerical inconsistencies
        t2 = inf;
    }
    else
      t2 = inf; /* +inf */

    /* the step is chosen as the minimum of t1 and t2 */
    t = std::min(t1, t2);

//    std::cout << "Step sizes: " << t << " (t1 = " << t1 << ", t2 = " << t2 << ") ";
    /* Step 2c: determine new S-pair and take step: */

    /* case (i): no step in primal or dual space */
    if (t >= inf)
    {
      /* QPP is infeasible */
      // FIXME: unbounded to raise
      q = iq;
      return inf;
    }

    /* case (ii): step in dual space */
    if (t2 >= inf)
    {
      /* set u = u +  t * [-r 1] and drop constraint l from the active set A */
      for (k = 0; k < iq; k++)
        u[k] -= t * r[k];
      u[iq] += t;
      iai[l] = l;
      delete_constraint(R, J, A, u, n, p, iq, l);
      goto l2a;
    }

    /* case (iii): step in primal and dual space */

    /* set x = x + t * z */
    x += t * z;
    /* update the solution value */
    f_value += t * np.dot(z) * (0.5 * t + u[iq]);
    /* u = u + t * [-r 1] */
    u.segment(0, iq) -= t * r.segment(0, iq);
    u[iq] += t;

    if (fabs(t - t2) < std::numeric_limits<double>::epsilon())
    {
      /* full step has taken */
      /* add constraint ip to the active set*/
      if (!add_constraint(R, J, d, iq, R_norm))
      {
        iaexcl[ip] = false;
        delete_constraint(R, J, A, u, n, p, iq, ip);
//        std::cout << "R" << std::endl << R << std::endl;
//        std::cout << "A " << A.transpose() << std::endl;
        for (i = 0; i < m; i++)
          iai[i] = i;

        A.segment(p, iq - p) = A_old.segment(p, iq - p);
        u.segment(p, iq - p) = u_old.segment(p, iq - p);

        for (i = p; i < iq; i++)
        {
            iai[A[i]] = -1;
        }

        x = x_old;

        goto l2; /* go to step 2 */
      }
      else
        iai[ip] = -1;
      goto l1;
    }

    /* a patial step has taken */
    /* drop constraint l */
    iai[l] = l;
    delete_constraint(R, J, A, u, n, p, iq, l);

    /* update s[ip] = CI^T * x + ci0 */
    s[ip] = CI.col(ip).dot(x) + ci0[ip];

    goto l2a;
}

inline void compute_d(Eigen::VectorXd& d, const Eigen::MatrixXd& J, const Eigen::VectorXd& np)
{
  /* compute d = H^T * np */
  d = J.transpose() * np;
}

inline void update_z(Eigen::VectorXd& z, const Eigen::MatrixXd& J, const Eigen::VectorXd& d, int iq)
{
    const int sz = d.size() - iq;
  /* setting of z = H * d */
  z = J.rightCols(sz) * d.bottomRows(sz);
}

inline void update_r(const Eigen::MatrixXd& R, Eigen::VectorXd& r, const Eigen::VectorXd& d, int iq)
{
  /* setting of r = R^-1 d */
    r.topRows(iq) = R.topLeftCorner(iq, iq).triangularView<Eigen::Upper>().solve(d.topRows(iq));
}

bool add_constraint(Eigen::MatrixXd& R, Eigen::MatrixXd& J, Eigen::VectorXd& d, int& iq, double& R_norm)
{
  int n = d.size();
  int i, j, k;
  double cc, ss, h, t1, t2, xny;
//  std::cout << "Add constraint " << iq << '/';

  /* we have to find the Givens rotation which will reduce the element
    d[j] to zero.
    if it is already zero we don't have to do anything, except of
    decreasing j */
  for (j = n - 1; j >= iq + 1; j--)
  {
    /* The Givens rotation is done with the matrix (cc cs, cs -cc).
    If cc is one, then element (j) of d is zero compared with element
    (j - 1). Hence we don't have to do anything.
    If cc is zero, then we just have to switch column (j) and column (j - 1)
    of J. Since we only switch columns in J, we have to be careful how we
    update d depending on the sign of gs.
    Otherwise we have to apply the Givens rotation to these columns.
    The i - 1 element of d has to be updated to h. */
    cc = d[j - 1];
    ss = d[j];
    h = sqrt(cc * cc + ss * ss);
    if (h < std::numeric_limits<double>::epsilon()) // h == 0
      continue;
    d[j] = 0.0;
    ss = ss / h;
    cc = cc / h;
    if (cc < 0.0)
    {
      cc = -cc;
      ss = -ss;
      d[j - 1] = -h;
    }
    else
      d[j - 1] = h;
    xny = ss / (1.0 + cc);

    const Eigen::VectorXd t1 = J.col(j - 1);
    const Eigen::VectorXd t2 = J.col(j);
    J.col(j - 1) = cc * t1 + ss * t2;
    J.col(j) = ss * t1 - cc * t2;
  }
  /* update the number of constraints added*/
  iq++;
  /* To update R we have to put the iq components of the d vector
    into column iq - 1 of R
    */
  R.block(0, iq - 1, iq, 1) = d.topRows(iq);

//  std::cout << iq << std::endl;

  if (fabs(d[iq - 1]) <= std::numeric_limits<double>::epsilon() * R_norm)
  {
    // problem degenerate
    return false;
  }
  R_norm = std::max<double>(R_norm, fabs(d[iq - 1]));
  return true;
}

void delete_constraint(Eigen::MatrixXd& R, Eigen::MatrixXd& J, Eigen::VectorXi& A, Eigen::VectorXd& u, int n, int p, int& iq, int l)
{
  int i, j, k, qq = -1; // just to prevent warnings from smart compilers
  double cc, ss, h, xny, t1, t2;

  /* Find the index qq for active constraint l to be removed */
  for (i = p; i < iq; i++)
    if (A[i] == l)
    {
      qq = i;
      break;
    }

  /* remove the constraint from the active set and the duals */
  for (i = qq; i < iq - 1; i++)
    {
      A[i] = A[i + 1];
      u[i] = u[i + 1];
      R.col(i) = R.col(i + 1);
    }

  A[iq - 1] = A[iq];
  u[iq - 1] = u[iq];
  A[iq] = 0;
  u[iq] = 0.0;
  R.col(iq - 1).setZero();
  /* constraint has been fully removed */
  iq--;

  if (iq == 0)
    return;

  for (j = qq; j < iq; j++)
  {
    cc = R(j, j);
    ss = R(j + 1, j);
    h = sqrt(cc * cc + ss * ss);
    if (h < std::numeric_limits<double>::epsilon()) // h == 0
      continue;
    cc = cc / h;
    ss = ss / h;
    R(j + 1, j) = 0.0;
    if (cc < 0.0)
    {
      R(j, j) = -h;
      cc = -cc;
      ss = -ss;
    }
    else
      R(j, j) = h;

    xny = ss / (1.0 + cc);

    Eigen::RowVectorXd t1 = R.block(j, j + 1, 1, iq - j - 1);
    Eigen::RowVectorXd t2 = R.block(j + 1, j + 1, 1, iq - j - 1);
    R.block(j, j + 1, 1, iq - j - 1) = cc * t1 + ss * t2;
    R.block(j + 1, j + 1, 1, iq - j - 1) = ss * t1 - cc * t2;

    t1 = J.col(j);
    t2 = J.col(j + 1);
    J.col(j) = cc * t1 + ss * t2;
    J.col(j + 1) = ss * t1 - cc * t2;
  }
}
