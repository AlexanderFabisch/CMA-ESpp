#pragma once

#include <string>
#include <ostream>

/**
 * @class Parameters
 * Holds all parameters that can be adjusted by the user.
 */
class Parameters
{
  friend class CMAES;
  friend std::ostream& operator<<(std::ostream& stream, const Parameters& parameters);

public:

  /* input parameter */
  int N; //!< problem dimension, must stay constant
  double* xstart;
  double* typicalX;
  bool typicalXcase;
  double* rgInitialStds;
  double* rgDiffMinChange;

  /* termination parameters */
  double stopMaxFunEvals;
  double facmaxeval;
  double stopMaxIter;

  struct { bool flg; double val; } stStopFitness;
  double stopTolFun;
  double stopTolFunHist;
  double stopTolX;
  double stopTolUpXFactor;

  /* internal evolution strategy parameters */
  int lambda; //!< population size
  int mu; //!< offspring size
  double mucov, mueff; //!< <- weights
  double* weights; //!< <- mu, -> mueff, mucov, ccov
  double damps; //!< <- cs, maxeval, lambda
  double cs; //!< -> damps, <- N
  double ccumcov; //!< <- N
  double ccov; //!< <- mucov, <- N
  double diagonalCov; //!< number of initial iterations
  struct { double modulo; double maxtime; } updateCmode;
  double facupdateCmode;

  enum Weights
  {
    UNINITIALIZED_WEIGHTS, LINEAR_WEIGHTS, EQUAL_WEIGHTS, LOG_WEIGHTS
  } weightMode;
  std::string resumefile;

  Parameters();
  Parameters(const Parameters& parameters);
  ~Parameters();
  Parameters& operator=(const Parameters& parameters);
  /*
   * @param dimension Dimension of the search space \f$N\f$. No default
   *                  available, must be defined here or you have to set the
   *                  member manually.
   * @param initialX Initial point in search space \f$x_0\f$, default (NULL) is
   *                 \f$(0.5,\ldots,0.5)^T + N(0, initialStdDev^2) \in R^N\f$.
   *                 This must be an array of size \f$N\f$.
   * @param initialStdDev Coordinatewise initial standard deviation of the
   *                      sample distribution (\f$\sigma \cdot \sqrt{C_{ii}} =
   *                      initialStdDev[i]\f$). The expected initial distance
   *                      between initialX and the optimum per coordinate should
   *                      be roughly initialStdDev. The entries should not
   *                      differ by several orders of magnitude. Default (NULL)
   *                      is \f$(0.3,\ldots,0.3)^T \in R^N\f$. This must be
   *                      an array of size \f$N\f$.
   * @param lambda Population size, number of sampled candidate solutions per
   *               generation. Default (0) is \f$4 + \lfloor 3\log{N} \rfloor\f$
   */
  void init(int dimension = 0, const double* initialX = 0,
      const double* initialStdDev = 0);

private:
  void assign(const Parameters& p);
  /**
   * Supplements default parameter values.
   */
  void supplementDefaults();
  /**
   * Initializes the offspring weights.
   */
  void setWeights(Weights mode);
};

std::ostream& operator<<(std::ostream& stream, const Parameters& parameters);

