/**
 * @file cmaes.h
 * @author Nikolaus Hansen, ported to C++ by Alexander Fabisch
 *
 * \mainpage
 * CMA-ES for non-linear function minimization.
 *
 * Copyright of C implementation by Nikolaus Hansen (e-mail:
 * hansen .AT. bionik.tu-berlin.de, hansen .AT. lri.fr), ported to C++ by
 * <a href="mailto:afabisch@googlemail.com"> Alexander Fabisch</a>.
 *
 * \section lgpl License
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License, version 2,
 * as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * \section purpose General Purpose
 *
 * The CMA-ES (Evolution Strategy with Covariance Matrix Adaptation) is a
 * robust search/optimization method. The goal is to minimize a given
 * objective function, \f$f: R^N \rightarrow R\f$. The CMA-ES should be
 * applied, if e.g. BFGS and/or conjugate gradient methods fail due to a
 * rugged search landscape (e.g. discontinuities, outliers, noise, local
 * optima, etc.). Learning the covariance matrix in the CMA-ES is similar
 * to learning the inverse Hessian matrix in a quasi-Newton method. On
 * smooth landscapes the CMA-ES is roughly ten times slower than BFGS,
 * assuming derivatives are not directly available. For up to \f$N=10\f$
 * parameters the simplex direct search method (Nelder & Mead) is
 * sometimes faster, but less robust than CMA-ES.  On considerably hard
 * problems the search (a single run) is expected to take between
 * \f$100\cdot N\f$ and \f$300\cdot N^2\f$ function evaluations. But you
 * might be lucky...
 *
 * \section application Application Remark
 *
 * The adaptation of the covariance matrix (e.g. by the CMA) is
 * equivalent to a general linear transformation of the problem
 * variables. Nevertheless, every problem specific knowledge about the
 * best problem transformation should be exploited before starting the
 * search procedure and an appropriate a priori transformation should be
 * applied to the problem. In particular a decision should be taken
 * whether variables, which are positive by nature, should be taken in
 * the log scale. A hard lower variable bound can also be realized by
 * taking the square. All variables should be re-scaled such that they
 * "live" in a similar search range width (for example, but not
 * necessarily between zero and one), such that the initial standard
 * deviation can be chosen the same for all variables.
 *
 *
 * \section links Links
 *  - http://www.lri.fr/~hansen/cmaesintro.html
 *  - http://www.lri.fr/~hansen/publications.html
 *
 * \section tut Tutorial
 * 	- http://www.lri.fr/~hansen/cmatutorial.pdf
 *
 * \section references References
 *
 * - Hansen, N, and S. Kern (2004).  Evaluating the CMA Evolution
 *   Strategy on Multimodal Test Functions. In: Eighth International
 *   Conference on Parallel Problem Solving from Nature PPSN VIII,
 *   Proceedings, pp. 282-291, Berlin: Springer
 *
 * - Hansen, N., S.D. Müller and P. Koumoutsakos (2003): Reducing the
 *   Time Complexity of the Derandomized Evolution Strategy with
 *   Covariance Matrix Adaptation (CMA-ES). Evolutionary Computation,
 *   11(1).
 *
 * - Hansen, N. and A. Ostermeier (2001). Completely Derandomized
 *   Self-Adaptation in Evolution Strategies. Evolutionary Computation,
 *   9(2), pp. 159-195.
 *
 * - Hansen, N. and A. Ostermeier (1996). Adapting arbitrary normal
 *   mutation distributions in evolution strategies: The covariance
 *   matrix adaptation. In Proceedings of the 1996 IEEE International
 *   Conference on Evolutionary Computation, pp. 312-317.
 */

#pragma once

#include "random.h"
#include "timings.h"
#include "parameters.h"
#include <string>

/**
 * @class CMAES
 * Evolution Strategies with Covariance Matrix Adaptation. The public interface
 * of the optimization algorithm.
 */
class CMAES
{
public:

  /**
   * Keys for get().
   */
  enum GetScalar
  {
    NoScalar = 0,
    AxisRatio = 1,
    Eval = 2, Evaluations = 2,
    FctValue = 3, FuncValue = 3, FunValue = 3, Fitness = 3,
    FBestEver = 4,
    Generation = 5, Iteration = 5,
    MaxEval = 6, MaxFunEvals = 6, StopMaxFunEvals = 6,
    MaxGen = 7, MaxIter = 7, StopMaxIter = 7,
    MaxAxisLength = 8,
    MinAxisLength = 9,
    MaxStdDev = 10,
    MinStdDev = 11,
    Dim = 12, Dimension = 12,
    Lambda = 13, SampleSize = 13, PopSize = 13,
    Sigma = 14
  };

  /**
   * Keys for getPtr().
   */
  enum GetVector
  {
    NoVector = 0,
    DiagC = 1,
    DiagD = 2,
    StdDev = 3,
    XBestEver = 4,
    XBest = 5,
    XMean = 6
  };

  /**
   * Keys for writeToFile().
   */
  enum WriteKey
  {
    WCNone = 0,
    WKResume = 1,
    WKXMean = 2,
    WKC = 4,
    WKAll = 8,
    WKFewInfo = 16,
    WKFew = 32,
    WKEval = 64,
    WKFitness = 128,
    WKFBestEver = 256,
    WKCGeneration = 512,
    WKSigma = 1024,
    WKLambda = 2048,
    WKB = 4096,
    WKXBest = 8192,
    WKClock = 16384,
    WKDim = 32768
  };

private:

  std::string version; //!< Implementation version.
  Random rand; //!< Random number generator.
  Parameters sp; //!< CMA-ES parameters.

  double sigma; //!< step size
  double* xmean; //!< mean x vector, "parent"
  double* rgxbestever; //!< best sample ever
  double** rgrgx; //!< x-vectors, lambda offspring
  int* index; //!< sorting index of sample pop.
  double* arFuncValueHist;

  double chiN;
  double** C; //!< lower triangular matrix: i>=j for C[i][j]
  double** B; //!< matrix with normalize eigenvectors in columns
  double* rgD; //!< axis lengths

  double* rgpc;
  double* rgps;
  double* rgxold;
  double* rgout;
  double* rgBDz; //!< for B*D*z
  double* rgdTmp; //!< temporary (random) vector used in different places
  double* rgFuncValue;
  double* publicFitness; //!< returned by init()

  double gen; //!< Generation number
  enum {INITIALIZED, SAMPLED, UPDATED} state; //!< algorithm state

  double maxdiagC; //!< repeatedly used for output
  double mindiagC;
  double maxEW;
  double minEW;

  bool eigensysIsUptodate;
  bool doCheckEigen; //!< control via signals.par
  double genOfEigensysUpdate;

  double dMaxSignifKond;
  double dLastMinEWgroesserNull;

  bool isResumeDone;

  time_t printtime;
  time_t writetime; //!< ideally should keep track for each output file
  time_t firstwritetime;
  time_t firstprinttime;

  std::string stopMessage; //!< A message that contains all matched stop criteria.

  std::string getTimeStr(void);

  /** 
   * Calculating eigenvalues and vectors.
   * @param rgtmp (input) N+1-dimensional vector for temporal use. 
   * @param diag (output) N eigenvalues. 
   * @param Q (output) Columns are normalized eigenvectors.
   */
  void eigen(double* diag, double** Q, double* rgtmp);

  /** 
   * Exhaustive test of the output of the eigendecomposition, needs O(n^3)
   * operations writes to error file.
   * @return number of detected inaccuracies
   */
  int checkEigen(double* diag, double** Q);

  /**
   * Symmetric tridiagonal QL algorithm, iterative.
   * Computes the eigensystem from a tridiagonal matrix in roughtly 3N^3 operations
   * code adapted from Java JAMA package, function tql2.
   * @param d input: Diagonale of tridiagonal matrix. output: eigenvalues.
   * @param e input: [1..n-1], off-diagonal, output from Householder
   * @param V input: matrix output of Householder. output: basis of
   *          eigenvectors, according to d
   */
  void ql(double* d, double* e, double** V);

  /**
   * Householder transformation of a symmetric matrix V into tridiagonal form.
   * Code slightly adapted from the Java JAMA package, function private tred2().
   * @param V input: symmetric nxn-matrix. output: orthogonal transformation
   *          matrix: tridiag matrix == V* V_in* V^t.
   * @param d output: diagonal
   * @param e output: [0..n-1], off diagonal (elements 1..n-1)
   */
  void householder(double** V, double* d, double* e);

  /**
   * Dirty index sort.
   */
  void sortIndex(const double* rgFunVal, int* index, int n);

  void adaptC2(const int hsig);

  /**
   * Treats minimal standard deviations and numeric problems. Increases sigma.
   */
  void testMinStdDevs(void);

  /**
   * Adds the mutation sigma*B*(D*z).
   * @param x Search space vector.
   * @param eps Mutation factor.
   */
  void addMutation(double* x, double eps = 1.0);

  /**
   * This hack reads key words from input key for data to be written to
   * a file, see file signals.par as input file. The length of the keys
   * is mostly fixed. If the key phrase does not match the expectation the
   * output might be strange.
   */
  void writeToStream(int key, std::ostream& file);

public:

  double countevals; //!< objective function evaluations
  Timing eigenTimings;

  /**
   * Releases the dynamically allocated memory, including that of the return
   * value of init().
   */
  ~CMAES();

  /**
   * Initializes the CMA-ES algorithm.
   * @param parameters The CMA-ES parameters.
   * @return Array of size lambda that can be used to assign fitness values and
   *         pass them to updateDistribution(). Not that after the desctructor
   *         was called, the array is deleted.
   */
  double* init(const Parameters& parameters);

  /**
   * Well, says hello.
   * @return eg. "(5,10)-CMA-ES(mu_eff=3.4), Ver="1.0alpha", dimension=9"
   */
  std::string sayHello();

  /**
   * Allows to restart with saved internal state (distribution) variables (use
   * writeToFile() for saving). Keyword "resume" followed by a filename in
   * initials.par invokes this function during initialization. Searches in
   * filename for the last occurrence of word "resume", followed by a dimension
   * number, and reads the subsequent values for xmean, evolution paths ps and
   * pc, sigma and covariance matrix.  Note that init() needs to be called
   * before calling resume_distribution() explicitely.  In the former all the
   * remaining (strategy-)parameters are set. It can be useful to edit the
   * written parameters, in particular to increase sigma, before resume.
   *
   * Not all internal state parameters are recovered. In particular generation
   * number and xbestever are not restored. For covariance matrices with large
   * condition numbers the writing precision of 6 digits is not sufficient and
   * resume will lead to poor result.
   * @param filename A file, that was written presumably by writeToFile().
   */
  void resumeDistribution(const std::string& filename);

  /**
   * The search space vectors have a special form: they are arrays with N+1
   * entries. Entry number -1 is the dimension of the search space N.
   * @return A pointer to a "population" of lambda N-dimensional multivariate
   * normally distributed samples.
   */
  double* const* samplePopulation();

  /**
   * Can be called after samplePopulation() to resample single solutions of the
   * population as often as desired. Useful to implement a box constraints
   * (boundary) handling.
   * @param i Index to an element of the returned value of samplePopulation().
   *          population[index] will be resampled where \f$0\leq i<\lambda\f$
   *          must hold.
   * @return A pointer to the resampled "population".
   */
  double* const* reSampleSingle(int i);

  /**
   * Can be called after samplePopulation() to resample single solutions. In
   * general, the function can be used to sample as many independent
   * mean+sigma*Normal(0,C) distributed vectors as desired.
   *
   * Input x can be a pointer to an element of the vector returned by
   * samplePopulation() but this is inconsistent with the const qualifier of the
   * returned value and therefore rather reSampleSingle() should be used.
   * @param x Solution vector that gets sampled a new value. If x == NULL new
   *          memory is allocated and must be released by the user using
   *          delete[].
   * @return A pointer to the resampled solution vector, equals input x for
   *         x != NULL on input.
   */
  double* sampleSingleInto(double* x);

  /**
   * Can be called after samplePopulation() to resample single solutions. In
   * general, the function can be used to sample as many independent
   * mean+sigma*Normal(0,C) distributed vectors as desired.
   * @param x Element of the return value of samplePopulation(), that is
   *          pop[0..\f$\lambda\f$]. This solution vector of the population gets
   *          sampled a new value.
   * @return A pointer to the resampled "population" member.
   */
  double const* reSampleSingleOld(double* x);

  /**
   * Used to reevaluate a slightly disturbed solution for an uncertaintly
   * measurement. In case if x == NULL on input, the memory of the returned x
   * must be released.
   * @param x Solution vector that gets sampled a new value. If x == NULL new
   *          memory is allocated and must be released by the user using
   *          delete[] x.
   * @param xmean Mean vector \f$\mu\f$ for perturbation.
   * @param eps Scale factor \f$\epsilon\f$ for perturbation:
   *            \f$x \sim \mu + \epsilon \sigma N(0,C)\f$.
   * @return A pointer to the perturbed solution vector, equals input x for
   *         x != NULL.
   */
  double* perturbSolutionInto(double* x, double const* xmean, double eps);

  /**
   * Core procedure of the CMA-ES algorithm. Sets a new mean value and estimates
   * the new covariance matrix and a new step size for the normal search
   * distribution.
   * @param fitnessValues An array of \f$\lambda\f$ function values.
   * @return Mean value of the new distribution.
   */
  double* updateDistribution(const double* fitnessValues);

  /**
   * Request a scalar parameter from CMA-ES.
   * @param key Key of the requested scalar.
   * @return The desired value.
   */
  double get(GetScalar key);

  /**
   * Request a vector parameter from CMA-ES.
   * @param key Key of the requested vector.
   * @return Pointer to the desired value array. Its content might be
   *         overwritten during the next call to any member functions other
   *         than get().
   */
  const double* getPtr(GetVector key);

  /**
   * Request a vector parameter from CMA-ES.
   * @param key Key of the requested vector.
   * @return Pointer to the desired value array with unlimited reading and
   *         writing access to its elements. The memory must be explicitly
   *         released using delete[].
   */
  double* getNew(GetVector key);

  /**
   * Request a vector parameter from CMA-ES.
   * @param key Key of the requested vector.
   * @param res Memory of size N == dimension, where the desired values are
   *            written into. For mem == NULL new memory is allocated as with
   *            calling getNew() and must be released by the user at some point.
   */
  double* getInto(GetVector key, double* res);

  /**
   * Some stopping criteria can be set in initials.par, with names starting
   * with stop... Internal stopping criteria include a maximal condition number
   * of about 10^15 for the covariance matrix and situations where the numerical
   * discretisation error in x-space becomes noticeably. You can get a message
   * that contains the matched stop criteria via getStopMessage().
   * @return Does any stop criterion match?
   */
  bool testForTermination();

  /**
   * A message that contains a detailed description of the matched stop
   * criteria.
   */
  std::string getStopMessage();

  /**
   * @param filename Output file name.
   * @param key Key of type WriteKey that indicates the content that should be
   *            written. You can combine multiple keys with |.
   */
  void writeToFile(int key, const std::string& filename);

  /**
   * Conducts the eigendecomposition of C into B and D such that
   * \f$C = B \cdot D \cdot D \cdot B^T\f$ and \f$B \cdot B^T = I\f$
   * and D diagonal and positive.
   * @param force For force == true the eigendecomposion is conducted even if
   *              eigenvector and values seem to be up to date.
   */
  void updateEigensystem(bool force);

  /**
   * Distribution mean could be changed before samplePopulation(). This might
   * lead to unexpected behaviour if done repeatedly.
   * @param xmean new mean, if it is NULL, it will be set to the current mean
   * @return new mean
   */
  double const* setMean(const double *xmean);
};
