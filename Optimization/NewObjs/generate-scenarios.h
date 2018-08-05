// generate scenarios of natural phosphorous inflows
// lognormal distribution, independent timesteps

#include <time.h>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;
using namespace std;

ublas::vector<double> generateSOW(int nDays, double mu, double sigma)
{
  ublas::vector<double> P_inflow(nDays);

  boost::mt19937 rng(time(NULL)); 
  boost::normal_distribution<double> nd(mu, sigma);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(rng, nd);

  for (int i = 0; i < nDays; i++)
    P_inflow(i) = var_nor();

  return P_inflow;
}
