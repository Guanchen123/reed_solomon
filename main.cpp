/*
 * -----------------------------------------------------------------------------
 * -----                               main.cpp                            -----
 * -----                          REED-SOLOMON CODES                       -----
 * -----------------------------------------------------------------------------
 *
 * File Description:
 *   This is the simulation file file for `reedSolomon` encoder/decoder
 *
 * Assumptions:
 *   None
 *
 * References:
 *   - http://downloads.bbc.co.uk/rd/pubs/whp/whp-pdf-files/WHP031.pdf
 *
 * Revision History
 *   Jun 02, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <string>
#include <sstream>
using namespace std;

#include "primitives.h" // primitive elements
#include "reedSolomon.h"

const int num_data = 15;            // number of different channel samples
const int num_trials_per_pt = 10000; // number of trials at each data point

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
/*
 * main
 * Description:
 *   main loop
 *
 * Arguments:
 *   none
 *
 * Return:
 *   Program exit status
 *
 * Operation:
 *   - initialize rand()'s seed. I do this because I like consistency in results
 *
 *   - pick default m,t values or prompt user to provide these values.
 *
 *   - For the purpose of this simulation, I determined the performance of the
 *     performance of the Reed-Solomon codes for two different choices of params
 *     + m=7, t=60
 *     + m-7, t=30
 *
 *   - error check whatever values of m,t are to be used
 *     + ensure m is less than number of bits in a system int. We want our
 *       Galois Field elements to always fit in an int.
 *     + ensure k (== n-2t == 2^m -1 -2t) > 0
 *
 */

int main(void)
{
  srand(time(0)); // initialize random seed

  unsigned int m, t;
  int k;

  // pick default m,t values
  m = 8;  // probably want values < 16 for top speed
  t = 16; // 30                   // remember n = 2^m-1, so pick t accordingly

  // crude error checking
  // check that m <= number of bits in an int
  if (m > sizeof(int) * 8)
  {
    cout << "m (== " << m << ") has to be <= int bit count of "
         << sizeof(int) * 8 << endl;
    exit(0);
  }

  // check that k > 0
  k = (pow(2, m) - 1 - 2 * t);
  if (k <= 0)
  {
    cout << "k (== n-2t == 2^m -1 -2t) = " << k << "is negative!!" << endl;
    exit(0);
  }

  // array holding Ps value at each data point
  double *Ps = new double[num_data];
  // array holding error rate at each data point
  double *Error_Rate = new double[num_data];
  double *SER = new double[num_data];
  // initial data point will have Probability of error == 0
  double Pss = 0.0;

  double EbN0_dB = 3.75;
  // perform a number of trials of generating
  int num_errors;
  int num_error_syms;

  // loop through data points, being sure to increment Ps by delta on each run
  // for (int i = 0; i < num_data; i++, Pss += delta)
  for (int i = 0; i < num_data; i++)
  {
    // keep users informed on progress status
    // cout << i << endl;

    EbN0_dB += 0.25;
    // perform a number of trials of generating
    num_errors = 0;
    num_error_syms = 0;

    for (int j = 0; j < num_trials_per_pt; j++)
    {
      reedSolomon rs(m, t);    // create the reed solomon object with m and t
      rs.gen_rand_msg();       // generate a random  message
      rs.encode();             // encode the given message
      rs.sim_channel(EbN0_dB); // pass encoded message through channel
      rs.decode();             // now decoded received message
      // rs.print_params();    // print some good stuff
      bool correctly_decoded = rs.compare(); // check if decoder worked
      num_error_syms += rs.comparesym();

      if (!correctly_decoded)
        num_errors++;
    }
    // printf("num_errors=%d\n",num_errors);
    // log Ps value and error rate for this data point
    Ps[i] = Pss;
    Error_Rate[i] = (double)num_errors / num_trials_per_pt;
    SER[i] = (double)num_error_syms / (num_trials_per_pt * k);

    // cout << "i=" << i << "\t\tEbN0_dB=" << EbN0_dB << "\t\tError_Rate=" << Error_Rate[i] << endl;

    cout << EbN0_dB << "\t\t" << num_error_syms << "\t\t" << Error_Rate[i] << "\t\t" << SER[i] << endl;
  }

  // remind user what parameters were used in simulation
  cout << "m: " << m << "\tt: " << t << endl;

  delete[] Error_Rate;
  delete[] SER;
  return 1;
}
