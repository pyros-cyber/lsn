#include "args.hxx"
#include "myStatFunc.h"
#include "varMC.h"

int main(int argc, char *argv[]) {
  args::ArgumentParser parser("This program computes the GS energy of a 1D "
                              "particle using a Metropolis algorithm");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::Group arguments(parser, "arguments", args::Group::Validators::DontCare,
                        args::Options::Global);
  args::ValueFlag<double> mu(arguments, "", "Variational parameter",
                             {'m', "mu"});
  args::ValueFlag<double> sigma(arguments, "", "Variational parameter",
                                {'s', "sigma"});
  args::ValueFlag<int> M(arguments, "", "Number of points to sample", {'M'});
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    cout << parser;
    return 0;
  } catch (args::ParseError e) {
    cerr << e.what() << std::endl;
    cerr << parser;
    return 1;
  } catch (args::ValidationError e) {
    cerr << e.what() << std::endl;
    cerr << parser;
    return 1;
  }

  // cout << "mu = " << args::get(mu) << " sigma = " << args::get(sigma)
  //     << " M = " << args::get(M) << endl;
  Particle p(args::get(mu), args::get(sigma), 0., 2.7);
  double _x, sum = 0., sum_sq = 0., err = 0.;
  int _M, N = 100;
  if (args::get(M) == 0) {
    _M = 100000;
  } else {
    _M = args::get(M);
  }
  int L = _M / N;
  vector<double> psi_sq;
  ofstream out_en("energy.dat");
  ofstream out_psi("psi_sq.dat");

  if (out_en.is_open() && out_psi.is_open()) {
    for (int i{}; i < N; ++i) {
      double en = 0., en_sq = 0.;
      for (int j{}; j < L; ++j) {
        p.MetropolisMove();
        _x = p.getX();
        en += p.LocEnergy(_x);
        psi_sq.push_back(_x);
      }
      en /= static_cast<double>(L);
      en_sq = pow(en, 2);
      sum = (sum * i + en) / (static_cast<double>(i + 1.));
      sum_sq = (sum_sq * i + en_sq) / (static_cast<double>(i + 1.));
      err = error(sum, sum_sq, i);
      out_en << i + 1 << " " << sum << " " << err << endl;
    }
    cout << "Acceptance rate: " << p.Acceptance() << endl;

    for (int i{}; i < L; ++i) {
      out_psi << psi_sq[i] << endl;
    }
  } else {
    cerr << "ERROR: can't open output files!" << endl;
  }
  return 0;
}
