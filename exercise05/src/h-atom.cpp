#include "args.hxx"
#include "hydrogen.h"
#include "metropolis.h"
#include "myStatFunc.h"
#include "random.h"

using namespace std;

int main(int argc, char *argv[]) {

  args::ArgumentParser parser(
      "This program estimates the value of <r> for the ground state and first "
      "excited state of the H2 atom");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::Group commands(parser, "commands");
  args::Command origin(commands, "origin",
                       "Start the simulation in the origin");
  args::Command far(commands, "far",
                    "Start the simulation far from the origin");
  args::Group arguments(parser, "arguments", args::Group::Validators::DontCare,
                        args::Options::Global);
  args::ValueFlag<string> state(
      arguments, "",
      "State of H2 atom, can be either 'ground' or 'first-excited'",
      {'s', "state"});
  args::ValueFlag<string> transition(
      arguments, "",
      "Transition probability to use in metropolis algorithm. Can either be "
      "'uniform' or 'gaussian'",
      {'t', "transition"});
  args::ValueFlag<double> step(arguments, "",
                               "Step of the Metropolis algorithm", {"step"});
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

  unsigned int M = 1e7, N = 100;
  int throws_per_block = M / N;
  if (origin) {
    HydrogenAtom hydro(0, 0, 0, args::get(state));
    Metropolis metro(args::get(step), hydro.x, hydro.y, hydro.z);
    vector<double> r(N, 0.);
    ofstream out_position("positions." + args::get(transition) + "." +
                          hydro.state + ".dat");
    if (out_position.is_open()) {
      for (int i{}; i < r.size(); ++i) {
        for (int j{}; j < throws_per_block; ++j) {
          out_position << metro.x_start << " " << metro.y_start << " "
                       << metro.z_start << endl;
          if (args::get(transition) == "uniform" && hydro.state == "ground") {
            metro.Run(
                [&metro](double coordinate) {
                  return metro.uniform_step(coordinate);
                },
                ground_state);
          } else if (args::get(transition) == "uniform" &&
                     hydro.state == "first-excited") {
            metro.Run(
                [&metro](double coordinate) {
                  return metro.uniform_step(coordinate);
                },
                first_excited_state);
          } else if (args::get(transition) == "gaussian" &&
                     hydro.state == "ground") {
            metro.Run(
                [&metro](double coordinate) {
                  return metro.gaussian_step(coordinate);
                },
                ground_state);
          } else if (args::get(transition) == "gaussian" &&
                     hydro.state == "first-excited") {
            metro.Run(
                [&metro](double coordinate) {
                  return metro.gaussian_step(coordinate);
                },
                first_excited_state);
          } else {
            cerr << "Unknown transition probability or hydrogen "
                    "state.\nExiting.\n";
            return 1;
          }

          r[i] += sqrt(metro.x_start * metro.x_start +
                       metro.y_start * metro.y_start +
                       metro.z_start * metro.z_start);
        }
        r[i] /= throws_per_block;
      }

      cout << "Metropolis algorithm acceptance: " << metro.acceptance() << endl;

      vector<double> radius_error = blocking_error(r);

      ofstream out("results." + args::get(transition) + "." + hydro.state +
                   ".origin.dat");
      if (out.is_open()) {
        for (int i{}; i < r.size(); ++i) {
          out << r[i] << " " << radius_error[i] << endl;
        }
      } else {
        cerr << "ERROR: unable to open output file" << endl;
        return 1;
      }
      out.close();
    } else {
      cerr << "ERROR: unable to open output file" << endl;
      return 1;
    }
  }
  if (far) {
    M = 10000;
    HydrogenAtom hydro(100, 100, 100, args::get(state));
    Metropolis metro(args::get(step), hydro.x, hydro.y, hydro.z);
    vector<double> r(M, 0.);
    for (int i{}; i < r.size(); ++i) {
      if (args::get(transition) == "uniform" && hydro.state == "ground") {
        metro.Run(
            [&metro](double coordinate) {
              return metro.uniform_step(coordinate);
            },
            ground_state);
      } else if (args::get(transition) == "uniform" &&
                 hydro.state == "first-excited") {
        metro.Run(
            [&metro](double coordinate) {
              return metro.uniform_step(coordinate);
            },
            first_excited_state);
      } else if (args::get(transition) == "gaussian" &&
                 hydro.state == "ground") {
        metro.Run(
            [&metro](double coordinate) {
              return metro.gaussian_step(coordinate);
            },
            ground_state);
      } else if (args::get(transition) == "gaussian" &&
                 hydro.state == "first-excited") {
        metro.Run(
            [&metro](double coordinate) {
              return metro.gaussian_step(coordinate);
            },
            first_excited_state);
      } else {
        cerr << "Unknown transition probability or hydrogen "
                "state.\nExiting.\n";
        return 1;
      }
      r[i] =
          sqrt(metro.x_start * metro.x_start + metro.y_start * metro.y_start +
               metro.z_start * metro.z_start);
    }
    cout << "Metropolis algorithm acceptance: " << metro.acceptance() << endl;

    ofstream out("results." + args::get(transition) + "." + hydro.state +
                 ".far.dat");
    if (out.is_open()) {
      for (int i{}; i < r.size(); ++i) {
        out << i << " " << r[i] << endl;
      }
    } else {
      cerr << "ERROR: unable to open output file" << endl;
      return 1;
    }
    out.close();
  }

  return 0;
}