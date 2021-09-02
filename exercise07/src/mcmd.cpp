#include "mcmd.h"

mcmd::mcmd(string simParameters, string initial_configuration,
           bool _instantenous)
    : rnd{"Primes", "seed.in"}, nbins{100} {
  instantenous = _instantenous;
  props = {"energy", "pressure"};
  for (auto &el : props) {
    walker[el] = 0.;
    block_average[el] = 0.;
    glob_average[el] = 0.;
    glob_average_sq[el] = 0.;
  }
  g_histo.resize(nbins);
  g_histo_block_ave.resize(nbins);
  g_histo_glob_ave.resize(nbins);
  g_histo_glob_ave_sq.resize(nbins);

  accepted = 0;
  attempted = 0;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl
       << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T "
       << endl
       << endl;
  cout << "The program uses Lennard-Jones units \n" << endl;
  cout << "Reading parameters from " << simParameters << endl;

  /* READING INPUT FILES AND INITIALIZING PARAMETERS */
  ifstream ReadInput(simParameters);
  if (ReadInput.fail()) {
    cerr << "ERROR: unable to open " << simParameters << " file.\nExiting.\n";
    exit(1);
  }
  ReadInput >> temp;
  beta = 1.0 / temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = static_cast<double>(npart) / rho;
  box = pow(vol, 1. / 3.);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;

  // Tail corrections for potential energy and pressure
  vtail = (8.0 * M_PI * rho) / (9.0 * pow(rcut, 9)) -
          (8.0 * M_PI * rho) / (3.0 * pow(rcut, 3));
  ptail = (32.0 * M_PI * rho) / (9.0 * pow(rcut, 9)) -
          (16.0 * M_PI * rho) / (3.0 * pow(rcut, 3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Tail correction for the pressure         = " << ptail << endl;

  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> nblk;

  cout << "The program perform Metropolis moves with uniform translations"
       << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

  r_range.resize(nbins + 1);
  getRadiusRange(r_range);

  x.resize(npart);
  y.resize(npart);
  z.resize(npart);

  ifstream ReadConf(initial_configuration);
  if (ReadConf.fail()) {
    cerr << "ERROR: unable to open " << initial_configuration
         << " file.\nExiting.\n";
    exit(1);
  }
  for (int i{}; i < npart; ++i) {
    double x_temp, y_temp, z_temp;
    ReadConf >> x_temp >> y_temp >> z_temp;
    x[i] = Pbc(x_temp * box);
    y[i] = Pbc(y_temp * box);
    z[i] = Pbc(z_temp * box);
  }
  ReadConf.close();

  /* OPENING OUTPUT STREAMS */
  Epot.open("results/epot.dat");
  Pres.open("results/pres.dat");
  Gerr.open("results/gerr.dat");
  Gave.open("results/gave.dat");
  Binn.open("results/binning.dat");

  if (instantenous) {
    ist_pot.open("results/instant_epot.dat");
    ist_pres.open("results/instant_pres.dat");
  }

  // Evaluate potential energy and pressure of the initial configuration
  Measure();
  // Print initial values for the potential energy and pressure
  cout << "Initial potential energy (with tail corrections) = "
       << walker.at("energy") / npart + vtail << endl;
  cout << "Virial                   (with tail corrections) = "
       << walker.at("pressure") / npart + ptail << endl;
  cout << "Pressure                 (with tail corrections) = "
       << rho * temp + (walker.at("pressure") + npart * ptail) / vol << endl
       << endl;
}

mcmd::~mcmd() {
  Epot.close();
  Pres.close();
  Gerr.close();
  Gave.close();
  if (instantenous) {
    ist_pot.close();
    ist_pres.close();
  }
}

void mcmd::getRadiusRange(vector<double> &vct) {
  double bin_size = (box * 0.5) / nbins;
  for (int i{}; i < vct.size(); ++i) {
    vct[i] = static_cast<double>(i) * bin_size;
  }
}

double mcmd::Boltzmann(double _x, double _y, double _z, unsigned int ip) {
  double ene = 0., dx = 0., dy = 0., dz = 0., dr = 0.;

  for (int i{}; i < npart; ++i) {
    if (i != ip) {
      // distance ip-i in pbc
      dx = Pbc(_x - x[i]);
      dy = Pbc(_y - y[i]);
      dz = Pbc(_z - z[i]);

      dr = sqrt(dx * dx + dy * dy + dz * dz);

      if (dr < rcut) {
        ene += 1. / pow(dr, 12) - 1. / pow(dr, 6);
      }
    }
  }
  return 4. * ene;
}

void mcmd::Reset(unsigned int iblk) {
  // Reset block averages
  // we choose to lookup the maps' elements with .at for safety:
  // If you access a key using the indexing operator [] that is
  // not currently a part of a map, then it automatically adds
  // a key for you!!
  for (auto &el : props) {
    block_average.at(el) = 0.;
  }
  fill(g_histo_block_ave.begin(), g_histo_block_ave.end(), 0.);
  attempted = 0.;
  accepted = 0.;
}

void mcmd::Move() {
  int o = 0;
  double p = 0., x_old = 0., y_old = 0., z_old = 0., x_new = 0., y_new = 0.,
         z_new = 0.;

  for (int i{}; i < npart; ++i) {
    o = static_cast<int>(rnd.Rannyu() * npart);

    x_old = x[o];
    y_old = y[o];
    z_old = z[o];

    x_new = Pbc(x[o] + delta * (rnd.Rannyu() - 0.5));
    y_new = Pbc(y[o] + delta * (rnd.Rannyu() - 0.5));
    z_new = Pbc(z[o] + delta * (rnd.Rannyu() - 0.5));

    p = exp(beta * (Boltzmann(x_old, y_old, z_old, o) -
                    Boltzmann(x_new, y_new, z_new, o)));
    if (p >= rnd.Rannyu()) {
      x[o] = x_new;
      y[o] = y_new;
      z[o] = z_new;

      accepted++;
    }
    attempted++;
  }
}
void mcmd::Measure(unsigned int istep) {
  double v = 0., w = 0., dx = 0., dy = 0., dz = 0., dr = 0.;
  unsigned int cf = 0;

  fill(g_histo.begin(), g_histo.end(), 0.);

  for (int i{}; i < npart - 1; ++i) {
    for (int j{i + 1}; j < npart; ++j) {
      // distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = sqrt(dx * dx + dy * dy + dz * dz);

      if (dr < box * 0.5) {
        cf = static_cast<unsigned int>((dr * nbins) / (box * 0.5));
        g_histo[cf] += 2;
      }

      if (dr < rcut) {
        v += 1. / pow(dr, 12) - 1. / pow(dr, 6);
        w += 1. / pow(dr, 12) - 0.5 / pow(dr, 6);
      }
    }
  }
  walker.at("energy") = 4. * v;
  walker.at("pressure") = 48. * w / 3.;

  if (instantenous) {
    double e = walker.at("energy") / static_cast<double>(npart) + vtail;
    double p =
        rho * temp +
        (walker.at("pressure") + ptail + static_cast<double>(npart)) / vol;
    ist_pot << istep << " " << e << endl;
    ist_pres << istep << " " << p << endl;
  }
}

void mcmd::Accumulate() {
  for (auto &el : props) {
    block_average.at(el) += walker.at(el);
  }
  for (int i{}; i < nbins; ++i) {
    g_histo_block_ave[i] += g_histo[i];
  }
}

void mcmd::Averages(unsigned int iblk) {

  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted / attempted << endl << endl;

  // Potential energy
  double stima_pot, err_pot;
  stima_pot =
      (block_average.at("energy") / nstep) / static_cast<double>(npart) + vtail;
  glob_average.at("energy") += stima_pot;
  glob_average_sq.at("energy") += stima_pot * stima_pot;
  err_pot =
      Error(glob_average.at("energy"), glob_average_sq.at("energy"), iblk);
  Epot << " " << iblk << " " << stima_pot << " "
       << glob_average.at("energy") / static_cast<double>(iblk) << " "
       << err_pot << endl;

  // Pressure
  double stima_pres, err_pres;
  stima_pres = rho * temp + ((block_average.at("pressure") / nstep) +
                             ptail * static_cast<double>(npart)) /
                                vol;
  glob_average.at("pressure") += stima_pres;
  glob_average_sq.at("pressure") += stima_pres * stima_pres;
  err_pres =
      Error(glob_average.at("pressure"), glob_average_sq.at("pressure"), iblk);
  Pres << " " << iblk << " " << stima_pres << " "
       << glob_average.at("pressure") / static_cast<double>(iblk) << " "
       << err_pres << endl;

  // Radial distribution
  double stima_g, norm_g, err_g;
  for (int i{}; i < nbins; ++i) {
    double bin_size = (box * 0.5) / nbins;
    norm_g = ((4. * M_PI) / 3) * npart * rho *
             (pow(r_range[i + 1], 3) - pow(r_range[i], 3));
    stima_g = (g_histo_block_ave[i] / nstep) / norm_g;
    g_histo_glob_ave[i] += stima_g;
    g_histo_glob_ave_sq[i] += stima_g * stima_g;
    err_g = Error(g_histo_glob_ave[i], g_histo_glob_ave_sq[i], iblk);
    Gave << g_histo_glob_ave[i] << " ";
    Gerr << err_g << " ";
    Binn << i * bin_size + bin_size / 2. << " ";
  }

  Gave << endl;
  Gerr << endl;
  Binn << endl;
}

void mcmd::ConfXYZ(unsigned int nconf) {
  ofstream WriteXYZ("frames/config_" + to_string(nconf) + ".xyz");

  if (WriteXYZ.is_open()) {
    WriteXYZ << npart << endl;
    WriteXYZ << "This is only a comment!" << endl;
    for (int i{}; i < npart; ++i) {
      WriteXYZ << "LJ  " << Pbc(x[i]) << "   " << Pbc(y[i]) << "   "
               << Pbc(z[i]) << endl;
    }
  } else {
    cerr << "ERROR: can't open output XYZ file." << endl;
    exit(1);
  }
  WriteXYZ.close();
}

void mcmd::ConfFinal() const {
  ofstream WriteConf("config.final");

  if (WriteConf.is_open()) {
    cout << "Print configuration to config.final" << endl;
    for (int i{}; i < npart; i++) {
      WriteConf << x[i] / box << " " << y[i] / box << " " << z[i] / box << endl;
    }
  } else {
    cerr << "ERROR: can't open conf.final" << endl;
    exit(1);
  }
  WriteConf.close();
}

void mcmd::Run() {
  unsigned int nconf = 1;
  // Simulation
  for (int iblk{1}; iblk <= nblk; ++iblk) {
    Reset(iblk); // Reset block averages
    for (int istep{1}; istep <= nstep; ++istep) {
      Move();
      Measure(istep);
      Accumulate(); // Update block averages

      if (istep % 10 == 0) {
        ConfXYZ(nconf);
        nconf++;
      }
    }
    Averages(iblk); // Print results for current block
  }
  ConfFinal(); // Write final configuration
  rnd.SaveSeed();
}
