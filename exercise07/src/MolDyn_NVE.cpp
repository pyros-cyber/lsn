#include "MolDyn_NVE.h"

void MolDyn_NVE::Input(string simParameters) {
  Epot.open("results/epot.dat");
  Ekin.open("results/ekin.dat");
  Temp.open("results/temp.dat");
  Etot.open("results/etot.dat");
  Press.open("results/press.dat");
  Gave.open("results/gave.dat");
  Gerr.open("results/gerr.dat");
  Binn.open("results/binning.dat");

  ifstream ReadInput(simParameters);
  if (ReadInput.fail()) {
    cerr << "ERROR: unable to open " << simParameters << " file" << endl;
    exit(1);
  }

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl
       << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  // reading simulation & physical parameters from the input file
  ReadInput >> temp;
  ReadInput >> npart;
  // resizing of vectors containing simulations' info
  x.resize(npart);
  y.resize(npart);
  z.resize(npart);
  xold.resize(npart);
  yold.resize(npart);
  zold.resize(npart);
  vx.resize(npart);
  vy.resize(npart);
  vz.resize(npart);
  fx.resize(npart);
  fy.resize(npart);
  fz.resize(npart);

  cout << "Number of particles = " << npart << endl;
  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = static_cast<double>(npart / rho);
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol, 1. / 3.);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method "
       << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;

  ReadInput >> measure_time_interval;
  cout << "Measures performed every " << measure_time_interval << " time steps."
       << endl;
  ReadInput >> n_blocks;
  std::cout << "Statistical error computed with " << n_blocks << " blocks."
            << endl
            << endl;

  // resizing of vectors to perform the blocking averages
  est_pot.resize(n_blocks);
  est_kin.resize(n_blocks);
  est_etot.resize(n_blocks);
  est_temp.resize(n_blocks);
  est_press.resize(n_blocks);
  // initializing them to zero
  fill(est_pot.begin(), est_pot.end(), 0.);
  fill(est_kin.begin(), est_kin.end(), 0.);
  fill(est_etot.begin(), est_etot.end(), 0.);
  fill(est_temp.begin(), est_temp.end(), 0.);
  fill(est_press.begin(), est_press.end(), 0.);

  ReadInput.close();

  // compute block size (for the blocking averages) and setting block index
  block_size = (nstep / measure_time_interval) / n_blocks;
  iblock = 0;
  imeasure = 0;
  
  bin_size = (box * 0.5) / nbins;
  g_histo.resize(nbins);
  for (auto &el : g_histo) {
    el.resize(n_blocks);
    fill(el.begin(), el.end(), 0.);
  }
  
  // r_range.resize(nbins + 1);
  // getRadiusRange(r_range);
}

MolDyn_NVE::MolDyn_NVE(string simParameters, string configFile)
    : rand{"Primes", "seed.in"}, nbins{100} {

  Input(simParameters);

  // Read initial configuration from the configuration file
  // we are starting from scratch our simulation in this CONSTRUCTOR
  cout << "Read initial configuration from file " + configFile << endl << endl;

  ifstream ReadConf(configFile);
  if (ReadConf.is_open()) {
    for (int i{}; i < npart; ++i) {
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
  } else {
    cerr << "ERROR: can't read configuration from: " << configFile
         << ". Stopping simulation." << endl;
    exit(1);
  }

  ReadConf.close();

  // Prepare initial velocities
  cout << "Prepare random velocities with Maxwell-Boltzmann distribution"
       << endl;
  double sumv[3] = {};
  for (int i{}; i < npart; ++i) {
    vx[i] = rand.Rannyu() - 0.5;
    vy[i] = rand.Rannyu() - 0.5;
    vz[i] = rand.Rannyu() - 0.5;

    sumv[0] += vx[i];
    sumv[1] += vy[i];
    sumv[2] += vz[i];
  }
  for (auto &elem : sumv) {
    elem /= static_cast<double>(npart);
  }
  // fs: velocity scale factor
  double sumv_sq = 0., fs;
  for (int i{}; i < npart; ++i) {
    vx[i] = vx[i] - sumv[0];
    vy[i] = vy[i] - sumv[1];
    vz[i] = vz[i] - sumv[2];

    sumv_sq += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
  }
  sumv_sq /= static_cast<double>(npart);
  fs = sqrt(3 * double(temp) / sumv_sq);

  for (int i{}; i < npart; ++i) {
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;
    xold[i] = Pbc(x[i] - vx[i] * delta);
    yold[i] = Pbc(y[i] - vy[i] * delta);
    zold[i] = Pbc(z[i] - vz[i] * delta);
  }
}

MolDyn_NVE::MolDyn_NVE(string simParameters, string configFile,
                       string oldConfigFile)
    : rand{"Primes", "seed.in"}, nbins{100} {

  Input(simParameters);

  // Read initial configuration from the configuration file (final in this
  // case!)
  cout << "Read initial configuration from file " + configFile << endl << endl;

  ifstream ReadConf(configFile);
  if (ReadConf.is_open()) {
    for (int i{}; i < npart; ++i) {
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] *= box;
      y[i] *= box;
      z[i] *= box;
    }
  } else {
    cerr << "ERROR: can't read configuration from: " << configFile
         << ". Stopping simulation." << endl;
    exit(1);
  }
  ReadConf.close();
  ReadConf.clear();

  // read previous configuration file (from old simulation)
  cout << "Read old configuration from file " + oldConfigFile << endl << endl;
  ReadConf.open(oldConfigFile);
  if (ReadConf.is_open()) {
    for (int i{}; i < npart; ++i) {
      ReadConf >> xold[i] >> yold[i] >> zold[i];
      xold[i] *= box;
      yold[i] *= box;
      zold[i] *= box;
    }
  } else {
    cerr << "ERROR: can't read configuration from: " << oldConfigFile
         << ". Stopping simulation." << endl;
    exit(1);
  }
  ReadConf.close();
  ReadConf.clear();

  // we compute the forces acting on the particles to use them in
  // the Verlet algorithm (take a look at the loop below:)
  // x(t+dt) = 2x(t) - x(t-dt) + F(t)dt^2/m (with PBCs)
  Move();

  // ESTIMATE for TEMPERATURE
  double t = 0.;
  for (int i{}; i < npart; ++i) {
    t += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
  }
  stima_temp = t / (npart * 3.);

  double fs = sqrt(temp / stima_temp);
  cout << "scaling factor = " << fs << endl;

  // updating velocities and positions
  for (int i{}; i < npart; ++i) {
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;

    double xi = xold[i];
    double yi = yold[i];
    double zi = zold[i];

    xold[i] = Pbc(x[i] - 2. * vx[i] * delta);
    yold[i] = Pbc(y[i] - 2. * vy[i] * delta);
    zold[i] = Pbc(z[i] - 2. * vz[i] * delta);

    x[i] = xi;
    y[i] = yi;
    z[i] = zi;
  }
}
MolDyn_NVE::~MolDyn_NVE() {
  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  Press.close();
  Gave.close();
  Gerr.close();
  Binn.close();
}
void MolDyn_NVE::ConfXYZ(int nconf) const {
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

void MolDyn_NVE::ConfFinal(string filename) const {
  ofstream WriteConf(filename);

  if (WriteConf.is_open()) {
    cout << "Print configuration to " + filename << endl;
    for (int i{}; i < npart; i++) {
      WriteConf << x[i] / box << " " << y[i] / box << " " << z[i] / box << endl;
    }
  } else {
    cerr << "ERROR: can't open " << filename << endl;
    exit(1);
  }
  WriteConf.close();
}

//void MolDyn_NVE::getRadiusRange(vector<double> &vct) {
//  double bin_size = (box * 0.5) / nbins;
//  for (int i{}; i < vct.size(); ++i) {
//    vct[i] = static_cast<double>(i) * bin_size;
//  }
//}

// Move particles with Verlet algorithm
void MolDyn_NVE::Move() {
  double xnew, ynew, znew;
  Force();
  // Verlet integration scheme
  for (int i{}; i < npart; ++i) {
    xnew = Pbc(2.0 * x[i] - xold[i] + fx[i] * pow(delta, 2));
    ynew = Pbc(2.0 * y[i] - yold[i] + fy[i] * pow(delta, 2));
    znew = Pbc(2.0 * z[i] - zold[i] + fz[i] * pow(delta, 2));

    vx[i] = Pbc(xnew - xold[i]) / (2.0 * delta);
    vy[i] = Pbc(ynew - yold[i]) / (2.0 * delta);
    vz[i] = Pbc(znew - zold[i]) / (2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
}

// Compute forces as -Grad_ip V(r)
void MolDyn_NVE::Force() {

  double d1, d2, d3, dr;
  double multiplier;
  fill(fx.begin(), fx.end(), 0.);
  fill(fy.begin(), fy.end(), 0.);
  fill(fz.begin(), fz.end(), 0.);

  for (int j{}; j < npart; ++j) {
    for (int i{}; i < npart; ++i) {
      if (i != j) {
        // distance j-i in pbc
        d1 = Pbc(x[j] - x[i]);
        d2 = Pbc(y[j] - y[i]);
        d3 = Pbc(z[j] - z[i]);
        dr = sqrt(d1 * d1 + d2 * d2 + d3 * d3);

        if (dr < rcut) {
          // -Grad_ip V(r)
          multiplier = (48.0 / pow(dr, 14) - 24.0 / pow(dr, 8));
          fx[j] += d1 * multiplier;
          fy[j] += d2 * multiplier;
          fz[j] += d3 * multiplier;
        }
      }
    }
  }
}

// Properties measurement
void MolDyn_NVE::Measure() {
  double v = 0., t = 0.;
  double dx, dy, dz, dr;
  unsigned int cf = 0;
  stima_press = 0.;
  // cycle over pairs of particles
  for (int i{}; i < npart - 1; ++i) {
    for (int j{i + 1}; j < npart; ++j) {
      // here I use old configurations [old = r(t)]
      // to be compatible with EKin which uses v(t)
      // => EPot should be computed with r(t)
      dx = Pbc(xold[i] - xold[j]);
      dy = Pbc(yold[i] - yold[j]);
      dz = Pbc(zold[i] - zold[j]);
      dr = sqrt(dx * dx + dy * dy + dz * dz);

      //if (dr < box * 0.5) {
      //  cf = static_cast<unsigned int>((dr * nbins) / (box * 0.5));
      //  //cout << "cf = " << cf << " iblock = " << iblock << endl;
      //  g_histo[cf][iblock] += 2;
      //}
      for (int k{}; k < nbins; ++k){
        if(dr > bin_size*k && dr < bin_size*(k+1)){
          g_histo[k][iblock] += 3. / (pow(bin_size+dr,3)-pow(dr,3)) / (2.*M_PI);
          break;
        }
      }
      if (dr < rcut) {
        // Potential energy
        v += 4.0 / pow(dr, 12) - 4.0 / pow(dr, 6);
        stima_press += 1. / pow(dr, 12) - 0.5 / pow(dr, 6);
      }
    }
  }

  // Kinetic energy
  for (int i{}; i < npart; ++i) {
    t += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
  }
  // Potential energy per particle
  stima_pot = static_cast<double>(v / npart);
  // Kinetic energy per particle
  stima_kin = static_cast<double>(t / npart);
  // Temperature
  stima_temp = static_cast<double>((2.0 / 3.0) * t / npart);
  // Total energy per particle
  stima_etot = static_cast<double>((t + v) / npart);
  stima_press = 16. * stima_press / (vol);
  stima_press += stima_temp * rho;

  if (Epot.is_open() && Ekin.is_open() && Temp.is_open() && Etot.is_open() &&
      Press.is_open()) {
    Epot << stima_pot << endl;
    Ekin << stima_kin << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Press << stima_press << endl;
  } else {
    cerr << "ERROR: unable to open output files." << endl;
    exit(1);
  }

  imeasure++;
  iblock = imeasure / block_size;
  est_pot[iblock] += stima_pot;
  est_kin[iblock] += stima_kin;
  est_temp[iblock] += stima_temp;
  est_etot[iblock] += stima_etot;
  est_press[iblock] += stima_press;
}

void MolDyn_NVE::BlockingResults() {
  for (auto &elem : est_pot) {
    elem /= block_size;
  }
  for (auto &elem : est_kin) {
    elem /= block_size;
  }
  for (auto &elem : est_temp) {
    elem /= block_size;
  }
  for (auto &elem : est_etot) {
    elem /= block_size;
  }
  for (auto &elem : est_press) {
    elem /= block_size;
  }
  for(auto &elem : g_histo) {
    for(auto &k : elem) {
      k /= block_size;
    }
  }
  //for (int i{}; i < g_histo.size(); ++i) {
  //  for (int j{}; j < g_histo[i].size(); ++j) {
  //    //cout << "r_range[" << j + 1 << "] = " << r_range[j + 1] << " r_range[" << j << "] = " << r_range[j] << endl;
  //    double norm_g = ((4. * M_PI) / 3) * npart * rho *
  //                    (pow(r_range[j + 1], 3) - pow(r_range[j], 3));
  //    //cout << "g_histo[" << i << "][" << j << "] = " << g_histo[i][j] << " norm_g = " << norm_g << " block_size = " << block_size << endl;
  //    g_histo[i][j] = g_histo[i][j] / block_size / norm_g;
  //    //cout << "g_histo[" << i << "][" << j << "] = " << g_histo[i][j] << endl;
  //  }
  //}

  vector<double> pot_err = blocking_error(est_pot);
  vector<double> kin_err = blocking_error(est_kin);
  vector<double> temp_err = blocking_error(est_temp);
  vector<double> etot_err = blocking_error(est_etot);
  vector<double> press_err = blocking_error(est_press);
  vector<vector<double>> g_histo_err(g_histo.size());
  for (int i{}; i < g_histo.size(); ++i) {
    g_histo_err[i] = blocking_error(g_histo[i]);
  }
  // for (int i{}; i < g_histo.size(); ++i) {
  //   for (int j{}; j < g_histo[i].size(); ++j) {
  //     cout << "g_histo[" << i << "][" << j << "] = " << g_histo[i][j] << " g_histo_err[" << i << "][" << j << "] = " << g_histo_err[i][j] << endl;
  //   }
  // }
  vector<double> g_norm(g_histo[0].size(), 0.);
  for(int j{}; j < g_norm.size(); ++j) {
    for(int i{}; i < g_histo.size(); ++i){
      g_norm[j] += g_histo[i][j];
    }
  }

  if (Gave.is_open() && Gerr.is_open() && Binn.is_open()){
    for(int i{}; i < g_histo[0].size(); ++i) {
        for(int j{}; j < g_histo.size(); ++j) {
          Gave << g_histo[j][i] / g_norm[i] << " ";
          Gerr << g_histo_err[j][i] / g_norm[i] << " ";
          Binn << j * bin_size + bin_size * 0.5 << " ";
        }
        Binn << endl;
        Gave << endl;
        Gerr << endl;
    }
  } else {
    cerr << "ERROR: unable to open output file." << endl;
    exit(1);
  }

  ofstream out("results/ave_epot.dat");
  if (out.is_open()) {
    for (int i{}; i < est_pot.size(); ++i)
      out << i << " " << est_pot[i] << " " << pot_err[i] << endl;
  } else {
    cerr << "ERROR: unable to open output file." << endl;
    exit(1);
  }
  out.close();

  out.open("results/ave_ekin.dat");
  if (out.is_open()) {
    for (int i{}; i < est_kin.size(); ++i)
      out << i << " " << est_kin[i] << " " << kin_err[i] << endl;
  } else {
    cerr << "ERROR: unable to open output file." << endl;
    exit(1);
  }
  out.close();

  out.open("results/ave_etot.dat");
  if (out.is_open()) {
    for (int i{}; i < est_etot.size(); ++i)
      out << i << " " << est_etot[i] << " " << etot_err[i] << endl;
  } else {
    cerr << "ERROR: unable to open output file." << endl;
    exit(1);
  }
  out.close();

  out.open("results/ave_temp.dat");
  if (out.is_open()) {
    for (int i{}; i < est_temp.size(); ++i)
      out << i << " " << est_temp[i] << " " << temp_err[i] << endl;
  } else {
    cerr << "ERROR: unable to open output file." << endl;
    exit(1);
  }
  out.close();

  out.open("results/ave_press.dat");
  if (out.is_open()) {
    for (int i{}; i < est_press.size(); ++i)
      out << i << " " << est_press[i] << " " << press_err[i] << endl;
  } else {
    cerr << "ERROR: unable to open output file." << endl;
    exit(1);
  }

  out.close();
}

void MolDyn_NVE::RunSimulation() {
  unsigned int nconf = 1;              // starting from configuration number 1
  ofstream out("frames/config_1.xyz"); // check if frames folder exists
  if (out.fail()) {
    cout << "\n\nERROR: Unable to write on folder 'frames', create it"
            "and run again!\n\n";
    exit(2);
  }
  out.close();
  // modified in order to print one xyz file less: the corresponding
  // configuration will be printed in the last lines to "old.0" and
  // "config.final".
  for (int istep{}; istep < nstep; ++istep) {
    Move();
    if (istep % iprint == 0) {
      // cout << "Number of time-steps: " << istep << endl;
    }
    if (istep % measure_time_interval == 0) {
      Measure();
      ConfXYZ(nconf);
      nconf++;
    }
  }
  cout << "Number of time-steps: " << nstep << endl << endl;
  ConfFinal("old.0");
  Move();
  Measure();
  ConfFinal("config.final");
  BlockingResults();
}
