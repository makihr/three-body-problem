#include <cmath>
#include <iostream>
#include <ps_types.hpp>
#include <iomanip>
#include <particle_simulator.hpp>
#include "kepler.hpp"

using namespace std;

struct Particle {
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64 mass;
};
struct OrbitalParam {
    PS::F64 ax;
    PS::F64 ecc;
    PS::F64 inc;
    PS::F64 OMG;
    PS::F64 omg;
    PS::F64 u;
};
/*
const double ax = 1.0;
const double ecc = 0.5;
const double inc = 0.1; //
const double OMG = 0.0; // 0-2pi
const double omg = 0.0;
const double u = 0.0;

const m0 = 1.0;
const m1 = 1.0;
PS::F64vec pos0, pos1, vel0, vel1;
void OrbParam2PosVel(pos0, pos1, vel0, vel1, m0, m1,
        ax, ecc, inc, OMG, omg, u);
template<typename Tptcl>
void make_hierarchical_system(const int nbody, Tptcl ptcl[]{



        }
*/

/*void integrate_euler(int i, double dt, double r[][3], double v[][3], double
a[][3])
{
        for (int k = 0; k < 3; k++)
        {
                r[i][k] += v[i][k] * dt;
                v[i][k] += a[i][k] * dt;
        }
}
*/
void acceleration(int n, double r[][3], double a[][3]) {
    const double m = 1;
    for (int i = 0; i < n; i++)
        for (int k = 0; k < 3; k++)
            a[i][k] = 0;

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double rji[3];
            for (int k = 0; k < 3; k++)
                rji[k] = r[j][k] - r[i][k];
            double r2 = 0;
            for (int k = 0; k < 3; k++) {
                r2 += rji[k] * rji[k];
            }

            double r3 = r2 * sqrt(r2);
            for (int k = 0; k < 3; k++) {
                a[i][k] += m * rji[k] / r3;
                a[j][k] -= m * rji[k] / r3;
            }
        }
    }
}

double t1[3], t2[3], t3[3], t4[3], t5[3], t6[3], t7[3];
double r12[3], v12[3],n12[3];

void accelerationPNS(int n, double r[][3], double v[][3], double a[][3]) {
    const double m = 1;
    double c = 100.0;

    for (int i = 0; i < n; i++)
        for (int k = 0; k < 3; k++)
            a[i][k] = 0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            double rji[3];
            for (int k = 0; k < 3; k++)
                rji[k] = r[j][k] - r[i][k];
            double r2 = 0;
            for (int k = 0; k < 3; k++) {
                r2 += rji[k] * rji[k];
            }
            double r = sqrt(r2);
            double r3 = r * r2;
            for (int k = 0; k < 3; k++) {
                n12[k] = -(rji[k]) / r;
                v12[k] = v[i][k] - v[j][k];
            }
            double n12v2_2 = 0.0;
            double v1_2 = 0.0;
            double v1v2 = 0.0;
            double v2_2 = 0.0;
            double n12v1 = 0.0;
            double n12v2 = 0.0;
            double r12r12 = 0.0;

            for (int k = 0; k < 3; k++) {
                n12v2_2 += n12[k] * v[j][k];
                v1_2 += v[i][k] * v[i][k];
                v1v2 += v[i][k] * v[j][k];
                v2_2 += v[j][k] * v[j][k];
                n12v1 = n12[k] * v[i][k];
                n12v2 = n12[k] * v[j][k];
            }

            n12v2_2 *= n12v2_2;
            for (int k = 0; k < 3; k++) {
                a[i][k] += ( - m * n12[k] / r2) + (1.0 / (c * c) * (((5 * m * m / r3 + 4 * m * m / r3 + m / r3 * (1.5 * n12v2_2 - v1_2 + 4.0 * v1v2 - 2.0 * v2_2)) * n12[k]) + m / r2 * (4 * n12v1 - 3 * n12v2) * v12[k]));
                a[j][k] -= ( - m * n12[k] / r2) + (1.0 / (c * c) * (((5 * m * m / r3 + 4 * m * m / r3 + m / r3 * (1.5 * n12v2_2 - v1_2 + 4.0 * v1v2 - 2.0 * v2_2)) * n12[k]) + m / r2 * (4 * n12v1 - 3 * n12v2) * v12[k]));
            }
        }
    }
}


void runge(int n, Particle ptcl[], double dt, double r[][3], double v[][3])
{
    double k1_r[n][3], k1_v[n][3];
    double k2_r[n][3], k2_v[n][3];
    double k3_r[n][3], k3_v[n][3];
    double k4_r[n][3], k4_v[n][3];
    double a[n][3];
    double r0[n][3], v0[n][3];
    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            r0[i][k] = r[i][k];
            v0[i][k] = v[i][k];
        }
    }

    //acceleration(n, r, a);
    accelerationPNS(n, r, v, a);
    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            k1_r[i][k] = v0[i][k];
            k1_v[i][k] = a[i][k];
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            r[i][k] = r0[i][k] + 0.5 * dt * k1_r[i][k];
            v[i][k] = v0[i][k] + 0.5 * dt * k1_v[i][k];
        }
    }
    //acceleration(n, r, a);
    accelerationPNS(n, r, v, a);
    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            k2_r[i][k] = v[i][k];
            k2_v[i][k] = a[i][k];
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            r[i][k] = r0[i][k] + 0.5 * dt * k2_r[i][k];
            v[i][k] = v0[i][k] + 0.5 * dt * k2_v[i][k];
        }
    }

    //acceleration(n, r, a);
    accelerationPNS(n, r, v, a);
    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            k3_r[i][k] = v[i][k];
            k3_v[i][k] = a[i][k];
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            r[i][k] = r0[i][k] + dt * k3_r[i][k];
            v[i][k] = v0[i][k] + dt * k3_v[i][k];
        }
    }

    //acceleration(n, r, a);
    accelerationPNS(n, r, v, a);
    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            k4_r[i][k] = v[i][k];
            k4_v[i][k] = a[i][k];
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            r[i][k] = r0[i][k] + (dt / 6) * (k1_r[i][k] + 2 * k2_r[i][k] + 2 * k3_r[i][k] + k4_r[i][k]);
            v[i][k] = v0[i][k] + (dt / 6) * (k1_v[i][k] + 2 * k2_v[i][k] + 2 * k3_v[i][k] + k4_v[i][k]);
        }
    }
}

void energy_calculation(int n, double r[][3], double v[][3], double& ekin, double& epot) {
    const double m = 1;
    ekin = epot = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double rji[3];
            for (int k = 0; k < 3; k++)
                rji[k] = r[j][k] - r[i][k];
            double r2 = 0;
            for (int k = 0; k < 3; k++)
                r2 += rji[k] * rji[k];
            double r = sqrt(r2);
            epot -= m * m / r;
        }
        for (int k = 0; k < 3; k++)
            ekin += 0.5 * m * v[i][k] * v[i][k];
    }
}


int main() {
    int n = 3;

    Particle ptcl[n];
    for (int i = 0; i < n; i++) {
        ptcl[i].mass = 1.0;
    }



    OrbitalParam orb[n - 1];
    orb[0].ax = 1.0;
    orb[0].ecc = 0.1;
    orb[0].inc = 0.0;
    orb[0].OMG = 0.0;
    orb[0].omg = 0.0;
    orb[0].u = 0.3;
    orb[1].ax = 3.0; // 7.5, 7.0, 6.5, ,,, 3.0
    orb[1].ecc = 0.5;
    orb[1].inc = 0.0;
    orb[1].OMG = 0.0;
    orb[1].omg = 0.0;
    orb[1].u = 1.0;
    if (orb[0].ecc == 0.0)
        orb[0].omg = 0.0;
    if (orb[1].ecc == 0.0)
        orb[1].omg = 0.0;
    if (orb[0].inc == 0.0)
        orb[0].OMG = 0.0;
    if (orb[1].inc == 0.0)
        orb[1].OMG = 0.0;
    MakeHierarchicalSystem(ptcl, orb, n);

    // �R�̖ڂ̎���
    double Tkep_out = 2.0 * M_PI * sqrt(orb[1].ax * orb[1].ax * orb[1].ax / 3.0);
    

    double r[n][3], v[n][3], a[n][3];
    for (int i = 0; i < n; i++) {
        r[i][0] = ptcl[i].pos.x;
        r[i][1] = ptcl[i].pos.y;
        r[i][2] = ptcl[i].pos.z;
        v[i][0] = ptcl[i].vel.x;
        v[i][1] = ptcl[i].vel.y;
        v[i][2] = ptcl[i].vel.z;
    }


    const double m = 1;
    double dt, t_end;

    cerr << "Please provide a value for the time step" << endl;
    cin >> dt;
    //cerr << "and for the duration of the run" << endl;
    //cin >> t_end;
    t_end = 10000 * Tkep_out;


    double dt_out = 0.01;
    double t_out = dt_out;
    double ekin = 0, epot = 0;

    //acceleration(n, r, a);
    accelerationPNS(n, r, v, a);

    energy_calculation(n, r, v, ekin, epot);
    double e_in = ekin + epot;
    cerr << "Initial totalenergyE_in=" << e_in << endl;
    dt_out = dt;
    t_out = dt_out;

    for (double t = 0; t < t_end; t += dt) {
        runge(n, ptcl, dt, r, v);
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < 3; k++)
                cout << r[i][k] << " ";
            for (int k = 0; k < 3; k++)
                cout << v[i][k] << " ";
        }
        cout << endl;
        t_out += dt_out;
    }

    energy_calculation(n, r, v, ekin, epot);
    double e_out = ekin + epot;

    cerr << "Final total energy E_out = " << e_out << endl;
    cerr << "absolute energy error: E_out- E_in = " << e_out - e_in << endl;
    cerr << "relative energy error: (E_out- E_in) / E_in = "
         << (e_out - e_in) / e_in << endl;
}