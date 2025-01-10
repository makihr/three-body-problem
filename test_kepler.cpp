#include <iostream>
#include <particle_simulator.hpp>
#include "kepler.hpp"

struct Ptcl {
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64 mass;
    Ptcl() : pos(0.0), vel(0.0), mass(0.0) {}
    Ptcl(const PS::F64vec& pos_, const PS::F64vec& vel_, const PS::F64 mass_) : pos(pos_), vel(vel_), mass(mass_) {}
};
struct OrbitalParam {
    PS::F64 ax;
    PS::F64 ecc;
    PS::F64 inc;
    PS::F64 OMG;
    PS::F64 omg;
    PS::F64 u;
};

int main() {
    constexpr double PI = 4.0 * atan(1.0);
    const PS::U32 seed = 0;
    std::mt19937 engine(seed);
    std::uniform_real_distribution<double> dist_0_p(0.0, PI);
    std::uniform_real_distribution<double> dist_0_2p(0.0, 2.0 * PI);
    std::uniform_real_distribution<double> dist_mp_m(-PI, PI);

    // Hierarchical system
    // the following settings are similar to Myllari et al. 2018
    constexpr int N = 3;
    Ptcl ptcl[N];
    ptcl[0].mass = 1.5;
    ptcl[1].mass = 0.5;
    ptcl[2].mass = 0.5;

    double ax_in_ar[] = {1.0};
    double e_out_ar[] = {0.5};
    double e_in_ar[] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9, 0.95, 0.99};
    const int n_ci = 100;  // cos(i)
    // const int n_ci = 1;
    double cosinc_ar[n_ci];
    cosinc_ar[0] = 1.0 - 1e-10;
    for (int i = 1; i < n_ci; i++) {
        cosinc_ar[i] = 1.0 - 2.0 * i / (n_ci - 1);
    }
    const int n_ax_out = 37;
    double ax_out_ar[n_ax_out];  // 10.0-2.0
    for (int i = 0; i < n_ax_out; i++) {
        ax_out_ar[i] = 10.0 - 0.25 * i;
    }
    for (auto e_in : e_in_ar) {
        for (auto ci : cosinc_ar) {
            for (auto ax_out : ax_out_ar) {
                for (auto e_out : e_out_ar) {
                    for (auto ax_in : ax_in_ar) {
                        // std::cout<<"ax_in= "<<ax_in<<" e_in= "<<e_in<<" ci= "<<ci<<" ax_out= "<<ax_out<<" e_out= "<<e_out<<std::endl;
                        OrbitalParam orb[N - 1];
                        orb[0].ax = ax_in;
                        orb[0].ecc = e_in;
                        orb[0].inc = 0.0;
                        orb[0].OMG = dist_0_2p(engine);
                        orb[0].omg = dist_0_2p(engine);
                        orb[0].u = get_eccentric_anomaly_from_mean_anomaly(dist_0_2p(engine), e_in);
                        orb[1].ax = ax_out;
                        orb[1].ecc = e_out;
                        orb[1].inc = acos(ci);
                        orb[1].OMG = dist_0_2p(engine);
                        orb[1].omg = dist_0_2p(engine);
                        orb[1].u = get_eccentric_anomaly_from_mean_anomaly(dist_0_2p(engine), e_out);
                        if (orb[0].ecc == 0.0) orb[0].omg = 0.0;
                        if (orb[1].ecc == 0.0) orb[1].omg = 0.0;
                        if (orb[0].inc == 0.0) orb[0].OMG = 0.0;
                        if (orb[1].inc == 0.0) orb[1].OMG = 0.0;
                        MakeHierarchicalSystem(ptcl, orb, N);
                        // for (int i = 0; i < N; i++) std::cout << "i= " << i << " " << ptcl[i].pos << " " << ptcl[i].vel << std::endl;
                        OrbitalParam orb_check[N - 1];
                        MakeOrbParamFromHierarchicalSystem(orb_check, ptcl, N);

                        for (int i = 0; i < N - 1; i++) {
                            /*
                            std::cout << "in) i= " << i << " ax= " << orb[i].ax << " ecc= " << orb[i].ecc << " inc= " << orb[i].inc
                                      << " OMG= " << orb[i].OMG << " omg= " << orb[i].omg << " u= " << orb[i].u << std::endl;
                            std::cout << "out)i= " << i << " ax= " << orb_check[i].ax << " ecc= " << orb_check[i].ecc << " inc= " << orb_check[i].inc
                                      << " OMG= " << orb_check[i].OMG << " omg= " << orb_check[i].omg << " u= " << orb_check[i].u << std::endl;
                                      */
                            auto diff_ax = orb[i].ax - orb_check[i].ax;
                            auto diff_ecc = orb[i].ecc - orb_check[i].ecc;
                            auto diff_inc = orb[i].inc - orb_check[i].inc;
                            auto diff_OMG = fmod((orb[i].OMG - orb_check[i].OMG), 2.0 * PI);
                            auto diff_omg = fmod((orb[i].omg - orb_check[i].omg), 2.0 * PI);
                            auto diff_u = fmod((orb[i].u - orb_check[i].u), 2.0 * PI);
                            std::cout<<"diff_ax= "<<diff_ax<<" diff_ecc= "<<diff_ecc<<" diff_inc= "<<diff_inc
                                     <<" diff_OMG= "<<diff_OMG<<" diff_omg= "<<diff_omg<<" diff_u= "<<diff_u<<std::endl;
                            assert(fabs(diff_ax) < 1.0e-11);
                            assert(fabs(diff_ecc) < 1.0e-11);
                            assert(fabs(diff_inc) < 1.0e-11);
                            assert(fabs(diff_OMG) < 1.0e-11);
                            assert(fabs(diff_omg) < 1.0e-11);
                            assert(fabs(diff_u) < 1.0e-11);                                     
                        }
                        for (int i = 0; i < N; i++) {
                            ptcl[i].pos = PS::F64vec(0.0);
                            ptcl[i].vel = PS::F64vec(0.0);
                        }
                    }
                }
            }
        }
    }

    /*
    OrbitalParam orb_in, orb_out;
    orb_in.ax = 1.0;
    orb_in.ecc = 0.1;
    orb_in.inc = 0.0;
    orb_in.OMG = 0.0;
    orb_in.omg = 0.0;
    orb_in.u = 0.0;
    orb_out.ax = 10.0;
    orb_out.ecc = 0.7;
    orb_out.inc = 1.0;
    orb_out.OMG = 0.5;
    orb_out.omg = 0.5;
    orb_out.u = 0.5;
    Ptcl ptcl2;
    ptcl2.mass = 1.0;
    MakeHierarchicalSystem(ptcl0, ptcl1, ptcl2, orb_in, orb_out);
    std::cout << "ptcl0.pos= " << ptcl0.pos << " ptcl0.vel= " << ptcl0.vel <<
    std::endl; std::cout << "ptcl1.pos= " << ptcl1.pos << " ptcl1.vel= " <<
    ptcl1.vel << std::endl; std::cout << "ptcl2.pos= " << ptcl2.pos << " ptcl2.vel=
    " << ptcl2.vel << std::endl;

    Ptcl ptcl_cm_in, ptcl_cm_out;
    MakeCMPtcl(ptcl_cm_in, ptcl0, ptcl1);
    std::cout << "ptcl_cm_in.pos= " << ptcl_cm_in.pos << " ptcl_cm_in.vel= " <<
    ptcl_cm_in.vel
              << std::endl;
    MakeCMPtcl(ptcl_cm_out, ptcl_cm_in, ptcl2);
    std::cout << "ptcl_cm_out.pos= " << ptcl_cm_out.pos << " ptcl_cm_out.vel= " <<
    ptcl_cm_out.vel
              << std::endl;

    OrbitalParam orb_in_2, orb_out_2;
    PosVel2OrbParam(orb_in_2, ptcl1, ptcl0);
    PosVel2OrbParam(orb_out_2, ptcl2, ptcl_cm_in);
    std::cout << "orb_in_2.ax= " << orb_in_2.ax << " orb_in_2.ecc= " << orb_in_2.ecc
              << " orb_in_2.inc= " << orb_in_2.inc << " orb_in_2.OMG= " <<
    orb_in_2.OMG
              << " orb_in_2.omg= " << orb_in_2.omg << " orb_in_2.u= " << orb_in_2.u
    << std::endl; std::cout << "orb_out_2.ax= " << orb_out_2.ax << " orb_out_2.ecc=
    " << orb_out_2.ecc
              << " orb_out_2.inc= " << orb_out_2.inc << " orb_out_2.OMG= " <<
    orb_out_2.OMG
              << " orb_out_2.omg= " << orb_out_2.omg << " orb_out_2.u= " <<
    orb_out_2.u
              << std::endl;
*/
    return 0;
}