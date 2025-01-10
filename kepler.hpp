#pragma once

#include "matrix3.hpp"
namespace PS = ParticleSimulator;

// 縮退について
// 軌道傾斜角が0度の場合、OMGは任意の値を取ることができるので，OMG=0とする．
// 離心率が0の場合、近点引数は任意の値を取ることができるので，omg=0とする．
#define DEGENERATE_OMG
#define DEGENERATE_omg

// a: semi-major axis
// l: mean anomaly
// e: eccentricity
// u: eccentric anomaly
// n: mean mortion

double keplereq(const double l, const double e, const double u) { return (u - e * sin(u) - l); }

double keplereq_dot(const double e, const double u) { return (1.0 - e * cos(u)); }

// return eccentric anomaly, u
double solve_keplereq_u(const double l, const double e) {
    // std::cout << "l= " << l << " e= " << e << std::endl;
    static double diff_u_crit = 1e-14;
    double u0 = l;
    double u1;
    int loop = 0;
    // std::cout<<"sin(u0)= "<<sin(u0)<<std::endl;
    // std::cout<<"cos(u0)= "<<cos(u0)<<std::endl;
    while (1) {
        loop++;
        double su0 = sin(u0);
        // double cu0 = sqrt(1.0 - su0 * su0);
        double cu0 = cos(u0);
        // u1 = u0 - keplereq(l, e, u0)/keplereq_dot(e, u0);
        u1 = u0 - ((u0 - e * su0 - l) / (1.0 - e * cu0));
        // if(loop % 1000000 == 1) std::cout<<"u0= "<<u0<<" u1= "<<u1<<" diff= "<<fabs(u1 - u0)<<std::endl;
        if (fabs(u1 - u0) < diff_u_crit) {
            return u1;
        } else {
            u0 = u1;
        }
    }
}

// get Pos and Vel after dt from Pos and Vel at current time
template <typename Tptcl>
void DriveKepler(Tptcl &ptcl0, Tptcl &ptcl1, const PS::F64 dt) {
    struct Orb {
        PS::F64 ax;
        PS::F64 ecc;
        PS::F64 inc;
        PS::F64 OMG;
        PS::F64 omg;
        PS::F64 u;
    } orb_now;
    PS::F64 tperi = PosVel2OrbParam(orb_now, ptcl0, ptcl1);
    const PS::F64 n = sqrt((ptcl0.mass + ptcl1.mass) / (orb_now.ax * orb_now.ax * orb_now.ax));
    const PS::F64 l_now = orb_now.u - orb_now.ecc * sin(orb_now.u);  // mean anomaly at now
    const PS::F64 l_new = n * dt + l_now;                            // mean anomaly at new time
    PS::F64 u_new = solve_keplereq_u(l_new, orb_now.ecc);            // eccentric anomaly at new time
    orb_now.u = u_new;
    OrbParam2PosVel(ptcl0, ptcl1, orb_now);
}

double get_eccentric_anomaly_from_mean_anomaly(const double l, const double e) {
    double u = solve_keplereq_u(l, e);
    return u;
}

/*
// not checked yet
double get_eccentric_anomaly_from_true_anomaly(const double f, const double e) {
    double sinu = sin(f) * sqrt(1.0 - e * e) / (1.0 + e * cos(f));
    double cosu = (cos(f) + e) / (1.0 + e * cos(f));
    double u = atan2(sinu, cosu);
    return u;
}

double get_mean_anomaly_from_eccentric_anomaly(const double u, const double e) {
    return (u - e * sin(u));
}

double get_mean_anomaly_from_true_anomaly(const double f, const double e) {
    double u = get_eccentric_anomaly_from_true_anomaly(f, e);
    return get_mean_anomaly_from_eccentric_anomaly(u, e);
}

double get_true_anomaly_from_eccentric_anomaly(const double u, const double e) {
    double sinf = sin(u) * sqrt(1.0 - e * e) / (1.0 - e * cos(u));
    double cosf = (cos(u) - e) / (1.0 - e * cos(u));
    double f = atan2(sinf, cosf);
    return f;
}

double get_true_anomaly_from_mean_anomaly(const double l, const double e) {
    double u = get_eccentric_anomaly_from_mean_anomaly(l, e);
    return get_true_anomaly_from_eccentric_anomaly(u, e);
}
*/

template <typename Tptcl>
Tptcl GetCMPtcl(const Tptcl &ptcl0, const Tptcl &ptcl1) {
    Tptcl ptcl_cm;
    ptcl_cm.mass = ptcl0.mass + ptcl1.mass;
    ptcl_cm.pos = (ptcl0.mass * ptcl0.pos + ptcl1.mass * ptcl1.pos) / ptcl_cm.mass;
    ptcl_cm.vel = (ptcl0.mass * ptcl0.vel + ptcl1.mass * ptcl1.vel) / ptcl_cm.mass;
    return ptcl_cm;
}

////////////////////////////////////////////////////////////////////////////////////////
//////////// CALCULATE PARTICLE POSITION AND VELOCITY FROM ORBITAL PARAMETERS //////////
// calculate relative (reduced) position and velocity
// pos_reduced: relative position defined as pos1-pos0
// vel_reduced: relative velocity defined as vel1-vel0
void OrbParam2PosVelReduced(PS::F64vec &pos_reduced, PS::F64vec &vel_reduced, const PS::F64 mass_tot, const double ax, const double ecc,
                            const double inc, const double OMG, const double omg, const double u) {
    double n = sqrt(mass_tot / (ax * ax * ax));
    double cosu = cos(u);
    double sinu = sin(u);
    double c0 = sqrt(1.0 - ecc * ecc);
    PS::F64vec pos_star(ax * (cosu - ecc), ax * c0 * sinu, 0.0);
    PS::F64vec vel_star(-ax * n * sinu / (1.0 - ecc * cosu), ax * n * c0 * cosu / (1.0 - ecc * cosu), 0.0);
    Matrix3<PS::F64> rot;
    auto OMG_tmp = OMG;
#if defined(DEGENERATE_OMG)
    OMG_tmp = inc == 0.0 ? 0.0 : OMG;
#endif
    auto omg_tmp = omg;
#if defined(DEGENERATE_OMG)
    omg_tmp = ecc == 0.0 ? 0.0 : omg;
#endif
    rot.rotation(inc, OMG_tmp, omg_tmp);
    pos_reduced = rot * pos_star;  // pos1-pos0
    vel_reduced = rot * vel_star;  // pos1-pos0
}
void OrbParam2PosVel(PS::F64vec &pos1, PS::F64vec &vel1, const PS::F64 mass1, const PS::F64vec &pos0, const PS::F64vec &vel0, const PS::F64 mass0,
                     const double ax, const double ecc, const double inc, const double OMG, const double omg, const double u) {
    double mass_tot = mass0 + mass1;
    PS::F64vec pos_reduced, vel_reduced;
    OrbParam2PosVelReduced(pos_reduced, vel_reduced, mass_tot, ax, ecc, inc, OMG, omg, u);
#if 1
    pos1 = pos0 + pos_reduced;
    vel1 = vel0 + vel_reduced;
#else
    PS::F64vec pos_cm = pos0 + (mass1 / m_tot) * pos_reduced;
    PS::F64vec vel_cm = vel0 + (mass1 / m_tot) * vel_reduced;
    pos1 = pos_cm + (m0 / m_tot) * pos_red;
    vel1 = vel_cm + (m0 / m_tot) * vel_red;
#endif
}
// ptcl1: target particle
// ptcl0: source particle whose position and velocity are known
template <typename Tptcl, typename Torb>
void OrbParam2PosVel(Tptcl &ptcl1, const Tptcl &ptcl0, const Torb &orb) {
    OrbParam2PosVel(ptcl1.pos, ptcl1.vel, ptcl1.mass, ptcl0.pos, ptcl0.vel, ptcl0.mass, orb.ax, orb.ecc, orb.inc, orb.OMG, orb.omg, orb.u);
}
// calculate positions and velocities of two particles from orbital parameters under given center of
// mass. ptcl0: source particle ptcl1: target particle ptcl_cm: center of mass orb: orbital
// parameters
template <typename Tptcl, typename Torb>
void OrbParam2PosVel(Tptcl &ptcl1, Tptcl &ptcl0, const Tptcl &ptcl_cm, const Torb &orb) {
    const auto mass_tot = ptcl0.mass + ptcl1.mass;
    PS::F64vec pos_reduced, vel_reduced;
    OrbParam2PosVelReduced(pos_reduced, vel_reduced, mass_tot, orb.ax, orb.ecc, orb.inc, orb.OMG, orb.omg, orb.u);
    ptcl0.pos = ptcl_cm.pos - (ptcl1.mass / mass_tot) * pos_reduced;
    ptcl0.vel = ptcl_cm.vel - (ptcl1.mass / mass_tot) * vel_reduced;
    ptcl1.pos = ptcl_cm.pos + (ptcl0.mass / mass_tot) * pos_reduced;
    ptcl1.vel = ptcl_cm.vel + (ptcl0.mass / mass_tot) * vel_reduced;
}

////////////////////////////////////////////////////////////////////////////////////////
//////////// CALCULATE ORBITAL PARAMETERS FROM PARTICLE POSITION AND VELOCITY //////////
// u is eccentric anomaly
// return the time passed from the peri-center passage
double PosVel2OrbParam(double &ax, double &ecc, double &inc, double &OMG, double &omg, double &u, const PS::F64vec &pos1, const PS::F64vec &vel1,
                       const PS::F64 mass1, const PS::F64vec &pos0, const PS::F64vec &vel0, const PS::F64 mass0) {
    double m_tot = mass0 + mass1;
    PS::F64vec pos_red = pos1 - pos0;
    PS::F64vec vel_red = vel1 - vel0;
    double r_sq = pos_red * pos_red;
    double r = sqrt(r_sq);
    double inv_dr = 1.0 / r;
    double v_sq = vel_red * vel_red;
    ax = 1.0 / (2.0 * inv_dr - v_sq / m_tot);
    PS::F64vec AM = pos_red ^ vel_red;
    //std::cout << "AM= " << AM << std::endl;
    // sinOMG : AM.x / h*sinI, cosOMG : -AM.y / h*sinI
    inc = atan2(sqrt(AM.x * AM.x + AM.y * AM.y), AM.z);
    OMG = atan2(AM.x, -AM.y);
    #if defined(DEGENERATE_OMG)
    if(fabs(inc) < 1e-15) OMG = 0.0;
    #endif
    PS::F64vec pos_bar, vel_bar;
    double cosOMG = cos(OMG);
    double sinOMG = sin(OMG);
    double cosinc = cos(inc);
    double sininc = sin(inc);
    pos_bar.x = pos_red.x * cosOMG + pos_red.y * sinOMG;
    pos_bar.y = (-pos_red.x * sinOMG + pos_red.y * cosOMG) * cosinc + pos_red.z * sininc;
    pos_bar.z = 0.0;
    vel_bar.x = vel_red.x * cosOMG + vel_red.y * sinOMG;
    vel_bar.y = (-vel_red.x * sinOMG + vel_red.y * cosOMG) * cosinc + vel_red.z * sininc;
    vel_bar.z = 0.0;
    double h = sqrt(AM * AM);
    double ecccosomg = h / m_tot * vel_bar.y - pos_bar.x * inv_dr;
    double eccsinomg = -h / m_tot * vel_bar.x - pos_bar.y * inv_dr;
    ecc = sqrt(ecccosomg * ecccosomg + eccsinomg * eccsinomg);
    omg = atan2(eccsinomg, ecccosomg);
    #if defined(DEGENERATE_omg)
    if(fabs(ecc) < 1e-15) omg = 0.0;
    #endif    
    double phi = atan2(pos_bar.y, pos_bar.x);  // f + omg (f: true anomaly)
    double f = phi - omg;
    double sinu = r * sin(f) / (ax * sqrt(1.0 - ecc * ecc));
    double cosu = (r * cos(f) / ax) + ecc;
    u = atan2(sinu, cosu);                  // eccentric anomaly
    double n = sqrt(m_tot / ax * ax * ax);  // mean mortion
    double l = u - ecc * sin(u);            // mean anomaly
    auto t_from_peri = l / n;
    OMG = OMG > 0.0 ? OMG : 2.0 * M_PI + OMG;
    omg = omg > 0.0 ? omg : 2.0 * M_PI + omg;
    u = u > 0.0 ? u : 2.0 * M_PI + u;
    return t_from_peri;  // time passed from the peri-center passage
}

// return the time passed from the peri-center passage
template <typename Tptcl, typename Torb>
double PosVel2OrbParam(Torb &orb, const Tptcl &ptcl1, const Tptcl &ptcl0) {
    return PosVel2OrbParam(orb.ax, orb.ecc, orb.inc, orb.OMG, orb.omg, orb.u, ptcl1.pos, ptcl1.vel, ptcl1.mass, ptcl0.pos, ptcl0.vel, ptcl0.mass);
}

/////////////////////////////////////////////////
//////////// MAKE HIERARCHICAL SYSTEMS //////////
// make hierarchical systems from orbital parameters
// the origin of the coordinate is the center of mass of the system
// ptcl[n_ptcl]
// orb[n_ptcl-1]
template <typename Tptcl, typename Torb>
void MakeHierarchicalSystem(Tptcl ptcl[], const Torb orb[], const int n_ptcl) {
    Tptcl ptcl_cm_in, ptcl_cm_out;
    ptcl_cm_in.mass = ptcl[0].mass;
    ptcl_cm_in.pos = PS::F64vec(0.0);
    ptcl_cm_in.vel = PS::F64vec(0.0);
    for (int i = 1; i < n_ptcl; i++) {
        OrbParam2PosVel(ptcl[i], ptcl_cm_in, orb[i - 1]);
        Tptcl ptcl_cm_out = GetCMPtcl(ptcl[i], ptcl_cm_in);
        for (int j = 0; j <= i; j++) {
            ptcl[j].pos -= ptcl_cm_out.pos;
            ptcl[j].vel -= ptcl_cm_out.vel;
        }
        ptcl_cm_in.mass += ptcl[i].mass;
    }
}

template <typename Tptcl, typename Torb>
void MakeOrbParamFromHierarchicalSystem(Torb orb[], const Tptcl ptcl[], const int n_ptcl) {
    Tptcl ptcl_cm = ptcl[0];
    for (int i = 0; i < n_ptcl - 1; i++) {
        // PosVel2OrbParam(orb[i], ptcl_cm, ptcl[i + 1]);
        PosVel2OrbParam(orb[i], ptcl[i + 1], ptcl_cm);
        ptcl_cm = GetCMPtcl(ptcl_cm, ptcl[i + 1]);
    }
}

/*
void PosVel2AxEccInc(double &ax, double &ecc, double &inc, const PS::F64vec &pos0,
                     const PS::F64vec &pos1, const PS::F64vec &vel0, const PS::F64vec &vel1,
                     const double mass0, const double mass1) {
    double m_tot = mass0 + mass1;
    PS::F64vec pos_red = pos1 - pos0;
    PS::F64vec vel_red = vel1 - vel0;
    double r_sq = pos_red * pos_red;
    double r = sqrt(r_sq);
    double inv_dr = 1.0 / r;
    double v_sq = vel_red * vel_red;
    ax = 1.0 / (2.0 * inv_dr - v_sq / m_tot);
    PS::F64vec AM = pos_red ^ vel_red;
    inc = atan2(sqrt(AM.x * AM.x + AM.y * AM.y), AM.z);
    PS::F64 OMG = atan2(AM.x, -AM.y);
    PS::F64vec pos_bar, vel_bar;
    double cosOMG = cos(OMG);
    double sinOMG = sin(OMG);
    double cosinc = cos(inc);
    double sininc = sin(inc);
    pos_bar.x = pos_red.x * cosOMG + pos_red.y * sinOMG;
    pos_bar.y = (-pos_red.x * sinOMG + pos_red.y * cosOMG) * cosinc + pos_red.z * sininc;
    pos_bar.z = 0.0;
    vel_bar.x = vel_red.x * cosOMG + vel_red.y * sinOMG;
    vel_bar.y = (-vel_red.x * sinOMG + vel_red.y * cosOMG) * cosinc + vel_red.z * sininc;
    vel_bar.z = 0.0;
    double h = sqrt(AM * AM);
    double ecccosomg = h / m_tot * vel_bar.y - pos_bar.x * inv_dr;
    double eccsinomg = -h / m_tot * vel_bar.x - pos_bar.y * inv_dr;
    ecc = sqrt(ecccosomg * ecccosomg + eccsinomg * eccsinomg);
}

// return semi major axis
double PosVel2Ax(const PS::F64vec &pos0, const PS::F64vec &pos1, const PS::F64vec &vel0,
                 const PS::F64vec &vel1, const double mass0, const double mass1) {
    // static const PS::F64 PI = 4.0 * atan(1.0);
    double m_tot = mass0 + mass1;
    // double m_red = (mass0*mass1)/m_tot;
    PS::F64vec pos_red = pos1 - pos0;
    PS::F64vec vel_red = vel1 - vel0;
    double r_sq = pos_red * pos_red;
    double r = sqrt(r_sq);
    double inv_dr = 1.0 / r;
    double v_sq = vel_red * vel_red;
    double ax = 1.0 / (2.0 * inv_dr - v_sq / m_tot);
    return ax;
}
*/