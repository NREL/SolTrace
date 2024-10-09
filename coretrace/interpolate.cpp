#include <math.h>
#include "interpolate.h"


LUdcmp::LUdcmp(MatDoub& a)
{
    n = (int)a.size();
    lu = a;
    aref = a;
    indx.resize(n);

    const double TINY = 1.0e-40;
    int i, imax, j, k;
    double big, temp;
    VectDoub vv(n);
    d = 1.0;
    for (i = 0; i < n; i++) {
        big = 0.0;
        for (j = 0; j < n; j++)
            if ((temp = fabs(lu.at(i).at(j))) > big) big = temp;
        if (big == 0.0) throw("Singular matrix in LUdcmp");
        vv[i] = 1.0 / big;
    }
    for (k = 0; k < n; k++) {
        big = 0.0;
        for (i = k; i < n; i++) {
            temp = vv[i] * fabs(lu.at(i).at(k));
            if (temp > big) {
                big = temp;
                imax = i;
            }
        }
        if (k != imax) {
            for (j = 0; j < n; j++) {
                temp = lu.at(imax).at(j);
                lu.at(imax).at(j) = lu.at(k).at(j);
                lu.at(k).at(j) = temp;
            }
            d = -d;
            vv[imax] = vv[k];
        }
        indx[k] = imax;
        if (lu.at(k).at(k) == 0.0) lu.at(k).at(k) = TINY;
        for (i = k + 1; i < n; i++) {
            temp = lu.at(i).at(k) /= lu.at(k).at(k);
            for (j = k + 1; j < n; j++)
                lu.at(i).at(j) -= temp * lu.at(k).at(j);
        }
    }
}

void LUdcmp::solve(VectDoub& b, VectDoub& x)
{
    int i, ii = 0, ip, j;
    double sum;
    if (b.size() != static_cast<size_t>(n) ||
        x.size() != static_cast<size_t>(n))
        throw("LUdcmp::solve bad sizes");
    for (i = 0; i < n; i++) x[i] = b[i];
    for (i = 0; i < n; i++) {
        ip = indx[i];
        sum = x[ip];
        x[ip] = x[i];
        if (ii != 0)
            for (j = ii - 1; j < i; j++) sum -= lu.at(i).at(j) * x[j];
        else if (sum != 0.0)
            ii = i + 1;
        x[i] = sum;
    }
    for (i = n - 1; i >= 0; i--) {
        sum = x[i];
        for (j = i + 1; j < n; j++) sum -= lu.at(i).at(j) * x[j];
        x[i] = sum / lu.at(i).at(i);
    }
}

void LUdcmp::solve(MatDoub& b, MatDoub& x)
{
    int i, j, m = (int)b.front().size();
    if (b.size() != static_cast<size_t>(n) ||
        x.size() != static_cast<size_t>(n) ||
        b.front().size() != x.front().size())
        throw("LUdcmp::solve bad sizes");
    VectDoub xx(n);
    for (j = 0; j < m; j++) {
        for (i = 0; i < n; i++) xx[i] = b.at(i).at(j);
        solve(xx, xx);
        for (i = 0; i < n; i++) x.at(i).at(j) = xx[i];
    }
}

void LUdcmp::inverse(MatDoub& ainv)
{
    int i, j;
    ainv.resize(n, VectDoub(n));
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) ainv.at(i).at(j) = 0.;
        ainv.at(i).at(i) = 1.;
    }
    solve(ainv, ainv);
}

double LUdcmp::det()
{
    double dd = d;
    for (int i = 0; i < n; i++) dd *= lu.at(i).at(i);
    return dd;
}

void LUdcmp::mprove(VectDoub& b, VectDoub& x)
{
    int i, j;
    VectDoub r(n);
    for (i = 0; i < n; i++) {
        double sdp = -b[i];
        for (j = 0; j < n; j++)
            sdp += (double)aref.at(i).at(j) * (double)x[j];
        r[i] = sdp;
    }
    solve(r, r);
    for (i = 0; i < n; i++) x[i] -= r[i];
}

//double Powvargram::SQR(const double a) { return a * a; };  // a squared
//
//Powvargram::Powvargram() {};
//
//Powvargram::Powvargram(MatDoub& x, VectDoub& y, const double beta, const double nug)
//{
//    bet = beta;
//    nugsq = nug * nug;
//
//    int i, j, k, npt = (int)x.size(), ndim = (int)x.front().size();
//    double rb, num = 0., denom = 0.;
//    for (i = 0; i < npt; i++) {
//        for (j = i + 1; j < npt; j++) {
//            rb = 0.;
//            for (k = 0; k < ndim; k++) rb += SQR(x.at(i).at(k) - x.at(j).at(k));
//            rb = pow(rb, 0.5 * beta);
//            num += rb * (0.5 * SQR(y[i] - y[j]) - nugsq);
//            denom += SQR(rb);
//        }
//    }
//    alph = num / denom;
//}

//double Powvargram::operator() (const double r) const
//{
//    return nugsq + alph * pow(r, bet);
//}

double GaussMarkov::SQR(const double a) { return a * a; };

double GaussMarkov::pvg_dist(const double r)
{
    return nugsq + alph * pow(r, bet);
}

void GaussMarkov::setup(const double beta, const double nug, const double* err)
{

    //Initialize powvargram parameter alpha
    bet = beta;
    nugsq = nug * nug;

    {
        int i, j, k, npt = (int)x.size(), ndim = (int)x.front().size();
        double rb, num = 0., denom = 0.;
        for (i = 0; i < npt; i++) {
            for (j = i + 1; j < npt; j++) {
                rb = 0.;
                for (k = 0; k < ndim; k++) rb += SQR(x.at(i).at(k) - x.at(j).at(k));
                rb = pow(rb, 0.5 * beta);
                num += rb * (0.5 * SQR(y[i] - y[j]) - nugsq);
                denom += SQR(rb);
            }
        }
        alph = num / denom;
    }


    //create local copies as needed
    npt = (int)x.size();
    ndim = (int)x.front().size();
    dstar.resize(npt + 1);
    vstar.resize(npt + 1);
    v.resize(npt + 1, VectDoub(npt + 1));
    yvi.resize(npt + 1);
    //y.resize(npt + 1);
    y.push_back(0.);
    //---------------

    int i, j;
    for (i = 0; i < npt; i++) {
        //y[i] = yy[i];
        for (j = i; j < npt; j++) {
            v.at(i).at(j) = v.at(j).at(i) = pvg_dist(rdist(&x.at(i), &x.at(j)));
        }
        v.at(i).at(npt) = v.at(npt).at(i) = 1.;
    }
    v.at(npt).at(npt) = y[npt] = 0.;
    if (err) for (i = 0; i < npt; i++) v.at(i).at(i) -= SQR(err[i]);
    vi = new LUdcmp(v);
    vi->solve(y, yvi);
}

GaussMarkov::GaussMarkov() {
    vi = 0;  //initialize null
}

GaussMarkov::~GaussMarkov() {
    if (vi != 0)
        delete vi;
}

double GaussMarkov::interp(VectDoub& xstar) {
    int i;
    for (i = 0; i < npt; i++) vstar[i] = pvg_dist(rdist(&xstar, &x.at(i)));
    vstar[npt] = 1.;
    lastval = 0.;
    for (i = 0; i <= npt; i++) lastval += yvi[i] * vstar[i];
    return lastval;
}

double GaussMarkov::interp(VectDoub& xstar, double& esterr) {
    lastval = interp(xstar);
    vi->solve(vstar, dstar);
    lasterr = 0;
    for (int i = 0; i <= npt; i++) lasterr += dstar[i] * vstar[i];
    esterr = lasterr = sqrt(fmax(0., lasterr));
    return lastval;
}

double GaussMarkov::rdist(VectDoub* x1, VectDoub* x2) {
    double d = 0.;
    for (int i = 0; i < ndim; i++) d += SQR(x1->at(i) - x2->at(i));
    return sqrt(d);
}
