/*
 * collatz_gmp.h
 * NonintegerCollatz: Collatz for non-integers
 *
 * Copyright (c) 2026 Bryan Franklin. All rights reserved.
 */
#ifndef COLLATZ_GMP_H
#define COLLATZ_GMP_H
#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <set>

static inline int init_v(std::vector<int>& v) {
    for(int i=v.size()-1; i>=0; --i) {
        v[i] = i;
    }
    return 1;
}

static const int D_MAX = 32;
static inline int next_v(std::vector<int>& v) {
    for(int i=0; i<(v.size()-1); ++i) {
        if( v[i] < v[i+1]-1 ) {
            v[i] += 1;
            return 1;
        } else {
            v[i] = i;
        }
    }
    if( v[v.size()-1] < D_MAX ) {
        v[v.size()-1] += 1;
        return 1;
    }

    return 0;
}

int init_vf(std::vector<double>& vf) {
    for(int i=vf.size()-1; i>=0; --i) {
        vf[i] = i;
    }
    return 1;
}


std::ostream& operator<<(std::ostream& ofs, std::vector<int> v) {
    ofs << "{";
    for(int i=0; i<v.size(); ++i) {
        ofs << v[i];
        if( i < v.size()-1 )
            ofs << " ";
    }
    ofs << "}";

    return ofs;
}


std::ostream& operator<<(std::ostream& ofs, std::vector<double> v) {
    ofs << "{";
    for(int i=0; i<v.size(); ++i) {
        ofs << v[i];
        if( i < v.size()-1 )
            ofs << " ";
    }
    ofs << "}";

    return ofs;
}


size_t lob_pos(const mpz_t& x) {
    // See section 5.15 of the GMP docs.
    return mpz_scan1(x, 0);
}


size_t hob_pos(const mpz_t& x) {
    // See section 5.15 of the GMP docs.
    return mpz_sizeinbase(x, 2)-1;
}


size_t hob_pos(const mpz_t& x, mpz_t &pos) {
    size_t length = hob_pos(x);
    mpz_set_si(pos, length);
    return length;
}


size_t gap(const mpz_t& x) {
    size_t length = hob_pos(x) - lob_pos(x);
    return length;
}


size_t gap(const mpz_t& x, mpz_t &gap_size) {
    size_t length = hob_pos(x) - lob_pos(x);
    mpz_set_si(gap_size, length);
    return length;
}


static inline std::string asBitString(mpz_t x_0, int bits=-1, int lowest_bit=0) {
    mpz_t x_j;
    mpz_init_set(x_j, x_0);

    if( mpz_cmp_ui(x_j, 0) == 0 ) return "0";
    if( bits<0 )    bits = hob_pos(x_j)+1;

    if( lowest_bit > 0 && lowest_bit < (hob_pos(x_j)) )
        mpz_tdiv_q_2exp(x_j, x_j, lowest_bit);  // right shift by lowest_bit

    std::string ret = "";
    while( mpz_cmp_ui(x_j, 0) > 0 && ret.length()<bits ) {
        int bit = mpz_tstbit(x_j, 0);
        ret = std::to_string(bit) + ret;
        mpz_tdiv_q_2exp(x_j, x_j, 1);  // right shift by 1 (one)
    }
    while( ret.length()<bits ) { ret = "0" + ret; }

    return ret;
}


static inline int v_to_z(const std::vector<int>& v, mpz_t &z) {
    mpz_init_set_ui(z, 0);
    mpz_t pow2, pow3, prod;
    mpz_init(pow2);
    mpz_init(pow3);
    mpz_init(prod);
    int m = v.size()-1;
    for(int i=0; i<v.size(); ++i) {
        // Compute z += 2^v[i] * 3^i
        mpz_ui_pow_ui(pow2, 2, v[i]);
        if( m > (i+1) )
            mpz_ui_pow_ui(pow3, 3, m-i-1);
        else
            mpz_set_ui(pow3, 1);    // last element is just 2^d
        mpz_mul(prod, pow2, pow3);
        mpz_add(z, z, prod);
    }

    return 0;
}


static inline int v_to_w(const std::vector<int>& v, mpz_t &w) {
    mpz_set_ui(w, 0);
    mpz_t bit;
    mpz_init(bit);
    for(int i=0; i<v.size(); ++i) {
        mpz_ui_pow_ui(bit, 2, v[i]);
        mpz_add(w, w, bit);
    }
    mpz_clear(bit);
    return 0;
}


static inline int x_to_w(const mpz_t x, mpz_t &w) {
    mpz_t x_i;
    mpz_t bit;

    if( mpz_cmp_ui(x, 0) <= 0 ) {
        // w isn't defined for x_0==0
        mpz_set_ui(w, -1);
        return -1;
    }

    mpz_set_ui(w, 0);
    mpz_init(bit);
    int divs = 0;
    mpz_init_set(x_i, x);
    // count divisions to get back to an odd
    while( mpz_divisible_ui_p(x_i, 2) ) {
        ++divs;
        mpz_tdiv_qr_ui(x_i, bit, x_i, 2);
    }
    mpz_ui_pow_ui(bit, 2, divs);
    mpz_add(w, w, bit);
    while( mpz_cmp_ui(x_i, 1) > 0 ) {
        // apply odd rule
        mpz_mul_ui(x_i, x_i, 3);
        mpz_add_ui(x_i, x_i, 1);

        // count divisions to get back to odd
        while( mpz_divisible_ui_p(x_i, 2) ) {
            ++divs;
            mpz_tdiv_qr_ui(x_i, bit, x_i, 2);
        }

        // add a bit to w
        mpz_ui_pow_ui(bit, 2, divs);
        mpz_add(w, w, bit);
    }
    mpz_clear(bit);

    return divs;
}


static inline int v_to_x(const std::vector<int> &v, mpq_t &xq, mpf_t &xf) {
    if( v.size() < 1 ) {
        std::cerr << "v must have at least one element." << std::endl;
        mpq_set_si(xq, -1, 1);
        mpf_set_d(xf, -1.0);
        return -1;
    }

    // attempt to compute x_0 from v vector
    mpz_t z, pow2, pow3, diff;
    mpz_init(z);
    mpz_init(pow2);
    mpz_init(pow3);
    mpz_init(diff);
    mpq_t num, den;
    mpq_init(num);
    mpq_init(den);

    v_to_z(v,z);
    int m = v.size()-1;
    int d = v[v.size()-1];

    mpz_ui_pow_ui(pow2, 2, d+1);
    mpz_sub(diff, pow2, z);

    mpz_ui_pow_ui(pow3, 3, m);

    mpq_set_z(num, diff);
    mpq_set_z(den, pow3);
    mpq_div(xq, num, den);
    mpf_set_q(xf, xq);
    mpq_clear(num);
    mpq_clear(den);

    mpz_clear(pow2);
    mpz_clear(pow3);
    mpz_clear(diff);
    mpz_clear(z);

    return 0;
}


static inline int w_to_v(const mpz_t &w, std::vector<int> &v) {

    int total_divs = 0;
    mpz_t w_i;
    mpz_init_set(w_i, w);
   
    v.clear();

    // deal with even inputs
    while( mpz_divisible_ui_p(w_i, 2) 
        && mpz_cmp_ui(w_i, 0) > 0 ) {
        //divide by 2
        mpz_divexact_ui(w_i, w_i, 2);
        ++total_divs;
    }
    // record initial number of divisions
    // (i.e., record the initial lob)
    v.push_back(total_divs);

    // start collatz loop (i.e., until x_i == 1)
    while( mpz_cmp_ui(w_i, 1) > 0 ) {

        // remove bit
        mpz_sub_ui(w_i, w_i, 1);

        // count even steps to get back to odd
        while( mpz_divisible_ui_p(w_i, 2) ) {
            // divide by 2
            mpz_divexact_ui(w_i, w_i, 2);
            ++total_divs;
        }

        // record number of divisions to return x_i to odd
        v.push_back(total_divs);
    }
    // final element of v will be d (of 2^d)

    return v.size()-1;
}


static inline int w_to_m(const mpz_t &w, mpz_t &m) {
    std::vector<int> v;

    w_to_v(w, v);
    int m_val = v.size()-1;
    mpz_set_ui(m, m_val);

    return m_val;
}


static inline int w_to_d(const mpz_t &w, mpz_t &d) {
    std::vector<int> v;

    w_to_v(w, v);
    int d_val = v[v.size()-1];
    mpz_set_ui(d, d_val);

    return d_val;
}


static inline int w_to_xq(const mpz_t &w, mpq_t &xq) {

    if( mpz_cmp_ui(w, 1) < 0 ) {
        std::cerr << "w must have at least one bit." << std::endl;
        mpq_set_si(xq, -1, 1);
        return -1;
    }

    // attempt to compute x_0 from v vector
    mpz_t z, pow2, pow3, diff;
    mpz_init(z);
    mpz_init(pow2);
    mpz_init(pow3);
    mpz_init(diff);
    mpq_t num, den;
    mpq_init(num);
    mpq_init(den);
    std::vector<int> v;

    w_to_v(w, v);
    v_to_z(v,z);
    int m = v.size()-1;
    int d = v[v.size()-1];

    // get 2^(d+1)
    mpz_ui_pow_ui(pow2, 2, d+1);
    mpz_sub(diff, pow2, z);

    // get 3^m
    mpz_ui_pow_ui(pow3, 3, m);

    // Compute: x_0 = (2^d - z) / 3^m
    mpq_set_z(num, diff);
    mpq_set_z(den, pow3);
    mpq_div(xq, num, den);

    mpq_clear(num);
    mpq_clear(den);
    mpz_clear(pow2);
    mpz_clear(pow3);
    mpz_clear(diff);
    mpz_clear(z);

    return 0;
}


static inline int w_to_xf(const mpz_t &w, mpf_t &xf) {
    if( mpz_cmp_ui(w, 1) < 0 ) {
        std::cerr << "w must have at least one bit." << std::endl;
        mpf_set_d(xf, -1.0);
        return -1;
    }

    // get exact quatient
    mpq_t xq;
    mpq_init(xq);
    w_to_xq(w, xq);

    // convert it to floating point
    mpf_set_q(xf, xq);
    mpq_clear(xq);

    return 0;
}


static inline int w_to_z(const mpz_t& w, mpz_t &z) {
    std::vector<int> v;
    w_to_v(w, v);
    v_to_z(v, z);

    char *str = mpz_get_str(NULL, 2, w);
    std::cout << __FUNCTION__ << " w->z: " << str << " -> " << z << std::endl;
    free(str); str=NULL;

    return 0;
}


int x_to_v(const mpz_t x, std::vector<int> &v) {
    // use collatz rules to build a v
    mpz_t x_i;
    mpz_t one;
    mpz_init(x_i);
    mpz_init(one);
    mpz_set_ui(one, 1);
    mpz_set(x_i, x);
    int total_divs = 0;
    v.clear();

    // deal with even inputs
    while( mpz_divisible_ui_p(x_i, 2) ) {
        //divide by 2
        mpz_divexact_ui(x_i, x_i, 2);
        ++total_divs;
    }
    // record initial number of divisions
    // (i.e., record the initial lob)
    v.push_back(total_divs);

    // start collatz loop (i.e., until x_i == 1)
    while( mpz_cmp_ui(x_i, 1) > 0 ) {

        // apply 3x+1 collatz rule for odds
        mpz_mul_ui(x_i, x_i, 3);
        mpz_add_ui(x_i, x_i, 1);

        // count even steps to get back to odd
        while( mpz_divisible_ui_p(x_i, 2) ) {
            // divide by 2
            mpz_divexact_ui(x_i, x_i, 2);
            ++total_divs;
        }

        // record number of divisions to return x_i to odd
        v.push_back(total_divs);
    }
    // final element of v will be d (of 2^d)

    return v.size()-1;
}


int x_to_vf(double x_0, std::vector<double> &vf, int max_d=10, int max_m=10) {
    for(int m=1; m<=max_m; ++m) {
        vf.resize(m);
        init_vf(vf);

        // set up and 'solve' multi-variate optimatization to find v using fnt
    }

    return 0;
}


int vf_to_x(std::vector<double> vf, double &x_0) {
    // compute x_0 directly from vf, m, and d.
    return 0;
}


int write_dataset4(std::string filename, int min, int max) {
    // open file
    std::ofstream ofs;
    ofs.open(filename, std::ofstream::out | std::ofstream::trunc);
    std::cout << "Writing " << filename << std::endl;

    int steps = (max - min) * 10;
    double x_step = (max - min) / (double)steps;
    double x_min = min;
    std::vector<int> v;
    mpz_t w, w2;
    mpz_t x_0;
    mpz_init(w);
    mpz_init(w2);
    mpz_init(x_0);
    int mults, divs;
    for(int i=min; i<max; ++i) {
        //double x_0 = x_min + x_step * i;
        mpz_set_si(x_0, i);
        x_to_v(x_0, v);
        v_to_w(v, w);
        x_to_w(x_0, w2);

        multsAndDivsToOne(i, mults, divs);

        double w_dbl = mpz_get_d(w);
        double w2_dbl = mpz_get_d(w2);
        ofs << i << " "
            << w_dbl << " "
            << divs << " "
            << mults << " "
            << w2_dbl << " "
            << std::endl;
    }
    mpz_clear(w);

    ofs.close();

    return 0;
}


int find_d(const mpz_t &z) {

    // find d (aka hob_pos(z))
    int k = 1;
    mpz_t pow2;
    mpz_init(pow2);
    mpz_ui_pow_ui(pow2, 2, k);
    while( mpz_cmp(pow2, z) < 0 ) {
        k *= 2;
        std::cout << "k=" << k << std::endl;
        mpz_ui_pow_ui(pow2, 2, k);
    }
    std::cout << "Looking for d, k=" << k << std::endl;
    int lower = k/2;
    int upper = k;
    while( upper > lower+1 ) {
        int mid = upper/2 + lower/2;
        std::cout << "lower=" << lower
            << "; mid=" << mid
            << "; upper=" << upper
            << std::endl;
        mpz_ui_pow_ui(pow2, 2, mid);
        int cmp = mpz_cmp(pow2, z);
        if( cmp > 0 ) {
            upper = mid;
        }
        if( cmp < 0 ) {
            lower = mid;
        }
    }
    k = lower;
    mpz_clear(pow2);
    std::cout << "d=" << k << std::endl;

    return k;
}


int get_candidates(const mpz_t &z_0, std::set<int> &candidates) {

    mpz_t sub, quot, rem;
    mpz_t pow2;
    mpz_t z_i;
    mpz_init(sub);
    mpz_init(quot);
    mpz_init(rem);
    mpz_init_set_ui(pow2, 0);
    mpz_init_set(z_i, z_0);
    candidates.clear();
    int k = -1;
    while( mpz_cmp(pow2, z_i) <= 0 ) {
        k += 1;

        if( k < 0 )    continue;
        mpz_ui_pow_ui(pow2, 2, k);
        if( mpz_cmp(pow2, z_i) == 0 ) {
            candidates.clear();
            candidates.insert(k);
            break;
        }
        if( mpz_cmp(pow2, z_i) > 0 )   continue;
        mpz_sub(sub, z_i, pow2);
        mpz_tdiv_qr_ui(quot, rem, sub, 3);
        if( mpz_cmp_ui(rem, 0) != 0 )  continue;

        candidates.insert(k);
    }

    mpz_clear(pow2);
    mpz_clear(sub);
    mpz_clear(quot);
    mpz_clear(rem);

    return candidates.size();
}


// TODO: add variables to help compute the branching factor
static long long sum_depths = 0;
static long long num_depths = 0;
static long long max_depth = 0;

static int z_to_v_helper(const mpz_t &z_0, std::vector<int> &v, int k_0, int
depth) {

    // update z_i
    // z_i = z_i - 2^k
    mpz_t z_i;
    mpz_t pow2;
    mpz_t sub, quot, rem;
    mpz_init(z_i);
    mpz_init(pow2);
    mpz_init(sub);
    mpz_init(quot);
    mpz_init(rem);

    // update z using lob = 2^k
    mpz_ui_pow_ui(pow2, 2, k_0);
    mpz_sub(sub, z_0, pow2);
    mpz_tdiv_qr_ui(z_i, rem, sub, 3);
    if( mpz_cmp_ui(rem, 0) != 0 ) {
        gmp_printf("%Zd not divisible by 3, rem=%Zd.\n", sub, rem);
        exit(1);
    }
    if( mpz_cmp_ui(z_i, 0) == 0 ) {
        std::cout << k_0 << " got to zero." << std::endl;
        return 1;
    }

    // find new candidates
    std::set<int> candidates;
    get_candidates(z_i, candidates);

    // if instead of making recursive calls, candidates were just added to an
    // overall queue that gets sorted by distance from some value, this could
    // be A* search.

    // TODO: perform look-ahead and sort candidates by resulting candidates

    // recursively check each candidate
    for(int k : candidates) {
        if( k >= k_0 )
            continue;

        //std::cout << "recursively checking candidate k=" << k << std::endl;
        if( z_to_v_helper(z_i, v, k, depth+1) == 1 ) {
            v.push_back(k);
            std::cout << k << " seems to work." << std::endl;
            return 1;
        }
    }

    mpz_clear(z_i);
    mpz_clear(pow2);
    mpz_clear(sub);
    mpz_clear(quot);
    mpz_clear(rem);

    sum_depths += depth;
    ++num_depths;

    if( depth > max_depth ) {
        max_depth = depth;
        std::cout << "Max depth: " << max_depth << std::endl;
    }
    if( (num_depths%(1<<20)) == 0 ) {
        double avg_depth = sum_depths / (double)num_depths;
        std::cout << "avg depth: " << avg_depth << std::endl;
    }

    return 0;
}


int z_to_v_direct(const mpz_t &z, std::vector<int> &v) {

    int d = find_d(z);

    v.clear();

    // z_i = z_i - 2^k
    mpz_t z_0;
    mpz_t pow2;
    mpz_t sub, quot, rem;
    mpz_init(z_0);
    mpz_init(pow2);
    mpz_init(sub);
    mpz_init(quot);
    mpz_init(rem);

    mpz_ui_pow_ui(pow2, 2, d);
    // y_0 = z - 2^d
    mpz_sub(z_0, z, pow2);
    if( mpz_cmp_ui(sub, 0) == 0 )  {
        v.clear();
    }

    for(int k=d-1; k>=0; --k) {
        std::cout << "checking top-level candidate k=" << k << std::endl;

        mpz_ui_pow_ui(pow2, 2, k);

        // compute largest power of 2 (2^k)
        //      such that (2^k < y_i) ^ (y_i - 2^k) % 3 == 0)
        if( mpz_cmp(pow2, z_0) >= 0 )   continue;

        mpz_sub(sub, z_0, pow2);
        if( mpz_cmp_ui(sub, 0) == 0 )  {
            v.clear();
            v.push_back(k);
        }

        mpz_tdiv_qr_ui(quot, rem, sub, 3);
        if( mpz_cmp_ui(rem, 0) != 0 )  continue;

        if( z_to_v_helper(z_0, v, k, 0) == 1 ) {
            // v is valid
            std::cout << "valid v: " << v << std::endl;
            v.push_back(k);
            std::cout << k << " seems to work." << std::endl;
            v.push_back(d);
            return 1;
        }
    }
    std::cout << "Failed to find a valid v." << std::endl;

    // success if y_i reaches 1
    if( mpz_cmp_ui(z_0, 1) == 0 )  return 1;

    mpz_clear(z_0);
    mpz_clear(pow2);
    mpz_clear(sub);
    mpz_clear(quot);
    mpz_clear(rem);

    // if y_i can't reach 1, return failure
    return 0;
}

#endif // COLLATZ_GMP_H
