/*
 * collatz_dbl.cpp
 * NonintegerCollatz: Collatz for non-integers
 *
 * Copyright (c) 2026 Bryan Franklin. All rights reserved.
 */
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdio.h>
#include "collatz.h"
#include "collatz_gmp.h"


// WARNING: This code is specific to x87's 80 bit extended_precision registers!
// see: https://en.wikipedia.org/wiki/Extended_precision

// to run:
// clear; CXXFLAGS=-O3\ -I/opt/local/include\ -L/opt/local/lib\ -lgmp\ -lgmpxx make collatz_dbl && time ./collatz_dbl && time ./plot_dbl.gplot


typedef union {
    long double x;
    unsigned char c[1];
    unsigned long long ll;
} ld2char_t;


int sign_bit(long double x) {
    ld2char_t ld;
    memset(&ld, '\0', sizeof(ld));
    ld.x = x;

    int bit = (ld.c[9] >> 7) & 1;;
    return bit;
}


unsigned long long significand(long double x) {
    ld2char_t ld;
    memset(&ld, '\0', sizeof(ld));
    ld.x = x;

    ld.c[7] |= 0x80;
    unsigned long long ret = ld.ll;

    return ret;
}


int exponent(long double x) {
    ld2char_t ld;
    memset(&ld, '\0', sizeof(ld));
    ld.x = x;

    int expon = ((ld.c[9] & 127) << 8) | ld.c[8];
    expon -= 16383;

    return expon;
}


int get_bit(unsigned char *buf, size_t n, int bit) {
    if( bit >= n )  return -1;

    return ((buf[(bit/8)]) >> (bit%8)) & 1;
}


int set_bit(unsigned char *buf, size_t n, int bit, int to) {
    if( bit >= n )  return -1;

    if( to == 0 ) {
        unsigned char ch = buf[bit/8];
        //std::cout << "Setting bit " << bit << " to 0; ch: " << ((int)ch) << std::endl;
        ch |= 1<<(bit%8);
        ch ^= 1<<(bit%8);
        buf[bit/8] = ch;
        //std::cout << "Set bit to 0; ch: " << ((int)ch) << std::endl;
    } else if( to == 1 ) {
        unsigned char ch = buf[bit/8];
        //std::cout << "Setting bit " << bit << " to 1; ch: " << ((int)ch) << std::endl;
        ch |= 1<<(bit%8);
        buf[bit/8] = ch;
        //std::cout << "Set bit to 1; ch: " << ((int)ch) << std::endl;
    } else {
        std::cerr << "Bit must be 0 or 1." << std::endl;
        return -1;
    }

    return -1;
}


long double set_exponent(long double x, int to) {
    ld2char_t ld;
    memset(&ld, '\0', sizeof(ld));
    ld.x = x;

    // clear old exponent
    ld.c[9] &= 0x80;
    ld.c[8] = 0;

    to += 16383;
    // set new bits
    ld.c[9] |= (to&0x7f00) >> 8;
    ld.c[8] |= to&0xff;
    
    return ld.x;
}


std::string long_double_as_hex(long double x) {
    ld2char_t ld;
    memset(&ld, '\0', sizeof(ld));
    ld.x = x;
    
    std::string ret = "";
    char asHex[4];
    for(int i=0; i<10; ++i) {
        snprintf(asHex, sizeof(asHex), "%02x", ld.c[i]);
        ret = std::string(asHex) + ret;
    }

    return "0x" + ret;
}


std::string long_double_as_binary(long double x) {
    ld2char_t ld;
    memset(&ld, '\0', sizeof(ld));
    ld.x = x;
    
    std::string ret = "";
    int num_bits = 80;
    for(int i=0; i<num_bits; ++i) {
        if( i==64 || i==78 )
            ret = " " + ret;
        int bit = get_bit(ld.c, sizeof(ld)*8, i);
        ret = std::to_string(bit) + ret;
    }

    return ret;
}


int lob_pos(long double x) {
    ld2char_t ld;
    memset(&ld, '\0', sizeof(ld));
    ld.x = x;
    
    // find lob bit
    int pos = -1;
    for(int i=0; i<113 && pos<0; ++i) {
        if( get_bit(ld.c, sizeof(ld)*8, i) == 1 ) {
            pos = i;
            break;
        }
    }
    return pos;
}


int gap(long double x) {
    // find lob bit
    int pos = lob_pos(x);
    // compute gap
    int gap = 63-pos;

    return gap;
}


long double lob_value(long double x) {
    ld2char_t ld;
    memset(&ld, '\0', sizeof(ld));
    ld.x = x;
    
    // find lob bit
    int pos = lob_pos(x);

    // compute number of bits to shift exponent
    int delta_exp = 63-pos;

    // extract exponent
    int expon = exponent(x);

    // update exponent
    //std::cout << "delta_exp: " << delta_exp << std::endl;
    expon -= delta_exp;

    // check that result is valid
    if( expon > 16383 || expon < -16383 )
        return -1.0;

    // place exponent back
    ld.x = set_exponent(1, expon);
    set_bit(ld.c, 63, 64, 1);

    // convert back to long double
    long double ret = ld.x;
    
    return ret;
}


long double collatz_dbl(const long double x) {
    long double lob = lob_value(x);
    if( lob < 0.0 ) return lob;

    if( lob <= 1 )
        return 3*x + lob;
    else
        return x/2;

    return -1.0;
}


void enumerate_dataset_double(std::ostream &fout, int max_x, double step_size) {
    long double x_dbl = -1.0;

    std::cout << "Enumerating up to " << max_x << " with step size " << step_size << "." << std::endl;
    long count = 0;
    long inv_step = (long)floor(1.0/step_size);

    for(long long i=0; i<max_x*inv_step; ++i) {
        long double x_0 = i * step_size;
        long m=0, d=0;
        if( (++count%1000) == 0 ) {
            std::cout << "." << std::flush;
            if( ((count/1000)%70) == 0 ) {
                std::cout << " " << std::setprecision(3)
                          << (100.0*count/(max_x/step_size)) << "%" 
                          << std::endl;
            }
        }

        // copy bits into significand of x_dbl;
        x_dbl = x_0;

        if( x_dbl == floor(x_dbl) )
            continue;

        // perform collatz operations to get m, d.
        while( x_dbl != 1.0 ) {
            long double lob = lob_value(x_dbl);
            if( lob <= 1 ) {
                ++m;
                x_dbl = 3*(x_dbl) + lob;
            } else {
                ++d;
                x_dbl = (x_dbl)/2;
            }
        }

        fout << std::setprecision(21)
            << -1.0 << " "
            << d << " "
            << m << " "
            << x_0 << " "
            << 1.0 << " "
            << std::endl;
    }
    std::cout << " 100.0%" << std::endl;
}


void sample_dataset_double(std::ostream &fout, int max_bits, int max_exp, int n_samples) {
    
    gmp_randclass rng(gmp_randinit_default);
    mpz_class rnd_value;
    mpz_t rnd_mpz;
    long double *x_ptr = NULL;
    x_ptr = (long double*)calloc(1, sizeof(*x_ptr));

    mpz_init(rnd_mpz);
    std::set<long double> samples;
    for(long n_bits=1; n_bits<=max_bits; ++n_bits) {
        if( n_samples > (1<<n_bits) ) {
            for(int i=0; i<(1<<n_bits); ++i) {
                long double x_0 = i / (long double)(1<<n_bits);
                // skip duplicates
                if( samples.count(x_0) > 0 )
                    continue;
                samples.insert(x_0);
            }
        } else {
            // choose samples randomly up to n_bits
            while( samples.size() < n_bits*n_samples ) {

                rnd_value = rng.get_z_bits(n_bits);
                mpz_set(rnd_mpz, rnd_value.get_mpz_t());

                // copy bits into significand of *x_ptr;
                *x_ptr = 1.0;
                for(int i=0; i<n_bits; ++i) {
                    int bit_value = mpz_tstbit(rnd_mpz, i);
                    set_bit((unsigned char*)x_ptr, 8*sizeof(*x_ptr), 62-i, bit_value);
                }
                long expon = exponent(*x_ptr);
                expon = lrand48() % (2*max_exp) - max_exp;
                *x_ptr = set_exponent(*x_ptr, expon);

                #if 0
                // skip integers
                if( *x_ptr == floor(*x_ptr) )
                    continue;
                #endif // 0

                // skip duplicates
                if( samples.count(*x_ptr) > 0 )
                    continue;

                samples.insert(*x_ptr);
            }
        }
    }
    std::cout << "total samples: " << samples.size() << std::endl;

    for(long double x_sample : samples) {
        long m=0, d=0;

        long double x_0 = 0.0;
        memset(&x_0, '\0', sizeof(x_0));
        x_0 = x_sample;
        long double x_i = x_0;

        // perform collatz operations to get m, d, and w.
        long double w = 0;
        std::vector<int> v;
        v.clear();
        while( gap(x_i) > 0 ) {
            // apply division-free recurrence
            long double lob = lob_value(x_i);
            x_i = 3*(x_i) + lob;
            w += lob;   // add lob to w

            // add lob position to v
            int lob_pos = round(log(lob)/log(2.0));
            v.push_back(lob_pos);
        }
        w += lob_value(x_i);    // add final lob (i.e. 2^d)
        d = 0;
        if( !v.empty() ) {
            d = v.back()-v[0];
        }
        m = v.size()-1;

        fout << std::setprecision(21)
            << w << " "
            << d << " "
            << m << " "
            << x_0 << " "
            << 1.0 << " "
            << std::endl;
    }
    mpz_clear(rnd_mpz);
    free(x_ptr); x_ptr=NULL;
}


int main() {

    std::fstream fout;
    fout.open("double_sampled.txt", std::fstream::out | std::fstream::trunc);
    sample_dataset_double(fout, 60, 76, 1000);
    fout.close();

    fout.open("double_enumerated.txt", std::fstream::out | std::fstream::trunc);
    enumerate_dataset_double(fout, 1e5, 1.0/32.0);
    fout.close();

    return 0;
}
