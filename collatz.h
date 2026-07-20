/*
 * collatz.h
 * NonintegerCollatz: Collatz for non-integers
 *
 * Copyright (c) 2025-2026 Bryan Franklin. All rights reserved.
 */
#ifndef COLLATZ_H
#define COLLATZ_H
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <vector>
#include <math.h>
#include <string.h>

//#define COLLATZ_TYPE unsigned long long int
#define COLLATZ_TYPE long long int
//#define COLLATZ_TYPE long double

static bool collatz_verbose = true;

COLLATZ_TYPE lobValue(COLLATZ_TYPE x) {
    return ((x | (x-1)) ^ x) + 1;
}
#define lob(x)  lobValue(x)

COLLATZ_TYPE hobValue(COLLATZ_TYPE x) {
    int i=0;

    while( x>1 ) { ++i; x/=2; }

    return (1<<i);
}
#define hob(x)  hobValue(x)

static inline COLLATZ_TYPE collatz(COLLATZ_TYPE x) {
    if( x < 1 ) return 0;   // too small

    if( (x%2) == 0 ) {
        return x/2;
    } else {
        if( x > (((COLLATZ_TYPE)1)<<((sizeof(COLLATZ_TYPE)*8)-2)) ) 
            return -1;  // too big
        return 3*x + 1;
    }
}


static inline int itersToOne(COLLATZ_TYPE x) {
    int ret = 0;
    while(x > 1 ) {
        x = collatz(x);
        ++ret;
    }
    return ret;
}


static inline int multsAndDivsToOne(COLLATZ_TYPE x, int &m, int &d) {
    int ret = 0;
    m = d = 0;
    while(x > 1 ) {
        if( (x%2) == 0 ) {
            ++d;
        } else {
            ++m;
        }
        x = collatz(x);
        ++ret;
    }
    return ret;
}


static inline COLLATZ_TYPE set_number(COLLATZ_TYPE ub, size_t offset, COLLATZ_TYPE lb) {
    return (ub<<offset)|lb;
}


static inline int getBit(COLLATZ_TYPE x, size_t k) {
    if( k < 0 || k > (8*sizeof(x)-1) ) return 0;
    return (int)((x>>k)&((COLLATZ_TYPE)1));
}


static inline std::string asBitString(COLLATZ_TYPE x_0, int bits=-1, int lowest_bit=0) {
    COLLATZ_TYPE x_j = x_0;

    if( x_0 == 0 ) return "0";
    if( bits<0 )    bits = 8*sizeof(x_0);

    if( lowest_bit > 0 && lowest_bit < (8*sizeof(x_j)) )
        x_j >>= lowest_bit;

    std::string ret = "";
    while( x_j > 0 && ret.length()<bits ) {
        COLLATZ_TYPE bit = x_j&1;
        ret = std::to_string(bit) + ret;
        x_j >>= 1;
    }
    while( ret.length()<bits ) { ret = "0" + ret; }

    return ret;
}


static inline std::string asBitStringLE(COLLATZ_TYPE x) {

    if( x == 0 ) return "0";

    std::string ret = "";
    while( x > 0 ) {
        COLLATZ_TYPE bit = x&1;
        ret = ret + std::to_string(bit);
        x >>= 1;
    }

    return ret;
}


std::string asShiftedBitString(COLLATZ_TYPE x, size_t zeros) {
    size_t max_cols = 130;
    std::string ret = asBitString(x);
    for(int i=0; i<zeros; ++i) {
        if( i==0 )
            ret += "o";
        else
            ret += "0";
    }
    for(int i=ret.length(); i<max_cols; ++i) {
        ret = " " + ret;
    }
    return ret;
}


static inline size_t countOnes(COLLATZ_TYPE x) {
    size_t ret = 0;
    while(x > 0 ) {
        x = (x & (x-1));
        ++ret;
    }
    return ret;
}


static inline int hobPos(COLLATZ_TYPE x) {

    int maxBit = sizeof(x) * 8 - 1;
    for(int i=maxBit-1; i>=0; --i) {
        if( getBit(x, i) != 0 ) {
            return i;
        }
    }

    return 0;
}


static inline int lobPos(COLLATZ_TYPE x) {

    int maxBit = sizeof(x) * 8 - 1;
    for(int i=0; i<maxBit; ++i) {
        if( getBit(x, i) != 0 ) {
            return i;
        }
    }

    return 0;
}


size_t gapSize(COLLATZ_TYPE x) {
    return hobPos(x)-lobPos(x);
}


size_t countZeros(COLLATZ_TYPE x) {
    size_t gap = gapSize(x);
    if( gap <= 1 )  return 0;   // 1 & 11 cases

    return countOnes(x) - 2;
}


double entropy(COLLATZ_TYPE x) {

    size_t ones = countOnes(x);  // include hob&lob ones
    size_t zeros = countZeros(x);
    double p_0 = -1.0;
    double p_1 = -1.0;
    if( ones > 0 || zeros > 0 ) {
        p_0 = zeros/(double)(ones+zeros);
        p_1 = ones/(double)(ones+zeros);
    }

    double sum = 0.0;
    if( p_0 > 0 )
        sum += zeros * p_0 * std::log2(p_0);
    if( p_1 > 0 )
        sum += ones * p_1 * std::log2(p_1);

    return -sum;
}


double internalEntropy(COLLATZ_TYPE x) {
    size_t hobPos_x = hobPos(x);
    size_t lobPos_x = lobPos(x);
    size_t gap = hobPos_x - lobPos_x;

    size_t bits = gap-2;  // ignore hob&lob bits
    size_t ones = countOnes(x)-2;  // ignore hob&lob ones
    size_t zeros = bits-ones;
    double p_0 = -1.0;
    double p_1 = -1.0;
    if( ones > 0 || zeros > 0 ) {
        p_0 = zeros/((double)bits);
        p_1 = ones/((double)bits);
    }
    std::cout << "p(b_i=0)=" << p_0 << "; p(b_i=1)=" << p_1 << std::endl;

    double sum = 0.0;
    if( zeros > 0 )
        sum += zeros * p_0 * std::log2(p_0);
    if( ones > 0 )
        sum += ones * p_1 * std::log2(p_1);

    return -sum;
}


COLLATZ_TYPE divFreeCollatz(COLLATZ_TYPE x) {
    return 3*x + (((COLLATZ_TYPE)1)<<lobPos(x));
}


static inline std::string bitSubstring(COLLATZ_TYPE x, size_t start, size_t len) {
    if( len > start+1 ) {
        return "Would extend beyond input.";
    }

    std::string ret = "";
    for(size_t j=start; j>(start-len); --j) {
        ret += ((x&(((COLLATZ_TYPE)1)<<j)) != 0 ) ? "1" : "0";
    }

    return ret;
}


static inline std::string topBits(COLLATZ_TYPE x, int numBits) {
    size_t i=hobPos(x);
    if( i < 0 )
        return "failed to find high order bit.";

    return bitSubstring(x, i, numBits);
}


static inline std::string bottomBits(COLLATZ_TYPE x, int numBits) {
    size_t i=lobPos(x);
    if( i < 0 )
        return "failed to find low order bit.";

    return bitSubstring(x, i+numBits-1, numBits);
}


static COLLATZ_TYPE find_x0_helper(COLLATZ_TYPE z, COLLATZ_TYPE x, int steps=-1) {
    if( collatz_verbose )
        std::cout << __FUNCTION__ << "(" << z << ", " << x << ", " << steps << ");" << std::endl;

    if( z < 0 ) {
        if( collatz_verbose ) {
            std::cout << "Negative z (" << z << ") not allowed." << std::endl;
        }
        return 0;
    }

    // check for base condition:
    // When z==0  the last bit has been subtracted from it.
    if( z==0 && x>=1 && steps<=0 ) {
        if( collatz_verbose ) {
            std::cout << "Found x=" << x << ", z==" << z << "." << std::endl;
        }
        return x;
    }

    if( steps == 0 ) {
        std::cout << "Ran out of steps." << std::endl;
        return 0;
    }

    // find values k such that (z-2^k)%3==0
    // report if multiple are found
    int hPos = hobPos(z);
    for(int i=hPos; i>=0; --i) {
        COLLATZ_TYPE ret = -1;
        COLLATZ_TYPE b_i = ((COLLATZ_TYPE)1)<<i;
        if( collatz_verbose ) {
            std::cout << "\ttrying b_i=" << b_i << "." << std::endl;
            std::cout << "\t\tz: " << z << std::endl;
            std::cout << "\t\tx: " << x << std::endl;
            std::cout << "\t\tz-b_i: " << (z-b_i) << std::endl;
            std::cout << "\t\t(z-b_i)%3: " << ((z-b_i)%3) << std::endl;
            std::cout << "\t\t(z-b_i)/3: " << ((z-b_i)/3.0) << std::endl;
        }

        if( b_i > z ) {
            if( collatz_verbose ) {
                std::cout << "b_i=" << b_i << " is too larger for z=" << z << "." << std::endl;
            }
            continue;
        }

        COLLATZ_TYPE zmb = z-b_i;
        COLLATZ_TYPE zmbm3 = zmb%3;
        if( zmbm3 != 0 ) {
            if( collatz_verbose ) {
                std::cout << "z-b_i=" << (zmb) << " is not divisible by 3, " << zmbm3 << "." << std::endl;
            }
            continue;
        }

        COLLATZ_TYPE zmbd3 = zmb/3;

        // make recursive call
        if( (ret=find_x0_helper( zmbd3, (x-b_i)/3, steps-1)) > 0 ) {
            if( collatz_verbose ) {
                std::cout << "b_i=" << b_i << " worked for " << ret << std::endl;
                std::cout << "\treturning " << ret << "." << std::endl;
            }

            // return on success
            return ret;
        }

        // report failure, if requested
        if( collatz_verbose ) {
            std::cout << "b_i=" << b_i << " failed (ret=" << ret << ")" << std::endl;
            std::cout << "2^" << i << " skipped." << std::endl;
        }
    }

    if( collatz_verbose ) {
        std::cout << __FUNCTION__ << ", ran out of bits to try.  Returning 0." << std::endl;
    }

    return 0;
}

// recursively find x from z
COLLATZ_TYPE find_x0(const COLLATZ_TYPE z, int steps=-1) {
    if( z <= 0 ) {
        if( collatz_verbose ) {
            std::cout << "z should be positive, not (" << z << ")." << std::endl;
        }
        return 0;
    }

    int hPos = hobPos(z);

    if( collatz_verbose ) {
        std::cout << "In " << __FUNCTION__ << "(" << z << "):" << std::endl;
        std::cout << "\thobPos(z): " << hPos << std::endl;
        std::cout << "\tsteps: " << steps << std::endl;
    }
    COLLATZ_TYPE ret = 0;
    for(int d=hPos+1; d>=0; --d) {
        COLLATZ_TYPE b = (((COLLATZ_TYPE)1)<<d);
        COLLATZ_TYPE y = z - b;
        if( collatz_verbose ) {
            std::cout << __FUNCTION__ << " trying:" << std::endl;
            std::cout << "\td: " << d << std::endl;
            std::cout << "\tb: " << b << std::endl;
            std::cout << "\ty: " << y << std::endl;
        }
        if( y < 0 ) {
            if( collatz_verbose ) {
                std::cout << "Skipping negative y." << std::endl;
            }
            continue;
        }
        if( (ret = find_x0_helper(y, b, steps)) > 0 ) {
            return ret;
        } else if( collatz_verbose ) {
            std::cout << __FUNCTION__ << " FAILed for b=" << b << " (ret="<< ret << ")." << std::endl;
        }
    }

    return 0;
}


// iteratively
COLLATZ_TYPE find_x0_iter(COLLATZ_TYPE z) {
    std::cout << __FUNCTION__ << "(" << z << ")" << std::endl;
    int hPos = hobPos(z);
    int lastBit = hPos;
    COLLATZ_TYPE x = 1;
    for(int i=hPos-1; i>=0; --i) {
        std::cout << "i=" << i << "; z=" << z << "; (1<<i)=" << (1<<i) << "; (z-(1<<i))=" << (z-(1<<i)) << " (0)" << std::endl;
        if( ((z-(1<<i))%3) == 0 && z>(1<<i) && z>0 ) {
            std::cout << "i=" << i << "; x=" << x << "; z=" << z << " (1)\t";
            int d = lastBit-i;
            std::cout << d << " division steps." << std::endl;
            x *= (1<<d);  // deal with division steps
            std::cout << "i=" << i << "; x=" << x << "; z=" << z << " (2)\t";
            if( z == 0 )
                return x;   // stop after final x/2 steps.
            if( ((x-1)%3) != 0 )
                std::cerr << "x-1 should be divisible by 3! (x-1=" << (x-1) << ")";
            x = (x-1)/3;            // deal with the 3x+1 step
            z = (z-(1<<i))/3;
            lastBit = i;
            std::cout << std::endl;
        }
    }
    if( z != 0 ) {
        std::cerr << "z=" << z << ", expcted to be zero." << std::endl;
    }

    x *= (1<<lastBit);
    std::cout << "x=" << x << std::endl;

    return x;
}

// second iterative attempt
COLLATZ_TYPE find_x0_iter2(COLLATZ_TYPE z) {
    std::cout << __FUNCTION__ << std::endl;
    if( z <= 0 ) {
        std::cout << "z (" << z << ") < 0, giving up." << std::endl;
        return 0;
    }
    COLLATZ_TYPE z_orig = z;

    // find, record, and remove hob to determine d.
    int lPos = lobPos(z);
    int hPos = hobPos(z);
    if( lPos == hPos ) {
        std::cout << "z (" << z << "), is powr of 2." << std::endl;
        return z;
    }
    int prevBit = hPos;
    COLLATZ_TYPE y = z - (((COLLATZ_TYPE)1)<<hPos);
    hPos = hobPos(y);

    // find encoded lob bits
    COLLATZ_TYPE x = 1;
    int mults = 0, divs = 0;
    for(int i=hPos-1; i>=0 && y>0; --i) {
        int j = prevBit-i;
        COLLATZ_TYPE b_i = ((COLLATZ_TYPE)1)<<i;

        std::cout << "y: " << y; // << std::endl;
        std::cout << "\tb_i: " << b_i; // << std::endl;
        std::cout << "\t(y-b_i): " << (y-b_i); // << std::endl;
        std::cout << "\t(y-b_i)%3: " << ((y-b_i)%3); // << std::endl;
        std::cout << "\t(((y-b_i)/3)%2): " << (((y-b_i)/3)%2) << std::endl;
        if( y < b_i ) continue;

        // update z and x
        if( b_i <= y
            && ((y-b_i)%3) == 0
            && (((y-b_i)/3)%2) == 1
            ) {

            divs += prevBit - i;
            std::cout << "divs: " << divs << std::endl;

            #if 1
            // update x
            std::cout << "D (by 2^" << j << ")" << std::endl;
            x <<= j;    // x *= 2^(j=prevBit-i)
            if( ((x-1)%3) != 0 )
                std::cout << "x should be divisible by 3! (x-1=" << (x-1) << ")" << std::endl;
            // 3x+1 step, in reverse, obviously.
            x -= 1;
            x /= 3;
            std::cout << "x: " << x << std::endl;
            #endif

            // update y
            y -= b_i;
            y /= 3;

            prevBit = i;
            ++mults;
        }
    }
    if( y != 0 ) {
        std::cout << "y: " << z << std::endl;
        ++mults;
    }

    // deal with any left-over bits
    divs += prevBit;
    if( prevBit > 0 )
        x *= ((COLLATZ_TYPE)1)<<prevBit;

    x = (pow(2, divs+1) - z_orig) / pow(3, mults);

    std::cout << __FUNCTION__ << ": x(z): " << x << "; mults: " << mults << "; divs: " << divs << "\t";

    return x;
}

std::vector<COLLATZ_TYPE> collatzPath(COLLATZ_TYPE x_0) {

    COLLATZ_TYPE x_i = x_0;
    std::vector<COLLATZ_TYPE> path;
    path.push_back(x_0);
    do {
        x_i = collatz(x_i);
        path.push_back(x_i);
    } while( x_i > 1 );

    return path;
}


std::string pathAsString(std::vector<COLLATZ_TYPE> &path) {
    std::string ret = "[ ";
    for(int i=0; i<path.size(); ++i) {
        if( i>0 )   ret += ", ";
        ret += std::to_string(path[i]);
    }
    ret += "]";
    return ret;
}


std::vector<COLLATZ_TYPE> commonPath(
    std::vector<COLLATZ_TYPE> path_a,
    std::vector<COLLATZ_TYPE> path_b) {
    std::vector<COLLATZ_TYPE> path_c;

    while( !path_a.empty()
        && !path_b.empty()
        && path_a.back() == path_b.back() ) {
        path_c.push_back(path_a.back());
        path_a.pop_back();
        path_b.pop_back();
    }

    return path_c;
}


COLLATZ_TYPE firstCommon(COLLATZ_TYPE a, COLLATZ_TYPE b) {

    std::vector<COLLATZ_TYPE> path_a = collatzPath(a);
    std::vector<COLLATZ_TYPE> path_b = collatzPath(b);

    if( collatz_verbose ) {
        std::cout << "path a: " << pathAsString(path_a) << std::endl;
        std::cout << "path b: " << pathAsString(path_b) << std::endl;
    }

    int pos_a = path_a.size()-1;
    int pos_b = path_b.size()-1;

    while( pos_a >= 0
            && pos_b >= 0 
            && path_a[pos_a] == path_b[pos_b] ) {
        if( collatz_verbose ) {
            std::cout << "\tpath_a[" << pos_a << "]: " << path_a[pos_a]
                << "; path_b[" << pos_b << "]: " << path_b[pos_b]
                << std::endl;
        }
        --pos_a;
        --pos_b;
    }

    if( collatz_verbose ) {
        std::cout << "pos_a: " << pos_a << " -> " << path_a[pos_a] << std::endl;
        std::cout << "pos_b: " << pos_b << " -> " << path_b[pos_b] << std::endl;
    }

    COLLATZ_TYPE ret = path_a[pos_a+1];
    if( ret != path_b[pos_b+1] )
        std::cerr << "Mismatch: path_a[" << (pos_a+1) << "]: " << ret
            << "; path_b[" << (pos_b+1) << "]: " << path_b[pos_b+1];

    if( collatz_verbose )
        std::cout << "ret: " << ret << std::endl;

    return ret;
}


typedef struct op_counts {
    size_t mults;
    size_t divs;
    double y_i;
} op_counts_t;

static inline int count_ops(COLLATZ_TYPE x_i, size_t &mults, size_t &divs, double &y_i) {
    static std::unordered_map<COLLATZ_TYPE, op_counts_t> cached;
    if( cached.count(x_i) > 0 ) {
        auto val = cached[x_i];
        mults = val.mults;
        divs = val.divs;
        y_i = val.y_i;
        return 0;
    }

    mults=0, divs=0;
    y_i = 0;

    while( x_i>1 ) {
        if( (x_i%2) == 0 ) {
            ++divs;
        } else {
            ++mults;
            y_i *= 3.0;
            y_i += std::pow(2.0, divs);
        }

        x_i = collatz(x_i);
    }

    cached[x_i] = {.mults=mults, .divs=divs, .y_i=y_i};

    return 0;
}


static inline size_t count_segments(COLLATZ_TYPE x_0, size_t min_gap = 2) {
    size_t ret = 0;

    COLLATZ_TYPE x_i = x_0;
    size_t zeros = 0;
    while( x_i > 1 ) {
        if( (x_i%2) == 0 ) {
            // zero bit
            ++zeros;
        } else {
            // one bit
            if( zeros>=min_gap ) {
                ++ret;
            }
            zeros = 0;
        }
        x_i /= 2;
    }

    return ret;
}


// a segment is based around the lob of he non-zero bits,
// and the zero bits above it.
// the segment containing the hob is defined to have 0 zeros.
// the 'offset' is the actual position of the lob for the segment.
typedef struct segment {
    size_t offset;
    COLLATZ_TYPE ones;
    size_t zeros;
} segment_t;


static inline std::vector<segment_t> get_segments(COLLATZ_TYPE x_0, size_t min_gap) {
    std::vector<segment_t> ret;

    COLLATZ_TYPE x_i = x_0;
    size_t hPos = hobPos(x_0);
    size_t lPos = lobPos(x_0);
    size_t zeros = 0;
    // first pass,
    size_t prev_one_pos = 0;
    for(size_t pos=lPos; pos<=hPos; ++pos) {
        int bit = getBit(x_i, pos);
        if( bit==0 ) {
            ++zeros;
        } else if (bit == 1 ) {
            if( zeros >= min_gap || ((pos-lPos) < min_gap && ret.size()==0) ) {
                //  find lob of each segment
                segment_t seg;
                memset(&seg, '\0', sizeof(seg));
                seg.offset = pos;   // record lob position for segment
                //  add entry to return vector
                ret.push_back(seg);
                if( ret.size() > 1 ) {
                    int prev = ret.size()-2;
                    // compute number of zeros since last non-zero.
                    ret[prev].zeros = pos-prev_one_pos-1;
                    // copy ones into segment
                    int num_ones = prev_one_pos-ret[prev].offset;
                    ret[prev].ones = ((1<<(num_ones+1))-1)
                                        & (x_i>>(ret[prev].offset));
                }
            }
            prev_one_pos = pos;
            zeros = 0;
        } else {
            std::cerr << "Unexpected bit value: " << bit << std::endl;
        }
    }

    // finish final segment
    if( zeros < min_gap ) { // zeros should be 0, as pos was hPos
        int last = ret.size()-1;
        // there are no zeros above the final segment
        ret[last].zeros = 0;
        // copy ones into segment
        ret[last].ones = x_i>>(ret[last].offset);
    }

    return ret;
}


static inline void print_segments(std::vector<segment_t> segs) {

    for(int i=0; i<segs.size(); ++i) {
        // find ones with 3 factors removed.
        COLLATZ_TYPE notThrees = segs[i].ones;
        int threeCount = 0;
        while( notThrees>0 && (notThrees%3) == 0 ) {
            notThrees /= 3;
            ++threeCount;
        }

        std::cout << "seg[" << i << "]: {"
            << segs[i].offset << ", "
            << asBitString(segs[i].ones) << ", "
            << segs[i].zeros << "} "
            << notThrees << " * 3^" << threeCount
            << std::endl;
    }
    std::cout << std::endl;
}


#endif // COLLATZ_H
