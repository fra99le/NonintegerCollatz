/*
 * collatz_w.cpp
 * NonintegerCollatz: Collatz for non-integers
 *
 * Copyright (c) 2026 Bryan Franklin. All rights reserved.
 */
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <unordered_set>
#include <vector>
#include "collatz.h"
#include "collatz_gmp.h"

// clear; c++ -I/opt/local/include -L/opt/local/lib -lgmp -lgmpxx collatz_w.cpp -o collatz_w && time ./collatz_w 33 && cat divfree_table.tex && time ./plot_w.gplot && time ./plot_dbl.gplot

typedef struct bit_info {
    int r, g, b, a;
} bit_info_t;


int write_divfree_table(const std::string filename, const std::string png, const long x_0) {

    if( (x_0%2) == 0 ) {
        std::cerr << "x_0 must be odd.\n";
        exit(1);
    }
    mpz_t x_i;
    mpz_t pow2, w;
    mpz_init_set_ui(x_i, x_0);
    mpz_init(pow2);
    mpz_init(w);
    int lobPos0 = lob_pos(x_i);
    std::cout << "lobPos0: " << lobPos0 << std::endl;

    // do initial run to gather basic stats
    int m = 0, d = 0;
    mpz_set_ui(w, 0);
    while( mpz_cmp_ui(x_i, 1) > 0 ) {
        if( mpz_odd_p(x_i) ) {
            mpz_mul_ui(x_i, x_i, 3);
            mpz_add_ui(x_i, x_i, 1);
            ++m;

            mpz_ui_pow_ui(pow2, 2, d);
            mpz_add(w, w, pow2);
        } else {
            mpz_divexact_ui(x_i, x_i, 2);
            ++d;
        }
    }
    mpz_ui_pow_ui(pow2, 2, d);
    mpz_add(w, w, pow2);
    std::cout << "x_0: " << x_0 << "; w: " << w << std::endl;
    std::cout << "m: " << m << "; d: " << d << std::endl;
    std::cout << "k: " << (m/(double)d) << std::endl;

    // build column format string
    std::string columnsStr = "|r|";
    std::string firstRow = "";
    firstRow += " & \\multicolumn{" + std::to_string(m+1) + "}{|c|}{Arithmetic Region ($AR_i$)} & & & \\\\\n";
    firstRow += "\\hline\n";
    firstRow += "$i$";
    int column = 0;
    for(int i=d+1; i>=0; --i) {
        if( mpz_tstbit(w, i) ) {
            firstRow += " & " + std::to_string(m-column);
            columnsStr += "r|";
            ++column;
        }
    }
    firstRow += " & $x_i$ & $v_i$ & $\\textit{lob}(x_i)$ \\\\";
    columnsStr += "r|r|r|r|";

    // create image to buffer bit histories
    bit_info_t bits[d+1][m+1];
    memset(&bits, '\0', sizeof(bits));
    mpz_init_set_ui(x_i, x_0);
    bit_info_t zero = {255, 255, 255, 255};
    bit_info_t one = {0, 0, 0, 255};
    bit_info_t hobBit = {255, 0, 0, 255};
    bit_info_t lobBit = {127, 0, 0, 255};

    mpz_set_ui(x_i, x_0);
    int iter = 0;
    int offset = 0;
    while( mpz_even_p(x_i) ) {
        mpz_divexact_ui(x_i, x_i, 2);
        bits[offset][0] = zero;
        ++offset;
    }
    std::vector<int> arWidths;
    int width = 0;
    while( true ) {

        // copy bits into image
        int hobPos = hob_pos(x_i);
        int lobPos = lob_pos(x_i);
        std::cout << "hob,lob: " << hobPos << ", " << lobPos << std::endl;
        for(int j=0; j<=hobPos; ++j) {
            int bit = mpz_tstbit(x_i, j);
            if( bit )
                bits[offset+j][iter] = one;
            else
                bits[offset+j][iter] = zero;
        }
        // mark hob with special values
        bits[offset+hobPos][iter] = hobBit;
        bits[offset+lobPos][iter] = lobBit;

        if( mpz_cmp_ui(x_i, 1) == 0 )   break;

        // apply 3x+1 rule
        mpz_mul_ui(x_i, x_i, 3);
        mpz_add_ui(x_i, x_i, 1);

        // handle 'division'
        while( mpz_even_p(x_i) ) {
            mpz_divexact_ui(x_i, x_i, 2);
            ++offset;
            ++width;
        }
        arWidths.push_back(width);
        width=0;

        ++iter;
    }

    // open file
    std::ofstream ofs;
    ofs.open(filename, std::ofstream::out | std::ofstream::trunc);
    std::cout << "Writing " << filename << std::endl;


    ofs << "\\begin{tabular}{" << columnsStr << "}" << std::endl;

    ofs << firstRow << std::endl;
    ofs << "\\hline" << std::endl;

    bit_info_t bit;
    std::string zeroStr = " ";
    int lobPos = -1;
    mpz_t x_val, w_val;
    mpz_init_set_ui(x_val, x_0);
    mpz_init_set_ui(w_val, 0);
    for(int j=0; j<=m; ++j) {
        zeroStr = " ";
        ofs << j << " & ";
        lobPos = (d+1);
        for(int i=d; i>=0; --i) {
            bit = bits[i][j];
            if( bit.a != 0 ) {
                if( bit.b != 0 ) {
                    ofs << "0";
                } else {
                    ofs << "1";
                    lobPos = i;
                }
            } else {
                ofs << zeroStr;
            }

            if( mpz_tstbit(w, i) ) {
                ofs << " & ";
                if( bit.a != 0 ) {
                    zeroStr = "0";
                } else {
                    zeroStr = " ";
                }
            }
        }
        ofs << x_val;
        mpz_ui_pow_ui(pow2, 2, lobPos);
        mpz_add(w_val, w_val, pow2);
        mpz_mul_ui(x_val, x_val, 3);    // 3*x_i
        mpz_add(x_val, x_val, pow2);    // add lob(x_i);
        ofs << " & " << lobPos;
        ofs << " & " << " $2^{" << lobPos << "}$ ";
        ofs << " \\\\" << std::endl;
    }
    ofs << "\\hline" << std::endl;

    // output w, alinged with x
    ofs << " w  & ";
    for(int i=d; i>=0; --i) {
        if( mpz_tstbit(w, i) )
            ofs << "1 & ";
        else
            ofs << "0";
    }
    ofs << " & & " << w_val << " \\\\" << std::endl;

    ofs << "\\end{tabular}" << std::endl;

    // close file
    ofs.close();

    // convert w back to x
    mpz_t w_gmp;
    mpq_t xq;
    mpq_init(xq);
    mpz_init_set(w_gmp, w_val);
    w_to_xq(w_gmp, xq);
    std::cout << "xq: " << xq << std::endl;

    mpz_clear(x_i);
    mpz_clear(pow2);
    mpz_clear(w_val);
    mpq_clear(xq);

    return 0;
}


int write_divfree_table_html(const std::string filename, const std::string png, const long x_0) {

    if( (x_0%2) == 0 ) {
        std::cerr << "x_0 must be odd.\n";
        exit(1);
    }
    mpz_t x_i;
    mpz_t pow2, w;
    mpz_init_set_ui(x_i, x_0);
    mpz_init(pow2);
    mpz_init(w);
    int lobPos0 = lob_pos(x_i);
    std::cout << "lobPos0: " << lobPos0 << std::endl;

    // do initial run to gather basic stats
    int m = 0, d = 0;
    mpz_set_ui(w, 0);
    while( mpz_cmp_ui(x_i, 1) > 0 ) {
        if( mpz_odd_p(x_i) ) {
            mpz_mul_ui(x_i, x_i, 3);
            mpz_add_ui(x_i, x_i, 1);
            ++m;

            mpz_ui_pow_ui(pow2, 2, d);
            mpz_add(w, w, pow2);
        } else {
            mpz_divexact_ui(x_i, x_i, 2);
            ++d;
        }
    }
    mpz_ui_pow_ui(pow2, 2, d);
    mpz_add(w, w, pow2);
    std::cout << "x_0: " << x_0 << "; w: " << w << std::endl;
    std::cout << "m: " << m << "; d: " << d << std::endl;
    std::cout << "k: " << (m/(double)d) << std::endl;

    // build column format string
    std::string firstRow = "";
    firstRow += " <td style=\"padding: 0px;\"> </td> <td style=\"padding: 0px;\" colspan=\"" + std::to_string(m+1) + "\">Arithmetic Region (AR<sub>i</sub>) </td><td style=\"padding: 0px;\"> &nbsp; </td><td style=\"padding: 0px;\"> &nbsp; </td><td style=\"padding: 0px;\"> &nbsp; </td>\n</tr>\n<tr style=\"border-bottom: 1px solid #000;\">\n";
    firstRow += "<td style=\"padding: 0px; text-align: left\"> i";
    int column = 0;
    for(int i=d+1; i>=0; --i) {
        if( mpz_tstbit(w, i) ) {
            firstRow += " </td><td style=\"padding: 0px; text-align: right;\"> " + std::to_string(m-column);
            ++column;
        }
    }
    firstRow += " </td><td style=\"padding: 0px; text-align: right;\"> x<sub>i</sub> </td>";
    firstRow += "<td style=\"padding: 0px; text-align: right;\"> v<sub>i</sub> </td>";
    firstRow += "<td style=\"padding: 0px; text-align: right;\"> lob(x<sub>i</sub>) </td>\n</tr>\n";


    // create image to buffer bit histories
    bit_info_t bits[d+1][m+1];
    memset(&bits, '\0', sizeof(bits));
    mpz_init_set_ui(x_i, x_0);
    bit_info_t zero = {255, 255, 255, 255};
    bit_info_t one = {0, 0, 0, 255};
    bit_info_t hobBit = {255, 0, 0, 255};
    bit_info_t lobBit = {127, 0, 0, 255};

    mpz_set_ui(x_i, x_0);
    int iter = 0;
    int offset = 0;
    while( mpz_even_p(x_i) ) {
        mpz_divexact_ui(x_i, x_i, 2);
        bits[offset][0] = zero;
        ++offset;
    }
    std::vector<int> arWidths;
    int width = 0;
    while( true ) {

        // copy bits into image
        int hobPos = hob_pos(x_i);
        int lobPos = lob_pos(x_i);
        std::cout << "hob,lob: " << hobPos << ", " << lobPos << std::endl;
        for(int j=0; j<=hobPos; ++j) {
            int bit = mpz_tstbit(x_i, j);
            if( bit )
                bits[offset+j][iter] = one;
            else
                bits[offset+j][iter] = zero;
        }
        // mark hob with special color
        bits[offset+hobPos][iter] = hobBit;
        bits[offset+lobPos][iter] = lobBit;

        if( mpz_cmp_ui(x_i, 1) == 0 )   break;

        // apply 3x+1 rule
        mpz_mul_ui(x_i, x_i, 3);
        mpz_add_ui(x_i, x_i, 1);

        // handle 'division'
        while( mpz_even_p(x_i) ) {
            mpz_divexact_ui(x_i, x_i, 2);
            ++offset;
            ++width;
        }
        arWidths.push_back(width);
        width=0;

        ++iter;
    }

    // open file
    std::ofstream ofs;
    ofs.open(filename, std::ofstream::out | std::ofstream::trunc);
    std::cout << "Writing " << filename << std::endl;


    ofs << "<table style=\"border-spacing: 1px;\"><tbody>" << std::endl;

    ofs << firstRow << std::endl;

    bit_info_t bit;
    std::string zeroStr = " ";
    int lobPos = -1;
    mpz_t x_val, w_val;
    mpz_init_set_ui(x_val, x_0);
    mpz_init_set_ui(w_val, 0);
    for(int j=0; j<=m; ++j) {
        zeroStr = " ";
        ofs << "<tr><td style=\"padding: 0px; text-align: left;\">" << j << " </td><td style=\"padding: 0px; text-align: right\"> ";
        lobPos = (d+1);
        for(int i=d; i>=0; --i) {
            bit = bits[i][j];
            if( bit.a != 0 ) {
                if( bit.b != 0 ) {
                    ofs << "0";
                } else {
                    ofs << "1";
                    lobPos = i;
                }
            } else {
                ofs << zeroStr;
            }

            if( mpz_tstbit(w, i) ) {
                ofs << " </td><td style=\"padding: 0px; text-align: right;\"> ";
                if( bit.a != 0 ) {
                    zeroStr = "0";
                } else {
                    zeroStr = " ";
                }
            }
        }
        ofs << x_val;
        mpz_ui_pow_ui(pow2, 2, lobPos);
        mpz_add(w_val, w_val, pow2);
        mpz_mul_ui(x_val, x_val, 3);    // 3*x_i
        mpz_add(x_val, x_val, pow2);    // add lob(x_i);
        ofs << " </td><td style=\"padding: 0px; text-align: right;\"> " << lobPos;
        ofs << " </td><td style=\"padding: 0px; text-align: right;\"> " << " 2<sup>" << lobPos << "</sup> </td> ";
        ofs << "</tr>" << std::endl;
        ofs << "<tr>" << std::endl;
    }

    // output w, alinged with x
    ofs << "<td style=\"padding: 0px; text-align: left;\">w</td><td style=\"padding: 0px; text-align: right\">";
    for(int i=d; i>=0; --i) {
        if( mpz_tstbit(w, i) )
            ofs << "1</td><td style=\"padding: 0px; text-align: right;\">";
        else
            ofs << "0";
    }
    ofs << "</td><td style=\"padding: 0px;\"> &nbsp; </td><td style=\"padding: 0px; text-align: right\">" << w_val << " </tr>" << std::endl;

    ofs << "</tbody></table>" << std::endl;

    // close file
    ofs.close();

    // convert w back to x
    mpz_t w_gmp;
    mpq_t xq;
    mpq_init(xq);
    mpz_init_set(w_gmp, w_val);
    w_to_xq(w_gmp, xq);
    std::cout << "xq: " << xq << std::endl;

    mpz_clear(x_i);
    mpz_clear(pow2);
    mpz_clear(w_val);
    mpq_clear(xq);

    return 0;
}


void enumerate_dataset_w(std::ostream &out, mpz_t &max_w, long long step_w) {
    mpz_t w_gmp, m, d;
    mpz_init(w_gmp);
    mpz_init(m);
    mpz_init(d);
    mpq_t xq;
    mpq_init(xq);
    std::cout << "max_w: " << max_w << "; step_w: " << step_w << std::endl;

    // sample regularly with step_w step size
    unsigned long int long_max_w = mpz_get_ui(max_w);
    for(mpz_set_ui(w_gmp, 1);
        mpz_cmp(w_gmp, max_w) < 0;
        mpz_add_ui(w_gmp, w_gmp, step_w)) {
        std::cout << "." << w_gmp << std::flush;
        w_to_xq(w_gmp, xq);

        w_to_m(w_gmp, m);
        w_to_d(w_gmp, d);

        if( mpz_cmp_ui(mpq_denref(xq), 1) != 0 ) {
            out << w_gmp << " "
                << d << " "
                << m << " "
                << mpq_numref(xq) << " "
                << mpq_denref(xq) << " "
                << std::endl;
        }
    }
    std::cout << std::endl << "done." << std::endl;

    mpz_clear(w_gmp);
    mpz_clear(m);
    mpz_clear(d);
    mpq_clear(xq);
}


void sample_dataset_w(std::ostream &out, const long max_bits, const long long n_samples) {
    mpz_t w_gmp, m, d;
    mpz_init(w_gmp);
    mpz_init(m);
    mpz_init(d);
    mpq_t xq;
    mpq_init(xq);
    // ramdomly sample for n_samples samples */
    // see: https://stackoverflow.com/a/11738909
    mpz_class value;
    gmp_randclass r(gmp_randinit_default);
    std::cout << value << std::endl;

    mpz_t rand_w;
    mpz_init(rand_w);
    std::set<std::string> samples;
    for(long n_bits=16; n_bits<=max_bits; ++n_bits) {
        samples.clear();
        while( samples.size() < n_samples 
                && samples.size() < 0.5*(1<<n_bits) ) {
            mpz_class w_class;
            w_class = r.get_z_bits(n_bits);

            size_t buf_size = mpz_sizeinbase(w_class.get_mpz_t(), 36) + 2;
            char *w_str = (char*)calloc(buf_size, sizeof(char));
            mpz_get_str(w_str, 36, w_class.get_mpz_t());
            std::cout << "w_str: " << w_str << std::endl;
            samples.insert(w_str);
        }
        for(std::string w_str : samples) {
            mpz_set_str(w_gmp, w_str.c_str(), 36);
            w_to_xq(w_gmp, xq);
            std::cout << "random w: " << w_gmp
                        << " (" << n_bits << " bits)" << std::endl;
            std::cout << "w -> x: " << xq << std::endl;

            w_to_m(w_gmp, m);
            w_to_d(w_gmp, d);
            std::cout << "m,d: " << m << ", " << d << std::endl;

            out << w_gmp << " "
                << d << " "
                << m << " "
                << mpq_numref(xq) << " "
                << mpq_denref(xq) << " "
                << std::endl;
        }
    }
    mpz_clear(w_gmp);
    mpz_clear(m);
    mpz_clear(d);
    mpz_clear(rand_w);
    mpq_clear(xq);
    std::cout << std::endl << "done." << std::endl;
}


void enumerate_dataset_x(std::ostream &out, int max_x, long step_size) {
    mpz_t w_gmp, m, d;
    mpz_t x_gmp;
    mpz_init(w_gmp);
    mpz_init(m);
    mpz_init(d);
    mpz_init(x_gmp);
    for(unsigned long int x=1; x<max_x; x+=step_size) {
        mpz_set_ui(x_gmp, x);
        x_to_w(x_gmp, w_gmp);
        
        w_to_m(w_gmp, m);
        w_to_d(w_gmp, d);

       std::cout << "x->w: " << x_gmp << " -> " << w_gmp << " (max_w: " << max_x << ")" << std::endl;
       std::cout << "m,d: " << m << ", " << d << std::endl;

        out << w_gmp << " "
            << d << " "
            << m << " "
            << x_gmp << " "
            << 1 << " "
            << std::endl;
    }

    mpz_clear(w_gmp);
    mpz_clear(m);
    mpz_clear(d);
    mpz_clear(x_gmp);
}


void sample_dataset_x(std::ostream &out, const int max_bits, const long n_samples) {
    mpz_class x_class;
    gmp_randclass r(gmp_randinit_default);
    mpz_t x_gmp;
    mpz_t w_gmp, m, d;
    mpz_init(w_gmp);
    mpz_init(m);
    mpz_init(d);
    mpz_init(x_gmp);
    std::cout << "up to " << max_bits << " bits, with " << n_samples << " per bits range." << std::endl;
    for(unsigned long int bits=16; bits<max_bits; ++bits) {
        long samples = n_samples;
        while( --samples > 0 ) {
            std::cout << bits << " / " << max_bits << " bits.  " << std::flush;
            std::cout << samples << " samples left.  " << std::flush;
            x_class = r.get_z_bits(bits);
            mpz_set(x_gmp, x_class.get_mpz_t());
            x_to_w(x_gmp, w_gmp);

            w_to_m(w_gmp, m);
            w_to_d(w_gmp, d);
            std::cout << "m,d: " << m << ", " << d << std::endl;

            std::cout << "x->w: " << x_gmp << " -> " << w_gmp << std::endl;

            out << w_gmp << " "
                << d << " "
                << m << " "
                << x_gmp << " "
                << 1 << " "
                << std::endl;
        }
    }

    mpz_clear(w_gmp);
    mpz_clear(m);
    mpz_clear(d);
    mpz_clear(x_gmp);
}


int main(int argc, char **argv) {
    #if 0
    if( argc < 1 ) {
        std::cerr << argv[0] << " example_x" << std::endl;
        return 1;
    }
    #endif // 0

    long long argv1 = atol(argv[1]);
    mpz_t gmp_arg;
    mpz_init(gmp_arg);
    write_divfree_table("divfree_table.tex", "divfree_table.png", argv1);
    write_divfree_table_html("MathML/01_divfree_table.html", "MathML/divfree_table.png", argv1);

    std::fstream ds;
    ds.open("divfree_w_sampled.txt", std::ofstream::out | std::ofstream::trunc);
    std::cout << "Writing sampling of w." << std::endl;
    sample_dataset_w(ds, 670, 500);
    ds.close();

    ds.open("divfree_x_sampled.txt", std::ofstream::out | std::ofstream::trunc);
    std::cout << "Writing sampling of x." << std::endl;
    sample_dataset_x(ds, 60, 10000);
    ds.close();

    std::cout << "Writing enumeration of x. " << std::endl;
    ds.open("divfree_x_enumerated.txt", std::ofstream::out | std::ofstream::trunc);
    enumerate_dataset_x(ds, 100000, 1);
    ds.close();

    std::cout << "Writing enumeration of w. " << std::endl;
    ds.open("divfree_w_enumerated.txt", std::ofstream::out | std::ofstream::trunc);
    mpz_init_set_ui(gmp_arg, 1000000);
    enumerate_dataset_w(ds, gmp_arg, 1);
    ds.close();

    mpz_clear(gmp_arg);

    return 0;
}
