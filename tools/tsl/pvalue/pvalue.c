/**
 *  C routine for computing p-values.  
 *
 *  Version 1.2
 *
 *  Vladimir Vacic, University of California, Riverside
 *  Lilia M. Iakoucheva, Rockefeller University, New York
 *  Predrag Radivojac, Indiana University, Bloomington
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "specfns.h"


/**
 *  
 */
double binom(int x, int N, double p)  {
    double ln_y, coef;

    coef = lgam(N + 1) - lgam(x + 1) - lgam(N - x + 1);
    ln_y = coef + x * log(p) + (N - x) * log(1 - p);
    return exp(ln_y);
}


/**
 *  Estimates p-value using the binomial test.
 */
double bin_p_value(int k1, int n1, int k2, int n2)  {
    double p, t, s, sum;
    double J[n2];
    int    i, j;

    p = 1.0 * (k1+k2) / (n1+n2);
    t = fabs(1.0*k1/n1 - 1.0*k2/n2);
    
    for(j=0; j<=n2; j++)
        J[j] = binom(j, n2, p);

    sum = 0;
    for(i=0; i<n1; i++)  {
        s = 0;
        for(j=0; j<=n2; j++)  
            if (fabs(1.0*i/n1 - 1.0*j/n2) >= t)
                s += J[j];
        
        sum += s * binom(i, n1, p);
    }

    return sum;
}


/** Estimates p-value using two sample t-test. */
double ttest_p_value(int k1, int n1, int k2, int n2)  {
    double mean1, mean2;
    double var1_mult, var2_mult;
    double svar, df, t;

    mean1 = 1.0 * k1 / n1;
    mean2 = 1.0 * k2 / n2;

    var1_mult = 1.0*k1*(1-mean1)*(1-mean1) + (n1-k1)*mean1*mean1;
    var2_mult = 1.0*k2*(1-mean2)*(1-mean2) + (n2-k2)*mean2*mean2;

    df   = n1 + n2 - 2;
    svar = (var1_mult + var2_mult) / df;
    t    = (mean1-mean2) / sqrt(svar*(1.0/n1 + 1.0/n2));

    return incbet(0.5*df, 0.5, df/(df+t*t));
}


void usage()  {
    printf("Usage: pvalue b|t COUNT_1_+ COUNT_1_- COUNTS_2_+ COUNT_2_- ...\n\n");

    printf("Calculates p-values for a group of symbol occurence counts in\n"); 
    printf("(+) positive/functional and (-) negative/control samples, under\n");
    printf("the null hypothesis that a symbol occurence shows no difference\n");
    printf("between the two samples.\n\n");

    printf("    b|t    Choice of a statistical test. 'b' for binomial test\n");
    printf("           or 't' for Student's t-test.\n\n");
}


int main(int argc, char** argv)  {
    int i, num_symbols;
    int *pos_counts, *neg_counts;
    int pos_total, neg_total;

    // Basic parameter checking
    if (1 == argc)  {
        usage(); exit(1);
    }
    if ('b'!=argv[1][0] && 't'!=argv[1][0])  {
        printf("ERROR: Statistical test parameter should be either 'b' or 't',\n");
        printf("but is '%s'.\n\n", argv[1]);
        usage(); exit(1);
    }
    if (argc < 6)  {
        printf("ERROR: At least 4 symbol counts should be provided: counts in\n");
        printf("positive and negative samples for at least two symbols.\n\n");
        usage(); exit(1);
    } 

    num_symbols = argc/2-1;

    pos_counts = (int*) calloc(num_symbols, sizeof(int));
    neg_counts = (int*) calloc(num_symbols, sizeof(int));

    pos_total = 0;
    neg_total = 0;

    for(i=0; i<num_symbols; i++)  {
        pos_counts[i] = atoi(argv[2*i+2]);
        neg_counts[i] = atoi(argv[2*i+3]);

        pos_total += pos_counts[i];
        neg_total += neg_counts[i];
    }

    if (0==pos_total)  {
        printf("ERROR: Total number of symbols in positives is 0.\n"); 
        usage(); exit(1);
    }

    if (0==neg_total)  {
        printf("ERROR: Total number of symbols in negatives is 0.\n"); 
        usage(); exit(1);
    }

    if ('b'==argv[1][0])  {
        for(i=0; i<num_symbols; i++)
            printf("%f\n", bin_p_value(pos_counts[i], pos_total, neg_counts[i], neg_total));
    }

    if ('t'==argv[1][0])  {
        for(i=0; i<num_symbols; i++)
            printf("%f\n", ttest_p_value(pos_counts[i], pos_total, neg_counts[i], neg_total));
    }
    
    free(pos_counts);
    free(neg_counts);

    return 0;
}

