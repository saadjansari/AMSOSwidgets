/* This routine adds contributions to the density-density correlation function for objects
 in a 2d reference frame defined by axes parallel and perpendicular to the object director.
 Periodic boundary conditions are assumed.
 
 Input: number of objects (n_objects)
 unit cell matrix (h)
 array of scaled object coordinates (s_object)
 array of object coordinates (r_object)
 array of object directors (u_object)
 number of bins in parallel direction (n_bins_par)
 number of bins in perpendicular direction (n_bins_perp)
 bin width in parallel direction (bin_width_par)
 bin width in perpendicular direction (bin_width_perp)
 lower limit of correlation function range in parallel direction (dr_par_min)
 lower limit of correlation function range in perpendicular direction (dr_perp_min)
 total pair correlation function array (pair_dist_par_perp)
 pair correlation function array for parallel pairs (pair_dist_par_perp_parallel)
 pair correlation function array for antiparallel pairs (pair_dist_par_perp_antiparallel)
 polar orientational pair correlation function array (pair_corr_par_perp_P1)

 Output: the correlation function arrays are modified on return. */

#include "bob.h"

void pair_correlation_par_perp_FN(int n_objects, double **h, double **s_object, double **r_object, double **u_object, double **u_object_lat,
                                  int n_bins_par, int n_bins_perp, double bin_width_par, double bin_width_perp,
                                  double dr_par_min, double dr_perp_min,
                                  double **pair_dist_par_perp, double **pair_dist_par_perp_parallel, double **pair_dist_par_perp_antiparallel,
                                  double **pair_corr_par_perp_P1, double **pair_corr_par_perp_P1_lat)
{
    int i_object, j_object, i, i_bin_par, i_bin_perp;
    double dot_product, dr2, dr_par, dr_perp, dot_product_lat;
    static int first_call = 1;
    static double *dr;
    
    /* Allocate memory for local arrays the first time the routine is called. */
    if (first_call) {
        dr = (double*) allocate_1d_array(3, sizeof(double));
        first_call = 0;
    }
    
    /* Loop over pairs of objects. */
    for (i_object = 0; i_object < n_objects - 1; ++i_object)
        for (j_object = i_object + 1; j_object < n_objects; ++j_object) {
            
            /* Calculate dot product of object directors. */
            dot_product = 0.0;
            for (i = 0; i < 3; ++i)
                dot_product += u_object[i_object][i] * u_object[j_object][i];
       
            /* Calculate dot product of lateral directors. */
            dot_product_lat = 0.0;
            for (i = 0; i < 3; ++i)
                dot_product_lat += u_object_lat[i_object][i] * u_object_lat[j_object][i];

            /* Calculate inter-object separation vector. */
            pair_separation(3, 3, h, r_object[i_object], s_object[i_object], r_object[j_object], s_object[j_object], dr);
            
            /* Calculate squared magnitude of pair separation. */
            dr2 = 0.0;
            for (i = 0; i < 3; ++i)
                dr2 += SQR(dr[i]);
            
            /* Calculate components of separation vector parallel and perpendicular to director of first object. */
            dr_par = 0.0;
            for (i = 0; i < 3; ++i)
                dr_par += dr[i] * u_object[i_object][i];  // the separation vector goes from i_object -> j_object
            dr_perp = sqrt(dr2 - SQR(dr_par));
            
            /* Add contributions to correlation functions. */
            i_bin_par = floor((dr_par - dr_par_min) / bin_width_par);
            i_bin_perp = floor((dr_perp - dr_perp_min) / bin_width_perp);
            if (i_bin_par >=0 && i_bin_par < n_bins_par && i_bin_perp >=0 && i_bin_perp < n_bins_perp) {
                pair_dist_par_perp[i_bin_par][i_bin_perp] += 1.0;
                pair_corr_par_perp_P1[i_bin_par][i_bin_perp] += dot_product;
                if (dot_product > 0)
                    pair_dist_par_perp_parallel[i_bin_par][i_bin_perp] += 1.0;
                else
                    pair_dist_par_perp_antiparallel[i_bin_par][i_bin_perp] += 1.0;
            }
            
            /* Calculate components of separation vector parallel and perpendicular to director of second object. */
            dr_par = 0.0;
            for (i = 0; i < 3; ++i)
                dr_par += - dr[i] * u_object[j_object][i];  // the sign of the separation vector is reversed here, to go from j_object -> i_object
            dr_perp = sqrt(dr2 - SQR(dr_par));
            
            /* Add contributions to correlation functions. */
            i_bin_par = floor((dr_par - dr_par_min) / bin_width_par);
            i_bin_perp = floor((dr_perp - dr_perp_min) / bin_width_perp);
            if (i_bin_par >=0 && i_bin_par < n_bins_par && i_bin_perp >=0 && i_bin_perp < n_bins_perp) {
                pair_dist_par_perp[i_bin_par][i_bin_perp] += 1.0;
                pair_corr_par_perp_P1[i_bin_par][i_bin_perp] += dot_product;
                pair_corr_par_perp_P1_lat[i_bin_par][i_bin_perp] += dot_product_lat;
                if (dot_product > 0)
                    pair_dist_par_perp_parallel[i_bin_par][i_bin_perp] += 1.0;
                else
                    pair_dist_par_perp_antiparallel[i_bin_par][i_bin_perp] += 1.0;
            }
        }
    
    return;
}
