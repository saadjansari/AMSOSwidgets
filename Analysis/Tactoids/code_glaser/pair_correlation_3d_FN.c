/* This routine adds contributions to the density-density correlation function between specific objects (atoms or bonds) in a local reference frame.
 
 Input: number of molecules (n_mols)
 number of atoms per molecule (n_atoms_per_mol)
 unit cell matrix (h)
 inverse unit cell matrix (h_inv)
 array of scaled atom coordinates (s_atom)
 array of atom coordinates (r_atom)
 array of molecular directors (u_mol)
 intramolecular labels of atoms defining local reference frame (ref_label_1, ref_label_2, ref_label_3, ref_label_4)
 first object type (object_type_1): 0 = atom, 1 = bond
 intramolecular atom labels defining first object (atom_label_1, atom_label_2)
 second object type (object_type_2): 0 = atom, 1 = bond
 intramolecular atom labels defining second object (atom_label_3, atom_label_4)
 array of dimensions of correlation function (n_bins)
 array of bin widths (bin_width)
 array of correlation function half ranges (half_range)
 total pair correlation function array (pair_dist_local)
 parallel pair correlation function array (pair_dist_local_parallel)
 antiparallel pair correlation function array (pair_dist_local_antiparallel)
 polar orientational pair correlation function array (pair_corr_local_P1)

 Output: the correlation function arrays are modified on return. */

#include "bob.h"

void pair_correlation_3d_FN(int n_mols, int n_atoms_per_mol, double **h, double **h_inv, double **s_atom, double **r_atom, double **u_mol,
                            int ref_label_1, int ref_label_2, int ref_label_3, int ref_label_4,
                            int object_type_1, int atom_label_1, int atom_label_2,
                            int object_type_2, int atom_label_3, int atom_label_4,
                            int *n_bins, double *bin_width, double *half_range,
                            double *pair_dist_local, double *pair_dist_local_parallel, double *pair_dist_local_antiparallel, double *pair_corr_local_P1)
{
    int i_mol, j_mol, ref_atom_1, ref_atom_2, ref_atom_3, ref_atom_4, atom_1, atom_2, atom_3, atom_4,
    i, j, k, i_linear, out_of_range;
    double pi, dot_product, object_dot_product, length, length2, mag, mag2;
    static int first_call = 1, *i_bin;
    static double **e_local, **v_local,
    *v_bond_1, *u_bond_1, *r_bond_1, *s_bond_1,
    *v_bond_2, *u_bond_2, *r_bond_2, *s_bond_2,
    *r_object_1, *s_object_1, *r_object_2, *s_object_2,
    *dr, *dr_local;
    
    /* Allocate memory for local arrays the first time the routine is called. */
    if (first_call) {
        i_bin = (int*) allocate_1d_array(3, sizeof(int));
        v_local = (double**) allocate_2d_array(3, 3, sizeof(double));
        e_local = (double**) allocate_2d_array(3, 3, sizeof(double));
        v_bond_1 = (double*) allocate_1d_array(3, sizeof(double));
        u_bond_1 = (double*) allocate_1d_array(3, sizeof(double));
        r_bond_1 = (double*) allocate_1d_array(3, sizeof(double));
        s_bond_1 = (double*) allocate_1d_array(3, sizeof(double));
        v_bond_2 = (double*) allocate_1d_array(3, sizeof(double));
        u_bond_2 = (double*) allocate_1d_array(3, sizeof(double));
        r_bond_2 = (double*) allocate_1d_array(3, sizeof(double));
        s_bond_2 = (double*) allocate_1d_array(3, sizeof(double));
        r_object_1 = (double*) allocate_1d_array(3, sizeof(double));
        s_object_1 = (double*) allocate_1d_array(3, sizeof(double));
        r_object_2 = (double*) allocate_1d_array(3, sizeof(double));
        s_object_2 = (double*) allocate_1d_array(3, sizeof(double));
        dr = (double*) allocate_1d_array(3, sizeof(double));
        dr_local = (double*) allocate_1d_array(3, sizeof(double));
        first_call = 0;
    }
    
    /* Compute constants. */
    pi = acos(-1.0);
    
    /* Loop over first molecule. */
    for (i_mol = 0; i_mol < n_mols; ++i_mol) {
        
        /* Get labels of atoms that define molecule-fixed reference frame. */
        ref_atom_1 = i_mol * n_atoms_per_mol + ref_label_1;  // the first two reference atoms define the z axis
        ref_atom_2 = i_mol * n_atoms_per_mol + ref_label_2;
        ref_atom_3 = i_mol * n_atoms_per_mol + ref_label_3;  // the second two reference atoms define the x axis
        ref_atom_4 = i_mol * n_atoms_per_mol + ref_label_4;
        
        /* Calculate unit vectors that define molecule-fixed reference frame. */
        bond_vector(3, 3, h, s_atom, r_atom, ref_atom_1, ref_atom_2,  // third basis vector is parallel to vector joining first pair of reference atoms
                    v_local[2], e_local[2], &length, &length2);
        bond_vector(3, 3, h, s_atom, r_atom, ref_atom_3, ref_atom_4,  // first vector joins second pair of reference atoms (must have component perpendicular to first basis vector)
                    v_local[0], e_local[0], &length, &length2);
        mag2 = 0.0;  // calculate and normalize cross product of two vectors (second basis vector)
        for (i = 0; i < 3; ++i) {
            j = (i + 1) % 3;
            k = (i + 2) % 3;
            v_local[1][i] = e_local[2][j] * e_local[0][k] - e_local[2][k] * e_local[0][j];
            mag2 += SQR(v_local[1][i]);
        }
        mag = sqrt(mag2);
        for (i = 0; i < 3; ++i)
            e_local[1][i] = v_local[1][i] / mag;
        for (i = 0; i < 3; ++i) {  // the cross product of the second and third basis vectors gives the first basis vector
            j = (i + 1) % 3;
            k = (i + 2) % 3;
            e_local[0][i] = e_local[1][j] * e_local[2][k] - e_local[1][k] * e_local[2][j];
        }
        
        /* Get labels of target atoms in first molecule. */
        atom_1 = i_mol * n_atoms_per_mol + atom_label_1;
        atom_2 = i_mol * n_atoms_per_mol + atom_label_2;
        
        /* Compute bond vector joining two target atoms in first molecule, and compute position of bond. */
        bond_vector(3, 3, h, s_atom, r_atom, atom_1, atom_2,
                    v_bond_1, u_bond_1, &length, &length2);
        bond_position(3, 3, h, h_inv, r_atom, atom_1,
                      v_bond_1, r_bond_1, s_bond_1);
        
        /* Determine position of object (bond or atom) in first molecule. */
        if (object_type_1 == 0) {  // if object_type_1 == 0, use position of first target atom
            for (i = 0; i < 3; ++i)
                s_object_1[i] = s_atom[atom_1][i];
        }
        else {  // if object_type_1 == 1, use position of bond
            for (i = 0; i < 3; ++i)
                s_object_1[i] = s_bond_1[i];
        }
        
        /* Loop over second molecule. */
        for (j_mol = 0; j_mol < n_mols; ++j_mol) {
            
            if (j_mol != i_mol) {  // if j_mol != i_mol, add contribution to intermolecular correlation function
                
                /* Get labels of target atoms in second molecule. */
                atom_3 = j_mol * n_atoms_per_mol + atom_label_3;
                atom_4 = j_mol * n_atoms_per_mol + atom_label_4;
                
                /* Compute bond vector joining two target atoms in second molecule, and compute position of bond. */
                bond_vector(3, 3, h, s_atom, r_atom, atom_3, atom_4,
                            v_bond_2, u_bond_2, &length, &length2);
                bond_position(3, 3, h, h_inv, r_atom, atom_3,
                              v_bond_2, r_bond_2, s_bond_2);
                
                /* Determine position of object (bond or atom) in second molecule. */
                if (object_type_2 == 0) {  // if object_type_2 == 0, use position of first target atom
                    for (i = 0; i < 3; ++i)
                        s_object_2[i] = s_atom[atom_3][i];
                }
                else {  // if object_type_2 == 1, use position of bond
                    for (i = 0; i < 3; ++i)
                        s_object_2[i] = s_bond_2[i];
                }
                
                /* Calculate inter-object separation vector in global reference frame. */
                pair_separation(3, 3, h, r_object_1, s_object_1, r_object_2, s_object_2, dr);
                
                /* Transform inter-object separation vector into local reference frame. */
                transform_vector(3, dr, e_local, dr_local);
                
                /* Calculate dot product of molecular directors. */
                dot_product = 0.0;
                for (i = 0; i < 3; ++i)
                    dot_product += u_mol[i_mol][i] * u_mol[j_mol][i];
                
                /* Calculate dot product of object directors. */
                object_dot_product = 0.0;
                for (i = 0; i < 3; ++i)
                    object_dot_product += u_bond_1[i] * u_bond_2[i];

                /* Add contributions to pair correlation functions. */
                out_of_range = 0;
                for (i = 0; i < 3; ++i) {
                    i_bin[i] = floor((dr_local[i] + half_range[i]) / bin_width[i]);
                    out_of_range = i_bin[i] < 0 || i_bin[i] >= n_bins[i];
                    if (out_of_range) break;
                }
                if (!out_of_range) {
                    i_linear = linear_index_row_major(2, 3, i_bin, n_bins);
                    pair_dist_local[i_linear] += 1.0;
                    pair_corr_local_P1[i_linear] += object_dot_product;
                    if (dot_product > 0.0)
                        pair_dist_local_parallel[i_linear] += 1.0;
                    else
                        pair_dist_local_antiparallel[i_linear] += 1.0;
                }
            }
        }
    }
    
    return;
}
