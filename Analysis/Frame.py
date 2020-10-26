import pdb
import numpy as np
import pandas as pd

from decorators import timer
from read_ascii_dat import read_dat_sylinder, read_dat_protein
from connected_components import get_nodes_in_clusters, get_edges_in_largest_cc
from calc_global_order import (calc_nematic_order, calc_polar_order,
                               calc_z_ordering)
from common_func import (calc_mean_pbc, calc_mean_separation,
                         calc_distance_pbc, unfold_coordinates)
from calc_tactoid_shape import calc_aspect_ratio
from calc_protein import calc_protein_energy

# A class to handle a single time frame in AMSOS


class Frame():
    def __init__(self, file_sylinder, file_protein, opts):

        self.file_sylinder = file_sylinder
        self.file_protein = file_protein
        self.opts = opts
        self.data = {}

    def read_frame(self):
        # get data frame for sylinder and protein files
        df_sylinder = read_dat_sylinder(self.file_sylinder)
        df_protein = read_dat_protein(self.file_protein)
        return df_sylinder, df_protein

    def analyze(self):

        # read frame
        df_sylinder, df_protein = self.read_frame()

        # Filament centers
        c = calc_mean_pbc(
            np.array(df_sylinder.pos0.tolist()),
            np.array(df_sylinder.pos1.tolist()),
            self.opts.boxsize)

        # Check existence of a cluster
        if self.opts.analyze_cluster:
            # Get largest connected component using information about xlinks
            cc, cc_bool = get_nodes_in_clusters(
                df_sylinder.gid.tolist(),
                df_protein.link0.tolist(),
                df_protein.link1.tolist(),
                min_size_ratio=0.1)
            df_sylinder_cc = df_sylinder[cc_bool]
            self.data['num_cluster'] = len(df_sylinder_cc)
        else:
            df_sylinder_cc = df_sylinder

        # Global order
        if self.opts.analyze_global_order:
            # Bulk and Cluster
            self.data['S_bulk'] = calc_nematic_order(
                np.array(df_sylinder.orientation.tolist()))
            self.data['P_bulk'] = calc_polar_order(
                np.array(df_sylinder.orientation.tolist()))
            self.data['S_cluster'] = calc_nematic_order(
                np.array(df_sylinder_cc.orientation.tolist()))
            self.data['P_cluster'] = calc_polar_order(
                np.array(df_sylinder_cc.orientation.tolist()))

        # Pair-pair separation
        if self.opts.analyze_pairpair_separation:

            # Filament centers
            c = calc_mean_pbc(
                np.array(df_sylinder.pos0.tolist()),
                np.array(df_sylinder.pos1.tolist()),
                self.opts.boxsize)
            self.data['mean_sep_bulk'] = list(
                calc_mean_separation(c, self.opts.boxsize))
            self.data['mean_sep_cluster'] = list(
                calc_mean_separation(c[cc, :], self.opts.boxsize))

        # Xlink
        if self.opts.analyze_xlinks:

            # Energy
            # xlink lengths
            xlink_lengths = np.linalg.norm(
                calc_distance_pbc(
                    np.array(df_protein.pos0.tolist()),
                    np.array(df_protein.pos1.tolist()),
                    self.opts.boxsize), axis=1)
            self.data['xlink_energy'] = list(
                calc_protein_energy(xlink_lengths, 0.05))

        if self.opts.analyze_aspect_ratio:
            self.data['tactoid_aspect_ratio'] = calc_aspect_ratio(
                unfold_coordinates(c[cc, :], self.opts.boxsize))

        if self.opts.analyze_z_ordering:
            self.data['z_order'] = calc_z_ordering(
                np.array(df_sylinder.orientation.tolist()))

        # if self.opts.analyze_local_order:
            # self.data['local_polar_order'] = calc_local_polar_order( np.array( df_sylinder.pos1.tolist() ),
            # self.opts.boxsize)

        # Length distribution inside vs outside cluster
        if self.opts.length_distribution:
            len_fils = calc_distance_pbc(
                np.array(df_sylinder.pos0.tolist()),
                np.array(df_sylinder.pos1.tolist()),
                self.opts.boxsize)
            self.data['length_mean_bulk'] = np.mean(len_fils)
            self.data['length_mean_cluster'] = np.mean(len_fils[cc_bool])
            self.data['length_mean_env'] = np.mean(
                len_fils[np.invert(cc_bool)])
