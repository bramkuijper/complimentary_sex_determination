#!/usr/bin/env python3

# run the actual simulation

import csd

# initialize the simulation
a = csd.CSD(
        pedigree_file_name="Lhet_adapted.csv"
        ,artificial=False # it is not an artificial experiment we try to generate
        ,generations=4 # maximum number of generations the experiment lasts
        ,nr_loci=1 # number of csd loci
        ,survival=1 # survival probability of diploid males
        ,linkage=0.5) # when arg_nr_loci > 1, linkage between the different csd loci [ 0: full linkage, 0.5 no linkage

# simulate the actual thing
a.simulate()

