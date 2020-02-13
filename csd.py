#!/usr/bin/env python3

### simulations that track sex ratios 
### over several generations of inbreeding
### experiments
### when organism has multiple and 
### potentially linked CSD loci

### written by Bram Kuijper


# import some useful python libraries
import sys
import os
from numpy import *
from random import sample
import time
import copy
import pandas as pd


# the CSD class contains all the functions we need for a simulation
# of an inbreeding experiment, while varying the number of loci + survival 
class CSD:

    # specification of class-global variables
    pedigree_data = [] # variable to store all the date from the pedigree file

    # the name of an artificial pedigree that will be generated 
    # when simulation generates an artificial experiment
    artificial_pedigree_file_name = "artificial_pedigree200.csv" 

    # variable indicating if we 
    # are indicating that we are running an artifical experiment
    artificial_experiment = True

    # differential survival of a diploid male
    diploid_male_survival = 1

    # the number of loci assumed to underly sex phenotype
    number_loci = 3 

    # the file where the outcome of the model is written to
    # the filename is now generated based on the parameters in the
    # __init()__ function below
    statsfile_name = ""

    # same, but then for family-level data
    famstatsfile_name = ""

    # the amount of linkage between the different loci
    # one could make it more complicated by also varying the linkage between
    # different pairs of loci, but this is not so informative (I think)
    linkage = 0.5

    replicate_runs = 50000# the required number of replicate simulations
    replicate = 0 # the current replicate
    recombination_vector = [] # the vector with recombination distances between loci
    mothers_P = [] # a vector filled with mommies from the parental generation

    # variable that counts the number of families
    # which produce no daughters at all
    count_no_daughters = 0
    
    # the number of generations that a simulation runs
    number_generations = 2 

    # the current generation that is running
    generation = 1 

    # the number of artificial lines; only important
    # when running in 'generator' mode
    number_artificial_lines = 10 
                                    
    # number to seed the random number generator
    random_seed = time.time(); 

    # specify the column names of the pedigree file
    ped_cn_generation = "generation"
    ped_cn_female_code = "female_code"
    ped_cn_diploids = "no_females"
    ped_cn_offspring_used = "amount_used_offspring"
    ped_cn_mother_code = "mother_id"
    ped_cn_haploids = "no_males"


    # constructor of the class (i.e., with which one
    # initializes things)
    # two arguments:
    #   artificial
    #       - simulation in 'generator' mode when artificial=True
    #            meaning that it generates an artificial pedigree
    #       - simulation in mornal mode when artificial=False
    #           it then accepts a given pedigree file and runs
    #           10000 replicate runs
    #   generations
    #       works only in 'generator' mode and determines
    #       the amount of generations for which a pedigree
    #       needs to be created
    def __init__(self
            ,pedigree_file_name
            ,pedigree_file_delimiter=","
            ,artificial=False
            ,generations=1
            ,nr_loci=1
            ,survival=1
            ,linkage=0.5):

        # set initializer arguments
        self.pedigree_file_name = pedigree_file_name
        self.artificial_experiment = artificial
        self.number_generations = generations
        self.number_loci = nr_loci
        self.diploid_male_survival = survival
        self.linkage = linkage
        self.csv_delimiter=pedigree_file_delimiter

        # generate filename
        self.statsfile_name = "csd_run_" \
            + str(self.number_loci) \
            + "loci_surv" \
            + str(self.diploid_male_survival) \
            + "_link_" + str(self.linkage)

        # generate the filename to store family-level data
        self.famstatsfile_name = self.statsfile_name + "fam"

        # get the pedigree data if it is not in 'generator' 
        # mode
        if not self.artificial_experiment:
            print("getting the experimental data")
            self.get_data()

        # initialize the recombination vector
        # specifying the linkage between different loci
        self.initialize_recomb_vec()
       
        # in case of a 'generator' run,
        # generate the artificial pedigree
        if self.artificial_experiment:
            self.generate_experiment()

    # function to read pedigree data from the pedigree file
    # this can be rather flexible,
    # as long as the data is structured as follows
    def get_data(self):

        # read the data
        pedigree_data_raw = pd.read_csv(
                filepath_or_buffer=self.pedigree_file_name
                ,sep=self.csv_delimiter

        # generate list of the breeders for each generation
        self.Nbreeders = [ 0 for i in range(0,self.number_generations) ]

        # initialize the data object that is going to store
        # our pedigree data
        self.pedigree_data = []

        # copy the pedigree data into that data object
        for row in dictreadert:

            # if the current generation in the current row
            # is larger than zero, copy this row of the dataset
            # to the pedigree dataset
            # as generation 0 is general outcrossing, rather than a mother-son
            # mating
            if int(row[self.ped_cn_generation]) > 0:
                self.pedigree_data.append(copy.deepcopy(row))
                
                if len(self.Nbreeders)  < int(row[self.ped_cn_generation]) + 1:
                    self.Nbreeders.append(0)

                # indicate that breeders will need to be produced
                # in this generation
                self.Nbreeders[int(row[self.ped_cn_generation])]+=1

                # also find the maximum number of generations 
                # that are needed by the simulation
                if int(row[self.ped_cn_generation]) > self.number_generations:
                    self.number_generations = int(row[self.ped_cn_generation])


    # more specific function for the rubicula data
    def get_data_rubicula(self):

        # get the csv file
        csvfile = open(self.pedigree_file_name,"r")
        
        # read the pedigree data from the file
        dictreadert = csv.DictReader(csvfile,delimiter=self.csv_delim)

        # initialize the data object that is going to store
        # our pedigree data
        self.pedigree_data = []

        # copy the pedigree data into that data object
        for row in dictreadert:

            the_row = copy.deepcopy(row)

            generation = 1 
            
            if the_row["Family.mother"] != "na":
                generation = 2

            if generation == 2:
                mother_id = the_row["Family.mother"]

                # okay offspring from a certain mother,
                # increment that mother's number of offspring
                # that is used for further breeding experiments
                for any_row in self.pedigree_data:
                    if any_row["female_code"] == mother_id:
                        any_row["amount_used_offspring"]+=1

            # get the data row in an appropriate format
            data_row = dict(
                    female_code=the_row["Replicate"], 
                    mother_id=the_row["Family.mother"], 
                    haploid_offspring=int(the_row["Haploidmale"]), 
                    diploid_offspring=int(the_row["Diploidmale"]) + float(the_row["Female"]), 
                    amount_used_offspring=0, 
                    generation=generation)

            # add the current row of data to the pedigree
            self.pedigree_data.append(data_row)


    # initialize the vector of recombination distances
    def initialize_recomb_vec(self):

        self.recombination_vector = [] 
    
        # fill the recombination vector with linkage values
        for i in range(0,self.number_loci):
            self.recombination_vector.append(self.linkage)

    # the main function used by a simulation 
    # in 'generator' mode
    # generate an artificial experiment
    # running over $generations
    def generate_experiment(self):
       
        # open the stats file to write sex ratio
        # data to
        self.statsfile = open(self.statsfile_name,"w")
        self.write_headers_stats()

        # loop through n generations
        self.generation = 1

        # for generation 1, initialize generation P
        self.artificial_P_generation()

        # and initialize this generation
        mothers = self.init_P_generation()

        # loop through the subsequent generations to set out
        # the pedigree further
        for self.generation in range(2, self.number_generations+1):
            
            # every generation generates new mothers
            # for the next generation
            mothers = self.run_single_generation(mothers)

        # after having run, close off the statistics file
        # and the pedigree data
        self.statsfile.close()
        self.write_artificial_pedigree()


    # in 'generator mode'
    # generates an artificial parental generation
    def artificial_P_generation(self):

        self.pedigree_data = []

        # generate some artificial mothers
        for i in range(1, self.number_artificial_lines+1):

            # own individual's id?
            # TODO

            # generate a mother's id for use in the pedigree
            mum_id = "A_" + str(i)

            # generate the 'contents' of this mother
            self.create_artificial_mother(mum_id)

    # when in 'generator' mode: create a mother 
    # that will have a certain number of offspring
    # of which a certain amount of daughters will be used for further
    # crosses 
    # note that this function does not yet do anything with 
    # genotypes; this will come later.
    def create_artificial_mother(self,the_id, mother_id):
       
        # get the number of diploid offspring generated by this mother
        # from a negative binomial distribution with a shape that 
        # approximates the offspring distribution from de Boer et al 2008
        number_diploids = random.negative_binomial(1, p=0.1, size=1.3)[0]

        # the number of daughters used for further crosses
        amount_used_offspring = 0 

        # we just need one individual per line
        if number_diploids > 0:
            amount_used_offspring = 1

        # add a row to the pedigree data with the data of this individual
        data_row = dict(
                female_code=the_id, 
                mother_id=mother_id, 
                diploid_offspring=number_diploids, 
                amount_used_offspring=amount_used_offspring, 
                generation=self.generation
                )
        self.pedigree_data.append(data_row)

        return data_row

    # find the parent with id parent_code
    # and add the number noffspring to its amount of used offspring
    def add_used_offspring(self,parent_code, noffspring):

        for row in self.pedigree_data:

            if self.ped_cn_female_code in row:
                if row[self.ped_cn_female_code] == parent_code:
                    row[self.ped_cn_offspring_used] += noffspring
                    break

    # wenjuan's data lacks a pedigree,
    # meaning that we have to simulate that in our model
    def generate_pedigree_wenjuan(self):

        code = 0
        current_generation = 1

        # define a list to store the parental ids
        # those will subsequently be used to assign to 
        # the offspring generation
        parental_codes = []

        # go through generation 1
        # and generate parental ids
        for row in self.pedigree_data:
            if int(row[self.ped_cn_generation]) == 1:
                row[self.ped_cn_female_code] = str(code)
                row[self.ped_cn_offspring_used] = 0
                code+=1
                parental_codes.append(str(code))

        # go through all the next generations
        # get the parental ids from the previous generations
        # and generate offspring ids from them
        for generation_i in range(2,self.number_generations):

            # there may be more offspring then there are currently parental codes
            # i.e., some parents may have given rise to multiple offspring breeders
            # add additional parental codes for that
            n_more_than_parents = self.Nbreeders[generation_i]-len(parental_codes)

            # add extra parents to the list of parental id's
            if n_more_than_parents > 0:
                parental_codes_new = parental_codes
                for i in range(0,n_more_than_parents):
                    random_parent_code = parental_codes[random.randint(len(parental_codes))]
                    parental_codes_new.append(random_parent_code)
            
                random.shuffle(parental_codes_new)

            elif n_more_than_parents < 0:
                random.shuffle(parental_codes)
                parental_codes = parental_codes[0:self.Nbreeders[generation_i]]

            parental_codes_tplus1 = []

            # now go through the list of parental codes
            # count the occurrence of each parental code
            # to assing the number of used offspring to each parent
            for code in list(set(parental_codes)):
                count = parental_codes.count(code)
                if count > 0:
                    self.add_used_offspring(code, count)       


            # loop through pedigree data
            for row in self.pedigree_data:

                # select the rows that fit the current generation
                if int(row[self.ped_cn_generation]) == generation_i:

                    # important, reset the number of offspring that 
                    # this parent generates
                    # to zero
                    row[self.ped_cn_offspring_used] = 0

                    if len(parental_codes) == 0:
                        print("parental codes exhausted")
                        sys.exit(1)
                    
                    # add a parental id to the current offspring
                    code_parent = parental_codes.pop()
                    code_current_off = code_parent + "_" + str(parental_codes.count(code_parent) + 1)
                    row[self.ped_cn_female_code] = code_current_off
                    row[self.ped_cn_mother_code] = code_parent

                    parental_codes_tplus1.append(code_current_off)

            parental_codes = parental_codes_tplus1

    # only used in 'generator' mode
    # write all the data from the artificially generated pedigree
    # to the pedigree file
    def write_artificial_pedigree(self):
        
        # open an artificial pedigree file 
        # and write the data-headers to it
        # it will be filled later
        self.artificial_pedigree_file = open(self.artificial_pedigree_file_name,"w")
        self.artificial_pedigree_file.write("female_code;diploid_offspring;amount_used_offspring;generation\n")

        # write a row of data to the file
        for row in self.pedigree_data:
            self.artificial_pedigree_file.write(row[self.ped_cn_female_code] + ";" + str(row[self.ped_cn_diploids]) + ";" + str(row[self.ped_cn_offspring_used]) + ";" + str(row[self.ped_cn_generation]) + "\n")

        # close the file after writing
        self.artificial_pedigree_file.close()

    # the main function when used in 'normal' 
    # mode
    # runs one instance of a simulation
    # inbreeding experiment
    def simulate(self):

        # open the statsfile to write the aggregrate data to
        # (e.g.,, overall sex ratio, total number of families, etc")
        self.statsfile = open(self.statsfile_name,"w")
        self.write_headers_stats()

        # open the statsfile to write the per-family data to
        # (e.g., number of diploids / family, number of diploid sons, etc)
        self.famstatsfile = open(self.famstatsfile_name,"w")

        # write headers to the family statistics file
        self.famstatsfile.write("replicate;" 
            + "generation;"
            + "female_id;"
            + "ndiploids;"
            + "ndiploid_sons;"
            + "prop_diploid_sons\n")
        
        # loop through all the replicate runs
        for self.replicate in range(0, self.replicate_runs):
            print("replicate run " 
                    + str(self.replicate) 
                    + " out of " 
                    + str(self.replicate_runs))

            self.generate_pedigree_wenjuan()

            # initialize the P generation females
            mothers = self.init_P_generation()

            # if there are no P generation mums,
            # we have a bug.
            if mothers == []:
                print("failing miserably!")
                sys.exit(1)

            # loop through all the generations
            for self.generation in range(1, self.number_generations):

                # each generation produces new mothers for the subsequent one 
                mothers = self.run_single_generation(mothers)

        # close output files
        self.statsfile.close()
        self.famstatsfile.close()

    # for a number of daughters, run brother sister matings
    # arguments:
    #       mother: a dictionary containing all the mother's data
    #       daughter_list: a list with daughters from that mother
    def brother_sister_matings(self, mother, daughter_list):

        # variable to store the new mothers after the 
        # the daughter brother mating
        new_mothers = []

        mother_id = mother[self.ped_cn_female_code]

        # get the list of possible daughters from this mother
        daughter_lookup_list = self.daughters_from_pedigree(mother_id)

        daughter_counter = 0

        # for each daughter, get her data from the file
        # and mate her with her brother
        for daughter in daughter_list:

            # generate a brother genome to fertilize her
            son_genome = self.get_maternal_gamete(mother["genomes"])

            # get the daughter's data in the data set
            daughter_data = daughter_lookup_list[daughter_counter]

            # if there is no daughter data...
            if daughter_data == None:

                # in 'generator mode' create a daughter
                if self.artificial_experiment:
                    daughter_data = self.create_artificial_mother(daughter_id, mother["id"]) 
                # if this is 'normal mode', something is wrong. All daughters that are used are 
                # present in the pedigree
                else:
                    print(
                            "cannot find daughter with id " 
                            + daughter_id 
                            + " for usage in daughter-son "
                            + "mating in generation " 
                            + str(self.generation)
                            )
                    sys.exit(1)

            # add the genomes of the daughter and her brother to the dataset
            # belonging the daughter
            daughter_data["genomes"] =  daughter["genomes"]
            daughter_data["genome_partner"] = son_genome

            # append this mated daughter to the list of daughters
            new_mothers.append(daughter_data)

            daughter_counter += 1

        return new_mothers

    # allow mothers to produce offspring
    def produce_offspring(self, mother):

        # make the diploid offspring belonging to this mother
        offspring_data = self.make_diploid_offspring(mother)

        assert(offspring_data["Noffspring"] > 0)

        self.famstatsfile.write(str(self.replicate) + ";"\
                + str(self.generation) + ";"\
                + mother["female_code"] + ";"\
                + str(offspring_data["Noffspring"]) + ";"\
                + str(offspring_data["Nsons"]) + ";"\
                + str(float(offspring_data["Nsons"])/float(offspring_data["Noffspring"])) + "\n")

        # get the number of diploids produced by this mother
        self.Ndiploids += offspring_data["Noffspring"]

        self.Nfamilies += 1

        if offspring_data["Noffspring"] > 0:
            self.Nmothers_with_diploids+=1
            if offspring_data["Nsons"] == offspring_data["Noffspring"]:
                self.count_no_daughters+=1

        self.Nsons+=offspring_data["Nsons"]

        if offspring_data["Noffspring"] > 0:
            self.prop_diploid_males += float(offspring_data["Nsons"]) / float(offspring_data["Noffspring"]) 


        if int(mother[self.ped_cn_haploids]) + offspring_data["Noffspring"] > 0:
            self.prop_sons_total += float(int(mother[self.ped_cn_haploids]) + offspring_data["Nsons"]) / float(int(mother[self.ped_cn_haploids]) + offspring_data["Noffspring"])

        self.all_daughters[mother[self.ped_cn_female_code]] = offspring_data["list_daughters"]

    # write the statistics for this generation
    def write_stats_this_generation(self):

        # take account of the sex ratio
        sr_this_generation = self.prop_diploid_males / float(self.Nfamilies)
        prop_sons_this_generation = self.prop_sons_total / float(self.Nfamilies)

        # write this in the statistics file
        self.statsfile.write(str(sr_this_generation) + ";"\
        + str(prop_sons_this_generation) + ";"\
        + str(self.Nsons) + ";"\
        + str(self.Ndiploids) + ";"\
        + str(self.Nfamilies) + ";"\
        + str(self.generation) + ";"\
        + str(self.replicate) + ";"\
        + str(self.number_loci) + "\n")


    # daughters from a particular mother are mated
    def mate_daughters(self, mother):

        # get the amount of offspring from this particular mother
        # that are used for further crosses
        # from the pedigree
        used_nr_offspring = int(mother[self.ped_cn_offspring_used])

        # get the daughters produced by the mother in the
        # previous round
        mother_id = mother[self.ped_cn_female_code]

        # if the number of daughters produced by this mother
        # is smaller than the number of offspring used for crosses
        # we have to reduce the latter number
        # 
        # this could occur because the number of daughters produced
        # is a stochastic variable, sampled from the distribution
        # of offspring numbers 
        if used_nr_offspring > len(self.all_daughters[mother_id]):
            used_nr_offspring = len(self.all_daughters[mother_id])

        # get a list of the required daughter genomes
        daughters_required = self.all_daughters[mother_id][0:used_nr_offspring]

        # get a list of daughters that we do not use now
        remaining_daughters = self.all_daughters[mother_id][used_nr_offspring:len(self.all_daughters[mother_id])]

        # put the unused daughters back on the stack
        if len(remaining_daughters) > 0:
            self.all_daughters[mother_id] = remaining_daughters
        else:
            del self.all_daughters[mother_id] # otherwise remove this mother's share

        # in case of artificial experiment:
        #
        # if the current mother has not produced any daughters, 
        # this line is extinct. It has to be replaced by a line from
        # another mother that has produced more than
        # 1 daughter
        if used_nr_offspring == 0 and self.artificial_experiment:

            lines_extinct += 1
            # first, reset the value of used offspring in the artificial
            # pedigree
            for row in self.pedigree_data:
                # check if id's match
                if row[self.ped_cn_female_code] == mother[self.ped_cn_female_code]:
                    row[self.ped_cn_offspring_used] = 0
                    break;

        # return mated daughters for this mother
        return self.brother_sister_matings(mother, daughters_required)

    # function will only be used in case we are simulating
    # an actual breeding experiment
    # since we want to make most optimal use of available resources
    # as soon as line goes extinct, it needs to be replaced by 
    # another one. This function does that.
    def replace_empty_positions(self):
        # condition below is only fulfilled in case of 
        # an artificial experiment
        # 
        # here, breeding positions are filled up again
        # that were empty due to line extinction
        # in case mothers did not reproduce (i.e., no offspring produced)
        # or only diploid sons
        while lines_extinct > 0:

            current_set_of_spots = lines_extinct

            moms_with_kids_left = self.all_daughters.keys()

            if len(moms_with_kids_left) < current_set_of_spots:
                current_set_of_spots = len(moms_with_kids_left)

            # not any more mothers
            if current_set_of_spots == 0:
                break

            random_moms = sample(moms_with_kids_left, current_set_of_spots)

            for mom in random_moms:
                current_daughter = self.all_daughters[mom][0:1]
                del self.all_daughters[mom][0:1]

                if len(self.all_daughters[mom]) == 0:
                    del self.all_daughters[mom]

                mother = self.female_from_pedigree(mom)
                
                for row in self.pedigree_data:
                    # check if id's match
                    if row[self.ped_cn_female_code] == mom:
                        row[self.ped_cn_offspring_used] = row[self.ped_cn_offspring_used]+1
                        break;

                brosismating = self.brother_sister_matings(mother, current_daughter)
                new_mothers.append('')
                new_mothers[-1:] = brosismating
                lines_extinct -=1
        
    # run a single generation
    # this is the workhorse of both the 'normal'
    # and the 'generator' model
    def run_single_generation(self, mothers):

        # during a single generation
        # - mothers make offspring
        # - write the stats with the number of diploids
        # - let daughters mate
        # - return the mated daughters as new mothers
        #   for the next generation

        # the average sex ratio is only taken over mothers
        # that actually do produce diploids
        # hence, count the number of diploids
        self.Nmothers_with_diploids = 0
        self.Ndiploids = 0
        self.Nsons = 0
        self.prop_diploid_males = 0
        self.Nfamilies = 0
        self.prop_sons_total = 0

        new_mothers = []
        self.all_daughters = {}

        # in case of artificial experiment,
        # variable that counts if there are lines extinct
        # and hence have to be replaced by lines from other
        # mothers
        lines_extinct = 0

        # first allow all mothers to produce offspring

        for mother in mothers:
            self.produce_offspring(mother)

        # write down the statistics (sex ratio etc)
        self.write_stats_this_generation()
          
        # end simulation if necessary
        if self.generation == self.number_generations:
            return

        
        # then allow all the produced daughters to mate
        #
        # in case of an artificial cross:
        #   when a mother did not produce any usable offspring
        #   pick second, or third offspring from other mothers
        for mother in mothers:
            new_mothers.append('')
            # add mothers to the stack
            new_mothers[-1:] = self.mate_daughters(mother)

        while lines_extinct > 0:
            self.replace_empty_positions()
#
                
        # allright, daughters have replaced previous generation
        # end of this generation
        return new_mothers

    def make_diploid_offspring(self, mother):

        daughters = []

        nr_sons = 0

        i = 0
        
        ndip = int(float(mother[self.ped_cn_diploids]))

        while i < ndip:

            gamete_mom = self.get_maternal_gamete(mother["genomes"])
            gamete_dad = mother["genome_partner"]

            if gamete_mom == gamete_dad:
                if random.uniform() < self.diploid_male_survival:
                    nr_sons +=1
                    i+=1
            else:
                daughters.append([gamete_mom,gamete_dad])
                i+=1

        # the daughters are then randomly shuffled 
        random.shuffle(daughters)

        # counter
        counter = 1
        newdt = []

        # then give every daughter an ID
        for dtr in daughters: 
            newdt.append(dict(id="na", mother_id = mother[self.ped_cn_female_code], genomes=dtr))
            counter+=1

        daughters = newdt

        nr_offspring = nr_sons + len(daughters)

        return(dict(Noffspring=nr_offspring,Nsons=nr_sons,list_daughters=daughters))

    # lookup a female based on her id in the dataset
    # return her data
    def female_from_pedigree(self, female_id):
        
        for data_row in self.pedigree_data:
            if data_row[self.ped_cn_female_code] == female_id:
                return data_row


    def daughters_from_pedigree(self, mother_id):

        daughters = []

        for data_row in self.pedigree_data:
            if self.ped_cn_mother_code in data_row and data_row[self.ped_cn_mother_code] == mother_id:
                daughters.append(data_row)

        # as a last check see if the number of used daughters 
        # matches with 
        mother = self.female_from_pedigree(mother_id)

        if mother[self.ped_cn_offspring_used] !=  len(daughters):

            print("error: number daughters found for mother with id " 
                    + mother_id 
                    + " is " 
                    + str(len(daughters)) 
                    + ", but the used amount of daughters is " 
                    + str(mother[self.ped_cn_offspring_used])
                    )
            sys.exit(1)


        return daughters

    # get all data of mothers from a particular generation
    def females_from_pedigree(self, generation):

        selection = []

        # loop through the data and get P parents
        for data_row in self.pedigree_data:

            # select the mothers from generation 1
            if int(data_row[self.ped_cn_generation]) == generation:
                selection.append(data_row)

        return(selection)

    # make a gamete out of mother's genomes
    def get_maternal_gamete(self, list_genomes):
       
        # random number indicating which genome copy
        # to start with
        strand = random.random_integers(0,1)

        # if there is only one locus,
        # just return a random haplotype
        if self.number_loci == 1:
            return(list_genomes[strand])

        haplotype = ""

        for locus in range(0, self.number_loci):

            haplotype += list_genomes[strand][locus] 

            # anything but the last locus may show recombination
            # which means that the strands flip
            if (locus < self.number_loci - 1 and random.uniform() < self.recombination_vector[locus]):
                strand = abs(strand-1) # flip between 0 and 1

        return haplotype

    # initialize the parental generation
    # allowing for a mother-son mating
    def init_P_generation(self):

        generation = 1

        # get the selection of the females that are used this
        # generation, from the pedigree
        current_females = self.females_from_pedigree(generation)

        if current_females == []:
            print("cannot initialize number of females in P generation")
            sys.exit(1)

        # loop through all mothers and give them
        # random fully heterozygous strings
        for P_mother in current_females:
            mother_genomes = self.generate_P_mother()
            P_mother["genomes"] = mother_genomes
            P_mother["genome_partner"] = self.get_maternal_gamete(mother_genomes) # make a son and make this the mum's partner

        return current_females

    def generate_P_mother(self):

        genome1 = ""
        genome2 = ""

        for i in range(0,self.number_loci):
            # first allele is 0 or 1
            allele1 = random.random_integers(0,1)
            # the other allele is thus always the opposite
            allele2 = abs(allele1-1)

            genome1+=str(allele1)
            genome2+=str(allele2)

        return([genome1,genome2])

    def write_headers_stats(self):
        self.statsfile.write(
                "prop_diploid_males;"
                + "prop_males;"
                + "nsons;"
                + "ndiploids;"
                + "nfamilies;"
                + "generation;"
                + "replicate;"
                + "number_loci\n")

