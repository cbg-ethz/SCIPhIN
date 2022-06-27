/**
 * SCIPhIN: Single-cell mutation identification via phylogenetic inference
 * <p>
 * Copyright (C) 2022 ETH Zurich, Jochen Singer
 * <p>
 * This file is part of SCIPhI.
 * <p>
 * SCIPhIN is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p>
 * SCIPhIN is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p>
 * You should have received a copy of the GNU General Public License
 * along with SCIPhI. If not, see <http://www.gnu.org/licenses/>.
 *
 * @author: Jochen Singer
 */
#define _SECURE_SCL 1

#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include <sstream>
#include <time.h>
#include <tuple>
#include <cmath>
#include <assert.h>

#include "sciphi_config.h"
#include "attachmentScores.h"
#include "mcmc.h"
#include "readData.h"
#include "output.h"
#include "version.h"

#include <boost/multi_array.hpp>
#include <boost/program_options.hpp>

template<typename TTreeType>
std::vector<double> setMoveProbs(Config<TTreeType> const &config);

// This function resets weights for the different move operations if the user decides to change the fraction of loops
// used to learn the parameters of the model.
// The weights determine how much time is spend in the different move/parameter estimation steps.
template<typename TTreeType>
void setMoveProbs(Config<TTreeType> &config, double newParamsEstimateRate) {

    config.moveProbs[0] = config.moveProbs[0] / (1.0 - config.moveProbs[3]);
    config.moveProbs[0] = config.moveProbs[0] * (1.0 - newParamsEstimateRate);
    config.moveProbs[1] = config.moveProbs[1] / (1.0 - config.moveProbs[3]);
    config.moveProbs[1] = config.moveProbs[1] * (1.0 - newParamsEstimateRate);
    config.moveProbs[2] = config.moveProbs[2] / (1.0 - config.moveProbs[3]);
    config.moveProbs[2] = config.moveProbs[2] * (1.0 - newParamsEstimateRate);
    config.moveProbs[3] = 1 - config.moveProbs[0] - config.moveProbs[1] - config.moveProbs[2];

    return;
}


// Function to read the parameters from the command line.
template<typename TTreeType>
int readParameters(Config<TTreeType> &config, int argc, char *argv[]) {

    boost::program_options::options_description generic("Generic options");
    generic.add_options()("help,h", "Print this help");

    // configuration options
    boost::program_options::options_description parseConfig("Configuration");
    parseConfig.add_options()
            ("im", boost::program_options::value<std::string>(&config.inFileName), 
             "input mpileup")
            ("il", boost::program_options::value<decltype(config.loadName)>(&config.loadName),
             "Directory from which to read intermediate results.")
            ("ol", boost::program_options::value<decltype(config.saveName)>(&config.saveName),
             "Directory to store intermediate results.")
            ("ex", boost::program_options::value<decltype(config.exclusionFileName)>(&config.exclusionFileName),
             "Filename of exclusion list (VCF format).")
            ("me", boost::program_options::value<decltype(config.mutationExclusionFileName)>(
                    &config.mutationExclusionFileName),
             "Filename of mutations to exclude during the sequencing error rate estimation (VCF format).")
            ("in", boost::program_options::value<decltype(config.cellInfo)>(&config.cellInfo),
             "Name of the BAM files used to create the mpileup.")
            (",o", boost::program_options::value<decltype(config.outFilePrefix)>(&config.outFilePrefix),
             "Prefix of output files.")
            (",l", boost::program_options::value<decltype(config.loops)>(&config.loops),
             "Maximal number of iterations per repetition. [1000000]")
            ("lz", boost::program_options::value<decltype(config.learnZygocity)>(&config.learnZygocity),
             "Set to 1 if zygosity should be learned. [0]")
            ("zyg", boost::program_options::value<double>(&std::get<0>(config.params[Config<TTreeType>::E_nu])),
             "Zygosity rate. [0]")
            ("sa", boost::program_options::value<decltype(config.sampling)>(&config.sampling), "Sampling step. "
             "If a value x different from 0 is specified, every x iteration, after the burn in phase, an index "
             "will be writen to disk to provide a posterior sampling. [0]")
            ("inc", boost::program_options::value<decltype(config.variantInclusionFileName)>(&config.variantInclusionFileName), 
             "File name of inclusion list (VCF format) containing Variants (CHROM, POS, REF, ALT) that should be included.")
            ("ls", boost::program_options::value<decltype(config.sampleLoops)>(&config.sampleLoops),
             "Number of sample iterations. [100000]")
            ("pr", boost::program_options::value<decltype(config.priorMutationRate)>(&config.priorMutationRate),
             "Prior mutation rate [0.0001].")
            ("uniq,u", boost::program_options::value<decltype(config.uniqTreshold)>(&config.uniqTreshold),
             "Filter mutations showing up to this number of cells from the VCF. [0]")
            (",e", boost::program_options::value<decltype(config.paramsEstimateRate)>(&config.paramsEstimateRate),
             "Paramter estimation rate, i.e. the fraction of loops used to estimate the different parameters. [0.2]")
            ("ur", boost::program_options::value<double>(&std::get<1>(config.dataUsageRate)),
             "Data usage rate increment steps. [0.1]")
            ("seed", boost::program_options::value<decltype(config.fixedSeed)>(&config.fixedSeed),
             "Seed for the random number generator. [42]")
            (",t", boost::program_options::value<decltype(config.scoreType)>(&config.scoreType),
             "Tree score type [m (max), s (sum)]. [s]")
            ("wildOverDis", boost::program_options::value<double>(&std::get<0>(config.params[0])),
             "Overdispersion for wild type. [100]")
            ("mutationOverDis", boost::program_options::value<double>(&std::get<0>(config.params[1])),
             "Overdispersion for mutant type. [2]")
            ("wildMean", boost::program_options::value<double>(&std::get<0>(config.params[2])),
             "Mean error rate. [0.001]")
            ("md", boost::program_options::value<decltype(config.windowSize)>(&config.windowSize),
             "Window size for maximum number of allowed mutations. [1]")
            ("sub", boost::program_options::value<decltype(config.sub)>(&config.sub), "PCR substitution rate. [0]")
            ("mmw", boost::program_options::value<decltype(config.maxMutPerWindow)>(&config.maxMutPerWindow),
             "Maximum number of mutations allowed per window. [1]")
            ("cwm",
             boost::program_options::value<decltype(config.numCellWithMutationMin)>(&config.numCellWithMutationMin),
             "Number of cells requiered to have a mutation in order to be called. [2]")
            ("ncf", boost::program_options::value<decltype(config.normalCellFilter)>(&config.normalCellFilter),
             "Normal cell filter. Currently there are three options: (0) Do not use the normal cells for filtering; "
             "(1) use a simple filtering scheme excluding mutations if the probability of being mutated is higher than "
             "not being mutated for any cell independently; (2) filter mutations where the probility that at least one "
             "cell is mutated is higher than no cell is mutated. Note that in contrast to (1) the cells are not "
             "independent and cells with no alternative support need to be explained via dropout events. [2]")
            ("mc", boost::program_options::value<decltype(config.minCoverage)>(&config.minCoverage),
             "Minimum coverage required per cell. [1]")
            ("mnp",
             boost::program_options::value<decltype(config.minNumCellsPassFilter)>(&config.minNumCellsPassFilter),
             "Number of cells which need to pass the filters.")
            ("ms", boost::program_options::value<decltype(config.minSupport)>(&config.minSupport),
             "Minimum number of reads required to support the alternative. [0]")
            ("mf", boost::program_options::value<decltype(config.minFreq)>(&config.minFreq),
             "Minimum required frequency of reads supporting the alternative per cell. [0]")
            ("mff", boost::program_options::value<decltype(config.meanFilter)>(&config.meanFilter),
             "Mean of acceptable variant allele frequency across all cells for a specific locus. [0.25]")
            ("bns", boost::program_options::value<decltype(config.maxSupInControlBulk)>(&config.maxSupInControlBulk),
             "Loci with up to this number of alternative supporting reads in the bulk control sample will be skiped. "
             "[2]")
            ("bnc", boost::program_options::value<decltype(config.minCovInControlBulk)>(&config.minCovInControlBulk),
             "Minimum required coverage of reads in the bulk control sample. [6]")
            ("mnc", boost::program_options::value<decltype(config.maxNumberNormalCellMutated)>(
                    &config.maxNumberNormalCellMutated), "Maximum number of control cells allowed to be mutated. [0]")
            ("unc", boost::program_options::value<decltype(config.useNormalCellsInTree)>(&config.useNormalCellsInTree),
             "Use normal cells for tree reconstruction. [false]")
            ("ll", boost::program_options::value<decltype(config.computeLossScore)>(&config.computeLossScore),
             "Compute the loss score = allow a mutation to be lost in a subtree.")
            ("lp", boost::program_options::value<decltype(config.computeParallelScore)>(&config.computeParallelScore),
             "Compute the parallel score = allow a mutation to be acquired twice independently in the tree.")
            ("llp", boost::program_options::value<decltype(config.lossScorePenalty)>(&config.lossScorePenalty),
             "Penalty when computing the loss score.")
            ("lpp", boost::program_options::value<decltype(config.parallelScorePenalty)>(&config.parallelScorePenalty),
             "Penalty when computing the parallel score.")
            ("ese", boost::program_options::value<decltype(config.estimateSeqErrorRate)>(&config.estimateSeqErrorRate),
             "Estimate the sequencing error rate. [1]")
            ("mlm", boost::program_options::value<decltype(config.ml_mode)>(&config.ml_mode),
             "Use the maximum likelihood mode instead of the MCMC approach. [false]")
            ("chi", boost::program_options::value<decltype(config.learnChi)>(&config.learnChi),
             "Learn the loss and parallel priors. [false]");

    // hidden options, i.e., input files
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
        ("smt", boost::program_options::value<decltype(config.mutToMaxName)>(&config.mutToMaxName), "Store to save mutations distribution of MAP tree.");

    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(parseConfig).add(hidden);
    boost::program_options::options_description visible("Allowed options");
    visible.add(generic).add(parseConfig);
    boost::program_options::variables_map global_options;

    /* parse program options */
    try {
        boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
                options(cmdline_options).
                run(), global_options);

        // show help options
        if (global_options.count("help")) {
            std::cout << visible << '\n';
            exit(EXIT_SUCCESS);
        }

        boost::program_options::notify(global_options);
    }
    catch (boost::program_options::error &e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        exit(EXIT_FAILURE);
    }

    // Reset the weights for the tree move operations and parameter estimation if the user resets the parameter
    // estimation rate.
    if (global_options.count("-e")) {
        setMoveProbs(config, config.paramsEstimateRate);
    }

    if (config.lossScorePenalty == -1)
    {
        std::cout << "Not computing loss score!" << std::endl;
        config.computeLossScore = false;
    }
    if (config.parallelScorePenalty == -1)
    {
        std::cout << "Not computing parallel score!" << std::endl;
        config.computeParallelScore = false;
    }

    config.bestName = config.outFilePrefix + "/best_index/";
    config.lastName = config.outFilePrefix + "/last_index/";
    config.samplingName = config.outFilePrefix + "/sample_index";

    // backup the initial parameters
    // Note: At all times a copy of all parameters and the tree structure is kept, such that in case the new tree in
    // the mcmc traversal is rejected the origianl value can be retained.
    std::get<1>(config.params[0]) = std::get<0>(config.params[0]);
    std::get<1>(config.params[1]) = std::get<0>(config.params[1]);
    std::get<1>(config.params[2]) = std::get<0>(config.params[2]);
    std::get<1>(config.params[3]) = std::get<0>(config.params[3]);
    std::get<1>(config.params[4]) = std::get<0>(config.params[4]);
    std::get<1>(config.params[5]) = std::get<0>(config.params[5]);

    return 0;
}

template<typename TTreeType>
void
learnChi(typename Config<TTreeType>::TGraph &bestTree,
        std::array<std::tuple<double, double>, 9> &bestParams,
        Config<TTreeType> &config,
        std::vector<std::vector<unsigned>> &sampleTrees) {
    
    bool ml_mode = config.ml_mode;
    unsigned sampleLoops = config.sampleLoops;
    config.lossScorePenalty = 100000;
    config.parallelScorePenalty = 100000;
    config.sampleLoops = 0;
    //config.ml_mode = true;
    
    runMCMC(bestTree,bestParams, config, sampleTrees);

    //config.ml_mode = true;
    //runMCMC(bestTree,bestParams, config, sampleTrees);

    //config.ml_mode = ml_mode;
    config.sampleLoops = sampleLoops;

    std::cout << "The loss penalty score is " << config.lossScorePenalty << " and the parallel penalty score is: " << config.parallelScorePenalty << std::endl;

}

int main(int argc, char *argv[]) {

    typedef Config<SampleTree> TConfig;
    TConfig config{};

    // read the command line arguments
    std::cout << "SCIPhIN v" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << std::endl;
    std::cout << "Reading the config file: ... " << std::flush;
    readParameters(config, argc, argv);
    std::cout << "done!" << std::endl;

    /* initialize the random number generator, either with a user defined seed, or a random number */
    if (config.fixedSeed == static_cast<unsigned>(-1)) {
        initRand();
    } else {
        srand(config.fixedSeed);
    }

    // if the user wants to read a stored index rather than start a new analysis extract the mutation data from the
    // index files
    if (config.loadName == "") {
        std::cout << "Reading the mpileup file: " << std::flush;
        getData(config);
        createInitialTree(config);
        std::cout << "done!" << std::endl;
    } else {
        std::cout << "Reading the stored results file: " << std::flush;
        readCellNames(config);

        //check if tree is provided
        std::string treePath = config.loadName + "/tree.gv";
        std::ifstream f(treePath.c_str());
        if (f.good())
        {
            readGraph(config);
        }
        else
        {
            std::cout << "Creating random tree since no tree was provided.";
            createInitialTree(config);
        }

        readNucInfo(config);
        config.getTmpAttachmentScore().resize(2 * config.getNumSamples() - 1);
        std::cout << "done!" << std::endl;
    }

    // After the data is read the log scores need to be initialised.
    computeLogScoresOP(config);

    // paramsCounter is a data structure to store the model parameters in each iteration of the posterior sampling.
    // It is used to provide the median, mean, and sd of each parameter
    // Initialize the vectors of paramsCounter
    config.paramsCounter.resize(0);

    if (config.mutToMaxName != "")
    {
        config.updateContainers(0);
        writeMutToMaxTree(config);
        return 0;
    }
    // As with the paramsCounter (see above) every mutation to cell probability during the posterior sampling is stored.
    config.initMutInSampleCounter();

    //optimal tree
    typename TConfig::TGraph optimalTree;
    // optimal parameters
    std::array<std::tuple<double, double>, 9> optimalParams;
    // list where tree samples are stored, if sampling based on posterior distribution is needed
    std::vector<std::vector<unsigned>> sampleTrees;

    if (config.learnChi){
        learnChi(optimalTree, optimalParams, config, sampleTrees);
    }

    // Find best scoring trees by MCMC
    std::vector<std::queue<double>> rootProbabilities(config.getNumMutations());
    runMCMC(optimalTree, optimalParams, config, sampleTrees, rootProbabilities);

    // Write the results to disk
    normalizeMutationCounts(config);
    printMutationProbability(config, config.outFilePrefix + ".probs");
    printRootProbability(config, config.outFilePrefix + ".rootProbs", rootProbabilities);
    writeParameters(config, config.outFilePrefix + ".params.txt");
    writeVCF(config, config.outFilePrefix + ".vcf");

    // Write the index of the final iteration to disk, such that it can be used to continue from.
    writeFinalIndex(config);

    // Get the best tree information
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, Vertex<SimpleTree>> TGraph;
    config.getTree() = optimalTree;
    config.params = optimalParams;
    computeLogScoresOP(config);

    // Write the mutation to node assignment to disk (note: this is the most likely assignment (MAP))
    printMutation2SampleAssignment(config, config.outFilePrefix + "_mut2Sample.tsv");

    // write a simplified tree, where the nodes are summarized if possible, to disk
    TGraph newTreeBest = simplifyTree(config);
    std::ofstream ofs2(config.outFilePrefix + ".gv");
    write_graphviz(ofs2,
            newTreeBest,
            my_label_writer(newTreeBest,
                    config.indexToPosition,
                    config.cellNames,
                    config.cellColours,
                    config.cellClusters));
    std::cout << "SCIPhIN complete" << std::endl;
}
