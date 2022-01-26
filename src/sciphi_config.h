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
#ifndef CONFIG_H
#define CONFIG_H

#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/dynamic_bitset.hpp>

#include <array>
#include <vector>
#include <tuple>
#include <random>
#include <utility>
#include <cassert>
#include <unordered_map>

#include "attachmentScores.h"
#include "passDownAttachmentScores.h"
#include "noise_counts.h"
#include "logScores.h"

// SCIPhIN is based on the SampleTree, which is a binary cell lineage tree that stores the cells as leafs and where the
// mutations are attached to the nodes.
struct SampleTree {
};

// The mutation Tree stors the mutations as leafs. This tree representation will be used for samples that consist of
// many cells but only few mutations. At the moment this functionality is not supported.
struct MutationTree {
};

// A simplified tree where the nodes are condensed, used for printing purposes.
struct SimpleTree {
};

template<typename TTreeType>
struct Vertex {
};

template<>
struct Vertex<MutationTree> {
    unsigned mutation;
};

template<>
struct Vertex<SampleTree> {
    int sample = -1;  // -1 indicates this is an inner node
};
template<>
struct Vertex<SimpleTree> {
    std::vector<unsigned> mutations;
    int sample = -1;
};

// This class is used to get statistics on the parameter learned
struct ParamsCounter {
    std::vector<double> mu;                 // drop out rate
    std::vector<double> nu;                 // zygousity rate
    std::vector<double> lambdaWildLoss;     // fraction of mutation losses
    std::vector<double> lambdaMutLoss;      // fraction of mutation losses
    std::vector<double> parallel;           // fraction of parallel mutations
    std::vector<double> wildAlpha;          // beta-binomial alpha parameter of the wild type (sequencing errors)
    std::vector<double> wildBeta;           // beta-binomial beta parameter of the wild type (sequencing errors)
    std::vector<double> mutAlpha;           // beta-binomial alpha parameter of the mutation model
    std::vector<double> mutBeta;            // beta-binomial alpha parameter of the mutation model

    ParamsCounter() {};

    unsigned resize(unsigned newSize) {
        mu.resize(newSize, 0);
        nu.resize(newSize, 0);
        lambdaWildLoss.resize(newSize, 0);
        lambdaMutLoss.resize(newSize, 0);
        parallel.resize(newSize, 0);
        wildAlpha.resize(newSize, 0);
        wildBeta.resize(newSize, 0);
        mutAlpha.resize(newSize, 0);
        mutBeta.resize(newSize, 0);
        return newSize;
    }
};

// This class combines all data structures such that they dont have to be passed between functions
// TODO: create several sub-classes and use OOP!
template<typename TTreeType>
class Config {

public:

    // Some type definitions
    typedef boost::adjacency_list<boost::vecS,
            boost::vecS,
            boost::bidirectionalS,
            Vertex<TTreeType>> TGraph;
    //typedef std::vector<bool>                                           TSampleNodes;
    typedef LogScores TLogScores;
    typedef AttachmentScores TAttachmentScores;
    typedef PassDownAttachmentScores TPassDownAttachmentScores;
    typedef std::vector<std::vector<std::tuple<unsigned, unsigned> > > TData;
    typedef std::tuple<double, double> TParamsTuple;
    typedef std::array<TParamsTuple, 9> TParams;
    typedef std::tuple<double, unsigned, unsigned> TLearningParamsTuple;
    typedef std::vector<TLearningParamsTuple> TLearningParams;
    typedef boost::dynamic_bitset<> TBitSet;

    // A practical enum for easy access to the parameters
    enum ParamType {
        E_wildOverDis = 0,
        E_mutationOverDis = 1,
        E_wildMean = 2,
        E_mu = 3,
        E_nu = 4,
        E_lambdaWildLoss = 5,
        E_lambdaMutLoss = 6,
        E_parallel = 7,
        E_mutationMean = 8
    };

    // Many small helper functions, mostly getter and setter
    void updateParamsCounter();

    TGraph &getTree();

    TGraph const &getTree() const;

    void setTree(TGraph const &newTree);

    TGraph &getTmpTree();

    void setTmpTree(TGraph const &newTree);

    unsigned getNumSamples() const;

    void setNumSamples(unsigned newNumSamples);

    ParamType getParamToOptimize() const;

    void setParamToOptimize(ParamType newParamToOptimize);

    unsigned getNumMutations() const;

    unsigned getNumAttachments() const;

    TLogScores &getLogScores();

    TLogScores const &getLogScores() const;

    void setLogScores(TLogScores newLogScores);

    TLogScores &getTmpLogScores();

    void setTmpLogScores(TLogScores &newTmpLogScores);

    TData &getData();

    TData const &getData() const;

    TData &getCompleteData();

    TData const &getCompleteData() const;

    TAttachmentScores &getTmpAttachmentScore();

    TAttachmentScores const &getTmpAttachmentScore() const;

    double getParam(ParamType param);

    double getParam(ParamType param) const;

    void setParam(ParamType param, double newParam);

    double getTmpParam(ParamType param);

    void setTmpParam(ParamType param, double newParam);

    double getSDParam(ParamType param);

    double getSDParam(ParamType param) const;

    void setSDParam(ParamType param, double newParam);

    unsigned getSDCountParam(ParamType param);

    unsigned getSDCountParam(ParamType param) const;

    void setSDCountParam(ParamType param, unsigned newParam);

    unsigned getSDTrialsParam(ParamType param);

    unsigned getSDTrialsParam(ParamType param) const;

    void setSDTrialsParam(ParamType param, unsigned newParam);

    unsigned getMoveTyp();

    void setMoveTyp(unsigned newMoveType);

    // This function is used to reset a newly drawn parameter if the proposed tree is rejected
    void resetParameters();

    // Print the learned parameters of the model
    void printParameters();

    // This function is used to increase the number of mutations used for the tree inference over time.
    // Especially for trees with many more mutations than cells this approach speeds up the MCMC scheme.
    bool updateContainers(unsigned currentLoop);

    // Initialize the counter of the mutation to node assignments
    void initMutInSampleCounter();

    std::default_random_engine &getGenerator();

    // THe members of the config class
    // The main tree structure - we make use of a pair to restore a tree if the proposed one is discarded
    std::pair<TGraph, TGraph> tree;

    // The parameters of the tree
    TParams params;

    // Substitution error, especially important for high coverage data sets to allow for early MDA errors
    double sub;

    // The parameters for learning the tree parameters
    TLearningParams learningParams;

    // The probabilities with which one of the move (including parameter estimation) is performed
    std::array<double, 4> moveProbs;

    // Random number generator
    std::default_random_engine generator;

    // Number of mcmc iteration
    unsigned loops;

    // Fixed, user specified seed
    unsigned fixedSeed;

    // Sum the probabilities of the mutation to node assignment or use the maximum
    char scoreType;

    // Rate with which the parameters and not the tree structure is learned
    double paramsEstimateRate;

    // The expected mutation rate
    double priorMutationRate;

    // The expected germline rate, important when it comes to single cell filtering based on single normal cells
    double priorGermlineRate;

    // Filter mutations showing up to this number of cells from the VCF.
    unsigned uniqTreshold;

    // Rate of data already used
    std::tuple<double, double> dataUsageRate;

    // Number of posterior sampling moves
    unsigned sampleLoops;
    // Positions used to estimate the sequencing error rate
    unsigned errorRateEstLoops;

    // The noise counts
    NoiseCounts noiseCounts;

    // The log score of the wild type, heterozygous, and homzygous model
    // We use a pair here to sote a copy in case the new tree is rejected
    std::pair<TLogScores, TLogScores> logScores;

    // The nucleotid counts currently in use
    TData data;

    // The complete nucleotide counts
    TData completeData;
    // TODO: delete this, it can be replaced by a static varaible
    TAttachmentScores _tmpAttachmentScore;

    // The number of samples
    unsigned numSamples;

    // The current move/parameter learning type
    unsigned moveType;

    // The current model parameter to learn
    ParamType paramToOptimize;

    // Prefix of all output files
    std::string outFilePrefix;

    // Infos, such as cell name and type (tumor, normal)
    std::string cellInfo;

    // List of positions to be excluded (e.g., dbSNP posisitions)
    std::string exclusionFileName;

    // List of positions to be excluded from the error learning, e.g., likely mutated positions in a panel data set
    std::string mutationExclusionFileName;

    // List of positions to be included
    std::string variantInclusionFileName;

    // The mpileup name
    std::string inFileName;

    // The index name to load
    std::string loadName;

    // The index name to store
    std::string saveName;

    // Name of the best index
    std::string bestName;
    
    // Name of the last index
    std::string lastName;

    // Name of the last index
    std::string samplingName;

    // Cell names
    std::vector<std::string> cellNames;

    // File to save mutations distribution of MAP tree.
    std::string mutToMaxName;

    // Colours of cells
    std::vector<std::string> cellColours;

    // Cluster id of cell
    std::vector<unsigned> cellClusters;

    // For each mutation the chromosome, the position on the chromosome, the reference, and the alternative allel
    std::vector<std::tuple<std::string, unsigned, char, char>> indexToPosition;

    // During the posterior sampling the mutation probabilities are summed up and normalized later
    std::vector<std::vector<TAttachmentScores::TAttachmentScore>> mutInSampleCounter;

    // Window size for maximum number of allowed mutations
    int windowSize;

    // Maximum number of mutations allowed per window
    unsigned maxMutPerWindow;

    // Minimal number of cell required to show the mutation
    unsigned numCellWithMutationMin;

    // Filtering scheme to apply when normal cells are present
    unsigned normalCellFilter;

    // Min coverage for a cell to be considered for candidate loci identification
    unsigned minCoverage;

    // Minimum number of cells that have to pass the filters
    unsigned minNumCellsPassFilter;

    // Minimum support required for a cell
    unsigned minSupport;

    // Minimum frequency required for a cell
    double minFreq;

    // Minimum coverage in bulk required for a position to be considered
    unsigned minCovInControlBulk;

    // Number of reads allowed to support the alternative in a control bulk sample
    unsigned maxSupInControlBulk;

    // Data structure to store the learned parameters of the model during the posterior sampling to obtain overview
    // statistics, sucha as median, mean, and sd
    ParamsCounter paramsCounter;

    // The score associated to noise (sequencing errors)
    double noiseScore;

    // Learn the homzygousity rate of the experiment
    bool learnZygocity;

    // Indicator if the computation of losing a mutation should be included
    bool computeLossScore;

    // Indicator if the computation of parallel a mutations should be included
    bool computeParallelScore;

    // Penalty when computing the loss score
    double lossScorePenalty;
    
    // Penalty when computing the lparallel score
    double parallelScorePenalty;

    // Indicator if the sequencing error rate should be learned
    bool estimateSeqErrorRate;

    // Mutations with a mean alternative frequency of this are  by means of learning two beta-binomials
    double meanFilter;

    // Minimum coverage required for normal cells
    unsigned minCovNormalCell;

    // Maximum number of normal cells allowed to support the alternative allele
    unsigned maxNumberNormalCellMutated;

    // Indicator to decide whether the normal cells should be included in the tree structur learning
    bool useNormalCellsInTree;
    
    unsigned sampling;


    Config() :
            params{{TParamsTuple{100.0, 100.0},             //overdispersion background
                           TParamsTuple{2, 2},              //overdispersion mutation
                           TParamsTuple{0.001, 0.001},      //sequencing error rate
                           TParamsTuple{0.9, 0.9},          //drop out rate
                           TParamsTuple{0, 0},              //zygousity rate
                           TParamsTuple{0, 0},              //loss rate
                           TParamsTuple{0, 0}}},            //parallele rate
            sub(0),
            learningParams{{TLearningParamsTuple{5.0, 0, 0},
                                   TLearningParamsTuple{0.1, 0, 0},
                                   TLearningParamsTuple{0.01, 0, 0},
                                   TLearningParamsTuple{0.01, 0, 0},
                                   TLearningParamsTuple{0.01, 0, 0},
                                   TLearningParamsTuple{0.01, 0, 0},
                                   TLearningParamsTuple{0.01, 0, 0},
                                   TLearningParamsTuple{0.01, 0, 0}}},
            moveProbs{{0.64, 0.16, 0, 0.2}},
            generator(42),
            loops(1000000),
            fixedSeed(42),
            scoreType('s'),
            paramsEstimateRate(0.2),
            priorMutationRate(0.0001),
            priorGermlineRate(0.001),
            uniqTreshold(0),
            dataUsageRate{0, 1},
            sampleLoops(100000),
            errorRateEstLoops(100000),
            outFilePrefix("sciphi"),
            windowSize(-1),
            maxMutPerWindow(1),
            numCellWithMutationMin(2),
            normalCellFilter(2),
            minCoverage(1),
            minNumCellsPassFilter(2),
            minSupport(3),
            minFreq(0),
            minCovInControlBulk(6),
            maxSupInControlBulk(2),
            noiseScore(0),
            learnZygocity(false),
            computeLossScore(false),
            computeParallelScore(false),
            lossScorePenalty(0),
            parallelScorePenalty(0),
            estimateSeqErrorRate(true),
            meanFilter(0.25),
            minCovNormalCell(5),
            maxNumberNormalCellMutated(0),
            useNormalCellsInTree(false) {};
};

template<typename TTreeType>
typename Config<TTreeType>::TGraph &
Config<TTreeType>::getTree() { return this->tree.first; }

template<typename TTreeType>
typename Config<TTreeType>::TGraph const &
Config<TTreeType>::getTree() const { return this->tree.first; }

template<typename TTreeType>
void
Config<TTreeType>::setTree(Config::TGraph const &newTree) { this->tree.first = newTree; }

// During the posterior sampling the learned values of the different parameters are stored using this function
template<typename TTreeType>
void
Config<TTreeType>::updateParamsCounter() {
    this->paramsCounter.wildAlpha.push_back(this->getParam(Config::E_wildOverDis) *
                                            this->getParam(Config::E_wildMean));
    this->paramsCounter.wildBeta.push_back(this->getParam(Config::E_wildOverDis) -
                                           this->getParam(Config::E_wildMean) *
                                           this->getParam(Config::E_wildOverDis));
    this->paramsCounter.mutAlpha.push_back(this->getParam(Config::E_mutationOverDis) *
                                           this->getParam(Config::E_mutationMean));
    this->paramsCounter.mutBeta.push_back(this->getParam(Config::E_mutationOverDis) -
                                          (0.5 - this->getParam(Config::E_wildMean)) *
                                          this->getParam(Config::E_mutationOverDis));
    this->paramsCounter.mu.push_back(this->getParam(Config::E_mu));
    this->paramsCounter.nu.push_back(this->getParam(Config::E_nu));
    this->paramsCounter.lambdaWildLoss.push_back(this->getParam(Config::E_lambdaWildLoss));
    this->paramsCounter.lambdaMutLoss.push_back(this->getParam(Config::E_lambdaMutLoss));
    this->paramsCounter.parallel.push_back(this->getParam(Config::E_parallel));
}

// Get the backup tree (in case the new tree is discarded the backup can be used without need for re-computation)
template<typename TTreeType>
typename Config<TTreeType>::TGraph &
Config<TTreeType>::getTmpTree() { return this->tree.second; }

// Set the backup tree (in case the new tree is discarded the backup can be used without need for re-computation)
template<typename TTreeType>
void
Config<TTreeType>::setTmpTree(Config::TGraph const &newTree) { this->tree.second = newTree; }

template<typename TTreeType>
unsigned
Config<TTreeType>::getNumSamples() const { return this->numSamples; }

template<typename TTreeType>
void
Config<TTreeType>::setNumSamples(unsigned newNumSamples) { this->numSamples = newNumSamples; }

template<typename TTreeType>
typename Config<TTreeType>::ParamType
Config<TTreeType>::getParamToOptimize() const { return this->paramToOptimize; }

template<typename TTreeType>
void
Config<TTreeType>::setParamToOptimize(
        Config<TTreeType>::ParamType newParamToOptimize) { this->paramToOptimize = newParamToOptimize; }

template<typename TTreeType>
unsigned
Config<TTreeType>::getNumMutations() const {
    assert(this->getData().size() > 0);
    return this->getData()[0].size();
}

template<typename TTreeType>
unsigned
Config<TTreeType>::getNumAttachments() const {
    assert(this->getData().size() > 0);
    return this->getData()[0].size();
}

template<typename TTreeType>
typename Config<TTreeType>::TLogScores &
Config<TTreeType>::getLogScores() { return this->logScores.first; }

template<typename TTreeType>
typename Config<TTreeType>::TLogScores const &
Config<TTreeType>::getLogScores() const { return this->logScores.first; }

template<typename TTreeType>
void
Config<TTreeType>::setLogScores(Config<TTreeType>::TLogScores newLogScores) { this->logScores.first = newLogScores; }

template<typename TTreeType>
typename Config<TTreeType>::TLogScores &
Config<TTreeType>::getTmpLogScores() { return this->logScores.second; }

template<typename TTreeType>
void
Config<TTreeType>::setTmpLogScores(
        Config<TTreeType>::TLogScores &newTmpLogScores) { this->logScores.second = newTmpLogScores; }

template<typename TTreeType>
typename Config<TTreeType>::TData &
Config<TTreeType>::getData() { return this->data; }

// get the currently used nucleotide count information
template<typename TTreeType>
typename Config<TTreeType>::TData const &
Config<TTreeType>::getData() const { return this->data; }

// get all nucleotide count information
template<typename TTreeType>
typename Config<TTreeType>::TData &
Config<TTreeType>::getCompleteData() { return this->completeData; }

template<typename TTreeType>
typename Config<TTreeType>::TData const &
Config<TTreeType>::getCompleteData() const { return this->completeData; }

// TODO: remove this as it can be realised with a static variable
template<typename TTreeType>
typename Config<TTreeType>::TAttachmentScores &
Config<TTreeType>::getTmpAttachmentScore() { return this->_tmpAttachmentScore; }

template<typename TTreeType>
typename Config<TTreeType>::TAttachmentScores const &
Config<TTreeType>::getTmpAttachmentScore() const { return this->_tmpAttachmentScore; }

template<typename TTreeType>
double
Config<TTreeType>::getParam(Config<TTreeType>::ParamType param) {
    // if the sequencing error rate is required compute it
    if (param == this->E_mutationMean)
        return 0.5 - this->getParam(this->E_wildMean);
    return std::get<0>(this->params[param]);
}

template<typename TTreeType>
double
Config<TTreeType>::getParam(Config<TTreeType>::ParamType param) const {
    if (param == this->E_mutationMean)
        return 0.5 - this->getParam(this->E_wildMean);
    return std::get<0>(this->params[param]);
}

template<typename TTreeType>
void
Config<TTreeType>::setParam(Config<TTreeType>::ParamType param, double newParam) {
    std::get<0>(this->params[param]) = newParam;
}

template<typename TTreeType>
double
Config<TTreeType>::getTmpParam(Config<TTreeType>::ParamType param) { return std::get<1>(this->params[param]); }

template<typename TTreeType>
void
Config<TTreeType>::setTmpParam(Config<TTreeType>::ParamType param, double newParam) {
    std::get<1>(this->params[param]) = newParam;
}

template<typename TTreeType>
double
Config<TTreeType>::getSDParam(Config<TTreeType>::ParamType param) { return std::get<0>(this->learningParams[param]); }

template<typename TTreeType>
double
Config<TTreeType>::getSDParam(Config<TTreeType>::ParamType param) const {
    return std::get<0>(this->learningParams[param]);
}

template<typename TTreeType>
void
Config<TTreeType>::setSDParam(Config<TTreeType>::ParamType param, double newParam) {
    std::get<0>(this->learningParams[param]) = newParam;
}

template<typename TTreeType>
unsigned
Config<TTreeType>::getSDCountParam(Config<TTreeType>::ParamType param) {
    return std::get<1>(this->learningParams[param]);
}

template<typename TTreeType>
unsigned
Config<TTreeType>::getSDCountParam(Config<TTreeType>::ParamType param) const {
    return std::get<1>(this->learningParams[param]);
}

template<typename TTreeType>
void
Config<TTreeType>::setSDCountParam(Config<TTreeType>::ParamType param, unsigned newParam) {
    std::get<1>(this->learningParams[param]) = newParam;
}

template<typename TTreeType>
unsigned
Config<TTreeType>::getSDTrialsParam(Config<TTreeType>::ParamType param) {
    return std::get<2>(this->learningParams[param]);
}

template<typename TTreeType>
unsigned
Config<TTreeType>::getSDTrialsParam(Config<TTreeType>::ParamType param) const {
    return std::get<2>(this->learningParams[param]);
}

template<typename TTreeType>
void
Config<TTreeType>::setSDTrialsParam(Config<TTreeType>::ParamType param, unsigned newParam) {
    std::get<2>(this->learningParams[param]) = newParam;
}

template<typename TTreeType>
unsigned
Config<TTreeType>::getMoveTyp() { return this->moveType; }

template<typename TTreeType>
void
Config<TTreeType>::setMoveTyp(unsigned newMoveType) { this->moveType = newMoveType; }

template<typename TTreeType>
std::default_random_engine &
Config<TTreeType>::getGenerator() { return this->generator; }

// This function is used to reset a newly drawn parameter if the proposed tree is rejected
template<typename TTreeType>
void
Config<TTreeType>::resetParameters() {
    switch (this->getParamToOptimize()) {
        case (E_wildOverDis) : {
            this->setParam(E_wildOverDis, this->getTmpParam(E_wildOverDis));
            break;
        }
        case (E_mutationOverDis) : {
            this->setParam(E_mutationOverDis, this->getTmpParam(E_mutationOverDis));
            break;
        }
        case (E_wildMean) : {
            this->setParam(E_wildMean, this->getTmpParam(E_wildMean));
            break;
        }
        case (E_mu) : {
            this->setParam(E_mu, this->getTmpParam(E_mu));
            break;
        }
        case (E_nu) : {
            this->setParam(E_nu, this->getTmpParam(E_nu));
            break;
        }
        case (E_lambdaWildLoss) : {
            this->setParam(E_lambdaWildLoss, this->getTmpParam(E_lambdaWildLoss));
            break;
        }
        case (E_lambdaMutLoss) : {
            this->setParam(E_lambdaMutLoss, this->getTmpParam(E_lambdaMutLoss));
            break;
        }
        case (E_parallel) : {
            this->setParam(E_parallel, this->getTmpParam(E_parallel));
            break;
        }
        case (E_mutationMean) : {
            assert(false);
            break;
        }
    }
    this->setLogScores(this->getTmpLogScores());
};

template<typename TTreeType>
void
Config<TTreeType>::initMutInSampleCounter() {
    this->mutInSampleCounter.resize(this->getNumSamples());
    for (size_t i = 0; i < this->mutInSampleCounter.size(); ++i) {
        this->mutInSampleCounter[i].resize(this->getCompleteData()[0].size());
        for (unsigned j = 0; j < this->mutInSampleCounter[i].size(); ++j) {
            this->mutInSampleCounter[i][j] = AttachmentScore();
        }
    }
}

// This function is used to increase the number of mutations used for the tree inference over time
template<typename TTreeType>
bool
Config<TTreeType>::updateContainers(unsigned currentLoop) {

    // check if more candidate mutations should be used for the tree inference
    if (currentLoop < this->loops &&
         static_cast<double>(currentLoop) / static_cast<double>(this->loops) >= std::get<0>(this->dataUsageRate)) {

        unsigned newDataSize = -1;
        // check if we are already in the posterior sampling phase
        if (currentLoop >= this->loops) {
            newDataSize = this->getCompleteData()[0].size();
        } else {
            // increase the data usage rate by a pre defined fraction
            std::get<0>(this->dataUsageRate) += std::get<1>(this->dataUsageRate);
            newDataSize = std::round(std::get<0>(this->dataUsageRate) * (this->getCompleteData()[0].size()));
        }

        std::cout << "newDataSize: " <<
                  newDataSize << " "
                  << std::get<0>(this->dataUsageRate) * (this->getCompleteData()[0].size()) << std::endl;

        // insert the new data points into the set of currently used ones
        size_t j = this->getNumMutations();
        for (; j < newDataSize; ++j) {
            for (size_t i = 0; i < this->getNumSamples(); ++i) {
                this->getData()[i].push_back(this->getCompleteData()[i][j]);
            }
        }

        // resize the log score vector
        this->getLogScores().resizeNumCells(this->getNumSamples());
        this->getTmpLogScores().resizeNumCells(this->getNumSamples());
        this->getLogScores().resizeNumMuts(this->getNumMutations());
        this->getTmpAttachmentScore().resize(2 * this->getNumSamples() - 1);

        // re-compute the log scores
        computeLogScoresOP(*this);
        this->setTmpLogScores(this->getLogScores());

        return true;
    }
    return false;
};

template<typename TTreeType>
void
Config<TTreeType>::printParameters() {

    std::cout << "num Samples:\t" << this->getNumSamples() << std::endl;
    std::cout << "total # mut:\t" << this->getCompleteData()[0].size() << "\tcurrently used:\t"
              << this->getNumMutations() << std::endl;
    std::cout << "seq error  - freq:    " << this->getParam(Config::E_wildMean) << " SD: "
              << this->getSDParam(Config::E_wildMean) << " count: "
              << this->getSDCountParam(Config::E_wildMean) << " trails: " << this->getSDTrialsParam(Config::E_wildMean)
              << std::endl;
    std::cout << "normal     - overDis: " << this->getParam(Config::E_wildOverDis) << " SD: "
              << this->getSDParam(Config::E_wildOverDis)
              << " count: " << this->getSDCountParam(Config::E_wildOverDis) << " trails: "
              << this->getSDTrialsParam(Config::E_wildOverDis) << std::endl;
    std::cout << "normal     - alpha:   " << this->getParam(Config::E_wildOverDis) * this->getParam(Config::E_wildMean)
              << " beta: " << this->getParam(Config::E_wildOverDis) -
                              this->getParam(Config::E_wildMean) * this->getParam(Config::E_wildOverDis) << std::endl;
    std::cout << "mutation   - overDis: " << this->getParam(Config::E_mutationOverDis) << " SD: "
              << this->getSDParam(Config::E_mutationOverDis)
              << " count: " << this->getSDCountParam(Config::E_mutationOverDis) << " trails: "
              << this->getSDTrialsParam(Config::E_mutationOverDis) << std::endl;
    std::cout << "mutation   - alpha: "
              << this->getParam(Config::E_mutationOverDis) * this->getParam(Config::E_mutationMean) << " beta: "
              << this->getParam(Config::E_mutationOverDis) -
                 (0.5 - this->getParam(Config::E_wildMean)) * this->getParam(Config::E_mutationOverDis) << std::endl;
    std::cout << "drop: " << 1 - this->getParam(Config::E_mu) << " SD: " << this->getSDParam(Config::E_mu) << " count: "
              << this->getSDCountParam(Config::E_mu) << " trails: " << this->getSDTrialsParam(Config::E_mu)
              << std::endl;
    std::cout << "zyg: " << this->getParam(Config::E_nu) << " SD: " << this->getSDParam(Config::E_nu) << " count: "
              << this->getSDCountParam(Config::E_nu) << " trails: " << this->getSDTrialsParam(Config::E_nu)
              << std::endl;
    std::cout << "lossWild: " << this->getParam(Config::E_lambdaWildLoss) << " SD: " << this->getSDParam(Config::E_lambdaWildLoss)
              << " count: " << this->getSDCountParam(Config::E_lambdaWildLoss) << " trails: "
              << this->getSDTrialsParam(Config::E_lambdaWildLoss) << std::endl;
    std::cout << "lossMut: " << this->getParam(Config::E_lambdaMutLoss) << " SD: " << this->getSDParam(Config::E_lambdaMutLoss)
              << " count: " << this->getSDCountParam(Config::E_lambdaMutLoss) << " trails: "
              << this->getSDTrialsParam(Config::E_lambdaMutLoss) << std::endl;
    std::cout << "para: " << this->getParam(Config::E_parallel) << " SD: " << this->getSDParam(Config::E_parallel)
              << " count: " << this->getSDCountParam(Config::E_parallel) << " trails: "
              << this->getSDTrialsParam(Config::E_parallel) << std::endl;
};

// Get the median, mean, and sd of all sampled parameters
std::array<double, 3> getStats(std::vector<double> &values) {
    std::array<double, 3> result;

    // Compute the median
    std::sort(values.begin(), values.end());
    result[0] = values[values.size() / 2];

    // Compute the sum
    double sum = 0;
    for (unsigned i = 0; i < values.size(); ++i) {
        sum += values[i];
    }
    result[1] = sum / static_cast<double>(values.size());

    // Compute the standard deviation
    double variance = 0;
    for (unsigned i = 0; i < values.size(); ++i) {
        variance += std::pow(values[i] - result[1], 2);
    }
    result[2] = std::sqrt(variance / static_cast<double>(values.size()));

    return result;
}

template<typename TTreeType>
void
writeParameters(Config<TTreeType> &config, std::string const &fileName) {
    std::ofstream outFile;
    outFile.open(fileName);
    outFile << "background frequency:\t" << config.getParam(Config<TTreeType>::E_wildMean) << std::endl;
    std::array<double, 3> stats = getStats(config.paramsCounter.wildAlpha);
    outFile << "wildAlpha:\t" << stats[0] << "\t" << stats[1] << "\t" << stats[2] << std::endl;
    stats = getStats(config.paramsCounter.wildBeta);
    outFile << "wildBeta:\t" << stats[0] << "\t" << stats[1] << "\t" << stats[2] << std::endl;
    outFile << "average mutation frequency:\t" << config.getParam(Config<TTreeType>::E_mutationMean) << std::endl;
    stats = getStats(config.paramsCounter.mutAlpha);
    outFile << "altAlpha:\t" << stats[0] << "\t" << stats[1] << "\t" << stats[2] << std::endl;
    stats = getStats(config.paramsCounter.mutBeta);
    outFile << "altBeta:\t" << stats[0] << "\t" << stats[1] << "\t" << stats[2] << std::endl;
    stats = getStats(config.paramsCounter.mu);
    outFile << "mu:\t" << stats[0] << "\t" << stats[1] << "\t" << stats[2] << std::endl;
    stats = getStats(config.paramsCounter.nu);
    outFile << "nu:\t" << stats[0] << "\t" << stats[1] << "\t" << stats[2] << std::endl;
    stats = getStats(config.paramsCounter.lambdaWildLoss);
    outFile << "lambdaWildLoss:\t" << stats[0]<< "\t" << stats[1] << "\t" << stats[2] << std::endl;
    stats = getStats(config.paramsCounter.lambdaMutLoss);
    outFile << "lambdaMutLoss:\t" << stats[0]<< "\t" << stats[1] << "\t" << stats[2] << std::endl;
    stats = getStats(config.paramsCounter.parallel);
    outFile << "parallel:\t" << stats[0]<< "\t" << stats[1] << "\t" << stats[2] << std::endl;
}

// this prints the graph as it is currently used
class my_label_writer_complete {
    public:

    my_label_writer_complete(boost::adjacency_list<boost::vecS,boost::vecS, boost::bidirectionalS, Vertex<SampleTree>> const & sampleTree_) :
        sampleTree(sampleTree_)
    {}

    template <class VertexOrEdge>
    void operator()(std::ostream& out, const VertexOrEdge& v) const {

        if (sampleTree[v].sample == -1)
        {
            out << "[label=\"" << v;
        }
        else
        {
            out << "[shape=box,label=\"" << sampleTree[v].sample + boost::num_vertices(sampleTree) / 2 - 1;
        }
        out << "\"]";
    }

    private:
    boost::adjacency_list<boost::vecS,boost::vecS, boost::bidirectionalS, Vertex<SampleTree>> const & sampleTree;
};

#endif
