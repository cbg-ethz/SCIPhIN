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
#ifndef MCMC_H
#define MCMC_H

#include "sciphi_config.h"
#include <vector>
#include <stdlib.h>
#include <cmath>
#include <string.h>
#include <float.h>
#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <limits.h>
#include "trees.h"
#include "rand.h"
#include "scoreTree.h"
#include "probabilities.h"

// This function keeps the coefficients between 0 and 1 by mirroring at 0 or 1
double adjusteCoefficient(double coEff)
{
    while(true)
    {
        if(coEff < 0.0)
        {
            coEff = std::abs(coEff);
        }
        else if (coEff > 1.0)
        {
            coEff =  1.0 - std::abs(1 - coEff);
        }
        else
        {
            break;
        }
    }
    return coEff;
}

// This function assures that the overdispersion of the heterozygous mutation cannot take a U form (alpha and beta
// are very close or bigger than 1)
// TODO: Can we omit this function since some data sets indeed show U shape nucleotide distributions
double adjustMutOverDis(double mutationOverDis)
{
    if (mutationOverDis < 2.0)
    {
        return 2.0 + (2.0 - mutationOverDis);
    }
    return mutationOverDis;
}

// Choose which parameter to optimize
template <typename TTreeType>
void
changeParameters(Config<TTreeType> & config)
{
    // store the old scores, in case the new parameter will be rejected
    config.setTmpLogScores(config.getLogScores());

    //choose which parameter to optimize
    unsigned numParams = 4 + 2 * config.computeLossScore + config.learnZygocity + config.computeParallelScore;
    unsigned paramIdx = rand() % numParams;
    if(paramIdx > 3) {
        if (paramIdx == 4 && !config.learnZygocity)
        {
            ++paramIdx;
        }
        if (paramIdx == 5 && !config.computeLossScore)
        {
            ++paramIdx;
        }
        if (paramIdx == 6 && !config.computeLossScore)
        {
            ++paramIdx;
        }
    }
    typename Config<TTreeType>::ParamType param = (typename Config<TTreeType>::ParamType)paramIdx;
    config.setParamToOptimize(param);


    // update the standard deviation of the chosen parameter every 1000 steps
    if (config.getSDTrialsParam(param) == 1000)
    {
        double sucRate = static_cast<double>(config.getSDCountParam(param))/static_cast<double>(config.getSDTrialsParam(param));
        config.setSDParam(param, config.getSDParam(param) * (std::exp((sucRate - 0.5)/10.0)));
        config.setSDCountParam(param, 0);
        config.setSDTrialsParam(param, 0);
    }
    config.setSDTrialsParam(param, config.getSDTrialsParam(param) + 1);

    // store old param value in case the new one will be rejected
    config.setTmpParam(param, config.getParam(param));
    double oldParamValue = config.getParam(param);

    // draw new parameter value
    std::normal_distribution<double> distribution(oldParamValue, config.getSDParam(param));
    double newParamValue = distribution(config.getGenerator());

    // If the zygocity, loss, or parallel rate is changed they might have to be normalized
    if (paramIdx > 3)
    {
        newParamValue = adjusteCoefficient(newParamValue);
        double sumParams = config.getParam(Config<TTreeType>::ParamType::E_nu)
                + config.getParam(Config<TTreeType>::ParamType::E_lambdaWildLoss)
                + config.getParam(Config<TTreeType>::ParamType::E_lambdaMutLoss)
                + config.getParam(Config<TTreeType>::ParamType::E_parallel);
        if (sumParams > 1)
        {
            for (unsigned i = 4; i < 8; ++i) {
                config.setParam(typename Config<TTreeType>::ParamType(i),
                        config.getParam(typename Config<TTreeType>::ParamType(i)) / sumParams);
            }
        }
    }
    else if (param == Config<TTreeType>::ParamType::E_mu)
    {
        newParamValue = adjusteCoefficient(newParamValue);
    }
    else if (param == Config<TTreeType>::ParamType::E_mutationOverDis)
    {
        newParamValue = adjustMutOverDis(newParamValue);
    }

    // if one of the wild type parameters is influenced the noise scores have to be re-computed
    config.setParam(param, std::abs(newParamValue));

    if (param == Config<TTreeType>::ParamType::E_wildOverDis || param == Config<TTreeType>::ParamType::E_wildMean)
    {
        computeNoiseScore(config);
    }

    // Re-compute the log scores with the new values
    computeLogScoresOP(config);
}

// Choose a tree changing move or the parameter adjustment
template <typename TTreeType>
void
proposeNextConfiguration(Config<TTreeType> & config)
{
	config.setMoveTyp(sampleRandomMove(config.moveProbs));      // pick the move type

	if(config.getMoveTyp() == 4) /* change one of the parameters */
    {
        changeParameters(config);
    }
    else
    {
        config.setTmpTree(config.getTree());
        if(config.getMoveTyp()==1){ /* prune and re-attach */
            pruneAndReAttach(config);
        }
        else if(config.getMoveTyp()==2){ /* swap two node labels  */
            swapNodeLabels(config);
        }
        else if(config.getMoveTyp()==3){  /*  swap two subtrees  */
            swapSubtrees(config);
        }
        else
        {
            std::cout << "This should not have happened!" << std::endl;
        }
    }
}

template <typename TTreeType>
void manageBestTrees(Config<TTreeType> & /*config*/,
                     double & /*bestTreeLogScore*/,
                     double /*currTreeLogScore*/,
                     TTreeType & /*bestTrees*/)
{}

// This function checks if the new score is better than any so far computed score and stores it to disk
inline
void
manageBestTrees(Config<SampleTree> & config,
                     double & bestTreeLogScore,
                     double currTreeLogScore,
                     Config<SampleTree>::TGraph & bestTree,
                     std::array<std::tuple<double, double>, 9> & bestParams)
{
    if(currTreeLogScore > bestTreeLogScore)
    {
        std::cout << "The new best score is: " <<  currTreeLogScore << std::endl;
        //save the current state to disk
        writeIndex(config, config.bestName);
        bestTreeLogScore = currTreeLogScore;
        bestTree = config.getTree();
        bestParams = config.params;

        //save the current best tree to disk
        boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, Vertex<SimpleTree>> newTreeBest = simplifyTree(config);
        std::ofstream ofs(config.outFilePrefix + ".gv");
        write_graphviz(ofs, newTreeBest, my_label_writer(newTreeBest, config.indexToPosition, config.cellNames, config.cellColours, config.cellClusters));
        ofs.close();
    }
}


// This function is actually doing the sampling of the mutations from the posterior distribution
inline
void
updateMutInSampleCounts(Config<SampleTree> & config)
{
    Config<SampleTree>::TAttachmentScores & attachmentScores = config.getTmpAttachmentScore();
    static Config<SampleTree>::TAttachmentScores attachmentSumScores;
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentScores;
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentSumScores;
    attachmentSumScores.resize(attachmentScores.size());
    passDownAttachmentScores.resize(attachmentScores.size());
    passDownAttachmentSumScores.resize(attachmentScores.size());
    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;

    // for all mutations
    for (size_t attachment = 0; attachment < config.getNumAttachments(); ++attachment) {
        // get the attachment scores
        scoreSum = AttachmentScore();
        getAllAttachmentScores(attachmentScores, scoreSum, config, attachment);

        // traverse the attachment scores to the leaves
        // scoreSum = AttachmentScore();
        long double sumParallel = 0;
        PassScoreToChildrenBFSVisitor visBFS(config,
                                             attachmentScores,
                                             attachmentSumScores,
                                             passDownAttachmentScores,
                                             passDownAttachmentSumScores,
                                             sumParallel,
                                             //scoreSum,
                                             attachment);
        breadth_first_search(config.getTree(), num_vertices(config.getTree()) - 1, visitor(visBFS));
        sumParallel -= std::log(2);

        static std::array<double, 6> res;
        res = scoreSum.computeMutTypeContribution(config, scoreSum, true);

        double logPHet = res[0] - res[5];
        double logPHom = res[1] - res[5];
        double logPLossWild = res[2] - res[5];
        double logPLossMut = res[3] - res[5];
        double logPPara = res[4] - res[5];

        // for all attachment positions
        for (unsigned i = 0; i < attachmentScores.size(); ++i)
        {
            if (config.getTree()[i].sample != -1)
            {
                // normalize probabilities
                attachmentSumScores[i].hetScore() -= scoreSum.hetScore();
                attachmentSumScores[i].homScore() -= scoreSum.homScore();
                attachmentSumScores[i].lossWildScore() -= scoreSum.lossWildScore();
                attachmentSumScores[i].lossAltRScore() -= scoreSum.lossAltRScore();
                passDownAttachmentSumScores[i].paralleleScore() -= sumParallel;

                // compute final score
                attachmentSumScores[i].hetScore() = attachmentSumScores[i].hetScore() + logPHet;
                attachmentSumScores[i].finalScore() =  attachmentSumScores[i].hetScore();
                if (config.learnZygocity)
                {
                    attachmentSumScores[i].homScore() = attachmentSumScores[i].homScore() + logPHom;
                    attachmentSumScores[i].finalScore() = addLogProb(attachmentSumScores[i].finalScore(),
                            attachmentSumScores[i].homScore());
                }
                if (config.computeLossScore)
                {
                    attachmentSumScores[i].lossWildScore() = attachmentSumScores[i].lossWildScore() + logPLossWild;
                    attachmentSumScores[i].finalScore() = addLogProb(attachmentSumScores[i].finalScore(),
                            attachmentSumScores[i].lossWildScore());

                    attachmentSumScores[i].lossAltRScore() = attachmentSumScores[i].lossAltRScore() + logPLossMut;
                    attachmentSumScores[i].finalScore() = addLogProb(attachmentSumScores[i].finalScore(),
                            attachmentSumScores[i].lossAltRScore());

                }
                if (config.computeParallelScore)
                {
                    attachmentSumScores[i].lcaRScore() =
                            passDownAttachmentSumScores[i].paralleleScore() + logPPara;
                    attachmentSumScores[i].finalScore() = addLogProb(attachmentSumScores[i].finalScore(),
                            attachmentSumScores[i].lcaRScore());
                }

                // update the counter
                config.mutInSampleCounter[config.getTree()[i].sample][attachment].addInRealSpace(attachmentSumScores[i]);
            }
        }
    }
}

template<typename TTreeType>
double
runMCMC(typename Config<TTreeType>::TGraph &bestTree,
        std::array<std::tuple<double, double>, 9> &bestParams,
        Config<TTreeType> &config,
        std::vector<std::vector<unsigned>> &sampleTrees)
{
    std::vector<std::queue<double>> rootProbabilities(config.getNumMutations());
    return runMCMC(bestTree,bestParams,config,sampleTrees, rootProbabilities);
}

template<typename TTreeType>
double
runMCMC(typename Config<TTreeType>::TGraph &bestTree,
        std::array<std::tuple<double, double>, 9> &bestParams,
        Config<TTreeType> &config,
        std::vector<std::vector<unsigned>> &sampleTrees,
        std::vector<std::queue<double>> &rootProbabilities)
    {

    double bestTreeLogScore = -DBL_MAX;         // initialize best tree score
    config.updateContainers(0);                 // update the container usage (define which mutations to use currently)
    computeNoiseScore(config);                  // compute the noise score
    bestTree = config.getTree();                // init the best tree
    bestParams = config.params;                 // init the best parameter values
    double currTreeLogScore = scoreTree(config);// get score of random tree

    // store the tree
    manageBestTrees(config,
                    bestTreeLogScore,
                    currTreeLogScore,
                    bestTree,
                    bestParams);


    double propTreeLogScore;                    // the new tree log score
    double random;                              // variable for random numbers

    rootProbabilities.resize(config.getNumMutations());
    std::vector<double> sumRootProb(config.getNumMutations(), 0);
    std::vector<bool> mutInRoot(config.getNumMutations());
    

    // these variables are needed for the ml mode
    unsigned numRejectedTreeMoves = 0;
    unsigned maxNumRejectedTreeMoves = 10 * config.numSamples;
    bool acceptAllTreeMoves = false;
    unsigned numAcceptAllTreeMoves = 0;
    unsigned maxNumAcceptAllTreeMoves = config.numSamples;
    unsigned restarts = 0;
    unsigned numMlRestarts = config.ml_mode ? 10 : 0;

    unsigned numMutInRoot = 0;

    // set the number of loops to something very high for the hill climbing mode
    if (config.ml_mode)
    {
        config.loops = 100000000;
    }

    for (int it = 0; it < config.loops + config.sampleLoops; it++) {        // run the iterations of the MCMC


        // This function handles the resizing of the corresponding containers
        if (config.updateContainers(it)) {
            bestTreeLogScore = -DBL_MAX;
            config.getTree() = bestTree;
            currTreeLogScore = scoreTree(config);
            manageBestTrees(config,
                            bestTreeLogScore,
                            currTreeLogScore,
                            bestTree,
                            bestParams);

        }

        // Print some output every 10000 iteration
        if (it % 10000 == 0) {
            config.printParameters();
            std::cout.precision(15);
            std::cout << "iterations: " << it << std::endl;
            std::cout << "score: " << currTreeLogScore << std::endl;
            std::cout << "bestScore: " << bestTreeLogScore << std::endl;
        }

        // Sample the next tree configuration
        proposeNextConfiguration(config);
        // Compute score of new tree
        propTreeLogScore = scoreTree(config);

        // if the maximum likelihood approach is chosen and the current change changed the tree only accept better trees
        random = (double) rand() / RAND_MAX;
        if (config.ml_mode) // && (config.getMoveTyp() != 4))
        {
            random = 1; // only accept better trees
            if (numRejectedTreeMoves == maxNumRejectedTreeMoves)
            {
                acceptAllTreeMoves = true; // once a maximum is hit we restart 10 times. the restart happens from a tree mutated from the current tree
                if (restarts == numMlRestarts) // max number of restarts is reached
                {
                    it = config.loops + config.sampleLoops;
                    config.sampleLoops = 1;
                    random = std::exp((propTreeLogScore - currTreeLogScore)) + 1;
                }
                else if (numAcceptAllTreeMoves < maxNumAcceptAllTreeMoves)
                {
                    random = 0;
                    ++numAcceptAllTreeMoves;
                }
                else if (numAcceptAllTreeMoves == maxNumAcceptAllTreeMoves)
                {
                    restarts += 1;
                    std::cout << "Restart number : " << restarts << std::endl;
                    it = -1; 
                    numAcceptAllTreeMoves = 0;
                    numRejectedTreeMoves = 0;
                    acceptAllTreeMoves = false;
                }

            }
        }

        // the proposed tree is accepted
        if (random < std::exp((propTreeLogScore - currTreeLogScore))) {
            // update counter of parameter estimation
            if (config.getMoveTyp() == 4)
            {
                config.setSDCountParam(config.getParamToOptimize(),
                                       config.getSDCountParam(config.getParamToOptimize()) + 1);
                if (!acceptAllTreeMoves)
                    numRejectedTreeMoves = 0;
            }
            currTreeLogScore = propTreeLogScore;

            manageBestTrees(config,
                            bestTreeLogScore,
                            currTreeLogScore,
                            bestTree,
                            bestParams);
        } 
        else // the proposed tree is rejected
        {
            if (config.getMoveTyp() == 4) {
                config.resetParameters();
            } else {
                config.getTree().swap(config.getTmpTree());

                if (!acceptAllTreeMoves)
                    ++numRejectedTreeMoves;
            }
        }
        
        if (config.learnChi)
        {
            getRootProbabilities(config, rootProbabilities, sumRootProb);
            if (it == 1000)
            {
                for (unsigned i = 0; i < config.getNumMutations(); ++i)
                {
                    mutInRoot[i] = (sumRootProb[i]/1000 >= 0.95) ? true : false;
                    numMutInRoot += mutInRoot[i];
                }
            }
            else if (it > 1000)
            {
                unsigned missed = 0;
                for (unsigned i = 0; i < config.getNumMutations(); ++i)
                {
                    if (mutInRoot[i])
                    {
                        if (sumRootProb[i]/1000 <= 0.95)
                        {
                            missed += 1;
                        }
                    }
                }
                if ((double)missed/(double)(numMutInRoot) > 0.95 || config.parallelScorePenalty < 1.0)
                {
                    config.lossScorePenalty *= 10;
                    config.parallelScorePenalty *= 10;
                    config.learnChi = false;
                    break;
                }
                if (it % 1000 == 0)
                {
                    config.lossScorePenalty /= 10;
                    config.parallelScorePenalty /= 10;
                }
                  
            }

        }


        // sample from the posterior distribution
        if (it >= config.loops) {
            if (config.ml_mode)
            {
                // Because we increased loops to a very large number
                // we will only get here (in the hill climbing mode) if
                // we complete all restarts
                // we therefore now take the global optimum
                config.getTree() = bestTree;
                config.params = bestParams;
                computeLogScoresOP(config);
                scoreTree(config);
            }
            getRootProbabilities(config, rootProbabilities, sumRootProb);
            updateMutInSampleCounts(config);
            config.updateParamsCounter();
            if (config.sampling != 0 && it%config.sampling == 0)
            {
                writeIndex(config, config.samplingName + std::to_string(it) + "/");
            }
        }
    }

    return bestTreeLogScore;
}

// Exponentiate and normalize the mutation counts
inline
void
normalizeMutationCounts(Config<SampleTree> & config)
{
    for (unsigned int j = 0; j < config.getNumMutations(); ++j)
    {
        for (unsigned int i = 0; i < config.getNumSamples(); ++i)
        {
            config.mutInSampleCounter[i][j].exp();
            config.mutInSampleCounter[i][j] /= static_cast<double>(config.sampleLoops);
		}
	}
}

#endif
