/**
 * SCIPhI: Single-cell mutation identification via phylogenetic inference
 * <p>
 * Copyright (C) 2018 ETH Zurich, Jochen Singer
 * <p>
 * This file is part of SCIPhI.
 * <p>
 * SCIPhI is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p>
 * SCIPhI is distributed in the hope that it will be useful,
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

template <typename TTreeType>
void 
printLogScores(Config<TTreeType> & config)
{
	for (unsigned int i = 0; i < config.getLogScores().size(); ++i)
	{
		for (unsigned int j = 0; j < config.getLogScores()[0].size(); ++j)
		{
            std::cout << std::get<0>(config.getLogScores()[i][j]) << " | " << std::get<1>(config.getLogScores()[i][j]) << "\t";
		}
        std::cout << std::endl;
	}
}

// This function keeps the coefficients between 0 and 1
// by mirroring at 0 or 1
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
            //coEff =  1.0 - std::abs(1 - coEff);
            coEff = 2.0 - coEff;
        }
        else
        {
            break;
        }
    }
    return coEff;
}

// This function assures that the overdispersion of the heterozygous mutation
// cannot take a U form (alpha and beta are very close or bigger than 1)
double adjustMutOverDis(double mutationOverDis)
{
    if (mutationOverDis < 2.0)
    {
        return 2.0 + (2.0 - mutationOverDis);
    }
    return mutationOverDis;
}

template <typename TTreeType>
void
changeParameters(Config<TTreeType> & config)
{
    // store the old scores, in case the new parameter will be rejected
    config.setTmpLogScores(config.getLogScores());

    //choose which parameter to optimize
    unsigned numParams = 4 + config.computeLossScore + config.learnZygocity + config.computeParallelScore;
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
    }
    typename Config<TTreeType>::ParamType param = (typename Config<TTreeType>::ParamType)paramIdx;
    config.setParamToOptimize(param);

    /*
    if (config.computeLossScore && config.learnZygocity) 
    {
        config.setParamToOptimize((typename Config<TTreeType>::ParamType)(rand() % 6));
    }
    else if (config.learnZygocity || config.computeLossScore)
    {
        config.setParamToOptimize((typename Config<TTreeType>::ParamType)(rand() % 5));
        if (config.getParamToOptimize() == Config<TTreeType>::nu && config.computeLossScore)
        {
            config.setParamToOptimize(Config<TTreeType>::lambda);
        }
    }
    else
    {
        config.setParamToOptimize((typename Config<TTreeType>::ParamType)(rand() % 4));
    }
    typename Config<TTreeType>::ParamType const & param = config.getParamToOptimize();
     */

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
    
    // draw new paramter value
    std::normal_distribution<double> distribution(oldParamValue, config.getSDParam(param));
    double newParamValue = distribution(config.getGenerator());

    if (paramIdx > 3)
    {
        newParamValue = adjusteCoefficient(newParamValue);
        double sumParams = config.getParam(Config<TTreeType>::ParamType::E_nu)
                + config.getParam(Config<TTreeType>::ParamType::E_lambda)
                + config.getParam(Config<TTreeType>::ParamType::E_parallel);
        if (sumParams > 1)
        {
            for (unsigned i = 4; i < 7; ++i) {
                config.setParam(typename Config<TTreeType>::ParamType(i),
                        config.getParam(typename Config<TTreeType>::ParamType(i)) / sumParams);
            }
        }
    }
    else if (param == Config<TTreeType>::ParamType::E_mutationOverDis)
    {
        newParamValue = adjustMutOverDis(newParamValue);
    }

    config.setParam(param, std::abs(newParamValue));
    if (param == Config<TTreeType>::ParamType::E_wildOverDis || param == Config<TTreeType>::ParamType::E_wildMean)
    {
        computeNoiseScore(config);
    }
    computeLogScoresOP(config);
}

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
                     TTreeType & /*bestTrees*/,
                     double /*minDist*/)
{}

inline
void
manageBestTrees(Config<SampleTree> & config,
                     double & bestTreeLogScore,
                     double currTreeLogScore,
                     std::vector<Config<SampleTree>::TGraph> & bestTrees,
                     std::array<std::tuple<double, double>, 7> & bestParams,
                     double /*tag*/)
{
    if(currTreeLogScore > bestTreeLogScore)
    {   
        std::cout << "The new best score is: " <<  currTreeLogScore << std::endl;
        //save the current state to disk
        writeIndex(config, config.bestName);
        bestTreeLogScore = currTreeLogScore;
        bestTrees[0] = config.getTree();
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
updateMutInSampleCounts(Config<SampleTree> & config, unsigned it)
{
    Config<SampleTree>::TAttachmentScores & attachmentScores = config.getTmpAttachmentScore();
    static Config<SampleTree>::TAttachmentScores attachmentSumScores;
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentScores;
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentSumScores;
    attachmentSumScores.resize(attachmentScores.size());
    passDownAttachmentScores.resize(attachmentScores.size());
    passDownAttachmentSumScores.resize(attachmentScores.size());
    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;


    {
        std::ofstream ofs(config.outFilePrefix +  "_r_" + std::to_string(it) + ".gv");
        write_graphviz(ofs, config.getTree(), my_label_writer_complete(config.getTree()));
        ofs.close();
    }

    {
        typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, Vertex<SimpleTree>> TGraph;
        TGraph newTreeBest = simplifyTree(config);
        std::ofstream ofs2(config.outFilePrefix + "_" + std::to_string(it) + ".gv");
        write_graphviz(ofs2, newTreeBest,
                       my_label_writer(newTreeBest, config.indexToPosition, config.cellNames, config.cellColours,
                                       config.cellClusters));
    }

    {
        writeIndex(config);
    }

    // for all mutations
    for (size_t attachment = 0; attachment < config.getNumAttachments(); ++attachment) {
        // get the attachment scores
        getAllAttachmentScores(attachmentScores, scoreSum, config, attachment);

        // traverse the attachment scores to the leaves
        scoreSum = AttachmentScore();
        PassScoreToChildrenBFSVisitor visBFS(config,
                                             attachmentScores,
                                             attachmentSumScores,
                                             passDownAttachmentScores,
                                             passDownAttachmentSumScores,
                                             scoreSum,
                                             attachment);
        breadth_first_search(config.getTree(), num_vertices(config.getTree()) - 1, visitor(visBFS));

        static std::array<double, 5> res;
        scoreSum.computeMutTypeContribution(config, scoreSum, res);

        double logPHet = res[0] - res[4];
        double logPHom = res[1] - res[4];
        double logPLoss = res[2] - res[4];
        double logPPara = res[3] - res[4];

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
                passDownAttachmentSumScores[i].paralleleScore() -= scoreSum.lcaRScore();

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
                    // TODO: overwriting lossWildScore in this way is error prone
                    attachmentSumScores[i].lossWildScore() =
                            addLogProbWeight(attachmentSumScores[i].lossWildScore(),
                                    attachmentSumScores[i].lossAltRScore(), 0.5) + logPLoss;
                    attachmentSumScores[i].finalScore() = addLogProb(attachmentSumScores[i].finalScore(),
                            attachmentSumScores[i].lossWildScore());
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

template <typename TTreeType>
double
runMCMC(std::vector<typename Config<TTreeType>::TGraph> & bestTrees, 
               std::array<std::tuple<double, double>, 7> & bestParams,
			   Config<TTreeType> & config,
			   std::vector<std::vector<unsigned>> & sampleTrees)
{


	int minDist = INT_MAX;                       // initialize distance to true tree if given
	double bestTreeLogScore = -DBL_MAX;          // initialize best tree score

	for(unsigned r=0; r<config.reps; r++){   // repeat the MCMC, start over with random tree each time, only best score and list of best trees is kept between repetitions
        
        config.updateContainers(0);
        computeNoiseScore(config);
        
        bestTrees.resize(1);
        bestTrees[0] = config.getTree();
        bestParams = config.params;

		double currTreeLogScore = scoreTree(config);       // get score of random tree

        manageBestTrees(config,
                bestTreeLogScore, 
                currTreeLogScore,
                bestTrees,
                bestParams,
                minDist);


        double propTreeLogScore;
        double random;

        for(unsigned it=0; it<config.loops + config.sampleLoops; it++){        // run the iterations of the MCMC

            // Mutations are added in batches of 10% of all mutations 10 times
            // This function handels the resizing of the corresponding containers
            if (config.updateContainers(it))
            {
                minDist = INT_MAX;                       // initialize distance to true tree if given
                bestTreeLogScore = -DBL_MAX;          // initialize best tree score
                config.getTree() = bestTrees[0];
		        currTreeLogScore = scoreTree(config);       // get score of random tree
                manageBestTrees(config,
                        bestTreeLogScore, 
                        currTreeLogScore,
                        bestTrees,
                        bestParams,
                        minDist);
            }
        	
            if(it % 10000 == 0)
			{ 
                config.printParameters();
                std::cout.precision(15);
                std::cout << "iterations: " << it << std::endl;
                std::cout << "score: " << currTreeLogScore << std::endl;
                std::cout << "bestScore: " << bestTreeLogScore << std::endl;
			}
            
            // sample the next tree configuration
            proposeNextConfiguration(config); 
            propTreeLogScore = scoreTree(config);


            random = (double) rand() / RAND_MAX;
            
            // the proposed tree is accepted
            if (random < std::exp((propTreeLogScore-currTreeLogScore)))
            {
                // update counter of parameter estimation
                if (config.getMoveTyp() == 4)
                    config.setSDCountParam(config.getParamToOptimize(), config.getSDCountParam(config.getParamToOptimize()) + 1);
                currTreeLogScore  = propTreeLogScore;

                manageBestTrees(config,
                        bestTreeLogScore, 
                        currTreeLogScore,
                        bestTrees,
                        bestParams,
                        minDist);
            }
            else // the proposed tree is rejected
            {
                if (config.getMoveTyp() == 4)
                {
                    config.resetParameters();
                }
                else
                {
                    config.getTree().swap(config.getTmpTree());
                }
            }
           
            // sample from the posterior distribution
            if (it >= config.loops)
            {
                updateMutInSampleCounts(config,it);
                config.updateParamsCounter();
                if (config.sampling != 0 && it%config.sampling == 0)
                {
                    writeIndex(config, config.samplingName + std::to_string(it) + "/");
                }
            }
        }
	} 
	
 	return bestTreeLogScore;
}

inline
void 
normalizeMutationCounts(Config<SampleTree> & config)
{
    for (unsigned int j = 0; j < config.getNumMutations(); ++j)
    {
        for (unsigned int i = 0; i < config.getNumSamples(); ++i)
        {
            config.mutInSampleCounter[i][j].exp();
			if (j < config.getNumMutations() - config.numUniqMuts)
            {
                config.mutInSampleCounter[i][j] /= static_cast<double>(config.sampleLoops);
            }
		}
	}
}

#endif
