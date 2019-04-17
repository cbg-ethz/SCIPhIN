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
#ifndef PROBABILITIES_H
#define PROBABILITIES_H

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
#include "rand.h"
#include "sciphi_config.h"
#include "trees.h"

#include <dlib/global_optimization.h>

#include <boost/math/special_functions/digamma.hpp>

// The computation of the beta-binomial function can be divided into three parts
// The three functions below each compute one of the three parts
inline
double
logBetaBinCountsTerm(double sup, double cov) {
    return std::lgamma(cov + 1.0)
           - std::lgamma(sup + 1.0)
           - std::lgamma(cov - sup + 1.0);
}

inline
double
logBetaBinMixedTerm(double sup, double cov, double mean, double overDis) {
    return std::lgamma(sup + mean * overDis)
           + std::lgamma(cov - sup + overDis * (1.0 - mean))
           - std::lgamma(cov + overDis);
}

inline
double
logBetaBinParamsTerm(double mean, double overDis) {
    return std::lgamma(overDis)
           - std::lgamma(mean * overDis)
           - std::lgamma(overDis * (1.0 - mean));
}

// Combine the three terms computed above
inline
double
logBetaBinPDF(double sup, double cov, double mean, double overDis) {
    if (cov == 0)
        return 0;

    return logBetaBinCountsTerm(sup, cov) +
           logBetaBinMixedTerm(sup, cov, mean, overDis) +
           logBetaBinParamsTerm(mean, overDis);
}

// Compute the tree terms above, however, the first term is actually not needed as it is basically a constant between
// iterations
inline
double
logBetaBinPDFOP(double sup, double cov, double mean, double overDis) {
    if (cov == 0)
        return 0;

    return //logBetaBinCountsTerm(sup, cov) + // we do not need the counts term as we are only interested in the score between iterations
            logBetaBinMixedTerm(sup, cov, mean, overDis) +
            logBetaBinParamsTerm(mean, overDis);
}

// Add two values in real space by first exponentiating
inline
double
addLogProb(double x, double y) {
    if (x == -INFINITY) {
        return y;
    }

    if (y == -INFINITY) {
        return x;
    }

    double maxScore;
    double minScore;
    if (x > y) {
        maxScore = x;
        minScore = y;
    } else {
        maxScore = y;
        minScore = x;
    }
    return std::log(1.0 + std::exp(minScore - maxScore)) + maxScore;
}

// Subtract two values in real space by first exponentiating
inline
double
subLogProb(double x, double y) {
    if (y == -INFINITY) // it is important to first check this case
    {                   // because it avoids - -INFINITY in the next case
        return x;
    }

    if (std::abs(x - y) <= 100 * std::numeric_limits<double>::epsilon()) {
        return -INFINITY;
    }

    double maxScore;
    double minScore;
    if (x > y) {
        maxScore = x;
        minScore = y;
    } else {
        maxScore = y;
        minScore = x;
    }
    return log(1.0 - exp(minScore - maxScore)) + maxScore;
}

// Add two weighted values in real space by first exponentiating
inline
double
addLogProbWeight(double x, double y, double nu) // = (1.0 - nu) * x + nu * y
{
    if (nu == 0.0) {
        return x;
    }

    if (nu == 1.0) {
        return y;
    }

    if (x == -INFINITY) {
        return y;
    }

    if (y == -INFINITY) {
        return x;
    }

    if (x > y) {
        return std::log((1.0 - nu) + nu * std::exp(y - x)) + x;
    }
    return std::log(nu + (1.0 - nu) * std::exp(x - y)) + y;
}

// Compute the actual plain wild type beta-binomial function
template<typename TTreeType>
double
computeRawWildLogScore(Config<TTreeType> const &config, double altCount, double coverage) {
    return logBetaBinPDF(altCount,
                         coverage,
                         config.getParam(Config<TTreeType>::E_wildMean),
                         config.getParam(Config<TTreeType>::E_wildOverDis));
}

// Compute the actual plain wild type beta-binomial function using the optimization above
template<typename TTreeType>
double
computeRawWildLogScoreOP(Config<TTreeType> const &config, double altCount, double coverage) {
    return logBetaBinPDFOP(altCount,
                           coverage,
                           config.getParam(Config<TTreeType>::E_wildMean),
                           config.getParam(Config<TTreeType>::E_wildOverDis));
}

// Compute the actual plain mutation type beta-binomial function
template<typename TTreeType>
double
computeRawMutLogScore(Config<TTreeType> const &config, double altCount, double coverage) {
    return logBetaBinPDF(altCount,
                         coverage,
                         config.getParam(Config<TTreeType>::E_mutationMean),
                         config.getParam(Config<TTreeType>::E_mutationOverDis));
}

// Compute the actual plain mutation type beta-binomial function using the optimization above
template<typename TTreeType>
double
computeRawMutLogScoreOP(Config<TTreeType> &config, double altCount, double coverage) {
    return logBetaBinPDFOP(altCount,
                           coverage,
                           config.getParam(Config<TTreeType>::E_mutationMean),
                           config.getParam(Config<TTreeType>::E_mutationOverDis));
}

template<typename TTreeType>
double
computeWildLogScoreOP(Config<TTreeType> &config, double altCount, double coverage) {
    double logWild = computeRawWildLogScoreOP(config, altCount, coverage);
    double logMut = computeRawMutLogScoreOP(config, altCount, coverage);
    return addLogProbWeight(logWild, logMut, config.sub);
}


template<typename TTreeType>
double
computeHomoLogScoreOP(Config<TTreeType> &config, double altCount, double coverage) {
    return logBetaBinPDFOP(coverage - altCount,
                           coverage,
                           config.getParam(Config<TTreeType>::E_wildMean),
                           config.getParam(Config<TTreeType>::E_wildOverDis));
}

// This function combines the wild type beta-bin, the hetero beta-bin and the homo beta-bin
template<typename TTreeType>
double
computeHeteroLogScoreOP(Config<TTreeType> &config, double altCount, double coverage) {
    double mu = config.getParam(Config<TTreeType>::E_mu);
    double oneMinusMu = 1.0 - mu;
    double logHomo = addLogProbWeight(computeRawWildLogScoreOP(config, altCount, coverage),
                                      computeHomoLogScoreOP(config, altCount, coverage), 0.5);
    return addLogProbWeight(computeRawMutLogScoreOP(config, altCount, coverage), logHomo, oneMinusMu);
}

template<typename TTreeType>
void computeWildLogScoresOP(Config<TTreeType> &config) {
    for (unsigned int i = 0; i < config.getLogScores().size(); ++i) {
        for (unsigned int j = 0; j < config.getLogScores()[0].size(); ++j) {
            std::get<0>(config.getLogScores()[i][j]) = computeWildLogScoreOP(config,
                                                                             std::get<1>(config.getData()[i][j]),
                                                                             std::get<0>(config.getData()[i][j]));
        }
    }
}

template<typename TTreeType>
void computeLogScoresOP(Config<TTreeType> &config) {
    for (unsigned int i = 0; i < config.getLogScores().numCells(); ++i) {
        for (unsigned int j = 0; j < config.getLogScores().numMuts(); ++j) {
            config.getLogScores().wtScore(i, j) = computeWildLogScoreOP(config, std::get<1>(config.getData()[i][j]),
                                                                        std::get<0>(config.getData()[i][j]));
            config.getLogScores().hetScore(i, j) = computeHeteroLogScoreOP(config, std::get<1>(config.getData()[i][j]),
                                                                           std::get<0>(config.getData()[i][j]));
            config.getLogScores().homScore(i, j) = computeHomoLogScoreOP(config, std::get<1>(config.getData()[i][j]),
                                                                         std::get<0>(config.getData()[i][j]));
        }
    }
}

inline
void
computeNoiseScore(Config<SampleTree> &config) {

    long double mean = config.getParam(Config<SampleTree>::E_wildMean);
    long double overDis = config.getParam(Config<SampleTree>::E_wildOverDis);
    long double alpha = mean * overDis;
    long double beta = overDis * (1.0 - mean);

    config.noiseScore = config.noiseCounts.numPos * logBetaBinParamsTerm(mean, overDis);
    for (unsigned i = 0; i < config.noiseCounts.sup.size(); ++i) {
        config.noiseScore += config.noiseCounts.sup[i].second * std::lgamma(config.noiseCounts.sup[i].first + alpha);
    }
    for (unsigned i = 0; i < config.noiseCounts.covMinusSup.size(); ++i) {
        config.noiseScore +=
                config.noiseCounts.covMinusSup[i].second * std::lgamma(config.noiseCounts.covMinusSup[i].first + beta);
    }
    for (unsigned i = 0; i < config.noiseCounts.cov.size(); ++i) {
        config.noiseScore -= config.noiseCounts.cov[i].second * std::lgamma(config.noiseCounts.cov[i].first + overDis);
    }
}

// This class is used to pass the score from the root towards the leaves
class PassScoreToChildrenBFSVisitor : public boost::default_bfs_visitor {
    Config<SampleTree> &config;
    Config<SampleTree>::TAttachmentScores &attachmentScores;
    Config<SampleTree>::TAttachmentScores &attachmentSumScores;
    Config<SampleTree>::TPassDownAttachmentScores &passDownAttachmentScores;
    Config<SampleTree>::TPassDownAttachmentScores &passDownAttachmentSumScores;
    Config<SampleTree>::TAttachmentScores::TAttachmentScore &sumScore;

public:

    PassScoreToChildrenBFSVisitor(Config<SampleTree> &config_,
                                  Config<SampleTree>::TAttachmentScores &attachmentScores_,
                                  Config<SampleTree>::TAttachmentScores &attachmentSumScores_,
                                  Config<SampleTree>::TPassDownAttachmentScores &passDownAttachmentScores_,
                                  Config<SampleTree>::TPassDownAttachmentScores &passDownAttachmentSumScores_,
                                  Config<SampleTree>::TAttachmentScores::TAttachmentScore &sumScore_,
                                  unsigned attachment_) :
            config(config_),
            attachmentScores(attachmentScores_),
            attachmentSumScores(attachmentSumScores_),
            passDownAttachmentScores(passDownAttachmentScores_),
            passDownAttachmentSumScores(passDownAttachmentSumScores_),
            sumScore(sumScore_) {}

    template<typename TVertex>
    void discover_vertex(TVertex v, typename Config<SampleTree>::TGraph const &g) const {
        (void) g;
        if (v == num_vertices(g) - 1) {
            return;
        }

        attachmentSumScores[v] = AttachmentScore();
        passDownAttachmentScores[v] = PassDownAttachmentScore();
        passDownAttachmentSumScores[v] = PassDownAttachmentScore();

        unsigned pN = source(*in_edges(v, g).first, g); // parent Node

        if (pN == num_vertices(config.getTree()) - 1) {
            sumScore.addInRealSpace(attachmentScores[v]);
            attachmentSumScores[v] = attachmentScores[v];
            return;
        }

        if (g[v].sample == -1) {
            auto it = out_edges(v, g).first;
            unsigned lN = target(*it, g);                   // left node id
            unsigned rN = target(*(it + 1), g);             // right node id
            unsigned sN = getSibling(g, v);                  // sibling node id

            passDownAttachmentScores.computeLogLossInCurrentInnerNode(g, attachmentScores, v);

            sumScore.addInRealSpace(attachmentScores[v]);

            attachmentSumScores[v].hetScore() =
                    addLogProb(attachmentScores[v].hetScore(),
                               attachmentSumScores[pN].hetScore());
            attachmentSumScores[v].homScore() =
                    addLogProb(attachmentScores[v].homScore(),
                               attachmentSumScores[pN].homScore());
            attachmentSumScores[v].lossAltRScore() =
                    addLogProb(attachmentScores[v].lossAltRScore(),
                               attachmentSumScores[pN].lossAltRScore());
            passDownAttachmentSumScores[v].lossAltInCurrentNodeRScore() =
                    addLogProb(passDownAttachmentScores[v].lossAltInCurrentNodeRScore(),
                               passDownAttachmentSumScores[pN].lossAltInCurrentNodeRScore());
            attachmentSumScores[v].lossWildScore() =
                    addLogProb(attachmentScores[v].lossWildScore(),
                               attachmentSumScores[pN].lossWildScore());

            unsigned ppN = source(*in_edges(pN, g).first, g);
            passDownAttachmentSumScores[v].sibNodeScore() = attachmentScores[getSibling(g, v)].hetSumScore();
            if (ppN != num_vertices(config.getTree()) - 1) {
                passDownAttachmentSumScores[v].sibNodeScore() = addLogProb(
                        passDownAttachmentSumScores[v].sibNodeScore(),
                        passDownAttachmentSumScores[pN].sibNodeScore());
            }
            passDownAttachmentScores[v].paralleleScore() =
                    attachmentScores[v].hetScore() + passDownAttachmentSumScores[v].sibNodeScore();
            if (g[sN].sample == -1) {
                passDownAttachmentScores[v].paralleleScore() = subLogProb(
                        passDownAttachmentScores[v].paralleleScore(),
                        attachmentScores[v].hetScore() +
                        addLogProb(attachmentScores[sN].hetScore(), attachmentScores[sN].childHetSumScore()));
            }
            if (ppN != num_vertices(config.getTree()) - 1) {
                if (g[getSibling(g, (TVertex) pN)].sample == -1) {
                    passDownAttachmentScores[v].paralleleScore() = subLogProb(
                            passDownAttachmentScores[v].paralleleScore(),
                            attachmentScores[v].hetScore() + attachmentScores[getSibling(g, (TVertex) pN)].hetScore());
                }
            }
            passDownAttachmentSumScores[v].paralleleScore() = addLogProb(
                    passDownAttachmentScores[v].paralleleScore(),
                    passDownAttachmentSumScores[pN].paralleleScore());
        } else {
            sumScore.hetScore() = addLogProb(sumScore.hetScore(), attachmentScores[v].hetScore());

            attachmentSumScores[v].hetScore() =
                    addLogProb(attachmentScores[v].hetScore(),
                               attachmentSumScores[pN].hetScore());
            attachmentSumScores[v].homScore() = attachmentSumScores[pN].homScore();
            attachmentSumScores[v].lossAltRScore() = subLogProb(
                    attachmentSumScores[pN].lossAltRScore(),
                    passDownAttachmentSumScores[pN].lossAltInCurrentNodeRScore());
            attachmentSumScores[v].lossWildScore() = attachmentSumScores[pN].lossWildScore();
            passDownAttachmentSumScores[v].paralleleScore() = passDownAttachmentSumScores[pN].paralleleScore();
        }
    }
};

/*
 * This function is used to optimize mean and overdispersion
 * for a given locus
 */
struct OptimizeBetaBinMeanOverDis {
    std::vector<std::pair<unsigned, unsigned>> const &counts;

    OptimizeBetaBinMeanOverDis(
            std::vector<std::pair<unsigned, unsigned>> const &counts_) :
            counts(counts_) {};

    double operator()(const dlib::matrix<double, 0, 1> &x) const {
        double result = 0;
        for (size_t cell = 0; cell < this->counts.size(); ++cell) {
            result += logBetaBinPDF(this->counts[cell].first, this->counts[cell].second, x(0), x(1));
        }
        return result;
    };
};

struct OptimizeBetaBinMeanOverDisDerivates {
    std::vector<std::pair<unsigned, unsigned>> const &counts;

    OptimizeBetaBinMeanOverDisDerivates(
            std::vector<std::pair<unsigned, unsigned>> const &counts_) :
            counts(counts_) {}

    dlib::matrix<double> operator()(const dlib::matrix<double, 0, 1> &x) const {
        double mean = x(0);
        double overDis = x(1);
        dlib::matrix<double, 0, 1> res = {0, 0};

        double temp = 0;
        unsigned counter = 0;
        for (size_t cell = 0; cell < this->counts.size(); ++cell) {
            unsigned k = this->counts[cell].first;
            unsigned n = this->counts[cell].second;
            temp += overDis * boost::math::digamma(k + mean * overDis)
                    - overDis * boost::math::digamma(n - k + overDis - overDis * mean);
            ++counter;
        }
        res(0) = counter * (-overDis * boost::math::digamma(mean * overDis) +
                            overDis * boost::math::digamma(overDis - overDis * mean)) + temp;

        temp = 0;
        for (size_t cell = 0; cell < this->counts.size(); ++cell) {
            unsigned k = this->counts[cell].first;
            unsigned n = this->counts[cell].second;
            temp += mean * boost::math::digamma(k + mean * overDis) +
                    (1.0 - mean) * boost::math::digamma(n - k + overDis - overDis * mean) -
                    boost::math::digamma(n + overDis);
        }
        res(1) = counter * (boost::math::digamma(overDis) - mean * boost::math::digamma(mean * overDis) -
                            (1.0 - mean) * boost::math::digamma(overDis - overDis * mean)) + temp;

        return res;
    }
};

struct OptimizeBetaBinOverDis {
    std::vector<std::pair<unsigned, unsigned>> const &counts;
    double meanFilter;

    OptimizeBetaBinOverDis(
            std::vector<std::pair<unsigned, unsigned>> const &counts_,
            double meanFilter_) :
            counts(counts_),
            meanFilter(meanFilter_) {};

    double operator()(const dlib::matrix<double, 0, 1> &x) const {
        double result = 0;
        for (size_t cell = 0; cell < this->counts.size(); ++cell) {
            result += logBetaBinPDF(this->counts[cell].first, this->counts[cell].second, meanFilter, x(0));
        }
        return result;
    };
};

struct OptimizeBetaBinOverDisDerivates {
    std::vector<std::pair<unsigned, unsigned>> const &counts;
    double meanFilter;

    OptimizeBetaBinOverDisDerivates(
            std::vector<std::pair<unsigned, unsigned>> const &counts_,
            double meanFilter_) :
            counts(counts_),
            meanFilter(meanFilter_) {};

    dlib::matrix<double> operator()(const dlib::matrix<double, 0, 1> &x) const {
        double mean = this->meanFilter;
        double overDis = x(0);
        dlib::matrix<double, 0, 1> res = {0};

        double temp = 0;
        unsigned counter = 0;
        for (size_t cell = 0; cell < this->counts.size(); ++cell) {
            unsigned k = this->counts[cell].first;
            unsigned n = this->counts[cell].second;
            temp += mean * boost::math::digamma(k + mean * overDis) +
                    (1.0 - mean) * boost::math::digamma(n - k + overDis - overDis * mean) -
                    boost::math::digamma(n + overDis);
            ++counter;
        }
        res(0) = counter * (boost::math::digamma(overDis) - mean * boost::math::digamma(mean * overDis) -
                            (1.0 - mean) * boost::math::digamma(overDis - overDis * mean)) + temp;

        return res;
    }
};

#endif
