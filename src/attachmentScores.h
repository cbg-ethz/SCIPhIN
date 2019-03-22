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

#ifndef ATTACHMENTSCORES_H
#define ATTACHMENTSCORES_H

#include <limits.h>

#include "sciphi_config.h"

// Forward declaration
double addLogProbWeight(double x, double y, double nu);

double addLogProb(double x, double y);

double subLogProb(double x, double y);

double addLog_nan_x(double result, double rlNScore);

/*
 * This class stores values associated with the attachment of mutations
 * to nodes in the phylogenetic tree.
 */
struct AttachmentScore {

    typedef std::array<double, 16> TAttachmentScore;

    TAttachmentScore attachmentScore;


    enum Type {
        E_hetScore = 0,                     // all cells below the node are heterozygous
        E_homScore = 1,                     // all cells below the node are homozygous
        E_hetSumScore = 2,                  // the sum of all hetScores below the node
        E_lossAltScore = 3,                 // a mutation is gained in the node and the mutated/altered allele is lost in any
                                            // node below
        E_lossAltRScore = 4,                // a mutation is gained in the node and the mutated/altered allele is lost in any
                                            // node below except for the child nodes
        E_lossWildScore = 5,                // a mutation is gained in the node and the reference allele is lost in any node below
                                            // except for the child nodes
        E_lcaScore = 6,                     // probability that the node is the lowest common ancestor of two independent mutation
                                            // events
        E_childHetSumScore = 7,             // the probability that one of the child nodes is mutated
        E_lcaRScore = 8,                   // probability that the node is the lowest common ancestor of two independent mutation
                                            // events with additional restrictions
        E_grantChildHetSumScore = 9,       // probability of a mutation below the current node, but not in itself
                                            // or the child nodes
        E_numAltPoss = 10,                  // number of different possibilities to loose the current mutation
        E_numAltRPoss = 11,                 // number of different possibilities to loose the current mutation (but not directly in
                                            // the child nodes
        E_numInnerNodes = 12,               // number of inner nodes below the current one
        E_numInnerChildNodes = 13,          // number of child nodes that are inner nodes (0, 1, or 2)
        E_numLcaRPoss = 14,                 // number of possibilities for a parallel mutation below the current node
        E_finalScore = 15                   // final combined score

    };

    // init everything to -INFINITY because we are doing the computation in 
    // log space and exp(-INFINITY) = 0
    AttachmentScore() :
            attachmentScore({{-INFINITY,
                                     -INFINITY,
                                     -INFINITY,
                                     -INFINITY,
                                     -INFINITY,
                                     -INFINITY,
                                     -INFINITY,
                                     -INFINITY,
                                     -INFINITY,
                                     -INFINITY,
                                     0,
                                     0,
                                     0,
                                     0,
                                     0,
                                     -INFINITY}}) {};

    // transform attachment scores into log
    void log() {
        for (unsigned i = 0; i < this->attachmentScore.size(); ++i) {
            this->attachmentScore[i] = std::log(this->attachmentScore[i]);
        }
    }

    // exponentiate attachment scores
    void exp() {
        for (unsigned i = 0; i < this->attachmentScore.size(); ++i) {
            this->attachmentScore[i] = std::exp(this->attachmentScore[i]);
        }
    }

    double &hetScore() {
        return this->attachmentScore[E_hetScore];
    }

    double const &hetScore() const {
        return this->attachmentScore[E_hetScore];
    }

    double &hetSumScore() {
        return this->attachmentScore[E_hetSumScore];
    }

    double const &hetSumScore() const {
        return this->attachmentScore[E_hetSumScore];
    }

    double &homScore() {
        return this->attachmentScore[E_homScore];
    }

    double const &homScore() const {
        return this->attachmentScore[E_homScore];
    }

    double &lossWildScore() {
        return this->attachmentScore[E_lossWildScore];
    }

    double const &lossWildScore() const {
        return this->attachmentScore[E_lossWildScore];
    }

    double &lossAltScore() {
        return this->attachmentScore[E_lossAltScore];
    }

    double const &lossAltScore() const {
        return this->attachmentScore[E_lossAltScore];
    }

    double &lossAltRScore() {
        return this->attachmentScore[E_lossAltRScore];
    }

    double const &lossAltRScore() const {
        return this->attachmentScore[E_lossAltRScore];
    }

    double &lcaScore() {
        return this->attachmentScore[E_lcaScore];
    }

    double const &lcaScore() const {
        return this->attachmentScore[E_lcaScore];
    }

    double const &lcaRScore() const {
        return this->attachmentScore[E_lcaRScore];
    }

    double &lcaRScore() {
        return const_cast<double &>(static_cast<const AttachmentScore*>(this)->lcaRScore());
    }

    double const &grantChildHetSumScore() const {
        return this->attachmentScore[E_grantChildHetSumScore];
    }

    double &grantChildHetSumScore() {
        return const_cast<double &>(static_cast<const AttachmentScore*>(this)->grantChildHetSumScore());
    }

    double &childHetSumScore() {
        return this->attachmentScore[E_childHetSumScore];
    }

    double const &childHetSumScore() const {
        return this->attachmentScore[E_childHetSumScore];
    }

    double const &numAltPoss() const {
        return this->attachmentScore[E_numAltPoss];
    }

    double &numAltPoss() {
        return const_cast<double &>(static_cast<const AttachmentScore*>(this)->numAltPoss());
    }

    double const &numAltRPoss() const {
        return this->attachmentScore[E_numAltRPoss];
    }

    double &numAltRPoss() {
        return const_cast<double &>(static_cast<const AttachmentScore*>(this)->numAltRPoss());
    }

    double const &numInnerNodes() const {
        return this->attachmentScore[E_numInnerNodes];
    }

    double &numInnerNodes() {
        return const_cast<double &>(static_cast<const AttachmentScore*>(this)->numInnerNodes());
    }

    double const &numInnerChildNodes() const {
        return this->attachmentScore[E_numInnerChildNodes];
    }

    double &numInnerChildNodes() {
        return const_cast<double &>(static_cast<const AttachmentScore*>(this)->numInnerChildNodes());
    }

    double const &numLcaRPoss() const {
        return this->attachmentScore[E_numLcaRPoss];
    }

    double &numLcaRPoss() {
        return const_cast<double &>(static_cast<const AttachmentScore*>(this)->numLcaRPoss());
    }

    double &finalScore() {
        return this->attachmentScore[E_finalScore];
    }

    double const &finalScore() const {
        return this->attachmentScore[E_finalScore];
    }

    /*AttachmentScore &operator+=(AttachmentScore const &rightSide) {
        this->hetScore() += rightSide.hetScore();
        this->homScore() += rightSide.homScore();
        // this check is necessary as there is not always a loss score available to be added
        if (!isnan(rightSide.lossWildScore())) {
            this->lossWildScore() += rightSide.lossWildScore();
        }
        if (!isnan(rightSide.lossAltScore())) {
            this->lossAltScore() += rightSide.lossAltScore();
        }
        this->finalScore() += rightSide.finalScore();

        return *this;
    }

    AttachmentScore &operator+=(AttachmentScore &rightSide) {
        this->hetScore() += rightSide.hetScore();
        this->homScore() += rightSide.homScore();

        // this check is necessary as there is not always a loss score available to be added
        if (!isnan(rightSide.lossWildScore())) {
            this->lossWildScore() += rightSide.lossWildScore();
        }
        if (!isnan(rightSide.lossAltScore())) {
            this->lossAltScore() += rightSide.lossAltScore();
        }
        this->finalScore() += rightSide.finalScore();

        return *this;
    }*/

    AttachmentScore &operator-=(AttachmentScore &rightSide) {
        for (unsigned i =0; i < this->attachmentScore.size(); ++i)
        {
            this->attachmentScore[i] -= rightSide.attachmentScore[i];
        }

        /*
        this->hetScore() -= rightSide.hetScore();
        this->homScore() -= rightSide.homScore();
        // this check is necessary as there is not always a loss score available to be added
        if (!isnan(rightSide.lossWildScore())) {
            this->lossWildScore() -= rightSide.lossWildScore();
        }
        if (!isnan(rightSide.lossAltScore())) {
            this->lossAltScore() -= rightSide.lossAltScore();
        }
        this->finalScore() -= rightSide.finalScore();
         */

        return *this;
    }

    // this function takes two scores in log space, exponentiates them and 
    // returns the log of the sum

    void addInRealSpace(AttachmentScore const &rightSide) {
        for (unsigned i = 0; i < this->attachmentScore.size(); ++i) {
            this->attachmentScore[i] = addLogProb(this->attachmentScore[i], rightSide.attachmentScore[i]);
        }

        /*
        this->hetScore() = addLogProb(this->hetScore(), rightSide.hetScore());
        this->homScore() = addLogProb(this->homScore(), rightSide.homScore());
        if (!isnan(rightSide.lossWildScore())) {
            this->lossWildScore() = addLogProb(this->lossWildScore(), rightSide.lossWildScore());
        }
        if (!isnan(rightSide.lossAltScore())) {
            this->lossAltScore() = addLogProb(this->lossAltScore(), rightSide.lossAltScore());
        }
        this->finalScore() = addLogProb(this->finalScore(), rightSide.finalScore());*/
    }

    AttachmentScore &operator/=(double rightSide) {
        for (unsigned i = 0; i < 8; ++i) {
            this->attachmentScore[i] /= rightSide;
        }
        return *this;
    }

    // this functions computes the final score of an attachment point
    void computeFinalScore(
            double nu,
            double lambda,
            unsigned numNodes,
            unsigned numMutPlacements,
            bool isLeaf,
            bool useLoss,
            bool useParallel) {
        // if the node is a leaf just return the weighted hetero score
        if (isLeaf) {
            this->finalScore() = this->hetScore() - std::log(numNodes * 2 + 1);
            return;
        }

        if (!useLoss) { // || isnan(this->lossWildScore())) {
            this->finalScore() = addLogProbWeight(this->hetScore() - std::log(numNodes * 2 + 1),
                                                  this->homScore() - std::log(numNodes), nu);
            return;
        }

        if (!useParallel) {
            double loss = addLogProbWeight(this->lossWildScore(), this->lossAltScore(), 0.5);
            this->finalScore() = std::log((1.0 - lambda - nu) / (numNodes * 2 + 1) * std::exp(this->hetScore()) +
                                          nu / numNodes * std::exp(this->homScore()) +
                                          lambda / numMutPlacements * std::exp(loss));
        }

        return;
    }

    void setMinusInfinity() {
        for (unsigned i = 0; i < this->attachmentScore.size(); ++i) {
            this->attachmentScore[i] = -INFINITY;
        }
    }
};

template<unsigned N>
std::ostream &operator<<(std::ostream &os, AttachmentScore const &obj) {
    os << obj.hetScore() << "|" << obj.homScore();
    if (N == 3) {
        return os << "|" << obj.finalScore();
    }
    os << "|" << obj.lossWildScore() << "|" << obj.lossAltScore() << "|" << obj.finalScore();
    return os;
}

std::ostream &operator<<(std::ostream &os, AttachmentScore const &obj) {
    os << "hetScore:         " << obj.hetScore() << "\n";
    os << "homScore:         " << obj.homScore() << "\n";
    os << "hetSumScore:      " << obj.hetSumScore() << "\n";
    os << "lossAltScore:     " << obj.lossAltScore() << "\n";
    os << "lossAltRScore:    " << obj.lossAltRScore() << "\n";
    os << "lossWildScore:    " << obj.lossWildScore() << "\n";
    os << "lcaScore:         " << obj.lcaScore() << "\n";
    os << "childHetSumScore: " << obj.childHetSumScore() << "\n";
    os << "lcaRScore:        " << obj.lcaRScore() << "\n";
    os << "numAltRPoss:      " << obj.numAltRPoss() << "\n";
    os << "finalScore:       " << obj.finalScore() << "\n";
    return os;
}

AttachmentScore operator+(AttachmentScore &leftSide, AttachmentScore &rightSide) {
    AttachmentScore result;
    for (unsigned i = 0; i < 8; ++i) {
        result.attachmentScore[i] = leftSide.attachmentScore[i] + rightSide.attachmentScore[i];
    }
    return result;
}

/*
 * This class combines many AttachmentScores into one object
 */
struct AttachmentScores {
    typedef AttachmentScore TAttachmentScore;
    typedef std::vector<AttachmentScore> TAttachmentScores;

    TAttachmentScores attachmentScores;

    AttachmentScore &operator[](unsigned pos) {
        return this->attachmentScores[pos];
    }

    AttachmentScore const &operator[](unsigned pos) const {
        return this->attachmentScores[pos];
    }

    unsigned size() {
        return this->attachmentScores.size();
    }

    void resize(unsigned newSize) {
        this->attachmentScores.resize(newSize);
    }

    double &hetScore(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].hetScore();
    }

    double const &hetScore(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].hetScore();
    }

    double &hetSumScore(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].hetSumScore();
    }

    double const &hetSumScore(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].hetSumScore();
    }

    double &homScore(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].homScore();
    }

    double const &homScore(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].homScore();
    }

    double &lossWildScore(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].lossWildScore();
    }

    double const &lossWildScore(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].lossWildScore();
    }

    double &lossAltScore(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].lossAltScore();
    }

    double const &lossAltScore(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].lossAltScore();
    }

    double &lossAltRScore(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].lossAltRScore();
    }

    double const &lossAltRScore(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].lossAltRScore();
    }

    double &lcaScore(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].lcaScore();
    }

    double const &lcaScore(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].lcaScore();
    }

    double &lcaRScore(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].lcaRScore();
    }

    double const &lcaRScore(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].lcaRScore();
    }

    double &grantChildHetSumScore(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].grantChildHetSumScore();
    }

    double const &grantChildHetSumScore(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].grantChildHetSumScore();
    }

    double &numAltPoss(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].numAltPoss();
    }

    double const &numAltPoss(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].numAltPoss();
    }

    double &numAltRPoss(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].numAltRPoss();
    }

    double const &numAltRPoss(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].numAltRPoss();
    }

    double &numInnerNodes(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].numInnerNodes();
    }

    double const &numInnerNodes(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].numInnerNodes();
    }

    double &numInnerChildNodes(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].numInnerChildNodes();
    }

    double const &numInnerChildNodes(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].numInnerChildNodes();
    }

    double &numLcaRPoss(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].numLcaRPoss();
    }

    double const &numLcaRPoss(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].numLcaRPoss();
    }

    double &childHetSumScore(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].childHetSumScore();
    }

    double const &childHetSumScore(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].childHetSumScore();
    }

    void computeLogHetScoreLeaf(unsigned attachPoint, double wtScore, double hetScore) {
        this->hetScore(attachPoint) = hetScore - wtScore;
    }

    void computeLogHetScoreInnerNode(unsigned attachPoint, double hetScoreLeft, double hetScoreRigth) {
        this->hetScore(attachPoint) = hetScoreLeft + hetScoreRigth;
    }

    void computeLogHomScoreLeaf(unsigned attachPoint, double wtScore, double homScore) {
        this->homScore(attachPoint) = homScore - wtScore;
    }

    void computeLogHomScoreInnerNode(unsigned attachPoint, double homScoreLeft, double homScoreRigth) {
        this->homScore(attachPoint) = homScoreLeft + homScoreRigth;
    }

    double computeLogLossScoreLeaf() {
        return -INFINITY;
    }

    void computeLogLossWildScoreLeaf(unsigned attachPoint) {
        this->lossWildScore(attachPoint) = computeLogLossScoreLeaf();
    }

    void computeLogLossAltScoreLeaf(unsigned attachPoint) {
        this->lossAltScore(attachPoint) = computeLogLossScoreLeaf();
    }

    template<typename TTree>
    void computeLogLossScoreInnerNode(
            TTree const &tree,
            unsigned attachPoint){

        auto it = out_edges(attachPoint, tree).first;
        unsigned lN = target(*it, tree);                // left node id
        unsigned rN = target(*(it + 1), tree);          // right node id

        double leftContAltLos = -INFINITY;
        // loss wild type chromosome
        double leftContWildLos = -INFINITY;
        // loss alternative chromosome
        double rightContAltLos = -INFINITY;
        // loss wild type chromosome
        double rightContWildLos = -INFINITY;

        bool computeValues = false;

        // check if child nodes are inner node
        if (tree[lN].sample == -1)
        {
            // loss alternative chromosome
            leftContAltLos = addLogProb(this->lossAltScore(lN), 0) + this->hetScore(rN);
            // loss wild type chromosome
            leftContWildLos = addLogProb(this->lossWildScore(lN), this->homScore(lN)) + this->hetScore(rN);

            this->numAltPoss(attachPoint) += this->numAltPoss(lN) + 1;

            this->numAltRPoss(attachPoint) += this->numAltPoss(lN);

            computeValues = true;

        }
        if (tree[rN].sample == -1)
        {
            // loss alternative chromosome
            rightContAltLos = addLogProb(this->lossAltScore(rN), 0) + this->hetScore(lN);
            // loss wild type chromosome
            rightContWildLos = addLogProb(this->lossWildScore(rN), this->homScore(rN)) + this->hetScore(lN);

            this->numAltPoss(attachPoint) += this->numAltPoss(rN) + 1;

            this->numAltRPoss(attachPoint) += this->numAltPoss(rN);

            computeValues = true;
        }

        if(computeValues == true)
        {
            // compute P_LA
            this->lossAltScore(attachPoint) = addLogProb(leftContAltLos, rightContAltLos);

            // compute P_LAR
            this->lossAltRScore(attachPoint) = addLogProb(this->lossAltScore(lN) + this->hetScore(rN),
                                                          this->hetScore(lN) + this->lossAltScore(rN));
            // compute P_LW
            this->lossWildScore(attachPoint) = addLogProb(leftContWildLos, rightContWildLos);
        }
    }

    template<typename TTree>
    void computeLogLcaScoreInnerNode(
            TTree const &tree,
            unsigned attachPoint){

        auto it = out_edges(attachPoint, tree).first;
        unsigned lN = target(*it, tree);                // left node id
        unsigned rN = target(*(it + 1), tree);          // right node id

        this->hetSumScore(attachPoint) = addLogProb(hetScore(attachPoint),
                addLogProb(hetSumScore(lN), hetSumScore(rN)));

        // check if child nodes are inner nodes
        if(tree[lN].sample == -1 && tree[rN].sample == -1) {
//            this->hetSumScore(attachPoint) = addLogProb(hetScore(attachPoint),
//                    addLogProb(hetSumScore(lN), hetSumScore(rN)));

            this->childHetSumScore(attachPoint) = addLogProb(this->hetScore(lN), this->hetScore(rN));
//            this->grantChildHetSumScore(attachPoint) = subLogProb(
//                    this->hetSumScore(attachPoint), addLogProb(
//                    this->hetScore(attachPoint), this->childHetSumScore(attachPoint)));
            this->lcaScore(attachPoint) = hetSumScore(lN) + hetSumScore(rN);

            double hetContLcaR = this->hetScore(lN) + this->hetScore(rN);
            double leftContLcaR = this->childHetSumScore(lN) + this->hetScore(rN);
            double rightContLcaR = this->hetScore(lN) + this->childHetSumScore(rN);
            this->lcaRScore(attachPoint) = subLogProb(
                    subLogProb(this->lcaScore(attachPoint), hetContLcaR),
                    addLogProb(leftContLcaR, rightContLcaR));

            this->numInnerNodes(attachPoint) = this->numInnerNodes(lN) + this->numInnerNodes(rN) + 2;
            this->numInnerChildNodes(attachPoint) = 2;
            this->numLcaRPoss(attachPoint) = (this->numInnerNodes(lN) + 1) * (this->numInnerNodes(rN) + 1)
                                             - 1
                                             - this->numInnerChildNodes(rN)
                                             - this->numInnerChildNodes(lN);
            return;
        }

        if(tree[lN].sample == -1)
        {
            this->numInnerNodes(attachPoint) = 1 + this->numInnerNodes(lN);
            this->numInnerChildNodes(attachPoint) = 1;
            return;
        }

        if(tree[rN].sample == -1)
        {
            this->numInnerNodes(attachPoint) = 1 + this->numInnerNodes(rN);
            this->numInnerChildNodes(attachPoint) = 1;
            return;
        }
    }


};

#endif
