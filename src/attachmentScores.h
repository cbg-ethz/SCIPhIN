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

double addLog_nan_x(double result, double rlNScore);

/*
 * This class stores values associated with the attachment of mutations
 * to nodes in the phylogenetic tree.
 */
struct AttachmentScore {

    typedef std::array<double, 9> TAttachmentScore;

    TAttachmentScore attachmentScore;


    enum Type {
        E_hetScore = 0,         // all cells below the node are heterozygous
        E_homScore = 1,         // all cells below the node are homozygous
        E_hetSumScore = 2,      // the sum of all hetScores below the node
        E_lossAltScore = 3,     // a mutation is gained in the node and the mutated/altered allele is lost in any
        // node below
                E_lossAltScoreR = 4,    // a mutation is gained in the node and the mutated/altered allele is lost in any
        // node below except for the child nodes
                E_lossWildScore = 5,     // a mutation is gained in the node and the reference allele is lost in any node below
        // except for the child nodes
                E_lcaScore = 6,         // probability that the node is the lowest common ancestor of two independent mutation
        // events
                E_lcaScoreR = 7,        // probability that the node is the lowest common ancestor of two independent mutation
        // events with additional restrictions
                E_finalScore = 8        // final combined score
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
                                     -INFINITY}}) {};

    // transform attachment scores into log
    void log() {
        for (unsigned i = 0; i < this->attachmentScore.size(); ++i) {
            this->attachmentScore[i] = std::log(this->attachmentScore[i]);
        }
    }

    // exponentiate attachment scores
    void exp() {
        for (unsigned i = 0; i < 7; ++i) {
            this->attachmentScore[i] = std::exp(this->attachmentScore[i]);
        }
    }

    double &hetScore() {
        return this->attachmentScore[E_hetScore];
    }

    double const &hetScore() const {
        return this->attachmentScore[E_hetScore];
    }

    double &sumHetScore() {
        return this->attachmentScore[E_hetSumScore];
    }

    double const &sumHetScore() const {
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

    double &lcaScore() {
        return this->attachmentScore[E_lcaScore];
    }

    double const &paralleleScore() const {
        return this->attachmentScore[E_lcaScore];
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

    AttachmentScore cellProbReal(AttachmentScore const &sumReal) {
        AttachmentScore result = {};
        result.hetScore() = this->hetScore() / sumReal.hetScore();
        result.homScore() = this->homScore() / sumReal.homScore();
        result.lossWildScore() = this->lossWildScore() / sumReal.lossWildScore();
        result.lossAltScore() = this->lossAltScore() / sumReal.lossAltScore();
        result.finalScore() = this->finalScore() / sumReal.finalScore();

        return result;
    }

    void setMinusInfinity() {
        this->hetScore() = -INFINITY;
        this->homScore() = -INFINITY;
        this->lossWildScore() = -INFINITY;
        this->lossAltScore() = -INFINITY;
        this->finalScore() = -INFINITY;
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
    os << obj.hetScore() << "|" <<
       obj.homScore() << "|" <<
       obj.lossWildScore() << "|" <<
                                  obj.lossAltScore() << "|" <<
       obj.paralleleScore() << "|" <<
       obj.finalScore();
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

    double &hetScore(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].hetScore();
    }

    double const &hetScore(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].hetScore();
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

    double &lossHomScore(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].lossAltScore();
    }

    double const &lossHomScore(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].lossAltScore();
    }

    unsigned size() {
        return this->attachmentScores.size();
    }

    void resize(unsigned newSize) {
        this->attachmentScores.resize(newSize);
    }

    /*
    AttachmentScore max() {
        AttachmentScore result = this->attachmentScores[0];

        for (unsigned i = 1; i < this->size(); ++i) {
            if (result.hetScore() < this->attachmentScores[i].hetScore()) {
                result.hetScore() = this->attachmentScores[i].hetScore();
            }
            if (result.homScore() < this->attachmentScores[i].homScore()) {
                result.homScore() = this->attachmentScores[i].homScore();
            }
            if (result.lossWildScore() < this->attachmentScores[i].lossWildScore()) {
                result.lossWildScore() = this->attachmentScores[i].lossWildScore();
            }
            if (result.lossAltScore() < this->attachmentScores[i].lossAltScore()) {
                result.lossAltScore() = this->attachmentScores[i].lossAltScore();
            }
        }
        return result;
    }
     */

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
        return std::nan("");
    }

    void computeLogLossWildScoreLeaf(unsigned attachPoint) {
        this->lossWildScore(attachPoint) = computeLogLossScoreLeaf();
    }

    void computeLogLossHomScoreLeaf(unsigned attachPoint) {
        this->lossHomScore(attachPoint) = computeLogLossScoreLeaf();
    }


    template<typename TTree>
    void computeLogLossHomScoreInnerNode(
            TTree const &tree,
            unsigned attachPoint) {
        double result = std::nan("");
        auto it = out_edges(attachPoint, tree).first;
        unsigned lN = target(*it, tree); // left node i
        unsigned rN = target(*(it + 1), tree); // right node id

        // check if left child is inner node

        bool lNisInnerNode = tree[lN].sample == -1;
        if (lNisInnerNode) {

            // check if left child of the left child is an inner node
            auto lIt = out_edges(lN, tree).first;
            unsigned llN = target(*lIt, tree);// left child of left node
            bool llNisInnerNode = tree[llN].sample == -1;
            if (llNisInnerNode) {
                result = this->hetScore(attachPoint)
                         - this->hetScore(llN)
                         + this->homScore(llN);
            }

            // check if right child of the left child is an inner node
            unsigned rlN = target(*(lIt + 1), tree); // right child of left node
            bool rlNisInnerNode = tree[rlN].sample == -1;
            if (rlNisInnerNode) {
                double rlNScore = this->hetScore(attachPoint)
                                  - this->hetScore(rlN)
                                  + this->homScore(rlN);

                result = addLog_nan_x(result, rlNScore);
            }

            //check if left child has a lossture score
            if (!std::isnan(this->lossHomScore(lN))) {
                result = addLogProb(result,
                                    this->lossHomScore(lN) + this->hetScore(rN));
            }
        }

        // check if left right is inner node
        bool rNisInnerNode = tree[rN].sample == -1;
        if (rNisInnerNode) {

            // check if left child of the left child is an inner node
            auto rIt = out_edges(rN, tree).first;
            unsigned lrN = target(*rIt, tree);// left child of left node
            bool lrNisInnerNode = tree[lrN].sample == -1;
            if (lrNisInnerNode) {
                double lrNScore = this->hetScore(attachPoint)
                                  - this->hetScore(lrN)
                                  + this->homScore(lrN);
                result = addLog_nan_x(result, lrNScore);
            }

            // check if right child of the left child is an inner node
            unsigned rrN = target(*(rIt + 1), tree); // right child of left node
            bool rrNisInnerNode = tree[rrN].sample == -1;
            if (rrNisInnerNode) {
                double rrNScore = this->hetScore(attachPoint)
                                  - this->hetScore(rrN)
                                  + this->homScore(rrN);
                result = addLog_nan_x(result, rrNScore);
            }

            //check if right child has a lossture score
            if (!std::isnan(this->lossHomScore(rN))) {
                result = addLogProb(result,
                                    this->lossHomScore(rN) + this->hetScore(lN));
            }
        }
        this->lossHomScore(attachPoint) = result;
    }

    template<typename TTree>
    void computeLogLossWildScoreInnerNode(
            TTree const &tree,
            unsigned attachPoint){

        double result = std::nan("");
        auto it = out_edges(attachPoint, tree).first;
        unsigned lN = target(*it, tree); // left node i
        unsigned rN = target(*(it + 1), tree); // right node id

        // check if left child is inner node

        bool lNisInnerNode = tree[lN].sample == -1;
        if (lNisInnerNode) {

            // check if left child of the left child is an inner node
            auto lIt = out_edges(lN, tree).first;
            unsigned llN = target(*lIt, tree);// left child of left node
            bool llNisInnerNode = tree[llN].sample == -1;
            if (llNisInnerNode) {
                result = this->hetScore(attachPoint)
                         - this->hetScore(llN);
            }

            // check if right child of the left child is an inner node
            unsigned rlN = target(*(lIt + 1), tree); // right child of left node
            bool rlNisInnerNode = tree[rlN].sample == -1;
            if (rlNisInnerNode) {
                double rlNScore = this->hetScore(attachPoint)
                                  - this->hetScore(rlN);

                result = addLog_nan_x(result, rlNScore);
            }

            //check if left child has a lossture score
            if (!std::isnan(this->lossWildScore(lN))) {
                result = addLogProb(result,
                                    this->lossWildScore(lN) + this->hetScore(rN));
            }
        }

        // check if left right is inner node
        bool rNisInnerNode = tree[rN].sample == -1;
        if (rNisInnerNode) {

            // check if left child of the left child is an inner node
            auto rIt = out_edges(rN, tree).first;
            unsigned lrN = target(*rIt, tree);// left child of left node
            bool lrNisInnerNode = tree[lrN].sample == -1;
            if (lrNisInnerNode) {
                double lrNScore = this->hetScore(attachPoint)
                                  - this->hetScore(lrN);
                result = addLog_nan_x(result, lrNScore);
            }

            // check if right child of the left child is an inner node
            unsigned rrN = target(*(rIt + 1), tree); // right child of left node
            bool rrNisInnerNode = tree[rrN].sample == -1;
            if (rrNisInnerNode) {
                double rrNScore = this->hetScore(attachPoint)
                                  - this->hetScore(rrN);
                result = addLog_nan_x(result, rrNScore);
            }

            //check if right child has a lossture score
            if (!std::isnan(this->lossWildScore(rN))) {
                result = addLogProb(result,
                                    this->lossWildScore(rN) + this->hetScore(lN));
            }
        }
        this->lossWildScore(attachPoint) = result;
    }
};

#endif
