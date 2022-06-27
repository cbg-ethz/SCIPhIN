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

#ifndef ATTACHMENTSCORES_H
#define ATTACHMENTSCORES_H

#include <limits.h>

#include "sciphi_config.h"

// Forward declaration
double addLogProbWeight(double x, double y, double nu);

double addLogProb(double x, double y);

double subLogProb(double x, double y);

template <typename TTreeType>
class Config;

/*
 * This class stores values associated with the attachment of mutations
 * to nodes in the phylogenetic tree.
 */

struct AttachmentScore {

    typedef std::array<double, 18> TAttachmentScore;

    TAttachmentScore attachmentScore;


    enum Type {
                // all cells below the node are heterozygous
                E_hetScore = 0,
                // all cells below the node are homozygous
                E_homScore = 1,
                // the sum of all hetScores below the node
                E_hetSumScore = 2,
                // a mutation is gained in the node and the mutated/altered allele is lost in any node below
                E_lossAltScore = 3,
                // a mutation is gained in the node and the mutated/altered allele is lost in any node below except
                // for the child nodes
                E_lossAltRScore = 4,
                // a mutation is gained in the node and the reference allele is lost in any node below
                E_lossWildScore = 5,
                // probability that the node is the lowest common ancestor of two independent mutation events
                E_lcaScore = 6,
                // the probability that one of the child nodes is mutated
                E_childHetSumScore = 7,
                // probability that the node is the lowest common ancestor of two independent mutation events with
                // additional restrictions
                E_lcaRScore = 8,
                // probability of a mutation below the current node, but not in itself or the child nodes
                E_grantChildHetSumScore = 9,
                // final combined score
                E_finalScore = 10,
                // number of different possibilities to loose the current mutation (but not directly in the child nodes
                E_numAltRPoss = 11,
                // number of inner nodes below the current one
                E_numInnerNodes = 12,
                // number of child nodes that are inner nodes (0, 1, or 2)
                E_numInnerChildNodes = 13,
                // number of possibilities for a parallel mutation below the current node
                E_numLcaRPoss = 14,
                // number of possibilities for a parallel mutation below the current node
                E_numFirstMutWildLoss = 15,
                E_numFirstMutAltLoss = 16,
                E_numFirstLca = 17
    };

    // init everything but the counter to -INFINITY because we are doing the computation in
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
                                     -INFINITY,
                                     0,
                                     0,
                                     0,
                                     0,
                                     0,
                                     0,
                                     0}}) {};

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

    double const &numFirstMutWildLoss() const {
        return this->attachmentScore[E_numFirstMutWildLoss];
    }

    double &numFirstMutWildLoss() {
        return const_cast<double &>(static_cast<const AttachmentScore*>(this)->numFirstMutWildLoss());
    }
    double const &numFirstMutAltLoss() const {
        return this->attachmentScore[E_numFirstMutAltLoss];
    }

    double &numFirstMutAltLoss() {
        return const_cast<double &>(static_cast<const AttachmentScore*>(this)->numFirstMutAltLoss());
    }
    double const &numFirstLca() const {
        return this->attachmentScore[E_numFirstLca];
    }

    double &numFirstLca() {
        return const_cast<double &>(static_cast<const AttachmentScore*>(this)->numFirstLca());
    }



    double &finalScore() {
        return this->attachmentScore[E_finalScore];
    }

    double const &finalScore() const {
        return this->attachmentScore[E_finalScore];
    }

    // this function takes two scores in log space, exponentiates them and returns the log of the sum
    void addInRealSpace(AttachmentScore const &rightSide) {
        for (unsigned i = 0; i < 11; ++i) {
            this->attachmentScore[i] = addLogProb(this->attachmentScore[i], rightSide.attachmentScore[i]);
        }
        for (unsigned i = 11; i < this->attachmentScore.size(); ++i) {
            this->attachmentScore[i] += rightSide.attachmentScore[i];
        }
    }

    AttachmentScore &operator/=(double rightSide) {
        for (unsigned i = 0; i < this->attachmentScore.size(); ++i) {
            this->attachmentScore[i] /= rightSide;
        }
        return *this;
    }

    // This function computes the contribution of the different mutation types. In other words: How much of the total
    // score is contributed by the heterozygous, the homozygous, the loss, and the parallel model. Each contribution
    // is the product of three parts: its likelihood, its model parameter value, and the inverse of the number of
    // different placement possibilities

    template <typename TTreeType>
    std::array<double, 6> computeMutTypeContribution(Config<TTreeType> const & config,
            AttachmentScore const & sumScore,
            bool isSumScore = false)
    {

        double logPHet = -INFINITY;
        double logPHom = -INFINITY;
        double logPLossWild = -INFINITY;
        double logPLossMut = -INFINITY;
        double logPPara = -INFINITY;
        double logPD = -INFINITY;

        // The weight of the heterozygous score depends on the weights of the other mutation types. We initialize it
        // to one and substrat the other contributions afterwards.
        double hetWeight = 1.0;

        // Compute the homozygous contribution
        unsigned numInnerNodes = config.getNumSamples() - 1;
        if (config.learnZygocity)
        {
            logPHom = this->homScore()
                      - std::log(numInnerNodes)
                      + std::log(config.getParam(Config<TTreeType>::E_nu));
            hetWeight -= config.getParam(Config<TTreeType>::E_nu);
        }

        // Compute the loss contribution. THe loss contribution consists of the model where the reference allele is
        // lost and where the alternative allele is lost
        if (config.computeLossScore)
        {
            logPLossWild= this->lossWildScore() 
                       - std::log(sumScore.numFirstMutWildLoss())
                       + std::log(config.getParam(Config<TTreeType>::E_lambdaWildLoss));
            if (!isSumScore){    
                    logPLossWild-= std::log(this->numInnerNodes());
            }
            hetWeight -= config.getParam(Config<TTreeType>::E_lambdaWildLoss);
            
            logPLossMut = this->lossAltRScore() 
                       - std::log(sumScore.numFirstMutAltLoss())
                       + std::log(config.getParam(Config<TTreeType>::E_lambdaMutLoss));
            if (!isSumScore){    
                      logPLossMut -= std::log(this->numAltRPoss());
            }
            hetWeight -= config.getParam(Config<TTreeType>::E_lambdaMutLoss);
        }

        // Compute the parallel contribution.
        if (config.computeParallelScore)
        {
            logPPara = this->lcaRScore()
                       - std::log(sumScore.numFirstLca())
                       + std::log(config.getParam(Config<TTreeType>::E_parallel));
            if (!isSumScore){
                 if (this->numLcaRPoss() > 0)
                 {
                    logPPara -= std::log(this->numLcaRPoss());
                 }
            }
            hetWeight -= config.getParam(Config<TTreeType>::E_parallel);
        }

        // Finally compute the heterozygous contribution.
        unsigned numNodes = config.getNumSamples() * 2 - 1;
        logPHet = this->hetScore()
                  - std::log(numNodes)
                  + std::log(hetWeight);

        logPD = addLogProb(addLogProb(addLogProb(addLogProb(logPHet,logPHom),logPLossWild),logPLossMut),logPPara);
        
        return {{logPHet, logPHom, logPLossWild, logPLossMut, logPPara, logPD}};
    }

    // This functions computes the final score of an attachment point
    template <typename TTreeType>
    void computeFinalScore(Config<TTreeType> const & config,
            AttachmentScore const & sumScore,
            bool isLeaf,
            bool isSumScore = false) {

        static std::array<double, 6> res;
        res = this->computeMutTypeContribution(config, sumScore, isSumScore);

        // If the node is a leaf just return the weighted hetero score
        if (isLeaf) {
            this->finalScore() = this->hetScore() - std::log(config.getNumSamples() * 2 - 1);
            return;
        }

        this->finalScore() = res[5];
        if (std::isnan(this->finalScore()))
        {
            this->finalScore() = -std::numeric_limits<double>::infinity();
        }
        return;
    }
};

std::ostream &operator<<(std::ostream &os, AttachmentScore const &obj) {
    os << "hetScore:         " << obj.hetScore() << "\n";
    os << "homScore:         " << obj.homScore() << "\n";
    os << "hetSumScore:      " << obj.hetSumScore() << "\n";
    os << "lossAltScore:     " << obj.lossAltScore() << "\n";
    os << "lossAltRScore:    " << obj.lossAltRScore() << "\n";
    os << "lossWildScore:    " << obj.lossWildScore() << "\n";
    os << "lcaScore:         " << obj.lcaScore() << "\n";
    os << "lcaRScore:         " << obj.lcaRScore() << "\n";
    os << "childHetSumScore: " << obj.childHetSumScore() << "\n";
    os << "numLcaRPoss:      " << obj.numLcaRPoss() << "\n";
    os << "numInnerNodes:  " << obj.numInnerNodes() << "\n";
    os << "numAltRPoss:      " << obj.numAltRPoss() << "\n";
    os << "finalScore:       " << obj.finalScore() << "\n";
    return os;
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

    double &numFirstMutWildLoss(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].numFirstMutWildLoss();
    }

    double const &numFirstMutWildLoss(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].numFirstMutWildLoss();
    }

    double &numFirstMutAltLoss(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].numFirstMutAltLoss();
    }

    double const &numFirstMutAltLoss(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].numFirstMutAltLoss();
    }

    double &numFirstLca(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].numFirstLca();
    }

    double const &numFirstLca(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].numFirstLca();
    }

    double &childHetSumScore(unsigned attachPoint) {
        return this->attachmentScores[attachPoint].childHetSumScore();
    }

    double const &childHetSumScore(unsigned attachPoint) const {
        return this->attachmentScores[attachPoint].childHetSumScore();
    }

    // We divide by the wild type score. In doing so we do not need to account for it in every computation because
    // the log of the quotient for the wild type score would be 0. However, when computing the tree score we need to
    // multiply the wild type scores.
    void computeLogHetScoreLeaf(unsigned attachPoint, double wtScore, double hetScore) {
        this->hetScore(attachPoint) = hetScore - wtScore;
    }

    void computeLogHetScoreInnerNode(unsigned attachPoint, double hetScoreLeft, double hetScoreRigth) {
        this->hetScore(attachPoint) = hetScoreLeft + hetScoreRigth;
    }

    // We divide by the wild type score. In doing so we do not need to account for it in every computation because
    // the log of the quotient for the wild type score would be 0. However, when computing the tree score we need to
    // multiply the wild type scores.
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


        double leftContAltLos = -INFINITY;              // left child node alternative allele lost contribution
        double leftContWildLos = -INFINITY;             // left child node wild type allele lost contribution
        double rightContAltLos = -INFINITY;             // right child node alternative allele lost contribution
        double rightContWildLos = -INFINITY;            // right child node wild type allele lost contribution

        bool computeValues = false;

        // Check if left child node is inner node and compute the corresponding contributions
        if (tree[lN].sample == -1)
        {
            // Loss alternative chromosome
            // Combine the likelihoods from the left child node, which consists of the likelihood that the allele was
            // lost below the left child node or directly in the left child node. The combined score has to be
            // multiplied by the probability that the right subtree is heterozygous
            leftContAltLos = addLogProb(this->lossAltScore(lN), 0) + this->hetScore(rN);

            // Loss wild type chromosome
            // Combine the likelihoods from the left child node, which consists of the likelihood that the allele was
            // lost below the left child node or directly in the left child node. The combined score has to be
            // multiplied by the probability that the right subtree is heterozygous
            leftContWildLos = addLogProb(this->lossWildScore(lN), this->homScore(lN)) + this->hetScore(rN);

            // Count the number of different possibilities
            this->numInnerNodes(attachPoint) = this->numInnerNodes(lN) + 1;
            this->numAltRPoss(attachPoint) = this->numInnerNodes(lN);

            computeValues = true;

        }

        // Check if right child node is inner node and compute the corresponding contributions
        if (tree[rN].sample == -1)
        {
            // Loss alternative chromosome
            // Combine the likelihoods from the right child node, which consists of the likelihood that the allele was
            // lost below the right child node or directly in the right child node. The combined score has to be
            // multiplied by the probability that the left subtree is heterozygous
            rightContAltLos = addLogProb(this->lossAltScore(rN), 0) + this->hetScore(lN);

            // Loss wild type chromosome
            // Combine the likelihoods from the right child node, which consists of the likelihood that the allele was
            // lost below the right child node or directly in the right child node. The combined score has to be
            // multiplied by the probability that the left subtree is heterozygous
            rightContWildLos = addLogProb(this->lossWildScore(rN), this->homScore(rN)) + this->hetScore(lN);

            // Count the number of different possibilities
            this->numInnerNodes(attachPoint) += this->numInnerNodes(rN) + 1;
            this->numAltRPoss(attachPoint) += this->numInnerNodes(rN);

            computeValues = true;
        }

        // In case one of the child nodes was an inner node
        if(computeValues == true)
        {
            // Compute the probability that the alternative allele was lost below the current node by adding the left
            // and right likelihoods.
            this->lossAltScore(attachPoint) = addLogProb(leftContAltLos, rightContAltLos);

            // Compute the probability that the alternative allele was lost below the current node but not in the
            // child nodes by adding the likelihood that the mutation was lost below the left or right child nodes
            this->lossAltRScore(attachPoint) = addLogProb(this->lossAltScore(lN) + this->hetScore(rN),
                                                          this->hetScore(lN) + this->lossAltScore(rN));

            // Compute the probability that the reference allele was lost below the current node by adding the left
            // and right likelihoods.
            this->lossWildScore(attachPoint) = addLogProb(leftContWildLos, rightContWildLos);
        }
        this->numFirstMutWildLoss(attachPoint) = isinf(this->lossWildScore(attachPoint)) ? 0 : 1;
        this->numFirstMutAltLoss(attachPoint) = isinf(this->lossAltRScore(attachPoint)) ? 0 : 1;
    }

    template<typename TTree>
    void computeLogLcaScoreInnerNode(
            TTree const &tree,
            unsigned attachPoint){

        auto it = out_edges(attachPoint, tree).first;
        unsigned lN = target(*it, tree);                // left node id
        unsigned rN = target(*(it + 1), tree);          // right node id

        // Compute the likelihood that any node in the current subtree is mutated (excluding leafs)
        this->hetSumScore(attachPoint) = addLogProb(hetScore(attachPoint),
                addLogProb(hetSumScore(lN), hetSumScore(rN)));

        // If the child nodes are leafs assign -INFINITY as parallel mutation in leafs are prohibited
        double lChildHetScore = tree[lN].sample == -1 ? this->hetScore(lN) : -INFINITY;
        double rChildHetScore = tree[rN].sample == -1 ? this->hetScore(rN) : -INFINITY;

        // Compute the likelihood that one or both children are mutated
        this->childHetSumScore(attachPoint) = addLogProb(lChildHetScore, rChildHetScore);

        // Compute the likelihood that there is one mutation in the left and one mutation in the right subtree
        this->lcaScore(attachPoint) = hetSumScore(lN) + hetSumScore(rN);

        // The final lowest common ancestor score is computed by substracting all prohibited cases from the
        // unrestricted score (lcaScore)
        // Both child nodes cannot be mutated (it would be more likely that the current node was mutated)
        double hetContLcaR = lChildHetScore + rChildHetScore;
        // The right child node cannot be mutated together with the children of the left child
        double leftContLcaR = this->childHetSumScore(lN) + rChildHetScore;
        // The left child node cannot be mutated together with the children of the right child
        double rightContLcaR = lChildHetScore + this->childHetSumScore(rN);
        // Substract the prohibited cases
        
        this->lcaRScore(attachPoint) = subLogProb(
                subLogProb(this->lcaScore(attachPoint), hetContLcaR),
                addLogProb(leftContLcaR, rightContLcaR));

        std::cout.precision(20);

        // Count the number of placing possibilities
        if(tree[lN].sample == -1 && tree[rN].sample == -1) {
            this->numInnerNodes(attachPoint) = this->numInnerNodes(lN) + this->numInnerNodes(rN) + 2;
            this->numInnerChildNodes(attachPoint) = 2;
            this->numLcaRPoss(attachPoint) = (this->numInnerNodes(lN) + 1) * (this->numInnerNodes(rN) + 1)
                                             - 1
                                             - this->numInnerChildNodes(rN)
                                             - this->numInnerChildNodes(lN);
            this->numFirstLca(attachPoint) += 1;
            return;

        }

        // Count the number of placing possibilities
        if(tree[lN].sample == -1)
        {
            this->numInnerNodes(attachPoint) = 1 + this->numInnerNodes(lN);
            this->numInnerChildNodes(attachPoint) = 1;
            this->numFirstLca(attachPoint) += 1;
            return;
        }

        // Count the number of placing possibilities
        if(tree[rN].sample == -1)
        {
            this->numInnerNodes(attachPoint) = 1 + this->numInnerNodes(rN);
            this->numInnerChildNodes(attachPoint) = 1;
            this->numFirstLca(attachPoint) += 1;
            return;
        }
    }
};

#endif
