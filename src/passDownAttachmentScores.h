/**
 * SCIPhI: Single-cell mutation identification via phylogenetic inference
 * <p>
 * Copyright (C) 2019 ETH Zurich, Jochen Singer
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

#ifndef PASSDOWNATTACHMENTSCORES_H
#define PASSDOWNATTACHMENTSCORES_H

/*
 * This class stores values associated with the attachment of mutations to nodes in the phylogenetic tree.
 * However, in contrast to AttachmentScore this function only stores values when it comes to the computation of the
 * probability of a certain cell to be affected by a certain mutation. In order to compute this value the previously
 * computed attachment scores are passed down from the root node to the child nodes.
 */
struct PassDownAttachmentScore {

    typedef std::array<double, 5> TPassDownAttachmentScore;

    TPassDownAttachmentScore passDownAttachmentScore;

    enum Type {
        // A mutation was gained in a node towards the root and is lost in the current node
                E_lossAltInCurrentNodeScore = 0,
        // A mutation was gained in a node towards the root, but not the parent node and is lost in the current node
                E_lossAltInCurrentNodeRScore = 1,
        // A mutation is present in any node of the subtree rooted at the sibling node
                E_sib = 2,
        // The current node and any other node in a different subtree is mutated
                E_parallel = 3
    };

    // init everything to -INFINITY because we are doing the computation in
    // log space and exp(-INFINITY) = 0
    PassDownAttachmentScore() :
            passDownAttachmentScore({{-INFINITY,
                                             -INFINITY,
                                             -INFINITY,
                                             -INFINITY,
                                             -INFINITY,
                                     }}) {};

    double const &lossAltInCurrentNodeScore() const {
        return this->passDownAttachmentScore[E_lossAltInCurrentNodeScore];
    }

    double &lossAltInCurrentNodeScore() {
        return const_cast<double &>(static_cast<const PassDownAttachmentScore *>(this)->lossAltInCurrentNodeScore());
    }

    double const &lossAltInCurrentNodeRScore() const {
        return this->passDownAttachmentScore[E_lossAltInCurrentNodeRScore];
    }

    double &lossAltInCurrentNodeRScore() {
        return const_cast<double &>(static_cast<const PassDownAttachmentScore *>(this)->lossAltInCurrentNodeRScore());
    }

    double const &paralleleScore() const {
        return this->passDownAttachmentScore[E_parallel];
    }

    double &paralleleScore() {
        return const_cast<double &>(static_cast<const PassDownAttachmentScore *>(this)->paralleleScore());
    }

    double const &sibNodeScore() const {
        return this->passDownAttachmentScore[E_sib];
    }

    double &sibNodeScore() {
        return const_cast<double &>(static_cast<const PassDownAttachmentScore *>(this)->sibNodeScore());
    }
};

/*
 * This class combines many AttachmentScores into one object
 */
struct PassDownAttachmentScores {
    typedef PassDownAttachmentScore TPassDownAttachmentScore;
    typedef std::vector<PassDownAttachmentScore> TPassDownAttachmentScores;

    TPassDownAttachmentScores passDownAttachmentScores;

    PassDownAttachmentScore &operator[](unsigned pos) {
        return this->passDownAttachmentScores[pos];
    }

    PassDownAttachmentScore const &operator[](unsigned pos) const {
        return this->passDownAttachmentScores[pos];
    }

    unsigned size() {
        return this->passDownAttachmentScores.size();
    }

    void resize(unsigned newSize) {
        this->passDownAttachmentScores.resize(newSize);
    }

    double const &lossAltInCurrentNodeScore(unsigned attachPoint) const {
        return this->passDownAttachmentScores[attachPoint].lossAltInCurrentNodeScore();
    }

    double &lossAltInCurrentNodeScore(unsigned attachPoint) {
        return const_cast<double &>(static_cast<const PassDownAttachmentScores *>(this)->lossAltInCurrentNodeScore(
                attachPoint));
    }

    double const &lossAltInCurrentNodeRScore(unsigned attachPoint) const {
        return this->passDownAttachmentScores[attachPoint].lossAltInCurrentNodeRScore();
    }

    double &lossAltInCurrentNodeRScore(unsigned attachPoint) {
        return const_cast<double &>(static_cast<const PassDownAttachmentScores *>(this)->lossAltInCurrentNodeRScore(
                attachPoint));
    }

    double const &paralleleScore(unsigned attachPoint) const {
        return this->passDownAttachmentScores[attachPoint].paralleleScore();
    }

    double &paralleleScore(unsigned attachPoint) {
        return const_cast<double &>(static_cast<const PassDownAttachmentScores *>(this)->paralleleScore(attachPoint));
    }

    double const &sibNodeScore(unsigned attachPoint) const {
        return this->passDownAttachmentScores[attachPoint].sibNodeScore();
    }

    double &sibNodeScore(unsigned attachPoint) {
        return const_cast<double &>(static_cast<const PassDownAttachmentScores *>(this)->sibNodeScore(attachPoint));
    }

    template<typename TTree>
    void computeLogLossInCurrentInnerNode(
            TTree const &tree,
            AttachmentScores &attachmentScores,
            unsigned attachPoint) {

        unsigned pN = source(*in_edges(attachPoint, tree).first, tree);     // parent node id
        unsigned sN = getSibling(tree, attachPoint);                        // sibling node id

        // If the current node is the root or a leaf do noting
        if (pN == num_vertices(tree) - 1 || tree[attachPoint].sample != -1) {
            return;
        }

        // The likelihood that the mutation was lost in the current node
        this->lossAltInCurrentNodeScore(attachPoint) = addLogProb(this->lossAltInCurrentNodeScore(pN)
                                                                  + attachmentScores.hetScore(sN),
                                                                  attachmentScores.hetScore(pN)
                                                                  - attachmentScores.hetScore(attachPoint));

        // The likelihood that the mutation was lost in the current node under the restricted model
        this->lossAltInCurrentNodeRScore(attachPoint) = subLogProb(this->lossAltInCurrentNodeScore(attachPoint),
                                                                   attachmentScores.hetScore(pN)
                                                                   - attachmentScores.hetScore(attachPoint));
    }
};


#endif //PASSDOWNATTACHMENTSCORES_H
