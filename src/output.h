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

#include "sciphi_config.h"
#include "version.h"

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <iomanip>
#include <ctime>
#include <math.h>

//#include "trees.h"

// Some forward declarations of functions required below

void writeToFile(std::string content, std::string fileName);

std::string getGraphVizFileContentNumbers(int *parents, int n);

void getGraphVizFileContentNames(std::string &output,
                                 Config<MutationTree> &config,
                                 std::vector<unsigned> const &parents,
                                 std::vector<std::string> const &geneNames,
                                 bool attachSamples,
                                 std::vector<std::vector<bool>> const &ancMatrix);

std::string getBestAttachmentString(std::vector<std::vector<bool>> const &ancMatrix,
                                    Config<MutationTree> &config,
                                    std::vector<std::string> geneNames);

void writeTree(Config<SampleTree> const &config);

// This function prints the probabilities of the different mutation models for each cell and mutation to disk.
inline
void
printMutationProbability(Config<SampleTree> const &config,
                         std::string const &outfileName) {

    std::ofstream outFile;
    outFile.open(outfileName);
    outFile << "chrom\tposition\t";
    for (unsigned int i = 0; i < config.getNumSamples(); ++i)
        outFile << config.cellNames[i] << "\t";
    outFile << "\n";
    for (unsigned int j = 0; j < config.getNumMutations(); ++j) {
        outFile << std::get<0>(config.indexToPosition[j]) << "\t" << std::get<1>(config.indexToPosition[j]) << "\t";
        for (unsigned int i = 0; i < config.getNumSamples(); ++i) {
            outFile << config.mutInSampleCounter[i][j].hetScore() << "|" <<
                    config.mutInSampleCounter[i][j].homScore() << "|" <<
                    config.mutInSampleCounter[i][j].lossWildScore() << "|" <<
                    config.mutInSampleCounter[i][j].lossAltRScore() << "|" << 
                    config.mutInSampleCounter[i][j].lcaRScore() << "|" <<
                    config.mutInSampleCounter[i][j].finalScore() << "\t";
        }
        outFile << std::endl;
    }
    outFile.close();
}

inline
void
printRootProbability(Config<SampleTree> const &config,
                         std::string const &outfileName,
                         std::vector<std::queue<double>> rootProbabilities) {

    std::ofstream outFile;
    outFile.open(outfileName);
    for (unsigned int j = 0; j < config.getNumMutations(); ++j) {
        outFile << std::get<0>(config.indexToPosition[j]) << "\t" << std::get<1>(config.indexToPosition[j]) <<"\t";
        while (rootProbabilities[j].size() > 1)
        {
            outFile << rootProbabilities[j].front() <<",";
            rootProbabilities[j].pop();
        }
        outFile << rootProbabilities[j].front() << std::endl;
    }
    outFile.close();
}

inline
void
writeVCFHeader(Config<SampleTree> const &config,
               std::ofstream &outFile) {
    outFile << "##fileformat=VCFv4.1\n";

    // get current date
    
    outFile << "##source=SCIPhIN " << "SCIPhIN v" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << "\n";
    outFile << "##FILTER=<ID=LowQual,Description=\"Low quality\">\n";
    outFile
            << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n";
    outFile << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for alt alleles\">\n";
    outFile
            << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n";
    outFile << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
    outFile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    outFile
            << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n";

    outFile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    if (config.cellInfo != "") {
        for (unsigned i = 0; i < config.getNumSamples(); ++i) {
            outFile << "\t" << config.cellNames[i];
        }
    } else {
        for (unsigned i = 0; i < config.getNumSamples(); ++i) {
            outFile << "\tcell" << i;
        }
    }
    outFile << "\n";
}

inline
bool
checkIfMutationIsUniq(Config<SampleTree> const &config,
                      std::vector<unsigned> &entries,
                      unsigned mutationPos) {
    unsigned numMutatedCells = 0;
    for (unsigned sample = 0; sample < config.getNumSamples(); ++sample) {
        if (config.mutInSampleCounter[sample][mutationPos].finalScore() >= 0.95) {
            ++numMutatedCells;
        }
    }
    return numMutatedCells <= config.uniqTreshold;
}

inline
void
writeVCFMultiEntry(std::ofstream &outFile,
                   Config<SampleTree> const &config,
                   std::vector<unsigned> &entries) {
    for (unsigned alt = 0; alt < entries.size(); ++alt) {
        // chrom
        outFile << std::get<0>(config.indexToPosition[entries[alt]]) << "\t";
        // pos
        outFile << std::get<1>(config.indexToPosition[entries[alt]]) << "\t";
        // id
        outFile << "*\t";
        // ref
        outFile << std::get<2>(config.indexToPosition[entries[alt]]) << "\t";
        // alt
        outFile << std::get<3>(config.indexToPosition[entries[alt]]) << "\t";
        // qual
        outFile << "*\t";
        // filter
        if (checkIfMutationIsUniq(config, entries, entries[alt])) {
            outFile << ".\t";
        } else {
            outFile << "PASS\t";
        }
        // info
        outFile << "DP=";
        unsigned dpAll = 0;
        for (unsigned sample = 0; sample < config.getNumSamples(); ++sample)
            dpAll += std::get<0>(config.getData()[sample][entries[alt]]);
        outFile << dpAll << "\t";
        // format
        outFile << "GT:AD:DP:GQ:PL\t";
        for (unsigned sample = 0; sample < config.getNumSamples(); ++sample) {
            if (config.mutInSampleCounter[sample][entries[alt]].finalScore() > 0.95) {
                outFile << "0/1:";
            } else {
                outFile << "0/0:";
            }
            outFile << std::get<1>(config.getData()[sample][entries[alt]]) << ":";
            outFile << std::get<0>(config.getData()[sample][entries[alt]]) << ":";
            outFile << static_cast<unsigned>(-10 * std::log10(
                    std::max(0.00000001, 1 - config.mutInSampleCounter[sample][entries[alt]].finalScore()))) << "\t";
        }
        outFile << "\n";
    }
}

inline
void
writeVCFEntry(std::ofstream &outFile,
              Config<SampleTree> const &config,
              std::vector<unsigned> &entries) {
    writeVCFMultiEntry(outFile, config, entries);
}


inline
void
writeVCF(Config<SampleTree> const &config,
         std::string const &outfileName) {
    // Create the output file
    std::ofstream outFile;
    outFile.open(outfileName);

    // print the header information
    writeVCFHeader(config, outFile);

    // create a map to sort the mutations according to chromosome and position
    std::map<std::pair<std::string, unsigned>, std::vector<unsigned>> posMap;
    for (unsigned i = 0; i < config.indexToPosition.size(); ++i) {
        std::pair<std::string, unsigned> newPos(std::get<0>(config.indexToPosition[i]),
                                                std::get<1>(config.indexToPosition[i]));
        auto it = posMap.find(newPos);
        if (it != posMap.end()) {
            it->second.push_back(i);
        } else {
            posMap[newPos] = {i};
        }
    }

    for (auto it = posMap.begin(); it != posMap.end(); ++it) {
        writeVCFEntry(outFile, config, it->second);
    }
}

void writeTree(Config<SampleTree> const & config, std::string const & dir)
{
    std::ofstream ofs(dir + "/tree.gv");
    write_graphviz(ofs, config.getTree(), my_label_writer_complete(config.getTree()));
    ofs.close();
}

void writeNucInfo(Config<SampleTree> const & config, std::string const & dir)
{
    std::string makeDir = "mkdir -p " + dir;
    std::system(makeDir.c_str());
    std::ofstream outFile;
    outFile.open(dir + "/nuc.tsv");
    
    outFile << "=numSamples=" << "\n";
    outFile << config.getNumSamples() << "\n";

    outFile << "=params=" << "\n";
    outFile << config.getParam(Config<SampleTree>::E_wildMean) << "\t"
            << config.getSDParam(Config<SampleTree>::E_wildMean) << "\t"
            << config.getSDCountParam(Config<SampleTree>::E_wildMean) << "\t"
            << config.getSDTrialsParam(Config<SampleTree>::E_wildMean) << "\n";
    outFile << config.getParam(Config<SampleTree>::E_wildOverDis) << "\t"
            << config.getSDParam(Config<SampleTree>::E_wildOverDis) << "\t"
            << config.getSDCountParam(Config<SampleTree>::E_wildOverDis) << "\t"
            << config.getSDTrialsParam(Config<SampleTree>::E_wildOverDis) << "\n";
    outFile << config.getParam(Config<SampleTree>::E_mutationOverDis) << "\t"
            << config.getSDParam(Config<SampleTree>::E_mutationOverDis) << "\t"
            << config.getSDCountParam(Config<SampleTree>::E_mutationOverDis) << "\t"
            << config.getSDTrialsParam(Config<SampleTree>::E_mutationOverDis) << "\n";
    outFile << config.getParam(Config<SampleTree>::E_mu) << "\t" << config.getSDParam(Config<SampleTree>::E_mu) << "\t"
            << config.getSDCountParam(Config<SampleTree>::E_mu) << "\t"
            << config.getSDTrialsParam(Config<SampleTree>::E_mu) << "\n";
    outFile << config.getParam(Config<SampleTree>::E_nu) << "\t" << config.getSDParam(Config<SampleTree>::E_nu) << "\t"
            << config.getSDCountParam(Config<SampleTree>::E_nu) << "\t"
            << config.getSDTrialsParam(Config<SampleTree>::E_nu) << "\n";
    outFile << config.getParam(Config<SampleTree>::E_lambdaWildLoss) << "\t" << config.getSDParam(Config<SampleTree>::E_lambdaWildLoss)
            << "\t" << config.getSDCountParam(Config<SampleTree>::E_lambdaWildLoss) << "\t"
            << config.getSDTrialsParam(Config<SampleTree>::E_lambdaWildLoss) << "\n";
    outFile << config.getParam(Config<SampleTree>::E_lambdaMutLoss) << "\t" << config.getSDParam(Config<SampleTree>::E_lambdaMutLoss)
            << "\t" << config.getSDCountParam(Config<SampleTree>::E_lambdaMutLoss) << "\t"
            << config.getSDTrialsParam(Config<SampleTree>::E_lambdaMutLoss) << "\n";
    outFile << config.getParam(Config<SampleTree>::E_parallel) << "\t"
            << config.getSDParam(Config<SampleTree>::E_parallel) << "\t"
            << config.getSDCountParam(Config<SampleTree>::E_parallel) << "\t"
            << config.getSDTrialsParam(Config<SampleTree>::E_parallel) << "\n";

    outFile << "=mutations=" << "\n";
    for (unsigned i = 0; i < config.indexToPosition.size(); ++i) {
        // chrom
        outFile << std::get<0>(config.indexToPosition[i]) << "\t";
        // pos
        outFile << std::get<1>(config.indexToPosition[i]) << "\t";
        // ref
        outFile << std::get<2>(config.indexToPosition[i]) << "\t";
        // alt
        outFile << std::get<3>(config.indexToPosition[i]);

        for (unsigned j = 0; j < config.getCompleteData().size(); ++j) {
            outFile << "\t" << std::get<0>(config.getCompleteData()[j][i]) << "\t"
                    << std::get<1>(config.getCompleteData()[j][i]);
        }
        outFile << "\n";
    }

    outFile << "=background=" << "\n";
    if (config.noiseCounts.cov.size() > 0) {
        outFile << config.noiseCounts.cov[0].first << "\t" << config.noiseCounts.cov[0].second;
        for (unsigned i = 1; i < config.noiseCounts.cov.size(); ++i) {
            outFile << "\t" << config.noiseCounts.cov[i].first << "\t" << config.noiseCounts.cov[i].second;
        }
        outFile << "\n";
    }

    if (config.noiseCounts.sup.size() > 0) {
        outFile << config.noiseCounts.sup[0].first << "\t" << config.noiseCounts.sup[0].second;
        for (unsigned i = 1; i < config.noiseCounts.sup.size(); ++i) {
            outFile << "\t" << config.noiseCounts.sup[i].first << "\t" << config.noiseCounts.sup[i].second;
        }
        outFile << "\n";
    }

    if (config.noiseCounts.covMinusSup.size() > 0) {
        outFile << config.noiseCounts.covMinusSup[0].first << "\t" << config.noiseCounts.covMinusSup[0].second;
        for (unsigned i = 1; i < config.noiseCounts.covMinusSup.size(); ++i) {
            outFile << "\t" << config.noiseCounts.covMinusSup[i].first << "\t"
                    << config.noiseCounts.covMinusSup[i].second;
        }
        outFile << "\n";
    }
    outFile.close();
}

void writeIndex(Config<SampleTree> const & config, std::string const & dir)
{
    writeNucInfo(config, dir);
    writeTree(config, dir);
}

void writeBestIndex(Config<SampleTree> const & config)
{
        writeNucInfo(config, config.bestName);
        writeTree(config, config.bestName);
}

void writeFinalIndex(Config<SampleTree> & config)
{
    writeIndex(config, config.lastName);
}

// Write the mutation to node assignment for the current tree to disk
template<typename TTreeType>
void writeMutToMaxTree(Config<TTreeType> &config) {
    typename Config<TTreeType>::TAttachmentScores::TAttachmentScore scoreSum;
    typename Config<TTreeType>::TAttachmentScores &attachmentScores = config.getTmpAttachmentScore();

    std::ofstream outFile;
    outFile.open(config.mutToMaxName);

    for (unsigned i = 0; i < config.getNumMutations(); ++i) {
        double score = getBestAttachmentScore(config, i);

        scoreSum = AttachmentScore();
        for (unsigned i = 0; i < attachmentScores.size(); ++i) {
            scoreSum.addInRealSpace(attachmentScores[i]);
        }
        for (unsigned i = 0; i < attachmentScores.size(); ++i) {
            outFile << exp(attachmentScores[i].finalScore() - scoreSum.finalScore()) << "\t";
        }
        outFile << std::endl;
    }
    outFile.close();
}
#endif /* OUTPUT_H_ */
