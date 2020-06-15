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
#ifndef READDATA_H
#define READDATA_H

#include <boost/algorithm/string.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/regex.hpp>

#include <cctype>
#include <set>
#include <queue>

#include "sciphi_config.h"
#include "probabilities.h"
#include "noise_counts.h"

#include <dlib/global_optimization.h>

// This functions is used to extract information about 
// sequencing errors. Only if a single cell shows a mutation
// the error rates will be collected.
inline
void updateSeqErrorStats(unsigned &seqErrors,
                         unsigned &seqErrorsCombCov,
                         std::vector<std::array<unsigned, 5>> const &counts) {

    for (unsigned short j = 0; j < 4; ++j) // loop over the nucleotides
    {
        for (unsigned i = 0; i < counts.size(); ++i) // loop over all cells
        {
            seqErrors += counts[i][j];
        }
    }
    for (unsigned i = 0; i < counts.size(); ++i) // update the coverage
    {
        seqErrorsCombCov += counts[i][4];
    }
}

// This function defines a code for the possible sequencing information, such as nucleotides, indels, ...
inline
unsigned charToIndex(char c) {
    switch (std::toupper(c)) {
        case ('A'):
            return 0;
        case ('C'):
            return 1;
        case ('G'):
            return 2;
        case ('T'):
            return 3;
        case ('-'):
        case ('+'):
            return 4;
        case ('^'):
            return 5;
        case ('$'):
            return 6;
        case ('.'):
        case (','):
            return 7;
        case ('*'):
            return 8;
    }

    if (('A' <= c && c <= 'Z') || ('a' <= c && c <= 'z')) {
        return 9;
    }

    std::cerr << "Unknown character \"" << c << "\" in pileup sequence!" << std::endl;
    std::exit(EXIT_FAILURE);

    return 10;
}

inline
char indexToChar(unsigned index) {
    switch (index) {
        case (0):
            return 'A';
        case (1):
            return 'C';
        case (2):
            return 'G';
        case (3):
            return 'T';
    }
    return 'N';
}

// This function is used to skip indels in the mpileup
unsigned skipIndels(std::string const &nucleotides,
                    unsigned currentPos) {
    // if not an indel return
    if (nucleotides[currentPos] != '-' && nucleotides[currentPos] != '+') {
        return currentPos;
    }

    unsigned numIndels = 0;
    unsigned i = currentPos + 1; // skip the '-' or '+'
    for (; i < nucleotides.size(); ++i) {
        // get the number of nucleotieds to skip (size of the indel)
        if (nucleotides[i] >= '0' && nucleotides[i] <= '9') {
            numIndels *= 10;
            numIndels += static_cast<int>(nucleotides[i]) - 48;
        }
            //break if new nucleotide starts
        else {
            break;
        }
    }
    return i + (numIndels - 1);
}

inline
void extractCellNucCountInformation(std::array<unsigned, 5> &counts,
                                    std::string const &nucleotides) {
    //counts = {{0, 0, 0, 0, 0}}; // (a,c,g,t, coverage)
    for (size_t j = 0; j < nucleotides.size(); ++j) // loop over the nucleotides of a cell in the pileup
    {
        unsigned index = charToIndex(nucleotides[j]);
        if (index < 4) {
            ++counts[index];
        } else if (index == 9) {
            --counts[4];
        } else if (index == 4) {
            j = skipIndels(nucleotides, j);
        } else if (index == 5) {
            ++j;
        } else if (index == 6) {
            continue;
        }

    }
}

// This function is used to collect the sequencing information of all cells
inline
void extractSeqInformation(std::array<unsigned, 5> &count,
                           std::vector<std::string> const &splitVec,
                           unsigned position) {
    count = {{0, 0, 0, 0, 0}}; // (a,c,g,t, coverage)
    count[4] = std::stoi(splitVec[position * 3 + 3]);
    extractCellNucCountInformation(count, splitVec[position * 3 + 4]);
}

// This function is used to collect the sequencing information of all cells
inline
void extractSeqInformation(std::vector<std::array<unsigned, 5>> &counts,
                           std::vector<std::string> const &splitVec,
                           std::vector<unsigned> const &positions) {
    for (unsigned i = 0; i < positions.size(); ++i) // loop over the cells
    {
        extractSeqInformation(counts[i], splitVec, positions[i]);
    }
}

// This function updates an array with sequencing error/noise information.
// The array is later used when computing the tree score to account for all
// the non mutated positions.
inline
void addNoiseCounts(std::vector<std::array<unsigned, 5>> &counts,
                    GappedNoiseCounts &gappedNoiseCounts) {
    for (unsigned cellId = 0; cellId < counts.size(); ++cellId) // loop over the cells
    {
        unsigned numErrors = 0;
        for (unsigned nucId = 0; nucId < 4; ++nucId) // check all 4 different nucleotides
        {
            numErrors += counts[cellId][nucId];
        }
        gappedNoiseCounts.add(numErrors, counts[cellId][4]);
    }
}

// TODO: this is deprecated and needs to be removed
// This function sets the nucleotide information for cells with not sufficient coverage to zero
template<typename TTreeType>
void applyCoverageFilterPerCell(std::vector<std::array<unsigned, 5>> &counts,
                                Config<TTreeType> const &config) {
    for (unsigned i = 0; i < counts.size(); ++i) {
        if (counts[i][4] < config.minCoverage) {
            for (unsigned j = 0; j < 5; ++j) {
                counts[i][j] = 0;
            }
        }
    }
}

// Check if the cell passes the support filter
bool passSuppFilter(unsigned altCount,
                    unsigned minSupport) {
    if (altCount >= minSupport) {
        return true;
    }
    return false;
}

// Check if the cell passes the frequency filter
bool passFreqFilter(double altCount,
                    double coverage,
                    double minFreq) {
    if (altCount / coverage >= minFreq) {
        return true;
    }
    return false;
}

// Check if the cell passes the coverage filter
bool passCovFilter(unsigned coverage,
                   unsigned minCov) {
    if (coverage >= minCov) {
        return true;
    }
    return false;
}

// Check if the three filters above are all passed
template<typename TTreeType>
bool applyFilterAcrossCells(std::vector<std::array<unsigned, 5>> const &counts,
                            Config<TTreeType> const &config,
                            unsigned nucId) {
    unsigned numCellsAboveThresh = 0;
    for (unsigned i = 0; i < counts.size(); ++i) {
        if (passCovFilter(counts[i][4], config.minCoverage) &&
            passFreqFilter(counts[i][nucId], counts[i][4], config.minFreq) &&
            passSuppFilter(counts[i][nucId], config.minSupport)) {
            ++numCellsAboveThresh;
            if (numCellsAboveThresh >= config.minNumCellsPassFilter) {
                return true;
            }
        }
    }
    return false;
}

// This function reads a list of positions to exclude from the analysis
inline
void
readExclusionList(std::set<std::tuple<std::string, std::string>> &exMap,
                  std::string const &fileName) {
    std::ifstream inputStream(fileName);
    std::string currLine;                    // line currently processed
    std::vector<std::string> splitVec;  // vector storing the words in currLine dived by a tab
    while (std::getline(inputStream, currLine)) {
        if (currLine[0] != '#') {
            boost::split(splitVec, currLine, boost::is_any_of("\t"));
            exMap.insert(std::tuple<std::string, std::string>(splitVec[0], std::to_string(std::stoi(splitVec[1]))));
        }
    }
}

<<<<<<< HEAD
std::set<std::tuple<std::string, std::string, std::string, std::string>>
readInclusionVCF(std::string const & fileName)
{
    std::set<std::tuple<std::string, std::string, std::string, std::string>> incMap;
    if(fileName.empty()) return incMap;

    std::cout << "Reading inclusion VCF " << fileName << " ..." << std::flush;
    std::ifstream inputStream(fileName);
    if (!inputStream){
        throw std::runtime_error("Could not open inclusionVCF (--inc) " + fileName);
    }
    std::string currLine;                    // line currently processed
    std::vector<std::string> splitVec;  // vector storing the words in currLine dived by a tab
    while (std::getline(inputStream, currLine))
    {
        if (currLine[0] != '#')
        {
            boost::split(splitVec, currLine, boost::is_any_of("\t"));
            incMap.insert(std::make_tuple(splitVec[0], splitVec[1], splitVec[3], splitVec[4]));
        }
    }
    std::cout << " read " << incMap.size() << " variants." << std::endl;

    return incMap;
}

template <typename TTreeType>
void
readCellNames(Config<TTreeType> &config) {
    std::ifstream inputStream(config.cellInfo);
    std::vector<std::string> splitVec;
    std::vector<std::string> splitVecEntry;

    unsigned counter = 0;
    std::string currLine;

    while (getline(inputStream, currLine)) {
        if (currLine != "") {
            boost::split(splitVec, currLine, boost::is_any_of("\t"));
            if (splitVec[1] == "CT" || (config.useNormalCellsInTree && splitVec[1] == "CN")) {
                boost::split(splitVecEntry, splitVec[0], boost::is_any_of("/"));
                config.cellNames.push_back(splitVecEntry.back());
                if (splitVec.size() > 2) {
                    config.cellColours.push_back(splitVec[2]);
                    if (splitVec.size() > 3) {
                        config.cellClusters.push_back(std::stoi(splitVec[3]));
                    } else {
                        config.cellClusters.push_back(1);
                    }
                } else {
                    config.cellColours.push_back("white");
                    config.cellClusters.push_back(1);
                }
            }
        }
    }
}

// Thus function reads cell information, such as experiment type (bulk, single cell), name and cell type
template<typename TTreeType>
void
readCellInformation(std::vector<unsigned> &tumorCellPos,
                    std::vector<unsigned> &normalCellPos,
                    unsigned &tumorBulkPos,
                    unsigned &normalBulkPos,
                    Config<TTreeType> &config) {

    readCellNames(config);

    tumorCellPos.resize(0);
    normalCellPos.resize(0);
    tumorBulkPos = UINT_MAX;
    normalBulkPos = UINT_MAX;

    std::ifstream inputStream(config.cellInfo);

    std::vector<std::string> splitVec;
    unsigned counter = 0;
    std::string currLine;

    while (getline(inputStream, currLine)) {
        if (currLine != "") {
            boost::split(splitVec, currLine, boost::is_any_of("\t"));
            if (splitVec[1] == "BN") {
                if (normalBulkPos != UINT_MAX) {
                    std::cout << "WARNING: Multiple bulk control files provided. Only using the first one!";
                } else {
                    normalBulkPos = counter;
                }
            }
            if (splitVec[1] == "BT") {
            }
            if (splitVec[1] == "CN") {
                normalCellPos.push_back(counter);
            }
            if (splitVec[1] == "CT" || (config.useNormalCellsInTree && splitVec[1] == "CN")) {
                tumorCellPos.push_back(counter);
            }
            ++counter;
        }
    }
}

// This function computes the log of n choose 2
inline
double logNChoose2(unsigned numMut) {
    return log(static_cast<double>(numMut)) + log((static_cast<double>(numMut) - 1.0) / 2.0);
}

// This function computes the log of n choose k
inline
double logNChooseK(unsigned n, unsigned k, double logNChoosekMinusOne) {
    if (k == 0)
        return 0;
    return logNChoosekMinusOne + log(static_cast<double>(n + 1 - k) / static_cast<double>(k));
}

// This function computes the log of n choose k
inline
double logNChooseK(unsigned n, unsigned k) {
    if (k == 0)
        return 0;
    return std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1);
}

// This function computes the log of the summation of several values.
// The input values are log values and exponentiated for the summation.
inline
double sumValuesInLogSpace(std::vector<double>::const_iterator itBegin, std::vector<double>::const_iterator itEnd) {
    std::vector<double>::const_iterator it = itBegin;
    double maxLogValue = *it;
    ++it;
    for (; it != itEnd; ++it) {
        if (maxLogValue < *it) {
            maxLogValue = *it;
        }
    }

    it = itBegin;
    double h1 = 0;
    for (; it != itEnd; ++it) {
        h1 += exp(*it - maxLogValue);
    }

    return log(h1) + maxLogValue;

}

// check if the reference base is known
bool isRefKnown(std::string const & n)
{
    if (n[0] == 'n' || n[0] == 'N')
        return false;
    return true;
}

template<typename TTreeType>
void estimateSeqErrorRate(Config<TTreeType> & config,
    std::set<std::tuple<std::string, std::string>> const & exMap,
    std::set<std::tuple<std::string, std::string>> const & errExMap,
    std::vector<unsigned> const & tumorCellPos,
    std::vector<unsigned> const & normalCellPos,
    std::ifstream & inputStream,
    std::stringstream & inFileHeadBuff)
{
    //estimate the error rate
    std::ifstream inputStream(config.inFileName);
    if (inputStream.fail()) {
        std::cerr << "The pileup file provided does not exist or you lack permission to access it." << std::endl;
        std::exit(1);
    }

    // helper variables
    std::string currLine;
    std::vector<std::string> splitVec;
    std::vector<std::array<unsigned, 5>> tumor_counts(tumorCellPos.size(), {{0, 0, 0, 0, 0}});
    std::vector<std::array<unsigned, 5>> normal_counts(normalCellPos.size(), {{0, 0, 0, 0, 0}});
    unsigned seqErrors = 0;
    unsigned seqErrorsCombCov = 0;
    for (size_t lineNumber = 0; lineNumber < config.errorRateEstLoops && getline(inputStream, currLine); ++lineNumber)
    {
        inFileHeadBuff << currLine << '\n';
        boost::split(splitVec, currLine, boost::is_any_of("\t"));
        
        if(!isRefKnown(splitVec[2]))
        {
            ++config.errorRateEstLoops;
            continue;
        }

        auto itEx = exMap.find(std::tuple<std::string, std::string>(splitVec[0], splitVec[1]));
        auto itErrEx = errExMap.find(std::tuple<std::string, std::string>(splitVec[0], splitVec[1]));
        if (itEx == exMap.end() && itErrEx == errExMap.end()) {
            // Extract the nucelotide information of the tumor cells
            extractSeqInformation(tumor_counts, splitVec, tumorCellPos);
            // Update the sequencing error rate and the coverage
            updateSeqErrorStats(seqErrors, seqErrorsCombCov, tumor_counts);
            // Extract the nucelotide information of the normal cells
            extractSeqInformation(normal_counts, splitVec, normalCellPos);
            // Update the sequencing error rate and the coverage
            updateSeqErrorStats(seqErrors, seqErrorsCombCov, normal_counts);
        } else {
            --lineNumber;
        }
    }
    ++seqErrors; // add pseudo counts

    // set the new estimated error rate
    if (seqErrorsCombCov > 0) {
        config.setParam(Config<TTreeType>::E_wildMean,
                        static_cast<double>(seqErrors) / static_cast<double>(seqErrorsCombCov));
        config.setTmpParam(Config<TTreeType>::E_wildMean, config.getParam(Config<TTreeType>::E_wildMean));
    }
}

// This function prints the current chromosome
void printCurrentChrom(std::string &currentChrom, std::string const &testChrom) {
    if (testChrom != currentChrom) {
        currentChrom = testChrom;
        std::cout << currentChrom << std::endl << std::flush;
    }
}

// This function adjusts for the probability to place the specified number of mutations into a tree, in case
// of tumor cells, or the probability of a random combination of the specified number of mutated cells.
double updateLogH1Temp(double logH1Temp, unsigned numCells, unsigned numMuts, bool tumorCells, double dropOut) {
    if (tumorCells) {
        return logH1Temp + 2 * logNChooseK(numCells, numMuts) - std::log(2 * numMuts - 1) -
               logNChooseK(2 * numCells, 2 * numMuts);
    }
    return logH1Temp + logNChooseK(numCells, numMuts) + std::log(std::pow(1.0 - dropOut, numMuts)) +
           std::log(std::pow(dropOut, numCells - numMuts)) - (1.0 - std::pow(dropOut, numCells));
}

// Check if H1 has a chance or not. This can be done of the peak number of mutated cells times the number of cells is
// stil smaller than H0
bool mustH0Win(double &logH1Max, double logH1Temp, double logNumCells, double logH0) {
    if (logH1Temp >= logH1Max) {
        logH1Max = logH1Temp;
        return false;
    } else {
        if (logH1Max + logNumCells < logH0) {
            return true;
        }
    }
    return false;
}

// Check of current position and next position are within minDist.
bool isInRange(int minDist,
               std::string const &currentChrom,
               std::string const &nextChrom,
               int currentPos,
               int nextPos) {
    if ((currentChrom == nextChrom) && (std::abs(currentPos - nextPos) <= minDist)) {
        return true;
    }
    return false;
}

// This function determines for each candidate locus how many other candidate loci are within a certain distance.
template<typename TTreeType, typename TData>
std::vector<unsigned> getProximity(Config<TTreeType> &config, TData const &data) {
    std::vector<unsigned> proximity(data.size(), 0);
    for (unsigned i = 0; i < data.size(); ++i) {
        int l = i;
        while (l > 0) {
            --l;
            if (isInRange(config.windowSize,
                          std::get<0>(std::get<0>(data[i])),
                          std::get<0>(std::get<0>(data[l])),
                          std::get<1>(std::get<0>(data[i])),
                          std::get<1>(std::get<0>(data[l])))) {
                ++proximity[i];
            } else {
                break;
            }
        }
        unsigned r = i;
        while (r + 1 < data.size()) {
            ++r;
            if (isInRange(config.windowSize,
                          std::get<0>(std::get<0>(data[i])),
                          std::get<0>(std::get<0>(data[r])),
                          std::get<1>(std::get<0>(data[i])),
                          std::get<1>(std::get<0>(data[r])))) {
                ++proximity[i];
            } else {
                break;
            }
        }
    }
    return proximity;
}

// Randomize the order of the candidate loci.
template<typename TEntry>
void randomizeMutPositions(std::vector<TEntry> &data) {
    for (size_t i = 0; i < data.size(); ++i) {
        swap(data[i], data[std::rand() % data.size()]);
    }
}

// Filter out candidate loci which have to many close by other candidates and store the result in the respective data
// structures. Location information is stored in indexToPosition and the nucleotide counts are stored in completeData.
template<typename TTreeType, typename TData>
void insertData(Config<TTreeType> &config, TData const &data) {
    std::string currentChrom = "";
    std::vector<unsigned> proximity = getProximity(config, data);

    TData dataFiltered;

    for (size_t i = 0; i < data.size(); ++i) {
        if (proximity[i] < config.maxMutPerWindow) {
            dataFiltered.push_back(data[i]);
        }
    }

    randomizeMutPositions(dataFiltered);

    for (size_t i = 0; i < dataFiltered.size(); ++i) {
        config.indexToPosition.push_back(std::get<0>(dataFiltered[i]));
        for (size_t j = 0; j < std::get<1>(dataFiltered[i]).size(); ++j) {
            config.getCompleteData()[j].push_back(std::get<1>(dataFiltered[i])[j]);
        }
    }
}

// Check if the mutation has evidence in the control bulk sample and ignore it in that case.
template<typename TTreeType>
bool
passNormalFilter(std::array<unsigned, 5> const &normalCounts, Config<TTreeType> const &config) {
    if (normalCounts[4] >= config.minCovInControlBulk) {
        for (unsigned i = 0; i < 4; ++i) {
            if (normalCounts[i] >= config.maxSupInControlBulk) {
                return false;
            }
        }
        return true;
    }
    return false;
}

//Check if at least one normal cell has enough coverage to make any claims
template<typename TTreeType>
bool
passNormalCellCoverage(std::vector<std::array<unsigned, 5>> const &normalCellCounts,
                       Config<TTreeType> const &config) {
    unsigned maxCov = 0;
    unsigned numMutated = 0;
    for (size_t i = 0; i < normalCellCounts.size(); ++i) {
        if (normalCellCounts[i][4] > maxCov) {
            maxCov = normalCellCounts[i][4];
        }
    }

    if (maxCov < config.minCovNormalCell) {
        return false;
    }

    return true;
}

// Check if the number of mutated control cell is smaller than the threshold.
template<typename TTreeType>
bool
passNormalCellFilter(std::vector<std::array<unsigned, 5>> const &normalCellCounts,
                     unsigned j,
                     Config<TTreeType> const &config) {
    unsigned maxCov = 0;
    unsigned numMutated = 0;
    for (size_t i = 0; i < normalCellCounts.size(); ++i) {
        if (computeRawWildLogScore(config, normalCellCounts[i][j], normalCellCounts[i][4]) <
            computeRawMutLogScore(config, normalCellCounts[i][j], normalCellCounts[i][4])) {
            ++numMutated;
        }
    }

    if (numMutated > config.maxNumberNormalCellMutated) {
        return false;
    }

    return true;
}

// This function computes the probability that 1, 2, 3, ..., all cells are mutated. It does so using a dynamic
// programming approach where the the probability that x cells are mutated uses the information computed when
// computing the probability that x-1 cells are mutated.
// In order to compute the probabilities efficiently two vectors (logProbs and tempLogProbs) are constantly updated
// and swapped.
template<typename TTreeType>
bool computeProbCellsAreMutated(
        Config<TTreeType> const &config,
        std::vector<long double> &logProbs,
        std::vector<long double> &tempLogProbs,
        std::vector<double> &logProbTempValues,
        std::vector<std::array<unsigned, 5>> &counts,
        std::vector<double> &cellsNotMutated,
        std::vector<double> &cellsMutated,
        unsigned currentChar,
        bool tumorCells) {

    unsigned numCells = counts.size();
    double logNumCells = log(numCells);

    // first compute for each cell the probability that it is mutated and the probability the cell is not mutated.
    // Also compute the probability that no cell is mutated.
    tempLogProbs[0] = 0.0;
    for (size_t i = 1; i < counts.size() + 1; ++i) {
        cellsNotMutated[i - 1] = computeRawWildLogScore(config,
                                                        counts[i - 1][currentChar],
                                                        counts[i - 1][4]);
        cellsMutated[i - 1] = computeRawMutLogScore(config,
                                                    counts[i - 1][currentChar],
                                                    counts[i - 1][4]);
        tempLogProbs[i] = tempLogProbs[i - 1] + cellsNotMutated[i - 1];
    }
    swap(logProbs, tempLogProbs);

    // the prior probability of a mutation depend on whether the cell is a normal or tumor cell
    double prior = tumorCells ? config.priorMutationRate : config.priorGermlineRate;

    // Adjust the probability of no mutation with the prior probability of no mutation
    double logH0 = logProbs.back() + log(1.0 - prior);

    logProbTempValues[0] = logH0;

    // compute the probabilitues of observing 1, 2, 3, ... mutations
    double logH1Max = -DBL_MAX; // the current best alternative score
    long double logNOverK = 0;  // helper to efficiently compute nChooseK
    size_t numMut = 1;          // number of mutations currently computed

    for (; numMut <= numCells; ++numMut) {
        double logProbAllPrevCellsMutated = logProbs[numMut - 1];
        double currentCellMutated = cellsMutated[numMut - 1];
        tempLogProbs[numMut] = logProbAllPrevCellsMutated + currentCellMutated;
        for (size_t i = numMut + 1; i < counts.size() + 1; ++i) {
            double previousCellNotMutated = logProbs[i - 1];
            currentCellMutated = cellsMutated[i - 1];
            double previousCellMutated = tempLogProbs[i - 1];
            double currentCellNotMutated = cellsNotMutated[i - 1];
            tempLogProbs[i] = addLogProb(previousCellNotMutated + currentCellMutated,
                                         previousCellMutated + currentCellNotMutated);

        }
        swap(logProbs, tempLogProbs);
        logNOverK = logNChooseK(numCells, numMut, logNOverK);
        double logH1Temp = logProbs.back() + log(prior) - logNOverK;
        logH1Temp = updateLogH1Temp(logH1Temp, numCells, numMut, tumorCells,
                                    (1.0 - config.getParam(Config<TTreeType>::E_mu)) / 2.0);

        // check whether the alternative hypothesis can win
        bool h0Wins = mustH0Win(logH1Max, logH1Temp, logNumCells, logH0);
        logProbTempValues[numMut] = logH1Temp;
        if (h0Wins) {
            return true;
        }
    }
    return false;
}


template<typename TTreeType>
bool readMpileupFile(Config<TTreeType> &config) {
    typedef std::tuple<unsigned, unsigned> TCountType;
    typedef std::tuple<std::string, unsigned, char, char> TPositionInfo;
    typedef std::tuple<TPositionInfo, std::vector<TCountType>, long double> TDataEntry;
    typedef std::vector<TDataEntry> TData;

    // read the cellular information
    std::vector<unsigned> tumorCellPos;
    std::vector<unsigned> normalCellPos;
    unsigned tumorBulkPos;
    unsigned normalBulkPos;
    readCellInformation(tumorCellPos, normalCellPos, tumorBulkPos, normalBulkPos, config);

    // the variables for number of samples
    unsigned numCells = tumorCellPos.size();
    config.setNumSamples(numCells);

    // resize the container holding the nucleotide counts of likely mutations
    config.getCompleteData().resize(numCells);
    config.getData().resize(numCells);

    // determine which positions to skip
    std::set<std::tuple<std::string, std::string>> exMap;
    readExclusionList(exMap, config.exclusionFileName);

    // determine which positions to skip during the error rate estimation
    std::set<std::tuple<std::string, std::string>> errExMap;
    readExclusionList(errExMap, config.mutationExclusionFileName);

    // determine which variants to keep
    auto incMap = readInclusionVCF(config.variantInclusionFileName);

    // init data
    TData data;

    std::array<unsigned, 5> normalBulkCounts;
    // vector to hold the nucleotide information {a,c,g,t,coverage}
    std::vector<std::array<unsigned, 5>> counts(numCells, {{0, 0, 0, 0, 0}});
    // vector to hold the nucleotide information {a,c,g,t,coverage}
    std::vector<std::array<unsigned, 5>> countsNormal(normalCellPos.size(), {{0, 0, 0, 0, 0}});
    // vector to hold the filtered nucleotide information {a,c,g,t,coverage}
    std::vector<std::array<unsigned, 5>> filteredCounts(numCells, {{0, 0, 0, 0, 0}});
    // helper vector to store nucleotide information
    std::vector<TCountType> tempCounts(numCells);
    // probabilities of observing 0, 1, 2, 3, 4 ... mutations
    std::vector<long double> logProbs(numCells + 1, 0);
    // probabilities of observing 0, 1, 2, 3, 4 ... mutations
    std::vector<long double> logProbsNormal(normalCellPos.size() + 1, 0);
    // helper array for probabilities of observing 0, 1, 2, 3, 4 ... mutations
    std::vector<long double> tempLogProbs(numCells + 1, 0);
    // helper array for probabilities of observing 0, 1, 2, 3, 4 ... mutations
    std::vector<long double> tempLogProbsNormal(normalCellPos.size() + 1);
    // helper array for probabilities of observing 0, 1, 2, 3, 4 ... mutations
    std::vector<double> logProbTempValues(numCells + 1);
    // helper array for probabilities of observing 0, 1, 2, 3, 4 ... mutations
    std::vector<double> logProbTempValuesNormal(normalCellPos.size() + 1);

    // helper arrays in which the probabilities of the cells to be mutated and to be bot mutated is stored.
    std::vector<double> cellsNotMutated(numCells);
    std::vector<double> cellsNotMutatedNormal(normalCellPos.size());
    std::vector<double> cellsMutated(numCells);
    std::vector<double> cellsMutatedNormal(normalCellPos.size());

    std::stringstream inFileHeadBuff;
    std::ifstream inputStream(config.inFileName, std::ifstream::in);
    if (!inputStream){
        throw std::runtime_error("Could not open input file " + config.inFileName);
    }

    // Before the mutation probabilities are computed estimate the sequencing error rate
    if (config.estimateSeqErrorRate)
    {
        estimateSeqErrorRate(config, exMap, errExMap, tumorCellPos, normalCellPos, inputStream, inFileHeadBuff);
    }

    std::ifstream inputStream(config.inFileName, std::ifstream::in);
    if (inputStream.fail()) {
        std::cerr << "The pileup file provided does not exist or you lack permission to access it." << std::endl;
        std::exit(1);
    }
    std::string currentChrom = "";
    std::string currLine;
    std::vector<std::string> splitVec;

    // Store the noise counts, aka. the nucleotide counts of the positions not mutated
    GappedNoiseCounts gappedNoiseCounts;

    // Print the parameter values to the console to be able to compare it to the later imputed values
    std::cout.precision(15);
    config.printParameters();

    while (getline(inFileHeadBuff, currLine) || getline(inputStream, currLine))
    {
        // solit the current line into easily accessible chunks
        boost::split(splitVec, currLine, boost::is_any_of("\t"));
        // print the progress
        printCurrentChrom(currentChrom, splitVec[0]);

        if(!isRefKnown(splitVec[2]))
        {
            continue;
        }

        // check if the current pos is to be excluded
        auto it = exMap.find(std::tuple<std::string, std::string>(splitVec[0], splitVec[1]));
        if (it == exMap.end()) {

            // Extract the nucleotide information from the normal bulk
            // If there is enough evidence that there is a germline mutation skip the position
            if (normalBulkPos != UINT_MAX) {
                extractSeqInformation(normalBulkCounts, splitVec, normalBulkPos);
                if (!passNormalFilter(normalBulkCounts, config)) {
                    continue;
                }
            }

            // Extract the nucleotide information from the normal cells and check if there is enough coverage in them
            if (normalCellPos.size() > 0) {
                extractSeqInformation(countsNormal, splitVec, normalCellPos);
                if (!passNormalCellCoverage(countsNormal, config)) {
                    continue;
                }
            }

            // Extract the nucleotide information from the tumor cells
            extractSeqInformation(counts, splitVec, tumorCellPos);

            filteredCounts = counts;
            applyCoverageFilterPerCell(filteredCounts, config); // TODO: this should not be used here

            bool positionMutated = false;
            for (unsigned short altAlleleIdx = 0; altAlleleIdx < 4; ++altAlleleIdx)
            {
                if(altAlleleIdx == charToIndex(splitVec[2][0]))
                {
                    continue;
                }
                bool h0Wins = !applyFilterAcrossCells(filteredCounts, config, altAlleleIdx);
              
                if(normalCellPos.size() > 0)
                {
                    // use the normal cell filter
                    if (config.normalCellFilter == 1)
                    {
                        if (!passNormalCellFilter(countsNormal, altAlleleIdx, config))
                        {
                            positionMutated = true;
                            h0Wins = true;
                            continue;
                        }
                    }
                    // Use the normal cell filter in which the probability if 1, 2, 3, ... cells being mutated is
                    // combined.
                    else if (config.normalCellFilter == 2)
                    {
                        bool h0WinsNormal = computeProbCellsAreMutated(config, 
                                logProbsNormal, 
                                tempLogProbsNormal, 
                                logProbTempValuesNormal, 
                                countsNormal,
                                cellsNotMutatedNormal,
                                cellsMutatedNormal,
                                altAlleleIdx,
                                false);

                        if (!h0WinsNormal)
                        {
                            // if the result is not very clear check if there is enough support from at least
                            // maxNumberNormalCellMutated cells
                            double logH0Normal = sumValuesInLogSpace(logProbTempValuesNormal.begin(),logProbTempValuesNormal.begin() + config.maxNumberNormalCellMutated + 1);
                            double logH1Normal = sumValuesInLogSpace(logProbTempValuesNormal.begin() + config.maxNumberNormalCellMutated + 1, logProbTempValuesNormal.end());

                            if (logH0Normal < logH1Normal) 
                            {
                                positionMutated = true;
                                h0Wins = true;
                                continue;
                            }
                        }
                    }
                }

                double logH0 = -DBL_MAX;
                double logH1 = -DBL_MAX;
                // h0Wins and positionMutated are true if there is enough evidence that the normal sample(s) is/are
                // mutated. positionMutated is set to true, such that the locus is not used for the noise counts
                if(!h0Wins)
                {
                    h0Wins = computeProbCellsAreMutated(config, 
                            logProbs, 
                            tempLogProbs, 
                            logProbTempValues, 
                            filteredCounts,
                            cellsNotMutated,
                            cellsMutated,
                            altAlleleIdx,
                            true);
                }

                // if h0Wins set the log probabilities of the null and alternative hypothesis appropriately.
                if (h0Wins) {
                    logH1 = -DBL_MAX;
                    logH0 = DBL_MAX;
                } else {
                    // if the result is not very clear check if there is enough support from at least
                    // numCellWithMutationMin cells
                    logH0 = sumValuesInLogSpace(logProbTempValues.begin(),
                                                logProbTempValues.begin() + config.numCellWithMutationMin);
                    logH1 = sumValuesInLogSpace(logProbTempValues.begin() + config.numCellWithMutationMin,
                                                logProbTempValues.end());
                }

                // If the probabilty that the locus is mutated is higher than it is not (taking all filter criteria
                // into account) check if the signal does not come from systematic sequencing error.
                if (logH1 > logH0 || incMap.count(std::make_tuple(splitVec[0], splitVec[1], splitVec[2], std::string(1, indexToChar(altAlleleIdx)))) != 0)
                {
                    positionMutated = true;

                    // Only look at loci which are likely mutated
                    static std::vector<std::pair<unsigned, unsigned>> testCounts;
                    testCounts.resize(0);
                    for (size_t cell = 0; cell < counts.size(); ++cell)
                    {
                        if(cellsNotMutated[cell] < cellsMutated[cell])
                        {
                            testCounts.push_back({counts[cell][altAlleleIdx], counts[cell][4]});
                        }
                    }

                    // In the followinf two beta binomials are learned.
                    // This first one learnes the overdispersion but the mean is fixed at a user defined value (0.25
                    // per default). This beta-binomial is representing the sequnecing error model
                    dlib::matrix<double, 0, 1> startingPointMeanOverDis = {0.05, 5.0};
                    dlib::matrix<double, 0, 1> dLibMinMeanOverDis = {0.001, 0.1};
                    dlib::matrix<double, 0, 1> dLibMaxMeanOverDis = {0.999, 10000.0};
                    OptimizeBetaBinMeanOverDis optBetaBinMeanOverDis(testCounts);
                    OptimizeBetaBinMeanOverDisDerivates optBetaBinDerMeanOverDis(testCounts);
                    double resultMeanOverDis = dlib::find_max_box_constrained(
                            dlib::bfgs_search_strategy(),  // Use BFGS search algorithm
                            // Stop when the change in rosen() is less than 1e-7
                            dlib::objective_delta_stop_strategy(1e-7),
                            optBetaBinMeanOverDis,
                            optBetaBinDerMeanOverDis,
                            startingPointMeanOverDis,
                            dLibMinMeanOverDis,
                            dLibMaxMeanOverDis);

                    // A second beta-binomial with free mean and overdispersion is learned. This one represents the
                    // mutation model
                    dlib::matrix<double, 0, 1> startingPointOverDis = {2.0};
                    dlib::matrix<double, 0, 1> dLibMinOverDis = {0.1};
                    dlib::matrix<double, 0, 1> dLibMaxOverDis = {10000.0};
                    OptimizeBetaBinOverDis optBetaBinOverDis(testCounts, config.meanFilter);
                    OptimizeBetaBinOverDisDerivates optBetaBinDerOverDis(testCounts, config.meanFilter);
                    double resultOverDis = dlib::find_max_box_constrained(
                            dlib::bfgs_search_strategy(),  // Use BFGS search algorithm
                            // Stop when the change in rosen() is less than 1e-7
                            dlib::objective_delta_stop_strategy(1e-7),
                            optBetaBinOverDis,
                            optBetaBinDerOverDis,
                            startingPointOverDis,
                            dLibMinOverDis,
                            dLibMaxOverDis);

                    // Compute the Chi-squared distributed p-Value
                    double pValue;
                    if (resultOverDis > resultMeanOverDis) // if true the results are within the optimization limit
                    {
                        pValue = 1;
                    } else {
                        double testStat = -2 * (resultOverDis - resultMeanOverDis);
                        boost::math::chi_squared mydist(1);
                        pValue = 1 - boost::math::cdf(mydist, testStat);
                    }

                    // Discard the mutation if its mean significantly deviates from the noise mean AND is smaller
                    // than the noise mean
                    if(pValue > 0.05 || startingPointMeanOverDis(0) >= config.meanFilter || incMap.count(std::make_tuple(splitVec[0], splitVec[1], splitVec[2], std::string(1, indexToChar(altAlleleIdx)))) != 0)
                    {
                        for (size_t cell = 0; cell < counts.size(); ++cell)
                        {
                            std::get<0>(tempCounts[cell]) = counts[cell][4];
                            std::get<1>(tempCounts[cell]) = counts[cell][altAlleleIdx];
                        }
                
                        data.push_back(TDataEntry(TPositionInfo(splitVec[0], std::stoi(splitVec[1]), splitVec[2][0], indexToChar(altAlleleIdx)), tempCounts, logH1/logH0));
                    }
                }
            }
            if (!positionMutated) {
                addNoiseCounts(counts, gappedNoiseCounts);
            }
        }
    }

    config.noiseCounts = NoiseCounts(gappedNoiseCounts);
    computeNoiseScore(config);

    insertData(config, data);
    writeNucInfo(config, config.bestName);

    return 0;
}

template<typename TTreeType>
void getData(Config<TTreeType> &config) {
    readMpileupFile(config);
}

template<typename TTreeType>
void readGraph(Config<TTreeType> &config) {
    std::ifstream inFile;
    std::get<1>(config.dataUsageRate) = 1.0;
    inFile.open(config.loadName + "/tree.gv");

    std::string line;
    std::vector<std::string> splitVec;
    std::vector<int> sampleIds;

    while (std::getline(inFile, line)) {
        boost::algorithm::split_regex(splitVec, line, boost::regex("->"));
        if (splitVec.size() > 1) {
            boost::algorithm::split_regex(splitVec, line, boost::regex("->|[ ;]"));
            boost::add_edge(std::stoi(splitVec[0]), std::stoi(splitVec[1]), config.getTree());
        }
        boost::algorithm::split_regex(splitVec, line, boost::regex("[\[\"]"));
        if (splitVec.size() == 4) {
            if (std::stoi(splitVec[0]) != sampleIds.size()) {
                std::cout << "Problem reading the tree." << std::endl;
                return;
            }
            sampleIds.push_back(std::stoi(splitVec[2]));
        }
    }

    unsigned numVertices = num_vertices(config.getTree());
    for (unsigned i = numVertices / 2 - 1; i < numVertices - 1; ++i) {
        config.getTree()[i].sample = sampleIds[i] - (numVertices / 2 - 1);
    }
}

template<typename TTreeType>
void readNucInfo(Config<TTreeType> &config) {
    std::ifstream inFile;
    std::get<1>(config.dataUsageRate) = 1.0;
    inFile.open(config.loadName + "/nuc.tsv");

    std::string line;
    std::vector<std::string> splitVec;
    std::getline(inFile, line);
    if (line == "=numSamples=") {
        std::getline(inFile, line);
        config.setNumSamples(std::stod(line));
    }

    std::getline(inFile, line);
    if (line == "=params=") {
        std::getline(inFile, line);
        boost::split(splitVec, line, boost::is_any_of("\t"));
        config.setParam(Config<SampleTree>::E_wildMean, std::stod(splitVec[0]));
        config.setSDParam(Config<SampleTree>::E_wildMean, std::stod(splitVec[1]));
        config.setSDCountParam(Config<SampleTree>::E_wildMean, std::stoi(splitVec[2]));
        config.setSDTrialsParam(Config<SampleTree>::E_wildMean, std::stoi(splitVec[3]));

        std::getline(inFile, line);
        boost::split(splitVec, line, boost::is_any_of("\t"));
        config.setParam(Config<SampleTree>::E_wildOverDis, std::stod(splitVec[0]));
        config.setSDParam(Config<SampleTree>::E_wildOverDis, std::stod(splitVec[1]));
        config.setSDCountParam(Config<SampleTree>::E_wildOverDis, std::stoi(splitVec[2]));
        config.setSDTrialsParam(Config<SampleTree>::E_wildOverDis, std::stoi(splitVec[3]));

        std::getline(inFile, line);
        boost::split(splitVec, line, boost::is_any_of("\t"));
        config.setParam(Config<SampleTree>::E_mutationOverDis, std::stod(splitVec[0]));
        config.setSDParam(Config<SampleTree>::E_mutationOverDis, std::stod(splitVec[1]));
        config.setSDCountParam(Config<SampleTree>::E_mutationOverDis, std::stoi(splitVec[2]));
        config.setSDTrialsParam(Config<SampleTree>::E_mutationOverDis, std::stoi(splitVec[3]));

        std::getline(inFile, line);
        boost::split(splitVec, line, boost::is_any_of("\t"));
        config.setParam(Config<SampleTree>::E_mu, std::stod(splitVec[0]));
        config.setSDParam(Config<SampleTree>::E_mu, std::stod(splitVec[1]));
        config.setSDCountParam(Config<SampleTree>::E_mu, std::stoi(splitVec[2]));
        config.setSDTrialsParam(Config<SampleTree>::E_mu, std::stoi(splitVec[3]));

        std::getline(inFile, line);
        boost::split(splitVec, line, boost::is_any_of("\t"));
        config.setParam(Config<SampleTree>::E_nu, std::stod(splitVec[0]));
        config.setSDParam(Config<SampleTree>::E_nu, std::stod(splitVec[1]));
        config.setSDCountParam(Config<SampleTree>::E_nu, std::stoi(splitVec[2]));
        config.setSDTrialsParam(Config<SampleTree>::E_nu, std::stoi(splitVec[3]));

        std::getline(inFile, line);
        boost::split(splitVec, line, boost::is_any_of("\t"));
        config.setParam(Config<SampleTree>::E_lambda, std::stod(splitVec[0]));
        config.setSDParam(Config<SampleTree>::E_lambda, std::stod(splitVec[1]));
        config.setSDCountParam(Config<SampleTree>::E_lambda, std::stoi(splitVec[2]));
        config.setSDTrialsParam(Config<SampleTree>::E_lambda, std::stoi(splitVec[3]));

        std::getline(inFile, line);
        boost::split(splitVec, line, boost::is_any_of("\t"));
        config.setParam(Config<SampleTree>::E_parallel, std::stod(splitVec[0]));
        config.setSDParam(Config<SampleTree>::E_parallel, std::stod(splitVec[1]));
        config.setSDCountParam(Config<SampleTree>::E_parallel, std::stoi(splitVec[2]));
        config.setSDTrialsParam(Config<SampleTree>::E_parallel, std::stoi(splitVec[3]));
    }

    std::getline(inFile, line);
    if (line == "=mutations=") {
        config.data.resize(config.getNumSamples());
        while (getline(inFile, line)) {
            if (line == "=background=") {
                break;
            }

            boost::split(splitVec, line, boost::is_any_of("\t"));
            config.indexToPosition.push_back({splitVec[0], std::stoul(splitVec[1]), splitVec[2][0], splitVec[3][0]});
            for (unsigned i = 4; i < splitVec.size(); i += 2) {
                config.data[(i - 4) / 2].push_back({std::stoul(splitVec[i]), std::stoul(splitVec[i + 1])});
            }
        }
        config.completeData = config.data;
    }

    {
        if (std::getline(inFile, line)) {
            boost::split(splitVec, line, boost::is_any_of("\t"));
            for (unsigned i = 0; i < splitVec.size(); i += 2) {
                config.noiseCounts.cov.push_back({std::stoul(splitVec[i]), std::stoul(splitVec[i + 1])});
                config.noiseCounts.numPos += config.noiseCounts.cov.back().second;
            }

            std::getline(inFile, line);
            boost::split(splitVec, line, boost::is_any_of("\t"));
            for (unsigned i = 0; i < splitVec.size(); i += 2) {
                config.noiseCounts.sup.push_back({std::stoul(splitVec[i]), std::stoul(splitVec[i + 1])});
            }

            std::getline(inFile, line);
            boost::split(splitVec, line, boost::is_any_of("\t"));
            for (unsigned i = 0; i < splitVec.size(); i += 2) {
                config.noiseCounts.covMinusSup.push_back({std::stoul(splitVec[i]), std::stoul(splitVec[i + 1])});
            }
        }
    }
}

#endif

