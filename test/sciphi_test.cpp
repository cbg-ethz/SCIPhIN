#define BOOST_TEST_MODULE SciphiTest
#define BOOST_TEST_DYN_LINK

#include <math.h>
#include <limits>

#include <boost/test/unit_test.hpp>

#include "src/sciphi_config.h"
#include "src/mcmc.h"
#include "src/scoreTree.h"
#include "src/readData.h"
#include "src/probabilities.h"

inline
void createInitialTreeTest(Config<SampleTree> & config)
{
    unsigned i = 1;
    for (; i < config.getNumSamples() - 1; ++i)
    {
        add_edge((i-1)/2, i, config.getTree());
        config.getTree()[i].sample = -1;
    }
    for (; i < 2 * config.getNumSamples() - 1; ++i)
    {
        add_edge((i-1)/2, i, config.getTree());
        config.getTree()[i].sample = i - (config.getNumSamples() - 1);
    }   

    // add artificial node
    add_edge(2 * config.getNumSamples() - 1, 0, config.getTree());

    //for (unsigned i = 0; i < 10 * config.getNumSamples(); ++i)
    //{
    //    pruneAndReAttach(config);
    //}

    config.setTmpTree(config.getTree());
}

void initTree1(Config<SampleTree> &config) {
    /*
    *                                          9
    *                                          |
    *                                          0
    *                           ___________________________
    *                          |                           |
    *                          1                           8
    *               _______________________                0.6
    *              |                       |               0.5
    *              2                       7               0.1
    *       _______________                0.4
    *      |               |               0.6
    *      3               6               0.2
    *   _______            0.3
    *  |       |           0.7
    *  4       5           0.3
    *  0.1     0.2
    *  0.9     0.8
    *  0.5     0.4
    */
    config.setNumSamples(5);

    add_edge(0, 1, config.getTree());
    config.getTree()[0].sample = -1;
    config.getTree()[1].sample = -1;
    add_edge(0, 8, config.getTree());
    config.getTree()[8].sample = 4;
    add_edge(1, 2, config.getTree());
    config.getTree()[3].sample = -1;
    add_edge(1, 7, config.getTree());
    config.getTree()[7].sample = 3;
    add_edge(2, 3, config.getTree());
    config.getTree()[2].sample = -1;
    add_edge(2, 6, config.getTree());
    config.getTree()[6].sample = 2;
    add_edge(3, 4, config.getTree());
    config.getTree()[4].sample = 0;
    add_edge(3, 5, config.getTree());
    config.getTree()[5].sample = 1;
    add_edge(9, 0, config.getTree());

    std::get<0>(config.logScores).resizeNumCells(5);
    std::get<0>(config.logScores).resizeNumMuts(1);
    std::get<0>(config.logScores).wtScore(0, 0) = std::log(0.1);
    std::get<0>(config.logScores).hetScore(0, 0) = std::log(0.9);
    std::get<0>(config.logScores).homScore(0, 0) = std::log(0.5);

    std::get<0>(config.logScores).wtScore(1, 0) = std::log(0.2);
    std::get<0>(config.logScores).hetScore(1, 0) = std::log(0.8);
    std::get<0>(config.logScores).homScore(1, 0) = std::log(0.4);

    std::get<0>(config.logScores).wtScore(2, 0) = std::log(0.3);
    std::get<0>(config.logScores).hetScore(2, 0) = std::log(0.7);
    std::get<0>(config.logScores).homScore(2, 0) = std::log(0.3);

    std::get<0>(config.logScores).wtScore(3, 0) = std::log(0.4);
    std::get<0>(config.logScores).hetScore(3, 0) = std::log(0.6);
    std::get<0>(config.logScores).homScore(3, 0) = std::log(0.2);

    std::get<0>(config.logScores).wtScore(4, 0) = std::log(0.6);
    std::get<0>(config.logScores).hetScore(4, 0) = std::log(0.5);
    std::get<0>(config.logScores).homScore(4, 0) = std::log(0.1);

    config.noiseScore = -1000;
    config.data.resize(5);
    config.data[0].resize(1);
    config._tmpAttachmentScore.resize(9);
    config.setParam(Config<SampleTree>::E_nu, 0.1);
    config.setParam(Config<SampleTree>::E_lambdaWildLoss, 0.2);
    config.setParam(Config<SampleTree>::E_lambdaMutLoss, 0.25);
    config.setParam(Config<SampleTree>::E_parallel, 0.3);
}

void initTree2(Config<SampleTree> &config) {
    /*
     *                                          15
     *                                          |
     *                                          0
     *                           ___________________________
     *                          |                           |
     *                          1                           2
     *               _______________________             _______
     *              |                       |           |       |
     *              3                       4           13      14
     *       _______________             _______        0.7     0.8
     *      |               |           |       |       0.3     0.2
     *      5               6           11      12      0.6     0.5
     *   _______         _______        0.5     0.6
     *  |       |       |       |       0.4     0.4
     *  7       8       9       10      0.8     0.7
     *  0.1     0.2     0.3     0.4
     *  0.9     0.8     0.7     0.5
     *  0.5     0.4     0.3     0.2
     */
    config.setNumSamples(8);

    add_edge(0, 1, config.getTree());
    config.getTree()[0].sample = -1;
    config.getTree()[1].sample = -1;
    add_edge(0, 2, config.getTree());
    config.getTree()[2].sample = -1;
    add_edge(1, 3, config.getTree());
    config.getTree()[3].sample = -1;
    add_edge(1, 4, config.getTree());
    config.getTree()[4].sample = -1;
    add_edge(3, 5, config.getTree());
    config.getTree()[5].sample = -1;
    add_edge(3, 6, config.getTree());
    config.getTree()[6].sample = -1;

    add_edge(5, 7, config.getTree());
    config.getTree()[7].sample = 0;
    add_edge(5, 8, config.getTree());
    config.getTree()[8].sample = 1;

    add_edge(6, 9, config.getTree());
    config.getTree()[9].sample = 2;
    add_edge(6, 10, config.getTree());
    config.getTree()[10].sample = 3;

    add_edge(4, 11, config.getTree());
    config.getTree()[11].sample = 4;
    add_edge(4, 12, config.getTree());
    config.getTree()[12].sample = 5;

    add_edge(2, 13, config.getTree());
    config.getTree()[13].sample = 6;
    add_edge(2, 14, config.getTree());
    config.getTree()[14].sample = 7;

    add_edge(15, 0, config.getTree());

    std::get<0>(config.logScores).resizeNumCells(8);
    std::get<0>(config.logScores).resizeNumMuts(1);
    std::get<0>(config.logScores).wtScore(0, 0) = std::log(0.1);
    std::get<0>(config.logScores).hetScore(0, 0) = std::log(0.9);
    std::get<0>(config.logScores).homScore(0, 0) = std::log(0.5);

    std::get<0>(config.logScores).wtScore(1, 0) = std::log(0.2);
    std::get<0>(config.logScores).hetScore(1, 0) = std::log(0.8);
    std::get<0>(config.logScores).homScore(1, 0) = std::log(0.4);

    std::get<0>(config.logScores).wtScore(2, 0) = std::log(0.3);
    std::get<0>(config.logScores).hetScore(2, 0) = std::log(0.7);
    std::get<0>(config.logScores).homScore(2, 0) = std::log(0.3);

    std::get<0>(config.logScores).wtScore(3, 0) = std::log(0.4);
    std::get<0>(config.logScores).hetScore(3, 0) = std::log(0.5);
    std::get<0>(config.logScores).homScore(3, 0) = std::log(0.2);

    std::get<0>(config.logScores).wtScore(4, 0) = std::log(0.5);
    std::get<0>(config.logScores).hetScore(4, 0) = std::log(0.4);
    std::get<0>(config.logScores).homScore(4, 0) = std::log(0.8);

    std::get<0>(config.logScores).wtScore(5, 0) = std::log(0.6);
    std::get<0>(config.logScores).hetScore(5, 0) = std::log(0.4);
    std::get<0>(config.logScores).homScore(5, 0) = std::log(0.7);

    std::get<0>(config.logScores).wtScore(6, 0) = std::log(0.7);
    std::get<0>(config.logScores).hetScore(6, 0) = std::log(0.3);
    std::get<0>(config.logScores).homScore(6, 0) = std::log(0.6);

    std::get<0>(config.logScores).wtScore(7, 0) = std::log(0.8);
    std::get<0>(config.logScores).hetScore(7, 0) = std::log(0.2);
    std::get<0>(config.logScores).homScore(7, 0) = std::log(0.5);

    config.noiseScore = -1000;
    config.data.resize(8);
    config.data[0].resize(1);
    config._tmpAttachmentScore.resize(15);
    config.setParam(Config<SampleTree>::E_nu, 0.1);
    config.setParam(Config<SampleTree>::E_lambdaWildLoss, 0.2);
    config.setParam(Config<SampleTree>::E_lambdaMutLoss, 0.25);
    config.setParam(Config<SampleTree>::E_parallel, 0.3);
}

void initTree2_parallel(Config<SampleTree> &config) {
    
    config.setNumSamples(16);

    add_edge(0, 1, config.getTree());
    config.getTree()[0].sample = -1;
    config.getTree()[1].sample = -1;
    add_edge(0, 2, config.getTree());
    config.getTree()[2].sample = -1;
    add_edge(1, 3, config.getTree());
    config.getTree()[3].sample = -1;
    add_edge(1, 4, config.getTree());
    config.getTree()[4].sample = -1;
    add_edge(2, 5, config.getTree());
    config.getTree()[5].sample = -1;
    add_edge(2, 6, config.getTree());
    config.getTree()[6].sample = -1;
    add_edge(3, 7, config.getTree());
    config.getTree()[7].sample = -1;
    add_edge(3, 8, config.getTree());
    config.getTree()[8].sample = -1;
    add_edge(4, 9, config.getTree());
    config.getTree()[9].sample = -1;
    add_edge(4, 10, config.getTree());
    config.getTree()[10].sample = -1;
    add_edge(5, 11, config.getTree());
    config.getTree()[11].sample = -1;
    add_edge(5, 12, config.getTree());
    config.getTree()[12].sample = -1;
    add_edge(6, 13, config.getTree());
    config.getTree()[13].sample = -1;
    add_edge(6, 14, config.getTree());
    config.getTree()[14].sample = -1;

    add_edge(7, 15, config.getTree());
    config.getTree()[15].sample = 0;
    add_edge(7, 16, config.getTree());
    config.getTree()[16].sample = 1;
    add_edge(8, 17, config.getTree());
    config.getTree()[17].sample = 2;
    add_edge(8, 18, config.getTree());
    config.getTree()[18].sample = 3;
    add_edge(9, 19, config.getTree());
    config.getTree()[19].sample = 4;
    add_edge(9, 20, config.getTree());
    config.getTree()[20].sample = 5;
    add_edge(10, 21, config.getTree());
    config.getTree()[21].sample = 6;
    add_edge(10, 22, config.getTree());
    config.getTree()[22].sample = 7;
    add_edge(11, 23, config.getTree());
    config.getTree()[23].sample = 8;
    add_edge(11, 24, config.getTree());
    config.getTree()[24].sample = 9;
    add_edge(12, 25, config.getTree());
    config.getTree()[25].sample = 10;
    add_edge(12, 26, config.getTree());
    config.getTree()[26].sample = 11;
    add_edge(13, 27, config.getTree());
    config.getTree()[27].sample = 12;
    add_edge(13, 28, config.getTree());
    config.getTree()[28].sample = 13;
    add_edge(14, 29, config.getTree());
    config.getTree()[29].sample = 14;
    add_edge(14, 30, config.getTree());
    config.getTree()[30].sample = 15;

    add_edge(31, 0, config.getTree());

    std::get<0>(config.logScores).resizeNumCells(16);
    std::get<0>(config.logScores).resizeNumMuts(1);
    for (int i = 0; i < 16; ++i)
    {
        std::get<0>(config.logScores).wtScore(i, 0) = std::log(0.1);
        std::get<0>(config.logScores).hetScore(i, 0) = std::log(0.5);
        std::get<0>(config.logScores).homScore(i, 0) = std::log(0.005);
    }

    config.noiseScore = -1000;
    config.data.resize(16);
    config.data[0].resize(1);
    config._tmpAttachmentScore.resize(31);
    config.setParam(Config<SampleTree>::E_nu, 0.1);
    config.setParam(Config<SampleTree>::E_lambdaWildLoss, 0.2);
    config.setParam(Config<SampleTree>::E_lambdaMutLoss, 0.25);
    config.setParam(Config<SampleTree>::E_parallel, 0.3);
}
void initTree2_1(Config<SampleTree> &config) {
    /*
     *                                          15
     *                                          |
     *                                          0
     *                           ___________________________
     *                          |                           |
     *                          1                           2
     *               _______________________             _______
     *              |                       |           |       |
     *              3                       4           13      14
     *       _______________             _______        1       1
     *      |               |           |       |       0.65     0.9
     *      5               6           11      12      0.933    1.2667
     *   _______         _______        1       1
     *  |       |       |       |       0.15     0.4
     *  7       8       9       10      0.267    0.6
     *  1       1       1       1
     *  -0.85   -0.6    -0.35   -0.1
     *  -1.67   -0.73   -0.4    -0.06667
     */
    config.setNumSamples(8);

    add_edge(0, 1, config.getTree());
    config.getTree()[0].sample = -1;
    config.getTree()[1].sample = -1;
    add_edge(0, 2, config.getTree());
    config.getTree()[2].sample = -1;
    add_edge(1, 3, config.getTree());
    config.getTree()[3].sample = -1;
    add_edge(1, 4, config.getTree());
    config.getTree()[4].sample = -1;
    add_edge(3, 5, config.getTree());
    config.getTree()[5].sample = -1;
    add_edge(3, 6, config.getTree());
    config.getTree()[6].sample = -1;

    add_edge(5, 7, config.getTree());
    config.getTree()[7].sample = 0;
    add_edge(5, 8, config.getTree());
    config.getTree()[8].sample = 1;

    add_edge(6, 9, config.getTree());
    config.getTree()[9].sample = 2;
    add_edge(6, 10, config.getTree());
    config.getTree()[10].sample = 3;

    add_edge(4, 11, config.getTree());
    config.getTree()[11].sample = 4;
    add_edge(4, 12, config.getTree());
    config.getTree()[12].sample = 5;

    add_edge(2, 13, config.getTree());
    config.getTree()[13].sample = 6;
    add_edge(2, 14, config.getTree());
    config.getTree()[14].sample = 7;

    add_edge(15, 0, config.getTree());

    std::get<0>(config.logScores).resizeNumCells(8);
    std::get<0>(config.logScores).resizeNumMuts(1);

    for (unsigned i = 0; i < 8; ++i)
    {
        std::get<0>(config.logScores).wtScore(i, 0) = 0;
        std::get<0>(config.logScores).hetScore(i, 0) = (i+1)/4.0 - 1.1;
        std::get<0>(config.logScores).homScore(i, 0) = (i+1)/3.0 - 1.4;
        std::cout << std::get<0>(config.logScores).hetScore(0, 0) << std::endl;
        std::cout << std::get<0>(config.logScores).homScore(0, 0) << std::endl;
    }

    config.noiseScore = -1000;
    config.data.resize(8);
    config.data[0].resize(1);
    config._tmpAttachmentScore.resize(15);
    config.setParam(Config<SampleTree>::E_nu, 0.15);
    config.setParam(Config<SampleTree>::E_lambdaWildLoss, 0.2);
    config.setParam(Config<SampleTree>::E_lambdaMutLoss, 0.25);
    config.setParam(Config<SampleTree>::E_parallel, 0.3);
}

void initTree3(Config<SampleTree> &config) {
    /*
    *                                          11
    *                                          |
    *                                          0
    *                           ___________________________
    *                          |                           |
    *                          1                           2
    *               _______________________           ____________
    *              |                       |         |            |
    *              3                       8         9            10
    *       _______________                0.4       0.5          0.6
    *      |               |               0.6       0.4          0.3
    *      4               7               0.2       0.1          0.7
    *   _______            0.3
    *  |       |           0.7
    *  5       6           0.3
    *  0.1     0.2
    *  0.9     0.8
    *  0.5     0.4
    */
    config.setNumSamples(6);

    add_edge(0, 1, config.getTree());
    config.getTree()[0].sample = -1;
    config.getTree()[1].sample = -1;
    add_edge(0, 2, config.getTree());
    config.getTree()[2].sample = -1;
    add_edge(1, 3, config.getTree());
    config.getTree()[3].sample = -1;
    add_edge(1, 8, config.getTree());
    config.getTree()[8].sample = 3;
    add_edge(3, 4, config.getTree());
    config.getTree()[4].sample = -1;
    add_edge(3, 7, config.getTree());
    config.getTree()[7].sample = 2;

    add_edge(4, 5, config.getTree());
    config.getTree()[5].sample = 0;
    add_edge(4, 6, config.getTree());
    config.getTree()[6].sample = 1;
    add_edge(2, 9, config.getTree());
    config.getTree()[9].sample = 4;
    add_edge(2, 10, config.getTree());
    config.getTree()[10].sample = 5;

    add_edge(11, 0, config.getTree());

    std::get<0>(config.logScores).resizeNumCells(6);
    std::get<0>(config.logScores).resizeNumMuts(1);
    std::get<0>(config.logScores).wtScore(0, 0) = std::log(0.1);
    std::get<0>(config.logScores).hetScore(0, 0) = std::log(0.9);
    std::get<0>(config.logScores).homScore(0, 0) = std::log(0.5);

    std::get<0>(config.logScores).wtScore(1, 0) = std::log(0.2);
    std::get<0>(config.logScores).hetScore(1, 0) = std::log(0.8);
    std::get<0>(config.logScores).homScore(1, 0) = std::log(0.4);

    std::get<0>(config.logScores).wtScore(2, 0) = std::log(0.3);
    std::get<0>(config.logScores).hetScore(2, 0) = std::log(0.7);
    std::get<0>(config.logScores).homScore(2, 0) = std::log(0.3);

    std::get<0>(config.logScores).wtScore(3, 0) = std::log(0.4);
    std::get<0>(config.logScores).hetScore(3, 0) = std::log(0.6);
    std::get<0>(config.logScores).homScore(3, 0) = std::log(0.2);

    std::get<0>(config.logScores).wtScore(4, 0) = std::log(0.5);
    std::get<0>(config.logScores).hetScore(4, 0) = std::log(0.4);
    std::get<0>(config.logScores).homScore(4, 0) = std::log(0.1);

    std::get<0>(config.logScores).wtScore(5, 0) = std::log(0.6);
    std::get<0>(config.logScores).hetScore(5, 0) = std::log(0.3);
    std::get<0>(config.logScores).homScore(5, 0) = std::log(0.7);

    config.noiseScore = -1000;
    config.data.resize(6);
    config.data[0].resize(1);
    config._tmpAttachmentScore.resize(11);
    config.setParam(Config<SampleTree>::E_nu, 0.1);
    config.setParam(Config<SampleTree>::E_lambdaWildLoss, 0.2);
    config.setParam(Config<SampleTree>::E_lambdaMutLoss, 0.25);
    config.setParam(Config<SampleTree>::E_parallel, 0.3);
}
void initTree3_1(Config<SampleTree> &config) {
    /*
    *                                          11
    *                                          |
    *                                          0
    *                           ___________________________
    *                          |                           |
    *                          1                           2
    *               _______________________           ____________
    *              |                       |         |            |
    *              3                       8         9            10
    *       _______________                0.4       0.5          0.6
    *      |               |               0.6       0.4          0.3
    *      4               7               0.2       0.1          0.7
    *   _______            0.3
    *  |       |           0.7
    *  5       6           0.3
    *  0.1     0.2
    *  0.9     0.8
    *  0.5     0.4
    */
    config.setNumSamples(6);

    add_edge(0, 1, config.getTree());
    config.getTree()[0].sample = -1;
    config.getTree()[1].sample = -1;
    add_edge(0, 2, config.getTree());
    config.getTree()[2].sample = -1;
    add_edge(1, 3, config.getTree());
    config.getTree()[3].sample = -1;
    add_edge(1, 8, config.getTree());
    config.getTree()[8].sample = 3;
    add_edge(3, 4, config.getTree());
    config.getTree()[4].sample = -1;
    add_edge(3, 7, config.getTree());
    config.getTree()[7].sample = 2;

    add_edge(4, 5, config.getTree());
    config.getTree()[5].sample = 0;
    add_edge(4, 6, config.getTree());
    config.getTree()[6].sample = 1;
    add_edge(2, 9, config.getTree());
    config.getTree()[9].sample = 4;
    add_edge(2, 10, config.getTree());
    config.getTree()[10].sample = 5;

    add_edge(11, 0, config.getTree());

    std::get<0>(config.logScores).resizeNumCells(6);
    std::get<0>(config.logScores).resizeNumMuts(1);

    for (unsigned i = 0; i < 6; ++i)
    {
        std::get<0>(config.logScores).wtScore(i, 0) = 0;
        std::get<0>(config.logScores).hetScore(i, 0) = (i+1)/4.0 - 1.1;
        std::get<0>(config.logScores).homScore(i, 0) = (i+1)/3.0 - 1.4;
        std::cout << std::get<0>(config.logScores).hetScore(0, 0) << std::endl;
        std::cout << std::get<0>(config.logScores).homScore(0, 0) << std::endl;
    }

    config.noiseScore = -1000;
    config.data.resize(6);
    config.data[0].resize(1);
    config._tmpAttachmentScore.resize(11);
    config.setParam(Config<SampleTree>::E_nu, 0.15);
    config.setParam(Config<SampleTree>::E_lambdaWildLoss, 0.2);
    config.setParam(Config<SampleTree>::E_lambdaMutLoss, 0.25);
    config.setParam(Config<SampleTree>::E_parallel, 0.3);
}

void initTree4(Config<SampleTree> &config) {
    /*
    *                                          11
    *                                          |
    *                                          0
    *                           ___________________________
    *                          |                           |
    *                          1                           10
    *               _______________________                0.3
    *              |                       |               0.2
    *              2                       3               0.6
    *       _______________          ____________
    *      |               |         |           |
    *      4               7         8           9
    *   _______            0.3       0.4         0.5
    *  |       |           0.7       0.5         0.4
    *  5       6           0.3       0.2         0.1
    *  0.1     0.2
    *  0.9     0.8
    *  0.5     0.4
    */
    config.setNumSamples(8);

    add_edge(0, 1, config.getTree());
    config.getTree()[0].sample = -1;
    config.getTree()[1].sample = -1;
    add_edge(0, 10, config.getTree());
    config.getTree()[10].sample = 5;
    add_edge(1, 2, config.getTree());
    config.getTree()[2].sample = -1;
    add_edge(1, 3, config.getTree());
    config.getTree()[3].sample = -1;
    add_edge(2, 4, config.getTree());
    config.getTree()[4].sample = -1;
    add_edge(2, 7, config.getTree());
    config.getTree()[7].sample = 2;

    add_edge(3, 8, config.getTree());
    config.getTree()[8].sample = 3;
    add_edge(3, 9, config.getTree());
    config.getTree()[9].sample = 4;
    add_edge(4, 5, config.getTree());
    config.getTree()[5].sample = 0;
    add_edge(4, 6, config.getTree());
    config.getTree()[6].sample = 1;

    add_edge(11, 0, config.getTree());

    std::get<0>(config.logScores).resizeNumCells(6);
    std::get<0>(config.logScores).resizeNumMuts(1);
    std::get<0>(config.logScores).wtScore(0, 0) = std::log(0.1);
    std::get<0>(config.logScores).hetScore(0, 0) = std::log(0.9);
    std::get<0>(config.logScores).homScore(0, 0) = std::log(0.5);

    std::get<0>(config.logScores).wtScore(1, 0) = std::log(0.2);
    std::get<0>(config.logScores).hetScore(1, 0) = std::log(0.8);
    std::get<0>(config.logScores).homScore(1, 0) = std::log(0.4);

    std::get<0>(config.logScores).wtScore(2, 0) = std::log(0.3);
    std::get<0>(config.logScores).hetScore(2, 0) = std::log(0.7);
    std::get<0>(config.logScores).homScore(2, 0) = std::log(0.3);

    std::get<0>(config.logScores).wtScore(3, 0) = std::log(0.4);
    std::get<0>(config.logScores).hetScore(3, 0) = std::log(0.5);
    std::get<0>(config.logScores).homScore(3, 0) = std::log(0.2);

    std::get<0>(config.logScores).wtScore(4, 0) = std::log(0.5);
    std::get<0>(config.logScores).hetScore(4, 0) = std::log(0.4);
    std::get<0>(config.logScores).homScore(4, 0) = std::log(0.1);

    std::get<0>(config.logScores).wtScore(5, 0) = std::log(0.3);
    std::get<0>(config.logScores).hetScore(5, 0) = std::log(0.2);
    std::get<0>(config.logScores).homScore(5, 0) = std::log(0.6);

    config.noiseScore = -1000;
    config.data.resize(6);
    config.data[0].resize(1);
    config._tmpAttachmentScore.resize(11);
    config.setParam(Config<SampleTree>::E_nu, 0.1);
    config.setParam(Config<SampleTree>::E_lambdaWildLoss, 0.2);
    config.setParam(Config<SampleTree>::E_lambdaMutLoss, 0.25);
    config.setParam(Config<SampleTree>::E_parallel, 0.3);
}

void initTree5(Config<SampleTree> &config) {
    /*
    *                                                        13
    *                                                        |
    *                                                        0
    *                                           ______________________________
    *                                          |                              |
    *                                          1                              12
    *                           ___________________________                   0.7
    *                          |                           |                  0.2
    *                          2                           3                  0.8
    *               _______________________           ____________
    *              |                       |         |            |
    *              4                       9         10            11
    *       _______________                0.4       0.5          0.6
    *      |               |               0.6       0.4          0.3
    *      5               8               0.2       0.1          0.7
    *   _______            0.3
    *  |       |           0.7
    *  6       7           0.3
    *  0.1     0.2
    *  0.9     0.8
    *  0.5     0.4
    */
    config.setNumSamples(7);

    add_edge(0, 1, config.getTree());
    config.getTree()[0].sample = -1;
    config.getTree()[1].sample = -1;
    add_edge(0, 12, config.getTree());
    config.getTree()[12].sample = 6;
    add_edge(1, 2, config.getTree());
    config.getTree()[2].sample = -1;
    add_edge(1, 3, config.getTree());
    config.getTree()[3].sample = -1;
    add_edge(2, 4, config.getTree());
    config.getTree()[4].sample = -1;
    add_edge(2, 9, config.getTree());
    config.getTree()[9].sample = 3;
    add_edge(4, 5, config.getTree());
    config.getTree()[5].sample = -1;
    add_edge(4, 8, config.getTree());
    config.getTree()[8].sample = 2;
    add_edge(5, 6, config.getTree());
    config.getTree()[6].sample = 0;
    add_edge(5, 7, config.getTree());
    config.getTree()[7].sample = 1;
    add_edge(3, 10, config.getTree());
    config.getTree()[10].sample = 4;
    add_edge(3, 11, config.getTree());
    config.getTree()[11].sample = 5;

    add_edge(13, 0, config.getTree());

    std::get<0>(config.logScores).resizeNumCells(7);
    std::get<0>(config.logScores).resizeNumMuts(1);
    std::get<0>(config.logScores).wtScore(0, 0) = std::log(0.1);
    std::get<0>(config.logScores).hetScore(0, 0) = std::log(0.9);
    std::get<0>(config.logScores).homScore(0, 0) = std::log(0.5);

    std::get<0>(config.logScores).wtScore(1, 0) = std::log(0.2);
    std::get<0>(config.logScores).hetScore(1, 0) = std::log(0.8);
    std::get<0>(config.logScores).homScore(1, 0) = std::log(0.4);

    std::get<0>(config.logScores).wtScore(2, 0) = std::log(0.3);
    std::get<0>(config.logScores).hetScore(2, 0) = std::log(0.7);
    std::get<0>(config.logScores).homScore(2, 0) = std::log(0.3);

    std::get<0>(config.logScores).wtScore(3, 0) = std::log(0.4);
    std::get<0>(config.logScores).hetScore(3, 0) = std::log(0.6);
    std::get<0>(config.logScores).homScore(3, 0) = std::log(0.2);

    std::get<0>(config.logScores).wtScore(4, 0) = std::log(0.5);
    std::get<0>(config.logScores).hetScore(4, 0) = std::log(0.4);
    std::get<0>(config.logScores).homScore(4, 0) = std::log(0.1);

    std::get<0>(config.logScores).wtScore(5, 0) = std::log(0.6);
    std::get<0>(config.logScores).hetScore(5, 0) = std::log(0.3);
    std::get<0>(config.logScores).homScore(5, 0) = std::log(0.7);

    std::get<0>(config.logScores).wtScore(6, 0) = std::log(0.7);
    std::get<0>(config.logScores).hetScore(6, 0) = std::log(0.2);
    std::get<0>(config.logScores).homScore(6, 0) = std::log(0.8);

    config.noiseScore = -1000;
    config.data.resize(7);
    config.data[0].resize(1);
    config._tmpAttachmentScore.resize(13);
    config.setParam(Config<SampleTree>::E_nu, 0.1);
    config.setParam(Config<SampleTree>::E_lambdaWildLoss, 0.2);
    config.setParam(Config<SampleTree>::E_lambdaMutLoss, 0.25);
    config.setParam(Config<SampleTree>::E_parallel, 0.3);
}

void initTree6(Config<SampleTree> &config) {
    /*
    *                                                        15
    *                                                        |
    *                                                        0
    *                                           ______________________________
    *                                          |                              |
    *                                          1                              2
    *                           ___________________________               ____________
    *                          |                           |             |            |
    *                          3                           4             13           14
    *               _______________________           ____________       0.7          0.8
    *              |                       |         |            |      0.2          0.1
    *              5                       10        11            12    0.8          0.9
    *       _______________                0.4       0.5          0.6
    *      |               |               0.6       0.4          0.3
    *      6               9               0.2       0.1          0.7
    *   _______            0.3
    *  |       |           0.7
    *  7       8           0.3
    *  0.1     0.2
    *  0.9     0.8
    *  0.5     0.4
    */
    config.setNumSamples(8);

    add_edge(0, 1, config.getTree());
    config.getTree()[0].sample = -1;
    config.getTree()[1].sample = -1;
    add_edge(0, 2, config.getTree());
    config.getTree()[2].sample = -1;
    add_edge(1, 3, config.getTree());
    config.getTree()[3].sample = -1;
    add_edge(1, 4, config.getTree());
    config.getTree()[4].sample = -1;
    add_edge(3, 5, config.getTree());
    config.getTree()[5].sample = -1;
    add_edge(3, 10, config.getTree());
    config.getTree()[10].sample = 3;
    add_edge(5, 6, config.getTree());
    config.getTree()[6].sample = -1;
    add_edge(5, 9, config.getTree());
    config.getTree()[9].sample = 2;
    add_edge(6, 7, config.getTree());
    config.getTree()[7].sample = 0;
    add_edge(6, 8, config.getTree());
    config.getTree()[8].sample = 1;
    add_edge(4, 11, config.getTree());
    config.getTree()[11].sample = 4;
    add_edge(4, 12, config.getTree());
    config.getTree()[12].sample = 5;
    add_edge(2, 13, config.getTree());
    config.getTree()[13].sample = 6;
    add_edge(2, 14, config.getTree());
    config.getTree()[14].sample = 7;

    add_edge(15, 0, config.getTree());

    std::get<0>(config.logScores).resizeNumCells(8);
    std::get<0>(config.logScores).resizeNumMuts(1);
    std::get<0>(config.logScores).wtScore(0, 0) = std::log(0.1);
    std::get<0>(config.logScores).hetScore(0, 0) = std::log(0.9);
    std::get<0>(config.logScores).homScore(0, 0) = std::log(0.5);

    std::get<0>(config.logScores).wtScore(1, 0) = std::log(0.2);
    std::get<0>(config.logScores).hetScore(1, 0) = std::log(0.8);
    std::get<0>(config.logScores).homScore(1, 0) = std::log(0.4);

    std::get<0>(config.logScores).wtScore(2, 0) = std::log(0.3);
    std::get<0>(config.logScores).hetScore(2, 0) = std::log(0.7);
    std::get<0>(config.logScores).homScore(2, 0) = std::log(0.3);

    std::get<0>(config.logScores).wtScore(3, 0) = std::log(0.4);
    std::get<0>(config.logScores).hetScore(3, 0) = std::log(0.6);
    std::get<0>(config.logScores).homScore(3, 0) = std::log(0.2);

    std::get<0>(config.logScores).wtScore(4, 0) = std::log(0.5);
    std::get<0>(config.logScores).hetScore(4, 0) = std::log(0.4);
    std::get<0>(config.logScores).homScore(4, 0) = std::log(0.1);

    std::get<0>(config.logScores).wtScore(5, 0) = std::log(0.6);
    std::get<0>(config.logScores).hetScore(5, 0) = std::log(0.3);
    std::get<0>(config.logScores).homScore(5, 0) = std::log(0.7);

    std::get<0>(config.logScores).wtScore(6, 0) = std::log(0.7);
    std::get<0>(config.logScores).hetScore(6, 0) = std::log(0.2);
    std::get<0>(config.logScores).homScore(6, 0) = std::log(0.8);

    std::get<0>(config.logScores).wtScore(7, 0) = std::log(0.8);
    std::get<0>(config.logScores).hetScore(7, 0) = std::log(0.1);
    std::get<0>(config.logScores).homScore(7, 0) = std::log(0.9);

    config.noiseScore = -1000;
    config.data.resize(8);
    config.data[0].resize(1);
    config._tmpAttachmentScore.resize(15);
    config.setParam(Config<SampleTree>::E_nu, 0.1);
    config.setParam(Config<SampleTree>::E_lambdaWildLoss, 0.2);
    config.setParam(Config<SampleTree>::E_lambdaMutLoss, 0.25);
    config.setParam(Config<SampleTree>::E_parallel, 0.3);
}

BOOST_AUTO_TEST_CASE(adjuste_coefficient) {
    BOOST_CHECK(adjusteCoefficient(0.6) == 0.6);
    BOOST_CHECK(adjusteCoefficient(1.6) - 0.4 <= 10.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(adjusteCoefficient(-0.6) == 0.6);
    BOOST_CHECK(adjusteCoefficient(-1.6) - 0.4 <= 10.0 * std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(sample_tree_init) {
    Config<SampleTree> config;
    config.setNumSamples(3);

    createInitialTreeTest(config);

    BOOST_CHECK_MESSAGE(num_vertices(config.getTree()) == 6, num_vertices(config.getTree()) << " != 5");

    BOOST_CHECK(source(*in_edges(1, config.getTree()).first, config.getTree()) == 0);
    BOOST_CHECK(source(*in_edges(2, config.getTree()).first, config.getTree()) == 0);
    BOOST_CHECK(source(*in_edges(3, config.getTree()).first, config.getTree()) == 1);
    BOOST_CHECK(source(*in_edges(4, config.getTree()).first, config.getTree()) == 1);

    BOOST_CHECK(config.getTree()[0].sample == -1);
    BOOST_CHECK(config.getTree()[1].sample == -1);
    BOOST_CHECK(config.getTree()[2].sample == 0);
    BOOST_CHECK(config.getTree()[3].sample == 1);
    BOOST_CHECK(config.getTree()[4].sample == 2);

    // artificial root
    BOOST_CHECK_MESSAGE(source(*in_edges(0, config.getTree()).first, config.getTree()) == 5,
                        source(*in_edges(0, config.getTree()).first, config.getTree()) << " != 5");
}

BOOST_AUTO_TEST_CASE(get_sibling) {

    for (unsigned i = 0; i < 10000; ++i) {
        Config<SampleTree> config;
        config.setNumSamples(7);

        createInitialTreeTest(config);

        unsigned newNonDescendant = getNewSibling(config, 1);

        BOOST_CHECK_MESSAGE(newNonDescendant != 0, newNonDescendant << " == 0");
        BOOST_CHECK_MESSAGE(newNonDescendant != 1, newNonDescendant << " == 1");
        BOOST_CHECK_MESSAGE(newNonDescendant != 3, newNonDescendant << " == 3");
        BOOST_CHECK_MESSAGE(newNonDescendant != 4, newNonDescendant << " == 4");
        BOOST_CHECK_MESSAGE(newNonDescendant != 7, newNonDescendant << " == 7");
        BOOST_CHECK_MESSAGE(newNonDescendant != 8, newNonDescendant << " == 8");
        BOOST_CHECK_MESSAGE(newNonDescendant != 9, newNonDescendant << " == 9");
        BOOST_CHECK_MESSAGE(newNonDescendant != 10, newNonDescendant << " == 10");
        BOOST_CHECK_MESSAGE(newNonDescendant != 13, newNonDescendant << " == 13");
    }
}

BOOST_AUTO_TEST_CASE(sample_tree_prune_and_reattach) {

    std::vector<unsigned> counter = {0, 0, 0, 0, 0};
    for (unsigned i = 0; i < 10000; ++i) {
        bool hit = false;
        Config<SampleTree> config;
        config.setNumSamples(3);
        createInitialTreeTest(config);

        pruneAndReAttach(config);

        if (source(*in_edges(0, config.getTree()).first, config.getTree()) == 1 &&
            source(*in_edges(1, config.getTree()).first, config.getTree()) == 5 &&
            source(*in_edges(2, config.getTree()).first, config.getTree()) == 0 &&
            source(*in_edges(3, config.getTree()).first, config.getTree()) == 0 &&
            source(*in_edges(4, config.getTree()).first, config.getTree()) == 1) {
            hit = true;
            ++counter[0];
        } else if (source(*in_edges(0, config.getTree()).first, config.getTree()) == 5 &&
                   source(*in_edges(1, config.getTree()).first, config.getTree()) == 0 &&
                   source(*in_edges(2, config.getTree()).first, config.getTree()) == 1 &&
                   source(*in_edges(3, config.getTree()).first, config.getTree()) == 0 &&
                   source(*in_edges(4, config.getTree()).first, config.getTree()) == 1) {
            hit = true;
            ++counter[1];
        } else if (source(*in_edges(0, config.getTree()).first, config.getTree()) == 5 &&
                   source(*in_edges(1, config.getTree()).first, config.getTree()) == 0 &&
                   source(*in_edges(2, config.getTree()).first, config.getTree()) == 0 &&
                   source(*in_edges(3, config.getTree()).first, config.getTree()) == 1 &&
                   source(*in_edges(4, config.getTree()).first, config.getTree()) == 1) {
            hit = true;
            ++counter[2];
        } else if (source(*in_edges(0, config.getTree()).first, config.getTree()) == 5 &&
                   source(*in_edges(1, config.getTree()).first, config.getTree()) == 0 &&
                   source(*in_edges(2, config.getTree()).first, config.getTree()) == 1 &&
                   source(*in_edges(3, config.getTree()).first, config.getTree()) == 1 &&
                   source(*in_edges(4, config.getTree()).first, config.getTree()) == 0) {
            hit = true;
            ++counter[3];
        } else if (source(*in_edges(0, config.getTree()).first, config.getTree()) == 1 &&
                   source(*in_edges(1, config.getTree()).first, config.getTree()) == 5 &&
                   source(*in_edges(2, config.getTree()).first, config.getTree()) == 0 &&
                   source(*in_edges(3, config.getTree()).first, config.getTree()) == 1 &&
                   source(*in_edges(4, config.getTree()).first, config.getTree()) == 0) {
            hit = true;
            ++counter[4];
        }

        BOOST_CHECK_MESSAGE(hit, "Invalid tree configuration.");
        if (!hit)
            break;
    }
    BOOST_CHECK_MESSAGE(counter[0] > 0, "1 configuration never reached");
    BOOST_CHECK_MESSAGE(counter[1] > 0, "2 configuration never reached");
    BOOST_CHECK_MESSAGE(counter[2] > 0, "3 configuration never reached");
    BOOST_CHECK_MESSAGE(counter[3] > 0, "4 configuration never reached");
    BOOST_CHECK_MESSAGE(counter[4] > 0, "5 configuration never reached");

}

BOOST_AUTO_TEST_CASE(sample_tree_swap_node_labels) {
    std::vector<unsigned> counter = {0, 0, 0};
    for (unsigned i = 0; i < 10000; ++i) {
        bool hit = false;
        Config<SampleTree> config;
        config.setNumSamples(3);
        createInitialTreeTest(config);

        BOOST_CHECK(source(*in_edges(0, config.getTree()).first, config.getTree()) == 5);
        BOOST_CHECK(source(*in_edges(1, config.getTree()).first, config.getTree()) == 0);
        BOOST_CHECK(source(*in_edges(2, config.getTree()).first, config.getTree()) == 0);
        BOOST_CHECK(source(*in_edges(3, config.getTree()).first, config.getTree()) == 1);
        BOOST_CHECK(source(*in_edges(4, config.getTree()).first, config.getTree()) == 1);

        swapNodeLabels(config);

        if (config.getTree()[0].sample == -1 &&
            config.getTree()[1].sample == -1 &&
            config.getTree()[2].sample == 1 &&
            config.getTree()[3].sample == 0 &&
            config.getTree()[4].sample == 2) {
            hit = true;
            ++counter[0];
        } else if (config.getTree()[0].sample == -1 &&
                   config.getTree()[1].sample == -1 &&
                   config.getTree()[2].sample == 0 &&
                   config.getTree()[3].sample == 2 &&
                   config.getTree()[4].sample == 1) {
            hit = true;
            ++counter[1];
        } else if (config.getTree()[0].sample == -1 &&
                   config.getTree()[1].sample == -1 &&
                   config.getTree()[2].sample == 2 &&
                   config.getTree()[3].sample == 1 &&
                   config.getTree()[4].sample == 0) {
            hit = true;
            ++counter[2];
        }

        BOOST_CHECK_MESSAGE(hit, "Invalid tree configuration.");
        if (!hit) {
            for (unsigned i = 0; i < num_vertices(config.getTree()); ++i)
                std::cout << "Sample of node " << i << " = " << config.getTree()[i].sample << std::endl;
            break;
        }
    }
    for (unsigned i = 0; i < counter.size(); ++i)
        BOOST_CHECK_MESSAGE(counter[i] > 0, i << ". configuration never reached");

}


BOOST_AUTO_TEST_CASE(sample_tree_get_mutation_of_nodes) {
    Config<SampleTree> config;
    config.setNumSamples(3);
    createInitialTreeTest(config);

    /*
            0
            |
        1_________2
        |
     3_______4

    */

    BOOST_CHECK(source(*in_edges(0, config.getTree()).first, config.getTree()) == 5);
    BOOST_CHECK(source(*in_edges(1, config.getTree()).first, config.getTree()) == 0);
    BOOST_CHECK(source(*in_edges(2, config.getTree()).first, config.getTree()) == 0);
    BOOST_CHECK(source(*in_edges(3, config.getTree()).first, config.getTree()) == 1);
    BOOST_CHECK(source(*in_edges(4, config.getTree()).first, config.getTree()) == 1);

    config.data.resize(3);
    config.data[0].resize(3, std::tuple<unsigned, unsigned>{10, 0});
    config.data[1].resize(3, std::tuple<unsigned, unsigned>{10, 0});
    config.data[2].resize(3, std::tuple<unsigned, unsigned>{10, 0});

    // cell 1 is node 2
    // cell 2 is node 3
    // cekk 3 is node 4
    // cell 2 and 3 have mutation 1
    std::get<1>(config.data[1][0]) = 10;
    std::get<1>(config.data[2][0]) = 10;

    // cell 2 has mutation 2
    std::get<1>(config.data[1][1]) = 10;

    // cell 2 has mutation 3
    std::get<1>(config.data[1][2]) = 10;

    config.getLogScores().resizeNumCells(config.getNumSamples());
    config.getTmpLogScores().resizeNumCells(config.getNumSamples());
    config.getLogScores().resizeNumMuts(config.getNumMutations());
    computeLogScoresOP(config);
    config._tmpAttachmentScore.resize(5);

    std::vector<std::vector<unsigned> > mutationsOfNodes = getMutationsOfNode(config);

    BOOST_CHECK_MESSAGE(mutationsOfNodes.size() == 6,
                        "mutationsOfNodes.size() = " << mutationsOfNodes.size() << " != 6");
    BOOST_CHECK_MESSAGE(mutationsOfNodes[0].size() == 0,
                        "mutationsOfNodes[0].size() = " << mutationsOfNodes[0].size() << " != 0");
    BOOST_CHECK_MESSAGE(mutationsOfNodes[1][0] == 0, "mutationsOfNodes[1][0] = " << mutationsOfNodes[1][0] << " != 0");
    BOOST_CHECK_MESSAGE(mutationsOfNodes[2].size() == 0,
                        "mutationsOfNodes[2].size() = " << mutationsOfNodes[2].size() << " != 0");
    BOOST_CHECK_MESSAGE(mutationsOfNodes[3][0] == 1, "mutationsOfNodes[3][0] = " << mutationsOfNodes[3][0] << " != 1");
    BOOST_CHECK_MESSAGE(mutationsOfNodes[3][1] == 2, "mutationsOfNodes[3][1] = " << mutationsOfNodes[3][1] << " != 2");
    BOOST_CHECK_MESSAGE(mutationsOfNodes[4].size() == 0,
                        "mutationsOfNodes[4].size() = " << mutationsOfNodes[4].size() << " != 0");
    BOOST_CHECK_MESSAGE(mutationsOfNodes[5].size() == 0,
                        "mutationsOfNodes[5].size() = " << mutationsOfNodes[5].size() << " != 0");
}


BOOST_AUTO_TEST_CASE(sample_tree_save_load_tree) {
    Config<SampleTree> config;
    initTree1(config);
    config.learnZygocity = true;
    config.computeLossScore = true;
    config.computeParallelScore = true;
    config.noiseScore = -1000;
    
    config.loadName = "test_tree";
    std::cout << "SCORE: " << scoreTree(config) << std::endl;

    std::cerr << "0" << std::endl;
    writeIndex(config, config.loadName);
    std::cerr << "1" << std::endl;

    readCellNames(config);
    std::cerr << "2" << std::endl;
    std::string treePath = config.loadName + "/tree.gv";
    std::cerr << "3" << std::endl;
    std::ifstream f(treePath.c_str());
    std::cerr << "4" << std::endl;
    std::cout<< config.bestName << std::endl;
    std::cerr << "5" << std::endl;
    readGraph(config);
    std::cerr << "6" << std::endl;
    readNucInfo(config);
    std::cerr << "7" << std::endl;

    std::cout << "SCORE: " << scoreTree(config) << std::endl;

    }

BOOST_AUTO_TEST_CASE(sample_tree_save_load_tree2) {
    Config<SampleTree> config;
    initTree2(config);
    config.learnZygocity = true;
    config.computeLossScore = true;
    config.computeParallelScore = true;
    config.noiseScore = -1000;
    
    config.loadName = "test_tree";
    std::cout << "SCORE: " << scoreTree(config) << std::endl;

    writeIndex(config, config.loadName);

    readCellNames(config);
    std::string treePath = config.loadName + "/tree.gv";
    std::ifstream f(treePath.c_str());
    std::cout<< config.bestName << std::endl;
    readGraph(config);
    readNucInfo(config);

    std::cout << "SCORE: " << scoreTree(config) << std::endl;

    }
BOOST_AUTO_TEST_CASE(sample_tree_save_load_tree5) {
    {
        Config<SampleTree> config;
        initTree3(config);
        config.learnZygocity = true;
        config.computeLossScore = true;
        config.computeParallelScore = true;
        config.noiseScore = -1000;
        
        config.loadName = "test_tree";
        std::cout << "SCORE3: " << scoreTree(config) << std::endl;

        writeIndex(config, config.loadName);

        readCellNames(config);
        std::string treePath = config.loadName + "/tree.gv";
        std::ifstream f(treePath.c_str());
        std::cout<< config.bestName << std::endl;
        readGraph(config);
        readNucInfo(config);

        std::cout << "SCORE3: " << scoreTree(config) << std::endl;
    }
    {
        Config<SampleTree> config;
        initTree4(config);
        config.learnZygocity = true;
        config.computeLossScore = true;
        config.computeParallelScore = true;
        config.noiseScore = -1000;
        
        config.loadName = "test_tree";
        std::cout << "SCORE4: " << scoreTree(config) << std::endl;

        writeIndex(config, config.loadName);

        readCellNames(config);
        std::string treePath = config.loadName + "/tree.gv";
        std::ifstream f(treePath.c_str());
        std::cout<< config.bestName << std::endl;
        readGraph(config);
        readNucInfo(config);

        std::cout << "SCORE4: " << scoreTree(config) << std::endl;
    }
    {
        Config<SampleTree> config;
        initTree5(config);
        config.learnZygocity = true;
        config.computeLossScore = true;
        config.computeParallelScore = true;
        config.noiseScore = -1000;
        
        config.loadName = "test_tree";
        std::cout << "SCORE5: " << scoreTree(config) << std::endl;

        writeIndex(config, config.loadName);

        readCellNames(config);
        std::string treePath = config.loadName + "/tree.gv";
        std::ifstream f(treePath.c_str());
        std::cout<< config.bestName << std::endl;
        readGraph(config);
        readNucInfo(config);

        std::cout << "SCORE5: " << scoreTree(config) << std::endl;
    }
    {
        Config<SampleTree> config;
        initTree6(config);
        config.learnZygocity = true;
        config.computeLossScore = true;
        config.computeParallelScore = true;
        config.noiseScore = -1000;
        
        config.loadName = "test_tree";
        std::cout << "SCORE6: " << scoreTree(config) << std::endl;

        writeIndex(config, config.loadName);

        readCellNames(config);
        std::string treePath = config.loadName + "/tree.gv";
        std::ifstream f(treePath.c_str());
        std::cout<< config.bestName << std::endl;
        readGraph(config);
        readNucInfo(config);

        std::cout << "SCORE6: " << scoreTree(config) << std::endl;
    }

    }

BOOST_AUTO_TEST_CASE(sample_tree_get_num_placements) {
    {
        Config<SampleTree> config;
        initTree1(config);
        config.computeLossScore = true;
        config.computeParallelScore = true;
        scoreTree(config);

        BOOST_CHECK(config._tmpAttachmentScore.numAltRPoss(4) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numAltRPoss(3) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numAltRPoss(2) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numAltRPoss(1) == 1);
        BOOST_CHECK(config._tmpAttachmentScore.numAltRPoss(0) == 2);

        BOOST_CHECK(config._tmpAttachmentScore.numInnerNodes(3) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numInnerNodes(3) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numInnerNodes(2) == 1);
        BOOST_CHECK(config._tmpAttachmentScore.numInnerNodes(1) == 2);
        BOOST_CHECK(config._tmpAttachmentScore.numInnerNodes(0) == 3);

        BOOST_CHECK(config._tmpAttachmentScore.numInnerChildNodes(4) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numInnerChildNodes(3) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numInnerChildNodes(2) == 1);
        BOOST_CHECK(config._tmpAttachmentScore.numInnerChildNodes(1) == 1);
        BOOST_CHECK(config._tmpAttachmentScore.numInnerChildNodes(0) == 1);

        BOOST_CHECK(config._tmpAttachmentScore.numLcaRPoss(4) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numLcaRPoss(3) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numLcaRPoss(2) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numLcaRPoss(1) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numLcaRPoss(0) == 0);
    }
    {
        Config<SampleTree> config;
        initTree2(config);
        config.computeLossScore = true;
        config.computeParallelScore = true;
        scoreTree(config);

        BOOST_CHECK(config._tmpAttachmentScore.numAltRPoss(7) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numAltRPoss(5) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numAltRPoss(3) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numAltRPoss(1) == 2);
        BOOST_CHECK(config._tmpAttachmentScore.numAltRPoss(0) == 4);

        BOOST_CHECK(config._tmpAttachmentScore.numInnerNodes(7) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numInnerNodes(5) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numInnerNodes(3) == 2);
        BOOST_CHECK(config._tmpAttachmentScore.numInnerNodes(1) == 4);
        BOOST_CHECK(config._tmpAttachmentScore.numInnerNodes(0) == 6);

        BOOST_CHECK(config._tmpAttachmentScore.numInnerChildNodes(7) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numInnerChildNodes(5) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numInnerChildNodes(3) == 2);
        BOOST_CHECK(config._tmpAttachmentScore.numInnerChildNodes(1) == 2);
        BOOST_CHECK(config._tmpAttachmentScore.numInnerChildNodes(0) == 2);

        BOOST_CHECK(config._tmpAttachmentScore.numLcaRPoss(7) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numLcaRPoss(5) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numLcaRPoss(3) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numLcaRPoss(1) == 0);
        BOOST_CHECK(config._tmpAttachmentScore.numLcaRPoss(0) == 2);

    }
}

BOOST_AUTO_TEST_CASE(sample_tree_get_all_attachment_scores) {
    {
        /*
        *
        *
        *              0
        *       _______________
        *      |               |
        *      1               2
        *   _______            0.3
        *  |       |           0.7
        *  3       4           0.3
        *  0.1     0.2
        *  0.9     0.8
        *  0.5     0.4
        */

        Config<SampleTree> config;
        config.setNumSamples(3);

        add_edge(0, 1, config.getTree());
        config.getTree()[0].sample = -1;
        config.getTree()[1].sample = -1;
        add_edge(0, 2, config.getTree());
        config.getTree()[2].sample = 2;
        add_edge(1, 3, config.getTree());
        config.getTree()[3].sample = 0;
        add_edge(1, 4, config.getTree());
        config.getTree()[4].sample = 1;

        add_edge(5, 0, config.getTree());

        std::get<0>(config.logScores).resizeNumCells(3);
        std::get<0>(config.logScores).resizeNumMuts(1);
        std::get<0>(config.logScores).wtScore(0, 0) = std::log(0.1);
        std::get<0>(config.logScores).hetScore(0, 0) = std::log(0.9);
        std::get<0>(config.logScores).homScore(0, 0) = std::log(0.5);

        std::get<0>(config.logScores).wtScore(1, 0) = std::log(0.2);
        std::get<0>(config.logScores).hetScore(1, 0) = std::log(0.8);
        std::get<0>(config.logScores).homScore(1, 0) = std::log(0.4);

        std::get<0>(config.logScores).wtScore(2, 0) = std::log(0.3);
        std::get<0>(config.logScores).hetScore(2, 0) = std::log(0.7);
        std::get<0>(config.logScores).homScore(2, 0) = std::log(0.3);

        config.noiseScore = -1000;
        config.data.resize(3);
        config.data[0].resize(1);
        config._tmpAttachmentScore.resize(5);

        config.computeLossScore = false;
        config.setParam(Config<SampleTree>::E_nu, 0);
        config.setParam(Config<SampleTree>::E_lambdaWildLoss, 0);
        config.setParam(Config<SampleTree>::E_lambdaMutLoss, 0);
        config.setParam(Config<SampleTree>::E_parallel, 0);

        scoreTree(config);

        double result = std::log(0.9 * 0.2 * 0.3);
        BOOST_CHECK_MESSAGE(
                std::abs(std::log(std::exp(config._tmpAttachmentScore.hetScore(3)) * 0.1 * 0.2 * 0.3) - result)
                <= std::numeric_limits<double>::epsilon(),
                std::log(std::exp(config._tmpAttachmentScore.hetScore(3)) * 0.1 * 0.2 * 0.3)
                        << " != " << result);

        result = std::log(0.9 * 0.8 * 0.7);
        BOOST_CHECK_MESSAGE(
                std::abs(std::log(std::exp(config._tmpAttachmentScore.hetScore(0)) * 0.1 * 0.2 * 0.3) - result)
                <= std::numeric_limits<double>::epsilon(),
                std::log(std::exp(config._tmpAttachmentScore.hetScore(0)) * 0.1 * 0.2 * 0.3)
                        << " != " << result);


        result = std::log(0.9 * 0.2 * 0.3 + 0.1 * 0.8 * 0.3 + 0.1 * 0.2 * 0.7 + 0.9 * 0.8 * 0.3 + 0.9 * 0.8 * 0.7);
        Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;
        getAllAttachmentScores(config._tmpAttachmentScore, scoreSum, config, 0);
        BOOST_CHECK_MESSAGE(std::abs(std::log(std::exp(scoreSum.hetScore()) * 0.1 * 0.2 * 0.3) - result)
                            <= 10 * std::numeric_limits<double>::epsilon(),
                            std::log(std::exp(scoreSum.hetScore()) * 0.1 * 0.2 * 0.3)
                                    << " != " << result);


        result = result - std::log(5) + config.noiseScore;
        BOOST_CHECK_MESSAGE(std::abs(sumScoreTree(config) - result)
                            <= std::numeric_limits<double>::epsilon(),
                            sumScoreTree(config) << " != " << result);
    }

}

BOOST_AUTO_TEST_CASE(sample_tree_compute_score_het_3) {
    Config<SampleTree> config;
    initTree3(config);
    scoreTree(config);

    double result = std::log(0.9) - std::log(0.1);
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.hetScore(5) - result)
                <= 10 * std::numeric_limits<double>::epsilon());

    result += std::log(0.8) - std::log(0.2);
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.hetScore(4) - result)
                <= 10 * std::numeric_limits<double>::epsilon());

    result += std::log(0.7) - std::log(0.3);
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.hetScore(3) - result)
                <= 10 * std::numeric_limits<double>::epsilon());

    result += std::log(0.6) - std::log(0.4);
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.hetScore(1) - result)
                <= 10 * std::numeric_limits<double>::epsilon());

    result += (std::log(0.4) - std::log(0.5)
               + std::log(0.3) - std::log(0.6));
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.hetScore(0) - result)
                <= 100 * std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(sample_tree_compute_score_hom_3) {
    Config<SampleTree> config;
    initTree3(config);
    scoreTree(config);

    double result = std::log(0.5) - std::log(0.1);
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.homScore(5) - result)
                <= 10 * std::numeric_limits<double>::epsilon());

    result += std::log(0.4) - std::log(0.2);
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.homScore(4) - result)
                <= 10 * std::numeric_limits<double>::epsilon());

    result += std::log(0.3) - std::log(0.3);
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.homScore(3) - result)
                <= 10 * std::numeric_limits<double>::epsilon());

    result += std::log(0.2) - std::log(0.4);
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.homScore(1) - result)
                <= 10 * std::numeric_limits<double>::epsilon());

    result += (std::log(0.1) - std::log(0.5)
               + std::log(0.7) - std::log(0.6));
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.homScore(0) - result)
                <= 100 * std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(sample_tree_compute_score_lossR_3) {
    Config<SampleTree> config;
    initTree3(config);
    config.computeLossScore = true;
    scoreTree(config);

    BOOST_CHECK(config._tmpAttachmentScore.lossWildScore(5) == -INFINITY);
    BOOST_CHECK(config._tmpAttachmentScore.lossWildScore(4) == -INFINITY);

    double result = std::log(0.5) - std::log(0.1) +
                    std::log(0.4) - std::log(0.2) +
                    std::log(0.7) - std::log(0.3);
    BOOST_CHECK(config._tmpAttachmentScore.lossWildScore(3) == result);

    result = addLogProb(std::log(0.5) - std::log(0.1) +
                        std::log(0.4) - std::log(0.2) +
                        std::log(0.7) - std::log(0.3) +
                        std::log(0.6) - std::log(0.4),
                        std::log(0.5) - std::log(0.1) +
                        std::log(0.4) - std::log(0.2) +
                        std::log(0.3) - std::log(0.3) +
                        std::log(0.6) - std::log(0.4));
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.lossWildScore(1) - result)
                <= 10 * std::numeric_limits<double>::epsilon());

    result = addLogProb(addLogProb(addLogProb(
            std::log(0.5) - std::log(0.1) +
            std::log(0.4) - std::log(0.2) +
            std::log(0.7) - std::log(0.3) +
            std::log(0.6) - std::log(0.4) +
            std::log(0.4) - std::log(0.5) +
            std::log(0.3) - std::log(0.6),
            std::log(0.5) - std::log(0.1) +
            std::log(0.4) - std::log(0.2) +
            std::log(0.3) - std::log(0.3) +
            std::log(0.6) - std::log(0.4) +
            std::log(0.4) - std::log(0.5) +
            std::log(0.3) - std::log(0.6)),
                                   std::log(0.5) - std::log(0.1) +
                                   std::log(0.4) - std::log(0.2) +
                                   std::log(0.3) - std::log(0.3) +
                                   std::log(0.2) - std::log(0.4) +
                                   std::log(0.4) - std::log(0.5) +
                                   std::log(0.3) - std::log(0.6)),
                        std::log(0.9) - std::log(0.1) +
                        std::log(0.8) - std::log(0.2) +
                        std::log(0.7) - std::log(0.3) +
                        std::log(0.6) - std::log(0.4) +
                        std::log(0.1) - std::log(0.5) +
                        std::log(0.7) - std::log(0.6));
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.lossWildScore(0) - result)
                <= 100 * std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(sample_tree_compute_score_lossA_3) {
    Config<SampleTree> config;
    initTree3(config);
    config.computeLossScore = true;
    scoreTree(config);

    BOOST_CHECK(config._tmpAttachmentScore.lossAltRScore(5) == -INFINITY);
    BOOST_CHECK(config._tmpAttachmentScore.lossAltRScore(4) == -INFINITY);
    BOOST_CHECK(config._tmpAttachmentScore.lossAltRScore(3) == -INFINITY);

    double result = std::log(0.7) - std::log(0.3) +
                    std::log(0.6) - std::log(0.4);
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.lossAltRScore(1) - result)
                <= 10 * std::numeric_limits<double>::epsilon());

    result = addLogProb(
            std::log(0.7) - std::log(0.3) +
            std::log(0.6) - std::log(0.4) +
            std::log(0.4) - std::log(0.5) +
            std::log(0.3) - std::log(0.6),
            std::log(0.6) - std::log(0.4) +
            std::log(0.4) - std::log(0.5) +
            std::log(0.3) - std::log(0.6));
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.lossAltRScore(0) - result)
                <= 100 * std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(sample_tree_compute_score_parallel_3) {
    Config<SampleTree> config;
    initTree3(config);
    config.computeLossScore = true;
    config.computeParallelScore = true;
    scoreTree(config);

    BOOST_CHECK(config._tmpAttachmentScore.lcaRScore(5) == -INFINITY);
    BOOST_CHECK(config._tmpAttachmentScore.lcaRScore(4) == -INFINITY);
    BOOST_CHECK(config._tmpAttachmentScore.lcaRScore(3) == -INFINITY);
    BOOST_CHECK(config._tmpAttachmentScore.lcaRScore(1) == -INFINITY);

    double result =
            std::log(0.9) - std::log(0.1) +
            std::log(0.8) - std::log(0.2) +
            std::log(0.4) - std::log(0.5) +
            std::log(0.3) - std::log(0.6);
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.lcaRScore(0) - result)
                <= 100 * std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(sample_tree_compute_score_parallel_4) {
    Config<SampleTree> config;
    initTree4(config);
    config.computeLossScore = true;
    config.computeParallelScore = true;
    scoreTree(config);

    for (unsigned i = 0; i <= 10; ++i) {
        BOOST_CHECK(config._tmpAttachmentScore.lcaRScore(i) == -INFINITY);
    }
}

BOOST_AUTO_TEST_CASE(sample_tree_compute_score_parallel_5) {
    Config<SampleTree> config;
    initTree5(config);
    config.computeLossScore = true;
    config.computeParallelScore = true;
    scoreTree(config);

    for (unsigned i = 0; i <= 10; ++i) {
        if (i != 1) {
            BOOST_CHECK(config._tmpAttachmentScore.lcaRScore(i) == -INFINITY);
        }
    }

    double result = std::log(0.9) - std::log(0.1) +
                    std::log(0.8) - std::log(0.2) +
                    std::log(0.3) - std::log(0.3) +
                    std::log(0.4) - std::log(0.4) +
                    std::log(0.4) - std::log(0.5) +
                    std::log(0.3) - std::log(0.6);
    std::cout << config._tmpAttachmentScore.lcaRScore(1) << " " << result << std::endl;

    BOOST_CHECK(std::abs(config._tmpAttachmentScore.lcaRScore(1) - result)
                <= 100 * std::numeric_limits<double>::epsilon());

    Config<SampleTree>::TAttachmentScores &attachmentScores = config.getTmpAttachmentScore();
    static Config<SampleTree>::TAttachmentScores attachmentSumScores;
    attachmentSumScores.resize(attachmentScores.size());
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentScores;
    passDownAttachmentScores.resize(attachmentScores.size());
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentSumScores;
    passDownAttachmentSumScores.resize(attachmentScores.size());

    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;
    long double sumParallel = 0;
    PassScoreToChildrenBFSVisitor visBFS(config,
                                         attachmentScores,
                                         attachmentSumScores,
                                         passDownAttachmentScores,
                                         passDownAttachmentSumScores,
                                         sumParallel,
                                         0);
    breadth_first_search(config.getTree(), num_vertices(config.getTree()) - 1, visitor(visBFS));

    BOOST_CHECK(std::abs(attachmentScores[1].lcaRScore() -
                         passDownAttachmentSumScores[5].paralleleScore())
                < 100.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(std::abs(passDownAttachmentSumScores[5].paralleleScore() -
                         passDownAttachmentSumScores[6].paralleleScore())
                < 100.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(std::abs(passDownAttachmentSumScores[5].paralleleScore() -
                         passDownAttachmentSumScores[10].paralleleScore())
                < 100.0 * std::numeric_limits<double>::epsilon());

    BOOST_CHECK(passDownAttachmentScores[2].paralleleScore() == -INFINITY);
    BOOST_CHECK(passDownAttachmentSumScores[2].paralleleScore() == -INFINITY);
    BOOST_CHECK(passDownAttachmentSumScores[8].paralleleScore() == -INFINITY);
    BOOST_CHECK(passDownAttachmentScores[12].paralleleScore() == -INFINITY);
    BOOST_CHECK(passDownAttachmentSumScores[12].paralleleScore() == -INFINITY);
}

BOOST_AUTO_TEST_CASE(sample_tree_compute_score_parallel_6) {
    Config<SampleTree> config;
    initTree6(config);
    config.computeLossScore = true;
    config.computeParallelScore = true;
    scoreTree(config);

    for (unsigned i = 0; i <= 10; ++i) {
        if (i != 1 && i != 0) {
            BOOST_CHECK(config._tmpAttachmentScore.lcaRScore(i) == -INFINITY);
        }
    }

    double result = std::log(0.9) - std::log(0.1) +
                    std::log(0.8) - std::log(0.2) +
                    std::log(0.3) - std::log(0.3) +
                    std::log(0.4) - std::log(0.4) +
                    std::log(0.4) - std::log(0.5) +
                    std::log(0.3) - std::log(0.6);
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.lcaRScore(1) - result)
                <= 1000 * std::numeric_limits<double>::epsilon());

    double result1 = std::log(0.9) - std::log(0.1) +
                     std::log(0.8) - std::log(0.2) +
                     std::log(0.3) - std::log(0.3) +
                     std::log(0.4) - std::log(0.4) +
                     std::log(0.5) - std::log(0.5) +
                     std::log(0.6) - std::log(0.6) +
                     std::log(0.2) - std::log(0.7) +
                     std::log(0.1) - std::log(0.8);
    double result2 = std::log(0.9) - std::log(0.1) +
                     std::log(0.8) - std::log(0.2) +
                     std::log(0.7) - std::log(0.3) +
                     std::log(0.4) - std::log(0.4) +
                     std::log(0.5) - std::log(0.5) +
                     std::log(0.6) - std::log(0.6) +
                     std::log(0.2) - std::log(0.7) +
                     std::log(0.1) - std::log(0.8);
    result = addLogProb(result1, result2);
    BOOST_CHECK(std::abs(config._tmpAttachmentScore.lcaRScore(0) - result)
                <= 100 * std::numeric_limits<double>::epsilon());

    Config<SampleTree>::TAttachmentScores &attachmentScores = config.getTmpAttachmentScore();
    static Config<SampleTree>::TAttachmentScores attachmentSumScores;
    attachmentSumScores.resize(attachmentScores.size());
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentScores;
    passDownAttachmentScores.resize(attachmentScores.size());
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentSumScores;
    passDownAttachmentSumScores.resize(attachmentScores.size());

    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;
    scoreSum = AttachmentScore();

    long double sumParallel = 0;
    PassScoreToChildrenBFSVisitor visBFS(config,
                                         attachmentScores,
                                         attachmentSumScores,
                                         passDownAttachmentScores,
                                         passDownAttachmentSumScores,
                                         sumParallel,
                                         0);
    breadth_first_search(config.getTree(), num_vertices(config.getTree()) - 1, visitor(visBFS));

    result = addLogProb(attachmentScores[0].lcaRScore(),
                        attachmentScores[1].lcaRScore());
    BOOST_CHECK(std::abs(result -
                         passDownAttachmentSumScores[6].paralleleScore())
                < 100.0 * std::numeric_limits<double>::epsilon());

    result1 = std::log(0.9) - std::log(0.1) +
                     std::log(0.8) - std::log(0.2) +
                     std::log(0.3) - std::log(0.3) +
                     std::log(0.4) - std::log(0.4) +
                     std::log(0.5) - std::log(0.5) +
                     std::log(0.6) - std::log(0.6) +
                     std::log(0.2) - std::log(0.7) +
                     std::log(0.1) - std::log(0.8);
    result2 = std::log(0.9) - std::log(0.1) +
                     std::log(0.8) - std::log(0.2) +
                     std::log(0.7) - std::log(0.3) +
                     std::log(0.4) - std::log(0.4) +
                     std::log(0.5) - std::log(0.5) +
                     std::log(0.6) - std::log(0.6) +
                     std::log(0.2) - std::log(0.7) +
                     std::log(0.1) - std::log(0.8);
    double result3 = std::log(0.9) - std::log(0.1) +
                     std::log(0.8) - std::log(0.2) +
                     std::log(0.3) - std::log(0.3) +
                     std::log(0.4) - std::log(0.4) +
                     std::log(0.4) - std::log(0.5) +
                     std::log(0.3) - std::log(0.6) +
                     std::log(0.7) - std::log(0.7) +
                     std::log(0.8) - std::log(0.8);

    result = addLogProb(addLogProb(result1, result2), result3);

    BOOST_CHECK(passDownAttachmentSumScores[6].paralleleScore() -
                result <=
                100 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(passDownAttachmentSumScores[5].paralleleScore() -
                result2 <=
                100 * std::numeric_limits<double>::epsilon());
//    BOOST_CHECK(result2 -
//                passDownAttachmentSumScores[9].paralleleScore() <=
//                100 * std::numeric_limits<double>::epsilon());

    std::cout << "S: << " << scoreSum.lcaRScore() << std::endl;
    std::cout << "0: " << config._tmpAttachmentScore.lcaRScore(0) << std::endl;
    std::cout << "2: " << passDownAttachmentScores[2].paralleleScore() << "\t"
              << passDownAttachmentSumScores[2].paralleleScore() << std::endl;
    std::cout << "13: " << passDownAttachmentScores[13].paralleleScore() << "\t"
              << passDownAttachmentSumScores[13].paralleleScore() << std::endl;
    std::cout << "14 " << passDownAttachmentScores[14].paralleleScore() << "\t"
              << passDownAttachmentSumScores[14].paralleleScore() << std::endl;
    std::cout << "5: " << passDownAttachmentScores[5].paralleleScore() << "\t"
              << passDownAttachmentSumScores[5].paralleleScore() << std::endl;
    std::cout << "9: " << passDownAttachmentScores[9].paralleleScore() << "\t"
              << passDownAttachmentSumScores[9].paralleleScore() << std::endl;
    std::cout << "6: " << passDownAttachmentScores[6].paralleleScore() << "\t"
              << passDownAttachmentSumScores[6].paralleleScore() << std::endl;
    std::cout << "7: " << passDownAttachmentScores[7].paralleleScore() << "\t"
              << passDownAttachmentSumScores[7].paralleleScore() << std::endl;
    std::cout << "8: " << passDownAttachmentScores[8].paralleleScore() << "\t"
              << passDownAttachmentSumScores[8].paralleleScore() << std::endl;


}


BOOST_AUTO_TEST_CASE(sample_tree_compute_score_lar_1) {
    Config<SampleTree> config;
    initTree1(config);
    config.computeLossScore = true;
    scoreTree(config);

    // test the case where the mutation is lost
    BOOST_CHECK_MESSAGE(config._tmpAttachmentScore.lossAltRScore(7) == -INFINITY,
                        "config._tmpAttachmentScore.lossAltRScore(7) != -inf - "
                                << config._tmpAttachmentScore.lossAltRScore(7));

    BOOST_CHECK_MESSAGE(config._tmpAttachmentScore.lossAltRScore(8) == -INFINITY,
                        "config._tmpAttachmentScore.lossAltRScore(8) != -inf - "
                                << config._tmpAttachmentScore.lossAltRScore(8));

    BOOST_CHECK_MESSAGE(config._tmpAttachmentScore.lossAltRScore(5) == -INFINITY,
                        "config._tmpAttachmentScore.lossAltRScore(5) != -inf - "
                                << config._tmpAttachmentScore.lossAltRScore(5));

    BOOST_CHECK_MESSAGE(config._tmpAttachmentScore.lossAltRScore(6) == -INFINITY,
                        "config._tmpAttachmentScore.lossAltRScore(6) != -inf - "
                                << config._tmpAttachmentScore.lossAltRScore(6));

    BOOST_CHECK_MESSAGE(config._tmpAttachmentScore.lossAltScore(6) == -INFINITY,
                        "config._tmpAttachmentScore.lossAltScore(6) != -inf - "
                                << config._tmpAttachmentScore.lossAltScore(6));

    BOOST_CHECK_MESSAGE(config._tmpAttachmentScore.lossAltRScore(3) == -INFINITY,
                        "config._tmpAttachmentScore.lossAltRScore(3) != -inf - "
                                << config._tmpAttachmentScore.lossAltRScore(3));

    double result = std::log(0.1) - std::log(0.1)
                    + std::log(0.2) - std::log(0.2)
                    + std::log(0.7) - std::log(0.3)
                    + std::log(0.6) - std::log(0.4);

    BOOST_CHECK_MESSAGE(std::abs(config._tmpAttachmentScore.lossAltRScore(1) - result)
                        <= 10 * std::numeric_limits<double>::epsilon(),
                        config._tmpAttachmentScore.lossAltRScore(1)
                                << " != " << result);

    double los5 = std::log(0.1) - std::log(0.1)
                  + std::log(0.2) - std::log(0.2)
                  + std::log(0.7) - std::log(0.3)
                  + std::log(0.6) - std::log(0.4)
                  + std::log(0.5) - std::log(0.6);

    double los3 = std::log(0.1) - std::log(0.1)
                  + std::log(0.2) - std::log(0.2)
                  + std::log(0.3) - std::log(0.3)
                  + std::log(0.6) - std::log(0.4)
                  + std::log(0.5) - std::log(0.6);

    result = addLogProb(los5, los3);

    BOOST_CHECK_MESSAGE(std::abs(config._tmpAttachmentScore.lossAltRScore(0) - result)
                        <= 10 * std::numeric_limits<double>::epsilon(),
                        config._tmpAttachmentScore.lossAltRScore(0)
                                << " != " << result);
}

BOOST_AUTO_TEST_CASE(sample_tree_compute_score_lossInCurrrentNode_1) {
    Config<SampleTree> config;
    initTree1(config);
    config.computeLossScore = true;
    scoreTree(config);

    Config<SampleTree>::TAttachmentScores &attachmentScores = config.getTmpAttachmentScore();
    static Config<SampleTree>::TAttachmentScores attachmentSumScores;
    attachmentSumScores.resize(attachmentScores.size());
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentScores;
    passDownAttachmentScores.resize(attachmentScores.size());
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentSumScores;
    passDownAttachmentSumScores.resize(attachmentScores.size());

    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;
    scoreSum = AttachmentScore();

    long double sumParallel = 0;
    PassScoreToChildrenBFSVisitor visBFS(config,
                                         attachmentScores,
                                         attachmentSumScores,
                                         passDownAttachmentScores,
                                         passDownAttachmentSumScores,
                                         sumParallel,
                                         0);
    breadth_first_search(config.getTree(), num_vertices(config.getTree()) - 1, visitor(visBFS));

    // test the case where the mutation is lost
    BOOST_CHECK_MESSAGE(passDownAttachmentScores.lossAltInCurrentNodeScore(0) == -INFINITY,
                        "passDownAttachmentScores.lossAltInCurrentNodeScore(0) != -inf - "
                                << passDownAttachmentScores.lossAltInCurrentNodeScore(0));

    BOOST_CHECK_MESSAGE(passDownAttachmentScores.lossAltInCurrentNodeScore(8) == -INFINITY,
                        "passDownAttachmentScores.lossAltInCurrentNodeScore(8) != -inf - "
                                << passDownAttachmentScores.lossAltInCurrentNodeScore(8));

    BOOST_CHECK_MESSAGE(passDownAttachmentScores.lossAltInCurrentNodeScore(7) == -INFINITY,
                        "passDownAttachmentScores.lossAltInCurrentNodeScore(7) != -inf - "
                                << passDownAttachmentScores.lossAltInCurrentNodeScore(7));

    BOOST_CHECK_MESSAGE(passDownAttachmentScores.lossAltInCurrentNodeScore(4) == -INFINITY,
                        "passDownAttachmentScores.lossAltInCurrentNodeScore(4) != -inf - "
                                << passDownAttachmentScores.lossAltInCurrentNodeScore(4));

    double result = std::log(0.1) - std::log(0.1)
                    + std::log(0.2) - std::log(0.2)
                    + std::log(0.3) - std::log(0.3)
                    + std::log(0.4) - std::log(0.4)
                    + std::log(0.5) - std::log(0.6);

    BOOST_CHECK_MESSAGE(std::abs(passDownAttachmentScores.lossAltInCurrentNodeScore(1) - result)
                        <= 10.0 * std::numeric_limits<double>::epsilon(),
                        passDownAttachmentScores.lossAltInCurrentNodeScore(1)
                                << " != " << result);
    double gain0 = std::log(0.1) - std::log(0.1)
                   + std::log(0.2) - std::log(0.2)
                   + std::log(0.3) - std::log(0.3)
                   + std::log(0.6) - std::log(0.4)
                   + std::log(0.5) - std::log(0.6);

    double gain1 = std::log(0.1) - std::log(0.1)
                   + std::log(0.2) - std::log(0.2)
                   + std::log(0.3) - std::log(0.3)
                   + std::log(0.6) - std::log(0.4);

    result = addLogProb(gain0, gain1);

    BOOST_CHECK_MESSAGE(std::abs(passDownAttachmentScores.lossAltInCurrentNodeScore(2) - result)
                        <= 10.0 * std::numeric_limits<double>::epsilon(),
                        passDownAttachmentScores.lossAltInCurrentNodeScore(2)
                                << " != " << result);

    gain0 = std::log(0.1) - std::log(0.1)
            + std::log(0.2) - std::log(0.2)
            + std::log(0.7) - std::log(0.3)
            + std::log(0.6) - std::log(0.4)
            + std::log(0.5) - std::log(0.6);

    gain1 = std::log(0.1) - std::log(0.1)
            + std::log(0.2) - std::log(0.2)
            + std::log(0.7) - std::log(0.3)
            + std::log(0.6) - std::log(0.4);

    result = addLogProb(gain0, gain1);

    double gain3 = std::log(0.1) - std::log(0.1)
                   + std::log(0.2) - std::log(0.2)
                   + std::log(0.7) - std::log(0.3);

    result = addLogProb(result, gain3);

    BOOST_CHECK_MESSAGE(std::abs(passDownAttachmentScores.lossAltInCurrentNodeScore(3) - result)
                        <= 10.0 * std::numeric_limits<double>::epsilon(),
                        passDownAttachmentScores.lossAltInCurrentNodeScore(3)
                                << " != " << result);
}

BOOST_AUTO_TEST_CASE(sample_tree_compute_score_lossInCurrrentNodeR_1) {
    Config<SampleTree> config;
    initTree1(config);
    config.computeLossScore = true;
    scoreTree(config);

    Config<SampleTree>::TAttachmentScores &attachmentScores = config.getTmpAttachmentScore();
    static Config<SampleTree>::TAttachmentScores attachmentSumScores;
    attachmentSumScores.resize(attachmentScores.size());
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentScores;
    passDownAttachmentScores.resize(attachmentScores.size());
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentSumScores;
    passDownAttachmentSumScores.resize(attachmentScores.size());

    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;
    scoreSum = AttachmentScore();
    long double sumParallel = 0;
    PassScoreToChildrenBFSVisitor visBFS(config,
                                         attachmentScores,
                                         attachmentSumScores,
                                         passDownAttachmentScores,
                                         passDownAttachmentSumScores,
                                         sumParallel,
                                         //scoreSum,
                                         0);
    breadth_first_search(config.getTree(), num_vertices(config.getTree()) - 1, visitor(visBFS));

    // test the case where the mutation is lost
    BOOST_CHECK_MESSAGE(passDownAttachmentScores.lossAltInCurrentNodeRScore(0) == -INFINITY,
                        "passDownAttachmentScores.lossAltInCurrentNodeRScore(0) != -inf - "
                                << passDownAttachmentScores.lossAltInCurrentNodeRScore(0));

    BOOST_CHECK_MESSAGE(passDownAttachmentScores.lossAltInCurrentNodeRScore(1) == -INFINITY,
                        "passDownAttachmentScores.lossAltInCurrentNodeRScore(1) != -inf - "
                                << passDownAttachmentScores.lossAltInCurrentNodeRScore(1));

    BOOST_CHECK_MESSAGE(passDownAttachmentScores.lossAltInCurrentNodeRScore(8) == -INFINITY,
                        "passDownAttachmentScores.lossAltInCurrentNodeRScore(8) != -inf - "
                                << passDownAttachmentScores.lossAltInCurrentNodeRScore(8));

    BOOST_CHECK_MESSAGE(passDownAttachmentScores.lossAltInCurrentNodeRScore(7) == -INFINITY,
                        "passDownAttachmentScores.lossAltInCurrentNodeRScore(7) != -inf - "
                                << passDownAttachmentScores.lossAltInCurrentNodeRScore(7));

    BOOST_CHECK_MESSAGE(passDownAttachmentScores.lossAltInCurrentNodeRScore(4) == -INFINITY,
                        "passDownAttachmentScores.lossAltInCurrentNodeRScore(4) != -inf - "
                                << passDownAttachmentScores.lossAltInCurrentNodeRScore(4));

    double gain0 = std::log(0.1) - std::log(0.1)
                   + std::log(0.2) - std::log(0.2)
                   + std::log(0.3) - std::log(0.3)
                   + std::log(0.6) - std::log(0.4)
                   + std::log(0.5) - std::log(0.6);

    double result = gain0;

    BOOST_CHECK_MESSAGE(std::abs(passDownAttachmentScores.lossAltInCurrentNodeRScore(2) - result)
                        <= 10.0 * std::numeric_limits<double>::epsilon(),
                        passDownAttachmentScores.lossAltInCurrentNodeRScore(2)
                                << " != " << result);

    gain0 = std::log(0.1) - std::log(0.1)
            + std::log(0.2) - std::log(0.2)
            + std::log(0.7) - std::log(0.3)
            + std::log(0.6) - std::log(0.4)
            + std::log(0.5) - std::log(0.6);

    double gain1 = std::log(0.1) - std::log(0.1)
                   + std::log(0.2) - std::log(0.2)
                   + std::log(0.7) - std::log(0.3)
                   + std::log(0.6) - std::log(0.4);

    result = addLogProb(gain0, gain1);

    BOOST_CHECK_MESSAGE(std::abs(passDownAttachmentScores.lossAltInCurrentNodeRScore(3) - result)
                        <= 10.0 * std::numeric_limits<double>::epsilon(),
                        passDownAttachmentScores.lossAltInCurrentNodeRScore(3)
                                << " != " << result);
}

BOOST_AUTO_TEST_CASE(sample_tree_het_prob_1) {
    Config<SampleTree> config;
    initTree1(config);
    config.computeLossScore = true;
    scoreTree(config);

    Config<SampleTree>::TAttachmentScores &attachmentScores = config.getTmpAttachmentScore();
    static Config<SampleTree>::TAttachmentScores attachmentSumScores;
    attachmentSumScores.resize(attachmentScores.size());
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentScores;
    passDownAttachmentScores.resize(attachmentScores.size());
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentSumScores;
    passDownAttachmentSumScores.resize(attachmentScores.size());

    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;
    scoreSum = AttachmentScore();
    long double sumParallel = 0;
    PassScoreToChildrenBFSVisitor visBFS(config,
                                         attachmentScores,
                                         attachmentSumScores,
                                         passDownAttachmentScores,
                                         passDownAttachmentSumScores,
                                         sumParallel,
                                         //scoreSum,
                                         0);
    breadth_first_search(config.getTree(), num_vertices(config.getTree()) - 1, visitor(visBFS));

    // test the case where the mutation is lost
    double allMutated = std::log(0.9) - std::log(0.1)
                        + std::log(0.8) - std::log(0.2)
                        + std::log(0.7) - std::log(0.3)
                        + std::log(0.6) - std::log(0.4)
                        + std::log(0.5) - std::log(0.6);

    double result = addLogProb(allMutated, std::log(0.5) - std::log(0.6));

    BOOST_CHECK_MESSAGE(std::abs(attachmentSumScores.hetScore(8) - result)
                        <= 10.0 * std::numeric_limits<double>::epsilon(),
                        attachmentSumScores.hetScore(8)
                                << " != " << result);

    double allBut2 = std::log(0.9) - std::log(0.1)
                     + std::log(0.8) - std::log(0.2)
                     + std::log(0.7) - std::log(0.3)
                     + std::log(0.6) - std::log(0.4)
                     + std::log(0.6) - std::log(0.6);

    result = addLogProb(allMutated, allBut2);
    result = addLogProb(result, std::log(0.6) - std::log(0.4));

    BOOST_CHECK_MESSAGE(std::abs(attachmentSumScores.hetScore(7) - result)
                        <= 10.0 * std::numeric_limits<double>::epsilon(),
                        attachmentSumScores.hetScore(7)
                                << " != " << result);

    double allBut2_4 = std::log(0.9) - std::log(0.1)
                       + std::log(0.8) - std::log(0.2)
                       + std::log(0.7) - std::log(0.3)
                       + std::log(0.4) - std::log(0.4)
                       + std::log(0.6) - std::log(0.6);

    double allBut2_4_6 = std::log(0.9) - std::log(0.1)
                         + std::log(0.8) - std::log(0.2)
                         + std::log(0.3) - std::log(0.3)
                         + std::log(0.4) - std::log(0.4)
                         + std::log(0.6) - std::log(0.6);

    double allBut2_4_6_8 = std::log(0.9) - std::log(0.1)
                           + std::log(0.2) - std::log(0.2)
                           + std::log(0.3) - std::log(0.3)
                           + std::log(0.4) - std::log(0.4)
                           + std::log(0.6) - std::log(0.6);

    result = addLogProb(allMutated, allBut2);
    result = addLogProb(result, allBut2_4);
    result = addLogProb(result, allBut2_4_6);
    result = addLogProb(result, allBut2_4_6_8);

    BOOST_CHECK_MESSAGE(std::abs(attachmentSumScores.hetScore(4) - result)
                        <= 10.0 * std::numeric_limits<double>::epsilon(),
                        attachmentSumScores.hetScore(4)
                                << " != " << result);
}

BOOST_AUTO_TEST_CASE(sample_tree_hom_prob_1) {
    Config<SampleTree> config;
    initTree1(config);
    config.computeLossScore = true;
    scoreTree(config);

    Config<SampleTree>::TAttachmentScores &attachmentScores = config.getTmpAttachmentScore();
    static Config<SampleTree>::TAttachmentScores attachmentSumScores;
    attachmentSumScores.resize(attachmentScores.size());
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentScores;
    passDownAttachmentScores.resize(attachmentScores.size());
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentSumScores;
    passDownAttachmentSumScores.resize(attachmentScores.size());

    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;
    scoreSum = AttachmentScore();
    long double sumParallel = 0;
    PassScoreToChildrenBFSVisitor visBFS(config,
                                         attachmentScores,
                                         attachmentSumScores,
                                         passDownAttachmentScores,
                                         passDownAttachmentSumScores,
                                         sumParallel,
                                         //scoreSum,
                                         0);
    breadth_first_search(config.getTree(), num_vertices(config.getTree()) - 1, visitor(visBFS));

    // test the case where the mutation is lost
    double allMutated = std::log(0.5) - std::log(0.1)
                        + std::log(0.4) - std::log(0.2)
                        + std::log(0.3) - std::log(0.3)
                        + std::log(0.2) - std::log(0.4)
                        + std::log(0.1) - std::log(0.6);

    double result = allMutated;

    BOOST_CHECK_MESSAGE(std::abs(attachmentSumScores.homScore(8) - result)
                        <= 10.0 * std::numeric_limits<double>::epsilon(),
                        attachmentSumScores.homScore(8)
                                << " != " << result);

    double allBut2 = std::log(0.5) - std::log(0.1)
                     + std::log(0.4) - std::log(0.2)
                     + std::log(0.3) - std::log(0.3)
                     + std::log(0.2) - std::log(0.4)
                     + std::log(0.6) - std::log(0.6);

    result = addLogProb(allMutated, allBut2);

    BOOST_CHECK_MESSAGE(std::abs(attachmentSumScores.homScore(7) - result)
                        <= 10.0 * std::numeric_limits<double>::epsilon(),
                        attachmentSumScores.homScore(7)
                                << " != " << result);

    double allBut2_4 = std::log(0.5) - std::log(0.1)
                       + std::log(0.4) - std::log(0.2)
                       + std::log(0.3) - std::log(0.3)
                       + std::log(0.4) - std::log(0.4)
                       + std::log(0.6) - std::log(0.6);

    double allBut2_4_6 = std::log(0.5) - std::log(0.1)
                         + std::log(0.4) - std::log(0.2)
                         + std::log(0.3) - std::log(0.3)
                         + std::log(0.4) - std::log(0.4)
                         + std::log(0.6) - std::log(0.6);

    //double allBut2_4_6_8 = std::log(0.5) - std::log(0.1)
    //                       + std::log(0.2) - std::log(0.2)
    //                       + std::log(0.3) - std::log(0.3)
    //                       + std::log(0.4) - std::log(0.4)
    //                       + std::log(0.6) - std::log(0.6);

    result = addLogProb(allMutated, allBut2);
    result = addLogProb(result, allBut2_4);
    result = addLogProb(result, allBut2_4_6);

    BOOST_CHECK_MESSAGE(std::abs(attachmentSumScores.homScore(4) - result)
                        <= 10.0 * std::numeric_limits<double>::epsilon(),
                        attachmentSumScores.homScore(4)
                                << " != " << result);
}

BOOST_AUTO_TEST_CASE(sample_tree_lcaR_prob_2) {
    Config<SampleTree> config;
    initTree2(config);
    config.computeLossScore = true;
    config.computeParallelScore = true;
    scoreTree(config);

    Config<SampleTree>::TAttachmentScores &attachmentScores = config.getTmpAttachmentScore();
    static Config<SampleTree>::TAttachmentScores attachmentSumScores;
    attachmentSumScores.resize(attachmentScores.size());
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentScores;
    passDownAttachmentScores.resize(attachmentScores.size());
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentSumScores;
    passDownAttachmentSumScores.resize(attachmentScores.size());

    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;
    scoreSum = AttachmentScore();
    long double sumParallel = 0;
    PassScoreToChildrenBFSVisitor visBFS(config,
                                         attachmentScores,
                                         attachmentSumScores,
                                         passDownAttachmentScores,
                                         passDownAttachmentSumScores,
                                         sumParallel,
                                         //scoreSum,
                                         0);
    breadth_first_search(config.getTree(), num_vertices(config.getTree()) - 1, visitor(visBFS));

    BOOST_CHECK(passDownAttachmentSumScores[0].sibNodeScore() == -INFINITY);
    BOOST_CHECK(passDownAttachmentSumScores[1].sibNodeScore() == config._tmpAttachmentScore[2].hetScore());
    BOOST_CHECK(passDownAttachmentSumScores[2].sibNodeScore() == config._tmpAttachmentScore[1].hetSumScore());
    BOOST_CHECK(passDownAttachmentSumScores[3].sibNodeScore() ==
                addLogProb(config._tmpAttachmentScore[2].hetScore(), config._tmpAttachmentScore[4].hetScore()));
    BOOST_CHECK(passDownAttachmentSumScores[4].sibNodeScore() ==
                addLogProb(config._tmpAttachmentScore[2].hetScore(), config._tmpAttachmentScore[3].hetSumScore()));
    BOOST_CHECK(std::abs(passDownAttachmentSumScores[5].sibNodeScore() -
                         addLogProb(config._tmpAttachmentScore[2].hetScore(),
                                    addLogProb(config._tmpAttachmentScore[4].hetScore(),
                                               config._tmpAttachmentScore[6].hetScore()))) <
                10.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(std::abs(passDownAttachmentSumScores[6].sibNodeScore() -
                         addLogProb(config._tmpAttachmentScore[2].hetScore(),
                                    addLogProb(config._tmpAttachmentScore[4].hetScore(),
                                               config._tmpAttachmentScore[5].hetScore()))) <
                10.0 * std::numeric_limits<double>::epsilon());

    BOOST_CHECK(std::abs(passDownAttachmentSumScores[2].paralleleScore() - attachmentScores[0].lcaRScore()
    ) < 10.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(passDownAttachmentScores[3].paralleleScore() == -INFINITY);
    BOOST_CHECK(passDownAttachmentSumScores[3].paralleleScore() == -INFINITY);

    double result = std::log(0.9) - std::log(0.1) +
                    std::log(0.8) - std::log(0.2) +
                    std::log(0.3) - std::log(0.7) +
                    std::log(0.2) - std::log(0.8);

    BOOST_CHECK(std::abs(passDownAttachmentSumScores[5].paralleleScore() - result)
                < 100.0 * std::numeric_limits<double>::epsilon());

    result = std::log(0.7) - std::log(0.3) +
             std::log(0.5) - std::log(0.4) +
             std::log(0.3) - std::log(0.7) +
             std::log(0.2) - std::log(0.8);

    // TODO, there is quite some imprecision involved
    BOOST_CHECK(std::abs(passDownAttachmentSumScores[6].paralleleScore() - result)
                < 1000.0 * std::numeric_limits<double>::epsilon());

    BOOST_CHECK(std::abs(attachmentScores[0].lcaRScore() -
                         addLogProb(passDownAttachmentSumScores[5].paralleleScore(),
                                    passDownAttachmentSumScores[6].paralleleScore()))
                < 100.0 * std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(sample_tree_lcaR_prob_3) {
    Config<SampleTree> config;
    initTree3(config);
    config.computeLossScore = true;
    config.computeParallelScore = true;
    scoreTree(config);

    Config<SampleTree>::TAttachmentScores &attachmentScores = config.getTmpAttachmentScore();
    static Config<SampleTree>::TAttachmentScores attachmentSumScores;
    attachmentSumScores.resize(attachmentScores.size());
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentScores;
    passDownAttachmentScores.resize(attachmentScores.size());
    static Config<SampleTree>::TPassDownAttachmentScores passDownAttachmentSumScores;
    passDownAttachmentSumScores.resize(attachmentScores.size());

    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;
    scoreSum = AttachmentScore();
    long double sumParallel = 0;
    PassScoreToChildrenBFSVisitor visBFS(config,
                                         attachmentScores,
                                         attachmentSumScores,
                                         passDownAttachmentScores,
                                         passDownAttachmentSumScores,
                                         sumParallel,
                                         //scoreSum,
                                         0);
    breadth_first_search(config.getTree(), num_vertices(config.getTree()) - 1, visitor(visBFS));

    BOOST_CHECK(std::abs(attachmentScores[0].lcaRScore() -
                         passDownAttachmentSumScores[4].paralleleScore())
                < 100.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(std::abs(passDownAttachmentSumScores[4].paralleleScore() -
                         passDownAttachmentSumScores[2].paralleleScore())
                < 100.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(std::abs(passDownAttachmentSumScores[5].paralleleScore() -
                         passDownAttachmentSumScores[2].paralleleScore())
                < 100.0 * std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(sample_tree_compute_mut_type_contribution1) {

    Config<SampleTree> config;
    initTree2(config);
    config.learnZygocity = false;
    config.computeLossScore = false;
    config.computeParallelScore = false;
    scoreTree(config);

    for (unsigned i = 0; i < 8; ++i) {
        std::get<0>(config.logScores).wtScore(i, 0) = std::log(1);
        std::get<0>(config.logScores).hetScore(i, 0) = std::log(1);
        std::get<0>(config.logScores).homScore(i, 0) = std::log(1);
    }

    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;
    scoreSum = AttachmentScore();

    typename Config<SampleTree>::TAttachmentScores &attachmentScores = config.getTmpAttachmentScore();
    getAllAttachmentScores(attachmentScores, scoreSum, config, 0);
    std::array<double, 6> res;

    res = scoreSum.computeMutTypeContribution(config, scoreSum);

    BOOST_CHECK(res[0] <= 10.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(res[1] == -INFINITY);
    BOOST_CHECK(res[2] == -INFINITY);
    BOOST_CHECK(res[3] == -INFINITY);
    BOOST_CHECK(res[4] == -INFINITY);
    BOOST_CHECK(res[5] <= 10.0 * std::numeric_limits<double>::epsilon());

    config.learnZygocity = true;
    getAllAttachmentScores(attachmentScores, scoreSum, config, 0);
    res = scoreSum.computeMutTypeContribution(config, scoreSum);

    // res[x] * #possibilities == #possibilities
    BOOST_CHECK(res[0] <= 10.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(res[1] <= 10.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(res[2] == -INFINITY);
    BOOST_CHECK(res[3] == -INFINITY);
    BOOST_CHECK(res[4] == -INFINITY);
    BOOST_CHECK(res[5] <= 10.0 * std::numeric_limits<double>::epsilon());

    config.computeLossScore = true;
    getAllAttachmentScores(attachmentScores, scoreSum, config, 0);
    res = scoreSum.computeMutTypeContribution(config, scoreSum);
    // res[x] * #possibilities == #possibilities
    BOOST_CHECK(res[0] <= 10.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(res[1] <= 10.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(res[2] <= 10.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(res[3] <= 10.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(res[4] == -INFINITY);
    BOOST_CHECK(res[5] <= 10.0 * std::numeric_limits<double>::epsilon());

    config.computeParallelScore = true;
    getAllAttachmentScores(attachmentScores, scoreSum, config, 0);
    res = scoreSum.computeMutTypeContribution(config, scoreSum);
    // res[x] * #possibilities == #possibilities
    BOOST_CHECK(res[0] <= 10.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(res[1] <= 10.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(res[2] <= 10.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(res[3] <= 10.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(res[4] <= 10.0 * std::numeric_limits<double>::epsilon());
    BOOST_CHECK(res[5] <= 10.0 * std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(sample_tree_compute_mut_type_contribution2_2) {

    Config<SampleTree> config;
    initTree2_1(config);
    config.learnZygocity = true;
    config.computeLossScore = true;
    config.computeParallelScore = true;
    scoreTree(config);

    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;

    typename Config<SampleTree>::TAttachmentScores &attachmentScores = config.getTmpAttachmentScore();
    getAllAttachmentScores(attachmentScores, scoreSum, config, 0);
    std::array<double, 6> res;

    for (unsigned i = 0; i < attachmentScores.size(); ++i)
    {
        std::cout << i << std::endl;
        std::cout << attachmentScores[i] << std::endl;
    }

    std::cout << scoreSum << std::endl;

    res = scoreSum.computeMutTypeContribution(config, scoreSum);
    for (unsigned i = 0; i < 6; ++i)
        std::cout << res[i] << " ";
    std::cout << std::endl;

}

BOOST_AUTO_TEST_CASE(sample_tree_compute_mut_type_contribution3_2) {

    /*
    *                                          11
    *                                          |
    *                                          0
    *                           ___________________________
    *                          |                           |
    *                          1                           2
    *               _______________________           ____________
    *              |                       |         |            |
    *              3                       8         9            10
    *       _______________                0.4       0.5          0.6
    *      |               |               0.6       0.4          0.3
    *      4               7               0.2       0.1          0.7
    *   _______            0.3
    *  |       |           0.7
    *  5       6           0.3
    *  0.1     0.2
    *  0.9     0.8
    *  0.5     0.4
    */
    Config<SampleTree> config;
    initTree3_1(config);
    config.learnZygocity = true;
    config.computeLossScore = true;
    config.computeParallelScore = true;
    scoreTree(config);

    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;

    typename Config<SampleTree>::TAttachmentScores &attachmentScores = config.getTmpAttachmentScore();
    getAllAttachmentScores(attachmentScores, scoreSum, config, 0);
    std::array<double, 6> res;

    for (unsigned i = 0; i < attachmentScores.size(); ++i)
    {
        std::cout << i << std::endl;
        std::cout << attachmentScores[i] << std::endl;
    }

    std::cout << scoreSum << std::endl;

    res = scoreSum.computeMutTypeContribution(config, scoreSum);
    for (unsigned i = 0; i < 6; ++i)
        std::cout << res[i] << " ";
    std::cout << std::endl;

    std::cout << "score: " << scoreTree(config) << std::endl;

}

BOOST_AUTO_TEST_CASE(sample_tree_compute_mut_type_het_2) {

    Config<SampleTree> config;
    initTree2(config);
    config.learnZygocity = true;
    config.computeLossScore = true;
    config.computeParallelScore = true;
    scoreTree(config);

    for (unsigned i = 0; i < 4; ++i) {
        std::get<0>(config.logScores).wtScore(i, 0) = std::log(0.01);
        std::get<0>(config.logScores).hetScore(i, 0) = std::log(1);
        std::get<0>(config.logScores).homScore(i, 0) = std::log(0.01);
    }
    for (unsigned i = 4; i < 8; ++i) {
        std::get<0>(config.logScores).wtScore(i, 0) = std::log(1);
        std::get<0>(config.logScores).hetScore(i, 0) = std::log(0.01);
        std::get<0>(config.logScores).homScore(i, 0) = std::log(0.01);
    }

    Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;
    scoreSum = AttachmentScore();

    typename Config<SampleTree>::TAttachmentScores &attachmentScores = config.getTmpAttachmentScore();
    getAllAttachmentScores(attachmentScores, scoreSum, config, 0);
    std::array<double, 6> res;

    //std::cout << config.getParam(Config<SampleTree>::E_lambda) << std::endl;
    //std::cout << config.getParam(Config<SampleTree>::E_parallel) << std::endl;

    double result = std::log(1) - std::log(0.01) +
                    std::log(1) - std::log(0.01) +
                    std::log(0.01) - std::log(0.01) +
                    std::log(0.01) - std::log(0.01) +
                    std::log(1) - std::log(1) +
                    std::log(1) - std::log(1) +
                    std::log(0.01) - std::log(1) +
                    std::log(0.01) - std::log(1);
    result = addLogProb(result, result);

    res = scoreSum.computeMutTypeContribution(config, scoreSum);

    std::cout << attachmentScores.lcaRScore(0) << " " << result << std::endl;
    BOOST_CHECK(attachmentScores.lcaRScore(0) - result <= 10.0 * std::numeric_limits<double>::epsilon());

    std::cout << std::exp(res[0] - res[4]) << " " <<
              std::exp(res[1] - res[4]) << " " <<
              std::exp(res[2] - res[4]) << " " <<
              std::exp(res[3] - res[4]) << " " <<
              std::exp(res[4] - res[4]) << std::endl;
}


BOOST_AUTO_TEST_CASE(sample_tree_simplify_tree) {
    Config<SampleTree> config;
    config.setNumSamples(3);
    createInitialTreeTest(config);

    config.data.resize(3);
    config.data[0].resize(3, std::tuple<unsigned, unsigned>{10, 0});
    config.data[1].resize(3, std::tuple<unsigned, unsigned>{10, 0});
    config.data[2].resize(3, std::tuple<unsigned, unsigned>{10, 0});

    // cell  1, 2 and 3 have mutation 1
    std::get<1>(config.data[0][0]) = 10;
    std::get<1>(config.data[1][0]) = 10;
    std::get<1>(config.data[2][0]) = 10;

    // cell 2 has mutation 2
    std::get<1>(config.data[1][1]) = 10;

    // cell 2 has mutation 3
    std::get<1>(config.data[1][2]) = 10;

    config.getLogScores().resizeNumCells(config.getNumSamples());
    config.getTmpLogScores().resizeNumCells(config.getNumSamples());
    config.getLogScores().resizeNumMuts(config.getNumMutations());
    config._tmpAttachmentScore.resize(5);
    computeLogScoresOP(config);


    typedef boost::adjacency_list<boost::vecS, boost::vecS,
            boost::bidirectionalS, Vertex<SimpleTree>> TGraph;
    TGraph newTree = simplifyTree(config);


    BOOST_CHECK_MESSAGE(num_vertices(newTree) == 4, "num_vertices(newTree) = " << num_vertices(newTree) << " != 4");
    BOOST_CHECK_MESSAGE(newTree[0].mutations[0] == 0,
                        "newTree[0].mutations[0] = " << newTree[0].mutations[0] << " != 0");
    BOOST_CHECK_MESSAGE(newTree[1].mutations[0] == 1,
                        "newTree[1].mutations[0] = " << newTree[1].mutations[0] << " != 1");
    BOOST_CHECK_MESSAGE(newTree[1].mutations[1] == 2,
                        "newTree[1].mutations[1] = " << newTree[1].mutations[1] << " != 2");
}

BOOST_AUTO_TEST_CASE(read_data_skip_indels) {
    std::string testNucs = ".,..,";
    BOOST_CHECK_MESSAGE(skipIndels(testNucs, 3) == 3, "skipIndels = " << skipIndels(testNucs, 3) << " != 3");
    testNucs = ".,.-3ACG.,";
    BOOST_CHECK_MESSAGE(skipIndels(testNucs, 3) == 7, "skipIndels = " << skipIndels(testNucs, 3) << " != 7");
    testNucs = ".,.+13ACGGGGGGGGGGG.,";
    BOOST_CHECK_MESSAGE(skipIndels(testNucs, 3) == 18, "skipIndels = " << skipIndels(testNucs, 3) << " != 18");
}

BOOST_AUTO_TEST_CASE(read_data_extract_seq_information) {
    std::vector<std::string> testLine = {"chr", "13", "REF", "5", "...,,^A", "IIIII", "5", "AC$GTN", "IIIII", "5",
                                         ".-3AACA.T.", "IIIII"};
    std::vector<std::array<unsigned, 5>> counts;
    counts.resize(3);
    counts[0] = {{0, 0, 0, 0, 0}};
    extractSeqInformation(counts[0], testLine, 0);
    BOOST_CHECK_MESSAGE(counts[0][0] == 0, " = counts[1][0] " << counts[0][0] << " != 0");
    BOOST_CHECK_MESSAGE(counts[0][1] == 0, " = counts[1][0] " << counts[0][1] << " != 0");
    BOOST_CHECK_MESSAGE(counts[0][2] == 0, " = counts[1][0] " << counts[0][2] << " != 0");
    BOOST_CHECK_MESSAGE(counts[0][3] == 0, " = counts[1][0] " << counts[0][3] << " != 0");
    BOOST_CHECK_MESSAGE(counts[0][4] == 5, " = counts[1][0] " << counts[0][4] << " != 5");

    counts[1] = {{0, 0, 0, 0, 0}};
    extractSeqInformation(counts[1], testLine, 1);
    BOOST_CHECK_MESSAGE(counts[1][0] == 1, " = counts[1][0] " << counts[1][0] << " != 1");
    BOOST_CHECK_MESSAGE(counts[1][1] == 1, " = counts[1][0] " << counts[1][1] << " != 1");
    BOOST_CHECK_MESSAGE(counts[1][2] == 1, " = counts[1][0] " << counts[1][2] << " != 1");
    BOOST_CHECK_MESSAGE(counts[1][3] == 1, " = counts[1][0] " << counts[1][3] << " != 1");
    BOOST_CHECK_MESSAGE(counts[1][4] == 4, " = counts[1][0] " << counts[1][4] << " != 4");

    counts[2] = {{0, 0, 0, 0, 0}};
    extractSeqInformation(counts[2], testLine, 2);
    BOOST_CHECK_MESSAGE(counts[2][0] == 1, " = counts[1][0] " << counts[2][0] << " != 1");
    BOOST_CHECK_MESSAGE(counts[2][1] == 0, " = counts[1][0] " << counts[2][1] << " != 1");
    BOOST_CHECK_MESSAGE(counts[2][2] == 0, " = counts[1][0] " << counts[2][2] << " != 1");
    BOOST_CHECK_MESSAGE(counts[2][3] == 1, " = counts[1][0] " << counts[2][3] << " != 1");
    BOOST_CHECK_MESSAGE(counts[2][4] == 5, " = counts[1][0] " << counts[2][4] << " != 5");
}

BOOST_AUTO_TEST_CASE(read_data_seq_error_stats) {
    std::vector<std::array<unsigned, 5>> counts;
    counts.resize(3);
    counts[0] = {{1, 0, 1, 0, 5}};
    counts[1] = {{1, 1, 0, 0, 4}};
    counts[2] = {{0, 0, 0, 3, 6}};

    unsigned errors = 0;
    unsigned cov = 0;
    updateSeqErrorStats(errors, cov, counts);

    BOOST_CHECK_MESSAGE(errors == 7, "errors: " << errors << " != 7");
    BOOST_CHECK_MESSAGE(cov == 15, "cov: " << cov << " != 15");
}

BOOST_AUTO_TEST_CASE(read_data_log_n_choose_k) {

    std::cout.precision(15);

    BOOST_CHECK_MESSAGE(std::round(exp(logNChooseK(10, 0, 0))) == 1, exp(logNChooseK(10, 0, 1)) << " != " << 1);
    BOOST_CHECK_MESSAGE(std::round(exp(logNChooseK(10, 0, log(2)))) == 1,
                        exp(logNChooseK(10, 0, log(2))) << " != " << 1);
    BOOST_CHECK_MESSAGE(std::round(exp(logNChooseK(10, 10, log(10)))) == 1,
                        exp(logNChooseK(10, 10, log(10))) << " != " << 1);
    BOOST_CHECK_MESSAGE(std::round(exp(logNChooseK(10, 6, log(252)))) == 210,
                        exp(logNChooseK(10, 6, log(252))) << " != " << 210);

    BOOST_CHECK_MESSAGE(std::round(exp(logNChoose2(2))) == 1, exp(logNChoose2(2)) << " != " << 1);
    BOOST_CHECK_MESSAGE(std::round(exp(logNChoose2(3))) == 3, exp(logNChoose2(3)) << " != " << 3);
    BOOST_CHECK_MESSAGE(std::round(exp(logNChoose2(10))) == 45, exp(logNChoose2(10)) << " != " << 45);
}

BOOST_AUTO_TEST_CASE(read_data_pass_cov_filter) {
    std::vector<std::string> testLine = {"chr", "13", "REF", "2", "..", "II", "2", "AC", "II", "5", ".-3AACA.T.",
                                         "IIIII"};
    std::vector<std::array<unsigned, 5>> counts;
    counts.resize(3);
    std::vector<unsigned> pos = {0, 1, 2};

    extractSeqInformation(counts, testLine, pos);

    BOOST_CHECK_MESSAGE(passCovFilter(counts[0][4], 3) == false,
                        "passCovFilter(counts[0][4], 3) = " << passCovFilter(counts[0][4], 3));
    BOOST_CHECK_MESSAGE(passCovFilter(counts[1][4], 3) == false,
                        "passCovFilter(counts[1][4], 3) = " << passCovFilter(counts[1][4], 3));
    BOOST_CHECK_MESSAGE(passCovFilter(counts[2][4], 3) == true,
                        "passCovFilter(counts[2][4], 3) = " << passCovFilter(counts[2][4], 3));
}

BOOST_AUTO_TEST_CASE(read_data_pass_supp_filter) {
    std::vector<std::string> testLine = {"chr", "13", "REF", "6", "..ACGG", "IIIIII", "2", "AC", "II", "5",
                                         ".-3AACAAT.",
                                         "IIIII"};
    std::vector<std::array<unsigned, 5>> counts;
    counts.resize(3);

    std::vector<unsigned> pos = {0, 1, 2};

    extractSeqInformation(counts, testLine, pos);

    BOOST_CHECK_MESSAGE(passSuppFilter(counts[0][0], 2) == false,
                        "passSuppFilter(counts[0][0], 2) = " << passSuppFilter(counts[0][0], 2));
    BOOST_CHECK_MESSAGE(passSuppFilter(counts[0][1], 1) == true,
                        "passSuppFilter(counts[0][1], 1) = " << passSuppFilter(counts[0][1], 1));
    BOOST_CHECK_MESSAGE(passSuppFilter(counts[0][2], 2) == true,
                        "passSuppFilter(counts[0][2], 2) = " << passSuppFilter(counts[0][2], 2));
}

BOOST_AUTO_TEST_CASE(read_data_pass_freq_filter) {
    std::vector<std::string> testLine = {"chr", "13", "REF", "6", "..ACGG", "IIIIII", "2", "AC", "II", "5",
                                         ".-3AACAAT.",
                                         "IIIII"};
    std::vector<std::array<unsigned, 5>> counts;
    counts.resize(3);

    std::vector<unsigned> pos = {0, 1, 2};

    extractSeqInformation(counts, testLine, pos);

    BOOST_CHECK_MESSAGE(passFreqFilter(counts[0][2], counts[0][4], 0.4) == false,
                        "passFreqFilter(counts[0][2], counts[0][4], 0.4) = "
                                << passFreqFilter(counts[0][2], counts[0][4], 0.4));
    BOOST_CHECK_MESSAGE(passFreqFilter(counts[0][2], counts[0][4], 0.33) == true,
                        "passFreqFilter(counts[0][2], counts[0][4], 0.4) = "
                                << passFreqFilter(counts[0][2], counts[0][4], 0.4));
}

BOOST_AUTO_TEST_CASE(read_data_apply_filter_across_cells) {
    Config<SampleTree> config;
    config.setNumSamples(3);
    config.minCoverage = 5;
    config.minFreq = 0.25;
    config.minSupport = 2;
    config.minNumCellsPassFilter = 2;

    std::vector<std::string> testLine = {"chr", "13", "REF", "4", "..AA", "IIII", "2", "AC", "II", "5", ".-3AACAAT.",
                                         "IIIII"};
    std::vector<std::array<unsigned, 5>> counts;
    counts.resize(3);
    std::vector<unsigned> pos = {0, 1, 2};

    extractSeqInformation(counts, testLine, pos);

    BOOST_CHECK_MESSAGE(applyFilterAcrossCells(counts, config, 0) == false,
                        "applyFilterAcrossCells(counts, config, 0) = " << applyFilterAcrossCells(counts, config, 0));

    testLine = {"chr", "13", "REF", "5", "...AA", "IIIII", "2", "AC", "II", "5", ".-3AACAAT.", "IIIII"};
    extractSeqInformation(counts, testLine, pos);

    BOOST_CHECK_MESSAGE(applyFilterAcrossCells(counts, config, 0) == true,
                        "applyFilterAcrossCells(counts, config, 0) = " << applyFilterAcrossCells(counts, config, 0));

    testLine = {"chr", "13", "REF", "10", ",,,,,....A", "IIIIIIIIII", "2", "AC", "II", "5", ".-3AACAAT.", "IIIII"};
    extractSeqInformation(counts, testLine, pos);

    BOOST_CHECK_MESSAGE(applyFilterAcrossCells(counts, config, 0) == false,
                        "applyFilterAcrossCells(counts, config, 0) = " << applyFilterAcrossCells(counts, config, 0));

    testLine = {"chr", "13", "REF", "5", "..A", "IIIII", "2", "AC", "II", "5", ".-3AACAAT.", "IIIII"};
    extractSeqInformation(counts, testLine, pos);

    BOOST_CHECK_MESSAGE(applyFilterAcrossCells(counts, config, 0) == false,
                        "applyFilterAcrossCells(counts, config, 0) = " << applyFilterAcrossCells(counts, config, 0));

}

BOOST_AUTO_TEST_CASE(logBetaBinPDF_test) {
    std::cout << logBetaBinPDF(200,400,0.01,100) << std::endl;
}

BOOST_AUTO_TEST_CASE(score_tree_big_tree) {
    {
        unsigned numSamples = 250;
        Config<SampleTree> config;
        config.setNumSamples(numSamples);
        config.noiseScore = -1000;
        config.data.resize(numSamples);
        config.data[0].resize(1);
        config._tmpAttachmentScore.resize(numSamples + numSamples - 1);
        config.setParam(Config<SampleTree>::E_nu, 0.1);
        config.setParam(Config<SampleTree>::E_lambdaWildLoss, 0.2);
        config.setParam(Config<SampleTree>::E_lambdaMutLoss, 0.25);
        config.setParam(Config<SampleTree>::E_parallel, 0.3);

        std::get<0>(config.logScores).resizeNumCells(numSamples);
        std::get<0>(config.logScores).resizeNumMuts(1);
        for (unsigned i = 0; i < numSamples; ++i)
        {
            std::get<0>(config.logScores).wtScore(i, 0) = std::log(0.1);
            std::get<0>(config.logScores).hetScore(i, 0) = std::log(0.5);
            std::get<0>(config.logScores).homScore(i, 0) = std::log(0.005);
        }

        config.learnZygocity = true;
        config.computeLossScore = true;
        config.computeParallelScore = true;
        createInitialTreeTest(config);
        std::cout << "1: " << scoreTree(config) << std::endl;
        
        Config<SampleTree>::TAttachmentScores & attachmentScores = config.getTmpAttachmentScore();
        Config<SampleTree>::TAttachmentScores::TAttachmentScore scoreSum;
        getAllAttachmentScores(attachmentScores, scoreSum, config, 0);
        scoreSum.computeFinalScore(config, scoreSum, false, true);

        static std::array<double, 6> res;
        res = scoreSum.computeMutTypeContribution(config, scoreSum, true);
        std::cout << "res: " << res[0] << " "  << res[1] << " " << res[2] << " " << res[3] << " " << res[4] << " " << res[5] << std::endl;
        
        //for (int i = 0; i < config._tmpAttachmentScore.size(); ++i)
        //    std::cout << i << ": " << config._tmpAttachmentScore[i] << std::endl;
        std::cout << "lossWildScore: ------" << std::endl;
        std::cout << "lossAltRScore: ------" << std::endl;
        std::cout << "lcaRScore: ------" << std::endl;
        std::cout << "1: " << scoreTree(config) << std::endl;
        std::string path = ".";
        writeTree(config, path);
    }
    {
        unsigned numSamples = 250;
        Config<SampleTree> config;
        config.setNumSamples(numSamples);
        config.noiseScore = -1000;
        config.data.resize(numSamples);
        config.data[0].resize(1);
        config._tmpAttachmentScore.resize(numSamples + numSamples - 1);
        config.setParam(Config<SampleTree>::E_nu, 0.1);
        config.setParam(Config<SampleTree>::E_lambdaWildLoss, 0.2);
        config.setParam(Config<SampleTree>::E_lambdaMutLoss, 0.25);
        config.setParam(Config<SampleTree>::E_parallel, 0.3);

        std::get<0>(config.logScores).resizeNumCells(numSamples);
        std::get<0>(config.logScores).resizeNumMuts(1);
        for (unsigned i = 0; i < numSamples; ++i)
        {
            std::get<0>(config.logScores).wtScore(i, 0) = std::log(0.005);
            std::get<0>(config.logScores).hetScore(i, 0) = std::log(0.5);
            std::get<0>(config.logScores).homScore(i, 0) = std::log(0.005);
        }

        config.learnZygocity = true;
        config.computeLossScore = true;
        config.computeParallelScore = true;
        createInitialTreeTest(config);
        std::cout << "1: " << scoreTree(config) << std::endl;
    }
    {
        unsigned numSamples = 250;
        Config<SampleTree> config;
        config.setNumSamples(numSamples);
        config.noiseScore = -1000;
        config.data.resize(numSamples);
        config.data[0].resize(1);
        config._tmpAttachmentScore.resize(numSamples + numSamples - 1);
        config.setParam(Config<SampleTree>::E_nu, 0.1);
        config.setParam(Config<SampleTree>::E_lambdaWildLoss, 0.05);
        config.setParam(Config<SampleTree>::E_lambdaMutLoss, 0.05);
        config.setParam(Config<SampleTree>::E_parallel, 0.05);

        std::get<0>(config.logScores).resizeNumCells(numSamples);
        std::get<0>(config.logScores).resizeNumMuts(1);
        for (unsigned i = 0; i < numSamples; ++i)
        {
            std::get<0>(config.logScores).wtScore(i, 0) = std::log(0.1);
            std::get<0>(config.logScores).hetScore(i, 0) = std::log(0.5);
            std::get<0>(config.logScores).homScore(i, 0) = std::log(0.005);
        }

        config.learnZygocity = true;
        config.computeLossScore = true;
        config.computeParallelScore = true;
        createInitialTreeTest(config);
        std::cout << "2: " << scoreTree(config) << std::endl;
    }
}
