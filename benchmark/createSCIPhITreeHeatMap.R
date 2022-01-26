# SCIPhI: Single-cell mutation identification via phylogenetic inference
#
# Copyright (C) 2018 ETH Zurich, Jochen Singer
#
# This file is part of SCIPhI.
#
# SCIPhIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SCIPhIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SCIPhI. If not, see <http://www.gnu.org/licenses/>.
# 
# @author: Jochen Singer

#! Rscript

library(ggplot2)
library(reshape)

args <- commandArgs(trailingOnly = TRUE)
inputName <- args[1]
outputName <- args[2]

header <- readLines(inputName, 1)
header <- unlist(strsplit(header, '\t'))[-1]
header <- gsub(header, pattern = '.bam', replacement = '')

df <- read.table(inputName, header = TRUE)
sciphiDist <- data.matrix(df)
meltData <- melt(t(sciphiDist[,-1]))

dict <- data.frame(label = header, lookup = names(df)[-1])
meltData$X1 <- dict$label[match(meltData$X1, dict$lookup)]


resultPlot <- ggplot(meltData, aes(X1,X2, fill=value)) +
  scale_fill_gradient(name = expression("P"[m]),
                    low = "white",
                    high = "steelblue",
                    limits=c(0,1)) +
  geom_raster() +
  scale_y_reverse() +
  ylab("Mutations") +
  xlab("Cells") +
  scale_x_discrete(limits=header, labels=header, breaks=header) +
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position="right",
        legend.direction="vertical",
        legend.box='horizontal',
        #plot.background = element_rect(fill = "transparent",colour = NA),
        #axis.text.x = element_text("",margin=margin(-20,0,0,0)),
        axis.title.x=element_text("",margin=margin(-3,0,0,0)),
        axis.text.x=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

ggsave(outputName)
