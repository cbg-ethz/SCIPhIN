# SCIPhI: Single-cell mutation identification via phylogenetic inference
#
# Copyright (C) 2018 ETH Zurich, Jochen Singer
#
# This file is part of SCIPhI.
#
# SCIPhI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SCIPhI is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SCIPhI. If not, see <http://www.gnu.org/licenses/>.
# 
# @author: Jochen Singer
library(ggplot2)

args <- commandArgs(TRUE)

inputName <- args[1]
options(scipen=999)
df <- read.table(inputName, header = TRUE)
df$drop <- as.factor(df$drop)

ggplot(data = df, aes(x = drop, y = recall, fill = tool)) +
  geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour = tool),  show.legend = FALSE) +
  geom_boxplot(outlier.size = NULL, outlier.shape = NA, alpha = 0.5) +
  xlab("Drop-out rate") +
  ylab("Recall") +
  scale_y_continuous(limits = c(0.5, 1)) +
  theme(legend.position = c(0.8, 0.2),
        text = element_text(size=25),
        legend.title=element_blank(),
        legend.text.align = 0,
        legend.key.size = unit(3., 'lines')) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=c("Monovar", expression(paste("SCI", Phi)), "VarScan3")) +
  scale_colour_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=c("Monovar", expression(paste("SCI", Phi)), "VarScan3"))
ggsave(paste(gsub(".txt","",inputName), "_rec.pdf", sep=""))

ggplot(data = df, aes(x = drop, y = precision, fill = tool)) +
  geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour =tool),  show.legend = FALSE) +
  geom_boxplot(outlier.size = NULL, outlier.shape = NA, alpha = 0.5) +
  xlab("Drop-out rate") +
  ylab("Precision") +
  scale_y_continuous(limits = c(0.8, 1)) +
  theme(legend.position = c(0.8, 0.2),
        text = element_text(size=25),
        legend.title=element_blank(),
        legend.text.align = 0,
        legend.key.size = unit(3., 'lines')) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=c("Monovar", expression(paste("SCI", Phi)), "VarScan3")) +
  scale_colour_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=c("Monovar", expression(paste("SCI", Phi)), "VarScan3"))
ggsave(paste(gsub(".txt","",inputName), "_pre.pdf", sep=""))

ggplot(data = df, aes(x = drop, y = f1, fill = tool)) +
  geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour = tool),  show.legend = FALSE) +
  geom_boxplot(outlier.size = NULL, outlier.shape = NA, alpha = 0.5) +
  xlab("Drop-out rate") +
  ylab("F1 score") +
  scale_y_continuous(limits = c(0.6, 1)) +
  theme(legend.position = c(0.2, 0.2),
        text = element_text(size=25),
        legend.title=element_blank(),
        legend.text.align = 0,
        legend.key.size = unit(3., 'lines')) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=c("Monovar", expression(paste("SCI", Phi)), "VarScan3")) +
  scale_colour_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=c("Monovar", expression(paste("SCI", Phi)), "VarScan2"))
ggsave(paste(gsub(".txt","",inputName), "_f1.pdf", sep=""))
