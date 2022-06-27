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
library(ggplot2)

args <- commandArgs(TRUE)

inputName <- args[1]
df <- read.table(inputName, header = TRUE)
df$loss <- as.factor(df$loss)
df$tool <- factor(df$tool, levels = c("Monovar", "SCIPhI", "SCIPhIN", "sciphi_max"))


ggplot(data = df, aes(x = loss, y = recall, fill = tool)) +
  geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour = tool), show.legend = FALSE) +
  geom_boxplot(outlier.size = NULL, outlier.shape = NA, alpha = 0.5) +
  xlab("Fraction of mutation losses") +
  ylab("Recall") +
  scale_y_continuous(limits = c(min(0.5, df$recall - 0.01), 1)) +

  theme(legend.title=element_blank(),
        legend.position = c(0.2, 0.2),
        text = element_text(size=25),
        legend.text.align = 0,
        legend.key.size = unit(3., 'lines')) +
  scale_color_manual(values = c("firebrick3", "lightsteelblue", "steelblue", "grey")) +
  scale_fill_manual(values = c("firebrick3", "lightsteelblue", "steelblue", "grey"), labels=c("Monovar", expression(paste("SCI", Phi)), expression(paste("SCI", Phi, "N")),expression(paste("SCI", Phi, "N_max"))))
ggsave(paste(gsub(".txt","",inputName), "_rec.pdf", sep=""))
  
ggplot(data = df, aes(x = loss, y = precision, fill = tool)) +
  geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour = tool), show.legend = FALSE) +
  geom_boxplot(outlier.size = NULL, outlier.shape = NA, alpha = 0.5) +
  xlab("Fraction of mutation losses") +
  ylab("Precision") +
  expand_limits(y=0.9) +
  theme(legend.title=element_blank(),
        legend.position = c(0.2, 0.2),
        text = element_text(size=25),
        legend.text.align = 0,
        legend.key.size = unit(3., 'lines')) +
  scale_color_manual(values = c("firebrick3", "lightsteelblue", "steelblue", "grey")) +
  scale_fill_manual(values = c("firebrick3", "lightsteelblue", "steelblue", "grey"), labels=c("Monovar", expression(paste("SCI", Phi)), expression(paste("SCI", Phi, "N")),expression(paste("SCI", Phi, "N_max"))))
ggsave(paste(gsub(".txt","",inputName), "_pre.pdf", sep=""))

ggplot(data = df, aes(x = loss, y = f1, fill = tool)) +
  geom_point(position = position_jitterdodge(jitter.width = 1), aes(colour = tool), show.legend = FALSE) +
  geom_boxplot(outlier.size = NULL, outlier.shape = NA, alpha = 0.5) +
  xlab("Fraction of mutation losses") +
  ylab("F1 score") +
  scale_y_continuous(limits = c(min(0.6, df$f1 - 0.01), 1)) +
  theme(legend.title=element_blank(),
        legend.position = c(0.2, 0.2),
        text = element_text(size=25),
        legend.text.align = 0,
        legend.key.size = unit(3., 'lines')) +
  scale_color_manual(values = c("firebrick3", "lightsteelblue", "steelblue", "grey")) +
  scale_fill_manual(values = c("firebrick3", "lightsteelblue", "steelblue", "grey"), labels=c("Monovar", expression(paste("SCI", Phi)), expression(paste("SCI", Phi, "N")),expression(paste("SCI", Phi, "N_max"))))
ggsave(paste(gsub(".txt","",inputName), "_f1.pdf", sep=""))
