library(dplyr)
library(gridExtra)
library(ggplot2)
library(optparse)

#######arguments
option_list <- list( 
  make_option(c("-i", "--input"), metavar="path", dest="ORFheaderdir",help="Specify the directory containing all the summary files generated",default=NULL)
)

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to draw the summary plots for contigs ORF prediction header file.
                               Usage: Rscript ~/pipeline/Github/CII_meta/ORF_header_summary.R -i ./"))
#######

setwd(opt$ORFheaderdir)
lengths_file <- "ORF_lengths.txt"
gc_conts_file <- "ORF_gc_conts.txt"
partials_file <- "ORF_partials.txt"
start_types_file <- "ORF_start_types.txt"


lengths_data <- read.table(lengths_file,  header = F , fill=TRUE, na.strings="", check.names=FALSE)
gc_conts_data <- read.table(gc_conts_file,  header = F , fill=TRUE, na.strings="", check.names=FALSE)
partials_data <- read.table(partials_file,  header = F , fill=TRUE, na.strings="", check.names=FALSE, colClasses=c("character"))
start_types_data <- read.table(start_types_file,  header = F , fill=TRUE, na.strings="", check.names=FALSE)

#sort(-lengths_data$V1)
head(lengths_data)
lengths_data2<-filter(lengths_data, V1 <= 2500)
bar_length <- ggplot(lengths_data2, aes(x=V1)) + 
  xlab("Length") + ylab("Density") + ggtitle("ORF Length Distribution") +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=100,
                 colour="black", fill="white") +
  theme(text = element_text(size = 12), axis.text.x = element_text(size = 12,vjust = 1,hjust = 1), plot.title = element_text(color="black", size=14, face="bold")) + 
  geom_density(alpha=.2, fill="#FF6666")

###
bar_gc <- ggplot(gc_conts_data, aes(x=V1)) + 
  xlab("GC content") + ylab("Density") + ggtitle("ORF GC Content Distribution") +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=0.05,
                 colour="black", fill="white") +
  theme(text = element_text(size = 12), axis.text.x = element_text(size = 12,vjust = 1,hjust = 1), plot.title = element_text(color="black", size=14, face="bold")) + 
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
#################

start_types_data2<-as.data.frame(table(start_types_data))
start_types_data3 <- start_types_data2 %>% mutate(Proportion = Freq/sum(Freq)*100)

bp_ST<- ggplot(start_types_data3, aes(x="", y=Proportion, fill=start_types_data))+ geom_bar(width = 1, colour="black", stat = "identity") + theme(panel.background = element_rect(fill = "white", colour = "white")) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ylab("Proportion (%)") + ggtitle("Pie Chart for ORF Start Types") +
  theme(text = element_text(size = 12), axis.text.x = element_text(size = 12,vjust = 1,hjust = 1), plot.title = element_text(color="black", size=14, face="bold"))
pie_ST <- bp_ST + coord_polar("y", start=0) + scale_fill_brewer(palette="Dark2")

###
partials_data2<-as.data.frame(table(partials_data))
partials_data3 <- partials_data2 %>% mutate(Proportion = Freq/sum(Freq)*100)

bp_P <- ggplot(partials_data3, aes(x="", y=Proportion, fill=partials_data))+ geom_bar(width = 1, colour="black", stat = "identity") + theme(panel.background = element_rect(fill = "white", colour = "white")) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ylab("Proportion (%)") + ggtitle("Pie Chart for ORF completeness") +
  theme(text = element_text(size = 12), axis.text.x = element_text(size = 12,vjust = 1,hjust = 1), plot.title = element_text(color="black", size=14, face="bold"))
pie_P <- bp_P + coord_polar("y", start=0)  + scale_fill_brewer(palette="Dark2")

merged<- grid.arrange(bar_length, bar_gc, pie_ST, pie_P, nrow = 2)
ggsave(plot = merged,filename = "ORF_summary.pdf", width = 12, height = 6,dpi = 300)
ggsave(plot = merged,filename = "ORF_summary.png", width = 12, height = 6,dpi = 300)

