library(optparse)

option_list <- list(
    make_option(c("-i","--input"), help="table for scatter plot"),
    make_option(c("-x", "--axisX")),
    make_option(c("-X", "--labelX"), default = NULL),
    make_option(c("-y", "--axisY")),
    make_option(c("-Y", "--labelY"), default = NULL),
    make_option(c("-c", "--category")),
    make_option(c("-o", "--outpath"), default="./scatter.pdf")
)

opt <- parse_args(OptionParser(option_list=option_list,description = "bar, stacked bar, pie plot of genes copy table"))


library(ggplot2)

input <-read.table(opt$input, quote="", na.strings="",
 comment.char="", check.names=F, stringsAsFactors=F, header = TRUE, sep = "\t")

p <- ggplot(input, aes_string(x = opt$axisX, y = opt$axisY, color=opt$category)) +
 geom_point(size = 0.5) +
 theme_classic() +
 guides(color=guide_legend(title = NULL)) +
 theme(panel.background=element_rect(color='black', fill='transparent'))

if(!is.null(opt$labelX)){
    p <- p + xlab(opt$labelX)
}

if(!is.null(opt$labelY)){
    p <- p + ylab(opt$labelY)
}

ggsave(plot = p, filename = opt$outpath, width = 12, height = 10)
