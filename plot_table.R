#!/usr/bin/env Rscript
library(optparse)

#######arguments
option_list <- list(
    make_option(c("-i", "--input"), help="Specify the path of gene copy table",default=NULL),
    make_option(c("-T", "--transpose"), default=FALSE, help="Transpose the input table"),
    make_option(c("-t", "--top"), help="Keep only the most abundant rows by compare the row sums.", default=NULL),
    make_option(c("-O", "--other"), default=FALSE, help="Calculate the sum of gene count except top."),
    make_option(c("-c", "--colname"), help="Specify colname in gene copy table for single column pie plot and barplot. If not specified, then row sums will be used",default=NULL),
    make_option(c("-b", "--bar"), help="output single column bar plot",default=NULL),
    make_option(c("-P", "--pie"), help="output single column pie plot",default=NULL),
    make_option(c("-l", "--percent"), help="Calculate percent for the pie labels",default=FALSE),
    make_option(c("-a", "--annotate"), help="Wheather to annotate the genes in plot",default=TRUE),
    make_option(c("-B", "--boxplot"), help="output boxplot",default=NULL),
    make_option(c("-H", "--heatmap"), help="output heatmap",default=NULL),
    make_option(c("-s", "--stackbar"), help="output stacked bar",default=NULL),
    make_option(c("-p", "--prefix"), help="The prefix of output files",default="")
)

opt <- parse_args(OptionParser(option_list=option_list,description = "bar, stacked bar, pie plot of genes copy table"))

library(stringr)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
df_map<-function(df, func){
  out<-c()
  for(i in 1:ncol(df)){
    out[i]<-func(df[, i])
  }
  return(out)
}

choose<-function(condition,choice1,choice2){
  if(condition){
    return(choice1)
  }else{
    return(choice2)
  }
}


pallet<-c(brewer.pal(8,"Set2")[-c(7,8)], rev(brewer.pal(12,"Paired")), brewer.pal(8,"Dark2"),brewer.pal(12,"Set3"),brewer.pal(8,"Accent"),brewer.pal(11,"Spectral"))

get_pheatmap_color <- function(annotation){
  ua <- sort(unique(as.character(annotation)))
  color_map <- c()
  color_map[ua] <- pallet[1:length(ua)]
  return(list(levels=color_map))
}


optl <- list()

for(idx in c("bar", "pie", "boxplot", "heatmap", "stackbar")){
  drt <- opt[[idx]]
  not_null <- !is.null(drt)
  optl[[idx]] <- not_null
  if(not_null){
    if(!dir.exists(drt)){dir.create(drt,recursive = T)}
    if(!endsWith(drt, '/')){drt<-opt[[idx]]<-paste(drt, '/', sep = "")}
    opt[[idx]]<-paste(opt[[idx]], opt$prefix, sep = "")
    TMP_DIR <- drt
  }
}



input <-read.table(opt$input, row.names = 1, quote="", na.strings="",
 comment.char="", check.names=F, stringsAsFactors=F, header = TRUE, sep = "\t")

# input <- na.omit(input)

gene_levels <- NULL

nc = ncol(input)
if(!is.numeric(input[, nc])){
    gene_levels <- input[, nc]
    gene_levels <- sapply(gene_levels, function(x){return(str_split(x, '; *')[[1]][1])})
    names(gene_levels) <- rownames(input)
    gene_levels[is.na(gene_levels)] <- "not ranked"
    gene_levels['Other'] <- 'Other'
    # unique_gene_levels <- unique(gene_levels)
}
input <- input[, df_map(input, is.numeric)]

global_annotation <- opt$annotate && !(is.null(gene_levels)||opt$transpose)

if(opt$transpose){
    input <- data.frame(t(input), check.names = F)
}

if(optl$bar||optl$pie){
    # print(input[opt$colname])
    if(is.null(opt$colname)){
        input_single <- rowSums(input)
    }else{
        input_single <- input[opt$colname][,1]
    }
    names(input_single) <- rownames(input)
    input_single <- input_single[input_single>0]
    input_single <- sort(input_single, decreasing=TRUE)
    if(!is.null(opt$top)){
        top<-as.numeric(opt$top)
        nr <- length(input_single)
        if(top<nr-1){
            single_data <- input_single[1:top]
            if(opt$other){
                single_data["Other"] <- sum(input_single[top+1:nr], na.rm=TRUE)
            }
        }else{
            single_data <- input_single
        }
    }else{
        single_data <- input_single
    }
    if(optl$pie){
        pdf(paste(opt$pie,"pie.pdf", sep = ""), width=18, height=10)
        if(opt$percent){
            percent <- round(single_data/sum(single_data)*100, 1)
            label <- paste(names(single_data), "(", percent, "% )")
        }else{
            label <- paste(names(single_data), "(", single_data, ")")
        }
        pie(single_data, main=choose(is.null(opt$colname),"",opt$colname), label=label, border="white", col=pallet)
        dev.off()
    }
    if(optl$bar){
        bar_data <- data.frame(genes=names(single_data), counts=single_data, check.names = F, stringsAsFactors=FALSE)
        b_wd <- 7
        b_ht <- nrow(bar_data) * 0.2 + 1.5
        if(global_annotation){
            bar_data$levels <- gene_levels[names(single_data)]
            bar_data <- bar_data[order(bar_data$levels), ]
            uni_levels <- unique(bar_data$levels)
            b_wd <- b_wd + (max(nchar(uni_levels))*0.05+0.3)*ceiling(length(uni_levels)/17)
        }
        b_ht<-ifelse(b_ht<50,b_ht,49.9)
        p<-ggplot(bar_data,aes(x=genes,y=counts))+
          choose(global_annotation, geom_bar(mapping = aes(fill=levels), stat = "identity",width = 0.7),
            geom_bar(stat = "identity", width = 0.7, fill=pallet[1])
          )+
          guides(fill=guide_legend(title = NULL))+
          scale_fill_manual(values = pallet)+
          scale_x_discrete(limits=bar_data$genes)+
          xlab("")+ylab("No. of genes")+ theme_classic() +
          scale_y_continuous(expand = c(0, 0)) + coord_flip()
        if(!is.null(opt$colname)){
            p = p+labs(title=opt$colname) + theme(plot.title = element_text(hjust = 0.5))
        }
        ggsave(filename = paste(opt$bar, "bar.pdf", sep = ""), plot = p , height = b_ht, width = b_wd)
    }
}

if(optl$heatmap||optl$stackbar||optl$boxplot){
    input_multi <- t(input)
    input_multi <- input_multi[, order(colSums(input_multi), decreasing = T)]

    if(!is.null(opt$top)){
        top<-as.numeric(opt$top)
        if(top<ncol(input_multi)-1){
          plot_data<-data.frame(input_multi[,1:top], check.names = F, check.rows = T)
          if(opt$other){
            plot_data$Other=apply(input_multi[,(top+1):ncol(input_multi)], 1, sum)
          }
        }else{
          plot_data<-data.frame(input_multi,check.names = F,check.rows = T)
        }
    }else{
        plot_data<-data.frame(input_multi, check.names = F, check.rows = T)
    }
    plot_data_save <- plot_data

    if(optl$heatmap){
        write.table(plot_data_save, paste(opt$heatmap,"heatmap.xls", sep = ""), sep = "\t", quote=FALSE, col.names=NA)
        annotation_col <- NA
        annotation_row <- NA
        annotation_color<-NA
        if(!is.null(gene_levels)&&opt$annotate){
            if(opt$transpose){
                annotation_row <- data.frame(levels=gene_levels[rownames(plot_data)])
                rownames(annotation_row) <- rownames(plot_data)
                annotation_color<-get_pheatmap_color(annotation_row[,1])
            }else{
                annotation_col <- data.frame(levels=gene_levels[colnames(plot_data)])
                rownames(annotation_col) <- colnames(plot_data)
                annotation_color<-get_pheatmap_color(annotation_col[,1])
            }
        }
        # print(annotation_color)
        pheatmap(plot_data, file=paste(opt$heatmap,"heatmap.pdf", sep = ""),
             annotation_row=annotation_row,
             annotation_col=annotation_col,
             # fontsize=15, 字号太大，会导致图出界
             border_color = "grey",
             # color = colorRampPalette(colors = cellcolors)(100),
             cluster_cols=TRUE,clustering_distance_cols="euclidean",
             cellwidth=20, cellheight=20,
             # cutree_cols=3,
             annotation_colors=annotation_color,
             cluster_rows=TRUE,clustering_distance_rows="euclidean"
        )
    }

    if(optl$stackbar||optl$boxplot){
        plot_data$id<-rownames(plot_data)
        plot_data<-plot_data[rev(colnames(plot_data))]
        plot_data<-melt(plot_data,id.vars = "id")
        if(optl$stackbar){
            p1<-(max(nchar(colnames(plot_data_save)))*0.05+0.3)*ceiling(ncol(plot_data_save)/17)+2.5
            bar_wd<-nrow(plot_data_save)*0.2+p1
            bar_wd<-ifelse(bar_wd<50,bar_wd,49.9)
            write.table(plot_data_save, paste(opt$stackbar,"stackbar.xls", sep = ""), sep = "\t", quote=FALSE, col.names=NA)
            p<-ggplot(plot_data,aes(x=id,y=value))+geom_bar(mapping = aes(fill=variable), stat = "identity",width = 0.7)+
              guides(fill=guide_legend(title = NULL))+
              scale_fill_manual(values = pallet)+
              # scale_x_discrete(limits=label_order)+
              xlab("")+ylab("No. of genes")+theme_classic()+
              theme(text = element_text(size = 10),
                    axis.text.x = element_text(angle = 90,size = 10,hjust = 1, vjust = 0.5)) +
              scale_y_continuous(
                # limits=c(0,101),
                expand = c(0, 0))
            ggsave(filename = paste(opt$stackbar, "stackbar.pdf", sep = ""), plot = p , height = 7, width = bar_wd)
        }
        if(optl$boxplot){
            boxplot_wd <- ncol(plot_data_save)*0.4 + 2.5
            if(global_annotation){
                plot_data$levels <- gene_levels[as.character(plot_data$variable)]
                plot_data <- plot_data[order(plot_data$levels), ]
                uni_levels <- unique(plot_data$levels)
                boxplot_wd <- boxplot_wd + (max(nchar(uni_levels))*0.05+0.3)*ceiling(length(uni_levels)/17)
            }
            boxplot_wd<-ifelse(boxplot_wd<50,boxplot_wd,49.9)
            write.table(plot_data_save, paste(opt$boxplot,"boxplot.xls", sep = ""), sep = "\t", quote=FALSE, col.names=NA)
            p<-ggplot(plot_data,aes(x=variable,y=value))+
              choose(global_annotation, geom_boxplot(mapping = aes(fill=levels)), geom_boxplot(fill=pallet[1]))+
              geom_jitter(position=position_jitter(0.17), size=1, alpha=0.3)+
              guides(fill=guide_legend(title = NULL))+
              scale_fill_manual(values = pallet)+
              choose(global_annotation, scale_x_discrete(limits=unique(plot_data$variable)), theme())+
              xlab("")+ylab("No. of genes")+theme_classic()+
              theme(text = element_text(size = 10),
                    axis.text.x = element_text(angle = 90,size = 10,hjust = 1, vjust = 0.5))
            ggsave(filename = paste(opt$boxplot, "boxplot.pdf", sep = ""), plot = p , height = 7, width = boxplot_wd)
        }
    }
}

# dev.cur()