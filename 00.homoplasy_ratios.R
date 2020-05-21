
############################################################
############################################################
####### ratio homoplmasic offsprings vs. non homo ##########
############################################################
############################################################

library(ape)
library(phangorn)
library(reshape2)

# folder_output <- "./dataset_test/"
folder_output <- "./dataset_3/" # ante Last dataset with NextStrain masking
folder_output <- "./dataset_4/" # ante Last dataset with DeMaio masking
folder_output <- "./dataset_5.1/" # last dataset with NextStrain masking
folder_output <- "./dataset_5.2/" #  Last dataset with DeMaio masking
folder_output <- "./dataset_6/" # SRA filter
folder_output <- "./dataset_7/" # SRA no low pos

setwd("/home/bacterio/Desktop/coronavirus_homoplasy_v2")
setwd(folder_output)

arbre_filtered_file <- list.files(path = ".", pattern = "^annotatedNewickTree_.*.tree_filtered")

# temp_tree <- read.tree(arbre_filtered_file)
# temp_tree$tip.label[temp_tree$tip.label!="NC_045512.2"] <- ""
# plot(temp_tree)

arbre_filtered <- read.tree(arbre_filtered_file)

# check format of node labels
# arbre_filtered$node.label[arbre_filtered$node.label != ""]

input <- read.csv("input_matrix", row.names=1,sep="\t",dec=",", stringsAsFactors=F,check.names=FALSE)

orphan_tips <- arbre_filtered$t[! arbre_filtered$t%in%colnames(input)]
# this object SHOULD be empty :
orphan_tips
# arbre_filtered <- drop.tip(arbre_filtered,orphan_tips)
# 
# jpeg(filename = "tree.jpg", width = 3500, height = 1600)
# plot(arbre_filtered,show.node.label=T,show.tip.label=F)
# dev.off()

#########################################################################
#########################################################################

# read the node labels and melt them to have a usefull format
annotation_nodes_useful_format <- melt(strsplit(arbre_filtered$node.label,"-"))
annotation_nodes_useful_format[,1] <- as.numeric(as.character(annotation_nodes_useful_format[,1])) # added 12 may

# get the architecture of the tree 
# this function gives, for each node and tip of the tree, the name of its descendants
architecture_arbre <- Descendants(arbre_filtered,type="all")

# FIRST CONDITION
# we don't want to study nodes that were annotated as carying an homoplasy AND have less than X descendants TIPS (equivalent to N nodes + tips)
# architecture_arbre[[312]] from 1 to 312 are tips
nb_offspring <- unlist(lapply(architecture_arbre, length))
names(nb_offspring) <- seq(from=1,to=length(nb_offspring),by=1)
min_offspring <- 3
enough_offsprings <- nb_offspring[nb_offspring>=min_offspring]

# STILL FIRST CONDITION
# i delete nodes that don't have enough offsprings from list of nodes to study homoplasies from
annotation_nodes_useful_format <- annotation_nodes_useful_format[(annotation_nodes_useful_format[,2]+length(arbre_filtered$tip.label))%in%names(enough_offsprings),]

# store all homoplasic SNP positions in a vector
homoplasies <- as.numeric(as.character(unique(annotation_nodes_useful_format[,1])))

# create object that will receive result
problems <- c()
raw_out_table <- matrix( nrow = 10000000, ncol = 3)
undef <- c()
i <- 0

# it's better to avoid for loops in R but here it improves readibility and we shouldn't study thousdands of homoplasies
for(homoplasy in homoplasies){  # homoplasy <- 1912
  corresponding_lines <- annotation_nodes_useful_format[annotation_nodes_useful_format[,1]==homoplasy,2]

  for(line in corresponding_lines){   # line <- 229
    node <- line+length(arbre_filtered$tip.label)
    tips <- architecture_arbre[[node]][architecture_arbre[[node]]<=(length(arbre_filtered$tip.label))]
    child_nodes <- architecture_arbre[[node]][architecture_arbre[[node]]>(length(arbre_filtered$tip.label))]
    # SECOND CONDITION 
    # we don't want the studied node to have any children node displaying the same homoplasy
    # so if we match a line of annotation_nodes_useful_format with the same homoplasy position on a tip descendant of the one we're studying, we discard it
    boolean_only_false_to_continue <- unlist(lapply(child_nodes, function(x) { annotation_nodes_useful_format[,1]==homoplasy & annotation_nodes_useful_format[,2]==x  }))
    
    if( sum(boolean_only_false_to_continue) == 0 ){
    
    NOT_homoplasy_result <- 0
    homoplasy_result <- 0
    
    alleles <- sapply(tips,function(x) input[as.character(homoplasy),arbre_filtered$t[x]])
    NOT_homoplasy_result <- sum(alleles=="ref")
    homoplasy_result <- sum(alleles=="not_ref")
    undef <- c(undef,sum(alleles=="undef"))

    # for a considered homoplasic node, it's weird not to have both "ref" and "not_ref" offsprings. If it's not the case, report positions:
    if(NOT_homoplasy_result == 0 | homoplasy_result == 0 ){
        problems <- c(problems,homoplasy);  print("problem")
    }else{
      i <- i+1
      raw_out_table[i,] <- c(homoplasy,homoplasy_result,NOT_homoplasy_result)
    }
    }
}
}

# Should be empty :
problems

object_backup <- raw_out_table
# raw_out_table <- object_backup

library(ggplot2)
nrow(raw_out_table)
colnames(raw_out_table) <- c("position","homoplasy_count","NOT_homoplasy_count")
raw_out_table <- as.data.frame(raw_out_table[!is.na(raw_out_table[,1]),])
###########################

# plot of all values of homoplasy allele count vs non-homoplasy allele count
out_table_long <- melt(raw_out_table, id=c("position"), measured=c("homoplasy_count","NOT_homoplasy_count"))

effectif <- nrow(raw_out_table)
jpeg(filename = paste("fig_1.non_filtered_count_histogram_min_offspring_",min_offspring,".jpg",sep=""), width = 1000, height = 700)
ggplot(out_table_long, mapping = aes (fill = variable, x = log(value))) + geom_density (alpha = .5) + ggtitle(paste("Homoplasy allele count vs non-homoplasy allele count\n non filtered dataset, n= ",effectif,sep=""))
dev.off()

out_table_ratios <- as.data.frame(cbind("position"=raw_out_table[,1],"homo/not_homo"=raw_out_table[,2]/raw_out_table[,3],"not_homo/homo"=raw_out_table[,3]/raw_out_table[,2]))
out_table_ratios_long <- melt(out_table_ratios, id=c("position"), measured=c("homo/not_homo","not_homo/homo"))

effectif <- sum(out_table_ratios_long[,2]=="homo/not_homo")
nb_unique_nodes <- length(unique(out_table_ratios_long[,1]))

jpeg(filename = paste("fig_2.homo_vs_nothomo_ratio_along_genome_min_offspring_",min_offspring,".jpg",sep=""), width = 1000, height = 700)
ggplot(out_table_ratios_long[out_table_ratios_long[,2]=="homo/not_homo",],
       aes(x=position, y=log(value), group=position))+
  geom_boxplot(outlier.shape = NA)+
  stat_summary(geom="point", fun.y=median)+
  geom_hline(yintercept = median(out_table_ratios_long[out_table_ratios_long[,2]=="homo/not_homo","value"]), linetype='dashed')+
  ggtitle(paste("log of the ratio Homoplasy allele count / non-homoplasy allele count\nnb of studied homoplasy positions = ",nb_unique_nodes,", nb of studied paired count values= ",effectif,"\ndashed line = mediane",sep=""))
dev.off()

jpeg(filename = paste("fig_3.nothomo_vs_homo_ratio_along_genome_min_offspring_",min_offspring,".jpg",sep=""), width = 1000, height = 700)
ggplot(out_table_ratios_long[out_table_ratios_long[,2]=="not_homo/homo",],
       aes(x=position, y=log(value), group=position))+
  geom_boxplot(outlier.shape = NA)+
  stat_summary(geom="point", fun.y=median)+
  geom_hline(yintercept = median(out_table_ratios_long[out_table_ratios_long[,2]=="not_homo/homo","value"]), linetype='dashed')+
  ggtitle(paste("log of the ratio non-Homoplasy allele count / homoplasy allele count\nnb of studied homoplasy positions = ",nb_unique_nodes,", nb of studied paired count values= ",effectif,"\ndashed line = mediane",sep=""))
dev.off()

jpeg(filename = paste("fig_4.non_filtered_ratio_histogram_min_offspring_",min_offspring,".jpg",sep=""), width = 1000, height = 700)
ggplot(out_table_ratios_long[out_table_ratios_long[,2]=="homo/not_homo",], mapping = aes (x = log(value))) + geom_density (alpha = .5) + ggtitle(paste("log of the ratio Homoplasy allele count / non-homoplasy allele count\n non filtered dataset, n= ",effectif,sep=""))
dev.off()


#####################################################################################
############ same plots with filtering of the data : nb of replicates ###############
############ eg. 5 nodes in the phylogeny detected as origin of a new homoplasy #####
#####################################################################################

# nb_rep <- 3
nb_rep <- 5
nb_rep <- 10

# i will calculate limits for the plots with nb_rep = 5 in order to get the exact same plots for both filters
out_table_restricted_5 <- raw_out_table[raw_out_table[,1]%in%names(table(raw_out_table[,1])[table(raw_out_table[,1])>=5]),]
ylimplot_min <- min(log(out_table_restricted_5[,2]/out_table_restricted_5[,3]))
ylimplot_max <- max(log(out_table_restricted_5[,2]/out_table_restricted_5[,3]))
ylimplot_min <- min(log10(out_table_restricted_5[,2]/out_table_restricted_5[,3]))
ylimplot_max <- max(log10(out_table_restricted_5[,2]/out_table_restricted_5[,3]))
ylimplot_min  <- c(-4)
ylimplot_max <- c(4)
# end of plot limit calculations

names(table(raw_out_table[,1])[table(raw_out_table[,1])>=nb_rep])

out_table_restricted <- raw_out_table[raw_out_table[,1]%in%names(table(raw_out_table[,1])[table(raw_out_table[,1])>=nb_rep]),]
str(out_table_restricted)

######################################################################################
############################# STATISTICAL ANALYSIS ###################################
######################################################################################
# Multiple test are perfoemed :
# paired t test #
# Shapiro test to check wether the differences follow a normal distrib
# Sign test that do not need the diff to follow a norm. distrib.
# it only allows to say if the distrib isequal or not to zero.
# The sign test is indifferent to magnitude of difference. It only cares if it is above or below M0. If you wish to incorporate the magnitude, you should look at the Signed Rank test (wilcox.test).

library(BSDA)

multiple_t_test <- sapply(sort(unique(out_table_restricted[,1])), function(x) {
  res <- t.test(out_table_restricted[out_table_restricted[,1] == x,2],out_table_restricted[out_table_restricted[,1] == x,3], paired=TRUE)
  # if p-value < 0.05, we reject equality hypothesis, so two distributions are not equal.
  # res$estimate is the mean value of differences between members of the pair
  # if negative, it means that x < y # in our case if negative, it means that homo < not_homo
  shapiro <- shapiro.test(out_table_restricted[out_table_restricted[,1] == x,2]-out_table_restricted[out_table_restricted[,1] == x,3])$p.value
  # if p-value < 0.05, distribution do NOT follow a normal distribution
  sign_test <- SIGN.test(out_table_restricted[out_table_restricted[,1] == x,2]-out_table_restricted[out_table_restricted[,1] == x,3], y = NULL, md = 0, alternative = "two.sided", conf.level = 0.95)$p.value
  sign_test_bonferroni <- p.adjust(sign_test, method = "bonferroni", n=length(unique(out_table_restricted[,1])))
  wilcox_test <- wilcox.test(out_table_restricted[out_table_restricted[,1] == x,2]-out_table_restricted[out_table_restricted[,1] == x,3], paired = F, mu=0, alternative = "two.sided")$p.value
  nb_replicates <- length(out_table_restricted[out_table_restricted[,1] == x,2])
  # null hypothesis is "median = 0". We do NOT reject NULL hypothesis if p-value > 0.05
  # so no significantly different from zero if p-value > 0.05
  return(c(position=x,t_test=res$p.value,mean_of_diff=res$estimate,shapiro=shapiro,sign_test=sign_test,sign_test_bonferroni=sign_test_bonferroni,wilcox_test=wilcox_test,nb_replicates=nb_replicates))
})


write.table(t(multiple_t_test), file = paste("T_test_min_",nb_rep,"_replicates_min_offspring_",min_offspring,".tsv",sep=""), append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"), fileEncoding = "")
# end paired t test

out_table_ratios_restricted_long <- as.data.frame(cbind(position=out_table_restricted[,1],"homo/not_homo"=out_table_restricted[,2]/out_table_restricted[,3]))

write.table(cbind(out_table_ratios_restricted_long,dataset=rep(folder_output,nrow(out_table_ratios_restricted_long))), file = paste("Homoplasy_RAW_pos_min_",nb_rep,"_replicates_min_offspring_",min_offspring,".tsv",sep=""), append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"), fileEncoding = "")

custom_x_titles <- c(0,sort(unique(out_table_ratios_restricted_long[,1])),29903)

ratios_mean_by_position <- aggregate(out_table_ratios_restricted_long[,2], list(Nucl_pos = out_table_ratios_restricted_long[,1]),mean)
ratios_median_by_position <- aggregate(out_table_ratios_restricted_long[,2], list(Nucl_pos = out_table_ratios_restricted_long[,1]),median)
ratios_sd_by_position <- aggregate(out_table_ratios_restricted_long[,2], list(Nucl_pos = out_table_ratios_restricted_long[,1]),sd)
ratios_mean_by_position <- cbind(ratios_mean_by_position,ratios_median_by_position[,2],ratios_sd_by_position[,2])
colnames(ratios_mean_by_position) <- c("Position","Mean_RoHo","Median_RoHo","Standard_deviation_RoHo")

write.table(ratios_mean_by_position, file = paste("Homoplasy_pos_min_",nb_rep,"_replicates_min_offspring_",min_offspring,".tsv",sep=""), append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"), fileEncoding = "")

table_effectifs <- table(out_table_ratios_restricted_long[,1])
table_effectifs <- as.data.frame(t(rbind(pos=as.numeric(names(table_effectifs)),effectif=as.numeric(unname(table_effectifs)))))

# jpeg(filename = paste("fig_1.homo_vs_nothomo_ratio_along_genome_min_",nb_rep,"_replicates_",filter,".jpg",sep=""), width = 1000, height = 500)
# svg(filename = paste("fig_5.RoHo_along_genome_min_",nb_rep,"_replicates_min_offspring_",min_offspring,".svg",sep=""), width = 10, height = 3.5)
svg(filename = paste("fig_5.RoHo_along_genome_min_",nb_rep,"_replicates_min_offspring_",min_offspring,".svg",sep=""), width = 10, height = 5)
ggplot(out_table_ratios_restricted_long, aes(x=position, y=log10(out_table_ratios_restricted_long[,2]), group=position))+
  coord_cartesian(xlim = c(0, 29903), ylim = c(ylimplot_min,ylimplot_max)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_vline(xintercept=table_effectifs$pos, linetype="dotted", color = "grey")+ # only for nrep=10
  # geom_boxplot(outlier.shape = NA)+
  # geom_boxplot( fill='#A4A4A4', color="darkred",width=200)+
  #geom_violin(trim=T, fill='#A4A4A4', color="darkred",width=400,position="dodge")+
  geom_violin(trim=T, fill='#A4A4A4', color="darkred",width=400,position=position_dodge(width = 90))+
  stat_summary(geom="point", fun.y=median)+
  scale_y_continuous(breaks=seq(from=trunc(ylimplot_min-1),to=trunc(ylimplot_max+1),by=1))+
  scale_x_continuous(breaks=custom_x_titles)+
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  # geom_text(data = homo_not_homo, aes(x=position, y=min(log(value)), group=position, label=position), size = 3, vjust = 0,angle = 90) +
stat_summary(geom = 'text', label = paste("n=",table_effectifs$effectif,sep=""), fun.y = function(x) max(x) + 0.5, vjust = 0,hjust = 0,angle = 45) +  # only for nrep=10
stat_summary(fun.y=mean, geom="point", shape=1, size=3, color="blue")

dev.off()

ylimplot_min  <- c(-3.5)
ylimplot_max <- c(3)

# if replicates = 5 :
svg(filename = paste("fig_5.RoHo_along_genome_min_",nb_rep,"_replicates_min_offspring_",min_offspring,".svg",sep=""), width = 10, height = 3.5)
ggplot(out_table_ratios_restricted_long, aes(x=position, y=log10(out_table_ratios_restricted_long[,2]), group=position))+
  coord_cartesian(xlim = c(0, 29903), ylim = c(ylimplot_min,ylimplot_max)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_violin(trim=T, fill='#A4A4A4', color="darkred")+
  stat_summary(geom="point", fun.y=median)+
  scale_y_continuous(breaks=seq(from=trunc(ylimplot_min-1),to=trunc(ylimplot_max+1),by=1))+
  scale_x_continuous(breaks=custom_x_titles)+
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y=mean, geom="point", shape=1, size=3, color="blue")

dev.off()

# here i can add a statistic comparison between the SRA RoHo index and the gisaid one
out_table_ratios_restricted_long

current_dataset_long <- cbind(out_table_ratios_restricted_long,dataset=rep(folder_output,nrow(out_table_ratios_restricted_long)))

dataset_6 <- read.csv("/home/bacterio/Desktop/coronavirus_homoplasy_v2/dataset_6/Homoplasy_RAW_pos_min_3_replicates_min_offspring_3.tsv", header=T,sep="\t",dec=",", stringsAsFactors=F,check.names=FALSE)
str(dataset_6)

dataset_6_matching <- dataset_6[dataset_6[,1]%in%current_dataset_long[,1],]
current_dataset_long_matching <- current_dataset_long[current_dataset_long[,1]%in%dataset_6[,1],]

two_experiment_dataset <- rbind(dataset_6_matching,current_dataset_long_matching)
two_experiment_dataset

two_experiment_dataset[,1] == 21137 & two_experiment_dataset[,3] == "./dataset_3/"

multiple_t_test <- sapply(sort(unique(two_experiment_dataset[,1])), function(x) {
  res <- t.test(as.numeric(two_experiment_dataset[two_experiment_dataset[,1] == x & two_experiment_dataset[,3] == "./dataset_3/",2]), as.numeric(two_experiment_dataset[two_experiment_dataset[,1] == x & two_experiment_dataset[,3] == "./dataset_6/",2]), paired=F)
  
  return(c(position=x,t_test=res$p.value,mean_of_diff=res$estimate))
})


# 
# jpeg(filename =  paste("fig_2.nothomo_vs_homo_ratio_along_genome_min_",nb_rep,"_replicates_",filter,".jpg",sep=""), width = 1000, height = 700)
# ggplot(out_table_ratios_restricted_long[out_table_ratios_restricted_long[,2]=="not_homo/homo",],
#        aes(x=position, y=log(value), group=position))+
#   geom_boxplot(outlier.shape = NA)+
#   stat_summary(geom="point", fun.y=median)+
#   ggtitle(paste("log of the ratio non-Homoplasy allele count / homoplasy allele count\nnb of studied homoplasy positions = ",nb_unique_nodes,", nb of studied paired count values= ",effectif,"\nmin replicates per homoplasy position = ",nb_rep,sep=""))+
#   scale_x_continuous(breaks=custom_x_titles)+
#   theme(axis.text.x = element_text(angle = 90))
# dev.off()
# 




