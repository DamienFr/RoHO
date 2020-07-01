
############################################################
############################################################
######### Ratio of Homoplasic Offspring (RoHO) #############
############################################################
############################################################

library(ape)
library(phangorn)
library(reshape2)
library(doParallel) # for parallelization
library(foreach) # for parallelization
library(ggplot2)

thread_number <- "auto" # set it to auto or to a number
# thread_number <- 6
deduplicated_dataset <- 0 # bolean, set to 1 (true) if dataset contains deduplicated sequences and if you have a table of re-duplication table
simulated_t_a_dataset <- 0 # set to 1 to run on simulated datasets, debug purposes # not compatible with deduplicated_dataset = 1
children_node_rule <- 1 # set to 1 to activate the filter that makes not study the nodes that have embedded homoplasic nodes 
min_offspring <- 4 # node annotated as carying an homoplasy and having less tips than "min_offspring" will be un-annotated
# they will not have their roho score computed and they will not count as embedded homoplasies neither.
tips_rule  <- 2

args <- commandArgs(TRUE)
folder_output <- paste("./",args[1],"/",sep="")

arbre_filtered_file <- list.files(path = ".", pattern = "^annotatedNewickTree_.*.tree_filtered")
arbre_filtered <- read.tree(arbre_filtered_file)

if(deduplicated_dataset){
  duplicated_indivs <- read.csv("./with_duplicated_isolates/iqtree_log_usefull3",sep="\t",dec=",", stringsAsFactors=F,check.names=FALSE)
}

input <- read.csv("input_matrix", row.names=1,sep="\t",dec=",", stringsAsFactors=F,check.names=FALSE)

if(deduplicated_dataset){
  # rapid test 
  duplicated_indivs[! duplicated_indivs[,1]%in%colnames(input),]
}
#########################################################################
#########################################################################

# read the node labels and melt them to have a usefull format
annotation_nodes_useful_format <- melt(strsplit(arbre_filtered$node.label,"-"))
annotation_nodes_useful_format[,1] <- as.numeric(as.character(annotation_nodes_useful_format[,1])) # added 12 may

# get the architecture of the tree 
# this function gives, for each node and tip of the tree, the name of its descendants
architecture_arbre <- Descendants(arbre_filtered,type="all")

# FIRST CONDITION
# we don't want to study nodes that were annotated as carying an homoplasy AND have less than X descendants TIPS
# architecture_arbre[[312]] from 1 to 312 are tips
nb_offspring <- unlist(lapply(architecture_arbre, length))
object <- table(nb_offspring)
names(nb_offspring) <- seq(from=1,to=length(nb_offspring),by=1)
min_offspring_tips_and_nodes <- (min_offspring * 2 )-2
enough_offsprings <- nb_offspring[nb_offspring>=min_offspring_tips_and_nodes]

# STILL FIRST CONDITION
# i delete nodes that don't have enough offsprings from list of nodes to study homoplasies from
annotation_nodes_useful_format <- annotation_nodes_useful_format[(annotation_nodes_useful_format[,2]+length(arbre_filtered$tip.label))%in%names(enough_offsprings),]

# store all homoplasic SNP positions in a vector
homoplasies <- as.numeric(as.character(unique(annotation_nodes_useful_format[,1])))

# for non-paralelized version, create object that will receive result
# i <- 0 ; raw_out_table <- matrix( nrow = 10000000, ncol = 3)

if (thread_number == "auto"){
  cores=detectCores()
  cl <- makeCluster(cores[1]-2) #not to overload your computer
}else{
  cl <- makeCluster(thread_number)
}

registerDoParallel(cl)

# for non-paralelized version, use  #for(homoplasy in homoplasies){ 
raw_out_table <- foreach(homoplasy=homoplasies, .combine=rbind) %dopar% {
  corresponding_lines <- annotation_nodes_useful_format[annotation_nodes_useful_format[,1]==homoplasy,2]
  temp_raw_out_table <- matrix( nrow = 0, ncol = 3)
  
  for(line in corresponding_lines){
    node <- line+length(arbre_filtered$tip.label)
    tips <- architecture_arbre[[node]][architecture_arbre[[node]]<=(length(arbre_filtered$tip.label))]
    
    child_nodes <- architecture_arbre[[node]][architecture_arbre[[node]]>(length(arbre_filtered$tip.label))]
    child_nodes <- child_nodes - length(arbre_filtered$tip.label)

    annotation_nodes_useful_format_considered_homoplasy <- annotation_nodes_useful_format[annotation_nodes_useful_format[,1]==homoplasy,]
    annotation_nodes_useful_format_considered_homoplasy_cons_child_nodes <- annotation_nodes_useful_format_considered_homoplasy[annotation_nodes_useful_format_considered_homoplasy[,2]%in%child_nodes,2]
    
    if( children_node_rule){
      count_children_node_with_homoplasy <- unlist(lapply(annotation_nodes_useful_format_considered_homoplasy_cons_child_nodes,function(x){
        tips <- architecture_arbre[[x+length(arbre_filtered$tip.label)]][architecture_arbre[[x+length(arbre_filtered$tip.label)]]<=(length(arbre_filtered$tip.label))]
        alleles <- sapply(tips,function(y) input[as.character(homoplasy),arbre_filtered$t[y]])
        
        if(deduplicated_dataset){
          
          weight <- sapply(tips,function(k) {z <- duplicated_indivs[duplicated_indivs[,1]==arbre_filtered$t[k],2]; if(length(z)>0){return(z) }else{return(0)}})
          weight <- weight+1
          NOT_homoplasy_count <- sum(weight[alleles=="ref"])
          homoplasy_count <- sum(weight[alleles=="not_ref"])
        }else{
          
          NOT_homoplasy_count <- sum(alleles=="ref")
          homoplasy_count <- sum(alleles=="not_ref")
          if(simulated_t_a_dataset){
            homoplasy_count <- sum(alleles=="t" | alleles=="T")
            NOT_homoplasy_count <- length(alleles)-homoplasy_count
          }
          
        }
        if( homoplasy_count >= tips_rule & NOT_homoplasy_count >= tips_rule ){return(1)}else{return(0)}
      } ) )
      
    }else{ count_children_node_with_homoplasy <- c(0) } # end children node rule
    if( sum(count_children_node_with_homoplasy) == 0 ){
      
      NOT_homoplasy_count <- 0
      homoplasy_count <- 0
      
      tips <- architecture_arbre[[node]][architecture_arbre[[node]]<=(length(arbre_filtered$tip.label))]
      alleles <- sapply(tips,function(x) input[as.character(homoplasy),arbre_filtered$t[x]])
      
      if(deduplicated_dataset){
        weight <- sapply(tips,function(x) {z <- duplicated_indivs[duplicated_indivs[,1]==arbre_filtered$t[x],2]; if(length(z)>0){return(z) }else{return(0)}})
        weight <- weight+1
        NOT_homoplasy_count <- sum(weight[alleles=="ref"])
        homoplasy_count <- sum(weight[alleles=="not_ref"])
      }else{
        NOT_homoplasy_count <- sum(alleles=="ref")
        homoplasy_count <- sum(alleles=="not_ref")
        if(simulated_t_a_dataset){
          homoplasy_count <- sum(alleles=="t" | alleles=="T")
          NOT_homoplasy_count <- length(alleles)-homoplasy_count
        }
      }
      
      if(homoplasy_count >= tips_rule & NOT_homoplasy_count >= tips_rule){
        # for non-paralelized version
        # i <- i+1 ; raw_out_table[i,] <- c(homoplasy,homoplasy_count,NOT_homoplasy_count)
        
        temp_raw_out_table <- rbind(temp_raw_out_table,c(homoplasy,homoplasy_count,NOT_homoplasy_count))
      }
    }
  }
  return(temp_raw_out_table)
}
stopCluster(cl)  

#########################################################################
#########################  Raw data saving  #############################
#########################################################################

save.image(file = paste("Workspace_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".RData",sep=""))
# load("Workspace_MinOffspring_4_embedded_1_MinTipsOfEachAllele_2.RData")
colnames(raw_out_table) <- c("position","homoplasy_count","NOT_homoplasy_count")
raw_out_table <- as.data.frame(raw_out_table[!is.na(raw_out_table[,1]),])

write.table(raw_out_table, file = paste("RAW_data_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".tsv",sep=""), append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"), fileEncoding = "")

#########################################################################
######################### Exploratory plots #############################
#########################################################################

out_table_long <- melt(raw_out_table, id=c("position"), measured=c("homoplasy_count","NOT_homoplasy_count"))
effectif <- nrow(raw_out_table)

# Density plot of values of homoplasy descendant counts and non-homoplasy descendant counts
jpeg(filename = paste("fig_1.non_filtered_count_histogram_min_offspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".jpg",sep=""), width = 1000, height = 700)
ggplot(out_table_long, mapping = aes (fill = variable, x = log(value))) + geom_density (alpha = .5) + ggtitle(paste("Homoplasy allele count vs non-homoplasy allele count\n non filtered dataset, n= ",effectif,sep=""))
dev.off()

out_table_ratios <- as.data.frame(cbind("position"=raw_out_table[,1],"homo/not_homo"=raw_out_table[,2]/raw_out_table[,3],"not_homo/homo"=raw_out_table[,3]/raw_out_table[,2]))
out_table_ratios_long <- melt(out_table_ratios, id=c("position"), measured=c("homo/not_homo","not_homo/homo"))

effectif <- sum(out_table_ratios_long[,2]=="homo/not_homo")
nb_unique_nodes <- length(unique(out_table_ratios_long[,1]))

# boxplot of all log10(RoHO) values along the genome. Unfiltered.
jpeg(filename = paste("fig_2.homo_vs_nothomo_ratio_along_genome_min_offspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".jpg",sep=""), width = 1000, height = 700)
ggplot(out_table_ratios_long[out_table_ratios_long[,2]=="homo/not_homo",],
       aes(x=position, y=log10(value), group=position))+
  geom_boxplot(outlier.shape = NA)+
  stat_summary(geom="point", fun.y=median)+
  geom_hline(yintercept = median(log10(out_table_ratios_long[out_table_ratios_long[,2]=="homo/not_homo","value"])), linetype='dashed')+
  ggtitle(paste("log of the ratio Homoplasy allele count / non-homoplasy allele count\nnb of studied homoplasy positions = ",nb_unique_nodes,", nb of studied paired count values= ",effectif,"\ndashed line = mediane",sep=""))
dev.off()

# Density plot of values of log10(RoHO)
jpeg(filename = paste("fig_4.non_filtered_ratio_histogram_min_offspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".jpg",sep=""), width = 1000, height = 700)
ggplot(out_table_ratios_long[out_table_ratios_long[,2]=="homo/not_homo",], mapping = aes (x = log(value))) + geom_density (alpha = .5) + ggtitle(paste("log of the ratio Homoplasy allele count / non-homoplasy allele count\n non filtered dataset, n= ",effectif,sep=""))
dev.off()


#####################################################################################
############ same plots with filtering of the data : nb of replicates ###############
############ eg. 5 nodes in the phylogeny detected as origin of a new homoplasy #####
#####################################################################################

nb_rep <- 3
# nb_rep <- 5

# will define limits for the plots in order to get the exact same plots for both filters
ylimplot_min  <- c(-3)
ylimplot_max <- c(3)

out_table_restricted <- raw_out_table[raw_out_table[,1]%in%names(table(raw_out_table[,1])[table(raw_out_table[,1])>=nb_rep]),]

######################################################################################
############################# STATISTICAL ANALYSIS ###################################
######################################################################################

library(BSDA)

multiple_t_test <- sapply(sort(unique(out_table_restricted[,1])), function(x) {
  res <- t.test(out_table_restricted[out_table_restricted[,1] == x,2],out_table_restricted[out_table_restricted[,1] == x,3], paired=TRUE)
  meanlog10roho <- mean(log10(out_table_restricted[out_table_restricted[,1] == x,2]/out_table_restricted[out_table_restricted[,1] == x,3]))
  medianlog10roho <- median(log10(out_table_restricted[out_table_restricted[,1] == x,2]/out_table_restricted[out_table_restricted[,1] == x,3]))
  sdlog10roho <- sd(log10(out_table_restricted[out_table_restricted[,1] == x,2]/out_table_restricted[out_table_restricted[,1] == x,3]))
  shapiro <- shapiro.test(out_table_restricted[out_table_restricted[,1] == x,2]-out_table_restricted[out_table_restricted[,1] == x,3])$p.value
  sign_test <- SIGN.test(out_table_restricted[out_table_restricted[,1] == x,2]-out_table_restricted[out_table_restricted[,1] == x,3], y = NULL, md = 0, alternative = "two.sided", conf.level = 0.95)$p.value
  sign_test_bonferroni <- p.adjust(sign_test, method = "bonferroni", n=length(unique(out_table_restricted[,1])))
  wilcox_test <- wilcox.test(out_table_restricted[out_table_restricted[,1] == x,2]-out_table_restricted[out_table_restricted[,1] == x,3], paired = F, mu=0, alternative = "two.sided")$p.value
  nb_replicates <- length(out_table_restricted[out_table_restricted[,1] == x,2])
  return(c(position=x,Mean_log10_RoHO=meanlog10roho,Median_log10_RoHO=medianlog10roho,SD_log10_RoHO=sdlog10roho,t_test=res$p.value,mean_of_diff=res$estimate,shapiro=shapiro,sign_test=sign_test,sign_test_bonferroni=sign_test_bonferroni,wilcox_test=wilcox_test,nb_replicates=nb_replicates))
})


write.table(t(multiple_t_test), file = paste("T_testMinReplicates_",nb_rep,"_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".tsv",sep=""), append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"), fileEncoding = "")

out_table_ratios_restricted_long <- as.data.frame(cbind(position=out_table_restricted[,1],"homo/not_homo"=out_table_restricted[,2]/out_table_restricted[,3]))

write.table(cbind(out_table_ratios_restricted_long,dataset=rep(folder_output,nrow(out_table_ratios_restricted_long))), file = paste("Homoplasy_RAW_posMinReplicates_",nb_rep,"_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".tsv",sep=""), append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T, qmethod = c("escape", "double"), fileEncoding = "")

custom_x_titles <- c(0,sort(unique(out_table_ratios_restricted_long[,1])),29903)

table_effectifs <- table(out_table_ratios_restricted_long[,1])
table_effectifs <- as.data.frame(t(rbind(pos=as.numeric(names(table_effectifs)),effectif=as.numeric(unname(table_effectifs)))))

if(simulated_t_a_dataset){
  fake_correspondance <- floor(runif(unique(out_table_ratios_restricted_long$position),min=100,max=29000))
  names(fake_correspondance) <- unique(out_table_ratios_restricted_long$position)
  out_table_ratios_restricted_long$position <- fake_correspondance[match(out_table_ratios_restricted_long$position,names(fake_correspondance))]
}

# for figure :
svg(filename = paste("fig_5.RoHo_MinReplicates_",nb_rep,"_MinOffspring_",min_offspring,"_embedded_",children_node_rule,"_MinTipsOfEachAllele_",tips_rule,".svg",sep=""), width = 10, height = 5)
ggplot(out_table_ratios_restricted_long, aes(x=position, y=log10(out_table_ratios_restricted_long[,2]), group=position))+
  coord_cartesian(xlim = c(0, 29903), ylim = c(ylimplot_min,ylimplot_max)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_vline(xintercept=table_effectifs$pos, linetype="dotted", color = "grey")+ 
  geom_violin(trim=T, fill='#A4A4A4', color="darkred",width=400,position=position_dodge(width = NULL)) +
  stat_summary(geom="point", fun.y=median)+
  scale_y_continuous(breaks=seq(from=trunc(ylimplot_min-1),to=trunc(ylimplot_max+1),by=1))+
  scale_x_continuous(breaks=custom_x_titles)+
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(geom = 'text', label = paste("n=",table_effectifs$effectif,sep=""), fun.y = function(x) max(x) + 0.5, vjust = 0,hjust = 0,angle = 45) +
  stat_summary(fun.y=mean, geom="point", shape=1, size=3, color="blue")
dev.off()
