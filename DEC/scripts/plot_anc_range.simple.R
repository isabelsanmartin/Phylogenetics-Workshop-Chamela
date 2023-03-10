source("scripts/plot_anc_range.util.R")

# file names
fp = "./" # edit to provide an absolute filepath
plot_fn = paste(fp, "output/simple.range.pdf",sep="")
tree_fn = paste(fp, "output/simple.ase.tre", sep="")
label_fn = paste(fp, "output/simple.state_labels.txt", sep="")
color_fn = paste(fp, "range_colors.n4.txt", sep="")

#install.packages("devtools", dependencies=TRUE)
#library(devtools)
#install_github("GuangchuangYu/ggtree")
#install_github("revbayes/RevGadgets")
#library(RevGadgets)

# get state labels and state colors
labs <- c("1"  = "A",   "2"  = "B", 
          "3"  = "C",   "4"  = "D", 
          "5"  = "AB",  "6"  = "AC", 
          "7"  = "BC",  "8"  = "AD", 
          "9"  = "BD",  "10" = "CD", 
          "11" = "ABC", "12" = "ABD", 
          "13" = "ACD", "14" = "BCD", 
          "15" = "ABCD")
ancstates <- processAncStates(tree_fn, state_labels = labs)

# plot the ancestral states
pp=plotAncStatesPie(t = ancstates, 
                     # Include cladogenetic events
                     cladogenetic = T,
                     # adjust tip labels 
                     tip_labels_offset = 0.1)


# get plot dimensions
x_phy = max(pp$data$x)       # get height of tree
x_label = 3.5                # choose space for tip labels
x_start = 45                  # choose starting age (greater than x_phy)
x0 = -(x_start - x_phy)      # determine starting pos for xlim
x1 = x_phy + x_label         # determine ending pos for xlim

# add axis
pp = pp + theme_tree2()
pp = pp + labs(x="Age (Ma)")

# change x coordinates
pp = pp + coord_cartesian(xlim=c(x0,x1), expand=TRUE)

# plot axis ticks
island_axis = sec_axis(~ ., breaks=x_phy-c(35.0, 17.0, 10.0, 3.5), labels=c("TEE","MMCO","LMC","BS") )
x_breaks = seq(0,x_start,1) + x0
x_labels = rev(seq(0,x_start,1))
pp = pp + scale_x_continuous(breaks=x_breaks, labels=x_labels, sec.axis=island_axis)

pp

# save 
ggsave(file=plot_fn, plot=pp, device="pdf", height=14, width=30, useDingbats=F)

