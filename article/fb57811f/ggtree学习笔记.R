library(ggplot2)
library(ggtree)
library(ape)

set.seed(2017)
tree <- rtree(4)
tree

library(tidytree)
x <- as_tibble(tree)
x
as.phylo(x)

d <- tibble(label = paste0('t', 1:4),
            trait = rnorm(4))

y <- full_join(x, d, by = 'label')
y
z = as.treedata(y)
as_tibble(z)

child(y, 7) # 选择父节点为5的子节点
parent(y, 2)
offspring(y, 5)


library(treeio)
beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
rst_file <- system.file("examples/rst", package="ggtree")
mlc_file <- system.file("examples/mlc", package="ggtree")
beast_tree <- read.beast(beast_file)
codeml_tree <- read.codeml(rst_file, mlc_file)

merged_tree <- merge_tree(beast_tree, codeml_tree)
merged_tree


library(dplyr)
df <- merged_tree %>% 
  as_tibble() %>%
  select(dN_vs_dS, dN, dS, rate) %>%
  subset(dN_vs_dS >=0 & dN_vs_dS <= 1.5) %>%
  tidyr::gather(type, value, dN_vs_dS:dS)
df$type[df$type == 'dN_vs_dS'] <- 'dN/dS'
df$type <- factor(df$type, levels=c("dN/dS", "dN", "dS"))
ggplot(df, aes(rate, value)) + geom_hex() + 
  facet_wrap(~type, scale='free_y') 


library(ape)
data(woodmouse)
d <- dist.dna(woodmouse)
tr <- nj(d)
bp <- boot.phylo(tr, 
                 woodmouse, 
                 function(x) nj(dist.dna(x)))

bp2 <- tibble(node=1:Nnode(tr) + # 计算父节点数
                Ntip(tr), # 计算tip数
              bootstrap = bp)
full_join(tr, bp2, by="node")


file <- system.file("extdata/BEAST", "beast_mcc.tree", package="treeio")
beast <- read.beast(file)
x <- tibble(label = as.phylo(beast)$tip.label, trait = rnorm(Ntip(beast)))
full_join(beast, x, by="label")

nwk <- '(((((((A:4,B:4):6,C:5):8,D:6):3,E:21):10,((F:4,G:12):14,H:8):13):13,((I:5,J:2):30,(K:11,L:11):2):17):4,M:56);'
tree <- read.tree(text=nwk)

groupClade(as_tibble(tree), c(17, 21))



set.seed(2017)
tr <- rtree(4)
x <- as_tibble(tr)
## the input nodes can be node ID or label
groupOTU(x, c('t1', 't4'), group_name = "fake_group")


p1 <- ggtree(merged_tree) + theme_tree2()
p2 <- ggtree(rescale_tree(merged_tree, 'dN')) + theme_tree2()
p3 <- ggtree(rescale_tree(merged_tree, 'rate')) + theme_tree2()

cowplot::plot_grid(p1, p2, p3, ncol=3, labels = LETTERS[1:3])


library(ggplot2)
f <- system.file("extdata/NHX", "phyldog.nhx", package="treeio")
nhx <- read.nhx(f)
to_drop <- c("Physonect_sp_@2066767",
             "Lychnagalma_utricularia@2253871",
             "Kephyes_ovata@2606431")
p1 <- ggtree(nhx) + geom_tiplab(aes(color = label %in% to_drop)) +
  scale_color_manual(values=c("black", "red")) + xlim(0, 0.8)

nhx_reduced <- drop.tip(as.phylo(nhx), to_drop)
p2 <- ggtree(nhx_reduced) + geom_tiplab() + xlim(0, 0.8)  
cowplot::plot_grid(p1, p2, ncol=2, labels = c("A", "B"))



beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)



p1 = ggtree(beast_tree) + 
  geom_tiplab() +  
  ggtitle('原始树') +
  xlim(0, 40) + theme_tree2()

tree2 = tree_subset(beast_tree, "A/Swine/HK/168/2012", levels_back=4)  
p2 <- ggtree(tree2, aes(color=group)) +
  ggtitle('取子集') +
  scale_color_manual(values = c("black", "red")) +
  geom_tiplab() +  xlim(0, 4) + theme_tree2() 

p3 <- ggtree(tree2, aes(color=group)) +
  geom_tiplab(hjust = -.1) + xlim(0, 5) + 
  geom_point(aes(fill = rate), shape = 21, size = 4) +
  ggtitle('用rate这个变量控制颜色') +
  scale_color_manual(values = c("black", "red"), guide = FALSE) +
  scale_fill_continuous(low = 'blue', high = 'red') +
  theme_tree2() + theme(legend.position = 'right')


p4 <- ggtree(tree2, aes(color=group), 
             root.position = as.phylo(tree2)$root.edge) +
  geom_tiplab() + xlim(18, 24) + 
  ggtitle('添加根节点但不显示') +
  scale_color_manual(values = c("black", "red")) +
  theme_tree2()

p5 <- ggtree(tree2, aes(color=group), 
             root.position = as.phylo(tree2)$root.edge) +
  geom_rootedge() + geom_tiplab() + xlim(0, 40) + 
  ggtitle('添加根节点且显示') +
  scale_color_manual(values = c("black", "red")) +
  theme_tree2()

plot_grid(p2, p3, p4, p5, ncol=2) %>%
  plot_grid(p1, ., ncol=2)







beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
tr <- read.beast(beast_file)

library(ggtree)
p = ggtree(tr, ladderize=F) + geom_tiplab() + ggtitle("ggtree")
edge=data.frame(as.phylo(tr)$edge, 
                edge_num=1:nrow(as.phylo(tr)$edge))
colnames(edge)=c("parent", "node", "edge_num")
p %<+% edge + geom_label(aes(x=branch, label=edge_num))


getNodeNum(tr)
Nnode2(tr)





clade <- tree_subset(beast_tree, node=121, levels_back=0)
clade2 <- tree_subset(beast_tree, node=121, levels_back=2)
p1 <- ggtree(clade) + 
  ggtitle('感兴趣的整个分支') +
  geom_tiplab() + xlim(0, 5)
p2 <- ggtree(clade2, aes(color=group)) + 
  ggtitle('感兴趣的整个分支 + 回退两个节点') +
  geom_tiplab() + 
  xlim(0, 8) + scale_color_manual(values=c("black", "red"))


library(ape)
library(tidytree)
library(treeio)
data(chiroptera)

# 如果不知道node的时候，就直接进行匹配
nodes <- grep("Plecotus", chiroptera$tip.label)
chiroptera <- groupOTU(chiroptera, nodes)

p3 <- ggtree(chiroptera, aes(colour = group)) + 
  ggtitle('整个进化树中选择感兴趣的整个分支') +
  scale_color_manual(values=c("black", "red")) +
  theme(legend.position = "none")

clade <- MRCA(chiroptera, nodes) # 最近的父节点
x <- tree_subset(chiroptera, clade, levels_back = 0)
p4 <- ggtree(x) + 
  ggtitle('感兴趣的特有分支') +
  geom_tiplab() + xlim(0, 5)

plot_grid(p1, p2, p3, p4,ncol=2)

library(ggplot2)
library(treeio)
library(ggtree)

nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)

ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()
ggtree(tree)

p1 = ggtree(tree, color="firebrick", size=2, linetype="dotted") +
  ggtitle('阶梯化')
p2 = ggtree(tree, color="firebrick", size=2, linetype="dotted", ladderize=FALSE) +
  ggtitle('非阶梯化')

cowplot::plot_grid(p1,p2)

ggtree(tree, branch.length="none")



library(ggtree)
set.seed(2017-02-16)
tree <- rtree(50)
p1 = ggtree(tree) +
  ggtitle('默认')
p2 = ggtree(tree, layout="roundrect")  +
  ggtitle('roundrect')
p3 = ggtree(tree, layout="slanted") +
  ggtitle('slanted')
p4 = ggtree(tree, layout="ellipse") +
  ggtitle('ellipse')
p5 = ggtree(tree, layout="circular") +
  ggtitle('circular')
p6 = ggtree(tree, layout="fan", open.angle=120) +
  ggtitle('fan')
p7 = ggtree(tree, layout="equal_angle") +
  ggtitle('equal_angle')
p8 = ggtree(tree, layout="daylight") +
  ggtitle('daylight')
p9 = ggtree(tree, branch.length='none') +
  ggtitle('none')
p10 = ggtree(tree, layout="ellipse", branch.length="none") +
  ggtitle('ellipse对齐')
p11 = ggtree(tree, branch.length='none', layout='circular') +
  ggtitle('circular对齐')
p12 = ggtree(tree, layout="daylight", branch.length = 'none') +
  ggtitle('daylight对齐')

cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, ncol = 4)


beast_file <- system.file("examples/MCC_FluA_H3.tree", 
                          package="ggtree")
beast_tree <- read.beast(beast_file)
ggtree(beast_tree, mrsd="2013-01-01") + theme_tree2()



p1 = ggtree(tree) + geom_treescale() + ggtitle('默认')
p2 = ggtree(tree) + geom_treescale(x=0, y=45, width=1, color='red') +
  ggtitle('设定位置、宽度、颜色')
p3 = ggtree(tree) + geom_treescale(fontsize=6, linesize=2, offset=1) +
  ggtitle('设定字体大小、线条大小、缩进')
p4 = ggtree(tree) + theme_tree2() +
  ggtitle('使用内置主题')

cowplot::plot_grid(p1,p2,p3,p4, ncol = 2)



p1 = ggtree(tree) + geom_point(aes(shape=isTip, color=isTip), size=3) +
  ggtitle('使用geom_point()函数')
p2 <- ggtree(tree) + geom_nodepoint(color="#b5e521", alpha=1/4, size=10) +
  geom_tippoint(color="#FDAC4F", shape=8, size=3) +
  ggtitle('使用两个函数')

cowplot::plot_grid(p1,p2, ncol = 2)


p1 = ggtree(tree) + geom_nodepoint(color="#b5e521", alpha=1/4, size=10) +
  geom_tippoint(color="#FDAC4F", shape=8, size=3) + 
  geom_tiplab(size=3, color="purple")

p2 = ggtree(tree, layout="circular") + 
  geom_tiplab(aes(angle=angle), color='blue')

p3 = ggtree(tree, branch.length = 'none') + 
  geom_tiplab(as_ylab=TRUE, color='firebrick')

cowplot::plot_grid(p1, p2,p3, ncol = 3)





## with root edge = 1
tree1 <- read.tree(text='((A:1,B:2):3,C:2):1;')
p1 = ggtree(tree1) + geom_tiplab() + geom_rootedge() +
  ggtitle('有根节点信息')

## without root edge
tree2 <- read.tree(text='((A:1,B:2):3,C:2);')
p2 = ggtree(tree2) + geom_tiplab() + geom_rootedge() +
  ggtitle('无根节点信息，默认无')

## setting root edge
tree2$root.edge <- 2
p3 = ggtree(tree2) + geom_tiplab() + geom_rootedge() +
  ggtitle('无根节点信息，添加信息')

## specify length of root edge for just plotting
## this will ignore tree$root.edge
p4 = ggtree(tree2) + geom_tiplab() + geom_rootedge(rootedge = 3) +
  ggtitle('无根节点信息，设置信息')

cowplot::plot_grid(p1,p2,p3,p4, ncol = 2)

ggtree(beast_tree, aes(color=rate)) +
  scale_color_continuous(low='darkgreen', high='red') +
  theme(legend.position="right")





anole.tree<-read.tree("http://www.phytools.org/eqg2015/data/anole.tre")
svl <- read.csv("http://www.phytools.org/eqg2015/data/svl.csv",
                row.names=1)
svl <- as.matrix(svl)[,1]
fit <- phytools::fastAnc(anole.tree,svl,vars=TRUE,CI=TRUE)

td <- data.frame(node = nodeid(anole.tree, names(svl)),
                 trait = svl)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)

d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree <- full_join(anole.tree, d, by = 'node')

p1 <- ggtree(tree, aes(color=trait), layout = 'circular', 
             ladderize = FALSE, continuous = TRUE, size=2) +
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
  geom_tiplab(hjust = -.1) + 
  xlim(0, 1.2) + 
  theme(legend.position = c(.05, .85)) 

p2 <- ggtree(tree, layout='circular', ladderize = FALSE, size=2.8) + 
  geom_tree(aes(color=trait), continuous=T, size=2) +  
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
  geom_tiplab(aes(color=trait), hjust = -.1) + 
  xlim(0, 1.2) + 
  theme(legend.position = c(.05, .85)) 

cowplot::plot_grid(p1, p2, ncol=2, labels=c("分支默认边框", "分支黑色边框"))    

library(treeio)
beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)
beast_tree

p1 <- ggtree(beast_tree, mrsd='2013-01-01') + theme_tree2() +
  labs(caption="时间序列")
p2 <- ggtree(beast_tree, branch.length='rate') + theme_tree2() +
  labs(caption="取代速率")

mlcfile <- system.file("extdata/PAML_Codeml", "mlc", package="treeio")
mlc_tree <- read.codeml_mlc(mlcfile)
p3 <- ggtree(mlc_tree) + theme_tree2() +
  labs(caption="单密码子核苷酸取代")
p4 <- ggtree(mlc_tree, branch.length='dN_vs_dS') + theme_tree2() +
  labs(caption="dN/dS")

cowplot::plot_grid(p1,p2,p3,p4, ncol = 2)


beast_tree2 <- rescale_tree(beast_tree, branch.length='rate')
ggtree(beast_tree2) + theme_tree2()








set.seed(2019)
x <- rtree(30)
p1 = ggtree(x, color="red") + theme_tree("steelblue")
p2 = ggtree(x, color="white") + theme_tree("black")

cowplot::plot_grid(p1,p2, ncol = 2)




## trees <- lapply(c(10, 20, 40), rtree)
## class(trees) <- "multiPhylo"
## ggtree(trees) + facet_wrap(~.id, scale="free") + geom_tiplab()

r8s <- read.r8s(system.file("extdata/r8s", "H3_r8s_output.log", package="treeio"))
ggtree(r8s) + facet_wrap( ~.id, scale="free") + theme_tree2()




library(ggtree)
treetext = "(((ADH2:0.1[&&NHX:S=human], ADH1:0.11[&&NHX:S=human]):
0.05 [&&NHX:S=primates:D=Y:B=100],ADHY:
0.1[&&NHX:S=nematode],ADHX:0.12 [&&NHX:S=insect]):
0.1[&&NHX:S=metazoa:D=N],(ADH4:0.09[&&NHX:S=yeast],
ADH3:0.13[&&NHX:S=yeast], ADH2:0.12[&&NHX:S=yeast],
ADH1:0.11[&&NHX:S=yeast]):0.1[&&NHX:S=Fungi])[&&NHX:D=N];"
tree <- read.nhx(textConnection(treetext))
ggtree(tree) + geom_tiplab() + 
  geom_label(aes(x=branch, label=S), fill='lightgreen') + 
  geom_label(aes(label=D), fill='steelblue') + 
  geom_text(aes(label=B), hjust=-.5)



set.seed(2015-12-21)
tree <- rtree(30)
p1 <- ggtree(tree) + xlim(NA, 6)

p2 = p1 + geom_cladelabel(node=45, label="test label") +
  geom_cladelabel(node=34, label="another clade")

p3 = p1 + geom_cladelabel(node=45, label="test label", align=TRUE,  offset = .2, color='red') +
  geom_cladelabel(node=34, label="another clade", align=TRUE, offset = .2, color='blue')

p4 = p1 + geom_cladelabel(node=45, label="test label", align=T, angle=270, hjust='center', offset.text=.5, barsize=1.5) +
  geom_cladelabel(node=34, label="another clade", align=T, angle=45, fontsize=8)

p5 = p1 + geom_cladelabel(node=34, label="another clade", align=T, geom='label', fill='lightblue')

cowplot::plot_grid(p2,p3,p4,p5,ncol = 2)



p1 = ggtree(tree, layout="daylight") + 
  geom_cladelabel(node=35, label="test label", angle=0, 
                  fontsize=8, offset=.5, vjust=.5)  + 
  geom_cladelabel(node=55, label='another clade', 
                  angle=-95, hjust=.5, fontsize=8)


p2 = ggtree(tree) + xlim(NA, 6) + 
  geom_tiplab() + 
  geom_strip('t10', 't30', barsize=2, color='red', 
             label="associated taxa", offset.text=.1) + 
  geom_strip('t1', 't18', barsize=2, color='blue', 
             label = "another label", offset.text=.1)

cowplot::plot_grid(p1,p2,ncol = 2)











nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)

p1 = ggtree(tree) + 
  geom_hilight(node=21, fill="steelblue", alpha=.6) +
  geom_hilight(node=17, fill="darkgreen", alpha=.6) 

p2 = ggtree(tree, layout="circular") + 
  geom_hilight(node=21, fill="steelblue", alpha=.6) +
  geom_hilight(node=23, fill="darkgreen", alpha=.6)


## type can be 'encircle' or 'rect'
ggtree(tree, layout="daylight", branch.length = 'none') + 
  geom_hilight(node=10) + 
  geom_hilight(node=16, fill='darkgreen', type="rect")


ggtree(tree) +
  geom_balance(node=16, fill='steelblue', color='white', alpha=0.6, extend=1) +
  geom_balance(node=19, fill='darkgreen', color='white', alpha=0.6, extend=1) 

## using external data
d <- data.frame(node=c(17, 21), type=c("A", "B"))
ggtree(tree) + geom_hilight(data=d, aes(node=node, fill=type))

## using data stored in tree object
x <- read.nhx(system.file("extdata/NHX/ADH.nhx", package="treeio"))
ggtree(x) + 
  geom_hilight(mapping=aes(subset = node %in% c(10, 12), fill = S)) +
  scale_fill_manual(values=c("steelblue", "darkgreen"))









p1 <- ggtree(tree) + geom_tiplab() + 
  geom_taxalink(taxa1='A', taxa2='E') + 
  geom_taxalink(taxa1='F', taxa2='K', 
                color='red', linetype = 'dashed',
                arrow=arrow(length=unit(0.02, "npc")))

p2 <- ggtree(tree, layout="circular") + 
  geom_taxalink(taxa1='A', taxa2='E', 
                color="grey",alpha=0.5, 
                offset=0.05,arrow=arrow(length=unit(0.01, "npc"))) + 
  geom_taxalink(taxa1='F', taxa2='K', 
                color='red', linetype = 'dashed', 
                alpha=0.5, offset=0.05,
                arrow=arrow(length=unit(0.01, "npc"))) +
  geom_taxalink(taxa1="L", taxa2="M", 
                color="blue", alpha=0.5, 
                offset=0.05,hratio=0.8, 
                arrow=arrow(length=unit(0.01, "npc"))) + 
  geom_tiplab()

# when the tree was created using reverse x, 
# we can set outward to FALSE, which will generate the inward curve lines.
p3 <- ggtree(tree, layout="inward_circular", xlim=150) +
  geom_taxalink(taxa1='A', taxa2='E', 
                color="grey", alpha=0.5, 
                offset=-0.2, 
                outward=FALSE,
                arrow=arrow(length=unit(0.01, "npc"))) +
  geom_taxalink(taxa1='F', taxa2='K', 
                color='red', linetype = 'dashed', 
                alpha=0.5, offset=-0.2,
                outward=FALSE,
                arrow=arrow(length=unit(0.01, "npc"))) +
  geom_taxalink(taxa1="L", taxa2="M", 
                color="blue", alpha=0.5, 
                offset=-0.2, 
                outward=FALSE,
                arrow=arrow(length=unit(0.01, "npc"))) +
  geom_tiplab(hjust=1) 

dat <- data.frame(from=c("A", "F", "L"), 
                  to=c("E", "K", "M"), 
                  h=c(1, 1, 0.1), 
                  type=c("t1", "t2", "t3"), 
                  s=c(2, 1, 2))
p4 <- ggtree(tree, layout="inward_circular", xlim=c(150, 0)) +
  geom_taxalink(data=dat, 
                mapping=aes(taxa1=from, 
                            taxa2=to, 
                            color=type, 
                            size=s), 
                ncp=10,
                offset=0.15) + 
  geom_tiplab(hjust=1) +
  scale_size_continuous(range=c(1,3))
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])


file <- system.file("extdata/MEGA7", "mtCDNA_timetree.nex", package = "treeio")
x <- read.mega(file)
p1 <- ggtree(x) + geom_range('reltime_0.95_CI', color='red', size=3, alpha=.3)
p2 <- ggtree(x) + geom_range('reltime_0.95_CI', color='red', size=3, alpha=.3, center='reltime')  
p3 <- p2 + scale_x_range() + theme_tree2()

cowplot::plot_grid(p1,p2,p3,ncol = 3)





file <- system.file("extdata/BEAST", "beast_mcc.tree", package="treeio")
beast <- read.beast(file)
ggtree(beast, aes(color=rate))  +
  geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=3) +
  scale_color_continuous(low="darkgreen", high="red") +
  theme(legend.position=c(.1, .8))













rstfile <- system.file("extdata/PAML_Codeml", "rst", 
                       package="treeio")
mlcfile <- system.file("extdata/PAML_Codeml", "mlc", 
                       package="treeio")

ml <- read.codeml(rstfile, mlcfile)
ggtree(ml, aes(color=dN_vs_dS), branch.length='dN_vs_dS') + 
  scale_color_continuous(name='dN/dS', limits=c(0, 1.5),
                         oob=scales::squish,
                         low='darkgreen', high='red') +
  geom_text(aes(x=branch, label=AA_subs), 
            vjust=-.5, color='steelblue', size=2) +
  theme_tree2(legend.position=c(.9, .3))







library(ggtree)
nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)
p1 = ggtree(tree) + geom_tiplab()
p2 = viewClade(p, MRCA(p, "I", "L"))

cowplot::plot_grid(p1,p2,ncol = 2, labels = c('原图','特定区域'))



tree2 <- groupClade(tree, c(17, 21))
p <- ggtree(tree2, aes(color=group)) + theme(legend.position='none') +
  scale_color_manual(values=c("black", "firebrick", "steelblue"))
scaleClade(p, node=17, scale=.1) 






p2 <- p + geom_tiplab()
node <- 21
collapse(p2, node, 'max') %>% expand(node)
collapse(p2, node, 'min') %>% expand(node)
collapse(p2, node, 'mixed') %>% expand(node)




library(ggsci)

data(iris)
rn <- paste0(iris[,5], "_", 1:150)
rownames(iris) <- rn
d_iris <- dist(iris[,-5], method="man")

c <- ape::bionj(d_iris)
grp <- list(setosa     = rn[1:50],
            versicolor = rn[51:100],
            virginica  = rn[101:150])

p_iris <- ggtree(tree_iris, layout = 'circular', branch.length='none')
groupOTU(p_iris, grp, 'group') + 
  aes(color=group) +
  scale_color_aaas() +
  theme(legend.position="right")




library(ggimage)
library(ggtree)

# 文件下载地址
# https://raw.githubusercontent.com/YuLab-SMU/treedata-book/master/data/tree_boots.nwk
# https://raw.githubusercontent.com/YuLab-SMU/treedata-book/master/data/tip_data.csv


x <- read.tree("tree_boots.nwk")
info <- read.csv("tip_data.csv")

p <- ggtree(x) %<+% info + xlim(-.1, 4)
p2 <- p + geom_tiplab(offset = .6, hjust = .5) +
  geom_tippoint(aes(shape = trophic_habit, color = trophic_habit, size = mass_in_kg)) + 
  theme(legend.position = "right") + scale_size_continuous(range = c(3, 10))

#https://raw.githubusercontent.com/YuLab-SMU/treedata-book/master/data/inode_data.csv
d2 <- read.csv("inode_data.csv")
p2 %<+% d2 + geom_label(aes(label = vernacularName.y, fill = posterior)) + 
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(3, "YlGnBu"))










library(ggtree)
## remote_folder <- paste0("https://raw.githubusercontent.com/katholt/",
##                         "plotTree/master/tree_example_april2015/")
remote_folder <- "data/tree_example_april2015/" 

## read the phylogenetic tree
tree <- read.tree(paste0(remote_folder, "tree.nwk"))

## read the sampling information data set
info <- read.csv(paste0(remote_folder,"info.csv"))

## read and process the allele table
snps<-read.csv(paste0(remote_folder, "alleles.csv"), header = F,
               row.names = 1, stringsAsFactor = F)
snps_strainCols <- snps[1,] 
snps<-snps[-1,] # drop strain names
colnames(snps) <- snps_strainCols

gapChar <- "?"
snp <- t(snps)
lsnp <- apply(snp, 1, function(x) {
  x != snp[1,] & x != gapChar & snp[1,] != gapChar
})
lsnp <- as.data.frame(lsnp)
lsnp$pos <- as.numeric(rownames(lsnp))
lsnp <- tidyr::gather(lsnp, name, value, -pos)
snp_data <- lsnp[lsnp$value, c("name", "pos")]

## read the trait data
bar_data <- read.csv(paste0(remote_folder, "bar.csv"))

## visualize the tree 
p <- ggtree(tree) 

## attach the sampling information data set 
## and add symbols colored by location
p <- p %<+% info + geom_tippoint(aes(color=location))

## visualize SNP and Trait data using dot and bar charts,
## and align them based on tree structure
p + geom_facet(panel = "SNP", data = snp_data, geom = geom_point, 
               mapping=aes(x = pos, color = location), shape = '|') +
  geom_facet(panel = "Trait", data = bar_data, geom = ggstance::geom_barh, 
             aes(x = dummy_bar_value, color = location, fill = location), 
             stat = "identity", width = .6) +
  theme_tree2(legend.position=c(.05, .85))




library(ggplot2)
library(ggtree)

set.seed(2019-10-31)
tr <- rtree(10)

d1 <- data.frame(
  # only some labels match
  label = c(tr$tip.label[sample(5, 5)], "A"),
  value = sample(1:6, 6))

d2 <- data.frame(
  label = rep(tr$tip.label, 5),
  category = rep(LETTERS[1:5], each=10),
  value = rnorm(50, 0, 3)) 

g <- ggtree(tr) + geom_tiplab(align=TRUE)

p1 <- ggplot(d1, aes(label, value)) + geom_col(aes(fill=label)) + 
  geom_text(aes(label=label, y= value+.1)) +
  coord_flip() + theme_tree2() + theme(legend.position='none')

p2 <- ggplot(d2, aes(x=category, y=label)) + 
  geom_tile(aes(fill=value)) + scale_fill_viridis_c() + 
  theme_tree2() 

cowplot::plot_grid(g, p2, p1, ncol=3) 

library(aplot)
p2 %>% insert_left(g) %>% insert_right(p1, width=.5) 




















library(ggimage)
library(ggtree)

nwk <- "((((bufonidae, dendrobatidae), ceratophryidae), (centrolenidae, leptodactylidae)), hylidae);"

x = read.tree(text = nwk)
ggtree(x) + xlim(NA, 7) + ylim(NA, 6.2) +
  geom_tiplab(aes(image=paste0("img/frogs/", label, '.jpg')), 
              geom="image", offset=2, align=2, size=.2)  + 
  geom_tiplab(geom='label', offset=1, hjust=.5) + 
  geom_image(x=.8, y=5.5, image="img/frogs/frog.jpg", size=.2)









library(ggtree)
newick <- "((Pongo_abelii,(Gorilla_gorilla_gorilla,(Pan_paniscus,Pan_troglodytes)Pan,Homo_sapiens)Homininae)Hominidae,Nomascus_leucogenys)Hominoidea;"

tree <- read.tree(text=newick)

d <- ggimage::phylopic_uid(tree$tip.label)
d$body_mass = c(52, 114, 47, 45, 58, 6)

p <- ggtree(tree) %<+% d + 
  geom_tiplab(aes(image=uid, colour=body_mass), geom="phylopic", offset=2.5) +
  geom_tiplab(aes(label=label), offset = .2) + xlim(NA, 7) +
  scale_color_viridis_c()








library(phytools)
library(treeio)
library(tidytree)
data(anoletree)
x <- getStates(anoletree,"tips")
tree <- anoletree

cols <- setNames(palette()[1:length(unique(x))],sort(unique(x)))
fitER <- ape::ace(x,tree,model="ER",type="discrete")
ancstats <- as.data.frame(fitER$lik.anc)
ancstats$node <- 1:tree$Nnode+Ntip(tree)

## cols parameter indicate which columns store stats
bars <- nodebar(ancstats, cols=1:6)
bars <- lapply(bars, function(g) g+scale_fill_manual(values = cols))

tree2 <- full_join(tree, data.frame(label = names(x), stat = x ), by = 'label')
p <- ggtree(tree2) + geom_tiplab() +
  geom_tippoint(aes(color = stat)) + 
  scale_color_manual(values = cols) +
  theme(legend.position = "right") + 
  xlim(NA, 8)
p + geom_inset(bars, width = .08, height = .05, x = "branch") 






library(phytools)
library(treeio)
library(tidytree)
data(anoletree)
x <- getStates(anoletree,"tips")
tree <- anoletree

cols <- setNames(palette()[1:length(unique(x))],sort(unique(x)))
fitER <- ape::ace(x,tree,model="ER",type="discrete")
ancstats <- as.data.frame(fitER$lik.anc)
ancstats$node <- 1:tree$Nnode+Ntip(tree)

pies <- nodepie(ancstats, cols = 1:6)
pies <- lapply(pies, function(g) g+scale_fill_manual(values = cols))

tree2 <- full_join(tree, data.frame(label = names(x), stat = x ), by = 'label')
p <- ggtree(tree2) + geom_tiplab() +
  geom_tippoint(aes(color = stat)) + 
  scale_color_manual(values = cols) +
  theme(legend.position = "right") + 
  xlim(NA, 8)

p + geom_inset(pies, width = .1, height = .1) 


pies_and_bars <- pies
i <- sample(length(pies), 20)
pies_and_bars[i] <- bars[i]
p + geom_inset(pies_and_bars, width=.08, height=.05)


library(ggplot2)
library(ggtree)
# install.packages('emojifont')

tt = '((snail,mushroom),(((sunflower,evergreen_tree),leaves),green_salad));'
tree = read.tree(text = tt)
d <- data.frame(label = c('snail','mushroom', 'sunflower',
                          'evergreen_tree','leaves', 'green_salad'),
                group = c('animal', 'fungi', 'flowering plant',
                          'conifers', 'ferns', 'mosses'))

ggtree(tree, linetype = "dashed", size=1, color='firebrick') %<+% d + 
  xlim(0, 4.5) + ylim(0.5, 6.5) +
  geom_tiplab(parse="emoji", size=15, vjust=.25) +
  geom_tiplab(aes(label = group), geom="label", x=3.5, hjust=.5)


p <- ggtree(tree, layout = "circular", size=1) +  
  geom_tiplab(parse="emoji", size=10, vjust=.25)
print(p)

## fan layout  
p2 <- open_tree(p, angle=200) 
print(p2)

p2 %>% rotate_tree(-90)



set.seed(123)
tr <- rtree(30)

ggtree(tr) + xlim(NA, 5.2) +
  geom_cladelabel(node=41, label="chicken", parse="emoji",
                  fontsize=12, align=TRUE, colour="firebrick") +
  geom_cladelabel(node=53, label="duck", parse="emoji",
                  fontsize=12, align=TRUE, colour="steelblue") +
  geom_cladelabel(node=48, label="family", parse="emoji",
                  fontsize=12, align=TRUE, colour="darkkhaki")











hc <- hclust(dist(mtcars))
hc







rm(list = ls())
library(ggtreeExtra)
library(ggtree)
library(phyloseq)
library(dplyr)

data("GlobalPatterns")
GP <- GlobalPatterns
GP <- prune_taxa(taxa_sums(GP) > 600, GP)
sample_data(GP)$human <- get_variable(GP, "SampleType") %in%
  c("Feces", "Skin")
mergedGP <- merge_samples(GP, "SampleType")
mergedGP <- rarefy_even_depth(mergedGP,rngseed=394582)
mergedGP <- tax_glom(mergedGP,"Order")

melt_simple <- psmelt(mergedGP) %>%
  filter(Abundance < 120) %>%
  select(OTU, val=Abundance)

p <- ggtree(mergedGP, layout="fan", open.angle=10) + 
  geom_tippoint(mapping=aes(color=Phylum), 
                size=1.5,
                show.legend=FALSE)
p <- rotate_tree(p, -90)

p <- p +
  geom_fruit(
    data=melt_simple, # 数据源
    geom=geom_boxplot, # 图形样式
    mapping = aes(
      y=OTU,
      x=val,
      group=label,
      fill=Phylum,
    ),
    size=.2,
    outlier.size=0.5,
    outlier.stroke=0.08,
    outlier.shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 3,
    ),
    grid.params=list()
  ) 

p <- p +
  scale_fill_discrete(
    name="Phyla",
    guide=guide_legend(keywidth=0.8, keyheight=0.8, ncol=1)
  ) +
  theme(
    legend.title=element_text(size=9), # The title of legend 
    legend.text=element_text(size=7) # The label text of legend, the sizes should be adjust with dpi.
  )
p













































































































































