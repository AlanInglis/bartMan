library(tidygraph)
library(ggraph)
library(dplyr)
library(gridExtra)
library(cowplot)
library(grid)
library(patchwork)

# create some data
nodes <- tibble(
  var = c("x4", "x1", NA, NA, NA),
  size = c(100, 65, 50,  15, 35)
)

edges <- tibble(
  from = c(1,2,2,1),
  to   = c(2,3,4,5)
)

# turn in tidygraph object
tg <- tbl_graph(nodes = nodes, edges = edges)
tg1 <- tbl_graph(nodes = nodes1, edges = edges1)

# plot width
g <- ggraph(tg, "partition", weight = size) +
  geom_node_tile(aes(fill = var)) +
  geom_node_label(aes(label = size, color = var)) +
  scale_y_reverse() +
  theme_void()+
  theme(legend.position = "none")
g

# width and height
g1 <- ggraph(tg1, "partition", weight = size) +
  geom_node_tile(aes(fill = var)) +
  geom_node_label(aes(label = size, color = var)) +
  scale_y_reverse() +
  theme_void()+
  theme(legend.position = "none")

g1

ggplot_build(g)$layout$panel_params[[1]]$y.range
ggplot_build(g1)$layout$panel_params[[1]]$y.range

g1 + ylim(-3.1, -0.9)

gridExtra::grid.arrange(g, g1, ncol = 2, heights = unit(0.5, "npc"))

grid.arrange(arrangeGrob(g, ncol=1, nrow = 1),
             arrangeGrob(g1, ncol=1, nrow = 1), heights=c(3,2))


plot_grid(g, g1, ncol = 2, rel_heights = c(3,2))

grid.draw(cbind(ggplotGrob(g), ggplotGrob(g1), size= "min"))

blank<-rectGrob(gp = gpar(col="white"))
grid.arrange(g, g1, blank, ncol=3, heights = c(1/2, 1/3, 1/9))

p = rectGrob()
grid.arrange(p, arrangeGrob(p,p, heights=c(3/4, 1/4), ncol=2),
             ncol=2)
grid.arrange(p,p, layout_matrix = cbind(c(1,1), c(2,3)))

g + (g1 *g)


ggraph(tg, "partition") +
  geom_node_tile(aes(fill = var, height = size/50)) +
  geom_node_label(aes(label = size, color = var)) +
  scale_y_reverse() +
  theme_void()+
  theme(legend.position = "none")

#
# ggplot(nodes, aes(var, size, height = size/4, width = size)) +
#   geom_tile(aes(fill = size))+
#   coord_fixed(ratio = 1)+
#   theme_void() +
#   theme(legend.position = 'none')
#
#
# ggplot(nodes, aes(var, group = size, fill = size)) +
#   geom_bar(position = 'stack',  width = 1) +
#   geom_text(aes(label = size), position = position_stack(vjust = 0.5), stat = 'count')
#
#
#
#
