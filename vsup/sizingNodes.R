library(tidygraph)
library(ggraph)
library(dplyr)

# create some data
nodes <- tibble(
  var = c("x4", "x1", NA, NA, NA),
  size = c(100, 65, 50,  35, 15)
)

edges <- tibble(
  from = c(1,2,2,1),
  to   = c(2,3,4,5)
)

# turn in tidygraph object
tg <- tbl_graph(nodes = nodes, edges = edges)

# plot width
ggraph(tg, "partition", weight = size) +
  geom_node_tile(aes(fill = var)) +
  geom_node_label(aes(label = size, color = var)) +
  scale_y_reverse() +
  theme_void()+
  theme(legend.position = "none")

# width and height
ggraph(tg, "partition", weight = size) +
  geom_node_tile(aes(fill = var, height = size/100)) +
  geom_node_label(aes(label = size, color = var)) +
  scale_y_reverse() +
  theme_void()+
  theme(legend.position = "none")

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
