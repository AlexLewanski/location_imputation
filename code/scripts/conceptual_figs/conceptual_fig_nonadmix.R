###################################################
###################################################
### CONCEPTUAL FIGURE FOR NON-ADMIXED INFERENCE ###
###################################################
###################################################

#####################
### SCRIPT SET-UP ###
#####################

#packages
library(here)
library(dplyr)
library(ggplot2)
library(plot3D)


nodes <- read.csv(here('simulation_output', 'conceptual_fig', 'nodes_nonadmix.csv'))
edges <- read.csv(here('simulation_output', 'conceptual_fig', 'edges_nonadmix.csv'))
span_info <- read.csv(here('simulation_output', 'conceptual_fig', 'span_info_nonadmix.csv'))


edge_df <- edges %>% 
  left_join(nodes, by = c('parent' = 'node_id', 'tree_ind')) %>% 
  rename(parent_time = time,
         parent_x = x,
         parent_y = y) %>% 
  left_join(nodes, by = c('child' = 'node_id', 'tree_ind')) %>% 
  rename(child_time = time,
         child_x = x,
         child_y = y)



phi <- 2
theta <- 220


infer_node <- 1


color_list <- list(
  c('#d85ed5', '#f3cef2', '#ac4baa'),
  c("#895129", "#e7dcd4", '#5f381c'),
  c('#48b46c', '#dcf4e4', '#2b6c40'),
  c("#ce380f", '#f5d7cf', '#90270a'),
  c("#254ad1", "#d4dcfa", '#162c7d')
)

color_list <- setNames(color_list, paste0('t', unique(edge_df$tree_ind)) )

for (TREE in unique(edge_df$tree_ind)) {
  tree_name <- paste0('tree_ind', TREE, '.png')
  png(here('plots', 'conceptual_figs', 'non_admix', tree_name), 
      width = 10, height = 10, units = 'cm', res = 500, bg = "transparent") 
  par(mar = c(0.1, 0.1, 0.1, 0.1), bg = "transparent")
  subset_edge <- edge_df[edge_df$tree_ind == TREE,]
  subset_nodes <- nodes[nodes$tree_ind == TREE,]
  
  rect3D(x0 = -2,x1 = 102,
         y0 = -2,y1 = 102,
         z0 = -5, z1 = NULL,
         col = adjustcolor('#e5e5e5',0.05), #adjustcolor(color_list[[paste0('t', TREE)]][2],0.05),
         border = NA, #color_list[[paste0('tree', TREE)]][2],
         bty='n',theta=theta,phi=phi,
         colkey=FALSE,xlab="",ylab="",zlab="",
         xlim=c(-2, 102),ylim=c(-2, 102),zlim=c(0, max(nodes$time) - 700),
         add=FALSE)
  
  title(main = paste0('Tree ', TREE), 
        line = -14, adj = 0.96, cex.main = 1,
        col.main = color_list[[paste0('t', TREE)]][1])
  
  for (i in c(1:nrow(subset_edge)) ) {
    new_df <- data.frame(x = unlist(subset_edge[i,c('parent_x', 'child_x')]), 
                         y = unlist(subset_edge[i,c('parent_y', 'child_y')]),
                         time = unlist(subset_edge[i,c('parent_time', 'child_time')]))
    
    lines3D(new_df$x, 
            new_df$y, 
            new_df$time,
            col = color_list[[paste0('t', TREE)]][1],
            lwd = 2,
            add = TRUE)
  }
  
  points3D(subset_nodes$x[subset_nodes$node_id != infer_node & subset_nodes$time == 0L], 
           subset_nodes$y[subset_nodes$node_id != infer_node & subset_nodes$time == 0L], 
           subset_nodes$time[subset_nodes$node_id != infer_node & subset_nodes$time == 0L],
           #col = "blue", 
           pch = 19, cex = 0.8, col = '#999999', #col = color_list[[paste0('t', TREE)]][1],
           xlab = "", ylab = "", zlab = "",
           #xlim = c(-10, 110), ylim = c(-10, 110),
           #theta = 40, phi = 10,
           bty = "n",
           add = TRUE)
  
  points3D(subset_nodes$x[subset_nodes$node_id != infer_node & subset_nodes$time > 0L], 
           subset_nodes$y[subset_nodes$node_id != infer_node & subset_nodes$time > 0L], 
           subset_nodes$time[subset_nodes$node_id != infer_node & subset_nodes$time > 0L],
           #col = "blue", 
           pch = 19, cex = 0.8, col = color_list[[paste0('t', TREE)]][1],
           xlab = "", ylab = "", zlab = "",
           #xlim = c(-10, 110), ylim = c(-10, 110),
           #theta = 40, phi = 10,
           bty = "n",
           add = TRUE)
  
  points3D(subset_nodes$x[subset_nodes$node_id == infer_node], 
           subset_nodes$y[subset_nodes$node_id == infer_node], 
           subset_nodes$time[subset_nodes$node_id == infer_node],
           #col = "blue", 
           pch = 19, cex = 1.5, col = 'black',
           xlab = "", ylab = "", zlab = "",
           #xlim = c(-10, 110), ylim = c(-10, 110),
           #theta = 40, phi = 10,
           bty = "n",
           add = TRUE)
  
  points3D(runif(1, 5, 70), 
           runif(1, 5, 70), 
           0,
           #col = "blue", 
           pch = 8, cex = 1.75, col = color_list[[paste0('t', TREE)]][3],
           xlab = "", ylab = "", zlab = "",
           #xlim = c(-10, 110), ylim = c(-10, 110),
           #theta = 40, phi = 10,
           bty = "n",
           add = TRUE)

  
  dev.off()
}


span_color <- sapply(color_list, function(x) x[1])
span_color[['n']] <- '#e5e5e5'

genome_span_plot <- span_info %>% 
  mutate(y1 = if_else(group == 'n', -1, -1.15),
         y2 = if_else(group == 'n', 1, 1.15)) %>% 
  ggplot() +
  geom_rect(aes(xmin = left, xmax = right, ymin = y1, ymax = y2, fill = group),
            color = 'white', lwd = 0.75) +
  geom_text(data = span_info %>% 
              filter(group != 'n') %>% 
              mutate(x_val = 0.5*(left + right),
                     group_num = gsub('t', '', group)),
            aes(x = x_val, y = -1.55, label = group_num, color = group),
            size = 3.75, fontface = "bold") +
  scale_fill_manual(values =  span_color) +
  scale_color_manual(values =  span_color) +
  theme_void() +
  theme(legend.position = 'none') +
  ylim(-1.8, 1.2)


ggsave(here('plots', 'conceptual_figs', 'non_admix', 'genome_span_fig.png'),
       genome_span_plot,
       width = 8, height = 0.5)


conditional_mean_df <- data.frame(
  tree = names(color_list),
  x = runif(5, 0, 100),
  y = runif(5, 0, 100)
)

sample_node_df <- nodes %>% 
  filter(time - 0 < 1e-10) %>% 
  group_by(node_id) %>% 
  slice_head(n = 1) %>% 
  ungroup()

marg_mean_map <- ggplot() +
  geom_point(data = sample_node_df %>% 
               filter(node_id != 1L),
             aes(x = x, y = y),
             size = 4, color = '#999999') +
  geom_point(data = sample_node_df %>% 
               filter(node_id == 1L),
             aes(x = x, y = y),
             size = 8, color = 'black') +
  geom_point(data = conditional_mean_df,
             aes(x = x, y = y, color = tree),
             shape = 8, size = 8.5) +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = '#e5e5e5'),
        panel.grid = element_blank(),
        axis.title = element_text(size = 14)) +
  scale_color_manual(values = span_color) +
  scale_y_reverse()

ggsave(here('plots', 'conceptual_figs', 'non_admix', 'marg_mean_map.png'),
       marg_mean_map,
       width = 5.3, height = 5)



#################################
### CODE NOT CURRENTLY IN USE ###
#################################

# nodes <- data.frame(
#   id = c("A","B","C","D"),
#   x = c(1,2,3,2),
#   y = c(3,3,2,1),
#   time = c(3,3,2,1)
# )
 
# # edges
# edges <- data.frame(
#   parent = c("A","B"),
#   child  = c("C","D")
# )

# edges_plot <- merge(edges, nodes, by.x="parent", by.y="id")
# edges_plot <- merge(edges_plot, nodes, by.x="child", by.y="id",
#                     suffixes = c("_p", "_c"))


# ggplot(nodes, aes(x, y, color = time)) +
#   geom_point() +
#   geom_segment(data = edges_plot,
#                aes(x = x_p, y = y_p,
#                    xend = x_c, yend = y_c),
#                inherit.aes = FALSE) +
#   scale_color_viridis_c() +
#   coord_fixed()


#library(plotly)


# fig <- plot_ly()
# 
# # add edges as line segments
# for(i in 1:nrow(edges_plot)) {
#   fig <- fig %>% add_trace(
#     type = "scatter3d",
#     mode = "lines",
#     x = c(edges_plot$x_p[i], edges_plot$x_c[i]),
#     y = c(edges_plot$y_p[i], edges_plot$y_c[i]),
#     z = c(edges_plot$time_p[i], edges_plot$time_c[i]),
#     line = list(color = 'gray'),
#     showlegend = FALSE
#   )
# }

# # add nodes
# fig <- fig %>% add_trace(
#   type = "scatter3d",
#   mode = "markers+text",
#   x = nodes$x,
#   y = nodes$y,
#   z = nodes$time,
#   text = nodes$id,
#   textposition = "top center",
#   marker = list(size = 4)
# )
#
# p <- plot_ly()
# p <- p %>%
#   add_trace(
#     data = nodes,
#     x = ~x, y = ~y, z = ~time,
#     type = "scatter3d",
#     mode = "markers",
#     marker = list(size = 5, color = "blue")
#   )
# 
# # Add line segments (each segment is drawn separately)
# for(i in 1:nrow(edge_df)) {
#   p <- p %>%
#     add_trace(
#       x = c(edge_df$parent_x[i], edge_df$child_x[i]),
#       y = c(edge_df$parent_y[i], edge_df$child_y[i]),
#       z = c(edge_df$parent_time[i], edge_df$child_time[i]),
#       type = "scatter3d",
#       mode = "lines",
#       line = list(width = 4, color = "blue"),
#       showlegend = FALSE
#     )
# }
#
#
# p %>%
#   layout(
#     scene = list(
#       xaxis = list(
#         showticklabels = FALSE,
#         showgrid = TRUE,     # keeps gridlines
#         zeroline = FALSE,
#         showbackground = TRUE,
#         title = '',
#         backgroundcolor = "rgba(230,230,230,0.5)"  # optional "floor"
#       ),
#       yaxis = list(
#         showticklabels = FALSE,
#         showgrid = TRUE,     # keeps gridlines
#         zeroline = FALSE,
#         showbackground = TRUE,
#         title = '',
#         #gridcolor = 'green',
#         backgroundcolor = "rgba(230,230,230,0.5)"  # optional "floor"
#       ),
#       zaxis = list(showgrid = FALSE,
#                    showticklabels = FALSE,
#                    title = '',
#                    label = FALSE # This changes the floor (x-y plane)
#       )
#     )
#   )
#
#
# node_df <- data.frame(node_time = c(0, 1, 2, 2, 2),
#            node_x = c(5, 6, 7, 9, 2),
#            node_y = c(6, 3, 7, 6.5, 8))

# line_df <- data.frame(
#   z_start,
#   z_end,
#   x_start,
#   x_end,
#   y_start
#   
# )
# node_df %>% 
#   ggplot() +
#   geom_point(aes(x = node_x, y = node_y)) +
#   xlim(0, 10) +
#   ylim(0, 10)

# 
# 
# fig <- plot_ly()
# 
# # Add points
# fig <- fig %>%
#   add_trace(
#     data = node_df,
#     x = ~node_x, y = ~node_y, z = ~node_time,
#     type = "scatter3d",
#     mode = "markers",
#     marker = list(size = 5, color = "blue")
#   )

# # df <- data.frame(
# #   id = 1:5,
# #   x = c(1, 2, 3, 4, 5),
# #   y = c(1, 3, 2, 5, 4),
# #   z = c(2, 1, 4, 3, 5)
# # )
# # 
# # # Define which points should be connected
# # # Each row is a segment: from (x,y,z) to (xend,yend,zend)
# # edges <- data.frame(
# #   x    = df$x[c(1, 2, 3)],
# #   y    = df$y[c(1, 2, 3)],
# #   z    = df$z[c(1, 2, 3)],
# #   xend = df$x[c(2, 3, 5)],
# #   yend = df$y[c(2, 3, 5)],
# #   zend = df$z[c(2, 3, 5)]
# # )
# 
# # Create plot
# p <- plot_ly()
# 
# # Add points
# p <- p %>%
#   add_trace(
#     data = df,
#     x = ~x, y = ~y, z = ~z,
#     type = "scatter3d",
#     mode = "markers",
#     marker = list(size = 5, color = "blue")
#   )
# 
# # Add line segments (each segment is drawn separately)
# for(i in 1:nrow(edges)) {
#   p <- p %>%
#     add_trace(
#       x = c(edges$x[i], edges$xend[i]),
#       y = c(edges$y[i], edges$yend[i]),
#       z = c(edges$z[i], edges$zend[i]),
#       type = "scatter3d",
#       mode = "lines",
#       line = list(width = 4, color = "gray"),
#       showlegend = FALSE
#     )
# }
# 
# 
# p %>%
#   layout(
#     scene = list(
#       xaxis = list(
#         showticklabels = FALSE,
#         showgrid = TRUE,     # keeps gridlines
#         zeroline = FALSE,
#         showbackground = TRUE,
#         backgroundcolor = "rgba(230,230,230,0.5)"  # optional "floor"
#       ),
#       yaxis = list(
#         showticklabels = FALSE,
#         showgrid = TRUE,     # keeps gridlines
#         zeroline = FALSE,
#         showbackground = TRUE,
#         backgroundcolor = "rgba(230,230,230,0.5)"  # optional "floor"
#       ),
#       zaxis = list(showgrid = FALSE,
#                    showticklabels = FALSE,
#                    label = FALSE # This changes the floor (x-y plane)
#     )
#   )
# )
