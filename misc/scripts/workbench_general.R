synthetic_gex <- synthetic_signal_matrix()

plot(synthetic_gex)

cpca_synthetic_data <- synthetic_c_pca_data()

plot(cpca_synthetic_data)

cpca_data <- cpca_synthetic_data

synthetic_gex$mat[1:5, 1:5]

synthetic_gex$group


plot_df <- as.data.table(synthetic_gex$mat, keep.rownames = "genes") %>%
  melt(id.vars = c("genes")) %>%
  merge(
    .,
    data.table::data.table(
      variable = names(synthetic_gex$group),
      grp = synthetic_gex$group
    ),
    by = "variable"
  )


ggplot(
  data = plot_df,
  mapping = aes(y = variable, x = genes, fill = value)
) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  ggtitle("Target data") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  ) +
  xlab("Features") +
  ylab("Samples") +
  geom_rug(
    data = unique(plot_df[, c("variable", "grp")]),
    mapping = aes(y = variable, color = grp),
    inherit.aes = FALSE,
    sides = "l",
    linewidth = 1.25,
    length = unit(0.0125, "npc")
  ) +
  labs(fill = "Value:", colour = "Group:")
