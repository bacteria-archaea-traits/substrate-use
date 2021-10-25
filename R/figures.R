
############
# Figure 1 #
############

# Figure combined from multiple plots generated in below loop. 
# Each plot represents a specific trait listed in 'plot_traits' below.

#Transform data
sub <- df %>%
  select(group_abbreviation,d1_mid,genome_size,temp_adjusted_maxgrowth,rRNA16S_genes,growth_tmp,hk_tot,stp) %>%
  mutate(d1_mid = log10(d1_mid),
         genome_size  = log10(genome_size),
         rRNA16S_genes = log10(rRNA16S_genes),
         hk_tot = sqrt(hk_tot),
         stp = sqrt(stp)) %>%
  mutate(group_abbreviation = as.factor(group_abbreviation))


# General defined features
plot_traits <- c("genome_size", "d1_mid", "temp_adjusted_maxgrowth", "rRNA16S_genes", "hk_tot", "stp", "growth_tmp")
plot_titles <- c("Log10 Genome size (Mbp)", "log10 Cell diameter", "Growth rate residuals", "log10 16S rRNA genes", "SQRT HK+HHK signalling", "SQRT signal trans. proteins", "Growth temperature (°C)")

# Define X,Y location of statistics for each plot
plot_letter_positions <- c(0,-0.8,-1.5,-0.2,0,0,0)
plot_count_positions <- c(0.08,-0.65,-1.3,-0.1,1,1.9,6.5)

# Crate loop to generate all files 
for(i in 1:length(plot_traits)) {
  
  trait <- plot_traits[i]
  title <- plot_titles[i]
  
  letter_position <- plot_letter_positions[i]
  count_position <- plot_count_positions[i]
  
  sub2 <- sub %>% filter(!is.na(get(trait))) %>% 
    group_by(group_abbreviation) %>%
    summarise(av = mean(get(trait)), sd = sd(get(trait), na.rm = TRUE), n = n()) %>% 
    ungroup()
  
  t <- sub %>% filter(!is.na(get(trait)))
  aov.test <- aov(get(trait)~group_abbreviation, t)
  t2 <- HSD.test(aov.test, trt = "group_abbreviation", group = TRUE, console = TRUE)$groups
  t2 <- cbind(rownames(t2), data.frame(t2, row.names=NULL))
  t2[,2] <- NULL
  names(t2)[1] <- "group_abbreviation"
  
  sub3 <- sub2 %>% inner_join(t2, by = "group_abbreviation") %>%
    mutate(group_abbreviation = fct_reorder(group_abbreviation, av, .desc = TRUE))
  
  overall_average <- mean(sub3$av)
  overall_sd <- sd(sub3$av)
  
  fig <- sub3 %>% filter(title == !!title) %>% 
    ggplot(aes(x = fct_reorder(group_abbreviation, av, .desc = TRUE), y = av)) +
    geom_hline(yintercept = overall_average, colour = "darkgrey", linetype = 2) +
    geom_ribbon(aes(ymin = overall_average - overall_sd, ymax = overall_average + overall_sd), fill = "grey", alpha = 0.5, group = 1) +
    geom_errorbar(aes(ymin = av - sd, ymax = av + sd), width = 0.3) +
    geom_point(aes(fill = group_abbreviation), shape = 21, colour = "black", size = 2.5, show.legend = FALSE) +
    scale_fill_manual(values = pathway_colours) +
    basic_layout +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )+
    geom_text(aes(label = groups, y = letter_position), size = 3, family = font_family) +
    geom_text(aes(label = n, y = count_position), size = 2.5, family = font_family) +
    labs(x = "", y = title)
  fig
  
  ggsave(filename = sprintf("Figure1_part_%s.tiff",trait), plot = fig, device = "tiff", path = save_path, units = "cm", width = 18, height = 6, dpi = 600, limitsize = TRUE) 
  
}


############
# Figure 2 #
############

# PCAs (with growth rate)

# Prepare data for PCAs
# Transform data according to data type
# Note: temp_adjusted_maxgrowth is the residuals of the log10 transformed growth rate 
# against growth temperature (see prep.R), and so should not be transformed again here.

org <- df %>%
  filter(!is.na(genome_size) & 
           !is.na(d1_mid) & 
           !is.na(rRNA16S_genes) & 
           !is.na(growth_tmp) & 
           !is.na(hk_tot) &
           !is.na(temp_adjusted_maxgrowth) & 
           !is.na(stp)) %>%
  mutate("Genome size" = log10(genome_size), 
         D1 = log10(d1_mid),
         RRN = sqrt(rRNA16S_genes),
         "HK tot" = sqrt(hk_tot), 
         "Growth tmp" = growth_tmp,
         "Max growth rate" = temp_adjusted_maxgrowth, 
         "Signal trans. proteins" = stp)

org <- org %>% select(group_abbreviation, species, `Genome size`, `HK tot`, `Growth tmp`, `Max growth rate`, D1, RRN, `Signal trans. proteins`)

# Define only substrate use groups with >= 7 species / rows, rest defined as 'Other'
dat <- org %>% group_by(group_abbreviation) %>% mutate(n = n()) %>% 
  ungroup() %>% 
  mutate(group_abbreviation = ifelse(n < 7, "Other", group_abbreviation)) %>% 
  select(-n)

# Set obs number to merge with pca later
dat$obs <-  seq(1,nrow(dat),1)

# Remove non-data columns for PCA
dat2 <- dat %>% select(-species, -group_abbreviation, -obs)

# Set rownames as number values
rownames(dat2) <- seq(1,nrow(dat2),1)

# Run PCA function
pca <- pca_format(dat2)

pca$obsnames <- as.numeric(as.character(pca$obsnames))

#Join back to original data via obs numbers to get group  names (more safe than just assuming same order)
pca <- pca %>% left_join(dat[,c("obs","group_abbreviation")], by = c("obsnames"="obs"))


# PC 1-2 #
##########

# Chose PC axes to output
plot_pcs <- c("PC1","PC2")

# Get % and vectors for specified PC axes (above)
props <- pca_prop_explained(dat2,plot_pcs)
vectors<- pca_vectors(dat2,plot_pcs)

# Plot data
# NOTE: To avoid overlap of data points geom_jitter is used for this plot
# This means the exact location of each dot will vary slightly between versions of this plot
a <- ggplot(pca, aes_string(x = plot_pcs[1], y = plot_pcs[2])) + 
  geom_jitter(aes(fill = group_abbreviation), shape = 21, colour = "black", alpha = 1, size = 2, width = 0.1, height = 0.1) +
  scale_x_continuous(limits = c(-5,5)) +
  scale_y_continuous(limits = c(-4,4)) +
  scale_fill_manual(values = pathway_colours_PCA) +
  basic_layout +
  theme(
    #legend.justification = c(0,1),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  ) +
  guides(fill = guide_legend(ncol = 3)) +
  labs(x = paste0(sprintf("%s (%s",plot_pcs[1], props[1]),"%)"), y = paste0(sprintf("%s (%s",plot_pcs[2], props[2]),"%)"),
       title = "", fill = "")

# Add PCA vectors
a <- a + 
  coord_equal() + 
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = v1, yend = v2), arrow = arrow(length = unit(0.1,"cm")), size = 0.3, alpha = 1, color = "black")

# Add vector names
a <- a + 
  geom_text(data = vectors, aes(x = v1*1.1, y = v2*1.1, label = varnames), size = 2, vjust = 0, hjust = 0, alpha = 1, colour = "black")


# PC 1-3 #
##########

# Chose PC axes to output
plot_pcs <- c("PC1","PC3")

# Get % and vectors for specified PC axes (above)
props <- pca_prop_explained(dat2,plot_pcs)
vectors<- pca_vectors(dat2,plot_pcs)

# Plot data
# NOTE: To avoid overlap of data points geom_jitter is used for this plot
# This means the exact location of each dot will vary slightly between versions of this plot
b <- ggplot(pca, aes_string(x = plot_pcs[1], y = plot_pcs[2])) + 
  geom_jitter(aes(fill = group_abbreviation), shape = 21, colour = "black", alpha = 1, size = 2, width = 0.1, height = 0.1) +
  scale_x_continuous(limits = c(-5,5)) +
  scale_y_continuous(limits = c(-4,4)) +
  scale_fill_manual(values = pathway_colours_PCA) +
  basic_layout +
  theme(
    #legend.justification = c(0,1),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  ) +
  guides(fill = guide_legend(ncol = 3)) +
  labs(x = paste0(sprintf("%s (%s",plot_pcs[1], props[1]),"%)"), y = paste0(sprintf("%s (%s",plot_pcs[2], props[2]),"%)"),
       title = "", fill = "")

# Add PCA vectors
b <- b + 
  coord_equal() + 
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = v1, yend = v2), arrow = arrow(length = unit(0.1,"cm")), size = 0.3, alpha = 1, color = "black")

# Add vector names
b <- b + 
  geom_text(data = vectors, aes(x = v1*1.1, y = v2*1.1, label = varnames), size = 2, vjust = 0, hjust = 0, alpha = 1, colour = "black")
b

#Combine both plots
# fig2_no_names <- ggarrange(a_no_names, b_no_names, ncol = 2, common.legend = TRUE, legend="bottom")
fig2 <- ggarrange(a, b, ncol = 2, common.legend = TRUE, legend="bottom")
  
# Save complete version
# ggsave(filename = "Figure2_clean.tiff", plot = fig2_no_names, device = "tiff", path = save_path, units = "cm", width = 18, height = 10, dpi = 600)
ggsave(filename = "Figure2.tiff", plot = fig2, device = "tiff", path = save_path, units = "cm", width = 18, height = 10, dpi = 600)



############
# Figure 3 #
############

# Phylogenetic analysis of relationships within resource groups

# NOTE: The output from below loop is also used for Figure 4.

# Count number of pairs of species within each phylogenetic level of each resource group
substrate_group <- unique(df$group_abbreviation)
levels <- c("genus","family","order","class","phylum","superkingdom")

results <- data.frame("substrate" = as.character(), "level" = as.character(), "pairs" = as.numeric(), "fraction" = as.numeric(), "weighted_sum" = as.numeric())

for(i in 1:length(substrate_group)) { 
  #For each substrate group
  
  this_substrate <- substrate_group[i]
  
  sub <- df %>% filter(group_abbreviation == this_substrate) %>% 
    select(group_abbreviation, species, genus, family, order, class, phylum, superkingdom) %>% 
    drop_na()
  
  already_counted <- 0
  
  total_pairs_possible <- nrow(sub)*(nrow(sub)-1)/2
  
  dat <- results[0,]
  
  for(a in 1:length(levels)) {
    #For each phylogenetic level count number of relationships 
  
    sub2 <- sub %>% select_at(levels[a]) %>% 
      na.omit() %>%
      group_by_at(levels[a]) %>% 
      summarise(n = n(), pairs = n*(n-1)/2) 
    
    this_sum <- sum(sub2$pairs)
    in_level <- this_sum - already_counted
    already_counted <- already_counted + in_level
    fraction <- in_level / total_pairs_possible
    
    if(levels[a] != "genus") {
      #Weighted sum is fraction * weight + sum of previous level
      weighted_sum <- fraction * a + dat$weighted_sum[nrow(dat)]
    } else {
      #Weighted sum for genus is just the fraction
      weighted_sum <- fraction
    }

    tmp <- data.frame("substrate" = as.character(this_substrate), "level" = as.character(levels[a]), "pairs" = as.numeric(in_level), "fraction" = as.numeric(fraction), "weighted_sum" = as.numeric(weighted_sum))
    #bind to data table for this substrate
    dat <- dat %>% bind_rows(tmp)

    #Add between superkingdom pairs to reach total pairs (=n_archaea * n_bacteria)
    if(levels[a] == "superkingdom") {
      
      a <- a+1 #a = 7 for betwe-superkingdom weighting
      
      if(length(sub2$superkingdom) == 2) {
        #There is more than one superkingdom
        in_level <- sub2$n[1]*sub2$n[2]
      } else {
        in_level <- 0
      }
      
      fraction <- in_level/total_pairs_possible
      weighted_sum <- fraction * a + dat$weighted_sum[nrow(dat)]
        
      tmp <- data.frame("substrate" = as.character(this_substrate), "level" = as.character("betw-superkingdom"), "pairs" = as.numeric(in_level), "fraction" = as.numeric(fraction), "weighted_sum" = as.numeric(weighted_sum))
      #bind to data table for this substrate
      dat <- dat %>% bind_rows(tmp)
      
    }
    
  }
  
  #Add this substrate results to main data output
  results <- results %>% bind_rows(dat)
}

results$pairs <- as.numeric(results$pairs)
results$fraction <- as.numeric(results$fraction)

#Add fake data for placeholders
row <- results[1:3,]
row[,] <- 0
row$substrate[1] <- "placeholder1"
row$substrate[2] <- "placeholder2"
row$substrate[3] <- "placeholder3"
row$level <- "genus"
results <- results %>% bind_rows(row)

results <- results %>% mutate(fraction = as.numeric(fraction))%>% 
  mutate(level = factor(level, levels = c("genus","family","order","class","phylum","superkingdom","betw-superkingdom")))

# Create histogram of phylogenetic distribution 
results$substrate <- factor(results$substrate, levels = plot_order)

# Create data frame with max weighted sums
dat_text <- results %>% filter(level == "betw-superkingdom")

p <- results %>% 
  ggplot(aes(x = level, y = fraction)) + 
  geom_col(fill = "#22638F", alpha = 1) + 
  scale_y_continuous(limits = c(0,1)) +
  basic_layout +
  theme(
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  geom_text(data = dat_text, mapping = aes(label = sprintf("ws = %s", signif(weighted_sum,3)), x = 2.5, y = 0.9), size = 2.5) +
  facet_wrap(~substrate, ncol = 6, nrow = 4, drop = FALSE) + 
  labs(x = "", y = "Scaled proportion")


# Remove placeholder plots

#panel-2-2
#panel-6-2
#panel-6-3

g <- ggplotGrob(p)
# get the grobs that must be removed
rm_grobs <- g$layout$name %in% c("panel-2-2", "panel-6-2", "panel-6-3", "strip-t-2-4", "strip-t-3-4", "strip-t-5-2")
# remove grobs
g$grobs[rm_grobs] <- NULL
g$layout <- g$layout[!rm_grobs, ]

## move axis closer to panel
g$layout[g$layout$name == "axis-b-2-4", c("t", "b")] <- c(19, 19)
g$layout[g$layout$name == "axis-b-3-4", c("t", "b")] <- c(19, 19)

grid.newpage()
grid.draw(g)

# Save figure
ggsave(filename = "Figure3.tiff", plot = g, device = "tiff", path = save_path, units = "cm", width = 18, height = 12, dpi = 300) 



#############
# Figure 4 #
#############

# NOTE: This figure is generated based on values generated under 'Figure 3'
# The above part of the code will have to be run first to generate this plot

#Get SDs for each trait across substrate groups
sds <- df %>%
  select(group_abbreviation,d1_mid,genome_size,temp_adjusted_maxgrowth,rRNA16S_genes,growth_tmp,hk_tot,stp) %>%
  mutate(d1_mid = log10(d1_mid),
         genome_size  = log10(genome_size),
         rRNA16S_genes = log10(rRNA16S_genes),
         hk_tot = sqrt(hk_tot),
         stp = sqrt(stp)) %>%
  mutate(group_abbreviation = as.factor(group_abbreviation))%>%
  rename_at(vars(plot_traits ), ~ plot_titles) %>%
  pivot_longer(2:8, names_to = "trait",values_to = "value") %>% 
  na.omit() %>%
  group_by(group_abbreviation,trait) %>% 
  summarise(sd = sd(value)) %>% 
  ungroup()

# Combine with maximum weighted sums (values from figure 3)
sds <- sds %>% left_join(results[results$level == "betw-superkingdom",c("substrate","weighted_sum")], by = c("group_abbreviation"="substrate"))

# Create table of P values from correlation analysis
cors_val <- data.frame("trait" = as.character(), "cor" = as.character(), "p-value" = as.character())
traits <- unique(sds$trait)  

for(i in 1:length(traits)) {
  
  this_trait <- traits[i]
  sub <- sds %>% filter(trait == this_trait)
  
  test <- cor.test(sub$weighted_sum, sub$sd, method = "pearson")
  row <- c(this_trait,signif(test$estimate,3),signif(test$p.value,3))
  names(row) <- names(cors_val)
  cors_val <- cors_val %>% bind_rows(row)
}

# Get max sd per trait for positioning of correlation text in figures
max <- sds %>% group_by(trait) %>% 
  filter(sd == max(sd)) %>% select(trait, sd)

dat_text <- cors_val %>% select(-p.value) %>% 
  left_join(max, by = "trait")

p <- sds %>% ggplot(aes(x = weighted_sum, y = sd)) + 
  geom_point(aes(fill = group_abbreviation), shape = 21, colour = "black", size = 2.5, show.legend = TRUE) +
  scale_fill_manual(values = pathway_colours) +
  scale_x_continuous(limits=c(3,6)) +
  basic_layout +
  theme(
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.title = element_blank()
    
  ) +
  geom_text(data = cors_val, mapping = aes(label = sprintf("cor = %s", cor), x = 3.7, y = 0.90*dat_text$sd), size = 2.5) +
  facet_wrap(~trait, scales = "free", ncol = 4, nrow = 2) + 
  labs(x = "Weighted average phylogenetic separation WS", y = "SD within substrate group")
p

ggsave(filename = "Figure4.tiff", plot = p, device = "tiff", path = save_path, units = "cm", width = 18, height = 12, dpi = 200, limitsize = TRUE) 


#########################
#########################
# Supplementary figures #
#########################
#########################

# Traits to include
traits <- c("genome_size","d1_mid","temp_adjusted_maxgrowth","rRNA16S_genes","growth_tmp","hk_tot","tcp_tot","ocp","stp")

df2 <- df 

#Add fake data for placeholders
row <- sub[1:3,]
row[,] <- NA
row$group_abbreviation[1] <- "placeholder1"
row$group_abbreviation[2] <- "placeholder2"
row$group_abbreviation[3] <- "placeholder3"

row$group_abbreviation[1] <- "placeholder1"
row$group_abbreviation[2] <- "placeholder2"
row$group_abbreviation[3] <- "placeholder3"

df2 <- df2 %>% bind_rows(row)
df2$group_abbreviation <- factor(df2$group_abbreviation, levels = plot_order)


############# 
# Figure S1 # genome size
#############

trait <- "genome_size"
title <- "Genome size (Mbp)"

sub <- df2 %>% filter(!is.na(.data[[trait]]))
sub_non <- full_df %>% filter(!is.na(!is.na(.data[[trait]])))

# Create data frame with total count within each group
dat_text <- sub %>% group_by(group_abbreviation) %>% summarise(n = n())

p <- sub %>% 
  ggplot(aes(x = .data[[trait]])) + 
  geom_histogram(data = sub_non, bins = 30, aes(x = .data[[trait]], y= ..ncount..), fill = "gray", alpha = 1) +
  geom_histogram(bins = 30, aes(x = .data[[trait]], y= ..ncount..), fill = "blue", alpha = 0.5) +
  scale_y_continuous() +
  scale_x_log10(limits = c(0.5,18)) +
  basic_layout +
  theme(
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7)
  ) +
  geom_text(data = dat_text, mapping = aes(label = sprintf("n = %s",n), x = 0.9, y = 0.9), size = 2.5) +
  facet_wrap(~group_abbreviation, ncol = 6, nrow = 4, drop = FALSE) + 
  labs(x = title, y = "Scaled proportion")
p

# Remove placeholder plots
g <- ggplotGrob(p)
# get the grobs that must be removed
rm_grobs <- g$layout$name %in% c("panel-2-2", "panel-6-2", "panel-6-3", "strip-t-2-4", "strip-t-3-4", "strip-t-5-2")
# remove grobs
g$grobs[rm_grobs] <- NULL
g$layout <- g$layout[!rm_grobs, ]

## move axis closer to panel
g$layout[g$layout$name == "axis-b-2-4", c("t", "b")] <- c(19, 19)
g$layout[g$layout$name == "axis-b-3-4", c("t", "b")] <- c(19, 19)

grid.newpage()
grid.draw(g)

ggsave(filename = "FigureS1.tiff", plot = g, device = "tiff", path = save_path, units = "cm", width = 18, height = 12, dpi = 200, limitsize = TRUE)
dev.off()


############# 
# Figure S2 # Cell diameter
#############

trait <- "d1_mid"
title <- diameter_axis_text

sub <- df2 %>% filter(!is.na(.data[[trait]]))
sub_non <- full_df %>% filter(!is.na(!is.na(.data[[trait]])))

# Create data frame with total count within each group
dat_text <- sub %>% group_by(group_abbreviation) %>% summarise(n = n())

p <- sub %>% 
  ggplot(aes(x = .data[[trait]])) + 
  geom_histogram(data = sub_non, bins = 25, aes(x = .data[[trait]], y= ..ncount..), fill = "gray", alpha = 1) +
  geom_histogram(bins = 25, aes(x = .data[[trait]], y= ..ncount..), fill = "blue", alpha = 0.5) +
  scale_y_continuous() +
  scale_x_log10(limits = c(0.1, 10)) +
  basic_layout +
  theme(
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7)
  ) +
  geom_text(data = dat_text, mapping = aes(label = sprintf("n = %s",n), x = 5, y = 0.9), size = 2.5) +
  facet_wrap(~group_abbreviation, ncol = 6, nrow = 4, drop = FALSE) + 
  labs(x = title, y = "Scaled proportion")
p

g <- ggplotGrob(p)
# get the grobs that must be removed
rm_grobs <- g$layout$name %in% c("panel-2-2", "panel-6-2", "panel-6-3", "strip-t-2-4", "strip-t-3-4", "strip-t-5-2")
# remove grobs
g$grobs[rm_grobs] <- NULL
g$layout <- g$layout[!rm_grobs, ]

## move axis closer to panel
g$layout[g$layout$name == "axis-b-2-4", c("t", "b")] <- c(19, 19)
g$layout[g$layout$name == "axis-b-3-4", c("t", "b")] <- c(19, 19)

grid.newpage()
grid.draw(g)

ggsave(filename = "FigureS2.tiff", plot = g, device = "tiff", path = save_path, units = "cm", width = 18, height = 12, dpi = 200, limitsize = TRUE)
dev.off()


############# 
# Figure S3 # Growth rate residuals
#############

trait <- "temp_adjusted_maxgrowth"
title <- "Growth rate residuals"

sub <- df2 %>% filter(!is.na(.data[[trait]]))
sub_non <- full_df %>% filter(!is.na(!is.na(.data[[trait]])))

# Create data frame with total count within each group
dat_text <- sub %>% group_by(group_abbreviation) %>% summarise(n = n())

p <- sub %>% 
  ggplot(aes(x = .data[[trait]])) + 
  geom_histogram(data = sub_non, bins = 20, aes(x = .data[[trait]], y= ..ncount..), fill = "gray", alpha = 1) +
  geom_histogram(bins = 20, aes(x = .data[[trait]], y= ..ncount..), fill = "blue", alpha = 0.5) +
  scale_y_continuous() +
  basic_layout +
  theme(
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7)
  ) +
  geom_text(data = dat_text, mapping = aes(label = sprintf("n = %s",n), x = 1.2, y = 0.9), size = 2.5) +
  facet_wrap(~group_abbreviation, ncol = 6, nrow = 4, drop = FALSE) + 
  labs(x = title, y = "Scaled proportion")
p


g <- ggplotGrob(p)
# get the grobs that must be removed
rm_grobs <- g$layout$name %in% c("panel-2-2", "panel-6-2", "panel-6-3", "strip-t-2-4", "strip-t-3-4", "strip-t-5-2")
# remove grobs
g$grobs[rm_grobs] <- NULL
g$layout <- g$layout[!rm_grobs, ]

## move axis closer to panel
g$layout[g$layout$name == "axis-b-2-4", c("t", "b")] <- c(19, 19)
g$layout[g$layout$name == "axis-b-3-4", c("t", "b")] <- c(19, 19)

grid.newpage()
grid.draw(g)

ggsave(filename = "FigureS3.tiff", plot = g, device = "tiff", path = save_path, units = "cm", width = 18, height = 12, dpi = 200, limitsize = TRUE)
dev.off()


############# 
# Figure S4 # 16S rRNA genes
#############

trait <- "rRNA16S_genes"
title <- "16S rRNA genes"

sub <- df2 %>% filter(!is.na(.data[[trait]]))
sub_non <- full_df %>% filter(!is.na(!is.na(.data[[trait]])))

# Create data frame with total count within each group
dat_text <- sub %>% group_by(group_abbreviation) %>% summarise(n = n())

p <- sub %>% 
  ggplot(aes(x = .data[[trait]])) + 
  geom_histogram(data = sub_non, bins = 10, aes(x = .data[[trait]], y= ..ncount..), fill = "gray", alpha = 1) +
  geom_histogram(bins = 10, aes(x = .data[[trait]], y= ..ncount..), fill = "blue", alpha = 0.5) +
  scale_y_continuous() +
  scale_x_log10() +
  basic_layout +
  theme(
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7)
  ) +
  geom_text(data = dat_text, mapping = aes(label = sprintf("n = %s",n), x = 12, y = 0.9), size = 2.5) +
  facet_wrap(~group_abbreviation, ncol = 6, nrow = 4, drop = FALSE) + 
  labs(x = title, y = "Scaled proportion")
p


g <- ggplotGrob(p)
# get the grobs that must be removed
rm_grobs <- g$layout$name %in% c("panel-2-2", "panel-6-2", "panel-6-3", "strip-t-2-4", "strip-t-3-4", "strip-t-5-2")
# remove grobs
g$grobs[rm_grobs] <- NULL
g$layout <- g$layout[!rm_grobs, ]

## move axis closer to panel
g$layout[g$layout$name == "axis-b-2-4", c("t", "b")] <- c(19, 19)
g$layout[g$layout$name == "axis-b-3-4", c("t", "b")] <- c(19, 19)

grid.newpage()
grid.draw(g)

ggsave(filename = "FigureS4.tiff", plot = g, device = "tiff", path = save_path, units = "cm", width = 18, height = 12, dpi = 200, limitsize = TRUE)
dev.off()


############# 
# Figure S5 # Hk+HKK signalling
#############

trait <- "hk_tot"
title <- "HK+HHK signalling"

sub <- df2 %>% filter(!is.na(.data[[trait]]))
sub_non <- full_df %>% filter(!is.na(!is.na(.data[[trait]])))

# Create data frame with total count within each group
dat_text <- sub %>% group_by(group_abbreviation) %>% summarise(n = n())

p <- sub %>% 
  ggplot(aes(x = .data[[trait]])) + 
  geom_histogram(data = sub_non, bins = 20, aes(x = .data[[trait]], y= ..ncount..), fill = "gray", alpha = 1) +
  geom_histogram(bins = 20, aes(x = .data[[trait]], y= ..ncount..), fill = "blue", alpha = 0.5) +
  scale_y_continuous() +
  scale_x_sqrt(breaks = c(10,50,100,200,300)) +
  basic_layout +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7)
  ) +
  geom_text(data = dat_text, mapping = aes(label = sprintf("n = %s",n), x = 200, y = 0.9), size = 2.5) +
  facet_wrap(~group_abbreviation, ncol = 6, nrow = 4, drop = FALSE) + 
  labs(x = title, y = "Scaled proportion")
p


g <- ggplotGrob(p)
# get the grobs that must be removed
rm_grobs <- g$layout$name %in% c("panel-2-2", "panel-6-2", "panel-6-3", "strip-t-2-4", "strip-t-3-4", "strip-t-5-2")
# remove grobs
g$grobs[rm_grobs] <- NULL
g$layout <- g$layout[!rm_grobs, ]

## move axis closer to panel
g$layout[g$layout$name == "axis-b-2-4", c("t", "b")] <- c(19, 19)
g$layout[g$layout$name == "axis-b-3-4", c("t", "b")] <- c(19, 19)

grid.newpage()
grid.draw(g)

ggsave(filename = "FigureS5.tiff", plot = g, device = "tiff", path = save_path, units = "cm", width = 18, height = 12, dpi = 200, limitsize = TRUE)
dev.off()


############# 
# Figure S6 # Signal transduction proteins
#############

trait <- "stp"
title <- "Signal transduction proteins"

sub <- df2 %>% filter(!is.na(.data[[trait]]))
sub_non <- full_df %>% filter(!is.na(!is.na(.data[[trait]])))

# Create data frame with total count within each group
dat_text <- sub %>% group_by(group_abbreviation) %>% summarise(n = n())

p <- sub %>% 
  ggplot(aes(x = .data[[trait]])) + 
  geom_histogram(data = sub_non, bins = 25, aes(x = .data[[trait]], y= ..ncount..), fill = "gray", alpha = 1) +
  geom_histogram(bins = 25, aes(x = .data[[trait]], y= ..ncount..), fill = "blue", alpha = 0.5) +
  scale_y_continuous() +
  scale_x_sqrt(breaks = c(50,250,500,1000,1500)) +
  basic_layout +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7)
  ) +
  geom_text(data = dat_text, mapping = aes(label = sprintf("n = %s",n), x = 1000, y = 0.9), size = 2.5) +
  facet_wrap(~group_abbreviation, ncol = 6, nrow = 4, drop = FALSE) + 
  labs(x = title, y = "Scaled proportion")
p


g <- ggplotGrob(p)
# get the grobs that must be removed
rm_grobs <- g$layout$name %in% c("panel-2-2", "panel-6-2", "panel-6-3", "strip-t-2-4", "strip-t-3-4", "strip-t-5-2")
# remove grobs
g$grobs[rm_grobs] <- NULL
g$layout <- g$layout[!rm_grobs, ]

## move axis closer to panel
g$layout[g$layout$name == "axis-b-2-4", c("t", "b")] <- c(19, 19)
g$layout[g$layout$name == "axis-b-3-4", c("t", "b")] <- c(19, 19)

grid.newpage()
grid.draw(g)

ggsave(filename = "FigureS6.tiff", plot = g, device = "tiff", path = save_path, units = "cm", width = 18, height = 12, dpi = 200, limitsize = TRUE)
dev.off()


############# 
# Figure S7 # Growth temperature
#############

trait <- "growth_tmp"
title <- "Growth temperature (C)"

sub <- df2 %>% filter(!is.na(.data[[trait]]))
sub_non <- full_df %>% filter(!is.na(!is.na(.data[[trait]])))

# Create data frame with total count within each group
dat_text <- sub %>% group_by(group_abbreviation) %>% summarise(n = n())

p <- sub %>% 
  ggplot(aes(x = .data[[trait]])) + 
  geom_histogram(data = sub_non, bins = 25, aes(x = .data[[trait]], y= ..ncount..), fill = "gray", alpha = 1) +
  geom_histogram(bins = 25, aes(x = .data[[trait]], y= ..ncount..), fill = "blue", alpha = 0.5) +
  scale_y_continuous() +
  basic_layout +
  theme(
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7)
  ) +
  geom_text(data = dat_text, mapping = aes(label = sprintf("n = %s",n), x = 90, y = 0.9), size = 2.5) +
  facet_wrap(~group_abbreviation, ncol = 6, nrow = 4, drop = FALSE) + 
  labs(x = title, y = "Scaled proportion")
p


g <- ggplotGrob(p)
# get the grobs that must be removed
rm_grobs <- g$layout$name %in% c("panel-2-2", "panel-6-2", "panel-6-3", "strip-t-2-4", "strip-t-3-4", "strip-t-5-2")
# remove grobs
g$grobs[rm_grobs] <- NULL
g$layout <- g$layout[!rm_grobs, ]

## move axis closer to panel
g$layout[g$layout$name == "axis-b-2-4", c("t", "b")] <- c(19, 19)
g$layout[g$layout$name == "axis-b-3-4", c("t", "b")] <- c(19, 19)

grid.newpage()
grid.draw(g)

ggsave(filename = "FigureS7.tiff", plot = g, device = "tiff", path = save_path, units = "cm", width = 18, height = 12, dpi = 200, limitsize = TRUE)
dev.off()


#############
# Figure S8 # - PCA (excluding growth rate)
#############

# Prepare data for PCAs
# Transform data according to data type
org <- df %>%
  filter(!is.na(genome_size) & 
           !is.na(d1_mid) & 
           !is.na(rRNA16S_genes) & 
           !is.na(growth_tmp) & 
           !is.na(hk_tot) & 
           !is.na(stp)) %>%
  mutate("Genome size" = log10(genome_size), 
         D1 = log10(d1_mid),
         RRN = sqrt(rRNA16S_genes),
         "HK tot" = sqrt(hk_tot), 
         "Growth tmp" = growth_tmp,
         "Signal trans. proteins" = stp)

org <- org %>% select(group_abbreviation, species, `Genome size`, `HK tot`, `Growth tmp`, D1, RRN, `Signal trans. proteins`)

# Define only substrate use groups with >= 7 species / rows, rest defined as 'Other'
dat <- org %>% group_by(group_abbreviation) %>% mutate(n = n()) %>% 
  ungroup() %>% 
  mutate(group_abbreviation = ifelse(n < 7, "Other", group_abbreviation)) %>% 
  select(-n)

# Set obs number to merge with pca later
dat$obs <-  seq(1,nrow(dat),1)

dat2 <- dat %>% select(-species, -group_abbreviation, -obs)

# Set rownames as number values
rownames(dat2) <- seq(1,nrow(dat2),1)

# Run PCA function to generate data for plots
pca <- pca_format(dat2)

pca$obsnames <- as.numeric(as.character(pca$obsnames))

#Join back to original data via obs numbers to get group  names (more safe than just assuming same order)
pca <- pca %>% left_join(dat[,c("obs","group_abbreviation")], by = c("obsnames"="obs"))

# PC 1-2 #
##########

# Chose PC axes to output
plot_pcs <- c("PC1","PC2")

# Get % and vectors for specified PC axes (above)
props <- pca_prop_explained(dat2,plot_pcs)
vectors<- pca_vectors(dat2,plot_pcs)

# Plot data
a <- ggplot(pca, aes_string(x = plot_pcs[1], y = plot_pcs[2])) + 
  geom_jitter(aes(fill = group_abbreviation), shape = 21, colour = "black", alpha = 1, size = 1.5, width = 0.1, height = 0.1, stroke = 0.5) +
  scale_x_continuous(limits = c(-6.5,6.5)) +
  scale_y_continuous(limits = c(-4.5,4.5)) +
  scale_fill_manual(values = pathway_colours_PCA) +
  basic_layout +
  theme(
    #legend.justification = c(0,1),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  ) +
  guides(fill = guide_legend(ncol = 3)) +
  labs(x = paste0(sprintf("%s (%s",plot_pcs[1], props[1]),"%)"), y = paste0(sprintf("%s (%s",plot_pcs[2], props[2]),"%)"),
       title = "", fill = "")

# Add PCA vectors
a <- a + 
  coord_equal() + 
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = v1, yend = v2), arrow = arrow(length = unit(0.1,"cm")), size = 0.3, alpha = 1, color = "red")

# Add vector names
a <- a + 
  geom_text(data = vectors, aes(x = v1*1.1, y = v2*1.1, label = varnames), size = 2, vjust = 0, hjust = 0, alpha = 1, colour = "red")
a

# PC 1-3 #
##########

# Chose PC axes to output
plot_pcs <- c("PC1","PC3")

# Get % and vectors for specified PC axes (above)
props <- pca_prop_explained(dat2,plot_pcs)
vectors<- pca_vectors(dat2,plot_pcs)

# Plot data
b <- ggplot(pca, aes_string(x = plot_pcs[1], y = plot_pcs[2])) + 
  geom_jitter(aes(fill = group_abbreviation), shape = 21, colour = "black", alpha = 1, size = 1.5, width = 0.1, height = 0.1, stroke = 0.5) +
  scale_x_continuous(limits = c(-6.5,6.5)) +
  scale_y_continuous(limits = c(-4.5,4.5)) +
  scale_fill_manual(values = pathway_colours_PCA) +
  basic_layout +
  theme(
    #legend.justification = c(0,1),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  ) +
  guides(fill = guide_legend(ncol = 3)) +
  labs(x = paste0(sprintf("%s (%s",plot_pcs[1], props[1]),"%)"), y = paste0(sprintf("%s (%s",plot_pcs[2], props[2]),"%)"),
       title = "", fill = "")

# Add PCA vectors
b <- b + 
  coord_equal() + 
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = v1, yend = v2), arrow = arrow(length = unit(0.1,"cm")), size = 0.3, alpha = 1, color = "red")
b

# Add vector names
b <- b + 
  geom_text(data = vectors, aes(x = v1*1.1, y = v2*1.1, label = varnames), size = 2, vjust = 0, hjust = 0, alpha = 1, colour = "red")
b

#Combine both plots
figS8 <- ggarrange(a, b, ncol = 2, common.legend = TRUE, legend="bottom")

# Save complete version
ggsave(filename = "FigureS8.tiff", plot = figS8, device = "tiff", path = save_path, units = "cm", width = 18, height = 10, dpi = 600)
