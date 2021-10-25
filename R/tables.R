
############
# Table 1 # - list of full category names, abbreviations and which terms have been grouped
############

# Get pathway_group names in final data frame and ADD nitrogen_fixation (for Other nitrogen fixation)
our_terms <- pw %>% filter(group_name %in% c(unique(df$group_name),"Nitrogen fix."))

#Get original list of pathway translations (to go all the way back to original terms)
conv <- read.csv("data/Madin et al.2020/renaming_pathways.csv")

#Get list of translated terms within selected pathway groups (existing in table of selected pathways 'pw')
pt <- our_terms %>% left_join(conv, by = c("pathway"="New"))

#Create table with each pathway and a list of old terms
pt2 <- pt %>% group_by(full_group_name) %>%
  mutate(original_terms = paste0(Original, collapse = ", ")) %>% 
  ungroup() %>%
  distinct(full_group_name, .keep_all = TRUE) %>% 
  select(full_group_name, group_abbreviation, original_terms) 

# Add description for cyanobacterial nitrogen fixation and Oxygenic photoautotrophy
# since these were created manually
pt2$original_terms[pt2$full_group_name == "Oxygenic photoautotrophy"] <- "Manually created: all phylum = Cyanobacteria"
pt2$original_terms[pt2$full_group_name == "Cyanobacterial nitrogen fixation"] <- "Manually created: all phylum = Cyanobacteria & pathway = nitrogen fixation"

#Move terms from general nitrogen fixation to 'other' nitrogen fixation
descr <- "excluding cyanobacteria"
pt2$original_terms[pt2$full_group_name == "Other nitrogen fixation"] <- paste(pt2$original_terms[pt2$full_group_name == "Nitrogen fixation"], descr, sep = "; ")

# Remove "nitrogen_fixation" from this list since this has been converted to "Other nitrogen fixation"
pt2 <- pt2 %>% filter(!(full_group_name == "Nitrogen fixation"))

# Rename columns
pt2 <- pt2 %>% rename("Full category name" = full_group_name, "Abbreviations" = group_abbreviation, "Terms included from data sources" = original_terms)

#Save
write.csv(pt2, "output/Table1.csv", row.names = FALSE)



#######
# Table 2: Summary of standard deviations between and within substrate-use groups
#######

#Transform data
dat <- df %>%
  select(group_name,d1_mid,genome_size,temp_adjusted_maxgrowth,rRNA16S_genes,growth_tmp,hk_tot,stp) %>%
  mutate(d1_mid = log10(d1_mid),
         genome_size  = log10(genome_size),
         rRNA16S_genes = log10(rRNA16S_genes),
         hk_tot = sqrt(hk_tot),
         stp = sqrt(stp))


# Gather data to enable simple calculations
dat2 <- dat %>% gather(key = "trait", value = "value", 2:8, na.rm = TRUE)

# Create table of means and SDs within each group and trait
within <- dat2 %>% group_by(group_name, trait) %>% 
  summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE))

# Calculate mean values of SDs for each trait
within_sds <- within %>% group_by(trait) %>% 
  summarise(sd_within = mean(sd, na.rm = TRUE)) %>% 
  mutate(sd_within = signif(sd_within, 3))

# Calculate SDS of across groups (on means of within)
across_sds <- within %>% group_by(trait) %>% 
  summarise(sd_across = sd(mean, na.rm = TRUE)) %>% 
  mutate(sd_across = signif(sd_across, 3))

fin <- across_sds %>% left_join(within_sds, by= "trait") %>% 
  mutate("ratio" = sd_within/sd_across) %>% 
  mutate(ratio = signif(ratio,3)) %>% 
  rename(`Ratio of within to between SDs` = ratio,
         `SD across group means` = sd_across,
         `Mean of SDs within groups` = sd_within)

write.csv(fin,"output/Table2.csv", row.names = FALSE)


#######
# Table 3: Maximum percentages of variation (R2 as percent) in each of seven quantitative
# traits explainable by substrate use groups and by majority habitat groups, plus percentage of
# variation that is overlapped
 #######

df2 <- df %>% rename("resource" = group_name)

table3 <- data.frame("Trait" = as.character(), "Substrate use groups" = as.character(), "Habitat groups" = as.character(), "Joint or overlapped explanation" = as.character())

# Genome size
vars <- df2 %>% select(genome_size, habitat, resource) %>% 
  na.omit() %>%
  mutate("genome_size_log10" = log10(genome_size)) %>% 
  select(-genome_size) 
# Run the hierarchical partitioning like this
mod <- hier.part(vars$genome_size_log10, vars[c("resource", "habitat")])

sub <- as.data.frame(t(c("Genome size", as.numeric(signif(mod$IJ$Total[1]*100,3)), as.numeric(signif(mod$IJ$Total[2]*100,3)), as.numeric(signif(mod$IJ$J[1]*100,3)))))
names(sub) <- names(table3)
table3 <- table3 %>% bind_rows(sub)


# d1_mid
vars <- df2 %>% select(d1_mid, habitat, resource) %>% 
  na.omit() %>%
  mutate("d1_mid_log10" = log10(d1_mid)) %>% 
  select(-d1_mid)
# Run the hierarchical partitioning like this
mod <- hier.part(vars$d1_mid_log10, vars[c("resource", "habitat")])

sub <- as.data.frame(t(c("Diameter", as.numeric(signif(mod$IJ$Total[1]*100,3)), as.numeric(signif(mod$IJ$Total[2]*100,3)), as.numeric(signif(mod$IJ$J[1]*100,3)))))
names(sub) <- names(table3)
table3 <- table3 %>% bind_rows(sub)


# temp_adjusted_maxgrowth
vars <- df2 %>% select(temp_adjusted_maxgrowth, habitat, resource) %>% 
  na.omit() 
# Run the hierarchical partitioning like this
mod <- hier.part(vars$temp_adjusted_maxgrowth, vars[c("resource", "habitat")])

sub <- as.data.frame(t(c("Temp adjusted max growth rate", as.numeric(signif(mod$IJ$Total[1]*100,3)), as.numeric(signif(mod$IJ$Total[2]*100,3)), as.numeric(signif(mod$IJ$J[1]*100,3)))))
names(sub) <- names(table3)
table3 <- table3 %>% bind_rows(sub)


# rRNA16S_genes
vars <- df2 %>% select(rRNA16S_genes, habitat, resource) %>% 
  na.omit() %>%
  mutate("rRNA16S_genes_log10" = log10(rRNA16S_genes)) %>% 
  select(-rRNA16S_genes)
# Run the hierarchical partitioning like this
mod <- hier.part(vars$rRNA16S_genes_log10, vars[c("resource", "habitat")])

sub <- as.data.frame(t(c("16S rRNA genes", as.numeric(signif(mod$IJ$Total[1]*100,3)), as.numeric(signif(mod$IJ$Total[2]*100,3)), as.numeric(signif(mod$IJ$J[1]*100,3)))))
names(sub) <- names(table3)
table3 <- table3 %>% bind_rows(sub)


# growth_tmp
vars <- df2 %>% select(growth_tmp, habitat, resource) %>% 
  na.omit() %>%
  mutate("growth_tmp_log10" = log10(growth_tmp)) %>% 
  select(-growth_tmp)
# Run the hierarchical partitioning like this
mod <- hier.part(vars$growth_tmp_log10, vars[c("resource", "habitat")])

sub <- as.data.frame(t(c("Growth temperature", as.numeric(signif(mod$IJ$Total[1]*100,3)), as.numeric(signif(mod$IJ$Total[2]*100,3)), as.numeric(signif(mod$IJ$J[1]*100,3)))))
names(sub) <- names(table3)
table3 <- table3 %>% bind_rows(sub)


# hk_tot
vars <- df2 %>% select(hk_tot, habitat, resource) %>% 
  na.omit() %>%
  mutate("hk_tot_sqrt" = sqrt(hk_tot)) %>% 
  select(-hk_tot)
# Run the hierarchical partitioning like this
mod <- hier.part(vars$hk_tot_sqrt, vars[c("resource", "habitat")])

sub <- as.data.frame(t(c("Histidine kinases", as.numeric(signif(mod$IJ$Total[1]*100,3)), as.numeric(signif(mod$IJ$Total[2]*100,3)), as.numeric(signif(mod$IJ$J[1]*100,3)))))
names(sub) <- names(table3)
table3 <- table3 %>% bind_rows(sub)


# stp
vars <- df2 %>% select(stp, habitat, resource) %>% 
  na.omit() %>%
  mutate("stp_sqrt" = sqrt(stp)) %>% 
  select(-stp)
# Run the hierarchical partitioning like this
mod <- hier.part(vars$stp_sqrt, vars[c("resource", "habitat")])

sub <- as.data.frame(t(c("Signal transduction proteins", as.numeric(signif(mod$IJ$Total[1]*100,3)), as.numeric(signif(mod$IJ$Total[2]*100,3)), as.numeric(signif(mod$IJ$J[1]*100,3)))))
names(sub) <- names(table3)
table3 <- table3 %>% bind_rows(sub)

#Fix column names
names(table3) <- gsub("\\."," ", names(table3))

#Save table
write.csv(table3, "output/Table3.csv", row.names = FALSE)



############
# Table S1 # - Count of species in pathway / habitat combinations
############

tableS1 <- df %>% count(habitat, full_group_name) %>% 
  pivot_wider(names_from = habitat, values_from = n) 

tableS1[is.na(tableS1)] <- 0

tableS1$total <- rowSums(tableS1[,2:10], na.rm = TRUE)

write.csv(tableS1, "output/TableS1.csv", row.names = FALSE)



############
# Table S2 # - Organisms excluded from plots S1 and S2
############

# Get organisms excluded using the x-axis restrictions in plot S1
tableS2 <- df %>% filter(genome_size < 0.5 | genome_size > 18) %>% 
  select(genus, species, genome_size) %>% 
  rename("value" = genome_size) %>% 
  mutate(excluded_by = "genome_size", range = "0.5 - 18")

# Get organisms excluded using the x-axis restrictions in plotcS2
tableS2 <- df %>% filter(d1_mid < 0.1 | d1_mid > 10) %>% 
  select(genus, species, d1_mid) %>% 
  rename("value" = d1_mid) %>% 
  mutate(excluded_by = "d1_mid", range = "0.1 - 10") %>%
  bind_rows(tableS2) 

# Save list of excluded organisms
write.csv(tableS2, "output/TableS2.csv", row.names = FALSE)
