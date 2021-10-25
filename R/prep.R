#########################
# Prepare main data frame
#########################

# Load main data file (from Madin et al. 2020)
df <- read.csv("data/Madin et al.2020/condensed_species_GTDB[NCBI_fill]_16102020.csv", as.is=TRUE)
# Ensure all rows have species designation
df <- df[!is.na(df$species),]


################################
# Define general figure layout #
################################

save_path <- "output/figures"

# Formats
text_size <- 11
text_color <- "black"
font_family <- "sans"
plot_line_color <- "black"
plot_line_width <- 0.5 #mm

# Define special axes
diameter_axis_text <- expression(paste("Cell diameter (",mu,m,")", sep = ""))

# Define theme standards 
basic_layout <- 
  theme_bw() + 
  theme(
    panel.border = element_rect(color = plot_line_color, size = plot_line_width),
    panel.grid = element_blank(),
    axis.title = element_text(size = text_size, family=font_family, color=text_color),
    axis.text = element_text(size = text_size, family=font_family, color=text_color),
    axis.line = element_line(color = plot_line_color, size = plot_line_width),
    axis.ticks = element_line(color = plot_line_color, size = plot_line_width),
    strip.background = element_blank(),
    strip.text = element_text(size = text_size)
  )


################################
# Add environment descriptions #
################################

en <- read.csv("data/Madin et al.2020/environments.csv", as.is=TRUE)

df <- df %>% left_join(en, by = c("isolation_source"="Type"))
rm(en)


######
# Merge with Mist transporter data
######

# Load mist data
mistdb <- read.csv("data/mistdb_extract_01022019.csv")

# Group transporters
df <- mistdb %>% group_by(species_tax_id) %>% 
  mutate(tcp_tot = sum(tcp.hk,tcp.hhk,tcp.rr,tcp.hrr,tcp.other, na.rm = TRUE), 
         hk_tot = sum(tcp.hk,tcp.hhk, na.rm = TRUE), 
         stp = sum(tcp_tot, ocp, na.rm = TRUE)) %>%
  select(species_tax_id,ocp, tcp_tot, stp, hk_tot, majormodes_total) %>% 
  right_join(df, by = "species_tax_id")

# Remove mistdb
rm(mistdb)


##########################
# Calculate mid diameter #
##########################

df$d1_mid <- ifelse(!is.na(df$d1_up), (df$d1_lo + df$d1_up)/2, df$d1_lo)
df$d2_mid <- ifelse(!is.na(df$d2_up), (df$d2_lo + df$d2_up)/2, df$d2_lo)

############################################
# Calculate growth rate from doubling time #
############################################

df$growth_rate <- log(2)/df$doubling_h

###############################
# Create new habitat category #
###############################

# First letter upper case for environment
df$Cobo.Simon.habitat <- str_to_title(df$Cobo.Simon.habitat, locale = "en")

#Change factor levels of habitat
#df$Cobo.Simon.habitat <- factor(df$Cobo.Simon.habitat, levels = c("Fresh","Marine","Soil","Rock","Host","Therm"))

# Use Cobo Simon habitat scheme
df$habitat <- as.character(df$Cobo.Simon.habitat)
# Set any isolation_source that does not fit into the Cobo Simon scheme to "Other"
df$habitat[!is.na(df$isolation_source) & is.na(df$habitat)] <- "Other"

# Regroup or Rename particular types of Cobo Simon habitats
df$habitat[df$habitat == "Rock"] <- "Other"
df$habitat[df$habitat == "Therm"] <- "Thermal"
df$habitat[grepl("plant|fungus|algae",df$isolation_source)] <- "Other"
df$habitat[grepl("feces|endotherm_surface|ectotherm_surface",df$isolation_source)] <- "Other"
df$habitat[df$isolation_source %in% c("host","host_animal")] <- "Other"

# Where habitat is not "Other" and isolation source is host-related, split into 'endo' and 'ecto'
df$habitat[grepl("endotherm",df$isolation_source) & !(df$habitat == "Other")] <- "Endotherm"
df$habitat[grepl("ectotherm",df$isolation_source) & !(df$habitat == "Other")] <- "Ectotherm"

# Set all species with no isolation_source/habitat information to habitat = "Other"
df$habitat[is.na(df$habitat)] <- "Other"

#Split fresh and marine into sediment and water column using isolation_source
df$habitat[!is.na(df$habitat) & df$habitat == "Fresh" & !is.na(df$isolation_source) & grepl('sediment',df$isolation_source)] <- "Fresh sediment"
df$habitat[!is.na(df$habitat) & df$habitat == "Marine" & !is.na(df$isolation_source) & grepl('sediment',df$isolation_source)] <- "Marine sediment"
df$habitat[!is.na(df$habitat) & df$habitat == "Fresh" & !is.na(df$isolation_source) & grepl('water',df$isolation_source)] <- "Fresh water"
df$habitat[!is.na(df$habitat) & df$habitat == "Marine" & !is.na(df$isolation_source) & grepl('water',df$isolation_source)] <- "Marine water"


####################################################
# Temperature adjusted growth rate using residuals #
####################################################
#NOTE: This is a form of temperature adjusted growth rates where residuals = adjusted growth rate

# Get data set for creating regression model on growth rate against growth temperature
sub <- df %>% filter(!is.na(growth_rate) & !is.na(growth_tmp))

#Ordinary least squares using lm() 
model <- lm(log10(growth_rate)~growth_tmp, data = sub)

#Use model to calculate residuals in original data frame
#measured growth rate minus estimated growth rate = residual
df <- df %>% mutate(
  temp_adjusted_maxgrowth = ifelse(!is.na(growth_rate), log10(growth_rate) - (model$coefficients['growth_tmp']*growth_tmp+model$coefficients['(Intercept)']), NA)) 

# Clean up
rm(model,sub)
   
###################################
# Save for plots on full data set #
###################################

full_df <- df %>% ungroup()

########################################################
# Restrict data frame to selected substrate use groups #
########################################################

# This creates a redundant data frame with each microbe 
# potentially represented by multiple rows

# Load list of substrate use groups we've decided to include in this study
pw <- read.csv("data/pathways_to_include.csv", as.is=TRUE)

# Extract all microbes in our data frame that fits into either of these pathways

all <- df[0,]

for(i in 1:nrow(pw)) {
  
  pathway <- pw$pathway[i]
  
  # Make sure to search for full terms 
  sub <- df[grepl(sprintf("\\b%s\\b",pathway),df$pathways),]
  
  sub$pathway <- pathway

  all <- all %>% bind_rows(sub)
  print(pathway);
}

# Add all Cyanobacteria as "oxygenic_photoautotrophy"
sub <- df[!is.na(df$phylum) & df$phylum == "Cyanobacteria",]
sub$pathway <- "oxygenic_photoautotrophy"

# Add to temporary data frame
all <- all %>% bind_rows(sub)

# Join naming info from table to include appropriate names and abbreviations
all <- all %>% inner_join(pw, by = "pathway")


###############################
# Deal with nitrogen fixation #
###############################

# Split out all cyanobacteria with N-fixation as cyanobacterial_N_fixation
sub <- all %>% filter(phylum == "Cyanobacteria" & pathway == "nitrogen_fixation")

sub$pathway <- "cyanobacterial_nitrogen_fixation"
sub$pathway_name <- pw$pathway_name[pw$pathway == "cyanobacterial_nitrogen_fixation"]
sub$group_name <- pw$group_name[pw$pathway == "cyanobacterial_nitrogen_fixation"]
sub$group_abbreviation <- pw$group_abbreviation[pw$pathway == "cyanobacterial_nitrogen_fixation"]
sub$full_group_name <- pw$full_group_name[pw$pathway == "cyanobacterial_nitrogen_fixation"]

# Add to temporary data frame
all <- all %>% bind_rows(sub)

# Remove all cyanobacteria from the group nitrogen_fixation
all <- all %>% filter(!(pathway == "nitrogen_fixation" & phylum == "Cyanobacteria"))

#Rename the remaining nitrogen fixers to "nitrogen_fixation_other"
all$pathway[!is.na(all$pathway) & all$pathway == "nitrogen_fixation"] <- "nitrogen_fixation_other"
all$pathway_name[!is.na(all$pathway) & all$pathway == "nitrogen_fixation_other"] <- pw$pathway_name[pw$pathway == "nitrogen_fixation_other"]
all$group_name[!is.na(all$pathway) & all$pathway == "nitrogen_fixation_other"] <- pw$group_name[pw$pathway == "nitrogen_fixation_other"]
all$group_abbreviation[!is.na(all$pathway) & all$pathway == "nitrogen_fixation_other"] <- pw$group_abbreviation[pw$pathway == "nitrogen_fixation_other"]
all$full_group_name[!is.na(all$pathway) & all$pathway == "nitrogen_fixation_other"] <- pw$full_group_name[pw$pathway == "nitrogen_fixation_other"]

#Rename 'all' back to 'df'
df <- all

#Ensure that each species is only included once in each pathway category 
df <- df %>% distinct(species, group_name, .keep_all = TRUE)
#(ex. a species with methylotrophy and methanol oxidation would end up in the same category 'Methylotrophy' twice)


############
# Finalise #
############

#Change factor levels of habitat
# df$habitat <- factor(df$habitat, levels = c("Fresh water","Fresh sediment","Marine water","Marine sediment","Soil","Thermal","Endotherm","Ectotherm","Intracellular","Other"))

#Convert genome size to Mbp in both data frame
df <- df %>% mutate(genome_size = genome_size/1000000) 
full_df <- full_df %>% mutate(genome_size = genome_size/1000000) 

#Remove any substrate use groups where nubmer of species n < 20
df <- df %>% group_by(group_name) %>% mutate(n = n()) %>% 
  filter(n >= 20) %>% 
  ungroup()


# Create plot order
plot_order <- c("Denitrific",
                "Fe red",
                "AsMnSeU red",
                "S-thio red",
                "Fumar red",
                "Acetogen",
                "H ox dark",
                "H2S-S-th ox",
                "Fe ox dark",
                "Nitrific", 
                "placeholder1",
                "Methylotrop",
                "Anox phot S",
                "Anox phot H",
                "Phot hetero",
                "O-gen phot",
                "Cyan N-fix",
                "Other N-fix",
                "CCLX deg", 
                "placeholder2", 
                "placeholder3",
                "CH4 gen",
                "Fermentat",
                "Aero hetero")


# Define colours for each substrate use group (this is fixed throughout)
pathway_colours <- c("#FF0000","#FF7F00","#FFD400","#6AFF00","#00EAFF","#0095FF","#0040FF","#AA00FF","#FF00AA",
                     "#EDB9B9","#E7E9B9","#B9EDE0","#DCB9ED","#8F2323","#8F6A23","#4F8F23","#23628F","#6B238F",
                     "#000000","#737373","#CCCCCC")

names(pathway_colours) <- unique(df$group_abbreviation)

#Add "Other" as gray for PCAs
pathway_colours_PCA <- c(pathway_colours,"#FFFFFF")
names(pathway_colours_PCA)[22] <- "Other"


#Check all groups
table(df$group_name)