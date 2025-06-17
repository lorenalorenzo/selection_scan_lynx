Comparative Genome Scans Reveal Species-Specific Selection in the Genus
Lynx
================
Lorena Lorenzo
17-Jun-2025

Here we load all the necessary R packages for plotting, selection scans,
and GO enrichment. We also set paths and some helper variables like
species codes, chromosome names, and plot colors.

``` r
#install.packages
library(tidyverse)
library(rehh)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(grid)
library(ggrepel)
library(writexl)
library(colorspace)
library(GOplot)
library(clusterProfiler)
library(gridExtra)
library(patchwork)
library(ggplot2)

#Set directories
path<- "/Users/lorenalorenzo/selection_scan_lynx"
files<- "/files/"
plots<- "/plots/"

#Name the variables
species<- c("lc", "ll", "lp", "lr")

names <- c("lc" = "Lynx canadensis",
            "lr" = "Lynx rufus",
            "ll" = "Lynx lynx",
            "lp" = "Lynx pardinus")

chr <- c( "A1", "A2", "A3", "B1", "B2", "B3", "B4", "C1", "C2", "D1", "D2", "D3", "D4", "E1", "E2", "E3", "F1", "F2")

colors <- c("lr"= "#FFB3B5",
            "lp"= "#BBCF85",
            "ll"="#61D8D6",
            "lc"="#D4BBFC")

#set X axis with chr_size dataframe (for plotting)
    axis_set <- read.table(paste0(path, files, "chr_size.txt"), sep="\t", col.names=c("chr", "size")) %>%
                 mutate(center = size/2)
    #Group by chr
    data_axis <- axis_set %>%
                  mutate(size_cum = lag(cumsum(as.numeric(size)), default = 0)) %>%
                  rowwise() %>%
                  mutate(center_cum = sum(center, size_cum))  
```

## 1. Genomic scans: saltiLASSI

After running the analyses, here we’re going to explore and plot the
results. I used a common treshold of 1% of highest statistic value to
define an outlier windows.

``` r
#Decide the threshold being used for this representation
percentage<- 0.99
  
for (sp in species)
{
# saltiLASSI results -----------------------------------------------------------

  #Read data
  lassi_scan<- read.table(paste0(path, files, sp, "_salti.lassip.hap.out"), sep="\t", header=T)
  
  #Group by chr
  data_cum <- lassi_scan %>% 
    group_by(chr) %>% 
    summarise(max_bp = max(pos)) %>%
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
    dplyr::select(chr, bp_add)
  
  data_lassi <- lassi_scan %>% 
    inner_join(data_cum, by = "chr") %>% 
    mutate(bp_cum = pos + bp_add) %>%
    mutate(dataset="saltiLASSI")

  #Change the "sp_L" column name to another in order to make it easier to manage
  colnames(data_lassi)[12]<- "statistic"

  ## define outliers ---------------------------------------------------------
  threshold <- sort(data_lassi$statistic)[round(length(data_lassi$statistic) * percentage)]
  lassi_outliers<- data_lassi %>% filter(statistic >= threshold)

  write.table(lassi_outliers, file=(paste0(path, files, sp, "_lassi_outliers")), sep="\t", row.names = FALSE, quote = FALSE)
   
  bed_out <- lassi_outliers  %>% dplyr::select(c(1,2,3,12))
  
  write.table(bed_out, file=(paste0(path, files, sp, "_lassi_outliers.bed")), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
 
  
## explore statistic distribution ---------------------------------------------------------
  dist_plot<- ggplot(data_lassi, aes(x= statistic)) +
         geom_density() + 
         labs (title=sp) +
        theme_minimal() +
        ggtitle (as.character(names[sp]))+
        geom_vline (xintercept= threshold, color="red")
  
  ggsave(filename = paste0(path, plots, sp, "_lassi_distribution.png"), plot= dist_plot, width= 12, height = 5 )

## plot the results ---------------------------------------------------------
    #customize the color palette
      # Get the darker version of the chosen species color
      darker_color <- darken(colors[sp], amount = 0.3)
      

  q <- ggplot() +
  geom_point(data_lassi, mapping=aes(x = bp_cum, y = statistic, 
                                     alpha = 0.5, color = as_factor(chr)), size= 2) +   
  geom_hline(yintercept= threshold, 
             linetype="dashed", show.legend = TRUE) +
  scale_x_continuous(label = data_axis$chr, breaks = data_axis$center_cum) +
  ylab ("Λ") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), legend.position = "none",
        axis.text = element_text(size = 18),  # Increase axis text size
        axis.title = element_text(size = 18),  # Increase axis title size
        #plot.title = element_text(size = 18),  # Increase plot title size
        panel.grid.minor = element_blank(), # Remove minor grid
        panel.grid.major = element_blank() # Remove major grid
        ) +  
        scale_color_manual(values=rep(c(as.character(colors[sp]), as.character(darker_color)), 9))
  
  assign(paste0(sp, "_plot"), q )
ggsave(filename = paste0(path, plots, sp, "_lassi_selection.png"), plot= q, width= 16, height = 3 )
  }  
```

### 1.1 Tryal with circular plot

In this section, I aim to experiment with circular plots to visualize
selection signals across chromosomes using tools like `circlize`

``` r
for (sp in species)
{

# saltiLASSI results -----------------------------------------------------------

# Read data
lassi_scan <- read.table(paste0(path, files, sp, "_salti.lassip.hap.out"), sep = "\t", header = TRUE)

# Group by chr
data_cum <- lassi_scan %>%
  group_by(chr) %>%
  summarise(max_bp = max(pos)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  dplyr::select(chr, bp_add)

data_lassi <- lassi_scan %>%
  inner_join(data_cum, by = "chr") %>%
  mutate(bp_cum = pos + bp_add) %>%
  mutate(dataset = "saltiLASSI")

# Rename the statistic column
colnames(data_lassi)[12] <- "statistic"

# Prepare data in genomic format for circlize
data_circ <- data_lassi %>%
  dplyr::select(chr, start, end, statistic)

# Set threshold
percentage <- 0.99
threshold <- sort(data_circ$statistic)[round(length(data_circ$statistic) * percentage)]

# Define chromosome color
chroms <- unique(data_circ$chr)
n_chr <- length(chroms)
species_color <- colors[sp]
chr_colors <- rep(c(species_color, darken(species_color, 0.3)), length.out = n_chr)
names(chr_colors) <- chroms

# Start image output
png(filename = paste0(path, plots, sp, "_lassi_circular_selection.png"),
    width = 1200, height = 1200, res = 150, bg = "transparent")

# Initialize circos plot
circos.clear()
circos.par(
  track.height = 0.45,  # taller track
  start.degree = 90,
  gap.degree = 2,
  canvas.xlim = c(-1.5, 1.5),  # more space for radial axis
  canvas.ylim = c(-1.5, 1.5),
  cell.padding = c(0.01, 0.01, 0.01, 0.01)
)

circos.initialize(factors = data_circ$chr, x = data_circ$start)

# Plot statistic values
circos.trackPlotRegion(
  factors = data_circ$chr,
  y = data_circ$statistic,
  ylim = c(0, max(data_circ$statistic)),
  track.height = 0.45,
  bg.border = NA,
  panel.fun = function(region, value, ...) {
    chr = CELL_META$sector.index
    x = data_circ$start[data_circ$chr == chr]
    y = data_circ$statistic[data_circ$chr == chr]

    # Points
    circos.points(x, y, col = chr_colors[chr], pch = 16, cex = 0.6)

    # Threshold line
    circos.lines(c(min(x), max(x)), c(threshold, threshold), col = "red", lty = 2)
  }
)

# Add chromosome names INSIDE the circle
circos.trackPlotRegion(
  ylim = c(0, 1),  # ✅ Add this line
  track.height = 0.05, bg.border = NA,
  panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ycenter - 0.1,
                labels = CELL_META$sector.index,
                facing = "bending.inside",
                niceFacing = TRUE,
                adj = c(0.5, 0.5),
                cex = 0.6)
  }
)

# Save the plot
dev.off()
}
```

## 2. Candidate regions and genes

Here I merge every overlapping outlier saltiLASSI window and calculate
the size of the merged regions. Then I extend the candidate regions up
and downstream 20kbp. Finally, I extract the genes overlapping with
these candidate regions. Candidate regions are defined without
extension, so extension is only considered for annotation.

``` bash
module load bedtools

species=(lc ll lp lr)

for sp in ${species[@]}
  do
    echo "$sp"
       #merge every overlapping outlier window in LASSI
       bedtools merge \
       -i ${sp}_lassi_outliers.bed \
       -c 1,4 \
       -o count,mean \
       > ${sp}_merged_tmp
       echo "$sp merged"
       
       #calculate the size
       awk '{print $1, $2, $3, $3-$2, $4, $5}' ${sp}_merged_tmp > ${sp}_merged2_tmp
       #print a header (column names)
       echo -e "chr start end size windows mean_stat" | cat - ${sp}_merged2_tmp | tr " " "\t" > ${sp}_candidate_regions
  done
  
  #Extend candidate regions up and downstream 20kbp
  for sp in ${species[@]}
  do    
      bedtools slop \
        -i  ${sp}_merged_tmp \
        -g chr_size.txt \
        -b 20000 \
        > ${sp}_extended_20kbp_tmp
        
      #print a header (column names)
       echo -e "chr start end windows" | cat - ${sp}_extended_20kbp_tmp | tr " " "\t" \
       > ${sp}_candidate_regions_extended_20kbp
       
  done  
  
  
  #Cross with annotation    
  for sp in ${species[@]}
  do    
     ###without extension 
        #1. Intersect with annotation file     
        bedtools intersect \
             -a ${sp}_candidate_regions \
             -b Felis_catus.Felis_catus_9.0.97.gff3 \
             -wa -wb \
             > ${sp}_candidate_regions_annotated 
            echo "$sp annotated"
        
        #2. Get only genes
          awk '$9=="gene"' ${sp}_candidate_regions_annotated  > ${sp}_genes_tmp
          echo "$sp genes"
        
        #3. Filter for interesting colums
           cut -f-1,2,3,4,5,6,10,11,15- ${sp}_genes_tmp | tr ' ' '\t' | cut -d';' -f1,2 | tr ';' '\t' | \
         awk '{if ($10 ~ /Name=[[:alnum:]]/) $10=$10; else $10="NA"; print $0}' | \
         awk '{if ($9 ~ /^ID=gene:/) {sub(/^ID=gene:/, "", $9)} if ($10 ~ /^Name=/) {sub(/^Name=/, "", $10)} print}' \
         > ${sp}_genes_filtered_tmp

        #4. Add a header
          echo -e "chr start end size windows mean_lassi gene_start gene_end ensembl_id gene_name" | \
          cat - ${sp}_genes_filtered_tmp > ${sp}_candidate_regions_genes
      
    ###with a 20kbp extension
        #1. Intersect with annotation file     
        bedtools intersect \
           -a ${sp}_candidate_regions_extended_20kbp \
           -b Felis_catus.Felis_catus_9.0.97.gff3 \
           -wa -wb \
           > ${sp}_candidate_regions_extended_20kbp_annotated 
          echo "$sp annotated"
          
        #2. Get only genes
          awk '$8=="gene"' ${sp}_candidate_regions_extended_20kbp_annotated   > ${sp}_genes_extend_tmp
          echo "$sp genes" 
      
        #3. Filter for interesting colums
          cut -f-1,2,3,4,5,9,10,14- ${sp}_genes_extend_tmp | tr ' ' '\t' | cut -d';' -f1,2 | tr ';' '\t' | \
            awk '{if ($9 ~ /Name=[[:alnum:]]/) $9=$9; else $9="NA"; print $0}' | \
            awk '{if ($8 ~ /^ID=gene:/) {sub(/^ID=gene:/, "", $8)} if ($9 ~ /^Name=/) {sub(/^Name=/, "", $9)} print}' \
            > ${sp}_genes_extend_filtered_tmp
      
        #4. Add a header
          echo -e "chr start end windows mean_lassi gene_start gene_end ensembl_id gene_name" | \
          cat - ${sp}_genes_extend_filtered_tmp > ${sp}_candidate_regions_20kbp_genes
   
   done 
  
  #remove temporal files
  rm *tmp
```

## 3. Enrichment analyses

Here I run enrichment analyses (ORA) for the candidate regions and genes
using the topGO package. I use the Felis catus Ensembl annotation from
ensembl.org to get all the genes and related GO terms.

``` r
##Dependencies
#BiocManager::install("biomaRt")
library(topGO)
library(biomaRt)

#####Get annotation from ensembl.org#####

#read from ensembl.org every felcat ensembl ID:
ensembl <- useMart("ensembl", dataset = "fcatus_gene_ensembl")  
###extract the GOterms for every ensembl_id
##ensembl_to_go <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "go_id"), ##mart = ensembl)

# Get attributes
attributes <- c("ensembl_gene_id", "external_gene_name", "go_id", "gene_biotype")

# Retrieve annotations
ensembl_to_go <- getBM(attributes = attributes, mart = ensembl) %>%
  filter(gene_biotype== "protein_coding")

#make a list with every GO terms per ENSEMBL_ ID. 
go_list <- split(ensembl_to_go$go_id, ensembl_to_go$ensembl_gene_id)

  #comment: Without filtering for coding genes, there is a total of 29550 ensembl_id corresponding to coding genes (19588) + non-coding genes (9468) + pseudogenes (494) according to https://www.ensembl.org/Felis_catus/Info/Annotation. We are getting 19564, 24 less than reported by the assembly (don't really know why)

for (sp in species)
{

#####read my gene data set#####  
    df_genes<- read.table(paste0(path, files, sp, "_candidate_regions_20kbp_genes"), 
                       sep=" ", header=T)
    genes <- df_genes %>%
              distinct(ensembl_id) %>%
              pull(ensembl_id) %>%
              as.character()  
    #this is a character object with the ensembl_ids ("ENSFCAG00000008235" "ENSFCAG00000008236" "ENSFCAG00000029529" ...)
    assign(paste0(sp, "_genes"), genes)

#cross the felcat annotation (ensembl_id with its GOterms) with my set of genes, to get a logical factor of true (if the gene is in the set) and false (if its not)
allgenes = factor(as.integer(names(go_list) %in% genes))
names(allgenes) <- names(go_list) #ensure names of allgenes are names of go_list (that is, the ensembl_id)

assign(paste0(sp, "_allgenes"), allgenes)


#####Creating the topGOdata object#####
godata <- new("topGOdata", ontology = "BP" , allGenes = allgenes, 
              #nodeSize= 10,
              annotationFun = annFUN.gene2GO, gene2GO = go_list)

assign(paste0(sp, "_BP_godata"), godata)
#####run the overrepresentation test#####
over_test <- runTest(godata, statistic = "fisher") #algorithm = "weight01" by default
assign(paste0(sp, "_test"), over_test)

#Get significant results without FDR correction
result_table <- GenTable(godata, Fisher=over_test, topNodes=over_test@geneData[2], numChar=1000) %>%  
  as_tibble() %>% 
  filter(Fisher<0.05 & Significant >1) 

  assign(paste0(sp, "_results"), results)

# Filter ensembl_to_go to include only ensembl_gene_id present in genes and select relevant columns
filtered_ensembl_to_go <- ensembl_to_go %>%
  filter(ensembl_gene_id %in% genes) %>%
  dplyr::select(ensembl_gene_id, external_gene_name, go_id)

# Group by GO.ID and summarize ensembl_gene_id and external_gene_name
grouped_ensembl_to_go <- filtered_ensembl_to_go %>%
  group_by(go_id) %>%
  summarize(
    ensembl_gene_id = paste(unique(ensembl_gene_id), collapse = ","),
    external_gene_name = paste(unique(external_gene_name), collapse = ",")
  )

# Merge result_table with grouped_ensembl_to_go to add ensembl_gene_id and external_gene_name
result_table <- result_table %>%
  left_join(grouped_ensembl_to_go, by = c("GO.ID" = "go_id"))

assign(paste0(sp, "_results"), result_table)
#######Save the enrichment results########
write.table(result_table, file=(paste0(path, files, sp, "_functional_enrichment.csv")), sep =  "\t", row.names = FALSE, quote = FALSE)
}
```

### 3.1 Go graph

``` r
library(ggplot2)
library(cowplot)  # for get_legend() and plot_grid
library(colorspace)

x_limits <- c(-1, 70)  # <- replace with actual range
x_breaks <- seq(0, 70, 10)  # <- define appropriate breaks



for (sp in species) {
   #Load and prepare data
##without filtering
data <- read.csv(paste0(path, files, sp, "_functional_enrichment.csv"), sep = "\t") %>%
  mutate(
    Fisher = as.numeric(Fisher),
    Significant = as.numeric(Significant),
    Fold_Enrichment = Significant / Expected,
    GO.ID = fct_reorder(GO.ID, Fold_Enrichment)  # <- this replaces Term
   )

##with filtering
#data <- get(paste0(sp, "_filtered"))

  # Define species-specific color gradient
  color_gradient <- c(
    darken(colors[sp], 0.5),
    darken(colors[sp], 0.2),
    colors[sp]
  )

  # Plot
  p <- ggplot(data, aes(Fold_Enrichment, GO.ID, size = Significant, fill = Fisher)) +
    geom_point(shape = 21, stroke = 0) +
    scale_size_continuous(range = c(2, 8), name = "Significant") +
    scale_fill_gradientn(colours = color_gradient, name = "Fisher p-value") +
    scale_x_continuous(
      #limits = x_limits, 
      breaks = x_breaks) +
    labs(x = "Fold Enrichment", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 8),
      axis.ticks.x = element_blank()
    )

  # Dynamically set height
  plot_height <- max(3, 0.15 * nrow(data))

  # Save
  ggsave(
    filename = paste0(path, plots. "GO_enrichment_", sp, "_GOID.png"),
    plot = p,
    width = 3,
    height = plot_height,
    dpi = 300
  )

}
####for extracting the legend (I already have it, no need to include in the for loop)
  # Create plot WITH legend (for extraction)
for (sp in species) {
  plot_with_legend <- base_plot +
    theme(legend.position = "right")

  # Extract legend
  legend <- cowplot::get_legend(
  plot_with_legend +
    guides(
      size = guide_legend(override.aes = list(stroke = 0.4, color = "black")),
      fill = guide_colorbar(order = 2)
    )
)


  # Save the legend
  ggsave(
    filename = paste0(path, plots, "GO_enrichment_legend_", sp, ".png"),
    plot = cowplot::plot_grid(legend),
    width = 3,
    height = 10,
    dpi = 300
  )
} 
```

### 3.2 Ancestor terms for GO IDs

``` r
library(GO.db)
library(AnnotationDbi)
library(dplyr)

# Function to retrieve broader ancestors
get_broader_terms <- function(go_id) {
  tryCatch({
    # Retrieve all ancestor terms
    ancestors <- AnnotationDbi::get(go_id, GOBPANCESTOR)
    if (is.null(ancestors)) return(NA)  # Handle missing ancestors
    ancestors
  }, error = function(e) NA)  # Handle errors
}

# List of species data frames and their names
species_data <- list(
  lc = lc_results,
  ll = ll_results,
  lp = lp_results,
  lr = lr_results
)

# Loop through each species
for (species_name in names(species_data)) {
  message("Processing species: ", species_name)
  
  # Extract GO terms for the species
  go_terms <- species_data[[species_name]]$GO.ID
  
  # Retrieve ancestors for each GO term
  ancestor_list <- lapply(go_terms, get_broader_terms)
  
  # Combine results into a data frame
  ancestor_df <- data.frame(
    GO.ID = rep(go_terms, sapply(ancestor_list, length)),
    Ancestor = unlist(ancestor_list),
    stringsAsFactors = FALSE
  )
  
  # Filter out input terms to keep only broader categories
  broader_categories <- ancestor_df %>%
    filter(!Ancestor %in% go_terms) %>%
    distinct()
  
  # Map ancestor descriptions using GO.db
  broader_categories <- broader_categories %>%
    mutate(Description = AnnotationDbi::select(GO.db, keys = Ancestor, columns = "TERM", keytype = "GOID")$TERM)
  
  # Save results to a CSV file
  output_file <- paste0(species_name, "_broader_go_categories.csv")
  write.csv(broader_categories, output_file, row.names = FALSE)
  
  message("Saved results to: ", output_file)
}
```

## 4. Assesing overlap between species candidate regions

``` bash

bedtools intersect -wa -wb \
       -a <(cut -f 1-3 lc_candidate_regions) \
       -b <(cut -f 1-3 ll_candidate_regions) <(cut -f 1-3 lp_candidate_regions) <(cut -f 1-3 lr_candidate_regions) \
       -names ll lp lr \
       -sorted \
       > lc_repeated_regions
       
bedtools intersect -wa -wb \
       -a <(cut -f 1-3 ll_candidate_regions) \
       -b <(cut -f 1-3 lc_candidate_regions) <(cut -f 1-3 lp_candidate_regions) <(cut -f 1-3 lr_candidate_regions) \
       -names lc lp lr \
       -sorted \
       > ll_repeated_regions
       
bedtools intersect -wa -wb \
       -a <(cut -f 1-3 lp_candidate_regions) \
       -b <(cut -f 1-3 lc_candidate_regions) <(cut -f 1-3 ll_candidate_regions) <(cut -f 1-3 lr_candidate_regions) \
       -names lc ll lr \
       -sorted \
       > lp_repeated_regions
       
bedtools intersect -wa -wb \
       -a <(cut -f 1-3 lr_candidate_regions) \
       -b <(cut -f 1-3 lc_candidate_regions) <(cut -f 1-3 ll_candidate_regions) <(cut -f 1-3 lp_candidate_regions) \
       -names lc ll lp \
       -sorted \
       > lr_repeated_regions
       
######Cross bed results with annotation file (gff3)    
    #match selected regions with genome annotation
species=(lc ll lp lr)

for sp in ${species[@]}
  do
    echo "$sp"
    bedtools intersect \
     -a ${sp}_repeated_regions \
     -b Felis_catus.Felis_catus_9.0.97.gff3 \
     -wa -wb \
      > ${sp}_annotated_tmp 
     echo "$sp annotated"
     #filter only genes
     awk '$10=="gene"' ${sp}_annotated_tmp  > ${sp}_gene_tmp
     echo "$sp gene filtered"   
     #get only interesting columns
    cut -f 1-8,11-12,16 ${sp}_gene_tmp > ${sp}_filtered_tmp
    #ensembl id and gene name as column data
    cut -d';' -f1,2 ${sp}_filtered_tmp | tr ';' '\t'  | awk '{if ($12 ~ /Name=[[:alnum:]]/) $12=$12; else $12="NA"; print $0}' > ${sp}_rep_genes_names_tmp 
    
    #cut extra-info from ens code and gene name columns
    paste <(cut -d" " -f-10 ${sp}_rep_genes_names_tmp) <(cut -d" " -f11 ${sp}_rep_genes_names_tmp | cut -d":" -f2) <(cut -d" " -f12 ${sp}_rep_genes_names_tmp | cut -d"=" -f2) |  tr "\t" " "  >> ${sp}_rep_columns_tmp
    
    #print a header (column names)
    echo -e "chr start end sp chr_2 start_2 end_2 gene_chr gene_start gene_end ensembl_id gene_name" | cat - ${sp}_rep_columns_tmp > ${sp}_repetitive_genomic_regions_annotated
      
        #remove temporal files
       rm *tmp
    done
```

### 4.1 Venn Diagram

In order to represent the repetitive selection between species in
relation to gene name, I will get the list of genes for each species and
compare those lists:

``` r
library(ggvenn)


# Initialize an empty list to store gene sets
gene_lists <- list()
ensembl_lists<- list()

# Read the gene data sets for each species and store in the list
for (sp in species) {
  df_genes <- read.table(paste0(path, files, sp, "_candidate_regions_20kbp_genes"), 
                         sep=" ", header=TRUE)
  ens_ids<- na.omit(df_genes$ensembl_id)
  ensembl_lists[[sp]]<- ens_ids
  genes <- na.omit(df_genes$gene_name)
  gene_lists[[sp]] <- genes
}

# Assign meaningful names to the gene sets
names(ensembl_lists) <- c("Lynx canadensis", "Lynx lynx", "Lynx pardinus", "Lynx rufus")

# Generate the Venn diagram using ggvenn
p <- ggvenn(
        ensembl_lists,
        fill_color = c(as.character(colors["lc"]), as.character(colors["ll"]), as.character(colors["lp"]), as.character(colors["lr"])),
        stroke_size = 0,
        show_percentage = FALSE,
        set_name_size = 6,
        text_size = 5,
)

ggsave(filename = paste0(path, plots, "genes_venn_diagram.pdf"), plot= p, width= 15, height = 10)
```

``` r
#check intersect btw GO terms
intersect(lc_results$GO.ID, ll_results$GO.ID)
```

We found intersection btw lynx canadensis and lynx lynx in the GO term:
methylation and btw l. canadensis and lynx pardinus in modulation by
host of viral process.

## 5. Customizing a supplementary table with the candidate regions info

``` r
for (sp in species) 
  {
  
# Read the files into data frames
df_regions <- read.table(paste0(path, files, sp, "_candidate_regions"), 
                         sep="\t", header=TRUE)
df_genes <- read.table(paste0(path, files, sp, "_candidate_regions_20kbp_genes"), 
                         sep=" ", header=TRUE)
df_genes <- df_genes %>%
  mutate(start = start + 20000,
         end = end - 20000) 

# Merge the data frames on the common columns
df_merged <- left_join(df_regions, df_genes, by=c('chr', 'start', 'end', 'windows', 'mean_lassi', 'max_lassi'))

# Concatenate gene_name and ensembl_id for the same genomic window
  table <- df_merged %>%
  group_by(chr, start, end, windows, size, mean_lassi, max_lassi) %>%
  summarise(gene_name = paste(na.omit(gene_name), collapse = ","),
            ensembl_id = paste(na.omit(ensembl_id), collapse = ","))   

#Save final candidate genomic table
#write.table(table, file=(paste0(path, files, sp, "_candidate_region.csv")), 
#              sep =  ",", row.names = FALSE, quote = FALSE)
#
#write_xlsx(table, path = paste0(path, files, sp, "_candidate_region.xlsx"))
#  
#
#assign(paste0(sp, "_table"), table)

# Create the histogram with custom x-axis breaks
w<- ggplot(table, aes(x = size / 1000)) +
  geom_histogram(binwidth = 50, color = "black") +
  scale_x_continuous(
    breaks = seq(floor(min(table$size) / 50000) * 50, 
                 ceiling(max(table$size) / 50000) * 50, 
                 by = 50),
    labels = scales::comma) +  # Optional: to format axis labels with commas
  labs(title = paste("Distribution of", names[sp], "candidate regions size "), 
       x = "Size (kb)", 
       y = "Frequency") +
  annotate("text", 
           x = max(table$size / 1000), 
           y = 12, 
           label = paste("Median:", round(median(table$size / 1000), 2), "kb"),
           hjust = 1,  # Align the text to the left of the x position
           color = "red")
assign(paste0(sp, "_size_hist"), w )
}
```
