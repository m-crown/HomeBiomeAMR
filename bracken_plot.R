library(stats) # v4.2.3
library(tidyverse) # v2.0.0
library(svglite) # v2.1.1
library(optparse)

plot_bracken <- function(file, top = NA, pal, otherColor){
  
  # file is the bracken data
  # top is an integer, trims the data to the number of top taxa (by median) in a given level for better plotting
  ## top = NA is the default and will plot all taxa in a given level
  # pal is a palette of hex colors
  # otherColor is the color used to represent taxa below median abundance threshold
  
  levdat <-
    dplyr::select(file,
                  -contains(c("_num","taxonomy_"))) |>
    tidyr::pivot_longer(cols = !name,
                        names_to = "sample",
                        values_to = "fraction") |>
    dplyr::mutate(sample = str_remove_all(sample, "_frac"))
  
  level <-
    dplyr::select(file,
                  contains("taxonomy_lvl")) |>
    unique()
  
  newpal <-
    paste0("#",
           stringr::str_split(pal, ",") |> unlist()
    )
  
  if (level == "K") {
    label <- "Kingdom"
  } else if (level == "P") {
    label <- "Phylum"
  } else if (level == "C") {
    label <- "Class"
  } else if (level == "O") {
    label <- "Order"
  } else if (level == "F") {
    label <- "Family"
  } else if (level == "G") {
    label <- "Genus"
  } else if (level == "S") {
    label <- "Species"
  } else if (level == "S1") {
    label <- "Strain"
  } else {
    message("Bracken taxonomy_lvl needs to be K, P, C, O, F, G, S, or S1")
  }
  
  if (!is.na(top) && top < length(unique(levdat$name))) {
    
    top_select <-
      stats::aggregate(x = levdat$fraction,
                       by = list(levdat$name),
                       FUN = mean) |>
      purrr::set_names(c("group","value")) |>
      dplyr::arrange(dplyr::desc(value)) |>
      dplyr::slice(1:top) |>
      dplyr::pull(var = group)
    
    levdat <- dplyr::filter(levdat, levdat$name %in% top_select)
    
    fill_unk <-
      stats::aggregate(x = levdat$fraction,
                       by = list(levdat$sample),
                       FUN = sum) |> 
      purrr::set_names(c("sample","sum")) |> 
      dplyr::mutate_if(is.numeric, round, 3)
    
    fill_unk$residual <- 1 - fill_unk$sum
    
    for (i in 1:nrow(fill_unk)) {
      
      levdat <- rbind(levdat, c("other", fill_unk[i,1], fill_unk[i,3]))
      
    }
    
    levdat <- transform(levdat, fraction = as.numeric(fraction))
    levdat$name <- factor(levdat$name, levels = unique(levdat$name))
    
    fill_vals <- c(rep_len(newpal, nrow(unique(levdat[1]))-1), otherColor)
  } else {
    
    fill_vals <- c(rep_len(newpal, nrow(unique(levdat[1]))))
  }
  
  print(nrow(levdat))

  ggplot2::ggplot(data = levdat) +
    ggplot2::geom_bar(mapping = ggplot2::aes(x = sample,
                                             y = fraction,
                                             fill = name),
                      position = "fill",
                      stat = "identity",
                      width = 0.75) + 
    ggplot2::scale_fill_manual(paste(label), values = fill_vals) + 
    ggplot2::theme_classic() +
    ggplot2::scale_x_discrete("") +
    ggplot2::scale_y_continuous("Relative Abundance") +
    ggplot2::theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 1, size = 12),
                   axis.text.y = element_text(size = 10),
                   axis.title.x = element_text(),
                   axis.title.y = element_text(size = 13),
                   legend.position = "right")
  
}

option_list = list(
  make_option(c("-d", "--dir"), type="character",
              help="directory of merged files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

notLogical <- cols(taxonomy_lvl = col_character())

family_file = readr::read_tsv(file = paste0(opt$dir, "merged_bracken_F.txt"),
                       col_names = TRUE,
                       col_types = notLogical) %>%
  select(-c(REANEG_frac, REANEG_num))
phylum_file = readr::read_tsv(file = paste0(opt$dir, "merged_bracken_P.txt"),
                              col_names = TRUE,
                              col_types = notLogical) %>%
  select(-c(REANEG_frac, REANEG_num))
class_file = readr::read_tsv(file = paste0(opt$dir, "merged_bracken_C.txt"),
                              col_names = TRUE,
                              col_types = notLogical) %>%
  select(-c(REANEG_frac, REANEG_num))
order_file = readr::read_tsv(file = paste0(opt$dir, "merged_bracken_O.txt"),
                              col_names = TRUE,
                              col_types = notLogical) %>%
  select(-c(REANEG_frac, REANEG_num))
genus_file = readr::read_tsv(file = paste0(opt$dir, "merged_bracken_G.txt"),
                              col_names = TRUE,
                              col_types = notLogical) %>%
  select(-c(REANEG_frac, REANEG_num))
species_file = readr::read_tsv(file = paste0(opt$dir, "merged_bracken_S.txt"),
                              col_names = TRUE,
                              col_types = notLogical) %>%
  select(-c(REANEG_frac, REANEG_num))

family_plot <- plot_bracken(file = family_file,
                   top = 10,
                   pal = "5c2751,ef798a,f7a9a8,00798c,6457a6,9dacff,76e5fc,a30000,ff7700,f5b841",
                   otherColor = "gray")



genus_plot <- plot_bracken(file = genus_file,
                            top = 20,
                            pal = "5c2751,ef798a,f7a9a8,00798c,6457a6,9dacff,76e5fc,a30000,ff7700,f5b841",
                            otherColor = "gray")

class_plot <- plot_bracken(file = class_file,
                            top = 10,
                            pal = "5c2751,ef798a,f7a9a8,00798c,6457a6,9dacff,76e5fc,a30000,ff7700,f5b841",
                            otherColor = "gray")



order_plot <- plot_bracken(file = order_file,
                            top = 10,
                            pal = "5c2751,ef798a,f7a9a8,00798c,6457a6,9dacff,76e5fc,a30000,ff7700,f5b841",
                            otherColor = "gray")



phylum_plot <- plot_bracken(file = phylum_file,
                            top = 10,
                            pal = "5c2751,ef798a,f7a9a8,00798c,6457a6,9dacff,76e5fc,a30000,ff7700,f5b841",
                            otherColor = "gray")



species_plot <- plot_bracken(file = species_file,
                            top = 20,
                            pal = "5c2751,ef798a,f7a9a8,00798c,6457a6,9dacff,76e5fc,a30000,ff7700,f5b841",
                            otherColor = "gray")

ggsave(filename = paste0(opt$dir,"bracken_plot_F.png"),
       plot = family_plot, units = "in",
       width = 10,
       height = 6)

ggsave(filename = paste0(opt$dir,"bracken_plot_G.png"),
                         plot = genus_plot, units = "in",
                         width = 10,
                         height = 6)
       
ggsave(filename = paste0(opt$dir,"bracken_plot_C.png"),
       plot = class_plot, units = "in",
       width = 10,
       height = 6)

ggsave(filename = paste0(opt$dir,"bracken_plot_O.png"),
       plot = order_plot, units = "in",
       width = 10,
       height = 6)

ggsave(filename = paste0(opt$dir,"bracken_plot_P.png"),
       plot = phylum_plot, units = "in",
       width = 10,
       height = 6)

ggsave(filename = paste0(opt$dir,"bracken_plot_S.png"),
       plot = species_plot, units = "in",
       width = 10,
       height = 6)


