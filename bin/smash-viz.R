#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(argparser, quietly = TRUE)))
suppressMessages(suppressWarnings(library(dplyr, quietly = TRUE)))
suppressMessages(suppressWarnings(library(ggplot2, quietly = TRUE)))
suppressMessages(suppressWarnings(library(purrr, quietly = TRUE)))
suppressMessages(suppressWarnings(library(stringr, quietly = TRUE)))

# Create a parser
p <- arg_parser("Produce a SMASH plot from results txt file")

# Add command line arguments
p <- add_argument(p, "--arms", help = "Definitions of chromosomal arms")
p <- add_argument(p, "--input", help = "Points to plot")
p <- add_argument(p, "--highlight", help = "Chromosomal arm(s) to highlight; if multiple, separate with comma")
p <- add_argument(p, "--color", help = "The color to use for highlighting", default = "blue")
p <- add_argument(p, "--title", help = "The title to use for the resulting plot")
p <- add_argument(p, "--output", help = "The file to write the plot to")

# Parse the command line arguments
argv <- parse_args(p)

if (!(argv$color %in% colors())) {
    stop(argv$color, " is not a valid color.")
}

arms <- read.csv(argv$arms, stringsAsFactors = FALSE)
df <- read.csv(argv$input, stringsAsFactors = FALSE, sep = "\t")

arms.clean <- arms %>%
    mutate(chrom = ifelse(chrom == "X", 23, chrom)) %>%
    mutate(chrom = ifelse(chrom == "Y", 24, chrom)) %>%
    mutate(chrom = as.numeric(chrom)) %>%
    arrange(chrom) %>%
    mutate(end = cumsum(as.numeric(length))) %>%
    mutate(chrom = ifelse(chrom == 23, "X", chrom)) %>%
    mutate(chrom = ifelse(chrom == 24, "Y", chrom)) %>%
    mutate(chrarm = paste0(chrom, arm)) %>%
    mutate(start = lag(end, default = 0)) %>%
    select(chrarm, start, end) %>%
    mutate(across(c(start, end), as.numeric))

df.clean <- df %>%
    filter(chr != "chrY") %>% # Highly repetitive so results are not high-quality
    mutate(across(c(abspos, ratio, seg_ratio), as.numeric)) %>%
    rowwise() %>%
    mutate(chrarm = arms.clean %>% rowwise() %>% filter(between(abspos, start, end)) %>% pluck("chrarm")) %>%
    ungroup() %>%
    mutate(
        ratio = ratio * 2,
        seg_ratio = seg_ratio * 2,
        abspos = abspos / 1e9
    ) %>%
    rowwise() %>%
    mutate(highlight = chrarm %in% unlist(str_split(argv$highlight, ","))) %>%
    ungroup()

chrom.breaks <- df.clean %>%
    group_by(chr) %>%
    summarise(minmax = mean(abspos)) %>%
    pluck("minmax") %>%
    unlist() %>%
    sort()

chrom.lines <- df.clean %>%
    group_by(chr) %>%
    summarise(minmax = c(ifelse(chr == "chrX", c(min(abspos), max(abspos)), min(abspos)))) %>%
    pluck("minmax") %>%
    unlist() %>%
    unique() %>%
    sort()

df.clean %>%
    mutate(seg_ratio = ifelse(abspos %in% chrom.lines, NA, seg_ratio)) %>% # Prevent drawing a trend line connecting chromosomes
    ggplot() +
    geom_point(
        data = df.clean %>% filter(!highlight),
        mapping = aes(x = abspos, y = ratio),
        color = ifelse(length(argv$highlight) > 0, "gray50", "black")
    ) +
    geom_point(
        data = df.clean %>% filter(highlight),
        mapping = aes(x = abspos, y = ratio),
        color = argv$color
    ) +
    geom_line(aes(x = abspos, y = seg_ratio), color = "red", size = 1) +
    geom_vline(
        data = data.frame(v = chrom.lines),
        mapping = aes(xintercept = v),
        linetype = "dashed"
    ) +
    labs(title = argv$title, x = "Genome Position (Gb)", y = "Copy Number") +
    scale_x_continuous(sec.axis = dup_axis(
        breaks = chrom.breaks,
        labels = c(1:22, "X"),
        name = NULL,
        guide = guide_axis(check.overlap = TRUE) # Drop labels on overlap
    )) +
    theme_classic() +
    theme(
        legend.position = "none",
        axis.title.x = element_text(hjust = 1),
        axis.title.y = element_text(hjust = 1, size = 18),
        axis.text = element_text(size = 14),
        axis.text.y = element_text(size = 18),
        title = element_text(size = 20),
        panel.border = element_rect(colour = "black", fill = NA)
    ) +
    ylim(0, 5)

ggsave(argv$output, width = 16, height = 9)
