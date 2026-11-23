library(ggplot2)
library(dplyr)
library(ggtext)

time_to_cd4_data <- read.csv("data/time_to_cd4_recovery.csv")

# ---- Palette & labels -------------------------------------------------------
stage_colours <- c(early = "#0F7C80", late = "#E0703A")
stage_labels <- c(early = "Early ART", late = "Late ART")
study_labels <- c(
    ding = "Ding et al. (2016)",
    le = "Le et al. (2013)"
)
base_font <- "sans"

# Last point of each trajectory: used for the direct "Early/Late ART" labels
# and the "reached 900 at X months" callouts.
end_points <- time_to_cd4_data |>
    group_by(study, stage) |>
    slice_max(month, n = 1) |>
    ungroup() |>
    mutate(
        label = paste0(
            "<b style='color:",
            stage_colours[stage],
            "'>",
            stage_labels[stage],
            "</b><br>",
            "<span style='color:grey40'>",
            month,
            " months</span>"
        )
    )

# ---- Plot -------------------------------------------------------------------
ggplot(
    data = time_to_cd4_data,
    aes(month, count, colour = stage)
) +
    # Recovery threshold reference
    geom_hline(
        yintercept = 900,
        linetype = "dashed",
        colour = "grey55",
        linewidth = 0.4
    ) +

    geom_line(linewidth = 1.3, lineend = "round") +
    geom_point(
        shape = 21,
        size = 3.2,
        stroke = 1.4,
        fill = "white"
    ) +

    # Direct labels at the end of each line (replaces the legend)
    geom_richtext(
        data = end_points,
        aes(label = label),
        hjust = 0,
        nudge_x = 0.25,
        size = 3.6,
        family = base_font,
        fill = "white",
        label.colour = NA,
        label.padding = unit(c(0.15, 0.2, 0.15, 0.2), "lines"),
        lineheight = 1.1
    ) +

    facet_wrap(
        ~study,
        ncol = 1,
        scales = "free_x",
        labeller = labeller(study = study_labels)
    ) +

    scale_colour_manual(values = stage_colours, guide = "none") +
    scale_x_continuous(
        breaks = scales::breaks_pretty(n = 6),
        expand = expansion(mult = c(0.01, 0.16))
    ) +
    scale_y_continuous(
        limits = c(0, 1000),
        breaks = seq(0, 1000, 200),
        expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(clip = "off") +

    labs(
        x = "Months since ART initiation",
        y = "CD4<sup>+</sup> cells/µL",
        title = "CD4<sup>+</sup> recovery trajectory",
        # subtitle = paste0(
        #     "Patients starting ART ",
        #     "<b style='color:", stage_colours["early"], "'>early</b>",
        #     " reach 900 cells/µL sooner than those starting ",
        #     "<b style='color:", stage_colours["late"], "'>late</b>"
        # ),
        caption = "Dashed line: 900 cells/µL recovery threshold"
    ) +

    theme_minimal(base_family = base_font, base_size = 13) +

    theme(
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(20, 30, 15, 20),

        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
        panel.spacing.y = unit(1.6, "lines"),

        axis.line.x = element_line(colour = "grey30", linewidth = 0.5),
        axis.ticks.x = element_line(colour = "grey30", linewidth = 0.5),
        axis.ticks.length.x = unit(4, "pt"),
        axis.text = element_text(colour = "grey35"),
        axis.title.x = element_text(colour = "grey25", margin = margin(t = 10)),
        axis.title.y = element_markdown(
            colour = "grey25",
            margin = margin(r = 10)
        ),

        strip.text = element_text(
            face = "bold",
            colour = "grey15",
            size = 13,
            hjust = 0,
            margin = margin(b = 6)
        ),

        plot.title.position = "panel",
        plot.title = element_markdown(
            face = "bold",
            size = 20,
            colour = "grey10",
            hjust = 0.5,
            margin = margin(b = 18)
        ),
        plot.subtitle = element_markdown(
            size = 13,
            colour = "grey35",
            margin = margin(b = 18)
        ),
        plot.caption = element_text(
            size = 10,
            colour = "grey50",
            margin = margin(t = 12)
        )
    )

ggsave(
    "Time to CD4 recovery.png",
    height = 8,
    width = 10,
    dpi = 300,
    bg = "white",
    path = "results/figures"
)
