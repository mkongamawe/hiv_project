library(tidyverse)
library(artyfarty)
library(ggtext)

time_to_cd4_data <- read_csv("data/time_to_cd4_recovery.csv")

ggplot(
    data = time_to_cd4_data,
    aes(month, count, colour = stage)
) +
    geom_point() +
    geom_line() +

    facet_wrap(
        ~study,
        nrow = 2,
        scales = "free_x",
        labeller = labeller(
            study = c(
                ding = "Ding et al. (2016)",
                le = "Le et al. (2013)"
            )
        )
    ) +

    labs(
        x = "Months since ART",
        y = "CD4+ cells/\u00b5L",
        title = "CD4+ recovery trajectory",
        colour = NULL
    ) +

    scale_colour_discrete(
        labels = c(
            early = "early ART",
            late = "late ART"
        )
    ) +

    theme_minimal(
        base_family = "Bookman Old Style"
    ) +

    theme(
        # Overall plot appearance
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),

        panel.grid.major = element_line(
            colour = "white",
            linewidth = 0.4
        ),

        panel.grid.minor = element_line(
            colour = "gray80",
            linewidth = 0.4,
            linetype = 2
        ),

        axis.line = element_line(
            colour = "black"
        ),

        # Strip text appearance
        strip.background = element_rect(
            fill = "oldlace"
        ),

        strip.text = element_markdown(
            face = "bold",
            family = "Bookman Old Style",
            colour = "black",
            size = 15
        ),

        # Legend at bottom, centred
        legend.position = "bottom",
        legend.justification = "center",

        # Plot title
        plot.title = element_markdown(
            family = "Bookman Old Style",
            face = "bold",
            hjust = 0.5,
            vjust = 0.5,
            size = 15
        )
    )

ggsave(
    "Time to CD4 recovery.png",
    height = 8,
    width = 10,
    path = "results/figures"
)
