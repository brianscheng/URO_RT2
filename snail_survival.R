# Snail survival: NC vs NH temperature regimes x 3 source populations
# Time origin = 1 Aug 2022 (first deployment); replacements enter late (Date_IN);
# survivors censored at sacrifice, 3 Dec 2024.
library(readxl); library(survival); library(coxphf); library(coxme)

d <- read_excel("RTdata_BSC_survival_analysis.xlsx")
t0  <- as.Date("2022-08-01"); end <- as.Date("2024-12-03")
d <- within(d, {
  entry <- as.numeric(as.Date(Date_IN) - t0)
  event <- as.integer(!is.na(Death_date))
  exit  <- ifelse(event == 1, as.numeric(as.Date(Death_date) - t0), as.numeric(end - t0))
  Treat <- factor(Treat, levels = c("NC", "NH"))
  Pop   <- factor(Pop,   levels = c("GA", "NC", "NH"))
})

# 1. Kaplan-Meier + log-rank
km <- survfit(Surv(entry, exit, event) ~ Treat + Pop, data = d)
plot(km, col = rep(c("#2a78d6", "#eb6834"), each = 3), lty = 1:3, ylim = c(0.6, 1))
survdiff(Surv(exit, event) ~ Treat + strata(Pop), data = d)  # note: survdiff ignores late entry;
                                                              # results nearly identical here (entry <= 22 d)
for (p in levels(d$Pop)) print(survdiff(Surv(exit, event) ~ Treat, data = subset(d, Pop == p)))

# 2. Cox, proportional hazards (shown for comparison) + PH test
m1 <- coxph(Surv(entry, exit, event) ~ Treat + Pop, data = d, cluster = Egg_tag)
summary(m1)
cox.zph(coxph(Surv(entry, exit, event) ~ Treat + Pop, data = d))  # expect strong violation for Treat

# 3. Piecewise Cox: separate regime effect in year 1 vs years 2+  (main model)
s <- survSplit(Surv(entry, exit, event) ~ ., data = d, cut = 365, episode = "period")
s$NH_yr1  <- (s$Treat == "NH") * (s$period == 1)
s$NH_yr2p <- (s$Treat == "NH") * (s$period == 2)
m2 <- coxph(Surv(entry, exit, event) ~ NH_yr1 + NH_yr2p + Pop, data = s)
summary(m2)
coxphf(Surv(entry, exit, event) ~ NH_yr1 + NH_yr2p + Pop, data = s)  # Firth: sparse NC-regime yr-1 deaths
coxme(Surv(entry, exit, event) ~ NH_yr1 + NH_yr2p + Pop + (1 | Egg_tag), data = s)  # capsule frailty

# 4. Population x regime (low power: 24 deaths)
m3 <- coxph(Surv(entry, exit, event) ~ Treat * Pop, data = d)
anova(coxph(Surv(entry, exit, event) ~ Treat + Pop, data = d), m3)

# 5. Sensitivity: original animals only (exclude replacements)
summary(coxph(Surv(entry, exit, event) ~ NH_yr1 + NH_yr2p + Pop, data = subset(s, Is_Replacement == "N")))

# 6. 
# Kaplan-Meier survival curves: NC vs NH temperature regime, by source population
# Recreates snail_KM_by_population.png (NH regime = blue, NC regime = red)
#
# Time origin = 1 Aug 2022 (first deployment). Replacement snails enter late
# (delayed entry at Date_IN). Survivors are censored at sacrifice, 3 Dec 2024.

library(readxl)
library(survival)
library(dplyr)
library(ggplot2)

# ---- settings ---------------------------------------------------------------
data_file <- "RTdata_BSC_survival_analysis.xlsx"
t0        <- as.Date("2022-08-01")          # time origin
end_date  <- as.Date("2024-12-03")          # sacrifice / censoring date
end_day   <- as.numeric(end_date - t0)

regime_cols <- c(NH = "#2a78d6",            # blue
                 NC = "#e34948")            # red
ink   <- "#3a3a38"
muted <- "#8a8a86"

panel_levels <- c("All populations", "GA population",
                  "NC population",   "NH population")

# ---- data -------------------------------------------------------------------
d <- read_excel(data_file) %>%
  mutate(
    entry = as.numeric(as.Date(Date_IN) - t0),
    event = as.integer(!is.na(Death_date)),
    exit  = ifelse(event == 1, as.numeric(as.Date(Death_date) - t0), end_day),
    Treat = factor(Treat, levels = c("NC", "NH"))
  )

# Stack the data so each snail appears in "All populations" and in its own panel
panels <- bind_rows(
  mutate(d, panel = "All populations"),
  mutate(d, panel = paste(Pop, "population"))
) %>%
  mutate(panel = factor(panel, levels = panel_levels))

# ---- Kaplan-Meier curves (delayed entry via counting-process Surv) ----------
km_curve <- function(dat) {
  fit <- survfit(Surv(entry, exit, event) ~ 1, data = dat)
  tibble(day  = c(0, fit$time, end_day),
         surv = c(1, fit$surv, tail(fit$surv, 1)))
}

curves <- panels %>%
  group_by(panel, Treat) %>%
  group_modify(~ km_curve(.x)) %>%
  ungroup() %>%
  mutate(date = t0 + day)

# ---- end-of-curve labels ("NH regime / 10/32 died") -------------------------
labels <- panels %>%
  group_by(panel, Treat) %>%
  summarise(n = n(), deaths = sum(event), .groups = "drop") %>%
  left_join(curves %>%
              group_by(panel, Treat) %>%
              summarise(s_end = last(surv), .groups = "drop"),
            by = c("panel", "Treat")) %>%
  group_by(panel) %>%
  # push the two labels apart when the curves end close together
  mutate(y = if (abs(diff(s_end)) < 0.05)
    mean(s_end) + ifelse(rank(-s_end, ties.method = "first") == 1,
                         0.027, -0.027)
    else s_end) %>%
  ungroup() %>%
  mutate(x = end_date + 15,
         label = paste0(Treat, " regime\n", deaths, "/", n, " died"))

year1_label <- tibble(panel = factor("All populations", levels = panel_levels),
                      x = t0 + 15, y = 0.61, label = "Year 1")

# ---- plot -------------------------------------------------------------------
p <- ggplot(curves, aes(date, surv, colour = Treat)) +
  annotate("rect", xmin = t0, xmax = t0 + 365, ymin = -Inf, ymax = Inf,
           fill = "#f0efec") +
  geom_step(linewidth = 0.9, direction = "hv") +
  geom_text(data = labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, size = 2.8,
            lineheight = 0.9, colour = ink) +
  geom_text(data = year1_label, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, size = 2.8, colour = muted) +
  facet_wrap(~ panel, nrow = 1) +
  scale_colour_manual(values = regime_cols,
                      breaks = c("NC", "NH"),
                      labels = c("NC temperature regime", "NH temperature regime"),
                      name = NULL) +
  scale_x_date(breaks = as.Date(c("2023-01-01", "2023-07-01",
                                  "2024-01-01", "2024-07-01")),
               date_labels = "%b\n%Y",
               limits = c(t0, end_date + 230),
               expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0.60, 1.00, 0.05)) +
  coord_cartesian(ylim = c(0.60, 1.01)) +
  labs(x = NULL, y = "Proportion surviving (Kaplan–Meier)") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position      = "top",
    legend.justification = "right",
    panel.grid.major.x   = element_blank(),
    panel.grid.minor     = element_blank(),
    panel.grid.major.y   = element_line(colour = "#e6e5e1", linewidth = 0.4),
    axis.line            = element_line(colour = muted, linewidth = 0.4),
    axis.ticks           = element_line(colour = muted, linewidth = 0.4),
    axis.text            = element_text(colour = muted),
    axis.title.y         = element_text(colour = ink),
    strip.text           = element_text(hjust = 0, size = 11, colour = ink),
    panel.spacing        = unit(1.2, "lines")
  )

print(p)
ggsave("snail_KM_by_population.png", p, width = 14, height = 4.4,
       dpi = 300, bg = "white")
