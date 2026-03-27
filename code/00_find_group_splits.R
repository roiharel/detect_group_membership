# Full fixed inter-group + within-group split analysis script
# - Robust column normalization (metadata using individual_local_identifier / group_id supported)
# - Detects within-group splits, numbers internal clusters largest->smallest
# - Produces per-group-per-date flags (split yes/no, missing data)
# - Robust plotting function that auto-detects/normalizes group column names and works when flags are missing
#
# Save and source this file, then call the functions in your session.
#
# Example:
source("code/intergroup_analysis_fixed.R")

cleaned_data <- readRDS("//10.126.19.90/EAS_shared/baboon/working/data/processed/2025/gps/v1_cleaned/gps_v1.RDS")

filtered_data <- cleaned_data %>%
  st_drop_geometry() %>%
  mutate(timestamp_dt = as.POSIXct(timestamp)) %>%
  #filter(year(timestamp_dt) == 2025) %>%
  #filter(month(timestamp_dt) == 8) %>%
  #filter(hour(timestamp_dt) > 15) %>%
  mutate(interval = floor_date(timestamp_dt, "hour")) %>%
  group_by(animal_id, interval) %>%
  slice_head(n = 1) %>%
  ungroup()


cluster_results <- filtered_data %>%
  group_by(interval) %>%
  group_modify(~ {
    coords <- .x %>% select(location.long, location.lat)
    db_result <- dbscan(coords, eps = 0.004, minPts = 1)
    .x %>% mutate(cluster_id = db_result$cluster)
  }) %>%
  ungroup()

cluster_centroids <- cluster_results %>%
  filter(cluster_id > 0) %>%  # Exclude noise points (cluster_id=0)
  group_by(interval, cluster_id) %>%
  summarize(
    centroid_lat = median(location.lat, na.rm = TRUE),
    centroid_lon = median(location.long, na.rm = TRUE),
    .groups = "drop"
  )

cluster_membership_data <- cluster_results %>%
  filter(cluster_id > 0) %>%
  left_join(cluster_centroids, by = c("interval", "cluster_id")) %>%
  rowwise() %>%
  mutate(distance_to_centroid = distHaversine(
    c(location.long, location.lat),
    c(centroid_lon, centroid_lat)
  )) %>%
  ungroup() %>%
  group_by(interval, cluster_id) %>%
  mutate(
    max_dist_in_cluster = max(distance_to_centroid, na.rm = TRUE),
    cluster_belonging = ifelse(
      max_dist_in_cluster == 0, 1,
      1 - (distance_to_centroid / max_dist_in_cluster)
    )
  ) %>%
  ungroup()


evening_locations <- cluster_membership_data %>%
  filter(hour(interval) == 16) %>%
  select(date = interval, animal_id = animal_id,
         group_id, cluster_id, location.lat, location.long) %>%
  mutate(date = as.Date(date))


prepped <- normalize_inputs(evening_locations, metadata,
             evening_animal_id_col = "animal_id",
             date_col = "date",
             cluster_col = "cluster_id",
             metadata_animal_id_col = "individual_local_identifier",
             metadata_group_col = "group_id")
splits_named <- identify_within_group_splits_numbered(prepped$evening_locations, prepped$animal_metadata)
flags <- group_day_split_flags(
  prepped$evening_locations,
  prepped$animal_metadata,
  metadata_group_col = "origin_group"   # already renamed by normalize_inputs
)
p <- plot_group_clusters_timeline(splits_named, flags, group_col = group_cols)
print(p)

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(plotly)
library(purrr)
library(forcats)
library(tibble)

# -------------------------
# UTIL: Normalize input column names
# -------------------------
normalize_inputs <- function(evening_locations, animal_metadata,
                             evening_animal_id_col = "animal_id",
                             date_col = "date",
                             cluster_col = "cluster_id",
                             metadata_animal_id_col = "animal_id",
                             metadata_group_col = "origin_group") {
  el <- as_tibble(evening_locations)
  meta <- as_tibble(animal_metadata)
  
  # rename evening_locations columns to standard names
  if (!"animal_id" %in% names(el)) {
    if (!(evening_animal_id_col %in% names(el))) stop("evening_locations missing animal id column: ", evening_animal_id_col)
    el <- el %>% rename(animal_id = !!sym(evening_animal_id_col))
  }
  if (!"date" %in% names(el)) {
    if (!(date_col %in% names(el))) stop("evening_locations missing date column: ", date_col)
    el <- el %>% rename(date = !!sym(date_col))
  }
  if (!"cluster_id" %in% names(el)) {
    if (!(cluster_col %in% names(el))) stop("evening_locations missing cluster id column: ", cluster_col)
    el <- el %>% rename(cluster_id = !!sym(cluster_col))
  }
  
  # rename metadata
  if (!"animal_id" %in% names(meta)) {
    if (!(metadata_animal_id_col %in% names(meta))) stop("animal_metadata missing animal id column: ", metadata_animal_id_col)
    meta <- meta %>% rename(animal_id = !!sym(metadata_animal_id_col))
  }
  if (!"origin_group" %in% names(meta)) {
    if (!(metadata_group_col %in% names(meta))) stop("animal_metadata missing group column: ", metadata_group_col)
    meta <- meta %>% rename(origin_group = !!sym(metadata_group_col))
  }
  
  # ensure date is Date
  if (!inherits(el$date, "Date")) {
    el <- el %>% mutate(date = as.Date(date))
    if (any(is.na(el$date))) stop("Could not convert evening_locations$date to Date. Provide Date or ISO strings.")
  }
  
  return(list(evening_locations = el, animal_metadata = meta))
}

# -------------------------
# WITHIN-GROUP SPLIT DETECTION (numbered)
# -------------------------
identify_within_group_splits_numbered <- function(evening_locations, animal_metadata,
                                                  continuity_gap_days = 1) {
  # Assumes evening_locations and animal_metadata already normalized (animal_id/date/cluster_id/origin_group)
  el <- as_tibble(evening_locations)
  meta <- as_tibble(animal_metadata)
  
  # Join metadata into evening locations if origin_group missing in el
  if (!"origin_group" %in% names(el)) {
    el <- el %>% left_join(meta %>% select(animal_id, origin_group), by = "animal_id")
  }
  
  # Make sure join succeeded
  if (!"origin_group" %in% names(el)) stop("origin_group not present after joining metadata. Check inputs.")
  
  # cluster-day details
  cluster_day_details <- el %>%
    group_by(date, origin_group, cluster_id) %>%
    summarise(
      n_in_cluster = n_distinct(animal_id),
      animal_ids = list(sort(unique(animal_id))),
      .groups = "drop"
    )
  
  # group-day summary (how many clusters per group per date)
  group_day_summary <- cluster_day_details %>%
    group_by(date, origin_group) %>%
    summarise(
      n_clusters = n(),
      total_group_on_date = sum(n_in_cluster),
      clusters = list(tibble(cluster_id = cluster_id, n_in_cluster = n_in_cluster, animal_ids = animal_ids)),
      .groups = "drop"
    )
  
  # keep only split dates
  split_dates <- group_day_summary %>% filter(n_clusters > 1) %>% arrange(origin_group, date)
  if (nrow(split_dates) == 0) {
    message("No within-group splits detected")
    return(list(timestamps_df = tibble(), summary_df = tibble()))
  }
  
  # group consecutive split dates into sessions
  split_sessions <- split_dates %>%
    group_by(origin_group) %>%
    arrange(date) %>%
    mutate(
      days_since_last = as.numeric(difftime(date, lag(date), units = "days")),
      is_continuous = is.na(days_since_last) | days_since_last <= continuity_gap_days,
      split_number = cumsum(!is_continuous) + 1,
      split_session = paste0(origin_group, "_split_", split_number)
    ) %>%
    ungroup()
  
  # expand cluster-level rows for timestamps
  timestamps_df <- cluster_day_details %>%
    inner_join(split_sessions %>% select(date, origin_group, n_clusters, total_group_on_date, split_number, split_session),
               by = c("date", "origin_group")) %>%
    mutate(prop_of_group = ifelse(total_group_on_date > 0, n_in_cluster / total_group_on_date, NA_real_)) %>%
    arrange(origin_group, split_number, date, cluster_id) %>%
    select(date, origin_group, split_session, split_number, n_clusters, cluster_id, n_in_cluster, prop_of_group, animal_ids, total_group_on_date)
  
  # compute mean cluster size per cluster across session for stable ranking
  cluster_mean_sizes <- timestamps_df %>%
    group_by(origin_group, split_number, split_session, cluster_id) %>%
    summarise(
      mean_size = mean(n_in_cluster, na.rm = TRUE),
      max_size = max(n_in_cluster, na.rm = TRUE),
      members = list(sort(unique(unlist(animal_ids)))),
      .groups = "drop"
    )
  
  cluster_ranked <- cluster_mean_sizes %>%
    group_by(origin_group, split_number, split_session) %>%
    arrange(desc(mean_size), cluster_id) %>%
    mutate(cluster_rank = row_number(), cluster_label = paste0(origin_group, "_", cluster_rank)) %>%
    ungroup()
  
  timestamps_df <- timestamps_df %>%
    left_join(cluster_ranked %>% select(origin_group, split_number, split_session, cluster_id, cluster_rank, cluster_label, mean_size),
              by = c("origin_group", "split_number", "split_session", "cluster_id")) %>%
    arrange(origin_group, split_number, cluster_rank, date)
  
  # summary per split session, including clusters tibble
  summary_df <- timestamps_df %>%
    group_by(origin_group, split_number, split_session) %>%
    summarise(
      start_date = min(date),
      end_date = max(date),
      duration_days = as.numeric(difftime(max(date), min(date), units = "days")) + 1,
      n_nights = n_distinct(date),
      avg_n_clusters = round(mean(n_clusters), 2),
      avg_cluster_size = round(mean(n_in_cluster), 2),
      max_cluster_size = max(n_in_cluster),
      min_cluster_size = min(n_in_cluster),
      total_unique_individuals = length(unique(unlist(animal_ids))),
      participants = list(sort(unique(unlist(animal_ids)))),
      clusters = list(
        timestamps_df %>%
          filter(origin_group == first(origin_group), split_number == first(split_number), split_session == first(split_session)) %>%
          distinct(cluster_id, cluster_label, cluster_rank, mean_size, n_in_cluster, animal_ids) %>%
          arrange(cluster_rank)
      ),
      .groups = "drop"
    ) %>%
    arrange(origin_group, start_date)
  
  # time since last split
  summary_df <- summary_df %>%
    group_by(origin_group) %>%
    arrange(start_date) %>%
    mutate(prev_end = lag(end_date), time_since_last_split_days = as.numeric(difftime(start_date, prev_end, units = "days"))) %>%
    ungroup() %>%
    select(origin_group, split_number, split_session, start_date, end_date, duration_days, n_nights, total_unique_individuals, participants, avg_n_clusters, avg_cluster_size, min_cluster_size, max_cluster_size, clusters, time_since_last_split_days)
  
  # Return
  return(list(timestamps_df = timestamps_df, summary_df = summary_df))
}

# -------------------------
# Per-group-per-date split flags (yes/no) and missing data
# -------------------------
group_day_split_flags <- function(evening_locations, animal_metadata,
                                  complete_dates = TRUE,
                                  evening_animal_id_col = "animal_id",
                                  date_col = "date",
                                  cluster_col = "cluster_id",
                                  metadata_animal_id_col = "animal_id",
                                  metadata_group_col = "origin_group") {
  # Normalize inputs
  norm <- normalize_inputs(evening_locations, animal_metadata,
                           evening_animal_id_col = evening_animal_id_col,
                           date_col = date_col,
                           cluster_col = cluster_col,
                           metadata_animal_id_col = metadata_animal_id_col,
                           metadata_group_col = metadata_group_col)
  el <- norm$evening_locations
  meta <- norm$animal_metadata
  
  # Ensure origin_group in el
  if (!"origin_group" %in% names(el)) {
    el <- el %>% left_join(meta %>% select(animal_id, origin_group), by = "animal_id")
  }
  if (!"origin_group" %in% names(el)) stop("origin_group not present after join. Check metadata.")
  
  group_day <- el %>%
    group_by(date, origin_group) %>%
    summarise(
      n_clusters = n_distinct(cluster_id),
      n_individuals = n_distinct(animal_id),
      clusters = list(sort(unique(cluster_id))),
      .groups = "drop"
    ) %>%
    mutate(split = if_else(n_clusters > 1, "yes", "no"))
  
  if (complete_dates) {
    all_dates <- seq.Date(min(el$date, na.rm = TRUE), max(el$date, na.rm = TRUE), by = "day")
    all_groups <- sort(unique(meta$origin_group))
    grid <- expand.grid(date = all_dates, origin_group = all_groups, stringsAsFactors = FALSE) %>% as_tibble()
    group_day <- grid %>%
      left_join(group_day, by = c("date", "origin_group")) %>%
      mutate(
        missing_data = if_else(is.na(n_clusters) & is.na(n_individuals), TRUE, FALSE),
        n_clusters = if_else(is.na(n_clusters), 0L, as.integer(n_clusters)),
        n_individuals = if_else(is.na(n_individuals), 0L, as.integer(n_individuals)),
        split = if_else(missing_data, NA_character_, split)
      )
  } else {
    group_day <- group_day %>% mutate(missing_data = FALSE)
  }
  
  return(group_day)
}

# -------------------------
# Flatten clusters to one-row-per-cluster per split session (sorted & renumbered)
# -------------------------
build_clusters_long_from_splits <- function(splits_named, group_col_name = "origin_group") {
  if (is.null(splits_named$summary_df) || nrow(splits_named$summary_df) == 0) stop("splits_named$summary_df missing or empty")
  
  flat <- splits_named$summary_df %>%
    # ensure column exists
    mutate("{group_col_name}" := origin_group) %>%
    select(all_of(group_col_name), split_number, split_session, clusters) %>%
    mutate(clusters_df = map(clusters, function(x){
      df <- as_tibble(x)
      if ("members" %in% names(df)) {
        df <- df %>% mutate(members = map(members, ~ if (is.factor(.x)) as.character(.x) else .x))
      } else {
        df <- df %>% mutate(members = map(seq_len(nrow(df)), ~ character(0)))
      }
      df %>% mutate(n_individuals = map_int(members, length))
    })) %>%
    select(-clusters) %>%
    unnest(clusters_df) %>%
    group_by(!!sym(group_col_name), split_number) %>%
    arrange(desc(n_individuals), cluster_id) %>%
    mutate(cluster_internal_number = row_number(),
           cluster_label_internal = paste0(!!sym(group_col_name), "_", cluster_internal_number)) %>%
    ungroup() %>%
    select(!!sym(group_col_name), split_number, split_session, cluster_label_internal, n_individuals, cluster_internal_number)
  
  return(flat)
}

# -------------------------
# Robust plotting: clusters per group over time (auto-detect group col; accept missing flags)
# -------------------------
plot_group_clusters_timeline <- function(splits_named,
                                         group_day_flags,
                                         group_col = NULL,
                                         show_single_cluster_days = TRUE,
                                         offset_scale = 0.18,
                                         date_limits = NULL) {
  
  # Extract dataframes
  timestamps <- as_tibble(splits_named$timestamps_df)
  flags <- as_tibble(group_day_flags)
  
  cat("DEBUG: timestamps nrows =", nrow(timestamps), "\n")
  cat("DEBUG: flags nrows =", nrow(flags), "\n")
  cat("DEBUG: timestamps cols =", paste(names(timestamps), collapse = ", "), "\n")
  cat("DEBUG: flags cols =", paste(names(flags), collapse = ", "), "\n")
  
  if (nrow(timestamps) == 0 && nrow(flags) == 0) {
    stop("Both timestamps and flags are empty")
  }
  
  # Find group column name (auto-detect)
  group_col_name <- NA_character_
  for (col in c("origin_group", "group_id", "group")) {
    if (col %in% names(timestamps) || col %in% names(flags)) {
      group_col_name <- col
      break
    }
  }
  
  if (is.na(group_col_name)) {
    stop("Cannot find group column.  Available in timestamps:  ", 
         paste(names(timestamps), collapse = ", "), 
         " | Available in flags: ", 
         paste(names(flags), collapse = ", "))
  }
  
  cat("DEBUG: Using group column name:", group_col_name, "\n")
  
  # Rename to standard name in both dataframes
  if (group_col_name != "origin_group") {
    if (group_col_name %in% names(timestamps)) {
      timestamps <- timestamps %>% rename(origin_group = !!sym(group_col_name))
    }
    if (group_col_name %in% names(flags)) {
      flags <- flags %>% rename(origin_group = !! sym(group_col_name))
    }
  }
  
  # Convert to character to avoid factor issues
  if ("origin_group" %in% names(timestamps)) {
    timestamps <- timestamps %>% mutate(origin_group = as.character(origin_group))
  }
  if ("origin_group" %in% names(flags)) {
    flags <- flags %>% mutate(origin_group = as.character(origin_group))
  }
  
  # Get all groups
  all_groups <- unique(c(
    if ("origin_group" %in% names(timestamps)) unique(timestamps$origin_group) else character(),
    if ("origin_group" %in% names(flags)) unique(flags$origin_group) else character()
  ))
  all_groups <- sort(na.omit(all_groups))
  
  cat("DEBUG: all_groups =", paste(all_groups, collapse = ", "), "\n")
  
  if (length(all_groups) == 0) {
    stop("No groups found in timestamps or flags")
  }
  
  # Create group index
  group_index <- tibble(origin_group = all_groups, group_y = seq_along(all_groups))
  cat("DEBUG: group_index created, nrows =", nrow(group_index), "\n")
  
  # Process timestamps
  if (nrow(timestamps) > 0) {
    cat("DEBUG: Processing timestamps...\n")
    timestamps <- timestamps %>%
      mutate(origin_group = as.character(origin_group)) %>%
      left_join(group_index, by = "origin_group")
    
    cat("DEBUG:  After join, timestamps has group_y?", "group_y" %in% names(timestamps), "\n")
    
    if (!"group_y" %in% names(timestamps)) {
      stop("group_y not found in timestamps after join.  Check that origin_group values match between timestamps and group_index.")
    }
    
    timestamps <- timestamps %>%
      mutate(
        cluster_rank = as.integer(coalesce(cluster_rank, 1)),
        n_in_cluster = as.integer(coalesce(n_in_cluster, 1))
      ) %>%
      group_by(origin_group, date) %>%
      mutate(n_clusters_day = first(n_clusters)) %>%
      ungroup() %>%
      mutate(
        y_offset = (cluster_rank - (coalesce(n_clusters_day, 1) + 1) / 2) * offset_scale,
        y = group_y + y_offset
      )
    
    cat("DEBUG: timestamps processed.  nrows =", nrow(timestamps), "\n")
  } else {
    timestamps <- tibble()
  }
  
  # Process flags
  if (nrow(flags) > 0 && show_single_cluster_days) {
    cat("DEBUG: Processing flags...\n")
    single_points <- flags %>%
      filter(! coalesce(missing_data, FALSE), n_clusters == 1) %>%
      mutate(origin_group = as.character(origin_group)) %>%
      left_join(group_index, by = "origin_group") %>%
      mutate(
        y = group_y,
        n_in_cluster = n_individuals
      )
    
    cat("DEBUG:  single_points created, nrows =", nrow(single_points), "\n")
  } else {
    single_points <- tibble()
  }
  
  # Build plot
  cat("DEBUG: Building ggplot...\n")
  
  p <- ggplot() +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "Within-group cluster timeline (sub-clusters shown as offsets)",
      subtitle = "Points = cluster presence; size = number of individuals; color = origin group",
      x = "Date",
      y = "Origin Group"
    ) +
    scale_y_continuous(
      breaks = group_index$group_y,
      labels = group_index$origin_group,
      expand = expansion(add = c(0.5, 0.5))
    ) +
    guides(size = guide_legend("Cluster size\n(n individuals)"))
  
  # Add split clusters
  if (nrow(timestamps) > 0) {
    cat("DEBUG: Adding split clusters to plot...\n")
    p <- p + geom_point(
      data = timestamps,
      aes(x = date, y = y, color = origin_group, size = pmax(1, n_in_cluster)),
      alpha = 0.8
    )
  }
  
  # Add single-cluster points
  if (nrow(single_points) > 0) {
    cat("DEBUG: Adding single-cluster points to plot...\n")
    p <- p + geom_point(
      data = single_points,
      aes(x = date, y = y, color = origin_group, size = pmax(1, n_in_cluster)),
      alpha = 0.5,
      shape = 16
    )
  }
  
  # Add group guide lines
  if (nrow(group_index) > 0) {
    p <- p + geom_hline(
      data = group_index,
      aes(yintercept = group_y),
      color = "grey90",
      size = 0.3
    )
  }
  
  # Colors
  if (! is.null(group_col) && is.character(group_col) && !is.null(names(group_col))) {
    if (all(all_groups %in% names(group_col))) {
      p <- p + scale_color_manual(values = group_col)
    } else {
      p <- p + scale_color_brewer(palette = "Set2")
    }
  } else {
    p <- p + scale_color_brewer(palette = "Set2")
  }
  
  # Date limits
  if (!is.null(date_limits) && length(date_limits) >= 2) {
    p <- p + xlim(as.Date(date_limits[1]), as.Date(date_limits[2]))
  }
  
  cat("DEBUG: Plot created successfully!\n")
  p
}

# -------------------------
# Example usage helper (convenience)
# -------------------------
run_example_pipeline <- function(evening_locations, metadata,
                                 evening_animal_id_col = "animal_id",
                                 date_col = "date",
                                 cluster_col = "cluster_id",
                                 metadata_animal_id_col = "individual_local_identifier",
                                 metadata_group_col = "group_id",
                                 continuity_gap_days = 1,
                                 group_col = NULL) {
  prepped <- normalize_inputs(evening_locations, metadata,
                              evening_animal_id_col = evening_animal_id_col,
                              date_col = date_col,
                              cluster_col = cluster_col,
                              metadata_animal_id_col = metadata_animal_id_col,
                              metadata_group_col = metadata_group_col)
  el <- prepped$evening_locations
  meta <- prepped$animal_metadata
  
  splits_named <- identify_within_group_splits_numbered(el, meta, continuity_gap_days = continuity_gap_days)
  flags <- group_day_split_flags(el, meta, complete_dates = TRUE)
  clusters_long <- if (is.data.frame(splits_named$summary_df) && nrow(splits_named$summary_df) > 0) build_clusters_long_from_splits(splits_named) else tibble()
  p <- plot_group_clusters_timeline(splits_named, flags, evening_locations = el, animal_metadata = meta, group_col = group_col)
  
  return(list(splits_named = splits_named, flags = flags, clusters_long = clusters_long, plot = p))
}


# Additional plot:  cluster count per group over time
# Shows a clear line/bar chart of how many clusters each group had on each date

plot_cluster_count_over_time <- function(group_day_flags,
                                         group_col = NULL,
                                         plot_type = "line") {
  
  # Extract relevant columns
  plot_data <- group_day_flags %>%
    filter(! missing_data) %>%
    select(date, origin_group, n_clusters) %>%
    # Additional plot:  cluster count per group over time
# Shows a clear line/bar chart of how many clusters each group had on each date

plot_cluster_count_over_time <- function(group_day_flags,
                                         group_col = NULL,
                                         plot_type = "line") {
  
  # Extract relevant columns
  plot_data <- group_day_flags %>%
    dplyr::filter(! missing_data) %>%
    dplyr::select(date, origin_group, n_clusters) %>%
    dplyr::mutate(origin_group = as.character(origin_group))
  
  if (nrow(plot_data) == 0) {
    stop("No non-missing data to plot")
  }
  
  # Ensure date is Date class
  plot_data <- plot_data %>%
    mutate(date = as.Date(date))
  
  # Base plot
  if (plot_type == "line") {
    p <- ggplot(plot_data, aes(x = date, y = n_clusters, color = origin_group, group = origin_group)) +
      geom_line(size = 1, alpha = 0.8) +
      geom_point(size = 2.5, alpha = 0.9) +
      labs(
        title = "Number of clusters per group over time",
        subtitle = "Shows when groups split into multiple clusters",
        x = "Date",
        y = "Number of clusters",
        color = "Origin Group"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(size = 13, face = "bold"),
        legend.position = "right"
      ) +
      scale_y_continuous(breaks = seq(0, max(plot_data$n_clusters) + 1, by = 1), limits = c(0, max(plot_data$n_clusters) + 0.5)) +
      scale_x_date(date_breaks = "7 days", date_labels = "%Y-%m-%d")
  } 
  else if (plot_type == "bar") {
    p <- ggplot(plot_data, aes(x = date, y = n_clusters, fill = origin_group)) +
      geom_col(position = "dodge", alpha = 0.8) +
      labs(
        title = "Number of clusters per group over time",
        subtitle = "Shows when groups split into multiple clusters",
        x = "Date",
        y = "Number of clusters",
        fill = "Origin Group"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        plot.title = element_text(size = 13, face = "bold"),
        legend.position = "right"
      ) +
      scale_y_continuous(breaks = seq(0, max(plot_data$n_clusters) + 1, by = 1)) +
      scale_x_date(date_breaks = "7 days", date_labels = "%Y-%m-%d")
  }
  else if (plot_type == "facet") {
    p <- ggplot(plot_data, aes(x = date, y = n_clusters, fill = origin_group)) +
      geom_col(alpha = 0.8) +
      facet_wrap(~origin_group, ncol = 2) +
      labs(
        title = "Number of clusters per group over time (faceted by group)",
        subtitle = "Each panel shows one group's cluster count over time",
        x = "Date",
        y = "Number of clusters",
        fill = "Origin Group"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(size = 13, face = "bold"),
        legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold")
      ) +
      scale_y_continuous(breaks = seq(0, max(plot_data$n_clusters) + 1, by = 1)) +
      scale_x_date(date_breaks = "14 days", date_labels = "%m-%d")
  }
  else {
    stop("plot_type must be 'line', 'bar', or 'facet'")
  }
  
  # Apply colors if provided
  if (!is.null(group_col) && is.character(group_col) && !is.null(names(group_col))) {
    if (plot_type == "line") {
      p <- p + scale_color_manual(values = group_col)
    } else {
      p <- p + scale_fill_manual(values = group_col)
    }
  } else {
    if (plot_type == "line") {
      p <- p + scale_color_brewer(palette = "Set2")
    } else {
      p <- p + scale_fill_brewer(palette = "Set2")
    }
  }
  
  # Add horizontal line at y=1 to mark single-cluster days
  p <- p + geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", size = 0.5, alpha = 0.5)
  
  p
}

# ==========================================
# EXAMPLE USAGE: 
# ==========================================
# Three versions of the same data:
#
# 1. Line plot (best for trends):
p1 <- plot_cluster_count_over_time(flags, group_col = NULL, plot_type = "line")
print(p1)
#
# 2. Bar plot (best for comparing counts):
p2 <- plot_cluster_count_over_time(flags, group_col = NULL, plot_type = "bar")
print(p2)
#
# 3. Faceted bar plot (best for individual group focus):
p3 <- plot_cluster_count_over_time(flags, group_col = NULL, plot_type = "facet")
print(p3)
#
# With custom colors:
# group_colors <- c("Bronze" = "#B87333", "Silver" = "#C0C0C0", "Gold" = "#FFD700")
# p1_colored <- plot_cluster_count_over_time(flags, group_col = group_colors, plot_type = "line")
# print(p1_colored)
#
# Make interactive:
# plotly::ggplotly(p1)
# plotly:: ggplotly(p2)
# plotly::ggplotly(p3)mutate(origin_group = character(origin_group))
  
  if (nrow(plot_data) == 0) {
    stop("No non-missing data to plot")
  }
  
  # Ensure date is Date class
  plot_data <- plot_data %>%
    mutate(date = as.Date(date))
  
  # Base plot
  if (plot_type == "line") {
    p <- ggplot(plot_data, aes(x = date, y = n_clusters, color = origin_group, group = origin_group)) +
      geom_line(size = 1, alpha = 0.8) +
      geom_point(size = 2.5, alpha = 0.9) +
      labs(
        title = "Number of clusters per group over time",
        subtitle = "Shows when groups split into multiple clusters",
        x = "Date",
        y = "Number of clusters",
        color = "Origin Group"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(size = 13, face = "bold"),
        legend.position = "right"
      ) +
      scale_y_continuous(breaks = seq(0, max(plot_data$n_clusters) + 1, by = 1), limits = c(0, max(plot_data$n_clusters) + 0.5)) +
      scale_x_date(date_breaks = "7 days", date_labels = "%Y-%m-%d")
  } 
  else if (plot_type == "bar") {
    p <- ggplot(plot_data, aes(x = date, y = n_clusters, fill = origin_group)) +
      geom_col(position = "dodge", alpha = 0.8) +
      labs(
        title = "Number of clusters per group over time",
        subtitle = "Shows when groups split into multiple clusters",
        x = "Date",
        y = "Number of clusters",
        fill = "Origin Group"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        plot.title = element_text(size = 13, face = "bold"),
        legend.position = "right"
      ) +
      scale_y_continuous(breaks = seq(0, max(plot_data$n_clusters) + 1, by = 1)) +
      scale_x_date(date_breaks = "7 days", date_labels = "%Y-%m-%d")
  }
  else if (plot_type == "facet") {
    p <- ggplot(plot_data, aes(x = date, y = n_clusters, fill = origin_group)) +
      geom_col(alpha = 0.8) +
      facet_wrap(~origin_group, ncol = 2) +
      labs(
        title = "Number of clusters per group over time (faceted by group)",
        subtitle = "Each panel shows one group's cluster count over time",
        x = "Date",
        y = "Number of clusters",
        fill = "Origin Group"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(size = 13, face = "bold"),
        legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold")
      ) +
      scale_y_continuous(breaks = seq(0, max(plot_data$n_clusters) + 1, by = 1)) +
      scale_x_date(date_breaks = "14 days", date_labels = "%m-%d")
  }
  else {
    stop("plot_type must be 'line', 'bar', or 'facet'")
  }
  
  # Apply colors if provided
  if (!is.null(group_col) && is.character(group_col) && !is.null(names(group_col))) {
    if (plot_type == "line") {
      p <- p + scale_color_manual(values = group_col)
    } else {
      p <- p + scale_fill_manual(values = group_col)
    }
  } else {
    if (plot_type == "line") {
      p <- p + scale_color_brewer(palette = "Set2")
    } else {
      p <- p + scale_fill_brewer(palette = "Set2")
    }
  }
  
  # Add horizontal line at y=1 to mark single-cluster days
  p <- p + geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", size = 0.5, alpha = 0.5)
  
  p
}

# ==========================================
# EXAMPLE USAGE: 
# ==========================================
# Three versions of the same data:
#
# 1. Line plot (best for trends):
p1 <- plot_cluster_count_over_time(flags, group_col = NULL, plot_type = "line")
print(p1)

# 2. Bar plot (best for comparing counts):
p2 <- plot_cluster_count_over_time(flags, group_col = NULL, plot_type = "bar")
print(p2)

# 3. Faceted bar plot (best for individual group focus):
p3 <- plot_cluster_count_over_time(flags, group_col = NULL, plot_type = "facet")
print(p3)

#With custom colors:
group_colors <- c("Bronze" = "#B87333", "Silver" = "#C0C0C0", "Gold" = "#FFD700")
p1_colored <- plot_cluster_count_over_time(flags, group_col = group_colors, plot_type = "line")
print(p1_colored)

#Make interactive:
plotly::ggplotly(p1)
plotly:: ggplotly(p2)
plotly::ggplotly(p3)


