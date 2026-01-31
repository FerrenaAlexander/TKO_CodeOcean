# =============================================================================
# Master Script to Run All Figure Scripts
# =============================================================================
# This script runs all 8 figure generation scripts in sequence
# =============================================================================

# Set working directory to script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Source configuration
source("config.R")

# Define scripts to run
scripts <- c(
  "codeocean_fig1_clean.R",
  "codeocean_fig2_clean.R",
  "codeocean_fig3_clean.R",
  "codeocean_fig4_clean.R",
  "codeocean_fig5_clean.R",
  "codeocean_fig6_clean.R",
  "codeocean_fig7_clean.R",
  "codeocean_fig8_clean.R"
)

# Track execution time and status
results <- data.frame(
  Script = scripts,
  Status = character(length(scripts)),
  Time_seconds = numeric(length(scripts)),
  stringsAsFactors = FALSE
)

cat("\n")
cat("================================================================\n")
cat("Running All Figure Scripts\n")
cat("================================================================\n")
cat("\n")

# Run each script
for (i in seq_along(scripts)) {
  script <- scripts[i]
  cat(sprintf("Running %s (%d/%d)...\n", script, i, length(scripts)))
  
  start_time <- Sys.time()
  
  tryCatch({
    source(script)
    results$Status[i] <- "SUCCESS"
  }, error = function(e) {
    results$Status[i] <- paste0("ERROR: ", e$message)
  })
  
  end_time <- Sys.time()
  results$Time_seconds[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat(sprintf("  Status: %s (%.1f seconds)\n\n", results$Status[i], results$Time_seconds[i]))
}

# Print summary
cat("\n")
cat("================================================================\n")
cat("Execution Summary\n")
cat("================================================================\n")
print(results)
cat("\n")

# Save summary
write.csv(results, file.path(BASE_OUTPUT_DIR, "execution_summary.csv"), row.names = FALSE)
cat("Summary saved to:", file.path(BASE_OUTPUT_DIR, "execution_summary.csv"), "\n\n")

# Calculate total time
total_time <- sum(results$Time_seconds)
cat(sprintf("Total execution time: %.1f seconds (%.1f minutes)\n", total_time, total_time/60))

# Count successes and failures
n_success <- sum(results$Status == "SUCCESS")
n_fail <- sum(results$Status != "SUCCESS")
cat(sprintf("Successful: %d/%d\n", n_success, length(scripts)))
cat(sprintf("Failed: %d/%d\n", n_fail, length(scripts)))

if (n_fail > 0) {
  cat("\nFailed scripts:\n")
  print(results[results$Status != "SUCCESS", ])
}

cat("\n")
cat("================================================================\n")
cat("All scripts completed!\n")
cat("================================================================\n")
