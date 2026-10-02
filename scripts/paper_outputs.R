## ============================================================================
## paper_outputs.R -- every figure and table in the paper, from the fitted chains
##
##   source("scripts/paper_outputs.R")
##   check_chains("wt")        # supplementary traces + convergence table
##   t1 <- table1()            # Table 1
##
## Reads the curated chains in outputs/final: strains/ holds <strain>_<n>.rds and
## simulations/ holds the joined simulation chains <scenario>_c_<iterations>.rds
## (join_chains.R). Nothing here fits anything -- fit_model.R does that -- so this
## script is safe to run repeatedly and everything it produces is derived from files
## already on disk.
##
## Sections
##   helpers       burn-in detection, chain file resolution
##   diagnostics   check_chains: traces, split-Rhat, ESS   (supplementary material)
##                 plot_chains: one quantity faceted across the P/Q/R simulation runs
##   draws         combine_chains: pooled post-burn-in draws for everything below
##   tables        table1
##   figures       figure3 (posterior predictive), figure4 (kinetics and intervals),
##                 figure_prior_posterior (supplementary: prior against posterior)
## ============================================================================

library(ggplot2)

## ---------------------------------------------------------------------------
## Opaque backgrounds on every saved figure.
##
## ggplot2's ggsave inherits the plot's background, and a figure assembled by
## ggpubr::ggarrange has a TRANSPARENT margin -- the strip carrying the row labels in
## Figure 3, for instance, came out at alpha 0.19. It looks correct on a white page and
## renders as a black bar in any dark-background viewer, and behaves unpredictably once
## placed in Word or converted to TIFF for submission. Nothing in the manuscript wants a
## transparent background.
##
## Overriding the default here rather than at each of the sixteen call sites keeps it
## uniform and applies to figures added later. An explicit bg = in a call still wins.
## ---------------------------------------------------------------------------
ggsave <- function(..., bg = "white") ggplot2::ggsave(..., bg = bg)

## ---------------------------------------------------------------------------
## Model form -- which fitted model every figure and table reads
##
##   use_form("thresh_tdist")   # the reported fits (default)
##   use_form("incidence")      # the hybrid: incidence indexed on crossing vl_thresh
##   save_paper_outputs(form = "incidence")
##
## Every function defaults to dir = STRAIN_DIR and writes through FIG_DIR / TAB_DIR, so
## switching the form redirects inputs and outputs together and the two sets never overwrite
## each other. The individual trajectory is identical under both forms -- growth from vl_thresh
## at rate a, peak, decay at rate b, onset T after crossing vl_thresh -- so Figure 4 and
## Tables 1, 3 and 4 need no form-specific code; only the predicted counts (Figure 3, the phase
## check) and the meaning of the incidence curve (Figure 5) differ, and those read the form
## from the fitted chain's x$form.
##
## The simulation-study outputs and the prior sensitivity analysis exist for thresh_tdist only.
## ---------------------------------------------------------------------------
## `wide` holds the prior-sensitivity refits -- the same data under widened priors -- which
## figure_prior_sensitivity() compares against `strains`. It belongs here rather than as a
## literal default on that function, because the two forms have separate sets and a default
## pointing at the thresh_tdist directory would silently compare an incidence fit against a
## thresh_tdist one, or, once that directory no longer exists, report "no widened fit" for
## every variant and draw the figure with only half its rows.
FORM_DIRS <- list(
  thresh_tdist = c(strains = "outputs/final/strains",
                   wide    = "outputs/final/strains_wide_priors",
                   figures = "outputs/figures", tables = "outputs/tables"),
  incidence    = c(strains = "outputs/final/strains_incidence",
                   ## Curated by curate_widened.R, which checks the 2x widening and the mode
                   ## and writes x$analysis_from onto each fit. Alpha is a splice of alpha_c
                   ## and alpha_c3 (join_branch.R), so it is longer than the other three.
                   wide    = "outputs/final/strains_widened",
                   figures = "outputs/figures/incidence", tables = "outputs/tables/incidence"))

use_form <- function(form = c("thresh_tdist", "incidence")) {
  form <- match.arg(form)
  d <- FORM_DIRS[[form]]
  assign("FORM",       form,           envir = globalenv())
  assign("STRAIN_DIR", d[["strains"]], envir = globalenv())
  assign("WIDE_DIR",   d[["wide"]],    envir = globalenv())
  assign("FIG_DIR",    d[["figures"]], envir = globalenv())
  assign("TAB_DIR",    d[["tables"]],  envir = globalenv())
  message(sprintf("form %s: chains from %s, figures to %s, tables to %s",
                  form, d[["strains"]], d[["figures"]], d[["tables"]]))
  invisible(form)
}

## First draw the analysis uses. Curated chains may record it as x$analysis_from -- the hybrid
## chains do, because their kernel was fixed at a known iteration and the draws before it come
## from an adapting proposal. Otherwise the log-posterior burn-in heuristic below.
analysis_start <- function(o) {
  if (!is.null(o$x$analysis_from)) as.integer(o$x$analysis_from)
  else as.integer(burnin(o, plot = FALSE))
}

## The model form a fitted chain used; objects written before x$form existed are thresh_tdist.
fit_form <- function(x) if (is.null(x$form)) "thresh_tdist" else x$form

## ---------------------------------------------------------------------------
## on_natural_scale -- undo the sampler's parameterisation
##
## Chains fitted with x$log_vthresh hold log(vl_thresh) in inc3, because sampling it on the log
## scale removes the hard bound at zero. Nothing outside the sampler wants that: Table 1 reports
## a viral load, fig4_draw compares inc3 against the peak, and the recovery machinery compares it
## against generating values stored on the natural scale. Converting once, wherever draws are
## read, keeps every downstream function ignorant of the parameterisation.
##
## Applied to a matrix of draws or a single named vector. Chains written before the
## reparameterisation have no flag and pass through untouched.
## ---------------------------------------------------------------------------
on_natural_scale <- function(d, x) {
  if (!isTRUE(x$log_vthresh)) return(d)
  if (is.null(dim(d))) { if ("inc3" %in% names(d)) d[["inc3"]] <- exp(d[["inc3"]]); return(d) }
  if ("inc3" %in% colnames(d)) d[, "inc3"] <- exp(d[, "inc3"])
  d
}

## Both kernels, compiled once per session. ct_sim_cached dispatches on form.
load_kernels <- function() {
  if (!exists("p_array_func")) source("scripts/model.R")
  if (!exists("ct_func2_cpp_thresh_tdist")) Rcpp::sourceCpp("scripts/simulate_onset.cpp")
  if (!exists("ct_func2_cpp_inc"))          Rcpp::sourceCpp("scripts/simulate_incidence.cpp")
}

burnin <- function(o, plot = TRUE, max_points = 4000) {
  fl <- max(which(!is.na(o$MCMC_posteriors)))
  lp <- o$MCMC_posteriors[1:fl, 1]
  
  ## first time the chain reaches the level it eventually settles at
  thresh <- median(lp[floor(fl / 2):fl])
  b      <- which(lp >= thresh)[1]
  
  if (plot) {
    ## thin for drawing only -- b is computed on the full chain
    keep <- unique(round(seq(1, fl, length.out = min(fl, max_points))))
    df   <- data.frame(iteration = keep, lp = lp[keep])
    
    p <- ggplot(df, aes(iteration, lp)) +
      annotate("rect", xmin = -Inf, xmax = b, ymin = -Inf, ymax = Inf,
               fill = "firebrick", alpha = 0.06) +
      geom_line(linewidth = 0.25, colour = "grey20") +
      geom_hline(yintercept = thresh, linetype = 3, colour = "grey45") +
      geom_vline(xintercept = b, colour = "firebrick", linewidth = 0.6) +
      scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.02))) +
      scale_y_continuous(labels = scales::comma) +
      labs(title = o$x$strain,
           subtitle = sprintf("burn-in %s of %s iterations (%.1f%%)",
                              format(b, big.mark = ","),
                              format(fl, big.mark = ","), 100 * b / fl),
           x = "Iteration", y = "Log posterior") +
      theme_bw(base_size = 12) +
      theme(panel.grid   = element_blank(),
            panel.border = element_blank(),
            axis.line    = element_line(colour = "black"),
            plot.title    = element_text(size = 13, face = "bold", colour = "black"),
            plot.subtitle = element_text(size = 11, colour = "grey30", margin = margin(b = 8)),
            axis.title.x  = element_text(size = 12, colour = "black", margin = margin(t = 8)),
            axis.title.y  = element_text(size = 12, colour = "black", margin = margin(r = 8)),
            axis.text     = element_text(size = 11, colour = "black"))
    print(p)
    attr(b, "plot") <- p
  }
  b
}

## ---------------------------------------------------------------------------
## chain_file / chain_tags -- resolve chains under either naming convention
##
## Two are in play. The runners write <strain>_<tag>_<iteration>.rds and replace the file
## as the chain grows, so the iteration has to be parsed and the highest taken. The curated
## set in outputs/final/strains is <strain>_<n>.rds with no iteration at all. Try the exact name
## first, then the iteration pattern.
## ---------------------------------------------------------------------------
## Resolve a `dir` argument that may be one directory for every strain, or a NAMED vector
## giving one per strain -- for when the variants are not all at the same stage, e.g. three
## curated in outputs/final/strains_incidence while omicron is still running and staged
## elsewhere. Returns a function of the strain. Unnamed multi-element vectors are rejected so
## a typo fails loudly rather than silently using one directory for everything.
dir_resolver <- function(dir, strains) {
  if (length(dir) <= 1L) return(function(s) dir)
  if (is.null(names(dir)) || !all(strains %in% names(dir)))
    stop("a per-strain `dir` must be named for every strain drawn: ",
         paste(setdiff(strains, names(dir)), collapse = ", "), " missing")
  miss <- dir[strains][!dir.exists(dir[strains])]
  if (length(miss)) stop("missing directory: ", paste(miss, collapse = ", "))
  function(s) unname(dir[[s]])
}

chain_file <- function(strain, tag, dir = STRAIN_DIR) {
  exact <- file.path(dir, sprintf("%s_%s.rds", strain, tag))
  if (file.exists(exact)) return(exact)
  f <- list.files(dir, sprintf("^%s_%s_[0-9]+\\.rds$", strain, tag), full.names = TRUE)
  if (!length(f)) return(NULL)
  f[which.max(as.numeric(sub(".*_([0-9]+)\\.rds$", "\\1", basename(f))))]
}

## Every chain for a strain in a directory, for when tags are not given. Anchored on the
## strain so "wt" cannot pick up another strain's files, and any trailing _<iteration> is
## stripped so several checkpoints of one chain collapse to a single tag.
chain_tags <- function(strain, dir = STRAIN_DIR) {
  f <- list.files(dir, sprintf("^%s_.*\\.rds$", strain))
  if (!length(f)) stop("no chains for ", strain, " in ", dir)
  sort(unique(sub("_[0-9]+$", "", sub("\\.rds$", "", sub(sprintf("^%s_", strain), "", f)))))
}

## ---------------------------------------------------------------------------
## check_chains -- do the parallel chains for one strain agree?
##
##   check_chains("wt")                                  # every chain in outputs/final/strains
##   check_chains("wt", pars = c("a_bar", "b_bar"))
##   check_chains("alpha", tags = c("Nig","Nij","Nik"), dir = "outputs")   # raw runs
##
## Prints a summary table, a convergence table and two figures, and returns them invisibly.
## The dashed vertical lines are each chain's burn-in, in that chain's own colour.
##
## The summary is the part that answers the mode question: `post` is the mean log
## posterior after burn-in, so chains sitting hundreds of units apart are in different
## basins, not merely unconverged.
## ---------------------------------------------------------------------------
check_chains <- function(strain,
                         tags       = NULL,
                         dir        = STRAIN_DIR,
                         pars       = c("a_bar", "b_bar", "vl_max_bar", "inc1"),
                         ## Convergence is tabulated for more parameters than are worth
                         ## plotting -- the sigmas and the nuisance parameters need an ESS
                         ## before any interval is quoted, but their traces say little.
                         diag_pars  = c("a_bar", "a_sigma", "b_bar", "b_sigma",
                                        "vl_max_bar", "vl_max_sigma",
                                        "inc1", "inc2", "inc3", "multiplier", "rdisp"),
                         max_points = 3000,
                         zoom       = 200,
                         ncol       = 2,
                         from       = c("burnin", "stationary", "start")) {

  from <- match.arg(from)
  if (is.null(tags)) tags <- chain_tags(strain, dir)

  ## Keep only the columns being drawn. Each checkpoint holds a copy of the data array;
  ## reading three and retaining them whole is needless.
  ch <- list()
  for (tag in tags) {
    f <- chain_file(strain, tag, dir)
    if (is.null(f)) { message("no chain '", tag, "' for ", strain, " in ", dir); next }
    o  <- readRDS(f)
    fl <- max(which(!is.na(o$MCMC_posteriors)))
    b  <- analysis_start(o)
    ch[[tag]] <- list(lp   = o$MCMC_posteriors[1:fl, 1],
                      out  = o$MCMC_output[1:fl,
                                           intersect(unique(c(pars, diag_pars)),
                                                     colnames(o$MCMC_output)),
                                           drop = FALSE],
                      fl   = fl,
                      b    = b,
                      ## ESS of the log posterior for this chain alone. A single number for
                      ## how much independent information 90,000 correlated draws carry --
                      ## and the one that makes a 6.5% acceptance rate concrete.
                      ess  = posterior::ess_bulk(o$MCMC_posteriors[b:fl, 1]),
                      ## Post-burn-in acceptance rate. This is what separates a chain that
                      ## mixes slowly because the posterior is awkward from one that mixes
                      ## slowly because its proposal is the wrong shape: a block MVN
                      ## proposal scaled 2.38^2/d should accept around 20-25%. Far below
                      ## that and the chain is standing still between rare large moves,
                      ## which is exactly what long smooth excursions look like.
                      acc  = mean(o$Acceptances[b:fl], na.rm = TRUE),
                      file = basename(f))
    rm(o); invisible(gc(verbose = FALSE))
  }
  if (!length(ch)) stop("no checkpoints found for ", strain)

  ## ---- summary -------------------------------------------------------------
  smry <- do.call(rbind, lapply(names(ch), function(nm) {
    z <- ch[[nm]]
    data.frame(chain = nm, file = z$file, iterations = z$fl, burnin = z$b,
               acc = round(z$acc, 3), ess_lp = round(z$ess),
               post = mean(z$lp[z$b:z$fl]), best = max(z$lp),
               stringsAsFactors = FALSE)
  }))
  smry$gap <- round(smry$post - max(smry$post), 1)   # 0 marks the best basin
  print(smry, row.names = FALSE, digits = 8)

  ## ---- convergence across chains -------------------------------------------
  ## Each chain is cut at its own burn-in, then all are truncated to a common length taken
  ## from the END of each, so a chain that settled late contributes only stationary draws
  ## and none is over-weighted. Split-Rhat above ~1.01 means the chains disagree relative to
  ## their own spread; ess_bulk is the effective sample size behind a posterior mean,
  ## ess_tail behind the interval edges.
  conv <- NULL
  if (length(ch) >= 2) {
    n  <- min(vapply(ch, function(z) z$fl - z$b + 1L, integer(1)))
    dk <- intersect(diag_pars, colnames(ch[[1]]$out))
    arr <- array(NA_real_, c(n, length(ch), length(dk) + 1L),
                 dimnames = list(NULL, names(ch), c("log posterior", dk)))
    for (j in seq_along(ch)) {
      z  <- ch[[j]]
      ix <- (z$fl - n + 1L):z$fl
      arr[, j, ] <- cbind(z$lp[ix], z$out[ix, dk, drop = FALSE])
    }
    conv <- as.data.frame(posterior::summarise_draws(
      posterior::as_draws_array(arr), "mean", "sd", "rhat", "ess_bulk", "ess_tail"))
    cat(sprintf("\nconvergence over the last %s post-burn-in draws of each chain\n",
                format(n, big.mark = ",")))
    print(conv, row.names = FALSE, digits = 4)
  }

  ## ---- shared pieces -------------------------------------------------------
  top   <- max(vapply(ch, function(z) max(z$lp), numeric(1)))
  first <- min(vapply(ch, `[[`, integer(1), "b"))

  ## Where the panels start. Drawing from iteration 1 puts the climb out of a flat
  ## incidence curve -- tens of thousands of log units -- on the same axis as the
  ## stationary phase, which compresses everything after burn-in into a horizontal line and
  ## forces y-labels so wide they collide. Starting at the earliest burn-in keeps the one
  ## interesting excursion (a chain that leaves late, as Neh does here) while putting the
  ## axis on the scale of the differences that matter. `from = "start"` restores the
  ## full view when the burn-in itself is what you want to look at.
  ##
  ## "stationary" starts at the LAST chain's burn-in, so every chain shown is past its own.
  ## This is the view to use when judging mixing: with "burnin" the y-scale of each panel is
  ## set by whichever chain made the largest pre-stationary excursion, which compresses the
  ## stationary phase and makes strains look more alike than they are.
  x0 <- switch(from,
               burnin     = first,
               stationary = max(vapply(ch, `[[`, integer(1), "b")),
               start      = 1L)

  ## A chain sitting in a different basin can be hundreds of units below the best. Clipping
  ## at a fixed `zoom` would push it off the panel and it would read as missing rather than
  ## as the finding. Widen so the worst chain is always on screen.
  zoom <- max(zoom, top - min(smry$post) + 20)

  bl  <- data.frame(chain = names(ch), b = vapply(ch, `[[`, integer(1), "b"))
  pal <- setNames(c("#1b6ca8", "#c0392b", "#1e8449", "#7d3c98",
                    "#b7791f", "#00838f")[seq_along(ch)], names(ch))

  ## Short x-labels ("50K"), because "100,000" is wide enough to be clipped at the panel
  ## edge on a default device.
  kfmt <- function(v) ifelse(abs(v) >= 1000, paste0(round(v / 1000), "K"), format(v))

  base <- theme_bw(base_size = 12) +
    theme(panel.grid    = element_blank(),
          panel.border  = element_blank(),
          axis.line     = element_line(colour = "black"),
          strip.background = element_blank(),
          strip.text    = element_text(size = 10, face = "bold", colour = "black",
                                       hjust = 0, margin = margin(b = 2)),
          legend.position = "top",
          legend.title  = element_blank(),
          legend.margin = margin(b = 0),
          plot.title    = element_text(size = 13, face = "bold", colour = "black"),
          plot.subtitle = element_text(size = 10, colour = "grey30", margin = margin(b = 6)),
          axis.text     = element_text(size = 9, colour = "black"))

  ## ---- traces --------------------------------------------------------------
  ## Burn-in is computed on the full chain above; thinning here is for drawing only.
  pk     <- intersect(pars, colnames(ch[[1]]$out))
  panels <- c("log posterior", pk)
  long <- do.call(rbind, lapply(names(ch), function(nm) {
    z    <- ch[[nm]]
    idx  <- x0:z$fl
    keep <- unique(round(seq(min(idx), z$fl, length.out = min(length(idx), max_points))))
    vals <- cbind(`log posterior` = z$lp[keep], z$out[keep, pk, drop = FALSE])
    data.frame(chain     = nm,
               iteration = keep,
               parameter = rep(colnames(vals), each = length(keep)),
               value     = as.vector(vals),
               stringsAsFactors = FALSE)
  }))
  long$parameter <- factor(long$parameter, levels = panels)

  ## The log posterior still spans the late-leaving chain's climb, which is worth seeing but
  ## would again dominate its own panel. Clip it to the same window as the zoom plot.
  ilp <- long$parameter == "log posterior"
  long$value[ilp] <- pmax(long$value[ilp], top - zoom)

  p_trace <- ggplot(long, aes(iteration, value, colour = chain)) +
    geom_vline(data = bl, aes(xintercept = b, colour = chain),
               linetype = 2, linewidth = 0.5, show.legend = FALSE) +
    ## Where chains agree they draw over one another and only the last is visible, which
    ## makes agreement indistinguishable from a missing chain. Alpha lets the overlap read
    ## as a blend instead.
    geom_line(linewidth = 0.3, alpha = 0.65) +
    ## Two columns by default. Stacked in one column the panels are short enough that a
    ## parameter moving across its whole plausible range looks like noise; the vertical
    ## resolution is what tells you whether the chains actually agree. Four parameters plus
    ## the log posterior leaves one empty cell -- pass three parameters for an exact 2x2.
    ## rec_labeller falls back to the plain name, which is what keeps the "log posterior"
    ## panel -- not a parameter, and not parseable as plotmath -- working alongside the rest.
    facet_wrap(~ parameter, scales = "free_y", ncol = ncol, strip.position = "top",
               labeller = rec_labeller()) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(labels = kfmt,
                       breaks = scales::pretty_breaks(n = if (ncol > 1) 4 else 6),
                       expand = expansion(mult = c(0.02, 0.03))) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = if (ncol > 1) 4 else 3),
                       labels = scales::comma) +
    labs(title = strain,
         subtitle = sprintf("from iteration %s; dashed lines mark each chain's burn-in",
                            format(x0, big.mark = ",")),
         x = "Iteration", y = NULL) +
    base

  ## ---- log posterior, on its own scale ---------------------------------------
  lpd <- long[ilp & long$iteration >= first, ]

  p_zoom <- ggplot(lpd, aes(iteration, value, colour = chain)) +
    geom_hline(yintercept = top, linetype = 3, colour = "grey45") +
    geom_vline(data = bl, aes(xintercept = b, colour = chain),
               linetype = 2, linewidth = 0.5, show.legend = FALSE) +
    geom_line(linewidth = 0.3, alpha = 0.65) +
    coord_cartesian(ylim = c(top - zoom, top)) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(labels = kfmt, breaks = scales::pretty_breaks(n = 6),
                       expand = expansion(mult = c(0.02, 0.03))) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = strain,
         subtitle = sprintf("post burn-in, within %s of the best value reached (dotted)",
                            format(round(zoom))),
         x = "Iteration", y = "Log posterior") +
    base

  print(p_trace); print(p_zoom)
  invisible(list(summary = smry, convergence = conv, trace = p_trace, zoom = p_zoom))
}

## ---------------------------------------------------------------------------
## plot_chains -- one quantity faceted across every run of a simulation family
##
##   plot_chains()                             # log-posterior, joined P/Q/R chains
##   plot_chains("a_sigma")                    # any MCMC_output column by name
##   plot_chains("b_sigma", letters = "R", chain = 1, dir = "outputs/old", thin = 1)
##
## Where check_chains compares parallel chains of one strain, this compares the same
## quantity across the runs of the P/Q/R simulation study (chain tag fixed, runs
## varying). Defaults to the joined chains (tag "c", from join_chains.R); the raw
## segments live in outputs/old under numeric tags. Uses the latest
## (highest-iteration) dump per run, so it is safe to call while runners are still
## replacing files. `param` is "posterior", "likelihood", or a column of MCMC_output;
## the preallocated NA tail of each matrix is dropped.
## ---------------------------------------------------------------------------
plot_chains <- function(param = "posterior", letters = c("P", "Q", "R"),
                        chain = "c", dir = "outputs/final/simulations",
                        thin = 10, ncol = 4) {
  files <- list.files(dir, pattern = sprintf("^[%s]\\d+_%s_\\d+\\.rds$",
                                             paste(letters, collapse = ""), chain))
  if (length(files) == 0) stop("No matching files in ", dir)
  run   <- sub("_\\d+\\.rds$", "", files)          # e.g. "P3_2"
  iters <- as.numeric(sub(".*_(\\d+)\\.rds$", "\\1", files))
  latest <- unlist(lapply(split(seq_along(files), run),
                          function(i) files[i][which.max(iters[i])]))

  df <- do.call(rbind, lapply(latest, function(f) {
    out <- readRDS(file.path(dir, f))
    v <- if (param == "posterior") {
      as.numeric(out$MCMC_posteriors)
    } else if (param == "likelihood") {
      as.numeric(out$MCMC_likelihoods)
    } else {
      if (!param %in% colnames(out$MCMC_output))
        stop("'", param, "' not found in ", f, ". Available: ",
             paste(colnames(out$MCMC_output), collapse = ", "))
      out$MCMC_output[, param]
    }
    v <- v[!is.na(v)]
    keep <- seq(1, length(v), by = thin)
    data.frame(run = sub("_\\d+\\.rds$", "", f),
               iteration = keep, value = v[keep])
  }))

  ggplot(df, aes(iteration, value)) +
    geom_line(linewidth = 0.2) +
    facet_wrap(~run, ncol = ncol, scales = "free") +
    labs(x = "Iteration", y = param) +
    theme_bw()
}

## ---------------------------------------------------------------------------
## combine_chains -- one chain's post-burn-in draws, with its continuations appended
##
##   combine_chains("delta", c("Nef", "Nhf"))
##
## Returns a plain matrix, all parameters. Use variant_diag() first to confirm the chains
## converged and agree; use this to build the pooled draws that feed the figures.
##
## `thin` keeps every nth draw. Delta's three combined chains run to ~870,000 rows across
## 43 columns, around 300 MB, and posterior summaries are unchanged by thinning to a tenth
## given ESS in the hundreds -- autocorrelation means those rows carry far less information
## than their number suggests.
## ---------------------------------------------------------------------------
combine_chains <- function(strain, tags, thin = 1L, dir = STRAIN_DIR) {
  out <- NULL; prev <- NULL
  for (k in seq_along(tags)) {
    f <- chain_file(strain, tags[k], dir)
    if (is.null(f)) stop("no chain '", tags[k], "' for ", strain, " in ", dir)
    o  <- readRDS(f)
    fl <- max(which(!is.na(o$MCMC_posteriors)))

    if (k == 1) {
      from <- analysis_start(o)
    } else {
      ## Same guard as variant_diag: a continuation only belongs to this chain if it ran
      ## the same kernel. Row 1 is its source's final state, so skip it.
      if (!isTRUE(all.equal(o$cov_matrix, prev, tolerance = 1e-10)))
        stop(tags[k], " uses a different proposal covariance than ", tags[k - 1])
      from <- 2L
    }

    out  <- rbind(out, on_natural_scale(o$MCMC_output[from:fl, , drop = FALSE], o$x))
    prev <- o$cov_matrix
    rm(o); invisible(gc(verbose = FALSE))
  }
  if (thin > 1L) out <- out[seq(1L, nrow(out), by = thin), , drop = FALSE]
  out
}

## ---------------------------------------------------------------------------
## table1 -- posterior medians and 95% credible intervals, all variants
##
##   t1 <- table1(); print(t1$raw); print(t1$derived)
##   write.csv(t1$raw, "outputs/table1.csv", row.names = FALSE)
##
## Replaces posterior_table() in results_figures.R, which reads through output_finder --
## unable to see the current filenames -- and summarises a single chain from iteration 1,
## burn-in included. This pools every chain for the variant, each cut at its own burn-in.
##
## `raw` is the parameters as fitted, which is what Table 1 reports. `derived` is the same
## posterior expressed in the units the Results discuss: growth and decay are shifted
## lognormals, so the population MEDIAN is `shift + exp(mu)` -- not the mean, which is
## `shift + exp(mu + sigma^2/2)` and runs well above the median when sigma is large. Using
## the mean is what made wild type's growth look implausible at a_sigma 1.33.
## ---------------------------------------------------------------------------
## The manuscript's row order and symbols, so output can be compared line by line with the
## submitted table. Equations 6 and 7 are the definitions; the caption states plainly that
## Table 1 carries population MEANS and SDs while Table 4 carries medians and IQRs, and the
## two differ substantially -- for wild type at a_sigma 1.5 the mean growth rate is 1.16
## against a median of 0.71.
T1_ROWS <- c(
  "a_mu_pop"     = "a_μ^pop",
  "a_sigma_pop"  = "a_σ^pop",
  "b_mu_pop"     = "b_μ^pop",
  "b_sigma_pop"  = "b_σ^pop",
  "vl_peak"      = "vl_peak",
  "vl_sigma"     = "vl_σ",
  "T_mu"         = "T_μ",
  "T_sigma"      = "T_σ",
  "v_thresh"     = "v_thresh")

table1 <- function(strains = c("wt", "alpha", "delta", "omicron"),
                   pars = c("a_bar", "a_sigma", "b_bar", "b_sigma",
                            "vl_max_bar", "vl_max_sigma", "inc1", "inc2", "inc3",
                            "multiplier", "rdisp"),
                   dir = STRAIN_DIR, thin = 10L, vl0 = 1.5, digits = 2,
                   latex = FALSE) {

  q <- c(0.5, 0.025, 0.975)
  ## One decimal setting for the whole table. The published Table 1 reports two everywhere,
  ## including the testing probabilities; this used to print those at three, so the two
  ## halves of one table disagreed on precision.
  fmt <- function(v, d) sprintf(paste0("%.", d, "f (%.", d, "f, %.", d, "f)"), v[1], v[2], v[3])

  ## Population mean and SD of a shifted lognormal, equations 6 and 7. The shift enters the
  ## mean and NOT the SD -- translating a distribution moves its centre without changing its
  ## spread. T_mu and T_sigma use the same expressions with no shift, as the manuscript
  ## specifies ("equation 6 with no-minimum term, and 7").
  pop_mean <- function(shift, mu, sg) shift + exp(mu + sg^2 / 2)
  pop_sd   <- function(mu, sg) sqrt((exp(sg^2) - 1) * exp(2 * mu + sg^2))

  ## Day-of-test probabilities from the discrete hazards test0..test5, exactly as
  ## p_array_func does it: read the hazards, write the probabilities somewhere else. The
  ## old posterior_table() wrote back into the same object, so from day 2 its survival term
  ## used already-converted probabilities rather than hazards -- worth checking the
  ## published Table 1 against this.
  test_probs <- function(d) {
    h  <- cbind(d[, paste0("test", 0:5), drop = FALSE], test6 = 1)
    pr <- h; surv <- rep(1, nrow(h))
    for (i in seq_len(ncol(h))) { pr[, i] <- surv * h[, i]; surv <- surv * (1 - h[, i]) }
    pr
  }

  keep <- NULL; raw <- NULL; der <- NULL; tst <- NULL; dnames <- NULL
  for (s in strains) {
    d <- do.call(rbind, lapply(chain_tags(s, dir),
                               function(t) combine_chains(s, t, thin = thin, dir = dir)))
    if (is.null(keep)) keep <- intersect(pars, colnames(d))
    raw <- cbind(raw, vapply(keep, function(p)
                    fmt(unname(quantile(d[, p], q)), digits), character(1)))

    ## Everything derived per draw, so each credible interval propagates the joint
    ## posterior rather than combining marginal endpoints -- growth, decay and peak are
    ## correlated, and duration depends on all three.
    ##
    ## Two levels of spread are in play and must not be confused. The POPULATION quantiles
    ## below describe variation between individuals at a given parameter value; the (2.5%,
    ## 97.5%) attached to each is posterior uncertainty about that quantity.
    ##
    ## Population percentiles rather than the SD. Growth and decay are shifted lognormals,
    ## so their SD is exp(mu + sigma^2/2) * sqrt(exp(sigma^2) - 1), which at wild type's
    ## a_sigma of 1.33 comes to 5.8 against a median of 1.58 -- a number describing the
    ## tail, not the population. The SD is reported too, for continuity with the published
    ## table, but the percentiles are the interpretable summary.
    lq  <- function(sh, mu, sg, p) sh + exp(mu + sg * qnorm(p))
    lsd <- function(mu, sg) exp(mu + sg^2 / 2) * sqrt(exp(sg^2) - 1)

    ab <- d[, "a_bar"]; asg <- d[, "a_sigma"]
    bb <- d[, "b_bar"]; bsg <- d[, "b_sigma"]
    vm <- d[, "vl_max_bar"]; vsg <- d[, "vl_max_sigma"]
    i1 <- d[, "inc1"];       i2  <- d[, "inc2"]; i3 <- d[, "inc3"]

    ## Names are the manuscript's own symbols, not the code's, so the printed table can be
    ## read against the paper without a translation step -- and so a renamed symbol has to be
    ## changed here too rather than drifting apart silently.
    ##
    ## vl_thresh is an ABSOLUTE viral load, on the same scale as the peak. The kernel
    ## compares it against the peak directly (simulate_onset.cpp:142, "if(vthresh > csv)")
    ## and measures the rise to peak as (csv - vthresh)/asv, both of which require it. Only
    ## the symbol changes here, from v_thresh; the definition in the manuscript is correct.
    ##
    ## Its posteriors -- 0.01 to 0.07 -- sit BELOW the 1.5 limit of detection, so the
    ## threshold is crossed before an individual becomes PCR-detectable. Symptoms are not:
    ## crossing the threshold starts a clock and onset follows T_i days later (eq 8), median
    ## 6.1 days for wild type, against about 0.94 days from threshold to detectability. Onset
    ## therefore falls ~5.2 days after detectability, which is what table4() reports as
    ## "Detectable -> symptoms", 5.10 (4.91-5.29). The estimate says immune activation begins
    ## at very low viral load, which is the mechanism the manuscript motivates at L136-137.
    ## The population MEDIAN is reported next to the mean for growth, decay and onset delay.
    ## All three are (shifted) lognormal, so their mean carries exp(sigma^2/2) and separates
    ## from the median as between-individual variation grows. At wild type's a_sigma of 1.33
    ## the mean growth rate is 3.11 against a median of 1.58 and sits at the 75th percentile
    ## of individuals -- so the mean alone makes wild type look like the fastest-growing
    ## variant when on the median the four agree (1.58, 1.60, 1.59, 1.20). Figure 4 and
    ## Table S6 already plot medians; printing only the mean here made the main text disagree
    ## with its own figure.
    ##
    ## The gap between the two is itself worth showing: it is a direct readout of the
    ## between-individual heterogeneity, which the mean and SD together obscure.
    ##
    ## No median row for vl_peak -- peak viral load is normal, so its mean and median are the
    ## same number and a second row would be noise.
    qs <- list(
      "a^pop_mu"     = pop_mean(0.5,  ab, asg),
      "a^pop_med"    = lq(0.5,  ab, asg, 0.5),
      "a^pop_sigma"  = pop_sd(ab, asg),
      "b^pop_mu"     = pop_mean(0.25, bb, bsg),
      "b^pop_med"    = lq(0.25, bb, bsg, 0.5),
      "b^pop_sigma"  = pop_sd(bb, bsg),
      "vlbar_peak"   = vm,
      "vl_sigma"     = vsg,
      ## Delay from threshold crossing to symptom onset, LogNormal(inc1, inc2) with no shift.
      "T^pop_mu"     = exp(i1 + i2^2 / 2),
      "T^pop_med"    = lq(0, i1, i2, 0.5),
      "T^pop_sigma"  = pop_sd(i1, i2),
      "vl_thresh"    = i3)

    if (is.null(dnames)) dnames <- names(qs)
    der <- cbind(der, vapply(qs, function(v) fmt(unname(quantile(v, q)), digits), character(1)))

    pr  <- test_probs(d)
    tst <- cbind(tst, vapply(colnames(pr), function(cn)
                     fmt(unname(quantile(pr[, cn], q)), digits), character(1)))
    rm(d, pr); invisible(gc(verbose = FALSE))
  }

  raw <- data.frame(parameter = keep, raw, stringsAsFactors = FALSE)
  der <- data.frame(quantity = dnames, der, stringsAsFactors = FALSE)
  tst <- data.frame(day = sprintf("P(symp -> test = %d)", 0:6), tst, stringsAsFactors = FALSE)
  colnames(raw) <- c("parameter", strains)
  colnames(der) <- c("parameter", strains)
  colnames(tst) <- c("parameter", strains)

  ## The manuscript's Table 1 is the kinetics/timing block stacked on the testing block, in
  ## that order. Returned assembled so it can be written out directly, with the two halves
  ## also available separately.
  t1 <- rbind(der, tst)
  rownames(raw) <- NULL; rownames(der) <- NULL; rownames(tst) <- NULL; rownames(t1) <- NULL

  tex <- table1_latex(t1, strains)
  if (latex) cat(tex, sep = "\n")
  list(table1 = t1, kinetics = der, testing = tst, raw = raw, latex = tex)
}

## LaTeX for Table 1, emitted so the manuscript is updated by pasting rather than by
## retyping sixty-four numbers. That matters more than convenience here: the submitted Table
## 1 still carries the pre-refit wild-type values (a^pop_mu 1.16, T_mu 12.84) while Tables 2
## and 3 carry the corrected ones, so the paper currently contradicts itself. Hand-copying is
## how that happened.
##
## Symbols follow NOTATION.md: vl_thresh rather than v_thresh, and the population
## superscripts placed as the manuscript sets them.
TABLE1_TEX <- c(
  "a^pop_mu"     = "$a^{\\mathrm{pop}}_{\\mu}$",
  "a^pop_med"    = "$a^{\\mathrm{pop}}_{\\mathrm{med}}$",
  "a^pop_sigma"  = "$a^{\\mathrm{pop}}_{\\sigma}$",
  "b^pop_mu"     = "$b^{\\mathrm{pop}}_{\\mu}$",
  "b^pop_med"    = "$b^{\\mathrm{pop}}_{\\mathrm{med}}$",
  "b^pop_sigma"  = "$b^{\\mathrm{pop}}_{\\sigma}$",
  "vlbar_peak"   = "$\\overline{vl}_{\\mathrm{peak}}$",
  "vl_sigma"     = "$vl_{\\sigma}$",
  "T^pop_mu"     = "$T^{\\mathrm{pop}}_{\\mu}$",
  "T^pop_med"    = "$T^{\\mathrm{pop}}_{\\mathrm{med}}$",
  "T^pop_sigma"  = "$T^{\\mathrm{pop}}_{\\sigma}$",
  "vl_thresh"    = "$vl_{\\mathrm{thresh}}$")

## Words as well as symbols, which is what R1.4a asks for ("describe parameters in tables, not
## just symbols"), and units, which the published table omitted entirely -- a growth rate of
## 3.11 means nothing without them.
##
## Each quantity gets an italic group header carrying the name and units, and its rows then
## name only the summary. That removes the triple-stacked a^pop_med, which was the most
## cramped label in the table, and puts the mean/median/SD triplet visibly together so the
## comparison between them can be read at a glance.
TABLE1_GROUP <- c(
  "a^pop_mu"    = "growth", "a^pop_med"    = "growth", "a^pop_sigma"  = "growth",
  "b^pop_mu"    = "decay",  "b^pop_med"    = "decay",  "b^pop_sigma"  = "decay",
  "vlbar_peak"  = "peak",   "vl_sigma"     = "peak",
  "T^pop_mu"    = "onset",  "T^pop_med"    = "onset",  "T^pop_sigma"  = "onset",
  "vl_thresh"   = "thresh")

TABLE1_GROUP_TEX <- c(
  growth = "Growth rate ($\\log_{10}$ copies mL$^{-1}$ day$^{-1}$)",
  decay  = "Decay rate ($\\log_{10}$ copies mL$^{-1}$ day$^{-1}$)",
  peak   = "Peak viral load ($\\log_{10}$ copies mL$^{-1}$)",
  onset  = "Delay from threshold crossing to symptom onset (days)",
  thresh = "Viral load at threshold crossing ($\\log_{10}$ copies mL$^{-1}$)",
  test   = "Probability of first test $k$ days after symptom onset")

## Within a group the row says only which summary it is, plus the symbol.
TABLE1_ROW_TEX <- c(
  "a^pop_mu"     = "population mean, $a^{\\mathrm{pop}}_{\\mu}$",
  "a^pop_med"    = "population median, $a^{\\mathrm{pop}}_{\\mathrm{med}}$",
  "a^pop_sigma"  = "population SD, $a^{\\mathrm{pop}}_{\\sigma}$",
  "b^pop_mu"     = "population mean, $b^{\\mathrm{pop}}_{\\mu}$",
  "b^pop_med"    = "population median, $b^{\\mathrm{pop}}_{\\mathrm{med}}$",
  "b^pop_sigma"  = "population SD, $b^{\\mathrm{pop}}_{\\sigma}$",
  "vlbar_peak"   = "population mean, $\\overline{vl}_{\\mathrm{peak}}$",
  "vl_sigma"     = "population SD, $vl_{\\sigma}$",
  "T^pop_mu"     = "population mean, $T^{\\mathrm{pop}}_{\\mu}$",
  "T^pop_med"    = "population median, $T^{\\mathrm{pop}}_{\\mathrm{med}}$",
  "T^pop_sigma"  = "population SD, $T^{\\mathrm{pop}}_{\\sigma}$",
  "vl_thresh"    = "$vl_{\\mathrm{thresh}}$")

table1_latex <- function(t1, strains, row_sep = "0.35ex") {
  lab <- c(wt = "Wild-type", alpha = "Alpha", delta = "Delta", omicron = "Omicron")

  ## Testing rows are not in the lookups: they are generated, and collapse from the old
  ## P(symp -> test = k) labels to "k = 0" under one group header.
  grp_of <- function(p) if (is.na(TABLE1_GROUP[p])) "test" else unname(TABLE1_GROUP[p])
  row_name <- function(p) {
    if (!is.na(TABLE1_ROW_TEX[p])) return(unname(TABLE1_ROW_TEX[p]))
    sprintf("$k = %s$", sub("^P\\(symp -> test = ([0-9])\\)$", "\\1", p))
  }
  ## Intervals are printed by fmt() as "m (lo, hi)"; LaTeX wants no space after the comma to
  ## match the published table, and the whole cell in text mode.
  cell <- function(x) gsub(", ", ",", x, fixed = TRUE)

  ## Group header rows span nothing -- they sit in the first column with the variant columns
  ## empty -- so the table stays a plain tabular and needs no extra package.
  body <- character(0); last <- ""
  for (i in seq_len(nrow(t1))) {
    p <- t1$parameter[i]
    g <- grp_of(p)
    if (g != last) {
      if (nchar(last)) body <- c(body, "\\addlinespace[0.6ex]")
      body <- c(body, paste0("\\textit{", unname(TABLE1_GROUP_TEX[g]), "}",
                             strrep(" &", length(strains)), " \\\\"))
      last <- g
    }
    body <- c(body,
              paste0("\\quad ", row_name(p), " & ",
                     paste(vapply(strains, function(s) cell(t1[[s]][i]), character(1)),
                           collapse = " & "), " \\\\"))
  }

  c("\\begin{table}[htbp]",
    "\\centering",
    "\\footnotesize",
    paste0("\\caption{\\textbf{Posterior medians and 95\\% credible intervals for model ",
           "parameters by SARS-CoV-2 variant.} For viral kinetics and symptom timing ",
           "parameters we report posterior summaries of the population parameters defined ",
           "in Equations 6 and 7, namely the population mean ($\\mu^{\\mathrm{pop}}$) and ",
           "population standard deviation ($\\sigma^{\\mathrm{pop}}$), together with the ",
           "population median ($\\mathrm{med}^{\\mathrm{pop}}$). Units are given in the ",
           "italicised heading for each quantity. Each entry is the posterior ",
           "median of that population quantity, with its 95\\% credible interval in ",
           "parentheses; the population summary describes variation between individuals, ",
           "while the interval describes posterior uncertainty about it. Growth and decay ",
           "rates are shifted lognormal, $a_i = 0.5 + \\mathrm{LogNormal}(\\eta_a, ",
           "\\omega_a)$ and $b_i = 0.25 + \\mathrm{LogNormal}(\\eta_b, \\omega_b)$, and the ",
           "delay to symptom onset is lognormal, so for these the population mean exceeds ",
           "the median by a factor that grows with between-individual variation. ",
           ## Quotes the reported (thresh_tdist) fit, so it is omitted for any other form.
           if (identical(FORM, "thresh_tdist"))
             paste0("The gap is ",
                    "largest for the wild-type growth rate, where $\\omega_a = 1.33$ places the mean ",
                    "at the 75th percentile of individuals: the mean separates the variants while the ",
                    "medians for wild type, Alpha and Delta are indistinguishable. ") else "",
           "We report both so ",
           "that the typical individual and the population average can each be read directly. ",
           "Peak viral load is normally distributed, so its mean and median coincide and only ",
           "one is given. The testing rows give the probability that an individual's first ",
           "test falls $k$ days after symptom onset, conditional on being tested within six ",
           "days.}"),
    "\\label{tab:1}",
    ## \addlinespace needs booktabs; \arraystretch alone would stretch the rules too. Left
    ## alignment on the parameter column, right on the estimates so the intervals line up.
    sprintf("\\setlength{\\extrarowheight}{%s}", row_sep),
    paste0("\\begin{tabular}{l", strrep("r", length(strains)), "}"),
    "\\hline",
    paste0("Quantity & ", paste(unname(lab[strains]), collapse = " & "), " \\\\"),
    "\\hline",
    body,
    "\\hline",
    "\\end{tabular}",
    "\\setlength{\\extrarowheight}{0pt}",
    "\\end{table}")
}

## ---------------------------------------------------------------------------
## figure_timescales -- the three time origins on one trajectory (R1.3)
##
##   figure_timescales()                        # wild type, posterior medians
##   figure_timescales("omicron")
##   figure_timescales(file = "outputs/figures/figS_timescales.png")
##
## The visual answer to "the treatment of time is unclear throughout". The manuscript uses
## three origins -- calendar day, detectability, and the crossing of vl_thresh -- and never
## draws them together, so a reader has to reconstruct the relationship from equations 3, 4
## and 8, which are pages apart and do not share a variable.
##
## Drawn from the POSTERIOR MEDIANS rather than sketched, so every interval on it is a
## quantity the paper reports and can be checked against Table 1 and Supplementary Table 6.
## A schematic with invented spacing would misplace the one thing the figure exists to show:
## that vl_thresh is crossed BEFORE detectability while symptom onset falls several days
## after it, and that peak viral load precedes onset.
## `vl_thresh` is drawn at an ILLUSTRATIVE value by default, not at its posterior median,
## and that is deliberate. The fitted parameter is pressed against its lower bound of zero --
## 59.5% of wild-type draws and 71.5% of delta draws lie within 0.01 of it -- so the data is
## not locating a threshold, it is pushing it as low as the prior permits. Drawing the
## schematic at 0.008 would assert a precision that is not there and would make the model
## look as though it requires a near-zero threshold, which it does not.
##
## Pass vl_thresh = "fitted" for the posterior median, or any number to place it by hand.
## Values ABOVE vl0 are handled: the threshold is then crossed after detectability and
## tau_det is negative, which is the case the notation has to survive.
figure_timescales <- function(strain = "wt", vl0 = 1.5, vl_thresh = 0.7,
                              dir = STRAIN_DIR,
                              thin = 50L, file = NULL, width = 8, height = 5) {

  d <- do.call(rbind, lapply(chain_tags(strain, dir),
                             function(t) combine_chains(strain, t, thin = thin, dir = dir)))
  md <- apply(d, 2, median)
  vt_fit <- median(d[, "inc3"]); vt_q <- quantile(d[, "inc3"], c(0.025, 0.975))
  rm(d); invisible(gc(verbose = FALSE))

  a  <- 0.5  + exp(md[["a_bar"]])
  b  <- 0.25 + exp(md[["b_bar"]])
  vp <- md[["vl_max_bar"]]
  vt <- if (identical(vl_thresh, "fitted")) vt_fit else as.numeric(vl_thresh)
  Ti <- exp(md[["inc1"]])                       # median of LogNormal(inc1, inc2)
  if (vt >= vp) stop("vl_thresh must lie below the peak (", round(vp, 2), ")")

  ## Everything measured from the threshold crossing, which is the origin equation 8 uses
  ## and the only one on which all four landmarks are positive.
  t_det  <- (vl0 - vt) / a
  t_peak <- (vp  - vt) / a
  t_end  <- t_peak + (vp - vl0) / b

  ## ---- compressed x axis --------------------------------------------------------
  ## Everything interesting happens in the first week: threshold, detectability, peak and
  ## onset fall within 6.2 days, while the decline runs to 22.5. Drawn linearly, the four
  ## landmarks are crushed into the left fifth of the panel and their labels collide.
  ##
  ## The decline after `brk` is therefore compressed by `squash`, and the axis is labelled
  ## with TRUE days throughout, with a break marker where the scale changes. The growth phase
  ## keeps full resolution, which is where the figure has to be read precisely.
  brk    <- min(Ti, t_peak) + 2
  squash <- 3.5
  tx  <- function(x) ifelse(x <= brk, x, brk + (x - brk) / squash)

  ## A fine sequence rather than three vertices, because the transform is piecewise and a
  ## straight segment drawn across the break would bend in the wrong place.
  xs <- sort(unique(c(seq(0, t_end, length.out = 400), t_peak, brk)))
  traj <- data.frame(
    x  = tx(xs),
    vl = ifelse(xs <= t_peak, vt + a * xs, vp - b * (xs - t_peak)))

  ## Landmarks, in trajectory order. Naming them here rather than inline keeps the label
  ## text and the geometry from drifting apart when the strain changes.
  ##
  ## Peak and onset are 1.2 days apart, so their labels are placed by hand with leader lines
  ## rather than stacked above the points: peak above, onset out to the right and below the
  ## declining limb, where the panel is empty.
  mk <- data.frame(
    x   = c(0, t_det, t_peak, Ti, t_end),
    vl  = c(vt, vl0, vp, NA, vl0),
    lab = c("crosses~italic(vl)[thresh]", "detectable~(italic(vl)[0])",
            "peak", "symptom~onset", "undetectable"),
    ## label offsets, in display units and viral-load units
    lx  = c(0, 0.4, 0, 2.6, 0),
    ly  = c(-0.9, 0.75, 0.75, -1.9, 0.75),
    hj  = c(0.5, 0, 0.5, 0, 1),
    stringsAsFactors = FALSE)
  ## Onset can fall after the peak, in which case its marker sits on the decline. Recompute
  ## rather than assuming the growth phase -- omicron's onset is earlier relative to peak.
  mk$vl[4] <- if (Ti <= t_peak) vt + a * Ti else vp - b * (Ti - t_peak)
  mk$x     <- tx(mk$x)
  mk$tx    <- mk$x + mk$lx
  mk$ty    <- mk$vl + mk$ly

  ## Two clock rows below the axis. y positions are in viral-load units because the panel
  ## shares one coordinate system; the axis simply does not label them.
  y1 <- -1.3; y2 <- -2.4
  arw <- grid::arrow(length = grid::unit(6, "pt"), type = "closed")

  ## True-day tick positions, thinned past the break so the compressed half is not crowded.
  ticks <- c(0, 2, 4, 6, seq(10, floor(t_end), by = 5))
  ticks <- ticks[ticks <= t_end]

  p <- ggplot() +
    ## the trajectory
    geom_line(data = traj, aes(x, vl), linewidth = 0.9, colour = STRAIN_COL[[strain]]) +
    geom_hline(yintercept = vl0, linetype = 3, colour = "grey35", linewidth = 0.4) +
    geom_hline(yintercept = vt,  linetype = 3, colour = "grey35", linewidth = 0.4) +
    ## drop lines from each landmark to the clock rows, so the timing reads vertically
    geom_segment(data = mk, aes(x = x, xend = x, y = vl, yend = y2),
                 linetype = 2, colour = "grey65", linewidth = 0.3) +
    ## leader lines from each label to its point, so labels can sit clear of the curve
    geom_segment(data = mk[mk$lx != 0 | mk$ly != 0, ],
                 aes(x = tx, xend = x, y = ty, yend = vl),
                 colour = "grey55", linewidth = 0.25) +
    geom_point(data = mk, aes(x, vl), size = 2.4, colour = "grey15") +
    geom_text(data = mk, aes(tx, ty, label = lab, hjust = hj), parse = TRUE,
              size = 3.1, colour = "grey15") +

    ## clock 1: tau, zero at detectability
    annotate("segment", x = tx(t_det), xend = tx(t_end) * 1.02, y = y1, yend = y1,
             colour = "grey20", linewidth = 0.5, arrow = arw) +
    annotate("point", x = tx(t_det), y = y1, size = 2, colour = "grey20") +
    annotate("text", x = tx(t_det), y = y1 + 0.45, label = "tau == 0", parse = TRUE,
             size = 3.1, hjust = 0.4) +
    annotate("text", x = tx(t_end) * 1.02, y = y1 + 0.45, hjust = 1, size = 3.1,
             label = "trajectory~time~tau", parse = TRUE) +

    ## clock 2: the delay T_i, threshold to onset
    annotate("segment", x = 0, xend = tx(Ti), y = y2, yend = y2,
             colour = "grey20", linewidth = 0.5, arrow = arw) +
    annotate("point", x = 0, y = y2, size = 2, colour = "grey20") +
    annotate("text", x = tx(Ti) / 2, y = y2 - 0.5, parse = TRUE, size = 3.1,
             label = sprintf("italic(T)[i] == %.1f~days", Ti)) +
    ## and the piece that reconciles the two. Sign-safe: tau_det runs threshold ->
    ## detectability, so it is drawn left to right when the threshold is below vl0 and right
    ## to left when it is above.
    annotate("segment", x = 0, xend = tx(t_det), y = y1, yend = y1,
             colour = "firebrick", linewidth = 0.7) +
    annotate("text", x = max(tx(t_det), 0), y = y1 - 0.5, hjust = 0, size = 3.1,
             colour = "firebrick", parse = TRUE,
             label = sprintf("tau[det] == %+.2f~days", t_det)) +

    ## the break in the x scale, marked so the compression cannot be read as real time
    annotate("text", x = tx(brk), y = vp + 0.5, label = "//", size = 3.4,
             colour = "grey45", angle = 20) +
    annotate("segment", x = tx(brk), xend = tx(brk), y = y2 - 0.7, yend = vp + 0.3,
             linetype = 3, colour = "grey70", linewidth = 0.3) +

    scale_y_continuous(breaks = seq(0, 8, 2)) +
    scale_x_continuous(breaks = tx(ticks), labels = ticks) +
    coord_cartesian(ylim = c(y2 - 1.1, vp + 1.3),
                    xlim = c(min(-1, tx(t_det) - 1), tx(t_end) * 1.10)) +
    labs(x = expression("days since viral load reaches "*italic(vl)[thresh]),
         y = expression(log[10]~RNA~copies~mL^-1),
         title = sprintf("Time scales in the model (%s)", STRAIN_LABEL[[strain]]),
         subtitle = sprintf(paste("kinetics at posterior medians; the threshold is drawn at",
                                  "%.2f for legibility (fitted median %.3f, 95%% CrI",
                                  "%.3f-%.3f, against a lower bound of 0).\nThe time axis is",
                                  "compressed %.1f-fold beyond day %.0f; tick labels are true",
                                  "days throughout."),
                            vt, vt_fit, vt_q[1], vt_q[2], squash, brk)) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text = element_text(colour = "black"))

  if (!is.null(file)) ggsave(file, p, width = width, height = height, dpi = 300)
  print(p)
  invisible(list(plot = p, landmarks = mk,
                 timings = c(tau_det = t_det, T_i = Ti, peak = t_peak, end = t_end)))
}

## ---------------------------------------------------------------------------
## figure_pairs -- posterior correlations between the fitted parameters (R2.7)
##
##   figure_pairs("alpha")                      # scatter matrix for one variant
##   figure_pairs(file = "outputs/figures/figS_pairs.png")
##   figure_corr()                              # correlation summary, all four variants
##
## R2.7 notes that Figure 4 shows marginals only and asks for the posterior covariance. The
## marginals genuinely do hide structure, and the specific thing to show is that SYMPTOM
## ONSET TIMING AND GROWTH RATE ARE JOINTLY CONSTRAINED: eta_T and a_bar correlate at -0.61
## to -0.92 across the four variants.
##
## That has a mechanism. The fitted vl_thresh sits at its lower bound of zero, so the time
## from threshold crossing to detectability is (vl_0 - vl_thresh)/a_i ~ 1.5/a_i, and the
## growth rate alone then determines where symptom onset falls along the trajectory. Data
## that pins the SHAPE of the trajectory more tightly than its POSITION in time will trade
## these two off, which is what the posterior shows.
##
## Note this is NOT the inc1-inc3 ridge seen in the simulation scenarios, where vl_thresh was
## free to move and absorbed the trade-off instead. In the variant fits inc1-inc3 is only
## -0.05 to -0.39. Same degeneracy, different pair, because the threshold is at its boundary.
## omega_T and vl_thresh are included because their joint behaviour is a finding, not a
## detail: viral load at symptom onset is vl_thresh + a*T while onset precedes the peak, so
## the two trade off along a ridge at fixed onset load. For omicron that shows as
## cor(vl_thresh, eta_T) near -0.7 while vl_thresh alone varies over a coefficient of
## variation of 68% and the combination over 2%. R2.7 asks for exactly this structure.
PAIR_PARS <- c("a_bar", "a_sigma", "b_bar", "b_sigma", "vl_max_bar", "vl_max_sigma",
               "inc1", "inc2", "inc3")
## Labels come from REC_SYM, the one symbol table, read at CALL time rather than assigned here.
## REC_SYM is defined further down the file and follows NOTATION_OMEGA; a copy taken at this
## point would be fixed at the omega spelling and would silently disagree with every other
## figure if that flag were ever set FALSE. Everything in PAIR_PARS has a REC_SYM entry, which
## pair_labels() asserts rather than leaving to produce NA strips.
pair_labels <- function(p) {
  miss <- setdiff(p, names(REC_SYM))
  if (length(miss))
    stop("no symbol in REC_SYM for: ", paste(miss, collapse = ", "), call. = FALSE)
  unname(REC_SYM[p])
}

pair_draws <- function(strain, pars = PAIR_PARS, thin = 20L,
                       dir = STRAIN_DIR) {
  d <- do.call(rbind, lapply(chain_tags(strain, dir),
                             function(t) combine_chains(strain, t, thin = thin, dir = dir)))
  d[, intersect(pars, colnames(d)), drop = FALSE]
}

figure_pairs <- function(strain = "alpha", pars = PAIR_PARS, thin = 20L, n_show = 2000L,
                         dir = STRAIN_DIR, file = NULL,
                         width = 8.5, height = 8) {
  m <- pair_draws(strain, pars, thin, dir)
  ## Thinned again for the scatter only: 2,000 points show the shape, and 40,000 overplot
  ## into a solid mass that hides exactly the structure the figure exists to show. The
  ## correlations printed in the upper triangle use ALL the draws.
  cm <- cor(m)
  idx <- if (nrow(m) > n_show) sort(sample.int(nrow(m), n_show)) else seq_len(nrow(m))

  ## Custom upper panel, passed as a bare function rather than through GGally::wrap(), which
  ## expects a named panel type. It reads the correlation out of `cm` -- computed above from
  ## ALL the draws -- so the printed number is not the correlation of the 2,000 plotted.
  upper_fn <- function(data, mapping, ...) {
    r <- cm[rlang::as_name(mapping$x), rlang::as_name(mapping$y)]
    strong <- abs(r) > 0.5
    ggplot() +
      annotate("text", 0, 0, label = sprintf("%+.2f", r), size = if (strong) 4 else 3.3,
               colour = if (strong) "firebrick" else "grey30",
               fontface = if (strong) "bold" else "plain") +
      theme_void()
  }

  lab <- pair_labels(colnames(m))
  p <- GGally::ggpairs(
    as.data.frame(m[idx, , drop = FALSE]),
    columnLabels = lab, labeller = "label_parsed",
    upper = list(continuous = upper_fn),
    lower = list(continuous = GGally::wrap("points", alpha = 0.15, size = 0.35,
                                           colour = STRAIN_COL[[strain]])),
    diag  = list(continuous = GGally::wrap("densityDiag", colour = STRAIN_COL[[strain]],
                                           fill = STRAIN_COL[[strain]], alpha = 0.3))) +
    labs(title = sprintf("Posterior correlations, %s", STRAIN_LABEL[[strain]]),
         subtitle = paste("correlations computed from all draws; scatter thinned to",
                          format(length(idx), big.mark = ","), "points")) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          axis.text = element_text(size = 6, colour = "black"))

  if (!is.null(file)) ggsave(file, p, width = width, height = height, dpi = 300)
  print(p)
  invisible(list(plot = p, cor = cm))
}

## The same structure across all four variants, as a correlation matrix rather than scatters
## -- four full pairs matrices would be unreadable, and the point is that the pattern repeats.
figure_corr <- function(strains = STRAINS, pars = PAIR_PARS, thin = 20L,
                        dir = STRAIN_DIR, file = NULL,
                        width = 9, height = 3) {
  d <- do.call(rbind, lapply(strains, function(s) {
    cm <- cor(pair_draws(s, pars, thin, dir))
    cm[upper.tri(cm, diag = TRUE)] <- NA
    z <- as.data.frame(as.table(cm))
    names(z) <- c("v1", "v2", "r")
    z$strain <- STRAIN_LABEL[[s]]
    z[!is.na(z$r), ]
  }))
  ## Rows and columns need SEPARATE level sets, and taking one for both is how eta_a went
  ## missing: below the diagonal the first parameter never appears as a row and the last never
  ## as a column, so a level set built from the rows alone excludes the first parameter -- and
  ## every correlation in which it appears as a COLUMN then falls outside the factor levels and
  ## is silently dropped. That removed eta_a's entire column, including its correlation with
  ## eta_T, which is one of the two the caption singles out.
  ##
  ## Dropping the empty first row and last column as well keeps the triangle compact; the
  ## alternative, one level set containing everything, renders a blank row and a blank column.
  allp <- intersect(pars, unique(c(as.character(d$v1), as.character(d$v2))))
  rows <- allp[-1]
  cols <- allp[-length(allp)]
  d$v1 <- factor(d$v1, levels = rows, labels = pair_labels(rows))
  d$v2 <- factor(d$v2, levels = rev(cols), labels = rev(pair_labels(cols)))
  if (anyNA(d$v1) || anyNA(d$v2))
    stop("correlation cells fell outside the axis levels -- ",
         sum(is.na(d$v1) | is.na(d$v2)), " of ", nrow(d), " would be dropped")
  d$strain <- factor(d$strain, levels = unname(STRAIN_LABEL[strains]))

  p <- ggplot(d, aes(v1, v2, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%+.2f", r)), size = 2.4,
              colour = ifelse(abs(d$r) > 0.6, "white", "grey20")) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_discrete(labels = scales::parse_format()) +
    facet_wrap(~strain, nrow = 1) +
    labs(x = NULL, y = NULL, title = "Posterior correlations by variant") +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold"),
          axis.text = element_text(colour = "black"),
          axis.ticks = element_blank())

  if (!is.null(file)) ggsave(file, p, width = width, height = height, dpi = 300)
  print(p)
  invisible(list(plot = p, table = d))
}

## ---------------------------------------------------------------------------
## save_paper_outputs -- regenerate every figure and table the manuscript uses
##
##   save_paper_outputs()                       # everything, into outputs/figures and tables
##   save_paper_outputs(what = "tables")        # tables only, seconds rather than minutes
##   save_paper_outputs(figures = c("4", "5"))  # just those two
##
## Figures land in outputs/figures and LaTeX in outputs/tables, named by their number in the
## manuscript rather than by what they plot, so a file maps to a place in the paper without
## anyone having to remember which is which.
##
## The names are the settled allocation, which is worth stating because it changed once:
##   Fig 7  precision against sample size          -- MAIN text, the R1.4c answer
##   Fig S  calibration, z against sample size     -- supplementary, NOT a recovery figure:
##          z is scale-free, so it has no expected trend with n and cannot show the benefit
##          of more data. It shows the posteriors are honest, which is a different claim.
##   Fig S  contraction against z                  -- supplementary
##
## Figures 1, 2 and 6 are NOT regenerated. Figure 1 is a hand-built schematic, Figure 6 is a
## standalone phase-bias illustration, and Figure 2's code was never saved -- it is analysis
## of the data rather than of the fits, so it needs rebuilding for the J4 export rather than
## for correctness.
## Directories come from use_form(). Re-sourcing keeps whichever form is already selected.
use_form(if (exists("FORM", envir = globalenv())) get("FORM", envir = globalenv()) else "thresh_tdist")

## `form` switches STRAIN_DIR, FIG_DIR and TAB_DIR for the duration of the call and restores
## the previous form afterwards. Under "incidence" the default figure set drops the
## simulation-study figures, which were fitted with thresh_tdist only, and adds the phase check
## (R2.5) and the convergence table, which have no reported equivalent for the hybrid yet.
save_paper_outputs <- function(what = c("both", "figures", "tables"),
                               figures = NULL, form = FORM,
                               seed = 20260907L) {
  what <- match.arg(what)
  old <- FORM; use_form(form); on.exit(use_form(old), add = TRUE)
  if (is.null(figures))
    figures <- if (form == "incidence")
      c("3", "4", "5", "S_prior_posterior", "S_correlations", "S_pairs", "S_phase_ppc")
    else c("3", "4", "5", "7", "S_calibration", "S_contraction", "S_prior_posterior",
           "S_correlations", "S_pairs")
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
  fp <- function(n) file.path(FIG_DIR, n)

  if (what %in% c("both", "tables")) {
    writeLines(table1()$latex,  file.path(TAB_DIR, "table1.tex"))
    writeLines(table2_latex(),  file.path(TAB_DIR, "supp_table4.tex"))
    writeLines(table3_latex(),  file.path(TAB_DIR, "supp_table5.tex"))
    ## table4() simulates populations per draw, so it is the one table that needs a seed to
    ## reproduce. Set here rather than inside it, so a caller who wants a different draw can
    ## still get one.
    set.seed(seed)
    writeLines(table4_latex(),  file.path(TAB_DIR, "supp_table6.tex"))
    if (form == "incidence")
      writeLines(capture.output(table_ess(latex = TRUE, show_rhat = TRUE)),
                 file.path(TAB_DIR, "supp_table_ess.tex"))
    message("tables written to ", TAB_DIR)
  }

  if (what %in% c("both", "figures")) {
    set.seed(seed)
    if ("3" %in% figures) figure3(file = fp("fig3_model_fit.png"))
    if ("4" %in% figures) figure4(file = fp("fig4_kinetics.png"))
    if ("5" %in% figures) figure5(file = fp("fig5_incidence.png"))
    if ("7" %in% figures) figure_precision(file = fp("fig7_precision.png"))
    if ("S_calibration"    %in% figures)
      figure_recovery_summary(file = fp("figS_calibration.png"))
    if ("S_contraction"    %in% figures)
      figure_contraction(file = fp("figS_contraction.png"))
    if ("S_prior_posterior" %in% figures)
      figure_prior_posterior(file = fp("figS_prior_posterior.png"))
    ## R2.7: one pairs matrix per variant, plus the cross-variant correlation summary. Four
    ## matrices because the reviewer asked for the posterior covariance and the structure is
    ## worth seeing per variant; the summary because the point is that it repeats.
    if ("S_correlations" %in% figures)
      figure_corr(file = fp("figS_correlations.png"))
    if ("S_pairs" %in% figures)
      for (s in STRAINS) figure_pairs(s, file = fp(sprintf("figS_pairs_%s.png", s)))
    if ("S_phase_ppc" %in% figures)
      figure_phase_ppc(file = fp("figS_phase_ppc.png"))
    message("figures written to ", FIG_DIR)
  }
  invisible(list(figures = FIG_DIR, tables = TAB_DIR))
}

## ---------------------------------------------------------------------------
## LaTeX for the three supplementary tables
##
##   writeLines(table2_latex(), "outputs/supp_table4.tex")   # prior mean + covariance
##   writeLines(table3_latex(), "outputs/supp_table5.tex")   # posteriors + prior column
##   writeLines(table4_latex(), "outputs/supp_table6.tex")   # population median + IQR
##
## NOTE THE NUMBERING. The manuscript's supplementary tables continue the main sequence, so
## table2() feeds Supplementary Table 4, table3() feeds Table 5 and table4() feeds Table 6.
## The generator names follow the order they were written, not the manuscript's. Getting
## this wrong pastes a covariance matrix over a posterior summary.
##
## All three were stale at resubmission: Table 4 carried the pre-switch longitudinal prior
## (eta_a -1.246111 against the current -1.1197317, see U7), and Tables 5 and 6 carried the
## wild-type spurious mode (phi_1 2.54, a 12.7-day onset; total detectable 30.37 days).
## ---------------------------------------------------------------------------

## The lognormal SCALE parameters are emitted as omega_a/omega_b rather than tau_a/tau_b,
## which frees tau to mean trajectory time (NOTATION.md section 2.1). Adopted 2026-09-07.
##
## sigma_a/sigma_b was the obvious alternative and was rejected: Table 1 reports
## a_sigma^pop and b_sigma^pop, so sigma_a (spread on the log scale) would sit beside
## a_sigma^pop (SD on the natural scale), two different quantities both reading as "sigma of
## a" -- and equation 7 converts between them, which makes the confusion likelier.
##
## Set FALSE to emit the manuscript's original tau_a/tau_b.
##
## Note this controls the LATEX only. The data frames returned by table2() and table3() keep
## tau_a/tau_b, because those are the names stored in long_post_6.rds and long_cov_6.rds and
## the console view should match the files it reads. tex_sym() is the single place the
## display name is decided.
NOTATION_OMEGA <- TRUE

## The manuscript is internally inconsistent here: Supplementary Table 4 writes the peak as
## \overline{vl}_peak and its SD as vl_sigma, while Table 5 writes \overline{v}_peak and
## sigma_vl for the same two parameters. Standardised on Table 4's form, which also matches
## Table 1 and equation 5.
tex_sym <- function() {
  ## Both renames move together. eta_T/omega_T exists precisely so that the onset delay uses
  ## the same location-scale convention as the growth and decay rates, so emitting eta_T
  ## alongside tau_a would defeat the purpose -- a reader would meet two conventions and no
  ## reason for either.
  sc <- if (NOTATION_OMEGA) c("\\omega_{a}", "\\omega_{b}") else c("\\tau_{a}", "\\tau_{b}")
  on <- if (NOTATION_OMEGA) c("\\eta_{T}", "\\omega_{T}") else c("\\phi_{1}", "\\phi_{2}")
  c(eta_a       = "$\\eta_{a}$",
    tau_a       = paste0("$", sc[1], "$"),
    eta_b       = "$\\eta_{b}$",
    tau_b       = paste0("$", sc[2], "$"),
    vl_peak_bar = "$\\overline{vl}_{\\mathrm{peak}}$",
    sigma_vl    = "$vl_{\\sigma}$",
    phi_1       = paste0("$", on[1], "$"),
    phi_2       = paste0("$", on[2], "$"),
    v_thresh    = "$vl_{\\mathrm{thresh}}$",
    xi          = "$\\xi$",
    delta_0     = "$\\delta_{0}$", delta_1 = "$\\delta_{1}$", delta_2 = "$\\delta_{2}$",
    delta_3     = "$\\delta_{3}$", delta_4 = "$\\delta_{4}$", delta_5 = "$\\delta_{5}$",
    M           = "$M$",
    mean        = "\\textbf{mean}")
}

## Intervals are printed as "m (lo, hi)" by table1/table3 and "m (lo-hi)" by table4. Both are
## emitted as "m (lo,hi)".
##
## An en-dash separator was tried first and abandoned: several Table 6 quantities are
## negative -- "Symptom onset -> peak" is -0.89 for wild type -- so the raw string is
## "-3.91--3.44", and turning the separator into an en-dash yields "---", an em-dash followed
## by a minus. A comma is unambiguous whatever the sign, and matches Tables 1 and 5.
tex_cell <- function(x) {
  x <- ifelse(is.na(x), "---", x)
  x <- gsub(", ", ",", x, fixed = TRUE)
  ## the separator only: a hyphen preceded by a digit. A leading minus is preceded by "(".
  gsub("(?<=[0-9])-", ",", x, perl = TRUE)
}

tex_wrap <- function(caption, label, colspec, header, body) {
  c("\\begin{table}[htbp]", "\\centering", "\\footnotesize",
    paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"),
    paste0("\\begin{tabular}{", colspec, "}"), "\\hline",
    header, "\\hline", body, "\\hline", "\\end{tabular}", "\\end{table}")
}

table2_latex <- function(t2 = table2(), digits = 7) {
  sym <- tex_sym()
  cols <- setdiff(names(t2), "row")
  num  <- function(v) ifelse(is.na(v), "---", formatC(v, format = "f", digits = digits))
  body <- vapply(seq_len(nrow(t2)), function(i) {
    nm <- t2$row[i]
    paste0(if (!is.na(sym[nm])) unname(sym[nm]) else nm, " & ",
           paste(num(unlist(t2[i, cols])), collapse = " & "), " \\\\")
  }, character(1))

  tex_wrap(
    paste0("\\textbf{Prior mean vector and covariance matrix for the viral kinetics ",
           "hyperparameters used in the main model.} These parameters were estimated by ",
           "fitting the within-host viral load model to longitudinal viral load ",
           "measurements. The resulting posterior mean vector and covariance matrix were ",
           "used to define a multivariate normal prior for the viral kinetics parameters ",
           "$\\theta_v$ in the main model. Given these parameters describe a lognormal ",
           "distribution governing the individual viral growth and decay rates, negative ",
           "values of $\\eta_a$ and $\\eta_b$ correspond to positive growth and decay rates ",
           "on the natural scale."),
    "tab:S4", paste0("l", strrep("r", length(cols))),
    paste0(" & ", paste(vapply(cols, function(c) unname(sym[c]), character(1)),
                        collapse = " & "), " \\\\"),
    body)
}

table3_latex <- function(t3 = table3()) {
  sym  <- tex_sym()
  ## The caption names the scale parameters in prose, so it has to follow NOTATION_OMEGA too
  ## -- hardcoding omega here would leave the caption contradicting the table body the moment
  ## anyone set the flag back.
  sc   <- if (NOTATION_OMEGA) c("\\omega_a", "\\omega_b") else c("\\tau_a", "\\tau_b")
  on   <- if (NOTATION_OMEGA) c("\\eta_T", "\\omega_T") else c("\\phi_1", "\\phi_2")
  cols <- setdiff(names(t3), "parameter")
  lab  <- c(prior = "Prior", wt = "Wild-type", alpha = "Alpha",
            delta = "Delta", omicron = "Omicron")
  body <- vapply(seq_len(nrow(t3)), function(i) {
    nm <- t3$parameter[i]
    paste0(if (!is.na(sym[nm])) unname(sym[nm]) else nm, " & ",
           paste(tex_cell(unlist(t3[i, cols])), collapse = " & "), " \\\\")
  }, character(1))

  tex_wrap(
    paste0("\\textbf{Posterior estimates for model parameters across SARS-CoV-2 variants.} ",
           "Values shown are posterior medians with 95\\% credible intervals in parentheses. ",
           "The prior column reports the corresponding prior distributions (median and 95\\% ",
           "interval) used in the model, sampled from the truncated distributions the ",
           "sampler evaluates. Parameters $\\eta_a$ and $", sc[1], "$ describe the location ",
           "and scale of the log-normal distribution governing the viral growth rate across ",
           "individuals, while $\\eta_b$ and $", sc[2], "$ describe the corresponding ",
           "distribution for viral decay rates. $\\overline{vl}_{\\mathrm{peak}}$ and ",
           "$vl_{\\sigma}$ denote the mean and standard deviation of the population ",
           "distribution of peak viral load. Parameters $", on[1], "$ and $", on[2], "$ ",
           "describe the location and scale of the lognormal distribution of the delay from ",
           "the time viral load first reaches the threshold $vl_{\\mathrm{thresh}}$ to ",
           "symptom onset. The parameter $\\xi$ denotes ",
           "the dispersion parameter of the observation model. Testing behaviour following ",
           "symptom onset is described by $\\delta_0, \\ldots, \\delta_5$, the probability of ",
           "being tested by day $k$ after symptom onset conditional on not having been ",
           "tested previously. $M$ denotes the spline multiplier governing overall incidence ",
           "scaling; no prior interval is reported for it because the sampler assigns it an ",
           "improper flat prior, the incidence gauge being fixed instead by pinning one ",
           "spline coefficient."),
    "tab:S5", paste0("l", strrep("l", length(cols))),
    paste0("Parameter & ", paste(unname(lab[cols]), collapse = " & "), " \\\\"),
    body)
}

table4_latex <- function(t4 = table4()) {
  cols <- c("lower_25", "median_50", "upper_75")
  ## Strain blocks, as the published table sets them: a bold spanning row then its rows.
  body <- character(0)
  for (s in unique(t4$strain)) {
    z <- t4[t4$strain == s, , drop = FALSE]
    body <- c(body, sprintf("\\multicolumn{4}{l}{\\textbf{%s}}\\\\", s),
              vapply(seq_len(nrow(z)), function(i)
                paste0("\\quad ", gsub("->", "$\\rightarrow$", z$parameter[i], fixed = TRUE),
                       " & ", paste(tex_cell(unlist(z[i, cols])), collapse = " & "), " \\\\"),
                character(1)))
  }

  tex_wrap(
    paste0("\\textbf{Posterior summaries of viral kinetics parameters by SARS-CoV-2 ",
           "variant.} For each parameter we report the population median and the 25th and ",
           "75th percentiles of the population distribution of individual-level values. ",
           "Entries show the posterior median estimate of each quantity, with 95\\% credible ",
           "intervals in parentheses. Note that these differ from Table~\\ref{tab:1}, which ",
           "reports population means and standard deviations; here we present population ",
           "medians and interquartile ranges."),
    "tab:S6", "llll",
    paste0("Strain / parameter & Lower 25$^{\\mathrm{th}}$ percentile & ",
           "50$^{\\mathrm{th}}$ percentile & Upper 25$^{\\mathrm{th}}$ percentile \\\\"),
    body)
}

## ---------------------------------------------------------------------------
## table2 -- the longitudinal prior: mean vector and covariance matrix
##
##   t2 <- table2(); print(t2, row.names = FALSE)
##   write.csv(table2(), "outputs/table2.csv", row.names = FALSE)
##
## Reads the two files model.R actually loads, so what this prints is the
## prior the sampler used. That matters here: the values in the SUBMITTED Table 2 do not come
## from these files. The manuscript reports eta_a = -1.246111 with variance 0.0573685932,
## where long_post_6.rds/long_cov_6.rds give -1.1197317 and 0.0478283162 -- the mean differs
## by 0.13 and the variance by 20%. The submitted table came from an earlier longitudinal fit
## (model.R still carries the commented-out paths it was switched away
## from), so Table 2 needs regenerating rather than copying forward.
##
## l1 is excluded. It is model parameter 7, one of the two in ignored_pars, never proposed and
## fixed at zero throughout -- which is why the published table shows six parameters and the
## stored objects carry seven.
## ---------------------------------------------------------------------------
T2_SYMBOLS <- c(a_bar = "eta_a", a_sigma = "tau_a", b_bar = "eta_b", b_sigma = "tau_b",
                vl_max_bar = "vl_peak_bar", vl_max_sigma = "sigma_vl")

table2 <- function(post_file = "longitudinal_fits/long_post_6.rds",
                   cov_file  = "longitudinal_fits/long_cov_6.rds",
                   pars = names(T2_SYMBOLS), mean_digits = 6, cov_digits = 10) {

  post <- readRDS(post_file)
  cv   <- as.matrix(readRDS(cov_file))

  miss <- setdiff(pars, intersect(names(post$mean), colnames(cv)))
  if (length(miss)) stop("not in both files: ", paste(miss, collapse = ", "))
  mu <- post$mean[pars]
  cv <- cv[pars, pars, drop = FALSE]

  ## The two files are separate records of one fit, so they must agree. post$sd should be the
  ## square root of the covariance diagonal; if it is not, one of them has been regenerated
  ## without the other and the prior the sampler builds is internally inconsistent.
  if (!is.null(post$sd)) {
    d <- abs(sqrt(diag(cv)) - post$sd[pars])
    if (any(d > 1e-8))
      warning("post$sd disagrees with sqrt(diag(cov)) for: ",
              paste(pars[d > 1e-8], collapse = ", "), call. = FALSE)
  }
  if (!isSymmetric(unname(cv))) warning("covariance is not symmetric", call. = FALSE)

  sym <- unname(T2_SYMBOLS[pars])
  out <- rbind(
    data.frame(row = "mean",
               setNames(as.list(round(as.numeric(mu), mean_digits)), sym),
               stringsAsFactors = FALSE, check.names = FALSE),
    data.frame(row = sym,
               setNames(as.data.frame(round(cv, cov_digits)), sym),
               stringsAsFactors = FALSE, check.names = FALSE))
  rownames(out) <- NULL
  out
}

## ---------------------------------------------------------------------------
## table3 -- fitted parameters on their own scale, against their priors
##
##   t3 <- table3(); print(t3, row.names = FALSE)
##
## The supplementary table of parameters AS FITTED -- the logmeans and logsds of the shifted
## lognormals rather than the population quantities Table 1 derives from them. This is where
## a reader checks whether an estimate is data-driven, so the prior column is the point of it.
##
## The multiplier M has no row-generating prior: prior() gives the multivariate block to
## parameters 1:8 and independent truncated normals to 9:(length-1), which stops one short of
## it. Its prior cell is therefore NA rather than invented. The submitted table prints
## 7.34 (0.46-23.16) there, which does not correspond to the N(1, 10) in the Supplementary
## and is worth checking before resubmission.
## ---------------------------------------------------------------------------
T3_ROWS <- c(a_bar = "eta_a", a_sigma = "tau_a", b_bar = "eta_b", b_sigma = "tau_b",
             vl_max_bar = "vl_peak_bar", vl_max_sigma = "sigma_vl",
             inc1 = "phi_1", inc2 = "phi_2", inc3 = "v_thresh", rdisp = "xi",
             test0 = "delta_0", test1 = "delta_1", test2 = "delta_2",
             test3 = "delta_3", test4 = "delta_4", test5 = "delta_5",
             multiplier = "M")

## prior_bounds truncates a prior for REPORTING only, where the stored bound is looser than the
## model's real support. rdisp is the case: lower_bound is -Inf, so its declared prior is
## N(1, 10) with 46% of its mass negative, and the cell came out as "0.99 (-18.54, 20.37)" -- a
## 95% interval for a dispersion parameter that starts below -18. The likelihood forbids that
## region, so truncating changes only a normalising constant that cancels in every acceptance
## ratio; see prior_draws(). Matches figure_prior_posterior, so the table and the figure cannot
## describe the same prior differently.
table3 <- function(strains = c("wt", "alpha", "delta", "omicron"),
                   dir = STRAIN_DIR, thin = 10L, n_prior = 50000L, digits = 2,
                   prior_bounds = c(rdisp = 0)) {

  q <- c(0.5, 0.025, 0.975)
  fmt <- function(v) sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
                             v[1], v[2], v[3])

  ## Prior column. prior_draws() samples exactly what prior() evaluates -- the TRUNCATED
  ## multivariate normal for the kinetics block, independent truncated normals for the rest
  ## -- so these are the priors the sampler actually used, not a re-derivation of them.
  pri_pars <- setdiff(names(T3_ROWS), "multiplier")
  pd <- prior_draws(strains[1], n_prior, dir, pars = pri_pars,
                    bounds = prior_bounds[intersect(names(prior_bounds), pri_pars)])

  ## The multiplier's cell is deliberately left empty, and the reason belongs in the caption.
  ##
  ## input_generator DECLARES a prior for it -- multiplier mean 1, SD 10, which enters
  ## prior_mean and prior_sds -- but prior() never applies it: it evaluates the truncated
  ## normals over parameters[9:last_par] with last_par = length(parameters) - 1, and the
  ## multiplier is the last parameter. So the sampler runs it under an improper flat prior.
  ##
  ## The submitted Table 5 reports 7.34 (0.46-23.16) here, which is the median and 95%
  ## interval of the DECLARED prior, N(1, 10) truncated at 0. That is a prior the likelihood
  ## never saw. Reporting it again would repeat the error, so the cell stays NA and the
  ## caption says the parameter is assigned a flat prior. The gauge degeneracy is broken by
  ## pinning a spline coefficient rather than by this prior, so nothing depends on it.
  mult_prior <- NA_character_

  keep <- NULL; cols <- NULL
  for (s in strains) {
    d <- do.call(rbind, lapply(chain_tags(s, dir),
                               function(t) combine_chains(s, t, thin = thin, dir = dir)))
    if (is.null(keep)) keep <- intersect(names(T3_ROWS), colnames(d))
    cols <- cbind(cols, vapply(keep, function(p)
                    fmt(unname(quantile(d[, p], q))), character(1)))
    rm(d); invisible(gc(verbose = FALSE))
  }

  prior_col <- vapply(keep, function(p)
    if (p %in% colnames(pd)) fmt(unname(quantile(pd[, p], q)))
    else if (p == "multiplier") mult_prior
    else NA_character_,
    character(1))

  out <- data.frame(parameter = unname(T3_ROWS[keep]), prior = prior_col, cols,
                    stringsAsFactors = FALSE)
  colnames(out) <- c("parameter", "prior", strains)
  rownames(out) <- NULL
  out
}

## ---------------------------------------------------------------------------
## table4 -- the supplementary table behind Figure 4
##
##   t4 <- table4(); print(t4)
##   write.csv(table4(vl0 = 2), "outputs/table4_100copies.csv", row.names = FALSE)
##
## Population MEDIAN and interquartile range of the individual-level distribution, each with
## its 95% credible interval -- the manuscript's caption is explicit that this differs from
## Table 1, which reports means and SDs.
##
## Growth, decay and peak come from closed forms (fig4_closed); the seven time intervals are
## SIMULATED (fig4_draw), because they are nonlinear functions of correlated random variables
## with no closed-form quantiles. That distinction matters numerically: plugging the median
## growth, decay and peak into the duration formula gives 31.2 d for wild type, while the
## median of the individual-level duration distribution is 30.4 d.
## ---------------------------------------------------------------------------
T4_ROWS <- c(growth         = "Growth rate",
             decay          = "Decay rate",
             thresh_to_symp = "Time threshold -> symptoms",
             peak_to_end    = "Peak -> undetectable",
             symp_to_end    = "Symptoms -> undetectable",
             symp_to_peak   = "Symptom onset -> peak",
             total          = "Total detectable time",
             det_to_peak    = "Detectable -> peak",
             det_to_symp    = "Detectable -> symptoms",
             peak           = "Peak viral load")

table4 <- function(strains = STRAINS, n_samples = 300, n_indiv = 20000, vl0 = 1.5,
                   dir = STRAIN_DIR, digits = 2) {
  fmt <- function(v) sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
                             v[2], v[1], v[3])

  do.call(rbind, lapply(strains, function(s) {
    message(s)
    dr <- sample_draws(s, n_samples, dir)$draws

    ## Per posterior draw: the population 25th, 50th and 75th percentiles of every quantity.
    per_draw <- lapply(seq_len(nrow(dr)), function(j) {
      p <- dr[j, ]
      d <- fig4_draw(p, n_indiv, vl0)
      if (is.null(d)) return(NULL)
      sim <- vapply(d, quantile, numeric(3), probs = c(0.25, 0.5, 0.75))
      rownames(sim) <- c("q25", "mid", "q75")
      cl <- fig4_closed(p)[c("q25", "mid", "q75"), , drop = FALSE]
      cbind(cl, sim)
    })

    drop <- vapply(per_draw, is.null, logical(1))
    if (all(drop)) stop("no usable draws for ", s)
    if (any(drop))
      message(sprintf("  %s: %d of %d draws dropped -- fewer than 100 detectable individuals",
                      s, sum(drop), length(per_draw)))
    per_draw <- per_draw[!drop]

    rows <- intersect(names(T4_ROWS), colnames(per_draw[[1]]))
    ## Across draws: the posterior median and 95% CrI of each population percentile.
    out <- do.call(rbind, lapply(rows, function(k) {
      grab <- function(r) vapply(per_draw, function(m) m[r, k], numeric(1))
      data.frame(strain    = STRAIN_LABEL[[s]],
                 parameter = unname(T4_ROWS[k]),
                 lower_25  = fmt(quantile(grab("q25"), c(0.025, 0.5, 0.975))),
                 median_50 = fmt(quantile(grab("mid"), c(0.025, 0.5, 0.975))),
                 upper_75  = fmt(quantile(grab("q75"), c(0.025, 0.5, 0.975))),
                 stringsAsFactors = FALSE)
    }))
    rm(dr, per_draw); invisible(gc(verbose = FALSE))
    out
  }))
}

## ============================================================================
## FIGURES
##
## These need the model itself -- p_array_func simulates the expected counts, which runs
## the C++ kernel. It is loaded lazily inside ppc_draws so that check_chains() and table1()
## work without compiling anything.
## ============================================================================

STRAINS      <- c("wt", "alpha", "delta", "omicron")
STRAIN_LABEL <- c(wt = "Wild type", alpha = "Alpha", delta = "Delta", omicron = "Omicron")
STRAIN_COL   <- c(wt = "#F8766D", alpha = "#7CAE00", delta = "#00BFC4", omicron = "#C77CFF")

## ---------------------------------------------------------------------------
## ppc_draws -- posterior predictive replicates for one variant
##
## For each posterior draw: simulate expected counts, draw negative binomial observations,
## collapse to the four marginals Figure 3 shows.
##
## The draw is made at WEEKLY resolution, matching the likelihood. The published figure
## drew daily and summed, and summing independent negative binomials gives variance
## M + xi*sum(mu_i^2) -- spreading the same mass over seven times as many cells shrinks the
## excess term about sevenfold. The published ribbons therefore understate the model's own
## predictive uncertainty, and points falling outside them read as misfit when they are
## partly an artefact. (U1)
##
## Draws are spread evenly across every chain, so the check reflects the whole posterior
## rather than one chain's corner of it.
## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
## sample_draws -- n posterior draws for a variant, spread evenly across its chains
##
## Returns the draws and the `x` list from the first chain -- identical across chains, since
## they share the data, the priors and every setting. Every figure draws through this, so
## each one reflects the whole posterior rather than one chain's corner of it.
## ---------------------------------------------------------------------------
sample_draws <- function(strain, n_samples, dir = STRAIN_DIR) {
  tags <- chain_tags(strain, dir)
  per  <- ceiling(n_samples / length(tags))
  dr   <- NULL; x <- NULL
  for (tg in tags) {
    o  <- readRDS(chain_file(strain, tg, dir))
    fl <- max(which(!is.na(o$MCMC_posteriors)))
    b  <- analysis_start(o)
    dr <- rbind(dr, on_natural_scale(o$MCMC_output[round(seq(b, fl, length.out = per)), , drop = FALSE], o$x))
    if (is.null(x)) x <- o$x
    rm(o); invisible(gc(verbose = FALSE))
  }
  list(draws = dr[seq_len(min(n_samples, nrow(dr))), , drop = FALSE], x = x)
}

ppc_draws <- function(strain, n_samples = 100, test_pop = NULL, ncores = 4,
                      dir = STRAIN_DIR) {

  load_kernels()

  sd_ <- sample_draws(strain, n_samples, dir)
  dr  <- sd_$draws; x <- sd_$x
  n_samples <- nrow(dr)
  if (is.null(test_pop)) test_pop <- x$test_pop

  obs_day <- x$data_array[-1, , , drop = FALSE]     # drop the undetected class
  obs_day[is.na(obs_day)] <- 0

  ct_cache_reset(ncores)
  for (i in seq_len(n_samples)) {
    p  <- dr[i, ]
    mu <- p_array_func(p, knots = x$knots, vk = 1, max_day = x$max_day,
                       population = x$population, data_array = x$data_array,
                       test_pop = test_pop, ncores = ncores, form = fit_form(x),
                       symp_delay_lim = x$symp_delay_lim, stoch = 0.5)[[1]]

    if (i == 1) {
      ## The NA pattern comes from start_j, not from the parameters, so compute it once.
      ## Weeks before the spline lead-in are dropped from draws and data alike, rather
      ## than zeroed and hidden behind a count filter.
      ok_day  <- !is.na(mu[1, , 1])
      n_wk    <- dim(weekly_aggregator(obs_day))[2]
      ok_wk   <- vapply(seq_len(n_wk),
                        function(w) all(ok_day[(7 * (w - 1) + 1):(7 * w)]), logical(1))
      obs     <- weekly_aggregator(obs_day)[, ok_wk, , drop = FALSE]
      wk_date <- as.Date(x$dates)[7 * (seq_len(n_wk) - 1) + 1][ok_wk]

      m_vl   <- matrix(NA_real_, n_samples, dim(obs)[1])
      m_wk   <- matrix(NA_real_, n_samples, dim(obs)[2])
      m_test <- matrix(NA_real_, n_samples, dim(obs)[3])
      m_2d   <- array(NA_real_, c(dim(obs)[1], dim(obs)[3], n_samples))
    }

    mu[is.na(mu)] <- 0
    mu_w <- weekly_aggregator(mu)[, ok_wk, , drop = FALSE]
    sim  <- array(rnbinom(length(mu_w), size = 1 / p[["rdisp"]], mu = mu_w), dim = dim(mu_w))

    m_vl[i, ]   <- apply(sim, 1, sum)
    m_wk[i, ]   <- apply(sim, 2, sum)
    m_test[i, ] <- apply(sim, 3, sum)
    m_2d[, , i] <- apply(sim, c(1, 3), sum)
    if (i %% 20 == 0) message("  draw ", i, "/", n_samples)
  }
  list(strain = strain, obs = obs, wk_date = wk_date,
       m_vl = m_vl, m_wk = m_wk, m_test = m_test, m_2d = m_2d)
}

## ---------------------------------------------------------------------------
## figure_phase_ppc -- posterior predictive check of the epidemic-phase pattern (R2.5)
##
##   figure_phase_ppc()
##   figure_phase_ppc(n_samples = 100, file = "outputs/figures/figS_phase_ppc.png")
##
## Reviewer 2 asks whether the fitted model can reproduce the growing/declining pattern of
## Figure 2 through posterior predictive simulation, and in particular Delta's departure from
## it. Figure 4 cannot answer that: it aggregates over the whole study period, so the phase
## contrast is averaged away.
##
## This is a stronger test than Figure 4. Phase bias is not fitted -- it emerges from
## integrating a time-varying incidence curve against the within-host trajectory -- so
## reproducing it is a check of the paper's central mechanism rather than of goodness of fit.
##
## The phase classification is the one used for the published Figure 2
## (archive/results_figures.R): a centred 7-day rolling mean of NATIONAL daily cases, growth
## taken as the ratio of that mean three days ahead to three days behind, and days labelled
## Growing above 1.1, Declining below 0.9, neutral otherwise.
##
## One deliberate difference from the archived phase_shift_plotter(): the lowest viral load row
## is dropped. It is the undetected class, which the likelihood drops and the predicted array
## never contains, so including it in an observed mean -- as the archived code did -- would
## compare unlike with unlike.
##
## The SAME day-to-phase map is applied to observed and predicted, so any difference between
## them is attributable to the model and not to a different split of the calendar.
##
## Means are computed from the expected counts rather than from negative binomial draws. The
## quantity is a mean over millions of tests, so observation noise contributes negligibly to
## it; the ribbon is posterior uncertainty in the expected mean viral load.
CASES_CSV <- "data/UK_covid_testing_data.csv"

phase_by_date <- function(csv = CASES_CSV, up = 1.1, down = 0.9) {
  cs <- read.csv(csv, stringsAsFactors = FALSE)
  cs <- cs[order(cs$date), c("date", "newCasesBySpecimenDate")]
  roll <- as.numeric(stats::filter(cs$newCasesBySpecimenDate, rep(1 / 7, 7), sides = 2))
  n <- length(roll)
  g <- c(roll[-seq_len(3)], rep(NA_real_, 3)) / c(rep(NA_real_, 3), roll[seq_len(n - 3)])
  setNames(ifelse(is.na(g), NA_character_,
                  ifelse(g > up, "Growing", ifelse(g < down, "Declining", "neutral"))),
           cs$date)
}

## phase_by_incidence -- the same growth rule applied to a variant's OWN fitted incidence
##
##   phase_by_incidence("delta")
##
## phase_by_date() classifies calendar days from the NATIONAL all-variant case series, which is
## what the published Figure 2 used. That is fine as a label for a calendar day and it keeps the
## two figures comparable, but it is not a statement about the variant being fitted, and the
## variant windows overlap: delta's runs to 2022-01-23, three months past the start of omicron's,
## so 31% of the days it labels "Growing" are days when national cases were climbing because of
## omicron while delta itself was in decline. Restricting delta to its pre-omicron days moves its
## observed growing-minus-declining contrast from -0.049 to +0.005 -- most of the "reversed"
## phase pattern is that contamination rather than a property of delta.
##
## Classifying by the variant's own fitted incidence removes the mixing, and is the sharper test
## of the mechanism besides: phase bias arises from convolving THAT incidence against the
## within-host trajectory, so this asks whether the model reproduces the bias implied by its own
## epidemic curve, rather than whether national cases happened to track the focal variant.
##
## ONE map per variant, from the posterior MEDIAN curve, not one per draw. The map has to be
## identical for the observed and predicted arrays or the comparison stops being like-for-like --
## a per-draw classification would re-split the calendar for every draw and the "observed" curve
## would move with the posterior.
##
## Onset incidence rather than the fitted crossing curve, because onset is what a case series
## counts and therefore what the national classifier approximates. Under thresh_tdist the two
## coincide; under the incidence form the crossing curve is convolved with the fitted onset
## delay, exactly as figure5() and p_array_func do it.
phase_by_incidence <- function(strain, dir = STRAIN_DIR, n_samples = 200L,
                               up = 1.1, down = 0.9) {
  if (!exists("infecteds_generator")) source("scripts/model.R")
  sd <- sample_draws(strain, n_samples, dir)
  x  <- sd$x; fm <- fit_form(x)
  L  <- if (fm == "incidence") x$symp_delay_lim else 0L

  inc <- t(vapply(seq_len(nrow(sd$draws)), function(j) {
    p <- sd$draws[j, ]
    I <- as.numeric(infecteds_generator(p, x$knots, x$population, x$max_day, form = fm))
    if (fm != "incidence") return(I)
    f <- stats::plnorm(0:L + 0.5, p[["inc1"]], p[["inc2"]]) -
         stats::plnorm(pmax(0:L - 0.5, 0), p[["inc1"]], p[["inc2"]])
    f <- f / sum(f)
    O <- rep(NA_real_, x$max_day)
    for (d in (L + 1L):x$max_day) O[d] <- sum(I[d - 0:L] * f)
    O
  }, numeric(x$max_day)))

  med  <- apply(inc, 2, median, na.rm = TRUE)
  roll <- as.numeric(stats::filter(med, rep(1 / 7, 7), sides = 2))
  n    <- length(roll)
  g    <- c(roll[-seq_len(3)], rep(NA_real_, 3)) / c(rep(NA_real_, 3), roll[seq_len(n - 3)])
  setNames(ifelse(is.na(g), NA_character_,
                  ifelse(g > up, "Growing", ifelse(g < down, "Declining", "neutral"))),
           x$dates)
}

## phase_by_counts -- the same growth rule applied to the VARIANT'S OWN OBSERVED case series
##
##   phase_by_counts("delta")
##
## The third classifier, and the one to prefer. The data array is (viral load x calendar day x
## test day), so summing over the first and third dimensions gives daily positive tests by
## specimen date for that variant alone -- the same quantity the national series holds, with the
## other variants removed. Nothing is fitted.
##
## Why this beats both alternatives:
##
##   * phase_by_date() mixes variants. Delta's window runs three months past the start of
##     omicron's, so 31% of the days it calls "Growing" for delta are days when national cases
##     were climbing because of omicron.
##
##   * phase_by_incidence() is partly self-referential. Phase bias arises from convolving the
##     fitted incidence against the fitted trajectory, so defining the phases from that same
##     fitted incidence and then asking whether the model reproduces the bias lets one half of
##     the model grade the other. It is the right diagnostic for asking what the model implies;
##     it is the wrong one for asking whether the model matches the world.
##
## The lowest viral load row is dropped, as everywhere else in this figure: it is the undetected
## class, and a case series counts positives.
##
## The day index is the SPECIMEN date in both this and the national series, so the two are the
## same kind of quantity and the 1.1 / 0.9 thresholds keep their meaning.
phase_by_counts <- function(strain, dir = STRAIN_DIR, up = 1.1, down = 0.9) {
  o <- readRDS(chain_file(strain, chain_tags(strain, dir)[1], dir))
  a <- o$x$data_array[-1, , , drop = FALSE]
  dt <- dimnames(o$x$data_array)[[2]]
  rm(o); invisible(gc(verbose = FALSE))
  a[is.na(a)] <- 0
  daily <- apply(a, 2, sum)

  roll <- as.numeric(stats::filter(daily, rep(1 / 7, 7), sides = 2))
  n    <- length(roll)
  g    <- c(roll[-seq_len(3)], rep(NA_real_, 3)) / c(rep(NA_real_, 3), roll[seq_len(n - 3)])
  setNames(ifelse(is.na(g), NA_character_,
                  ifelse(g > up, "Growing", ifelse(g < down, "Declining", "neutral"))), dt)
}

## Mean viral load by test day, within each phase, from a (vl x day x testday) count array.
mean_vl_by_phase <- function(arr, vl, ph, phases = c("Growing", "Declining")) {
  do.call(rbind, lapply(phases, function(f) {
    j <- which(ph == f)
    if (!length(j)) return(NULL)
    a <- arr[, j, , drop = FALSE]
    num <- apply(a, 3, function(m) sum(vl * rowSums(m)))
    den <- apply(a, 3, sum)
    data.frame(phase = f, delay = seq_len(dim(arr)[3]) - 1L,
               mean_vl = num / den, row.names = NULL)
  }))
}

## `classify` picks the phase map:
##   "counts"    -- the variant's own observed daily positives. The default, and the only one
##                  that is both variant-specific and independent of the fit.
##   "national"  -- phase_by_date(), the published Figure 2 rule. Keep it when the two figures
##                  have to agree; it mixes variants during transitions.
##   "incidence" -- phase_by_incidence(), the variant's own FITTED onset curve. Answers what the
##                  model implies, not whether it matches the data.
## `up` and `down` are the growth thresholds, exposed because the classification is a judgement
## and the figure should be reproducible at other cuts. 1.1 / 0.9 reproduces the published rule.
## Tightening them keeps only days whose direction is unambiguous: the observed contrast rises
## for every variant as they tighten, which is the signature of marginal days diluting it, so a
## stricter cut is a lower-variance estimate of the phase bias at the cost of fewer days.
figure_phase_ppc <- function(strains = STRAINS, n_samples = 50, ncores = 4,
                             classify = c("counts", "national", "incidence"),
                             up = 1.1, down = 0.9, titles = FALSE,
                             dir = STRAIN_DIR, file = NULL) {
  load_kernels()
  classify <- match.arg(classify)
  dir_of <- dir_resolver(dir, strains)

  pmap <- if (classify == "national") phase_by_date(up = up, down = down) else NULL
  obs_all <- pred_all <- NULL
  for (s in strains) {
    message(s, " <- ", dir_of(s), " [", classify, "]")
    sd_ <- sample_draws(s, n_samples, dir_of(s))
    dr  <- sd_$draws; x <- sd_$x

    obs_day <- x$data_array[-1, , , drop = FALSE]
    obs_day[is.na(obs_day)] <- 0
    vl <- as.numeric(sub("^vl_", "", dimnames(x$data_array)[[1]][-1]))

    ## The variant-specific map is built from its own chain, so it is derived per strain rather
    ## than once outside the loop.
    smap <- switch(classify,
                   national  = pmap,
                   counts    = phase_by_counts(s, dir_of(s), up = up, down = down),
                   incidence = phase_by_incidence(s, dir_of(s), up = up, down = down))
    ph <- unname(smap[dimnames(x$data_array)[[2]]])

    ct_cache_reset(ncores)
    for (i in seq_len(nrow(dr))) {
      p  <- dr[i, ]
      mu <- p_array_func(p, knots = x$knots, vk = 1, max_day = x$max_day,
                         population = x$population, data_array = x$data_array,
                         test_pop = x$test_pop, ncores = ncores, form = fit_form(x),
                         symp_delay_lim = x$symp_delay_lim, stoch = 0.5)[[1]]
      if (i == 1) {
        ## Days before the spline lead-in carry no prediction; drop them from both sides.
        ok <- !is.na(mu[1, , 1])
        o  <- mean_vl_by_phase(obs_day[, ok, , drop = FALSE], vl, ph[ok])
        o$strain <- unname(STRAIN_LABEL[s]); obs_all <- rbind(obs_all, o)
      }
      mu[is.na(mu)] <- 0
      q <- mean_vl_by_phase(mu[, ok, , drop = FALSE], vl, ph[ok])
      q$strain <- unname(STRAIN_LABEL[s]); q$draw <- i
      pred_all <- rbind(pred_all, q)
      if (i %% 20 == 0) message("  draw ", i, "/", nrow(dr))
    }
  }

  pr <- do.call(rbind, lapply(split(pred_all, list(pred_all$strain, pred_all$phase,
                                                   pred_all$delay), drop = TRUE),
    function(z) data.frame(strain = z$strain[1], phase = z$phase[1], delay = z$delay[1],
                           mid = median(z$mean_vl),
                           lo  = quantile(z$mean_vl, 0.025, names = FALSE),
                           hi  = quantile(z$mean_vl, 0.975, names = FALSE),
                           row.names = NULL)))

  lv <- unname(STRAIN_LABEL[strains])
  for (d in list(obs_all, pr)) d$strain <- factor(d$strain, levels = lv)
  obs_all$strain <- factor(obs_all$strain, levels = lv)
  pr$strain      <- factor(pr$strain,      levels = lv)

  ## STRAIN_COL is keyed on the code, the data on the label.
  pal <- setNames(unname(STRAIN_COL[strains]), unname(STRAIN_LABEL[strains]))

  p <- ggplot(pr, aes(delay)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = strain, group = phase), alpha = 0.25) +
    geom_line(aes(y = mid, colour = strain, linetype = phase), linewidth = 0.7) +
    geom_line(data = obs_all, aes(y = mean_vl, linetype = phase), colour = "grey25",
              linewidth = 0.4) +
    geom_point(data = obs_all, aes(y = mean_vl, shape = phase), size = 1.7,
               colour = "grey15", fill = "white") +
    scale_colour_manual(values = pal, guide = "none") +
    scale_fill_manual(values = pal, guide = "none") +
    scale_linetype_manual(values = c(Growing = 1, Declining = 2), name = NULL) +
    scale_shape_manual(values = c(Growing = 21, Declining = 4), name = NULL) +
    facet_wrap(~strain, nrow = 1) +
    labs(x = "Time from symptom onset (days)",
         y = expression(paste("Mean viral load (", log[10], " copies mL"^-1, ")")),
         ## titles = FALSE by default: this is a manuscript figure and the caption carries the
         ## description. That matters more than usual here, because the subtitle named the phase
         ## CLASSIFIER -- the three versions of this figure look alike and differ only in how
         ## the calendar was split. With the titles off, the caption has to say which classifier
         ## was used, or the figure is ambiguous in a way the reader cannot detect.
         ## Pass titles = TRUE for a standalone diagnostic copy, where the labelling is the
         ## only thing telling two otherwise identical files apart.
         title = if (!titles) NULL else
           "Posterior predictive check of the epidemic-phase pattern",
         subtitle = if (!titles) NULL else if (fit_form(x) == "incidence")
           paste0("grey = observed, coloured = posterior predictive median and 95% interval.",
                  "\nPhases from ",
                  switch(classify,
                         national  = "national all-variant case counts",
                         counts    = "each variant's own observed daily positive tests",
                         incidence = "each variant's own fitted onset incidence"),
                  sprintf(": 7-day mean, growth over 6 days, >%.2f growing, <%.2f declining.",
                          up, down),
                  "\nIncidence indexed on crossing the viral load threshold: the delay to onset",
                  " is convolved\nagainst calendar-time incidence, so predicted composition can",
                  " differ between phases.")
         else paste("grey = observed, coloured = posterior predictive median and 95%",
                    "interval.\nThe two predictive curves coincide exactly: conditional",
                    "on delay the model's viral load\ndistribution does not depend on",
                    "epidemic phase.")) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold"),
          axis.text = element_text(colour = "black"),
          ## Below the panels. With the title gone the legend was the topmost element, which put
          ## a floating key above the variant names and made the strip labels read as if they
          ## belonged to it. Underneath it sits next to the caption, where a reader looks for
          ## the key anyway. Kept at the top when titles are on, so the diagnostic copy still
          ## reads title -> subtitle -> key -> panels.
          legend.position = if (titles) "top" else "bottom",
          legend.margin = margin(t = 2, b = 0))

  ## The four lines of title and subtitle took about 0.9 in. Reclaimed as panel height rather
  ## than shrinking the canvas, since this is now a main-text figure.
  if (!is.null(file)) ggsave(file, p, width = 9, height = if (titles) 3.2 else 2.9)
  print(p)
  invisible(list(plot = p, observed = obs_all, predicted = pr))
}

## Column headings replace the per-panel y-axis titles: every panel counts positive tests,
## so "Count" repeated sixteen times said nothing the x-axis did not. The panel letter is
## folded into the heading rather than added by ggarrange, so the top row carries one text
## element instead of two. Only the top row is titled.
## Kept short deliberately. At width 8 the first three columns are about 1.7 inches, which
## holds roughly 20 characters of 9pt bold -- longer titles are silently clipped, and the
## clipping is easy to miss because panel D is 1.4x wider and its title survives. Each
## panel's x-axis already names the quantity, so the title only has to say which marginal
## it is. Units go in the caption.
PANEL_TITLES  <- c("A  Viral load", "B  Tests per week",
                   "C  Testing delay", "D  Load by delay")
PANEL_MAX_DAY <- 5   # panel D shows days 0-5, which fills a 2x3 grid exactly
VL_LAB        <- expression(paste("Viral load (", log[10], ")"))

fit_panels <- function(d, titles = FALSE, xlabs = FALSE) {
  ## x-axis titles only on the bottom row -- the four columns repeat identically down the
  ## figure. Tick labels stay on every row, since panel B's date range differs by variant.
  xl  <- function(s) if (xlabs) s else NULL
  col <- STRAIN_COL[[d$strain]]

  base <- function(p) p + theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black"),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 9, face = "bold", hjust = 0,
                                    margin = margin(b = 4))) +
    geom_point(aes(y = Count), colour = col, size = 0.9) +
    geom_line(aes(y = median), colour = col) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = col, alpha = 0.3)

  ## A 2.5% quantile of zero is legitimate -- the model can predict no counts in a cell --
  ## but log10(0) is -Inf, which ggplot drops with a warning and leaves the ribbon torn.
  ## Clamped to half a count, which is below anything observable and so cannot mislead.
  mk <- function(m, x) {
    q <- matrixStats::colQuantiles(m, probs = c(.025, .5, .975))
    q[, 1] <- pmax(q[, 1], 0.5)
    setNames(data.frame(x = x, q), c("x", "lower", "median", "upper"))
  }

  ## The array's first dimension is labelled "vl_1.5", "vl_2", ... -- as.numeric() on those
  ## returns NA, which silently blanks every panel indexed on viral load.
  vl_vals <- as.numeric(sub("^vl_", "", dimnames(d$obs)[[1]]))

  vl <- cbind(mk(d$m_vl, vl_vals),                      Count = apply(d$obs, 1, sum, na.rm = TRUE))
  dy <- cbind(mk(d$m_wk, d$wk_date),                    Count = apply(d$obs, 2, sum, na.rm = TRUE))
  td <- cbind(mk(d$m_test, seq_len(dim(d$obs)[3]) - 1), Count = apply(d$obs, 3, sum, na.rm = TRUE))
  dy <- dy[dy$Count > 0, ]

  q2   <- apply(d$m_2d, c(1, 2), quantile, probs = c(.025, .5, .975))
  twod <- expand.grid(vl = vl_vals, td = seq_len(dim(d$obs)[3]) - 1)
  twod <- cbind(twod, lower = pmax(c(q2[1, , ]), 0.5), median = c(q2[2, , ]),
                upper = c(q2[3, , ]),
                Count = c(apply(d$obs, c(1, 3), sum, na.rm = TRUE)))

  ## Days 0-5 fill a 2x3 grid exactly; all seven leave two empty cells and a stranded
  ## panel. Day 6 carries the fewest observations, so it is the least informative to omit
  ## -- but say so in the caption, since the model is fitted to all seven.
  twod       <- twod[twod$td <= PANEL_MAX_DAY, ]
  twod$panel <- factor(paste("Day", twod$td))
  ttl        <- if (titles) PANEL_TITLES else rep("", 4)

  list(
    base(ggplot(vl, aes(x))) +
      scale_y_log10(labels = scales::label_number()) +
      labs(title = ttl[1], x = xl(VL_LAB)),
    ## breaks_pretty adapts to the range: a fixed "3 months" gives omicron's ~13-week
    ## window no break at all, and crowds alpha's into an unreadable run.
    base(ggplot(dy, aes(x))) +
      scale_y_log10(labels = scales::label_number()) +
      scale_x_date(labels = scales::date_format("%b-%y"),
                   breaks = scales::breaks_pretty(n = 3)) +
      labs(title = ttl[2], x = xl("Week commencing")),
    base(ggplot(td, aes(x))) +
      scale_y_continuous(labels = scales::label_number()) +
      labs(title = ttl[3], x = xl("Day of test after onset")),
    base(ggplot(twod, aes(x = vl))) +
      facet_wrap(~panel, ncol = 3) +
      scale_y_log10(labels = scales::label_number(), breaks = 10^seq(0, 6, 2)) +
      theme(strip.background = element_blank(),
            strip.text = element_text(colour = "black", hjust = 0, size = 8,
                                      margin = margin(b = 2)),
            panel.spacing = grid::unit(4, "pt")) +
      labs(title = ttl[4], x = xl(VL_LAB))
  )
}

## ---------------------------------------------------------------------------
## figure3 -- posterior predictive check, one row per variant
##
##   figure3(n_samples = 100, file = "images/figure3.png")
## ---------------------------------------------------------------------------
figure3 <- function(strains = STRAINS, n_samples = 100, file = NULL,
                    dir = STRAIN_DIR, ncores = 4) {
  dir_of <- dir_resolver(dir, strains)
  rows <- lapply(seq_along(strains), function(k) {
    s <- strains[k]; message(s, " <- ", dir_of(s))
    r <- ggpubr::ggarrange(
      plotlist = fit_panels(ppc_draws(s, n_samples, ncores = ncores, dir = dir_of(s)),
                            titles = (k == 1), xlabs = (k == length(strains))),
      nrow = 1, ncol = 4, widths = c(1, 1, 1, 1.4))
    ggpubr::annotate_figure(r, left = ggpubr::text_grob(STRAIN_LABEL[[s]], face = "bold",
                                                        rot = 90, size = 10))
  })
  ## The titled row and the labelled row each need extra height so every row's plot area
  ## comes out the same.
  h <- rep(1, length(strains)); h[c(1, length(h))] <- h[c(1, length(h))] + 0.08
  p <- ggpubr::ggarrange(plotlist = rows, nrow = length(strains), heights = h)
  if (!is.null(file)) ggsave(file, p, width = 8, height = 1.8 * length(strains),
                             limitsize = FALSE)
  p
}

## ---------------------------------------------------------------------------
## FIGURE 4 -- kinetics parameters and derived intervals
## ---------------------------------------------------------------------------
FIG4 <- list(
  growth         = "A  Growth rate",
  decay          = "B  Decay rate",
  peak           = "C  Peak viral load",
  det_to_peak    = "D  Detectable to peak",
  peak_to_end    = "E  Peak to undetectable",
  total          = "F  Total detectable time",
  det_to_symp    = "G  Detectable to symptoms",
  symp_to_end    = "H  Symptoms to undetectable",
  thresh_to_symp = "I  Threshold to symptoms")

## One posterior draw -> one simulated population -> the six duration quantities.
##
## vl0 is the detection threshold: 1.5 log10 is the paper's, 2.0 is the 100 copies/mL that
## R1.2a asks for. Making it an argument also settles the vi=1 / vi=2 divergence between
## the two archived copies of vls_calculator -- whichever was intended, it is now explicit
## and the same value flows through every panel.
fig4_draw <- function(p, n_indiv, vl0) {
  a  <- 0.5  + exp(rnorm(n_indiv, p[["a_bar"]], p[["a_sigma"]]))
  b  <- 0.25 + exp(rnorm(n_indiv, p[["b_bar"]], p[["b_sigma"]]))
  vp <- rnorm(n_indiv, p[["vl_max_bar"]], p[["vl_max_sigma"]])
  Ti <- exp(rnorm(n_indiv, p[["inc1"]], p[["inc2"]]))
  vt <- p[["inc3"]]

  ## Never detectable above the threshold: excluded, as the kernel does.
  ##
  ## A parameter set where almost nobody clears the threshold yields quantiles computed
  ## from a handful of individuals, or none at all. That never happens under the posterior
  ## but is common under the prior, which is diffuse enough to propose peaks below the
  ## symptom threshold. NULL rather than a row of NaN, so the caller can count how many
  ## draws were unusable instead of silently propagating them.
  ok <- vp > vt & vp > vl0
  if (sum(ok) < 100) return(NULL)
  a <- a[ok]; b <- b[ok]; vp <- vp[ok]; Ti <- Ti[ok]

  rise <- vp - vl0
  dtp  <- rise / a
  pte  <- rise / b
  dts  <- (vt - vl0) / a + Ti

  data.frame(det_to_peak = dtp, peak_to_end = pte, total = dtp + pte,
             det_to_symp = dts, symp_to_end = dtp + pte - dts, thresh_to_symp = Ti,
             ## Table 4's "Symptom onset -> peak". Negative when symptoms follow the peak,
             ## which is the usual case for wild type (published median -1.93 d), so it must
             ## not be clamped at zero.
             symp_to_peak = dtp - dts)
}

## Growth rate, decay rate and peak viral load have closed-form population quantiles, so
## they are computed exactly rather than simulated -- which also avoids a real bias: a and
## b are shifted lognormals, and at a_sigma near 1.3 the SAMPLE mean is itself strongly
## right-skewed, so its median across posterior draws sits well below the true mean.
##
## The central point is the population MEDIAN, not the mean. For a lognormal with sigma
## above 1 the mean lies close to the upper quartile -- at wild type's a_sigma of 1.33 it
## falls near the 75th percentile -- so plotting it against an interquartile whisker would
## put the point at the whisker's right-hand end and misrepresent the typical individual.
## A median with quartiles is one coherent family of summary.
fig4_closed <- function(p) {
  lq <- function(shift, mu, sg, q) shift + exp(mu + sg * qnorm(q))
  cbind(
    growth = c(mid = 0.5  + exp(p[["a_bar"]]),
               q25 = lq(0.5,  p[["a_bar"]], p[["a_sigma"]], 0.25),
               q75 = lq(0.5,  p[["a_bar"]], p[["a_sigma"]], 0.75)),
    decay  = c(mid = 0.25 + exp(p[["b_bar"]]),
               q25 = lq(0.25, p[["b_bar"]], p[["b_sigma"]], 0.25),
               q75 = lq(0.25, p[["b_bar"]], p[["b_sigma"]], 0.75)),
    peak   = c(mid = p[["vl_max_bar"]],
               q25 = qnorm(0.25, p[["vl_max_bar"]], p[["vl_max_sigma"]]),
               q75 = qnorm(0.75, p[["vl_max_bar"]], p[["vl_max_sigma"]])))
}

## Per draw: the population median and IQR. Across draws: the 95% credible interval of each.
## Two levels of uncertainty, kept distinct -- the whisker is between-individual variation,
## the shaded rectangles are posterior uncertainty about where that variation sits.
## ---------------------------------------------------------------------------
## prior_draws -- draws from the prior, on the same footing as posterior draws
##
## The prior has two parts, matching prior() in model.R exactly:
##   parameters 1:8 (less the ignored l1/l_sigma) -- a TRUNCATED multivariate normal from
##     the longitudinal fit, so the correlations between growth, decay and peak are kept
##   inc1-inc3 -- independent truncated normals, mean and sd from prior_cov_final's diagonal
##
## Truncated, not plain, MVN. The archived prior_vs_posterior sampled the untruncated one
## while the sampler evaluates a truncated density, so its prior curves were subtly wrong
## -- most visibly for a_sigma and b_sigma, which are bounded below at zero and whose prior
## mass would otherwise extend into negative values the model can never visit.
## ---------------------------------------------------------------------------
## `bounds` overrides x$lower_bound for the named parameters, for the case where the stored
## bound is looser than the model's actual support. rdisp is the one: its lower_bound is -Inf,
## so the prior as specified is N(1, 10) with about 46% of its mass on negative values -- which
## the LIKELIHOOD then excludes, since the negative binomial needs a positive size. Sampling
## the unbounded prior therefore draws a region the posterior can never occupy and makes the
## panel span -18 to +20. Passing c(rdisp = 0) draws the prior restricted to the support the
## data actually see, which is what the posterior is being compared against.
##
## This changes nothing about the fit. Truncating a prior at a point the likelihood already
## forbids alters only its normalising constant, which cancels in every acceptance ratio.
prior_draws <- function(strain, n_samples, dir = STRAIN_DIR, pars = NULL, bounds = NULL) {
  o  <- readRDS(chain_file(strain, chain_tags(strain, dir)[1], dir))
  x  <- o$x; rm(o); invisible(gc(verbose = FALSE))

  if (length(bounds)) {
    bad <- setdiff(names(bounds), names(x$lower_bound))
    if (length(bad)) stop("no such parameter to bound: ", paste(bad, collapse = ", "))
    ## Tighten only. A bound looser than the model's own would draw a prior the sampler never
    ## had, which is a different object from the one this figure claims to show.
    loose <- names(bounds)[bounds < x$lower_bound[names(bounds)]]
    if (length(loose))
      stop("bounds may only tighten; ", paste(loose, collapse = ", "), " already bounded higher")
    x$lower_bound[names(bounds)] <- bounds
  }

  nm  <- names(x$parameters)
  vl  <- setdiff(nm[1:8], nm[x$ignored_pars])
  if (is.null(pars)) pars <- c(vl, "inc1", "inc2", "inc3")

  ## prior() gives the multivariate block to parameters 1:8 and independent truncated
  ## normals to 9:(length-1). The final parameter -- the multiplier -- falls outside both,
  ## so it has no proper prior and cannot be sampled from one.
  last <- length(nm) - 1L
  ind  <- setdiff(pars, vl)
  bad  <- setdiff(ind, nm[9:last])
  if (length(bad))
    stop("no prior is defined for: ", paste(bad, collapse = ", "))

  m <- matrix(NA_real_, n_samples, length(pars), dimnames = list(NULL, pars))
  kv <- intersect(pars, vl)
  if (length(kv))
    m[, kv] <- TruncatedNormal::rtmvnorm(n_samples, mu = x$prior_mean[kv],
                                         sigma = x$prior_cov_final[kv, kv],
                                         lb = x$lower_bound[kv],
                                         ub = rep(Inf, length(kv)))
  sds <- sqrt(diag(x$prior_cov_final))
  for (p in ind)
    m[, p] <- truncnorm::rtruncnorm(n_samples, a = x$lower_bound[[p]],
                                    mean = x$prior_mean[[p]], sd = sds[[p]])
  ## The bounds travel with the draws. Anything that estimates a density from these samples
  ## needs them: a Gaussian kernel put on a sample that piles up against a hard bound spreads
  ## mass across it and reports a peak INSIDE the support where the true density peaks at the
  ## bound itself. inc3 is the case that matters -- its prior mean is exactly its lower bound,
  ## so the prior is a half-normal whose mode is at zero, and an uncorrected estimate draws it
  ## with a mode near 0.4 and a dip to zero at the boundary.
  attr(m, "lower_bound") <- setNames(as.numeric(x$lower_bound[pars]), pars)
  m
}

## `draws` overrides the posterior sample, which is how the prior is put through identical
## machinery rather than summarised by a parallel code path that could differ.
fig4_summary <- function(strain, n_samples = 300, n_indiv = 20000, vl0 = 1.5,
                         dir = STRAIN_DIR, draws = NULL) {
  dr <- if (is.null(draws)) sample_draws(strain, n_samples, dir)$draws else draws

  ## The durations are ratios of the kinetics parameters and have no closed form, so they
  ## are simulated; n_indiv is large because they inherit the same heavy tails.
  s <- lapply(seq_len(nrow(dr)), function(j) {
    p <- dr[j, ]
    d <- fig4_draw(p, n_indiv, vl0)
    if (is.null(d)) return(NULL)
    cbind(fig4_closed(p),
          rbind(mid = apply(d, 2, median),
                q25 = apply(d, 2, quantile, 0.25),
                q75 = apply(d, 2, quantile, 0.75)))
  })

  drop <- vapply(s, is.null, logical(1))
  if (all(drop))
    stop("no usable draws for ", strain, ": every parameter set left fewer than 100 ",
         "individuals detectable above ", vl0)
  if (any(drop))
    message(sprintf("  %s: %d of %d draws dropped -- fewer than 100 individuals detectable",
                    strain, sum(drop), length(s)))
  s <- s[!drop]

  grab <- function(row) t(vapply(s, function(m) m[row, ], numeric(ncol(s[[1]]))))
  q    <- function(m) t(apply(m, 2, quantile, c(0.025, 0.5, 0.975)))

  mn <- q(grab("mid")); lo <- q(grab("q25")); hi <- q(grab("q75"))
  data.frame(strain = strain, param = rownames(mn),
             ## population median and the 95% CrI on it
             mid = mn[, 2], mid_lo = mn[, 1], mid_hi = mn[, 3],
             ## posterior median of each quartile (the whisker) and the 95% CrI on each
             iqr_lo = lo[, 2], iqr_hi = hi[, 2],
             q25_lo = lo[, 1], q25_hi = lo[, 3],
             q75_lo = hi[, 1], q75_hi = hi[, 3], row.names = NULL)
}

## ---------------------------------------------------------------------------
## figure4_compare -- Figure 4 with two sets of fits on the same axes
##
##   figure4_compare(file = "outputs/figures/fig4_compare.png")   # the two model forms
##   figure4_compare(vl0 = 2)
##   figure4_priors(file = "outputs/figures/incidence/fig4_priors.png")   # the two priors
##
## Each variant gets two rows, offset above and below its label. Everything else -- the point
## at the population median, the interquartile whisker, the three credible-interval bands --
## is encoded exactly as in figure4(), so the two can be read side by side.
##
## Two comparisons use this. The model-form one holds the prior fixed and varies the indexing;
## the prior one holds the model fixed and varies the prior width. `prior_mode` says which,
## and it is a hard check rather than a label:
##
##   identical  the two sets must share a prior exactly, or the comparison would confound the
##              thing being varied with a prior that had drifted between directories. One grey
##              prior row is drawn.
##   scaled     the two priors must differ by one common factor on every parameter. Both prior
##              rows are drawn, in their specification's colour, so the reader can see how far
##              the prior moved against how far the posterior did -- which is the entire point
##              of a prior-sensitivity figure.
##   none       no prior row.
##
## Offsets on a continuous axis rather than position_dodge, because the bands are geom_tile
## with an explicit height and dodging would move the tiles and the whiskers by different
## amounts -- the same reason figure_prior_sensitivity() offsets by hand.
## ---------------------------------------------------------------------------
FORM_LABEL <- c(thresh_tdist = "Onset-indexed (reported)", incidence = "Incidence-indexed")
FORM_COL   <- setNames(c("grey25", "#D95F02"), FORM_LABEL)

PRIOR_LABEL <- c(std = "Original priors", wide = "Widened priors")

figure4_compare <- function(strains = STRAINS, vl0 = 1.5, n_samples = 300, n_indiv = 20000,
                            file = NULL, show_prior = TRUE, n_prior = 1000, offset = 0.19,
                            width = 8.4, height = 8.4,
                            dirs = NULL, labels = FORM_LABEL,
                            prior_mode = c("identical", "scaled", "none")) {
  prior_mode <- match.arg(prior_mode)
  if (length(labels) != 2L) stop("`labels` must name exactly two sets of fits")
  forms <- names(labels)
  if (is.null(dirs))
    dirs <- vapply(forms, function(f) FORM_DIRS[[f]][["strains"]], character(1))
  if (is.null(names(dirs))) names(dirs) <- forms
  dirs <- dirs[forms]
  miss  <- dirs[!dir.exists(dirs)]
  if (length(miss)) stop("no chains in ", paste(miss, collapse = ", "))

  df <- do.call(rbind, lapply(forms, function(f) {
    d <- do.call(rbind, lapply(strains, function(s) {
      message(f, " ", s)
      fig4_summary(s, n_samples = n_samples, n_indiv = n_indiv, vl0 = vl0, dir = dirs[[f]])
    }))
    d$model <- unname(labels[f]); d
  }))

  pri <- NULL
  if (show_prior && prior_mode != "none") {
    KP <- c("a_bar", "a_sigma", "b_bar", "b_sigma", "vl_max_bar", "vl_max_sigma",
            "inc1", "inc2", "inc3")
    px <- lapply(forms, function(f) {
      o <- readRDS(chain_file(strains[1], chain_tags(strains[1], dirs[[f]])[1], dirs[[f]]))
      r <- list(mean = o$x$prior_mean[KP], sd = o$x$prior_sds[KP])
      rm(o); invisible(gc(verbose = FALSE)); r
    })
    if (any(abs(px[[1]]$mean - px[[2]]$mean) > 1e-10))
      stop("the two sets of fits used different prior MEANS -- the comparison would confound ",
           "the quantity being varied with a relocated prior")
    if (prior_mode == "identical") {
      if (any(abs(px[[1]]$sd - px[[2]]$sd) > 1e-10))
        stop("the two sets of fits used different priors -- a model comparison would ",
             "confound the two; use prior_mode = \"scaled\" if that is intended")
      pri <- fig4_summary(strains[1], n_indiv = n_indiv, vl0 = vl0, dir = dirs[[1]],
                          draws = prior_draws(strains[1], n_prior, dirs[[1]]))
      pri$strain <- "Prior"; pri$model <- "Prior"
    } else {
      ## One common factor, or the panels are not comparing a single manipulation.
      r <- px[[2]]$sd / px[[1]]$sd
      if (diff(range(r)) > 1e-6)
        stop("the prior widening is not a single factor across parameters (",
             sprintf("%.3f-%.3f", min(r), max(r)), ") -- the figure would imply one")
      message(sprintf("prior standard deviations widened %.2fx on every parameter", r[[1]]))
      pri <- do.call(rbind, lapply(forms, function(f) {
        z <- fig4_summary(strains[1], n_indiv = n_indiv, vl0 = vl0, dir = dirs[[f]],
                          draws = prior_draws(strains[1], n_prior, dirs[[f]]))
        z$strain <- "Prior"; z$model <- unname(labels[f]); z
      }))
    }
    df <- rbind(df, pri)
  }

  df <- df[df$param %in% names(FIG4), , drop = FALSE]
  LAB <- c(STRAIN_LABEL, Prior = "Prior")
  ord <- c(rev(unname(STRAIN_LABEL[strains])), if (!is.null(pri)) "Prior")
  df$label <- factor(LAB[df$strain], levels = ord)
  df$param <- factor(df$param, levels = names(FIG4), labels = unlist(FIG4))
  df$model <- factor(df$model, levels = c(unname(labels), "Prior"))
  ## The first set sits above its variant's label, the second below. A single shared prior row
  ## has no specification and sits on the label; under prior_mode = "scaled" the two prior rows
  ## carry their own specification and are offset like everything else.
  df$y <- as.numeric(df$label) +
          ifelse(df$model == labels[[1]], offset,
          ifelse(df$model == labels[[2]], -offset, 0))

  cols <- c(setNames(unname(FORM_COL), unname(labels)), Prior = "grey55")
  BAND_PANELS <- levels(df$param)[1:3]

  build <- function(d, legend) {
    d <- droplevels(d)
    lv <- levels(d$label)
    ggplot(d, aes(y = y, colour = model, fill = model)) +
      geom_tile(aes(x = (mid_lo + mid_hi) / 2, width = mid_hi - mid_lo),
                height = 0.32, alpha = 0.3, colour = NA) +
      geom_tile(aes(x = (q25_lo + q25_hi) / 2, width = q25_hi - q25_lo),
                height = 0.32, alpha = 0.3, colour = NA) +
      geom_tile(aes(x = (q75_lo + q75_hi) / 2, width = q75_hi - q75_lo),
                height = 0.32, alpha = 0.3, colour = NA) +
      geom_errorbar(aes(xmin = iqr_lo, xmax = iqr_hi), orientation = "y",
                    width = 0.22, linewidth = 0.5) +
      geom_point(aes(x = mid), size = 1.5) +
      scale_y_continuous(breaks = seq_along(lv), labels = lv,
                         limits = c(0.5, length(lv) + 0.5)) +
      scale_colour_manual(values = cols, name = NULL, drop = TRUE) +
      scale_fill_manual(values = cols, name = NULL, guide = "none") +
      facet_wrap(~param, scales = "free_x", ncol = 3) +
      labs(x = NULL, y = NULL) +
      theme_bw(base_size = 9) +
      theme(panel.grid = element_blank(), panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            strip.background = element_blank(),
            strip.text = element_text(hjust = 0, face = "bold", size = 9),
            axis.text = element_text(colour = "black"),
            legend.position = if (legend) "top" else "none")
  }

  top <- df[df$param %in% BAND_PANELS, ]
  ## Filter the lower block on the STRAIN, not the model: under prior_mode = "scaled" the prior
  ## rows carry a specification label rather than "Prior", so a model-based filter would leave
  ## them in and the derived durations would acquire a prior row they never had.
  bot <- df[!(df$param %in% BAND_PANELS) & df$strain != "Prior", ]
  bot$label <- droplevels(bot$label)
  bot$y <- as.numeric(bot$label) + ifelse(bot$model == labels[[1]], offset, -offset)

  p <- if (nrow(bot)) ggpubr::ggarrange(build(top, TRUE), build(bot, FALSE), nrow = 2,
                                        heights = c(1, 1.6)) else build(top, TRUE)
  if (!is.null(file)) ggsave(file, p, width = width, height = height)
  p
}

## figure4_priors -- the same figure comparing the two PRIOR specifications under one model.
## Both directories follow the active form, so this cannot silently compare across forms.
figure4_priors <- function(strains = STRAINS, vl0 = 1.5, n_samples = 300, n_indiv = 20000,
                           file = NULL, ...) {
  figure4_compare(strains = strains, vl0 = vl0, n_samples = n_samples, n_indiv = n_indiv,
                  file = file, dirs = c(std = STRAIN_DIR, wide = WIDE_DIR),
                  labels = PRIOR_LABEL, prior_mode = "scaled", ...)
}

## ---------------------------------------------------------------------------
## figure4 -- population medians and interquartile ranges, all variants
##
##   figure4(file = "images/figure4.png")                          # 1.5 log10
##   figure4(vl0 = 2, file = "images/figure4_100copies.png")       # R1.2a
## ---------------------------------------------------------------------------
## n_prior is larger than n_samples because prior draws cost nothing -- no C++, no chain
## reads -- and the band is a quantile of a quantile, which is noisy at small n.
## `dir` may be one directory for every strain, or a named vector giving one per strain --
## see dir_resolver().
figure4 <- function(strains = STRAINS, vl0 = 1.5, n_samples = 300, n_indiv = 20000,
                    file = NULL, dir = STRAIN_DIR, show_prior = TRUE,
                    n_prior = 1000) {
  dir_of <- dir_resolver(dir, strains)

  df <- do.call(rbind, lapply(strains, function(s) {
    message(s, " <- ", dir_of(s))
    fig4_summary(s, n_samples = n_samples, n_indiv = n_indiv, vl0 = vl0, dir = dir_of(s))
  }))

  ## The prior's own interquartile range, through the same machinery. Drawn as a grey band
  ## behind the whiskers so a posterior that has learned nothing is immediately visible --
  ## it sits on top of the grey. Panel B needs this: Alpha's and Omicron's b_sigma has a
  ## posterior/prior SD ratio of 0.94 and 0.96, so their wide decay whiskers are the prior
  ## showing through, and without the band the figure presents them exactly as it presents
  ## wild type's and Delta's, which are estimated.
  ##
  ## Computed ONCE and drawn across the full panel height, because the prior is the same
  ## for every variant. Drawing it per row instead made four bands that differed only by
  ## Monte Carlo noise in a quantile-of-a-quantile, which reads as a real difference
  ## between variants. Asserted rather than assumed, since a strain added later with its
  ## own prior would silently break it.
  pri <- NULL
  if (show_prior) {
    pm <- sapply(strains, function(s) {
      o <- readRDS(chain_file(s, chain_tags(s, dir_of(s))[1], dir_of(s)))
      o$x$prior_mean[c("a_bar","a_sigma","b_bar","b_sigma","vl_max_bar","vl_max_sigma",
                       "inc1","inc2","inc3")] })
    if (any(abs(pm - pm[, 1]) > 1e-10))
      stop("priors differ between strains -- draw the band per row rather than once")
    pri <- fig4_summary(strains[1], n_indiv = n_indiv, vl0 = vl0, dir = dir_of(strains[1]),
                        draws = prior_draws(strains[1], n_prior, dir_of(strains[1])))
  }

  ## The prior enters as a further row with identical encoding -- point, whisker,
  ## rectangles -- rather than as a background band. The same visual grammar makes the
  ## comparison direct: a variant whose row looks like the prior row has learned nothing
  ## about that quantity. The rectangles on the prior row are prior uncertainty about the
  ## population quantiles, exactly as the coloured ones are posterior uncertainty.
  if (!is.null(pri)) { pri$strain <- "Prior"; df <- rbind(df, pri) }

  ## Keep only the quantities this figure draws.
  ##
  ## fig4_draw also returns symp_to_peak, which exists for Supplementary Table 6's "Symptom
  ## onset -> peak" row and has no entry in FIG4. Without this filter, factor() maps it to NA
  ## and ggplot renders an unlabelled tenth panel -- a whole extra row of the figure. Adding
  ## a quantity to fig4_draw should not silently add a panel here.
  df <- df[df$param %in% names(FIG4), , drop = FALSE]

  ## Variants run top to bottom in the order given, with the prior ABOVE them. ggplot puts
  ## the first level at the bottom, so the variant labels are reversed and "Prior" appended
  ## last. The colour values are NAMED rather than positional -- supplying them in strains
  ## order against reversed levels shifts every variant onto the next one's colour.
  lab  <- c(unname(STRAIN_LABEL[strains]), if (!is.null(pri)) "Prior")
  cols <- setNames(c(unname(STRAIN_COL[strains]), if (!is.null(pri)) "grey45"), lab)
  LAB  <- c(STRAIN_LABEL, Prior = "Prior")

  ord <- c(rev(unname(STRAIN_LABEL[strains])), if (!is.null(pri)) "Prior")
  df$strain <- factor(LAB[df$strain], levels = ord)
  df$param  <- factor(df$param, levels = names(FIG4), labels = unlist(FIG4))

  ## Credible-interval bands: every variant, and the prior on panels A-C only.
  ##
  ## Those three are the model's own parameters, where the prior is defined rather than
  ## implied, and including its bands costs nothing -- measured, they leave the decay and
  ## peak axes unchanged, and the growth axis actually needs LESS range for the prior than
  ## for the variants.
  ##
  ## The prior is omitted from the six derived durations entirely -- no bands, no whisker,
  ## no point. Its bands there are ruinous: they stretch panel G nearly sixfold (0.3 to
  ## 36.7 days against the variants' 2.6 to 8.8) and panel H runs to -15 days, since
  ## symp_to_end goes negative whenever the symptom delay exceeds the detectable period.
  ## Showing the whisker without the bands would mean two different encodings of the same
  ## row within one figure, which is worse than a clean absence: the prior is defined on
  ## the parameters, not on quantities derived from them, so those panels leave its row
  ## blank and the caption says why.
  ## Built as two stacked plots rather than one nine-panel facet. facet_wrap shares its
  ## discrete y scale across every panel, so a single figure would carry an empty "Prior"
  ## row through D-I. Splitting lets droplevels() remove that level where the prior is not
  ## shown, and the two blocks are reassembled so the result looks like one figure.
  BAND_PANELS <- levels(df$param)[1:3]

  build <- function(d) {
    d <- droplevels(d)
    ggplot(d, aes(y = strain, colour = strain, fill = strain)) +
      ## Three bands of equal height: the 95% CrI of the population median, and of each
      ## quartile. geom_tile rather than geom_rect so the vertical position comes from the
      ## same discrete mapping as the whiskers -- computing ymin and ymax by hand from
      ## as.numeric(strain) silently drew every band one row from its own whisker.
      geom_tile(aes(x = (mid_lo + mid_hi) / 2, width = mid_hi - mid_lo),
                height = 0.6, alpha = 0.3, colour = NA) +
      geom_tile(aes(x = (q25_lo + q25_hi) / 2, width = q25_hi - q25_lo),
                height = 0.6, alpha = 0.3, colour = NA) +
      geom_tile(aes(x = (q75_lo + q75_hi) / 2, width = q75_hi - q75_lo),
                height = 0.6, alpha = 0.3, colour = NA) +
      ## the whisker itself: posterior median interquartile range, with end caps
      geom_errorbar(aes(xmin = iqr_lo, xmax = iqr_hi), orientation = "y",
                    width = 0.35, linewidth = 0.5) +
      geom_point(aes(x = mid), size = 1.6) +
      scale_colour_manual(values = cols, guide = "none") +
      scale_fill_manual(values   = cols, guide = "none") +
      facet_wrap(~param, scales = "free_x", ncol = 3) +
      labs(x = NULL, y = NULL) +
      theme_bw(base_size = 9) +
      theme(panel.grid = element_blank(), panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            strip.background = element_blank(),
            strip.text = element_text(hjust = 0, face = "bold", size = 9),
            axis.text = element_text(colour = "black"))
  }

  top <- df[df$param %in% BAND_PANELS, ]
  bot <- df[!(df$param %in% BAND_PANELS) & df$strain != "Prior", ]

  ## Heights in proportion to the rows each block draws -- one facet row of five categories
  ## against two facet rows of four -- so every row comes out the same height.
  p <- if (nrow(bot)) ggpubr::ggarrange(build(top), build(bot), nrow = 2,
                                        heights = c(1, 1.75)) else build(top)

  if (!is.null(file)) ggsave(file, p, width = 8, height = 7.2)
  p
}

## ---------------------------------------------------------------------------
## figure_prior_posterior -- prior against posterior, every fitted parameter
##
##   figure_prior_posterior(file = "images/figureS_prior_posterior.png")
##   figure_prior_posterior(pars = c(KIN_PARS, paste0("test", 0:5)))
##
## The fuller answer to R1.4a. Figure 4's prior row shows what the prior implies about the
## POPULATION; this shows the prior on the parameters themselves, which is where a reader
## checks whether an estimate is data-driven.
##
## Prior and posterior are both drawn from samples and passed through the same density
## estimator, so any difference between the curves is a difference between the
## distributions rather than between two ways of drawing them.
##
## The multiplier is absent by design: prior() covers parameters 9:(length-1), which stops
## one short of it, so the multiplier has no proper prior and there is nothing to compare
## a posterior against. Worth stating in the caption rather than leaving a reader to
## wonder why the incidence scaling is missing.
## ---------------------------------------------------------------------------
KIN_PARS <- c("a_bar", "a_sigma", "b_bar", "b_sigma", "vl_max_bar", "vl_max_sigma")
PP_PARS  <- c(KIN_PARS, "inc1", "inc2", "inc3", "rdisp")

## ---------------------------------------------------------------------------
## figure5 -- inferred daily incidence of symptom onset, by variant
##
##   figure5(file = "images/figure5.png")
##   f5 <- figure5(); f5$data          # the underlying series, for a supplementary table
##
## Equation 2: I_symp(t) = N_pop * exp(B(theta_s, t)) * M, evaluated per posterior draw by
## infecteds_generator() -- the same function the likelihood uses, so the curve plotted is
## the curve fitted rather than a reconstruction of it.
##
## Each variant is drawn over its own observation window. Day i of a variant's spline maps to
## x$dates[i] one-for-one, since max_day is set to n_days and dates are the data array's
## calendar column names.
##
## NOTE ON LABELLING (R1.3 / action F). This is incidence of SYMPTOM ONSET, not of infection.
## The figure caption in the submitted manuscript says so, but the Results text at L199 calls
## it "the incidence of infections giving rise to pillar 2 cases". One of the two needs to
## change, and the figure is the one that is right.
##
## Both ends are trimmed by default, for different reasons.
##
## LEAD-IN: the first 35 days. p_array_func only starts filling at start_j = 36, so no data
## constrains them and the curve there is prior plus spline smoothness.
##
## TAIL: two separate defects, so the cut is the earlier of the two. NO knots are padded
## after the data -- every variant's array ends on its last day with counts -- so neither of
## these is a designed lead-out.
##
##   (a) EXTRAPOLATION. splinefun is built on control points at c(0, knots * max_day), so
##       days after the final knot are extrapolated rather than interpolated and the cubic
##       runs off in whatever direction the last two control points imply. Untrimmed,
##       omicron falls to ~8,000 in mid-January then climbs back above 25,000 by the 23rd.
##
##       How far this reaches is INCIDENTAL. wt (217 days) and alpha (329) both happen to
##       leave a 7-day remainder over 14-day knots; omicron's knots are built as
##       seq(0, n_days - 7, 7) and so stop a week short by construction; delta (308 = 22 x 14)
##       leaves nothing at all. "Days past the last knot" is therefore not a rule, which is
##       why (b) has to be applied as well.
##
##   (b) PARTIAL OBSERVATION. p_array_func fills observation day j from
##       infecteds[j:(j - max_test_day)], so incidence on day s is seen by observations on
##       days s ... s+6 -- and only those that exist. The last day's incidence is used once,
##       by that day's day-0 tests; the day before it twice; and so on. The final
##       max_test_day days are therefore progressively less informed even where the spline
##       is still interpolating, which is exactly delta's situation: its last knot lands on
##       its last day, so (a) never fires and (b) is the only thing wrong with its tail.
##
## Both cuts are derived from the stored x -- knots and n_test_days -- rather than
## hard-coded, so they adapt to each variant's own knot spacing and testing window.
## ---------------------------------------------------------------------------
figure5 <- function(strains = STRAINS, n_samples = 300, file = NULL,
                    dir = STRAIN_DIR, trim_leadin = TRUE, trim_tail = TRUE,
                    log_y = FALSE, width = 9, height = 4,
                    quantity = c("onset", "fitted")) {

  ## Same lazy load as ppc_draws, so this file still works without the C++ kernel compiled.
  ## infecteds_generator is pure R -- it evaluates the spline and scales by the population --
  ## so nothing is built here.
  if (!exists("infecteds_generator")) source("scripts/model.R")
  quantity <- match.arg(quantity)

  ## WHICH CURVE. Under thresh_tdist the spline IS incidence of symptom onset, so both
  ## choices give the same curve. Under the incidence form the spline is incidence of CROSSING
  ## vl_thresh, and onset follows it by the fitted LogNormal(inc1, inc2) delay. "onset" (the
  ## default) convolves the two -- O(s) = sum_tau I(s - tau) f(tau), with f discretised and
  ## truncated at symp_delay_lim exactly as p_array_func does -- so Figure 5 plots the same
  ## quantity under either model and the two versions can be compared directly. "fitted"
  ## plots the spline itself.
  form_used <- NULL

  series <- do.call(rbind, lapply(strains, function(s) {
    message(s)
    sd <- sample_draws(s, n_samples, dir)
    x  <- sd$x
    fm <- fit_form(x); form_used <<- fm
    L  <- if (fm == "incidence") x$symp_delay_lim else 0L
    convolve_onset <- fm == "incidence" && quantity == "onset"

    inc <- t(vapply(seq_len(nrow(sd$draws)), function(j) {
      p <- sd$draws[j, ]
      I <- as.numeric(infecteds_generator(p, x$knots, x$population, x$max_day, form = fm))
      if (!convolve_onset) return(I)
      f <- stats::plnorm(0:L + 0.5, p[["inc1"]], p[["inc2"]]) -
           stats::plnorm(pmax(0:L - 0.5, 0), p[["inc1"]], p[["inc2"]])
      f <- f / sum(f)
      O <- rep(NA_real_, x$max_day)
      for (d in (L + 1L):x$max_day) O[d] <- sum(I[d - 0:L] * f)
      O
    }, numeric(x$max_day)))

    ## Lead-in: p_array_func fills from start_j = 36. Onset on day s needs crossings back to
    ## s - L, which day 36 always has (L = 20).
    from <- if (trim_leadin) 36L else if (convolve_onset) L + 1L else 1L

    ## (a) last day the spline interpolates rather than extrapolates. floor() because
    ##     knots * max_day need not land on an integer day. For the convolved onset curve the
    ##     extrapolated days feed only its final day, so the same cut applies.
    last_interp <- floor(max(x$knots) * x$max_day)
    ## (b) last day whose incidence is seen by a full set of test days 0..max_test_day. A
    ##     crossing is observed only after its onset delay as well, so the fitted crossing
    ##     curve under the incidence form loses a further L days.
    last_full   <- x$max_day - (x$n_test_days - 1L) -
                   if (fm == "incidence" && !convolve_onset) L else 0L
    to <- if (trim_tail) min(last_interp, last_full) else x$max_day

    if (to <= from) stop(s, ": nothing left after trimming (", from, " to ", to, ")")
    if (trim_tail && to < x$max_day)
      message(sprintf("  %s: kept to day %d of %d (dropped %d) -- extrapolation past day %d, partial observation past day %d",
                      s, to, x$max_day, x$max_day - to, last_interp, last_full))
    keep <- from:to
    q <- apply(inc[, keep, drop = FALSE], 2, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)

    data.frame(strain = STRAIN_LABEL[[s]],
               date   = as.Date(x$dates[keep]),
               lo = q[1, ], mid = q[2, ], hi = q[3, ], row.names = NULL)
  }))
  series$strain <- factor(series$strain, levels = unname(STRAIN_LABEL[strains]))
  cols <- setNames(unname(STRAIN_COL[strains]), unname(STRAIN_LABEL[strains]))

  p <- ggplot(series, aes(date, mid, colour = strain, fill = strain)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.25) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_fill_manual(values = cols, name = NULL) +
    scale_x_date(date_labels = "%b-'%y", date_breaks = "3 months") +
    labs(x = NULL, y = if (identical(form_used, "incidence") && quantity == "fitted")
                         expression("Number of individuals crossing "*vl[thresh]*", "*I(t))
                       else expression("Number of individuals developing symptoms, "*I^symp*(t))) +
    theme_bw(base_size = 11) +
    ## No gridlines, matching every other figure here. theme_bw draws both major and minor by
    ## default, and this used to blank only the minor ones.
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black"),
          legend.position = "right")

  p <- p + if (log_y) scale_y_log10(labels = scales::comma)
           else scale_y_continuous(labels = scales::comma,
                                   expand = expansion(mult = c(0, 0.05)))

  if (!is.null(file)) ggsave(file, p, width = width, height = height, dpi = 300)
  print(p)
  invisible(list(plot = p, data = series))
}

## rdisp is drawn with its prior truncated at zero, via `prior_bounds`. Its stored lower_bound
## is -Inf, so the prior as specified puts about 46% of its mass on negative values; the
## likelihood excludes those, because a negative binomial needs a positive size, so the region
## is not one the posterior can occupy and drawing it only stretched the panel to span -18 to
## +20. Truncating shows the prior over the support the data actually see. It is a display
## decision about which half of a prior to draw, not a change to the fit -- see prior_draws().
##
## Ten parameters also fill the 5 x 2 grid exactly, so no panel is left carrying its own x axis.
figure_prior_posterior <- function(strains = STRAINS,
                                   pars = PP_PARS, n_prior = 20000,
                                   thin = 20L, file = NULL, dir = STRAIN_DIR,
                                   ncol = 5, prior_bounds = c(rdisp = 0)) {

  smp <- lapply(strains, function(s)
    do.call(rbind, lapply(chain_tags(s, dir),
                          function(t) combine_chains(s, t, thin = thin, dir = dir))))
  names(smp) <- strains
  ## Only bound what is actually being drawn, so a caller narrowing `pars` does not trip the
  ## "no such parameter" check on a bound for a panel they left out.
  pb  <- prior_bounds[intersect(names(prior_bounds), pars)]
  pr  <- prior_draws(strains[1], n_prior, dir, pars = pars, bounds = pb)  # same for all variants
  lev <- intersect(pars, colnames(smp[[1]]))

  ## Densities are evaluated on an explicit grid rather than left to geom_density, because
  ## the display window has to be chosen and geom_density would set it from whichever
  ## sample is widest. For inc3 the prior is 364 times wider than the posterior, so the
  ## posterior becomes an invisible spike at the origin of a panel scaled to the prior.
  ##
  ## density(from, to) still estimates from the whole sample and only evaluates on that
  ## grid, so heights stay comparable between panels -- clipping the SAMPLES instead would
  ## renormalise each curve over its window and make a truncated prior look taller than it
  ## is.
  ##
  ## The window spans the posteriors plus the prior's central 90%. Including the prior's
  ## tails would restore the spike; excluding the prior entirely would lose the comparison
  ## for a_bar and vl_max_bar, where prior and posterior barely overlap.
  ## Truncation points, for the boundary correction below.
  LB <- attr(pr, "lower_bound")

  curves <- do.call(rbind, lapply(lev, function(p) {
    pv  <- pr[, p]
    pos <- lapply(smp, function(d) d[, p])
    lo  <- min(unlist(lapply(pos, min)), quantile(pv, 0.05))
    hi  <- max(unlist(lapply(pos, max)), quantile(pv, 0.95))
    pad <- 0.04 * (hi - lo); lo <- lo - pad; hi <- hi + pad
    b   <- if (!is.null(LB) && p %in% names(LB)) LB[[p]] else -Inf

    ## Reflection at a finite lower bound. A Gaussian kernel has support on the whole line, so
    ## a sample bounded at b loses the mass the kernels place below b and the estimate sags
    ## toward the boundary. Mirroring the sample about b and doubling the result on the valid
    ## side restores it, which is exact for a density whose derivative vanishes at b and a good
    ## correction otherwise.
    ##
    ## Without it inc3's prior -- a half-normal, mean and bound both zero -- is drawn with its
    ## mode near 0.4 and falling to the boundary, the opposite shape from the truth. The
    ## bandwidth is taken from the ORIGINAL sample and passed in, because the mirrored sample is
    ## symmetric about b and would otherwise widen bw.nrd0's spread estimate and oversmooth.
    ##
    ## For a sample sitting many bandwidths clear of its bound the correction is numerically a
    ## no-op, so it is applied to prior and posterior alike rather than guessed at per curve.
    dens <- function(v, g) {
      if (is.finite(b)) {
        bw <- stats::bw.nrd0(v)
        d  <- stats::density(c(v, 2 * b - v), bw = bw, from = lo, to = hi, n = 512)
        d$y <- ifelse(d$x >= b, 2 * d$y, 0)
      } else {
        d <- stats::density(v, from = lo, to = hi, n = 512)
      }
      data.frame(parameter = p, group = g, x = d$x, y = d$y, stringsAsFactors = FALSE)
    }
    rbind(dens(pv, "Prior"),
          do.call(rbind, lapply(strains, function(s) dens(pos[[s]], STRAIN_LABEL[[s]]))))
  }))

  ## Every curve scaled to peak at 1. Without it the comparison is unreadable wherever the
  ## data is much more informative than the prior: inc3's posterior is 364 times narrower,
  ## so on a shared vertical scale the prior is a line indistinguishable from the axis and
  ## looks like a plotting failure rather than a result. Height carries no information here
  ## anyway -- the y axis is unlabelled and what matters is where each distribution sits
  ## and how wide it is.
  curves$y <- ave(curves$y, curves$parameter, curves$group,
                  FUN = function(v) if (max(v) > 0) v / max(v) else v)

  curves$parameter <- factor(curves$parameter, levels = lev)
  post <- curves[curves$group != "Prior", ]
  pri  <- curves[curves$group == "Prior", ]
  post$strain <- factor(post$group, levels = unname(STRAIN_LABEL[strains]))
  cols <- setNames(unname(STRAIN_COL[strains]), unname(STRAIN_LABEL[strains]))

  p <- ggplot(post, aes(x = x, y = y, colour = strain)) +
    ## Prior first and filled, so it reads as the backdrop a posterior is judged against.
    ## Where it is flat across the panel, the prior is uninformative over the range the
    ## posterior occupies -- which is the point.
    geom_area(data = pri, aes(x = x, y = y), inherit.aes = FALSE,
              fill = "grey70", colour = "grey35", alpha = 0.5, linewidth = 0.3) +
    geom_line(linewidth = 0.55) +
    facet_wrap(~parameter, scales = "free", ncol = ncol, labeller = rec_labeller()) +
    scale_colour_manual(values = cols, name = NULL) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold", size = 9),
          axis.text = element_text(colour = "black"),
          legend.position = "top", legend.margin = margin(b = 0))

  if (!is.null(file))
    ## Width scales with the column count so a panel keeps about two inches whatever the
    ## layout. Fixed at 8 it gave 2.0 in at the old four columns and 1.6 in at five, which is
    ## too narrow for a density curve to carry its own x axis labels.
    ggsave(file, p, width = 2 * ncol, height = 1.9 * ceiling(length(lev) / ncol) + 0.6)
  p
}

## ---------------------------------------------------------------------------
## RECOVERY -- the simulation study: posterior against prior against truth
##
##   figure_recovery("M")                                  # one family
##   figure_recovery("M", file = "images/recovery_M.png")
##   r <- figure_recovery("N"); r$table                     # coverage, parameter by rung
##   recovery_table(c("M0","M1","M2","M3"))                 # the numbers, no plot
##
## This is E1 and R1.4c. Each panel is one parameter; each row is one rung of the incidence
## ladder, labelled by the number of positive tests it actually contains, so the narrowing of
## the intervals reads directly against sample size. The red dashed line is the generating
## value and the blue row is the prior, which is what a posterior has to beat to have learned
## anything.
##
## Ported from prior_posterior_plot() in simulation_study.R, with three changes: the window
## is the last N iterations rather than a fraction of the chain, coverage is computed and
## coloured rather than left to the eye, and the generating values are cross-checked against
## data/synthetic_truth.rds instead of trusted from the chain alone.
## ---------------------------------------------------------------------------
REC_PARS   <- c("a_bar", "a_sigma", "b_bar", "b_sigma", "vl_max_bar", "vl_max_sigma",
                "inc1", "inc2", "inc3", "rdisp")
TRUTH_FILE <- "data/synthetic_truth.rds"

## Manuscript symbols for the recovery figures, as plotmath. Code names (a_bar, inc1, rdisp)
## are internal and mean nothing to a reader; these are the names used in the Methods and in
## Table 1. Follows NOTATION_OMEGA so the scale parameters and the onset pair change with the
## rest of the paper rather than drifting out of step -- see NOTATION.md.
REC_SYM <- c(
  a_bar        = "eta[a]",
  a_sigma      = if (NOTATION_OMEGA) "omega[a]" else "tau[a]",
  b_bar        = "eta[b]",
  b_sigma      = if (NOTATION_OMEGA) "omega[b]" else "tau[b]",
  vl_max_bar   = "bar(vl)[peak]",
  vl_max_sigma = "vl[sigma]",
  inc1         = if (NOTATION_OMEGA) "eta[T]"   else "phi[1]",
  inc2         = if (NOTATION_OMEGA) "omega[T]" else "phi[2]",
  inc3         = "vl[thresh]",
  rdisp        = "xi")

## Panel letters shared by figure_scenarios() and figure_recovery_grid(), so a reader can
## carry a column of the recovery grid back to the trajectories that produced it. Defined once
## because the two figures must not drift: both order their panels less-variable profile first
## and then by onset regime, which is the order below.
## The A-F generation is named so that the family IS its panel, which is why its entries are
## the identity. M-R and S-X keep the mapping they were built with.
SCEN_PANEL <- c(P = "A", Q = "B", R = "C", M = "D", N = "E", O = "F",
                V = "A", W = "B", X = "C", S = "D", T = "E", U = "F",
                A = "A", B = "B", C = "C", D = "D", E = "E", F = "F")

## Parsed labels for the parameter facet, plain text for everything else.
##
## Anything absent from REC_SYM falls back to its own name, NOT to plotmath. The previous
## as_labeller(REC_SYM, label_parsed) passed unmatched names straight to parse(), which is fine
## for a valid R symbol like "multiplier" or "theta1" but a hard error for any label carrying a
## space -- check_chains() facets on "log posterior", and recovery_table() documents
## pars = c(REC_PARS, "thetas"). The fallback costs nothing and makes the labeller safe to
## attach to every parameter facet in the file rather than only the ones known to be REC_PARS.
## The mapping has to go through as_labeller() with label_parsed as its default, exactly as
## before. Handing labeller() a bare function instead does NOT work: it re-wraps anything
## lacking the "labeller" class with as_labeller(x, label_value), so a returned call prints as
## "eta[a]" and a returned expression as "expression(eta[a])".
##
## The only change from the original as_labeller(REC_SYM, label_parsed) is the fallback. A
## named character vector passes an unmatched label straight through to parse(), which is fine
## for a valid R symbol but a hard error on anything containing a space -- check_chains() facets
## on "log posterior", and recovery_table() documents pars = c(REC_PARS, "thetas").
## Backtick-quoting makes any name parse to a symbol carrying its own spelling, so it renders
## unchanged and the labeller is safe to attach to every parameter facet in the file.
rec_labeller <- function() {
  src <- function(s)
    ifelse(s %in% names(REC_SYM), unname(REC_SYM[s]), paste0("`", s, "`"))
  ggplot2::labeller(parameter = ggplot2::as_labeller(src, ggplot2::label_parsed),
                    .default = ggplot2::label_value)
}

## facet_nested spans an outer strip across the columns it groups, which plain facet_grid
## repeats once per column. Fall back cleanly if ggh4x is absent.
facet_span <- function(f, ...) {
  if (requireNamespace("ggh4x", quietly = TRUE)) ggh4x::facet_nested(f, ...)
  else ggplot2::facet_grid(f, ...)
}

## The last n_last filled draws of every chain for a scenario, pooled.
##
## A trailing window rather than everything-after-burn-in, because the kernel re-adapts
## at every resume and once more mid-continuation, and draws either side of an adaptation
## come from different proposals. Taking the tail keeps the sample inside one kernel. The
## warning fires when it does not.
##
## The default dir holds the JOINED chains (join_chains.R), one file per scenario under
## tag "c", so chain_tags() finds exactly one chain and pooling is a non-issue. The raw
## segments (_1 -> _2 -> _3, each resuming from the previous one's final state) are in
## outputs/old; if reading those, pass a single tag -- pooling segments would mix pre-
## and post-resume draws and count the same chain several times.
recovery_draws <- function(scenario, n_last = 10000L, dir = "outputs/final/simulations",
                           tags = NULL) {
  if (is.null(tags)) tags <- chain_tags(scenario, dir)
  dr <- NULL; x <- NULL
  for (tg in tags) {
    f <- chain_file(scenario, tg, dir)
    ## chain_file returns NULL when nothing matches, and readRDS(NULL) reports only "bad
    ## 'file' argument", which says nothing about which scenario or directory was wrong.
    ## The usual cause is the default dir: the ladder scenarios live in
    ## outputs/final/simulations, the variant-sim fits in outputs.
    if (is.null(f))
      stop("no chain '", tg, "' for ", scenario, " in ", dir,
           " -- check the `dir` argument", call. = FALSE)
    o  <- readRDS(f)
    fl <- max(which(!is.na(o$MCMC_posteriors)))
    ## A curated chain records where its analysis window opens, derived from its own
    ## cov_start: the draw after the last covariance adaptation that actually happened. Prefer
    ## it to n_last, which counts back from the end and cannot know where the kernel was last
    ## re-estimated -- for a joined chain it can reach back across the junction into a segment
    ## that ran a different proposal, or into pre-convergence draws.
    ##
    ## n_last = NULL means "work it out", and sim_window() is consulted only here, AFTER
    ## analysis_from has been ruled out. Calling it up front -- which recovery_table used to do
    ## -- makes a curated directory emit one "no segments found in outputs/old" warning per
    ## scenario for a number that is then thrown away.
    st <- if (!is.null(o$x$analysis_from)) as.integer(o$x$analysis_from)
          else max(1L, fl - as.integer(if (is.null(n_last)) sim_window(scenario)
                                       else n_last) + 1L)

    ## Guard against a window that reaches back across a covariance change -- but ONLY when
    ## the chain has not recorded its own window start.
    ##
    ## x$cov_start is the schedule input_generator was built with, in the numbering of the
    ## chain that first ran it. A joined chain inherits it unchanged, so for B3/D3/E3/F3 --
    ## a 95,001-iteration first segment plus a continuation -- it still reads
    ## 80,000/100,000/120,000/140,000 even though the 100k, 120k and 140k adaptations could
    ## never have fired: the segment ended first. Comparing analysis_from against it then
    ## reports a re-adaptation that never happened. Verified from the draws: the median
    ## accepted step and the acceptance rate are flat across 140,000 in all four, and change
    ## exactly once, at the junction that analysis_from already marks.
    ##
    ## analysis_from is written by the curation step, which knows where the junction is.
    ## Where it exists it is authoritative and cov_start is stale, so only the fallback path
    ## -- a window counted back from the end, which genuinely cannot know -- is checked.
    if (is.null(o$x$analysis_from)) {
      ad <- suppressWarnings(max(o$x$cov_start[o$x$cov_start <= fl]))
      if (is.finite(ad) && st < ad)
        warning(sprintf("%s chain %s: window starts at %s but the proposal was re-adapted at %s",
                        scenario, tg, format(st, big.mark = ","), format(ad, big.mark = ",")),
                call. = FALSE)
    }

    dr <- rbind(dr, on_natural_scale(o$MCMC_output[st:fl, , drop = FALSE], o$x))
    if (is.null(x)) x <- o$x
    rm(o); invisible(gc(verbose = FALSE))
  }
  list(draws = dr, x = x, n_chains = length(tags))
}

## Generating values for a scenario. The chain carries them (input_generator attaches
## x$truth), but so does data/synthetic_truth.rds, and the two are independent records of the
## same thing -- one written when the chain was built, one extracted from the definitions in
## simulation_study.R. Disagreement means the scenario was redefined after the array was
## generated, which would invalidate every recovery claim silently. So check.
scenario_truth <- function(x, scenario, pars) {
  tr <- x$truth
  if (is.null(tr) || all(is.na(tr[pars])))
    stop("no generating values in the chain for ", scenario, " -- fitted to real data?")
  if (file.exists(TRUTH_FILE)) {
    f <- readRDS(TRUTH_FILE)[[scenario]]
    if (!is.null(f)) {
      p <- intersect(pars, names(f))
      d <- abs(as.numeric(tr[p]) - as.numeric(f[p]))
      if (any(d > 1e-8, na.rm = TRUE))
        warning(scenario, ": chain and ", TRUTH_FILE, " disagree on ",
                paste(p[d > 1e-8], collapse = ", "), call. = FALSE)
    }
  }
  setNames(as.numeric(tr[pars]), pars)
}

## Truth, posterior quantiles and coverage, one row per scenario x parameter. This is the
## table the paper reports: with one fit per scenario, recovery is assessed marginally --
## does each 95% interval contain the known value -- and NOT as a coverage rate, which would
## need replicates.
## `pars` accepts the token "thetas", which expands PER SCENARIO to whatever spline
## coefficients that scenario has. They cannot be listed literally: the number of knots
## differs by variant -- 16 thetas for wtsim, 24 for alphasim, 23 for deltasim, 13 for
## omicronsim -- so a single fixed vector cannot serve them all.
##
##   recovery_table(SC, pars = c(REC_PARS, "thetas", "multiplier"), ...)
##
## Rows for the pinned spline coefficient are flagged rather than dropped. It is held at its
## generating value for the whole run and never proposed, so its interval is degenerate and
## its "recovery" is circular -- but it is worth seeing that it sat where it was meant to.
recovery_table <- function(scenarios, pars = REC_PARS, n_last = NULL,
                           dir = "outputs/final/simulations", tags = NULL) {
  do.call(rbind, lapply(scenarios, function(s) {
    ## n_last passes through untouched, NULL included: recovery_draws decides the window, and
    ## only it knows whether the chain records an analysis_from that makes sim_window()
    ## irrelevant. n_draws is then read back off the draws actually returned rather than from
    ## the request, so the column cannot disagree with what was summarised.
    d  <- recovery_draws(s, n_last, dir, tags)
    nl <- nrow(d$draws) %/% d$n_chains

    nm <- colnames(d$draws)
    pp <- unlist(lapply(pars, function(p)
      if (identical(p, "thetas")) grep("^theta[0-9]+$", nm, value = TRUE) else p))
    miss <- setdiff(pp, nm)
    if (length(miss))
      stop(s, ": no such parameter(s) -- ", paste(miss, collapse = ", "), call. = FALSE)

    pin <- names(d$x$parameters)[d$x$ignored_spline_par]
    tr <- scenario_truth(d$x, s, pp)
    q  <- t(apply(d$draws[, pp, drop = FALSE], 2, quantile, c(0.025, 0.5, 0.975)))
    sd_p <- apply(d$draws[, pp, drop = FALSE], 2, sd)
    data.frame(scenario = s,
               n_draws  = nl,
               n_tests  = sum(d$x$data_array),
               parameter = pp, truth = tr[pp],
               lo = q[, 1], mid = q[, 2], hi = q[, 3],
               covered = tr[pp] >= q[, 1] & tr[pp] <= q[, 3],
               ## held at truth and never proposed, so coverage here means nothing
               pinned = pp %in% pin,
               ## how far off in posterior SDs -- separates a near miss from a real failure.
               ## NA rather than Inf for the pinned coefficient, whose SD is zero.
               z     = ifelse(sd_p > 0, (q[, 2] - tr[pp]) / sd_p, NA_real_),
               width = q[, 3] - q[, 1],
               row.names = NULL)
  }))
}

## `scenarios` overrides family/rungs, for targets that are not a sample-size ladder --
## the <variant>sim fits in particular, whose rows are four different variants rather than
## four sample sizes. `row_label` follows: "n_tests" labels each row by its number of
## positive tests, which is the point of a ladder; "scenario" uses the name, which is the
## point when the rows are not nested.
figure_recovery <- function(family = "M", rungs = 0:3, scenarios = NULL, pars = REC_PARS,
                            n_last = 10000L, max_prior_ratio = 25, ncol = 5,
                            row_label = c("n_tests", "scenario"),
                            dir = "outputs/final/simulations", tags = NULL, file = NULL) {

  row_label <- match.arg(row_label)
  scen <- if (is.null(scenarios)) paste0(family, rungs) else scenarios
  post <- recovery_table(scen, pars, n_last, dir, tags)

  ## Take the parameter set from what recovery_table actually returned rather than from
  ## `pars`, so the "thetas" token works here too. Scenarios can carry different numbers of
  ## spline coefficients, in which case this is their union and a scenario missing one simply
  ## has no bar in that panel.
  pars <- unique(as.character(post$parameter))

  ## Marginal priors. The six kinetics parameters have a truncated MVN prior, but each
  ## marginal is a truncated normal with the same mean and SD, and only the marginal is
  ## being drawn here. The richest rung is the reference because input_generator derives
  ## the test0-test5 priors from the data array, so they are best determined there.
  ref <- recovery_draws(scen[length(scen)], 1L, dir, tags)$x
  ## Parameters absent from this scenario, and the multiplier, which prior() never evaluates
  ## (U8), give NA and simply have no prior bar.
  pm  <- ref$prior_mean[pars]; ps <- ref$prior_sds[pars]
  lb  <- ref$lower_bound[pars]; lb[is.na(lb)] <- -Inf
  ok  <- !is.na(pm) & !is.na(ps) & ps > 0
  prior <- data.frame(
    scenario = "Prior", parameter = pars,
    lo  = ifelse(ok, truncnorm::qtruncnorm(0.025, a = lb, mean = pm, sd = ps), NA_real_),
    mid = ifelse(ok, pm, NA_real_),
    hi  = ifelse(ok, truncnorm::qtruncnorm(0.975, a = lb, mean = pm, sd = ps), NA_real_),
    covered = NA, row.names = NULL)

  ## A prior far wider than every posterior flattens the panel onto a point. Drop it there
  ## and label the panel, rather than losing the posterior comparison to keep a bar that
  ## runs off the axis anyway.
  pw   <- tapply(post$hi - post$lo, post$parameter, max)
  wide <- prior$parameter[(prior$hi - prior$lo) / pmax(pw[prior$parameter], 1e-12) > max_prior_ratio]
  prior[prior$parameter %in% wide, c("lo", "mid", "hi")] <- NA

  lab <- if (row_label == "n_tests")
    vapply(split(post$n_tests, post$scenario)[scen],
           function(v) format(v[1], big.mark = ","), character(1))
  else setNames(scen, scen)

  df <- rbind(prior, post[, names(prior)])
  df$row       <- ifelse(df$scenario == "Prior", "Prior", lab[df$scenario])
  df$row       <- factor(df$row, levels = rev(c("Prior", unname(lab))))
  df$parameter <- factor(df$parameter, levels = pars)
  ## grey covers the truth, red misses it, blue is the prior
  df$kind <- ifelse(df$scenario == "Prior", "Prior",
                    ifelse(df$covered, "Covers truth", "Misses truth"))

  ## Kinetics are held fixed across the rungs of a family, so one truth line per panel is
  ## right for them. The THETAS are not: the ladder adds log(10) to every spline coefficient
  ## per rung, so their generating value differs on every row, and a single line would draw
  ## rung 0's value across all four -- three of them wrong by up to 6.9 log units. Constant
  ## parameters keep the line; varying ones get a tick on each row instead.
  ##
  ## The `covered` column was never affected: recovery_table compares each scenario against
  ## its own x$truth. This is a display fix only.
  spl   <- split(post$truth, post$parameter)[pars]
  const <- vapply(spl, function(v) diff(range(v)) < 1e-8, logical(1))

  tru <- data.frame(parameter = factor(pars[const], levels = pars),
                    truth = vapply(spl[pars[const]], `[`, numeric(1), 1), row.names = NULL)
  tru_row <- post[!post$parameter %in% pars[const], c("scenario", "parameter", "truth")]
  if (nrow(tru_row)) {
    tru_row$row       <- factor(lab[tru_row$scenario], levels = levels(df$row))
    tru_row$parameter <- factor(tru_row$parameter, levels = pars)
  }

  note <- if (length(wide))
    data.frame(parameter = factor(wide, levels = pars),
               row = factor("Prior", levels(df$row)),
               label = "prior wider than axis", row.names = NULL)

  n_miss <- sum(!post$covered)
  p <- ggplot(df, aes(y = row)) +
    geom_vline(data = tru, aes(xintercept = truth), colour = "firebrick",
               linetype = 2, linewidth = 0.4) +
    geom_linerange(aes(xmin = lo, xmax = hi, colour = kind), linewidth = 1.1, na.rm = TRUE) +
    geom_point(aes(x = mid, colour = kind), size = 1.7, na.rm = TRUE) +
    ## drawn after the intervals so a tick sitting inside one is still visible
    { if (nrow(tru_row)) geom_point(data = tru_row, aes(x = truth, y = row),
                                    inherit.aes = FALSE, shape = 124, size = 3.2,
                                    colour = "firebrick") } +
    { if (!is.null(note)) geom_text(data = note, aes(x = -Inf, label = label),
                                    hjust = -0.05, size = 2.5, colour = "steelblue") } +
    scale_colour_manual(values = c(`Covers truth` = "grey15",
                                   `Misses truth` = "firebrick",
                                   Prior          = "steelblue"), name = NULL) +
    facet_wrap(~parameter, scales = "free_x", ncol = ncol, labeller = rec_labeller()) +
    labs(x = NULL, y = if (row_label == "n_tests") "positive tests" else NULL,
         title = if (is.null(scenarios))
                   sprintf("Family %s: posterior against prior and generating value", family)
                 else "Posterior against prior and generating value",
         subtitle = sprintf("last %s iterations  |  red dashed = truth  |  %d of %d intervals miss",
                            format(n_last, big.mark = ","), n_miss, nrow(post))) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold"),
          axis.text = element_text(colour = "black"),
          axis.text.y = element_text(size = 7),
          legend.position = "top")

  if (!is.null(file))
    ggsave(file, p, width = 10, height = 1.9 * ceiling(length(pars) / ncol) + 1.2)
  print(p)
  invisible(list(plot = p, table = post))
}

## ---------------------------------------------------------------------------
## figure_recovery_grid -- every scenario in one figure: families down, parameters across
##
##   figure_recovery_grid()
##   figure_recovery_grid(file = "outputs/figures/figS_recovery_grid.png")
##
## figure_recovery() shows one family at a time, which makes it hard to see whether a
## parameter fails everywhere or only under one onset regime. This lays all six families out
## as rows against the same ten parameter columns, so a column that is red down its whole
## length is a parameter the method does not recover, while one red only in the P/Q/R half is
## a property of that kinetics profile.
##
## Row strips carry both groupings: kinetics profile (M/N/O more variable, P/Q/R less) and
## onset regime, since those are the two things the design varies apart from sample size.
##
## Rows within each cell are the four ladder rungs, labelled by their NOMINAL target rather
## than their realised count -- the realised counts differ between families (983 for M0,
## 1,007 for P0 and so on), and putting each on its own axis level would give thirty-odd
## levels that no longer line up across rows. figure_recovery() keeps the realised counts.
## layout = "tall" transposes it: parameters down, families across, intervals drawn
## vertically. Ten columns of parameters cannot be read on a portrait page -- each panel gets
## about 0.7 in -- whereas six columns of families on ten rows fits A4 portrait comfortably.
## The transpose is not cosmetic: facet_grid frees the x scale by column only, so with
## parameters as ROWS the value has to move to the y axis for each parameter to keep its own
## scale.
figure_recovery_grid <- function(families = c("M", "N", "O", "P", "Q", "R"), rungs = 0:3,
                                 pars = REC_PARS, n_last = 10000L, max_prior_ratio = 25,
                                 layout = c("wide", "tall"),
                                 prior_panel = c("inline", "side", "none"),
                                 symmetric_y = TRUE,
                                 dir = "outputs/final/simulations", tags = NULL,
                                 file = NULL) {
  layout <- match.arg(layout)
  prior_panel <- match.arg(prior_panel)
  ## "side" needs the patchwork composition built below for the tall layout only; in the wide
  ## layout the prior bar sits in the same free-x column as the posteriors and there is nowhere
  ## to put a separate panel. Fail loudly rather than silently dropping the prior.
  if (layout == "wide" && prior_panel == "side")
    stop("prior_panel = 'side' is only available for layout = 'tall'")

  rung_lab <- c("1,000", "10,000", "100,000", "1,000,000")[rungs + 1L]
  keep <- c("scenario", "parameter", "lo", "mid", "hi", "covered", "truth")

  post <- list(); prior <- list()
  for (f in families) {
    scen <- paste0(f, rungs)
    z <- recovery_table(scen, pars, n_last, dir, tags)
    z$family <- f
    z$row    <- rung_lab[match(z$scenario, scen)]
    post[[f]] <- z[, c(keep, "family", "row")]

    ## Marginal priors, taken from the richest rung as figure_recovery does. Recomputed per
    ## family rather than once: the test0-test5 priors are derived from the data array, so
    ## they are not guaranteed identical between families even though the kinetics priors are.
    ref <- recovery_draws(scen[length(scen)], 1L, dir, tags)$x
    pm <- ref$prior_mean[pars]; ps <- ref$prior_sds[pars]
    lb <- ref$lower_bound[pars]; lb[is.na(lb)] <- -Inf
    ok <- !is.na(pm) & !is.na(ps) & ps > 0
    prior[[f]] <- data.frame(
      scenario = "Prior", parameter = pars,
      lo   = ifelse(ok, truncnorm::qtruncnorm(0.025, a = lb, mean = pm, sd = ps), NA_real_),
      mid  = ifelse(ok, pm, NA_real_),
      hi   = ifelse(ok, truncnorm::qtruncnorm(0.975, a = lb, mean = pm, sd = ps), NA_real_),
      covered = NA, truth = NA_real_, family = f, row = "Prior", row.names = NULL)
  }
  post  <- do.call(rbind, post)
  prior <- do.call(rbind, prior)

  ## The wide-prior test is GLOBAL per parameter, not per family: facet_grid frees the x
  ## scale by column, so one runaway prior bar would flatten that parameter in every row.
  pw   <- tapply(post$hi - post$lo, post$parameter, max, na.rm = TRUE)
  wide <- unique(prior$parameter[(prior$hi - prior$lo) /
                                 pmax(pw[prior$parameter], 1e-12) > max_prior_ratio])
  prior[prior$parameter %in% wide, c("lo", "mid", "hi")] <- NA

  df <- rbind(prior, post)
  df$row       <- factor(df$row, levels = rev(c("Prior", rung_lab)))
  df$parameter <- factor(df$parameter, levels = pars)
  df$kind      <- ifelse(df$scenario == "Prior", "Prior",
                  ifelse(is.na(df$covered), "Prior",
                  ifelse(df$covered, "Covers truth", "Misses truth")))
  ## Columns are identified by (profile, onset regime), not by the scenario letter: M and P
  ## are both "mainly pre-peak" and differ only in kinetics profile, which the outer strip
  ## already says. The letters are internal names and carry nothing for a reader.
  ##
  ## Less-variable first, so the panels run in increasing order of between-individual
  ## variance from left to right. This re-levels LOCALLY rather than editing PROF_LEVELS,
  ## because PROF_LTY is keyed on that order and figure_precision draws its linetypes from it.
  ##
  ## The "1:"/"2:" prefixes are dropped here: with the order reversed they would read
  ## "2: less variable" then "1: more variable", which looks like a mistake. Restore them by
  ## deleting the sub() calls if the text refers to the profiles by number.
  prof_clean <- function(f) {
    strip <- function(v) sub("^[0-9]+:\\s*", "", v)
    factor(strip(unname(PROF_LABEL[f])), levels = strip(rev(PROF_LEVELS)))
  }
  ## Column identity carries the panel letter as well as the onset regime, so each column has
  ## a distinct label and can be matched to a panel of figure_scenarios(). Six levels rather
  ## than onset's three; the nesting under prof still holds because each letter belongs to
  ## exactly one profile.
  ## Restricted to the families actually drawn. SCEN_PANEL carries the incidence-indexed
  ## families (S-X) alongside the original ones and gives each pair the same letter -- V is
  ## P's twin -- so taking every entry would generate duplicate factor levels and fail.
  panel_lab <- function(f) {
    ord <- names(sort(SCEN_PANEL[unique(as.character(f))]))
    factor(sprintf("%s: %s", SCEN_PANEL[f], unname(ONSET_LABEL[f])),
           levels = sprintf("%s: %s", SCEN_PANEL[ord], unname(ONSET_LABEL[ord])))
  }
  df$prof  <- prof_clean(df$family)
  df$onset <- panel_lab(df$family)

  ## One truth line per panel. Every parameter in REC_PARS is held constant across the rungs
  ## of a family -- the ladder shifts only the spline coefficients, which are not in
  ## REC_PARS -- so unlike figure_recovery() no per-row ticks are needed. Truth DOES differ
  ## between families, which is the point of having them as separate rows.
  tru <- unique(post[!is.na(post$truth), c("family", "parameter", "truth")])
  tru$parameter <- factor(tru$parameter, levels = pars)
  tru$prof  <- prof_clean(tru$family)
  tru$onset <- panel_lab(tru$family)

  n_miss <- sum(!post$covered, na.rm = TRUE)
  n_tot  <- sum(!is.na(post$covered))
  cols <- c(`Covers truth` = "grey15", `Misses truth` = "firebrick", Prior = "steelblue")
  sub  <- sprintf("%s  |  red dashed = truth  |  %d of %d intervals miss",
                  if (is.null(n_last)) "each scenario's final single-kernel segment"
                  else sprintf("last %s iterations", format(n_last, big.mark = ",")),
                  n_miss, n_tot)
  ## The title has to follow what is actually drawn: prior_panel = "none" removes the prior,
  ## and a title promising one would be wrong.
  ttl  <- if (prior_panel == "none")
            "Simulation recovery: posterior against generating value"
          else "Simulation recovery: posterior against prior and generating value"

  ## Dropping the prior rows is layout-independent, so it happens before the branch. Only the
  ## "side" panel, which exists in the tall layout alone, needs pr1 afterwards.
  ##
  ## The prior is identical in all six columns -- the kinetics block comes from the same
  ## longitudinal fit and the inc/rdisp priors are fixed constants -- so drawing it six times
  ## spends a sixth of the space on one repeated bar. Worse, it is much wider than any
  ## posterior for several parameters, so with free scales it sets the row's range and squashes
  ## the posteriors: eta_T's row spans -1 to 3 because of the prior while every posterior sits
  ## between 1 and 2.2.
  pr1 <- NULL
  if (prior_panel != "inline") {
    pr1 <- prior[prior$family == families[1], ]
    others <- prior[prior$family != families[1], ]
    chk <- merge(pr1[, c("parameter", "lo", "hi")], others[, c("parameter", "lo", "hi")],
                 by = "parameter", suffixes = c(".1", ".o"))
    if (any(abs(chk$lo.1 - chk$lo.o) > 1e-8 | abs(chk$hi.1 - chk$hi.o) > 1e-8, na.rm = TRUE))
      warning("priors are not identical across families; the side panel shows ",
              families[1], "'s only", call. = FALSE)
    df <- df[df$scenario != "Prior", , drop = FALSE]
    df$row <- droplevels(df$row)
  }

  base <- function(sz) list(
    scale_colour_manual(values = cols, name = NULL),
    theme_bw(base_size = sz),
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_blank(),
          legend.position = "top",
          legend.margin = margin(0, 0, 0, 0),
          axis.text = element_text(colour = "black")))

  if (layout == "wide") {
    p <- ggplot(df, aes(y = row)) +
      geom_vline(data = tru, aes(xintercept = truth), colour = "firebrick",
                 linetype = 2, linewidth = 0.4) +
      geom_linerange(aes(xmin = lo, xmax = hi, colour = kind), linewidth = 0.9, na.rm = TRUE) +
      geom_point(aes(x = mid, colour = kind), size = 1.3, na.rm = TRUE) +
      facet_span(prof + onset ~ parameter, scales = "free_x", labeller = rec_labeller()) +
      labs(x = NULL, y = "positive tests (nominal)", title = ttl, subtitle = sub) +
      base(9) +
      theme(strip.text.x = element_text(hjust = 0, face = "bold"),
            strip.text.y = element_text(angle = 0, hjust = 0, size = 7),
            axis.text = element_text(colour = "black", size = 6))
    w <- 2.0 * length(pars) + 2.2; h <- 1.35 * length(families) + 1.4

  } else {
    ## prior_panel = "side" moves the prior to a single narrow panel on the right, letting each
    ## row rescale to its posteriors. The cost is that the side panel then carries its OWN y
    ## scale, so the prior is no longer directly comparable by eye -- it has to be read off
    ## its own axis. That trade must be stated in the caption.
    ## prior_panel = "none" drops it entirely; "inline" is the original behaviour.

    ## Short x labels: five ticks share about 1.1 in per column on A4 portrait, so the full
    ## "1,000,000" cannot be drawn even rotated.
    lev <- levels(df$row)
    short <- setNames(sub("^1(,000)+$", "", lev), lev)
    short[lev == "Prior"]     <- "prior"
    short[lev == "1,000"]     <- "1e3"
    short[lev == "10,000"]    <- "1e4"
    short[lev == "100,000"]   <- "1e5"
    short[lev == "1,000,000"] <- "1e6"

    p <- ggplot(df, aes(x = row)) +
      geom_hline(data = tru, aes(yintercept = truth), colour = "firebrick",
                 linetype = 2, linewidth = 0.35) +
      geom_linerange(aes(ymin = lo, ymax = hi, colour = kind), linewidth = 0.75, na.rm = TRUE) +
      geom_point(aes(y = mid, colour = kind), size = 0.9, na.rm = TRUE) +
      scale_x_discrete(limits = rev(lev), labels = short) +
      ## prof spans its three onset regimes as a single upper strip, so the two kinetics
      ## profiles read as groups rather than as six unrelated columns.
      facet_span(parameter ~ prof + onset, scales = "free_y", switch = "y",
                 labeller = rec_labeller()) +
      labs(x = NULL, y = NULL, title = ttl, subtitle = sub) +
      base(7) +
      theme(strip.text.x = element_text(face = "bold", size = 6.5),
            strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 6.5),
            strip.placement = "outside",
            axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = 5),
            panel.spacing.x = unit(3, "pt"),
            panel.spacing.y = unit(2, "pt"),
            plot.title = element_text(size = 9),
            plot.subtitle = element_text(size = 7))
    w <- 8.27; h <- 11.69                       # A4 portrait

    ## Symmetric y: each parameter's row is centred on the midpoint of its generating values
    ## and extends equally far above and below. Without it a single distant miss sets one end
    ## of the range and the row reads as though the estimates were systematically off in that
    ## direction -- vl_peak's row runs 7.2 to 8.1 purely because of one 1e3 interval, so a
    ## truth at 8.0 sits almost at the top of its own row.
    ##
    ## Implemented with an invisible geom_blank layer rather than scale limits, because the
    ## scales are free per row and only a layer can expand them panel by panel. The layer
    ## carries no column facet variable, so ggplot replicates it across every column of its
    ## row, which is exactly the intended effect.
    if (symmetric_y) {
      ctr <- tapply(tru$truth, tru$parameter, function(v) (min(v) + max(v)) / 2)
      allv <- rbind(data.frame(parameter = df$parameter,  v = df$lo),
                    data.frame(parameter = df$parameter,  v = df$hi),
                    data.frame(parameter = tru$parameter, v = tru$truth))
      allv <- allv[is.finite(allv$v), ]
      sp   <- split(allv$v, allv$parameter)
      lim  <- do.call(rbind, lapply(names(sp), function(nm) {
        c0 <- ctr[[nm]]
        if (!is.finite(c0) || !length(sp[[nm]])) return(NULL)
        h <- max(abs(sp[[nm]] - c0))
        data.frame(parameter = nm, y = c(c0 - h, c0 + h), stringsAsFactors = FALSE)
      }))
      lim$parameter <- factor(lim$parameter, levels = pars)
      lim$row       <- factor(levels(df$row)[1], levels = levels(df$row))
      p <- p + geom_blank(data = lim, aes(x = row, y = y), inherit.aes = FALSE)
    }

    if (prior_panel == "side") {
      pr1$pp <- factor("prior")
      pr1$parameter <- factor(pr1$parameter, levels = pars)
      q <- ggplot(pr1, aes(x = pp)) +
        geom_linerange(aes(ymin = lo, ymax = hi), colour = "steelblue",
                       linewidth = 0.75, na.rm = TRUE) +
        geom_point(aes(y = mid), colour = "steelblue", size = 0.9, na.rm = TRUE) +
        ## Its own y scale, drawn on the right so the two sets of numbers do not collide.
        scale_y_continuous(position = "right") +
        facet_grid(parameter ~ pp, scales = "free_y") +
        labs(x = NULL, y = NULL) +
        theme_bw(base_size = 7) +
        theme(panel.grid = element_blank(), panel.border = element_blank(),
              axis.line = element_line(colour = "black"),
              strip.background = element_blank(),
              strip.text.x = element_text(face = "bold", size = 6.5),
              strip.text.y = element_blank(),
              axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 5, colour = "black"),
              panel.spacing.y = unit(2, "pt"))
      if (!requireNamespace("patchwork", quietly = TRUE))
        stop("prior_panel = 'side' needs the patchwork package")
      ## Widths in proportion to the number of x positions each plot carries, so panel
      ## widths stay comparable rather than the prior panel being stretched.
      p <- patchwork::wrap_plots(p, q, nrow = 1,
                                 widths = c(length(families) * length(lev), 1.6))
    }
  }

  if (!is.null(file)) ggsave(file, p, width = w, height = h, limitsize = FALSE)
  print(p)
  invisible(list(plot = p, table = post))
}

## ---------------------------------------------------------------------------
## figure_prior_sensitivity -- what widening the kinetics priors does to the posterior
##
##   figure_prior_sensitivity()
##   figure_prior_sensitivity(file = "outputs/figures/figS_prior_sensitivity.png")
##
## The real-data counterpart to R1.4b. Each panel is one parameter; each row is the prior or
## one variant's posterior, carrying BOTH prior specifications as a pair of offset bars whose
## colour says which. If the data dominates, a variant's two bars sit on top of each other
## however far apart the two bars on the prior row are -- which is the claim the figure exists
## to test, and pairing them on one row is what makes it readable at a glance.
##
## The specification, not prior-versus-posterior, gets the colour, because the comparison the
## figure is making is within a row. Prior and posterior are drawn identically -- both are a
## median and a 95% interval -- and told apart by the row label and a rule beneath it.
##
## The widened fits carry their own analysis window as x$analysis_from, because two of them
## re-adapted the proposal covariance mid-run and only the draws after that point come from a
## single kernel. Reading it off the file rather than hard-coding it here means the figure
## cannot silently disagree with the decision recorded alongside the chains.
## inc3 (vl_thresh) is deliberately absent. Its posterior sits against its lower bound of zero
## in every variant, so the panel would show four intervals pinned to the axis and a prior
## forty times wider than the plotted range -- a display of a boundary, not of prior
## sensitivity. It is discussed where the parameter is defined instead. Pass it explicitly via
## `pars` to put it back.
##
## At ncol = 3 the six widened parameters fill the first two rows and the two unwidened
## symptom-timing parameters sit alone on the third, which is the grouping the figure wants.
PRIOR_SENS_PARS <- c("a_bar", "a_sigma", "b_bar", "b_sigma", "vl_max_bar", "vl_max_sigma",
                     "inc1", "inc2")

PRIOR_SENS_COL <- c("Original priors" = "grey15", "Widened priors" = "#D95F02")

## Post-burn-in draws for the requested columns only. Reading whole chains and subsetting
## afterwards costs a few hundred MB per variant here.
.draws_of <- function(f, pars, thin, from = NULL) {
  o  <- readRDS(f)
  fl <- max(which(!is.na(o$MCMC_posteriors)))
  b  <- if (!is.null(from)) from else
        analysis_start(o)
  ix <- seq.int(b, fl, by = thin)
  d  <- on_natural_scale(o$MCMC_output[ix, intersect(pars, colnames(o$MCMC_output)), drop = FALSE], o$x)
  rm(o); invisible(gc(verbose = FALSE))
  d
}

## Both directories follow the active form, so the figure cannot compare across forms.
figure_prior_sensitivity <- function(strains  = STRAINS,
                                     pars     = PRIOR_SENS_PARS,
                                     dir_std  = STRAIN_DIR,
                                     dir_wide = WIDE_DIR,
                                     n_prior  = 50000L, thin = 10L,
                                     max_prior_ratio = 25, ncol = 3, file = NULL) {

  q3 <- function(v) stats::quantile(v, c(0.025, 0.5, 0.975), na.rm = TRUE)
  row_of <- function(d, lab, spec) {
    do.call(rbind, lapply(intersect(pars, colnames(d)), function(p) {
      z <- q3(d[, p])
      data.frame(parameter = p, row = lab, spec = spec,
                 lo = z[[1]], mid = z[[2]], hi = z[[3]], row.names = NULL)
    }))
  }
  SPEC <- names(PRIOR_SENS_COL)

  out <- list()

  ## Priors. The kinetics block comes from one longitudinal fit and the inc priors are fixed
  ## constants, so the prior is common to all four variants and one row suffices for both
  ## specifications. Taken from wild type's target in each directory.
  out[["p0"]] <- row_of(prior_draws(strains[1], n_prior, dir_std,  pars), "Prior", SPEC[1])
  out[["p1"]] <- row_of(prior_draws(strains[1], n_prior, dir_wide, pars), "Prior", SPEC[2])

  for (s in strains) {
    lab <- unname(STRAIN_LABEL[s])

    ## Standard fits: three chains, pooled.
    d <- do.call(rbind, lapply(chain_tags(s, dir_std), function(tg)
      .draws_of(chain_file(s, tg, dir_std), pars, thin)))
    out[[paste0(s, "_std")]] <- row_of(d, lab, SPEC[1])
    rm(d); invisible(gc(verbose = FALSE))

    ## Widened fit: one chain, over its own recorded window.
    fw <- list.files(dir_wide, sprintf("^%s_c_[0-9]+\\.rds$", s), full.names = TRUE)
    if (!length(fw)) { message("no widened fit for ", s); next }
    out[[paste0(s, "_wide")]] <- row_of(.draws_of(fw[1], pars, thin), lab, SPEC[2])
  }

  df <- do.call(rbind, out)

  ## A prior far wider than every posterior flattens the panel onto a point. vl_thresh is the
  ## case here: its prior spans 0 to 12 while every posterior sits below 0.4.
  isp  <- df$row == "Prior"
  pw   <- tapply(df$hi[!isp] - df$lo[!isp], df$parameter[!isp], max)
  wide <- unique(df$parameter[isp &
                 (df$hi - df$lo) / pmax(pw[df$parameter], 1e-12) > max_prior_ratio])
  df[isp & df$parameter %in% wide, c("lo", "mid", "hi")] <- NA

  lv <- c("Prior", unname(STRAIN_LABEL[strains]))
  df$row       <- factor(df$row, levels = lv)
  df$parameter <- factor(df$parameter, levels = pars)
  df$spec      <- factor(df$spec, levels = SPEC)

  ## The two specifications share a row, so they are offset by hand on a continuous axis
  ## rather than dodged. position_dodge() dodges along x and would have to be coaxed into the
  ## horizontal orientation; an explicit offset is one line and cannot be misread. Row 1 is
  ## drawn at the top, so the level order is reversed into the numeric position.
  df$y <- length(lv) - as.integer(df$row) + 1 +
          ifelse(df$spec == SPEC[1], 0.16, -0.16)

  ## Two things a reader would otherwise have to infer: which panels have a prior too wide to
  ## draw, and which parameters were not widened at all (inc1-inc3, whose prior/posterior SD
  ## ratios are 18-424 -- the likelihood already dominates them). On those the two prior rows
  ## coincide by construction, and saying so is better than leaving it to be noticed.
  ## Compare interval WIDTHS, not endpoints: a relative test on the endpoints would call
  ## vl_peak unchanged, because its prior sits near 7.5 and doubling its SD moves each end by
  ## under 3% of that. The widths, 0.39 against 0.78, are unambiguous.
  same <- vapply(pars, function(p) {
    z <- df[df$parameter == p & isp, ]
    if (nrow(z) != 2L) return(FALSE)
    if (all(is.na(z$lo))) return(TRUE)
    w <- z$hi - z$lo
    abs(w[2] / w[1] - 1) < 0.1
  }, logical(1))

  lab <- ifelse(pars %in% wide,
                ifelse(same[pars], "prior off scale; not widened", "prior wider than axis"),
                ifelse(same[pars], "prior not widened", NA_character_))
  note <- if (any(!is.na(lab)))
    data.frame(parameter = factor(pars[!is.na(lab)], levels = pars),
               y = length(lv), label = lab[!is.na(lab)], row.names = NULL)

  p <- ggplot(df, aes(y = y, colour = spec)) +
    ## Prior and posteriors are drawn alike -- they are the same kind of quantity, a median and
    ## a 95% interval, and the row label says which is which. The rule separates the prior row
    ## from the variants so the eye does not read it as a fifth variant.
    geom_hline(yintercept = length(lv) - 0.5, linewidth = 0.25, colour = "grey70") +
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 1.0, na.rm = TRUE) +
    geom_point(aes(x = mid), size = 1.5, na.rm = TRUE) +
    { if (!is.null(note)) geom_text(data = note, aes(x = -Inf, y = y, label = label),
                                    inherit.aes = FALSE, hjust = -0.05, vjust = -1.1,
                                    size = 2.4, colour = "grey35") } +
    scale_y_continuous(breaks = length(lv):1, labels = lv,
                       expand = expansion(add = c(0.5, 0.9))) +
    scale_colour_manual(values = PRIOR_SENS_COL, name = NULL, drop = FALSE) +
    facet_wrap(~parameter, scales = "free_x", ncol = ncol, labeller = rec_labeller()) +
    labs(x = NULL, y = NULL,
         title = "Sensitivity of the posterior to the width of the kinetics priors",
         subtitle = paste("kinetics prior standard deviations doubled; points are medians,",
                          "bars 95% intervals.",
                          "Top row is the prior, remaining rows each variant's posterior")) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.line.x = element_line(colour = "black"),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold"),
          axis.text = element_text(colour = "black"),
          legend.position = "top")

  if (!is.null(file))
    ggsave(file, p, width = 9, height = 1.35 * ceiling(length(pars) / ncol) + 1.4)
  print(p)
  invisible(list(plot = p, table = df))
}

## ---------------------------------------------------------------------------
## table_ess -- convergence diagnostics per parameter and variant (R1.m5, U4)
##
##   table_ess()                       # data frame, long
##   table_ess(latex = TRUE)           # LaTeX for the Supplementary
##
## Same construction as check_chains(): each chain cut at its own burn-in, then all truncated
## to a common length taken from the END, so a chain that settled late contributes only
## stationary draws and none is over-weighted. Reusing that path matters -- a table that
## disagreed with the printed diagnostics by a different windowing rule would be worse than no
## table.
##
## The spline coefficients are summarised by their WORST value rather than listed. There are
## 20 to 30 of them per variant, they are nuisance parameters, and reporting only the kinetics
## would invite the obvious question of what was left out. One row answers it.
ESS_ROWS <- c("log posterior" = "log posterior",
              a_bar = "$\\eta_a$", a_sigma = "$\\omega_a$",
              b_bar = "$\\eta_b$", b_sigma = "$\\omega_b$",
              vl_max_bar = "$\\overline{vl}_{\\mathrm{peak}}$",
              vl_max_sigma = "$vl_{\\sigma}$",
              inc1 = "$\\eta_T$", inc2 = "$\\omega_T$",
              inc3 = "$vl_{\\mathrm{thresh}}$",
              rdisp = "$\\xi$", multiplier = "multiplier",
              spline = "spline coefficients (worst)")

## show_tail = FALSE reports bulk effective sample size only, which is the conventional
## minimum alongside R-hat and is what R1.m5 actually asks about. The tail values are not
## dropped silently: the caption gives the smallest one and says where it occurs, because the
## two weakest numbers in the study are tail values and a reviewer recomputing them should
## find them acknowledged rather than absent.
## show_rhat = TRUE gives each variant a pair of columns, split R-hat and bulk ESS, instead of
## ESS alone with the R-hat bounds relegated to the caption. That is the form to use when the
## table is the evidence for R1.m5, which is a question about R-hat: a reader asking whether the
## 1.01 threshold is met should be able to see every value rather than a maximum per variant.
## rhat_type picks WHICH Rhat, and the captions take their wording from the same choice, so a
## table can no longer print one statistic while naming another. The three differ materially
## here: over the four variants the maxima are 1.0070 (classic), 1.0108 (split) and 1.0073
## (rank-normalised split), so plain split-Rhat is the only one that breaches 1.01 -- which is
## what the earlier caption claimed not to happen while calling the statistic "split".
RHAT_LABEL <- c(
  classic = "$\\hat{R}$ statistic of \\citet{gelman1992inference}",
  split   = "split $\\hat{R}$ statistic",
  rank    = "rank-normalised split $\\hat{R}$ statistic \\citep{vehtari2021rank}")

table_ess <- function(strains = c("wt", "alpha", "delta", "omicron"),
                      dir = STRAIN_DIR, latex = FALSE, digits = 0,
                      show_tail = FALSE, show_rhat = FALSE,
                      rhat_type = c("rank", "split", "classic")) {
  rhat_type <- match.arg(rhat_type)
  rhat_fn <- switch(rhat_type,
    classic = function(m) posterior::rhat_basic(m, split = FALSE),
    split   = function(m) posterior::rhat_basic(m, split = TRUE),
    rank    = function(m) posterior::rhat(m))

  one <- function(s) {
    tags <- chain_tags(s, dir)
    ch <- lapply(tags, function(tg) {
      o  <- readRDS(chain_file(s, tg, dir))
      fl <- max(which(!is.na(o$MCMC_posteriors)))
      b  <- analysis_start(o)
      r  <- list(lp = o$MCMC_posteriors[1:fl, 1], out = o$MCMC_output[1:fl, , drop = FALSE],
                 fl = fl, b = b,
                 ## The middle spline knot is pinned at its prior mean to identify the incidence
                 ## scale against `multiplier` (model.R:461), so it is not a
                 ## sampled parameter and must not enter the diagnostics. Alpha's chain 3 holds it
                 ## 3.6e-07 from the other two, which is zero within-chain variance with non-zero
                 ## between-chain variance -- precisely rhat = Inf, and it reached supp_table_ess.
                 pin = if (is.null(o$x$ignored_spline_par)) character(0)
                       else names(o$x$parameters)[o$x$ignored_spline_par])
      rm(o); invisible(gc(verbose = FALSE)); r
    })
    n   <- min(vapply(ch, function(z) z$fl - z$b + 1L, integer(1)))
    nm  <- colnames(ch[[1]]$out)
    pin <- unique(unlist(lapply(ch, `[[`, "pin")))
    nch <- length(ch)
    dk  <- setdiff(c(intersect(setdiff(names(ESS_ROWS), c("log posterior", "spline")), nm),
                     grep("^theta", nm, value = TRUE)), pin)
    arr <- array(NA_real_, c(n, length(ch), length(dk) + 1L),
                 dimnames = list(NULL, seq_along(ch), c("log posterior", dk)))
    for (j in seq_along(ch)) {
      ix <- (ch[[j]]$fl - n + 1L):ch[[j]]$fl
      arr[, j, ] <- cbind(ch[[j]]$lp[ix], ch[[j]]$out[ix, dk, drop = FALSE])
    }
    rm(ch); invisible(gc(verbose = FALSE))
    ## Computed per parameter rather than via summarise_draws("rhat"), which is hardwired to the
    ## rank-normalised version and so could not honour rhat_type.
    d <- do.call(rbind, lapply(dimnames(arr)[[3]], function(p) {
      m <- arr[, , p, drop = TRUE]
      data.frame(variable = p, rhat = rhat_fn(m),
                 ess_bulk = posterior::ess_bulk(m), ess_tail = posterior::ess_tail(m),
                 stringsAsFactors = FALSE)
    }))

    ## Collapse the spline block to its worst row on each diagnostic.
    th <- grepl("^theta", d$variable)
    if (any(th))
      d <- rbind(d[!th, ],
                 data.frame(variable = "spline",
                            rhat     = max(d$rhat[th],     na.rm = TRUE),
                            ess_bulk = min(d$ess_bulk[th], na.rm = TRUE),
                            ess_tail = min(d$ess_tail[th], na.rm = TRUE)))
    d$strain <- s
    d$n_draws <- n
    d$n_chains <- nch
    d
  }

  d <- do.call(rbind, lapply(strains, one))
  d$variable <- factor(d$variable, levels = names(ESS_ROWS))
  d <- d[order(d$variable, match(d$strain, strains)), ]
  if (!latex) return(d)

  nd  <- tapply(d$n_draws, d$strain, max)[strains]
  lab <- unname(STRAIN_LABEL[strains])

  if (show_rhat) {
    ## Two columns per variant. Rows keep ESS_ROWS order; parameters absent from a fit (none at
    ## present, but the function is called on other directories) drop out rather than print NA.
    grab <- function(col, fmt)
      vapply(names(ESS_ROWS), function(p)
        paste(vapply(strains, function(s) {
          z <- d[d$strain == s & as.character(d$variable) == p, col]
          if (!length(z) || is.na(z)) "---" else sprintf(fmt, z)
        }, character(1)), collapse = " & "), character(1))

    rr <- grab("rhat", "%.3f")
    ee <- vapply(names(ESS_ROWS), function(p)
      paste(vapply(strains, function(s) {
        z <- d[d$strain == s & as.character(d$variable) == p, "ess_bulk"]
        if (!length(z) || is.na(z)) "---" else format(round(z), big.mark = ",", trim = TRUE)
      }, character(1)), collapse = " & "), character(1))

    keep <- names(ESS_ROWS)[!grepl("^-+$", gsub("[ &]", "", rr))]
    body <- vapply(keep, function(p) {
      r <- strsplit(rr[[p]], " & ", fixed = TRUE)[[1]]
      e <- strsplit(ee[[p]], " & ", fixed = TRUE)[[1]]
      paste0(unname(ESS_ROWS[p]), " & ",
             paste(sprintf("%s & %s", r, e), collapse = " & "), " \\\\")
    }, character(1))

    hdr <- c(paste0(" & ", paste(sprintf("\\multicolumn{2}{c}{%s}", lab), collapse = " & "),
                    " \\\\"),
             paste(sprintf("\\cmidrule(lr){%d-%d}", 2 * seq_along(strains),
                           2 * seq_along(strains) + 1), collapse = ""),
             paste0("Parameter & ",
                    paste(rep("$\\hat{R}$ & ESS", length(strains)), collapse = " & "), " \\\\"))

    ## Both the chain count and the Rhat claim are computed, not asserted. Omicron has four
    ## chains where the others have three, and the hardcoded "no parameter reaches 1.01" was
    ## false in the generated table it described (alpha omega_T 1.014, spline Inf).
    nc  <- tapply(d$n_chains, d$strain, max)[strains]
    rmx <- max(d$rhat, na.rm = TRUE)
    wr  <- d[which.max(d$rhat), ]
    cap <- paste0(
      "\\textbf{Convergence diagnostics for the variant fits.} The ",
      unname(RHAT_LABEL[rhat_type]), " and bulk ",
      "effective sample size (ESS) for every sampled parameter, the log posterior, and the ",
      "worst of the spline coefficients, computed across the ",
      if (length(unique(nc)) == 1L) sprintf("%d chains fitted for each variant", nc[[1]])
      else paste0("chains fitted for each variant (",
                  paste(sprintf("%s %d", lab, nc), collapse = ", "), ")"),
      ". The middle spline coefficient is excluded: it is pinned at its prior mean to identify ",
      "the incidence scale and is not sampled. Burn-in was discarded and the chains truncated ",
      "to a common length taken from the end of each (",
      paste(sprintf("%s %s", lab, trimws(format(nd, big.mark = ","))), collapse = ", "),
      " draws per chain respectively). ",
      if (rmx < 1.01) sprintf(
        "No parameter of any variant reaches $\\hat{R}=1.01$; the largest is %.4f (%s, %s). ",
        rmx, unname(STRAIN_LABEL[wr$strain]), unname(ESS_ROWS[as.character(wr$variable)]))
      else sprintf("The largest $\\hat{R}$ is %.4f (%s, %s). ", rmx,
                   unname(STRAIN_LABEL[wr$strain]),
                   unname(ESS_ROWS[as.character(wr$variable)])),
      "Requires the booktabs package for \\texttt{\\textbackslash cmidrule}.")

    cat(tex_wrap(cap, "tab:ess",
                 paste0("l", strrep(" rr", length(strains))), hdr, body), sep = "\n")
    return(invisible(d))
  }

  ## Wide: rows are parameters, one column per variant.
  cell <- if (show_tail)
    sprintf(paste0("%.", digits, "f / %.", digits, "f"), d$ess_bulk, d$ess_tail)
  else sprintf(paste0("%.", digits, "f"), d$ess_bulk)
  w <- reshape(data.frame(variable = d$variable, strain = d$strain, cell = cell),
               idvar = "variable", timevar = "strain", direction = "wide")
  names(w) <- sub("^cell\\.", "", names(w))
  body <- paste(apply(w, 1, function(r)
    paste0(unname(ESS_ROWS[as.character(r[["variable"]])]), " & ",
           paste(r[strains], collapse = " & "), " \\\\")), collapse = "\n")

  rmax <- tapply(d$rhat, d$strain, max, na.rm = TRUE)[strains]
  nc2  <- tapply(d$n_chains, d$strain, max)[strains]
  wt_min <- d[which.min(d$ess_tail), ]
  cap <- paste0(
    "\\textbf{Effective sample size by parameter and variant.} ",
    if (show_tail) "Each cell gives the bulk and tail effective sample size, computed "
    else "Each cell gives the bulk effective sample size, the number of independent draws the posterior mean is worth, computed ",
    "across the ",
    if (length(unique(nc2)) == 1L) sprintf("%d chains of each variant", nc2[[1]])
    else paste0("chains of each variant (", paste(sprintf("%s %d", lab, nc2), collapse = ", "),
                ")"),
    " after discarding burn-in and truncating all ",
    "chains to a common length taken from the end of each (",
    paste(sprintf("%s %s", lab, trimws(format(nd, big.mark = ","))), collapse = ", "),
    " draws per chain respectively). The middle spline coefficient is pinned at its prior mean ",
    "to identify the incidence scale and is excluded; the remaining spline coefficients are ",
    "nuisance parameters and are summarised by their worst value. The ",
    unname(RHAT_LABEL[rhat_type]), " did not exceed ",
    paste(sprintf("%.3f for %s", rmax, lab), collapse = ", "), ".",
    if (!show_tail) paste0(
      " Tail effective sample size, which governs the precision of the interval edges rather ",
      "than of the mean, was lowest for ",
      unname(ESS_ROWS[as.character(wt_min$variable)]), " under ",
      unname(STRAIN_LABEL[wt_min$strain]), " at ", round(wt_min$ess_tail),
      "; every other value exceeded ",
      round(min(d$ess_tail[-which.min(d$ess_tail)]) / 100) * 100, "."))

  cat(tex_wrap(cap, "tab:ess", paste0("l", strrep("r", length(strains))),
               paste0("Parameter & ", paste(lab, collapse = " & "), " \\\\"), body))
  invisible(d)
}

## ---------------------------------------------------------------------------
## figure_modes -- one primary-mode and one secondary-mode chain per variant
##
##   figure_modes()
##   figure_modes(file = file.path(FIG_DIR, "figS_modes.png"))
##
## The evidence behind the multimodality paragraph. The census counts how many chains reached
## each mode; this shows what reaching one looks like, which is the part a reader cannot take
## on trust from a proportion.
##
## Alpha and Omicron are the only two variants with a secondary mode, and they are NOT the same
## phenomenon: Alpha's is slower growth with a long onset delay, Omicron's is faster growth with
## a short delay and a tenfold higher threshold viral load. Putting them side by side is the
## point of the figure -- the shared signature is the inflated omega_a row, and every other row
## moves in opposite directions.
##
## Each panel gets its own y scale via facet_grid2(independent = "y"). plain facet_grid frees
## scales per ROW, which would force Alpha and Omicron onto one axis; their log posteriors are
## -23,700 and -6,200, so the two columns would be unreadable.
## ---------------------------------------------------------------------------
MODE_COL <- c(Primary = "#1F5FA9", Secondary = "#C4342B")

## Chains chosen as the longest representative of each mode, from the census over outputs/ and
## outputs/old. Alpha's secondary chains (i31/i33/i34) agree to within 0.02 on omega_a, so i31
## stands for all three; Omicron's secondary mode was reached by i33 alone.
MODE_CHAINS <- list(
  Alpha   = c(Primary   = "outputs/alpha_i32_195000.rds",
              Secondary = "outputs/old/alpha_i31_155000.rds"),
  Omicron = c(Primary   = "outputs/omicron_i31_190000.rds",
              Secondary = "outputs/omicron_i33_195000.rds"))

figure_modes <- function(chains = MODE_CHAINS, pars = REC_PARS, max_points = 2500L,
                         lp_drop = 6, file = NULL) {

  d <- do.call(rbind, lapply(names(chains), function(st)
    do.call(rbind, lapply(names(chains[[st]]), function(md) {
      p <- chains[[st]][[md]]
      if (!file.exists(p)) stop("no such chain file: ", p, call. = FALSE)
      o <- readRDS(p); x <- o$x
      ## A trace from a different form or truncation would be a different model, and the whole
      ## claim is that these two chains target the same posterior.
      if (is.null(x$form) || x$form != "incidence" ||
          is.null(x$symp_delay_lim) || as.integer(x$symp_delay_lim) != 28L)
        stop(p, " is not an incidence/sdl-28 fit", call. = FALSE)
      fl   <- max(which(!is.na(o$MCMC_posteriors[, 1])))
      keep <- unique(round(seq(1, fl, length.out = min(fl, max_points))))
      m    <- on_natural_scale(o$MCMC_output[keep, , drop = FALSE], x)
      miss <- setdiff(pars, colnames(m))
      if (length(miss))
        stop("parameters absent from ", p, ": ", paste(miss, collapse = ", "), call. = FALSE)
      out <- rbind(
        data.frame(strain = st, mode = md, iteration = keep, parameter = "log posterior",
                   value = o$MCMC_posteriors[keep, 1], stringsAsFactors = FALSE),
        do.call(rbind, lapply(pars, function(q)
          data.frame(strain = st, mode = md, iteration = keep, parameter = q,
                     value = m[, q], stringsAsFactors = FALSE))))
      rm(o, x, m); invisible(gc(verbose = FALSE))
      out
    }))))

  ## The log posterior climbs hundreds of units during burn-in and the gap between modes is a
  ## few tens, so an unclipped axis compresses the thing the figure exists to show into the top
  ## few percent of the panel. The floor is set relative to the gap itself -- lp_drop times it
  ## below the lower mode -- so it adapts to each variant instead of being a magic number.
  ##
  ## Off-scale draws are dropped to NA, NOT clamped to the floor. Clamping draws them as a flat
  ## line along the bottom of the panel, which looks exactly like a chain stuck at a low value
  ## and is the opposite of what is happening; NA leaves the trace to enter the panel when it
  ## climbs into range, which is what a reader already expects a truncated axis to mean.
  lp <- d[d$parameter == "log posterior", ]
  set <- do.call(rbind, lapply(split(lp, list(lp$strain, lp$mode), drop = TRUE), function(k)
    data.frame(strain = k$strain[1], mode = k$mode[1],
               level = median(k$value[k$iteration >= 0.9 * max(k$iteration)]),
               stringsAsFactors = FALSE)))
  n_clip <- 0L
  for (st in names(chains)) {
    lv  <- set$level[set$strain == st]
    flr <- min(lv) - lp_drop * max(abs(diff(lv)), 1)
    i   <- d$parameter == "log posterior" & d$strain == st & d$value < flr
    n_clip <- n_clip + sum(i)
    d$value[i] <- NA_real_
  }
  if (n_clip)
    message(sprintf("log posterior axis truncated %g x the mode gap below the lower mode: %d of %d drawn points fall below it and are not shown",
                    lp_drop, n_clip, sum(d$parameter == "log posterior")))
  for (st in names(chains)) {
    lv <- set$level[set$strain == st]
    message(sprintf("  %-8s modes at %.1f and %.1f, gap %.1f log posterior units",
                    st, max(lv), min(lv), max(lv) - min(lv)))
  }

  d$parameter <- factor(d$parameter, levels = c("log posterior", pars))
  d$mode      <- factor(d$mode,   levels = names(MODE_COL))
  d$strain    <- factor(d$strain, levels = names(chains))

  p <- ggplot(d, aes(iteration / 1000, value, colour = mode)) +
    geom_line(linewidth = 0.18) +
    scale_colour_manual(values = MODE_COL, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    ## facet_span routes to facet_nested, which has no `independent` argument -- this needs
    ## facet_grid2 specifically. Degrade to a shared row scale rather than failing outright.
    (if (requireNamespace("ggh4x", quietly = TRUE))
       ggh4x::facet_grid2(parameter ~ strain, scales = "free_y", independent = "y",
                          labeller = rec_labeller())
     else {
       warning("ggh4x absent: Alpha and Omicron share a y scale per row, which makes the ",
               "log posterior row unreadable", call. = FALSE)
       facet_grid(parameter ~ strain, scales = "free_y", labeller = rec_labeller())
     }) +
    labs(x = "Iteration (thousands)", y = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey94", colour = NA),
          strip.text.y = element_text(angle = 0))

  if (!is.null(file))
    ggsave(file, p, width = 7.5, height = 1.05 * (length(pars) + 1) + 1, dpi = 300,
           limitsize = FALSE)
  p
}

## ---------------------------------------------------------------------------
## figure_traces -- every reported chain, all four variants, one grid
##
##   figure_traces()                       # the analysis window, which is what the paper uses
##   figure_traces(window = "all")         # the whole history, adaptation boundary marked
##   figure_traces(file = file.path(FIG_DIR, "figS_traces.png"))
##
## The companion to the R-hat and ESS tables: those reduce mixing to two numbers per parameter,
## this shows what the numbers are summarising. Chains are coloured within a variant, so the
## question each column answers is whether its chains are on top of one another.
##
## The default window is x$analysis_from -- the draw after the final covariance adaptation --
## and NOT the whole chain. Every posterior summary in the paper comes from that window, so it
## is the window whose stationarity matters; including burn-in would also force the log
## posterior axis to span hundreds of units and flatten the part being judged. window = "all"
## shows the full history with the window start marked, for when the question is how the chain
## got there rather than whether it has arrived.
## ---------------------------------------------------------------------------
## `alpha` matters more than it looks: chains are drawn in order, so without transparency the
## last one painted hides the rest and a column of four well-mixed chains reads as one chain
## changing colour partway through -- which is the opposite of what the figure is for.
figure_traces <- function(strains = STRAINS, pars = REC_PARS, dir = STRAIN_DIR,
                          window = c("analysis", "all"), max_points = 2000L,
                          alpha = 0.55, file = NULL) {
  window <- match.arg(window)

  starts <- NULL
  d <- do.call(rbind, lapply(strains, function(st) {
    tags <- chain_tags(st, dir)
    do.call(rbind, lapply(seq_along(tags), function(i) {
      o <- readRDS(chain_file(st, tags[i], dir)); x <- o$x
      fl <- max(which(!is.na(o$MCMC_posteriors[, 1])))
      b  <- analysis_start(o)
      starts <<- rbind(starts, data.frame(strain = st, chain = sprintf("%d", i),
                                          b = b, stringsAsFactors = FALSE))
      rng  <- if (window == "analysis") b:fl else 1:fl
      keep <- unique(round(seq(min(rng), max(rng), length.out = min(length(rng), max_points))))
      m    <- on_natural_scale(o$MCMC_output[keep, , drop = FALSE], x)
      miss <- setdiff(pars, colnames(m))
      if (length(miss))
        stop("parameters absent from ", st, " chain ", tags[i], ": ",
             paste(miss, collapse = ", "), call. = FALSE)
      out <- rbind(
        data.frame(strain = st, chain = sprintf("%d", i), iteration = keep,
                   parameter = "log posterior", value = o$MCMC_posteriors[keep, 1],
                   stringsAsFactors = FALSE),
        do.call(rbind, lapply(pars, function(q)
          data.frame(strain = st, chain = sprintf("%d", i), iteration = keep, parameter = q,
                     value = m[, q], stringsAsFactors = FALSE))))
      rm(o, x, m); invisible(gc(verbose = FALSE))
      out
    }))
  }))

  ## Only needed for window = "all": with burn-in included the log posterior spans hundreds of
  ## units and the stationary part collapses to a sliver. Dropped to NA rather than clamped,
  ## so an off-scale excursion leaves a gap instead of a false flat line.
  n_drop <- 0L
  if (window == "all") {
    lp <- d[d$parameter == "log posterior", ]
    for (st in strains) {
      v   <- lp$value[lp$strain == st]
      flr <- min(v, na.rm = TRUE)
      q   <- stats::quantile(v, c(0.5, 1), na.rm = TRUE)
      flr <- q[[1]] - 6 * (q[[2]] - q[[1]])
      i   <- d$parameter == "log posterior" & d$strain == st & d$value < flr
      n_drop <- n_drop + sum(i, na.rm = TRUE)
      d$value[i] <- NA_real_
    }
    if (n_drop)
      message(sprintf("%d of %d drawn log posterior points fall below the truncated axis and are not shown",
                      n_drop, sum(d$parameter == "log posterior")))
  }

  for (st in strains) {
    k <- d[d$strain == st & d$parameter == "log posterior", ]
    message(sprintf("  %-8s %d chains, %s draws each after thinning, iterations %s-%s",
                    st, length(unique(k$chain)),
                    format(sum(k$chain == k$chain[1]), big.mark = ","),
                    format(min(k$iteration), big.mark = ","),
                    format(max(k$iteration), big.mark = ",")))
  }

  d$parameter <- factor(d$parameter, levels = c("log posterior", pars))
  d$strain    <- factor(d$strain, levels = strains, labels = unname(STRAIN_LABEL[strains]))
  starts$strain <- factor(starts$strain, levels = strains,
                          labels = unname(STRAIN_LABEL[strains]))

  p <- ggplot(d, aes(iteration / 1000, value, colour = chain)) +
    ## Only for window = "all": everything left of a chain's line is discarded before any
    ## posterior summary, so a reader can see which part of the trace the paper actually used.
    ## Same convention as check_chains -- dashed, in the chain's own colour.
    (if (window == "all")
       geom_vline(data = starts, aes(xintercept = b / 1000, colour = chain),
                  linetype = 2, linewidth = 0.3, show.legend = FALSE)
     else NULL) +
    geom_line(linewidth = 0.15, alpha = alpha, na.rm = TRUE) +
    scale_colour_brewer(palette = "Dark2", name = "Chain") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    (if (requireNamespace("ggh4x", quietly = TRUE))
       ggh4x::facet_grid2(parameter ~ strain, scales = "free", independent = "all",
                          labeller = rec_labeller())
     else {
       warning("ggh4x absent: variants share a y scale per row, which makes the log ",
               "posterior row unreadable", call. = FALSE)
       facet_grid(parameter ~ strain, scales = "free_y", labeller = rec_labeller())
     }) +
    labs(x = "Iteration (thousands)", y = NULL) +
    guides(colour = guide_legend(override.aes = list(linewidth = 1.2))) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey94", colour = NA),
          strip.text.y = element_text(angle = 0))

  if (!is.null(file))
    ggsave(file, p, width = 2.6 * length(strains) + 0.8,
           height = 1.05 * (length(pars) + 1) + 1, dpi = 300, limitsize = FALSE)
  p
}

## ---------------------------------------------------------------------------
## figure_scenarios -- what the six simulation scenarios actually look like
##
##   figure_scenarios()
##   figure_scenarios(file = "outputs/figures/figS_scenarios.png")
##
## Ported from plot_trajectories() in simulation_study.R, which is an exploratory script with
## stray top-level calls and reads the superseded _mult files. E3/R1.5 rests on every figure
## regenerating from this script, so a supplementary figure sourced from there would undercut
## the reproducibility claim.
##
## SIX panels, not 24. ladder() shifts only the spline coefficients, so the four rungs of a
## family share their kinetics exactly and M0 draws the same trajectories as M3; the distinct
## parameter sets are two kinetics profiles by three onset regimes.
##
## Changes from the original:
##   - only the onset-aligned frame is drawn. The original also aligned on peak viral load,
##     which is a statement about biology; this figure is about what the DATA can see, and
##     twelve panels would bury that.
##   - the horizontal reference is the limit of detection vl_0 = 1.5 rather than the plotting
##     floor, so the line means something in the model.
##   - percentages of individuals whose onset precedes peak are returned rather than printed,
##     so the caption can quote them without re-running.
##
## The shaded band is the 0-6 day testing window. It is the point of the figure: everything
## the likelihood sees lies inside it, and the whole later decline is extrapolated. That is
## the visual form of the clearance-phase caveat in the Discussion and of R1.m2.
figure_scenarios <- function(families = c("P", "Q", "R", "M", "N", "O"),
                             n = 200, n_pct = 200000L,
                             xlim = c(-13, 28),
                             vl0 = 1.5,
                             seed = 1, file = NULL) {

  truth <- readRDS(TRUTH_FILE)
  set.seed(seed)

  ## Individual-level draws for one scenario. Used twice: n for the lines, n_pct for the
  ## percentage, because the percentage is quoted in the caption and must not move with the
  ## seed. At n = 200 its standard error is about 3.5 points, which is enough to shift a
  ## reported "76%" to "70%" on a re-run.
  draw <- function(p, m) {
    a       <- 0.5  + exp(rnorm(m, p[["a_bar"]], p[["a_sigma"]]))
    b       <- 0.25 + exp(rnorm(m, p[["b_bar"]], p[["b_sigma"]]))
    vpeak   <- rnorm(m, p[["vl_max_bar"]], p[["vl_max_sigma"]])
    tsv     <- exp(rnorm(m, p[["inc1"]], p[["inc2"]]))   # threshold crossing -> onset, T_i
    ## Time is measured from SYMPTOM ONSET, so t1 is onset -> peak and is negative for
    ## individuals whose symptoms begin after their peak.
    list(a = a, b = b, vpeak = vpeak, tsv = tsv, vthresh = p[["inc3"]],
         t1 = (vpeak - p[["inc3"]]) / a - tsv)
  }

  one <- function(f) {
    ## Any rung would do -- kinetics are identical across the ladder -- so take rung 3.
    p <- truth[[paste0(f, "3")]]
    if (is.null(p)) stop("no generating values for ", f, "3 in ", TRUTH_FILE)

    d <- draw(p, n)

    ## Time is measured from PEAK viral load. In that frame the trajectory is a plain tent,
    ## vl = vpeak + a*s rising and vpeak - b*s falling, because vl_thresh + a*(t1 + T)
    ## collapses to vpeak by construction. Symptom onset sits at s = -t1, so an individual
    ## whose onset precedes their peak has their onset marked to the LEFT of zero -- which is
    ## what the column label asserts, and what the onset-aligned version obscured by putting
    ## every onset at zero and letting the peak move instead.
    ##
    ## Drawn from EXACT VERTICES rather than sampled on a shared time grid. The trajectory is
    ## piecewise linear, so three points describe it completely, and the endpoints can be put
    ## exactly on the limit of detection instead of wherever the grid happened to fall. That
    ## fixes two things at once: every line now starts and ends precisely at vl_0 so the
    ## bottom of the panel is flat, and the figure drops from ~250,000 line vertices to a few
    ## thousand.
    s_up <- (vl0 - d$vpeak) / d$a          # detectable from here (negative)
    s_dn <- (d$vpeak - vl0) / d$b          # undetectable after here (positive)
    vlf  <- function(s, i) ifelse(s <= 0, d$vpeak[i] + d$a[i] * s, d$vpeak[i] - d$b[i] * s)

    traj <- do.call(rbind, lapply(seq_len(n), function(i) {
      s <- c(s_up[i], 0, s_dn[i])
      data.frame(family = f, id = i, t = s, vl = vlf(s, i))
    }))

    ## The 0-6 day post-onset testing window is individual-specific once time is measured from
    ## peak, so it cannot be one shaded rectangle. The observed stretch of each trajectory is
    ## drawn over the top instead, which shows the same thing per individual rather than on
    ## average -- and makes plain that in the post-peak scenarios the growth phase is never
    ## observed at all.
    obs <- do.call(rbind, lapply(seq_len(n), function(i) {
      lo <- max(s_up[i], -d$t1[i]); hi <- min(s_dn[i], -d$t1[i] + 6)
      if (!(hi > lo)) return(NULL)
      s <- sort(unique(c(lo, if (lo < 0 && 0 < hi) 0, hi)))
      data.frame(family = f, id = i, t = s, vl = vlf(s, i))
    }))

    ## Onset markers are dropped where the individual is already undetectable at onset. Their
    ## viral load at that moment is far below the axis -- symptoms beginning three weeks after
    ## peak imply a value of -70 or worse -- and, more to the point, such an individual
    ## returns no positive test and contributes nothing to the data the figure is about. The
    ## threshold-crossing markers are NOT clipped: vl_thresh lies below the limit of detection
    ## in two scenarios, and that the threshold is crossed before detectability is a genuine
    ## feature of the model rather than something to hide.
    vl_on <- ifelse(d$t1 >= 0, d$vpeak - d$a * d$t1, d$vpeak + d$b * d$t1)
    vl_on[vl_on < vl0] <- NA
    marks <- data.frame(
      family = f, id = seq_len(n),
      t_cross  = -d$tsv - d$t1, vl_cross = d$vthresh,
      t_onset  = -d$t1,
      vl_onset = vl_on)
    list(traj = traj, obs = obs, marks = marks,
         pct = 100 * mean(draw(p, n_pct)$t1 >= 0))
  }

  z     <- lapply(families, one); names(z) <- families
  traj  <- do.call(rbind, lapply(z, `[[`, "traj"))
  obs   <- do.call(rbind, lapply(z, `[[`, "obs"))
  marks <- do.call(rbind, lapply(z, `[[`, "marks"))
  pct   <- vapply(z, `[[`, numeric(1), "pct")

  ## Less variable on top, matching figure_recovery_grid's left-to-right ordering.
  lab <- function(d) {
    d$prof  <- factor(sub("^[0-9]+:\\s*", "", unname(PROF_LABEL[d$family])),
                      levels = sub("^[0-9]+:\\s*", "", rev(PROF_LEVELS)))
    d$onset <- onset_factor(d$family)
    d
  }
  traj <- lab(traj); obs <- lab(obs); marks <- lab(marks)

  ## Panel letters, matching the recovery grid's columns, with the percentage of individuals
  ## whose onset precedes their peak. The percentage is what the column strip is asserting,
  ## so stating it in the panel makes the label checkable rather than something the reader has
  ## to take on trust while looking at peaks drawn on the other side of zero.
  pan <- unique(traj[, c("family", "prof", "onset")])
  pan$label <- sprintf("%s   %.0f%% onset before peak",
                       unname(SCEN_PANEL[pan$family]), pct[pan$family])

  p <- ggplot(traj, aes(t, vl, group = id)) +
    geom_line(alpha = 0.08, linewidth = 0.25, colour = "grey25") +
    ## The part of each trajectory the likelihood actually sees.
    geom_line(data = obs, alpha = 0.22, linewidth = 0.45, colour = "steelblue") +
    geom_point(data = marks, aes(t_cross, vl_cross), inherit.aes = FALSE,
               colour = "darkgreen", size = 0.7, alpha = 0.4) +
    geom_point(data = marks, aes(t_onset, vl_onset), inherit.aes = FALSE,
               colour = "purple3", size = 1.0, alpha = 0.55) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
    geom_hline(yintercept = vl0, linetype = 3, colour = "grey40", linewidth = 0.3) +
    ## geom_label, not geom_text: the peak line and the tallest trajectories run straight
    ## through the annotation otherwise.
    geom_label(data = pan, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
               hjust = -0.03, vjust = 1.15, fontface = "bold", size = 3,
               fill = "white", linewidth = 0, label.padding = unit(2, "pt")) +
    ## coord_cartesian rather than scale limits: threshold crossings for the slowest-growing
    ## individuals fall 20 to 30 days before onset, and letting them set the axis leaves the
    ## trajectories squeezed into the right-hand third with stray green points floating to
    ## the left of where any line is drawn.
    coord_cartesian(xlim = xlim) +
    facet_grid(prof ~ onset) +
    labs(x = "days since peak viral load", y = expression(log[10]~"viral load"),
         title = "Simulated within-host trajectories, by scenario",
         subtitle = paste("purple = symptom onset  |  green = threshold crossing  |",
                          "blue = the 0-6 days after onset, within which tests are retained")) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text = element_text(colour = "black"))

  if (!is.null(file)) ggsave(file, p, width = 9, height = 5.6)
  print(p)
  invisible(list(plot = p, pct_onset_before_peak = round(pct, 1)))
}

## ---------------------------------------------------------------------------
## figure_precision -- posterior precision against sample size (R1.4c)
##
##   figure_precision()                      # all six families, saved chains
##   figure_precision(c("M","N","O"))
##   figure_precision(file = "figS_precision.png")
##
## The direct answer to "downsampling: 10^3/10^4/10^5/10^6 -- precision vs sample size".
## Posterior SD against the number of positive tests, both on log scales, one panel per
## parameter and one line per family.
##
## The dashed reference has slope -1/2. A likelihood with n independent observations gives
## posterior SD proportional to n^-1/2, so a parameter tracking the reference is being
## learned from the data at the classical rate; one flattening above it has hit the floor
## the PRIOR imposes, and more data will not help. Reading which parameters do which is the
## substance of the answer, not the overall downward trend, which is uninformative on its own.
##
## The reference is anchored at each parameter's richest rung rather than fitted, so it is a
## statement of what n^-1/2 scaling WOULD look like, not a fit to be judged by its residuals.
##
## SIM_DIR/SIM_TAG default to the joined chains. `n_last` must lie inside the final segment,
## which is the only one that ran a fixed kernel: join_chains.R does not record the segment
## boundaries, and the shortest final segment across the 24 scenarios is about 40,000 draws.
SIM_DIR <- "outputs/final/simulations"
SIM_TAG <- "c"
SIM_SEG_DIR <- "outputs/old"   # the un-joined segments join_chains.R read
SIM_BURNIN  <- 10000L

## How many draws at the end of a joined chain come from ONE fixed proposal.
##
## join_chains.R stitches the segments but records no boundary, and only the final segment
## ran a single kernel throughout -- earlier ones contain adaptations. The boundary is
## recoverable from the segment filenames: <scenario>_<tag>_<iterations>.rds, where the
## highest tag is the last segment and its iteration count is the number of rows it
## contributed (row 1 of a continuation duplicates the source's final state and is dropped).
##
## This matters more than it sounds. A fixed 40,000-draw window was used previously because
## that is the SHORTEST final segment across the 24 scenarios, which threw away half the
## draws for M/N/O, whose final segments run to 90,000-105,000. Three scenarios appeared to
## have ESS below 100 as a result; at their true window they are 114, 216 and 130.
final_segment <- function(scenario, seg_dir = SIM_SEG_DIR) {
  f <- list.files(seg_dir, sprintf("^%s_[0-9]+_[0-9]+\\.rds$", scenario))
  if (!length(f)) return(NA_integer_)
  tg <- as.integer(sub(sprintf("^%s_([0-9]+)_[0-9]+\\.rds$", scenario), "\\1", f))
  it <- as.integer(sub(".*_([0-9]+)\\.rds$", "\\1", f))
  max(it[tg == max(tg)])
}

## The analysis window: the final segment less a burn-in, since a resumed chain starts from
## its source's last state under a NEW proposal and takes a while to settle into it.
## build_segments -- record, for every scenario, how many draws come from one fixed kernel
##
##   build_segments()                    # fill in anything missing, keep what is there
##   build_segments(overwrite = TRUE)    # re-derive everything from the segment filenames
##   build_segments(set = c(R3 = 55000)) # set one by hand
##
## The scan of outputs/old is a fallback, not a foundation: those are working files that may
## be archived or cleared, and if they vanish every scenario silently reverts to the 40,000
## default -- which is the bug that made three scenarios look like they had ESS below 100.
## This writes the numbers down once, in the outputs directory, next to the chains they
## describe. Existing entries are kept unless overwrite = TRUE, so an appended run recorded
## by append_chain.R is never clobbered by a re-scan.
build_segments <- function(scenarios = as.vector(t(outer(c("M","N","O","P","Q","R"), 0:3, paste0))),
                           set = NULL, overwrite = FALSE,
                           dir = SIM_DIR, seg_dir = SIM_SEG_DIR) {
  seg <- if (file.exists(SEG_FILE) && !overwrite) readRDS(SEG_FILE) else integer(0)
  for (s in scenarios) {
    if (!overwrite && s %in% names(seg)) next
    n <- final_segment(s, seg_dir)
    if (is.na(n)) { warning("no segments found for ", s, call. = FALSE); next }
    seg[[s]] <- as.integer(n)
  }
  if (!is.null(set)) for (nm in names(set)) seg[[nm]] <- as.integer(set[[nm]])
  seg <- seg[order(names(seg))]
  dir.create(dirname(SEG_FILE), recursive = TRUE, showWarnings = FALSE)
  saveRDS(seg, SEG_FILE)
  message(sprintf("%d scenarios written to %s", length(seg), SEG_FILE))
  invisible(seg)
}

## What each scenario actually contributes, for checking and for the supplementary.
segments_table <- function(scenarios = NULL, dir = SIM_DIR, tag = SIM_TAG,
                           burnin = SIM_BURNIN) {
  seg <- if (file.exists(SEG_FILE)) readRDS(SEG_FILE) else integer(0)
  if (is.null(scenarios)) scenarios <- names(seg)
  do.call(rbind, lapply(scenarios, function(s) {
    f <- chain_file(s, tag, dir)
    it <- if (is.null(f)) NA_integer_
          else as.integer(sub(".*_([0-9]+)\\.rds$", "\\1", basename(f)))
    data.frame(scenario = s, chain_iterations = it,
               final_segment = if (s %in% names(seg)) seg[[s]] else NA_integer_,
               burnin = burnin, window = sim_window(s), row.names = NULL)
  }))
}

## final_segments.rds, written by append_chain.R, takes precedence over the filename scan.
## It has to: a chain that has had a resumed run appended has a final segment that
## outputs/old cannot know about, since its tags there are numeric and the appended run
## carries a letter. Read fresh each call rather than cached, so appending a run mid-session
## takes effect without re-sourcing.
SEG_FILE <- file.path(SIM_DIR, "final_segments.rds")

sim_window <- function(scenario, seg_dir = SIM_SEG_DIR, burnin = SIM_BURNIN,
                       fallback = 40000L) {
  if (file.exists(SEG_FILE)) {
    s <- readRDS(SEG_FILE)
    ## `s` is a named integer vector, and `s[["absent"]]` throws on an atomic vector rather
    ## than returning NULL as it would for a list -- so test membership first.
    if (scenario %in% names(s) && is.finite(s[[scenario]]))
      return(max(as.integer(s[[scenario]]) - burnin, 5000L))
  }
  n <- final_segment(scenario, seg_dir)
  if (is.na(n)) {
    warning("no segments found for ", scenario, " in ", seg_dir,
            "; falling back to ", format(fallback, big.mark = ","), " draws", call. = FALSE)
    return(fallback)
  }
  max(n - burnin, 5000L)
}

## The internal family letters mean nothing to a reader, so every figure that shows them is
## relabelled by what the scenario actually is. The design is a 2 x 3 grid, and the two
## dimensions get their own aesthetic rather than being flattened into six colours: onset
## regime is the COLOUR and kinetics profile is the LINETYPE. A reader can then follow one
## dimension while holding the other fixed, which is what both comparisons in the Results
## require -- the heterogeneity result compares profiles at matched onset, and the
## identifiability result compares onset regimes at matched profile.
##
## Profile 1 (M/N/O) is the more variable between individuals, a_sigma 1.10 and b_sigma 0.85;
## profile 2 (P/Q/R) the less variable, at 0.75 and 0.70. Dashed for the more variable, since
## it is the one whose posteriors are looser.
## The subject of these labels has to be stated. "mainly pre-peak" is elliptical for "symptom
## onset mainly precedes the peak", but every figure that uses it is aligned on ONSET, so what
## the reader sees is where the PEAK falls -- to the right of zero in the pre-peak scenarios.
## A label reading "pre-peak" beside peaks drawn on the right reads as a contradiction. Naming
## onset explicitly costs two words and removes it.
## S-X are the same scenarios regenerated under incidence indexing (make_scenarios_incidence.R):
## S<-M, T<-N, U<-O, V<-P, W<-Q, X<-R, with the generating parameters unchanged, so each
## inherits its source's onset position and variability profile exactly.
ONSET_LABEL <- c(M = "onset before peak", N = "onset near peak", O = "onset after peak",
                 P = "onset before peak", Q = "onset near peak", R = "onset after peak",
                 S = "onset before peak", T = "onset near peak", U = "onset after peak",
                 V = "onset before peak", W = "onset near peak", X = "onset after peak",
                 A = "onset before peak", B = "onset near peak", C = "onset after peak",
                 D = "onset before peak", E = "onset near peak", F = "onset after peak")
ONSET_LEVELS <- c("onset before peak", "onset near peak", "onset after peak")
ONSET_COL    <- setNames(c("#1B9E77", "#7570B3", "#D95F02"), ONSET_LEVELS)

PROF_LABEL  <- c(M = "1: more variable", N = "1: more variable", O = "1: more variable",
                 P = "2: less variable", Q = "2: less variable", R = "2: less variable",
                 S = "1: more variable", T = "1: more variable", U = "1: more variable",
                 V = "2: less variable", W = "2: less variable", X = "2: less variable",
                 A = "2: less variable", B = "2: less variable", C = "2: less variable",
                 D = "1: more variable", E = "1: more variable", F = "1: more variable")
PROF_LEVELS <- c("1: more variable", "2: less variable")
PROF_LTY    <- setNames(c("22", "solid"), PROF_LEVELS)

onset_factor <- function(f) factor(unname(ONSET_LABEL[f]), levels = ONSET_LEVELS)
prof_factor  <- function(f) factor(unname(PROF_LABEL[f]),  levels = PROF_LEVELS)

## n_last = NULL uses each scenario's own final single-kernel segment, less SIM_BURNIN.
## Pass a number to force one window across all scenarios.
precision_table <- function(families = c("M","N","O","P","Q","R"), rungs = 0:3,
                            pars = REC_PARS, n_last = NULL,
                            dir = SIM_DIR, tag = SIM_TAG) {
  scen <- as.vector(t(outer(families, rungs, paste0)))
  do.call(rbind, lapply(scen, function(s) {
    f <- chain_file(s, tag, dir)
    if (is.null(f)) { message("no chain for ", s, " in ", dir); return(NULL) }
    o  <- readRDS(f)
    fl <- max(which(!is.na(o$MCMC_posteriors)))
    ## The same rule recovery_draws follows: a curated chain records analysis_from, the draw
    ## after the last covariance adaptation that actually happened, and that beats any window
    ## counted back from the end. Matching the two matters here in particular, because
    ## recovery_summary_table pairs this function's post_sd with recovery_table's interval --
    ## if they read different draws, the contraction and the coverage on the same point
    ## describe different posteriors. Read before sim_window() is consulted so a curated
    ## directory with no final_segments.rds does not warn about segments it never had.
    st <- if (!is.null(o$x$analysis_from)) as.integer(o$x$analysis_from)
          else max(1L, fl - (if (is.null(n_last)) sim_window(s) else as.integer(n_last)) + 1L)
    nl <- fl - st + 1L
    dr <- on_natural_scale(o$MCMC_output[st:fl, pars, drop = FALSE], o$x)
    n  <- sum(o$x$data_array)
    tr <- o$x$truth
    ## Prior SD for the same parameters, so precision can be read against what was assumed
    ## rather than only against the other rungs.
    psd <- o$x$prior_sds[pars]
    rm(o); invisible(gc(verbose = FALSE))
    ## ESS travels with the SD because it says whether the SD can be believed. inc1 and inc3
    ## are correlated at -0.96 in the late-onset families -- they trade off almost exactly,
    ## since one sets the threshold and the other the delay from crossing it -- and the
    ## sampler crawls along that ridge. R3's inc3 reaches ess_tail 31 from 40,000 draws, at
    ## which point its posterior SD is too noisy to plot as though it were measured.
    data.frame(scenario = s, family = substr(s, 1, 1), rung = as.integer(substr(s, 2, 2)),
               n_draws = nl, n_tests = n, parameter = pars,
               post_sd = apply(dr, 2, sd),
               width   = apply(dr, 2, function(v) diff(quantile(v, c(0.025, 0.975)))),
               ess     = apply(dr, 2, posterior::ess_tail),
               prior_sd = as.numeric(psd),
               truth = as.numeric(tr[pars]), row.names = NULL)
  }))
}

## `show_prior` defaults to FALSE, and that is a display decision worth recording. The prior
## SD is the natural reference for "how much was learned", but drawing it here destroys the
## figure: rdisp's prior SD is 10 against a posterior SD near 0.001, so that panel spans four
## decades and the sample-size trend -- the entire subject of R1.4c -- flattens to a line.
## figure_contraction carries the comparison against the prior properly, on a bounded scale
## where every parameter is commensurable. This figure answers precision against n, and the
## free y-axes should be set by the posteriors alone.
## rdisp is kept here. This figure plots posterior SD, which makes no reference to the prior, so
## rdisp's deliberately vague prior does not distort it: its SD falls from about 0.035 to 0.003
## across the ladder, the steepest slope of any parameter (-0.34), and that is a real result.
## Ten parameters also fill the 5 x 2 grid exactly, so no panel is left carrying its own x axis.
figure_precision <- function(families = c("M","N","O","P","Q","R"), rungs = 0:3,
                             pars = REC_PARS, n_last = 40000L, ncol = 5,
                             ess_min = 100, show_prior = FALSE,
                             dir = SIM_DIR, tag = SIM_TAG, file = NULL) {

  d <- precision_table(families, rungs, pars, n_last, dir, tag)
  d$parameter <- factor(d$parameter, levels = pars)
  d$scen      <- d$family
  d$onset     <- onset_factor(d$family)
  d$profile   <- prof_factor(d$family)
  ## Hollow points where the SD is not reliably estimated. Drawing them solid would assert a
  ## precision the chain cannot support, and it is exactly these points that produce the
  ## apparent rise of posterior SD with sample size in the late-onset families.
  d$reliable <- d$ess >= ess_min

  ## Reference line per panel: anchored on the richest rung, falling as n^-1/2.
  anch <- do.call(rbind, lapply(split(d, d$parameter), function(z) {
    z <- z[z$n_tests == max(z$n_tests), ][1, ]
    data.frame(parameter = z$parameter, n0 = z$n_tests, sd0 = z$post_sd)
  }))
  rng <- range(d$n_tests)
  ref <- do.call(rbind, lapply(seq_len(nrow(anch)), function(i)
    data.frame(parameter = anch$parameter[i],
               n_tests   = rng,
               post_sd   = anch$sd0[i] * sqrt(anch$n0[i] / rng))))

  ## Where the prior sits. A posterior SD level with it means nothing was learned.
  pri <- unique(d[, c("parameter", "prior_sd")])
  pri <- pri[is.finite(pri$prior_sd) & pri$prior_sd > 0, ]
  if (!show_prior) pri <- pri[0, ]

  p <- ggplot(d, aes(n_tests, post_sd)) +
    geom_line(data = ref, aes(group = parameter), linetype = 2,
              colour = "grey45", linewidth = 0.4) +
    { if (nrow(pri)) geom_hline(data = pri, aes(yintercept = prior_sd),
                                colour = "steelblue", linetype = 3, linewidth = 0.4) } +
    geom_line(aes(colour = onset, linetype = profile, group = scen), linewidth = 0.7) +
    geom_point(aes(colour = onset, shape = reliable), size = 1.6, fill = "white") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
    scale_colour_manual(values = ONSET_COL, name = "Symptom onset") +
    scale_linetype_manual(values = PROF_LTY, name = "Kinetics profile") +
    scale_x_log10(breaks = c(1e3, 1e4, 1e5, 1e6),
                  labels = c("1k", "10k", "100k", "1M")) +
    scale_y_log10() +
    facet_wrap(~parameter, scales = "free_y", ncol = ncol, labeller = rec_labeller()) +
    labs(x = "positive tests", y = "posterior SD", colour = NULL,
         title = "Posterior precision against sample size",
         subtitle = sprintf(paste("dashed grey = n^-1/2 reference anchored at the richest",
                                  "rung%s;\nhollow points have ESS < %d and their SD is not",
                                  "reliably estimated"),
                            if (show_prior) "; dotted blue = prior SD" else "", ess_min)) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold"),
          axis.text = element_text(colour = "black"),
          legend.position = "top")

  if (!is.null(file))
    ggsave(file, p, width = 10, height = 1.9 * ceiling(length(pars) / ncol) + 1.2)
  print(p)
  invisible(list(plot = p, table = d))
}

## ---------------------------------------------------------------------------
## figure_recovery_summary -- the whole simulation study in one panel (E1, main text)
##
##   figure_recovery_summary()
##   figure_recovery_summary(pars = c("a_bar","b_bar","vl_max_bar","inc1","b_sigma"))
##
## Six families x four sample sizes x ten parameters is 240 recovery checks, which no
## interval plot can show at once. This plots instead WHERE THE TRUTH SITS INSIDE ITS OWN
## 95% interval, rescaled so 0 is the lower bound and 1 the upper:
##
##     position = (truth - lo) / (hi - lo)
##
## Inside [0, 1] the interval covers the truth; outside it misses, and how far outside says
## how badly. 0.5 is dead centre. The quantity is dimensionless, so every parameter shares
## one axis and the panels are directly comparable -- which the raw intervals are not, since
## they differ by orders of magnitude in width.
##
## Read it as: do points sit inside the band, and do they tighten toward 0.5 as sample size
## grows. Both questions at once, which is what E1 asks.
recovery_summary_table <- function(families = c("M","N","O","P","Q","R"), rungs = 0:3,
                                   pars = REC_PARS, n_last = NULL,
                                   dir = SIM_DIR, tag = SIM_TAG) {
  scen <- as.vector(t(outer(families, rungs, paste0)))
  d <- recovery_table(scen, pars, n_last, dir, tag)
  d$family   <- substr(d$scenario, 1, 1)
  d$position <- (d$truth - d$lo) / (d$hi - d$lo)

  ## Precision, and how much of it came from the data rather than the prior. Contraction is
  ## 1 - posterior variance / prior variance: 1 means the data determined the parameter, 0
  ## means the posterior is the prior and nothing was learned. Pairing it with z is what
  ## stops a vague posterior from being scored as a success -- covering the truth with an
  ## interval as wide as the prior is not recovery.
  pr <- precision_table(families, rungs, pars, n_last, dir, tag)
  k  <- match(paste(d$scenario, d$parameter), paste(pr$scenario, pr$parameter))
  d$post_sd     <- pr$post_sd[k]
  d$prior_sd    <- pr$prior_sd[k]
  d$ess         <- pr$ess[k]
  d$contraction <- 1 - (d$post_sd / d$prior_sd)^2
  d$z           <- (d$mid - d$truth) / d$post_sd
  d
}

## The standardised error z = (posterior mean - truth) / posterior SD against sample size.
##
## Chosen over plotting where the truth sits inside its interval, which was the first
## attempt. That quantity throws the WIDTH away, so a parameter whose interval is wide
## because the prior is carrying it -- b_bar and b_sigma, precisely the two known to be
## weakly identified -- scores best. It also has no useful reference value: under a
## calibrated posterior the truth is UNIFORM within its own interval at every sample size,
## so there is nothing for the points to converge to.
##
## z has both properties the other lacks. It is dimensionless, so every parameter shares an
## axis; and a confidently-wrong posterior is penalised rather than rewarded, because
## shrinking the SD without moving toward the truth pushes z out. Under a well-calibrated
## fit z is approximately standard normal, so the +-1.96 band is where 95% of points belong
## and anything beyond it is the interval missing.
figure_recovery_summary <- function(families = c("M","N","O","P","Q","R"), rungs = 0:3,
                                    pars = REC_PARS, n_last = 40000L, ncol = 5,
                                    ess_min = 100, dir = SIM_DIR, tag = SIM_TAG, file = NULL) {

  d <- recovery_summary_table(families, rungs, pars, n_last, dir, tag)
  d$parameter <- factor(d$parameter, levels = pars)
  d$scen      <- d$family
  d$onset     <- onset_factor(d$family)
  d$profile   <- prof_factor(d$family)
  d$reliable  <- d$ess >= ess_min
  lim <- c(min(-4, quantile(d$z, 0.005)), max(4, quantile(d$z, 0.995)))
  d$shown <- pmin(pmax(d$z, lim[1]), lim[2])

  n_miss <- sum(!d$covered)
  p <- ggplot(d, aes(n_tests, shown)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.96, ymax = 1.96,
             fill = "grey85", alpha = 0.6) +
    geom_hline(yintercept = 0, colour = "grey40", linetype = 2, linewidth = 0.3) +
    geom_line(aes(colour = onset, linetype = profile, group = scen),
              linewidth = 0.6, alpha = 0.8) +
    geom_point(aes(colour = onset, shape = reliable), size = 1.7, fill = "white") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
    scale_colour_manual(values = ONSET_COL, name = "Symptom onset") +
    scale_linetype_manual(values = PROF_LTY, name = "Kinetics profile") +
    scale_x_log10(breaks = c(1e3, 1e4, 1e5, 1e6), labels = c("1k", "10k", "100k", "1M")) +
    coord_cartesian(ylim = lim) +
    facet_wrap(~parameter, scales = "free_x", ncol = ncol, labeller = rec_labeller()) +
    labs(x = "positive tests", y = "(posterior mean - truth) / posterior SD", colour = NULL,
         title = "Recovery of known parameter values across 24 simulated datasets",
         subtitle = sprintf(paste("shaded band = +-1.96 SD, where 95%% of points belong;",
                                  "hollow points have ESS < %d\n%d of %d 95%% intervals",
                                  "contain the generating value"),
                            ess_min, nrow(d) - n_miss, nrow(d))) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold"),
          axis.text = element_text(colour = "black"),
          legend.position = "top")

  if (!is.null(file))
    ggsave(file, p, width = 10, height = 1.9 * ceiling(length(pars) / ncol) + 1.4)
  print(p)
  invisible(list(plot = p, table = d))
}

## figure_contraction -- the whole study as one panel: did we learn it, and did we get it right
##
## x is posterior contraction, 1 - post_var/prior_var: 0 means the posterior is the prior and
## the data said nothing, 1 means the data determined the parameter. y is the standardised
## error. Every one of the 240 checks is a point.
##
## y_axis = "z" is that original figure. The two axes then separate failures that a coverage
## count conflates: a point at low contraction and z near 0 has not been recovered but left at
## the prior, which happens to sit near the truth, while a point at high contraction and large
## |z| is confidently wrong.
##
## y_axis = "onset" drops z and puts the symptom-onset regime on y instead. Two reasons to
## prefer it, and they are the reasons this is now the default:
##
##   * z is already shown, better, by figure_recovery_grid(): an interval drawn against its
##     generating value says everything z says and also shows the width and the direction of
##     the miss. Spending an axis on it here duplicates that figure and dilutes this one.
##   * onset regime was encoded NOWHERE in the z version. The study has three design
##     dimensions -- sample size, onset regime, variability profile -- and colour and shape
##     carried only two of them, so the scenarios that exist to test whether late onset damages
##     identifiability could not be read off the figure at all.
##
## What survives is the one thing no other figure carries. figure_recovery_grid() no longer
## draws the prior, and figure_precision() plots posterior SD with no prior reference, so
## contraction is the only place in the paper that says how much of the posterior's precision
## came from the data rather than from the longitudinal prior. That is the answer to R1.4b, and
## a figure that answers it alone is stronger than one that half-answers it alongside a
## duplicate of R1.4c.
## rdisp is drawn, and its panel needs a caveat the others do not. Its prior is N(1, 10) with no
## lower bound -- deliberately vague, and in fact placing about 46% of its mass on negative
## values the likelihood then excludes -- against a posterior SD near 0.003. Contraction is
## therefore 1.000 in all 24 cells, which measures the width of that prior rather than anything
## the data did. Kept because ten parameters fill the 5 x 2 grid exactly and because a reader
## comparing this figure with the recovery grid should find the same parameter set in both, but
## the caption has to say what its column means. Pass pars = setdiff(REC_PARS, "rdisp") to drop
## it.
figure_contraction <- function(families = c("M","N","O","P","Q","R"), rungs = 0:3,
                               pars = REC_PARS,
                               n_last = NULL, ess_min = 100, ncol = 5,
                               y_axis = c("onset", "z"),
                               dir = SIM_DIR, tag = SIM_TAG, file = NULL) {
  y_axis <- match.arg(y_axis)

  d <- recovery_summary_table(families, rungs, pars, n_last, dir, tag)
  d <- d[is.finite(d$contraction) & is.finite(d$z), ]
  d$parameter <- factor(d$parameter, levels = pars)
  d$reliable  <- d$ess >= ess_min
  d$profile   <- prof_factor(d$family)
  d$onset     <- onset_factor(d$family)

  ## Kinetics profile is the SHAPE here -- triangles for the more variable profile, circles
  ## for the less variable -- because this figure's main finding is the difference between
  ## them: b_sigma sits at 0.27 contraction under one and 0.60 under the other. Colour is
  ## already carrying sample size and linetype has nothing to attach to on a scatter.
  ##
  ## Reliability rides on the same scale rather than a second aesthetic: the unreliable
  ## variants are the hollow forms of each shape, giving four levels whose labels say what
  ## they are. Only nine of 240 points are hollow, all in one scenario family.
  d$mark <- factor(
    ifelse(d$reliable, as.character(d$profile),
           paste0(as.character(d$profile), ", ESS < ", ess_min)),
    levels = c(PROF_LEVELS[1], paste0(PROF_LEVELS[1], ", ESS < ", ess_min),
               PROF_LEVELS[2], paste0(PROF_LEVELS[2], ", ESS < ", ess_min)))
  mark_shape <- setNames(c(17, 2, 16, 1), levels(d$mark))
  ## Sample size as a four-level factor rather than a size aesthetic: within a panel the
  ## story is the left-to-right march of contraction as data accumulates, and four labelled
  ## colours read that far better than a continuous area scale.
  d$size_lab <- factor(c("1k", "10k", "100k", "1M")[as.integer(substr(d$scenario, 2, 2)) + 1L],
                       levels = c("1k", "10k", "100k", "1M"))
  lim <- c(min(-4, quantile(d$z, 0.005)), max(4, quantile(d$z, 0.995)))
  d$shown <- pmin(pmax(d$z, lim[1]), lim[2])

  ## The x scale is deliberately NOT free, in either variant. Contraction is already
  ## dimensionless and on a common 0-1 scale, and comparing panels against each other is the
  ## entire point: b_sigma sitting at 0.25 while inc1 sits at 0.98 is the finding.
  ##
  ## The lower limit runs slightly below zero on purpose. b_sigma comes out at -0.015 and
  ## -0.038 in the two smallest less-variable datasets: posterior SD marginally ABOVE the prior
  ## SD, with ESS above 1,600, so this is not a mixing artefact -- the posterior IS the prior
  ## and the data added nothing measurable. Clipping at zero would hide the clearest instance
  ## of the paper's own limitation.
  x_common <- list(
    scale_x_continuous(limits = c(-0.05, 1), breaks = c(0, 0.5, 1)),
    scale_colour_viridis_d(option = "D", end = 0.9, name = "Positive tests"),
    facet_wrap(~parameter, ncol = ncol, labeller = rec_labeller()),
    theme_bw(base_size = 10),
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold"),
          axis.text = element_text(colour = "black"),
          legend.position = "top"))

  if (y_axis == "onset") {
    ## Onset reversed so the three regimes read top to bottom in the order the trajectory
    ## visits them -- onset before the peak at the top -- which matches every other figure that
    ## lists them and the column order of the recovery grid.
    d$onset <- factor(d$onset, levels = rev(levels(d$onset)))
    ## "onset before peak" three times down a ten-panel figure spends about an inch of width on
    ## a word the axis title already says. Stripped to the distinguishing part, with the
    ## subject moved to the axis title, which is where a reader looks for it once.
    onset_short <- setNames(sub("^onset ", "", levels(d$onset)), levels(d$onset))
    ## Profiles are dodged within each onset row rather than left to overplot. At 1M tests most
    ## parameters sit between 0.95 and 1.00, so eight points per row would collapse onto each
    ## other and the profile contrast -- the thing worth seeing in the b_sigma panel -- would
    ## be invisible.
    p <- ggplot(d, aes(contraction, onset)) +
      geom_vline(xintercept = 0, colour = "grey40", linetype = 2, linewidth = 0.3) +
      geom_point(aes(colour = size_lab, shape = profile), size = 1.9, alpha = 0.9,
                 position = position_dodge(width = 0.7)) +
      scale_shape_manual(values = c(17, 16), name = "Kinetics profile") +
      scale_y_discrete(labels = onset_short) +
      x_common +
      labs(x = "Posterior contraction", y = "Symptom onset",
           title = "How much the data determined each parameter",
           subtitle = paste("each point is one parameter in one simulated dataset;",
                            "0 = posterior no narrower than the prior")) +
      theme(axis.text.y = element_text(colour = "black", size = 8))
    h <- 1.55 * ceiling(length(pars) / ncol) + 1.4

  } else {
    p <- ggplot(d, aes(contraction, shown)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.96, ymax = 1.96,
               fill = "grey88", alpha = 0.6) +
      geom_hline(yintercept = 0, colour = "grey40", linetype = 2, linewidth = 0.3) +
      geom_point(aes(colour = size_lab, shape = mark), size = 1.9, alpha = 0.9) +
      scale_shape_manual(values = mark_shape, name = "Kinetics profile") +
      coord_cartesian(ylim = lim) +
      x_common +
      ## Bare axis names, with both definitions carried by the figure caption instead.
      labs(x = "Posterior contraction",
           y = "Posterior Z-score",
           title = "What was learned, and whether it was right",
           subtitle = paste("each point is one parameter in one simulated dataset;",
                            "shaded band = +-1.96 posterior SD, where 95% of points belong"))
    h <- 1.9 * ceiling(length(pars) / ncol) + 1.4
  }

  if (!is.null(file)) ggsave(file, p, width = 10, height = h)
  print(p)
  invisible(list(plot = p, table = d))
}

## ---------------------------------------------------------------------------
## recovery_ess -- effective sample size for the simulation-study scenarios
##
##   recovery_ess("M")                          # M0-M3, current segment, whole chain
##   recovery_ess(c("M","N","O"))               # all twelve
##   recovery_ess("M", n_last = 20000)          # the window an interval is quoted from
##
## Why this and not check_chains(): these scenarios have ONE chain each -- run in segments
## (_1 -> _2 -> _3) and joined by join_chains.R into the tag "c" files read here.
## check_chains would read segments as parallel chains and report a split-Rhat between
## draws that are sequential by construction, which is not a convergence diagnostic at
## all. With one chain there is no Rhat to compute, so ESS is the whole diagnostic.
##
## posterior::ess_bulk is rank-normalised and split -- it halves the chain internally and so
## still detects drift within a single run. ess_bulk is the sample size behind a posterior
## mean, ess_tail behind the 2.5%/97.5% edges. Recovery is judged on intervals, so ess_tail
## is the binding number here; below ~100 per chain the interval edges are not resolved and
## a "miss" may be Monte Carlo noise rather than a real failure.
## n_last = NULL uses each scenario's own analysis window from sim_window(), so these numbers
## describe exactly the draws the figures and tables summarise. It used to fall back to
## burnin(), which finds the first iteration reaching the chain's median log posterior --
## right for a single-kernel run, wrong here, because it reaches back across resume points
## into segments that used a different proposal.
recovery_ess <- function(families = c("M", "N", "O", "P", "Q", "R"), rungs = 0:3,
                         pars = REC_PARS, n_last = NULL,
                         dir = "outputs/final/simulations", tag = "c") {

  scen <- if (all(nchar(families) == 1)) as.vector(t(outer(families, rungs, paste0))) else families

  out <- do.call(rbind, lapply(scen, function(s) {
    f <- chain_file(s, tag, dir)
    if (is.null(f)) { message("no segment '", tag, "' for ", s, " in ", dir); return(NULL) }
    nl <- if (is.null(n_last)) sim_window(s) else as.integer(n_last)
    o  <- readRDS(f)
    fl <- max(which(!is.na(o$MCMC_posteriors)))
    st <- max(1L, fl - nl + 1L)
    dr <- on_natural_scale(o$MCMC_output[st:fl, intersect(pars, colnames(o$MCMC_output)), drop = FALSE], o$x)
    lp <- o$MCMC_posteriors[st:fl, 1]
    acc <- mean(o$Acceptances[st:fl], na.rm = TRUE)
    rm(o); invisible(gc(verbose = FALSE))

    e <- do.call(rbind, lapply(colnames(dr), function(p)
      data.frame(scenario = s, parameter = p,
                 bulk = round(posterior::ess_bulk(dr[, p])),
                 tail = round(posterior::ess_tail(dr[, p])), row.names = NULL)))
    rbind(data.frame(scenario = s, parameter = "log posterior",
                     bulk = round(posterior::ess_bulk(lp)),
                     tail = round(posterior::ess_tail(lp))), e)
  }))
  if (is.null(out)) stop("no chains found")

  cat(sprintf("chain _%s | window: %s\n\n", tag,
              if (is.null(n_last)) "each scenario's own final segment, see segments_table()"
              else sprintf("%s draws, fixed", format(n_last, big.mark = ","))))
  ## one row per scenario, one column per parameter -- the shape you actually scan
  w <- reshape(out[, c("scenario", "parameter", "tail")], direction = "wide",
               idvar = "scenario", timevar = "parameter")
  names(w) <- sub("^tail\\.", "", names(w))
  w$draws <- vapply(w$scenario, function(s) sim_window(s), integer(1))
  w <- w[, c("scenario", "draws", setdiff(names(w), c("scenario", "draws")))]
  print(w, row.names = FALSE)

  lo <- out[out$tail < 100, ]
  cat(sprintf("\nESS_tail below 100: %d of %d\n", nrow(lo), nrow(out)))
  if (nrow(lo)) print(lo[order(lo$tail), ], row.names = FALSE)
  cat("\n(ess_tail above; ess_bulk in the returned data frame)\n")
  invisible(out)
}
