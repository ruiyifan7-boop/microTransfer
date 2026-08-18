\name{prediction}
\alias{balanced_accuracy}
\alias{macro_f1}
\alias{fit_nearest_centroid}
\alias{predict_nearest_centroid}
\alias{validate_prediction}
\alias{prediction_permutation_test}
\alias{label_noise_sensitivity}
\title{Held-Out-Cohort Prediction and Label-Noise Diagnostics}
\usage{
balanced_accuracy(truth, predicted)

macro_f1(truth, predicted)

fit_nearest_centroid(x, y, ...)

predict_nearest_centroid(model, newx, ...)

validate_prediction(
  counts, metadata, sample_col, study_col, outcome_col,
  features = NULL, transform = "clr", pseudocount = 0.5,
  min_prevalence = 0.1, min_total = 20,
  fit_fun = fit_nearest_centroid,
  predict_fun = predict_nearest_centroid, ...
)

prediction_permutation_test(
  counts, metadata, sample_col, study_col, outcome_col,
  features = NULL, transform = "clr", pseudocount = 0.5,
  min_prevalence = 0.1, min_total = 20,
  fit_fun = fit_nearest_centroid,
  predict_fun = predict_nearest_centroid,
  permutations = 999L, seed = 1L, ...
)

label_noise_sensitivity(
  counts, metadata, sample_col, study_col, outcome_col,
  features = NULL, transform = "clr", pseudocount = 0.5,
  min_prevalence = 0.1, min_total = 20,
  fit_fun = fit_nearest_centroid,
  predict_fun = predict_nearest_centroid,
  noise = c(0, 0.05, 0.1, 0.15, 0.2),
  replicates = 100L, seed = 1L, ...
)
}
\arguments{
  \item{truth}{Observed labels.}
  \item{predicted}{Predicted labels.}
  \item{x}{Training matrix.}
  \item{y}{Training labels.}
  \item{model}{Nearest-centroid model.}
  \item{newx}{Numeric test matrix.}
  \item{counts}{Sample-by-feature counts.}
  \item{metadata}{Sample metadata.}
  \item{sample_col}{Sample identifier column.}
  \item{study_col}{Cohort column.}
  \item{outcome_col}{Categorical outcome column.}
  \item{features}{Optional prespecified features.}
  \item{transform}{Feature transformation.}
  \item{pseudocount}{CLR pseudocount.}
  \item{min_prevalence}{Training prevalence filter.}
  \item{min_total}{Training total-count filter.}
  \item{fit_fun}{Function accepting \code{x} and \code{y}.}
  \item{predict_fun}{Function accepting a model and \code{newx}.}
  \item{noise}{Target label-noise fractions.}
  \item{replicates}{Replicates per nonzero noise level.}
  \item{permutations}{Number of within-cohort label permutations.}
  \item{seed}{Random seed.}
  \item{...}{Additional model arguments.}
}
\value{Performance metrics, a model, or a prediction validation object.}
\description{
Evaluates categorical prediction in held-out cohorts and quantifies
sensitivity to class-count-preserving training-label perturbations.
}
