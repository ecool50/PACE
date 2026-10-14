// pybind11 bindings for the PACE core.
//
// The core is R-free by construction, so the same translation units back both
// the R package and this module. Arrays cross the boundary as numpy arrays in
// Fortran order, which is what the core expects and what R already uses, so
// neither side copies.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "decomposition.hpp"
#include "neighbourhood.hpp"
#include "count_stats.hpp"
#include "hyperparameters.hpp"
#include "irls_driver.hpp"
#include "preprocess.hpp"
#include "sparse_product.hpp"
#include "random_design.hpp"
#include "special.hpp"
#include "mash.hpp"

namespace py = pybind11;

namespace {

using ColumnMajor = py::array_t<double, py::array::f_style | py::array::forcecast>;

void raise_if_failed(const pace::Status& status, const char* what) {
  if (status.is_ok()) return;
  throw std::runtime_error(std::string(what) + ": " + status.message);
}

pace::Span<const double> span_of(const ColumnMajor& array) {
  return pace::Span<const double>(array.data(), static_cast<std::int64_t>(array.size()));
}

// Rcpp's NumericMatrix arrives zero-filled and parts of the core rely on that,
// writing only the entries they compute. pybind11 does NOT zero a new array, so
// zero it here -- otherwise untouched entries are whatever was on the heap.
pace::Span<double> out_span_of(ColumnMajor& array) {
  double* data = array.mutable_data();
  std::fill(data, data + array.size(), 0.0);
  return pace::Span<double>(data, static_cast<std::int64_t>(array.size()));
}

// Canonical mash for one slice, end to end. `bhat` and `shat` are
// n_effects x n_conditions; `v` may be empty for the identity, which is what
// PACE publishes with.
py::dict mash(const ColumnMajor& bhat, const ColumnMajor& shat, const ColumnMajor& v,
              double null_weight, double grid_mult, std::int64_t max_iterations,
              double tolerance, const std::vector<std::string>& condition_names) {
  if (bhat.ndim() != 2 || shat.ndim() != 2) throw std::runtime_error("bhat and shat must be 2-D");
  if (bhat.shape(0) != shat.shape(0) || bhat.shape(1) != shat.shape(1)) {
    throw std::runtime_error("bhat and shat must have the same shape");
  }
  const std::int64_t n_effects = bhat.shape(0);
  const std::int64_t n_conditions = bhat.shape(1);

  const std::vector<pace::CovarianceComponent> covariances =
      pace::canonical_covariances(n_conditions, condition_names);

  std::vector<double> grid;
  raise_if_failed(pace::autoselect_grid(span_of(bhat), span_of(shat), n_effects, n_conditions,
                                        grid_mult, &grid),
                  "grid");

  const std::int64_t n_components = pace::n_mixture_components(
      static_cast<std::int64_t>(covariances.size()), static_cast<std::int64_t>(grid.size()), true);

  ColumnMajor loglik({n_effects, n_components});
  raise_if_failed(pace::mixture_loglik(span_of(bhat), span_of(shat), span_of(v), n_effects,
                                       n_conditions, covariances, grid, true, 1,
                                       out_span_of(loglik)),
                  "likelihood");

  std::vector<double> weights;
  std::int64_t iterations = 0;
  raise_if_failed(pace::optimize_mixture_weights(span_of(loglik), n_effects, n_components,
                                                 null_weight, max_iterations, tolerance, &weights,
                                                 &iterations),
                  "mixture weights");

  ColumnMajor posterior_mean({n_effects, n_conditions});
  ColumnMajor posterior_sd({n_effects, n_conditions});
  ColumnMajor lfsr({n_effects, n_conditions});
  raise_if_failed(pace::posterior_summaries(span_of(bhat), span_of(shat), span_of(v), n_effects,
                                            n_conditions, covariances, grid, true, weights,
                                            span_of(loglik), out_span_of(posterior_mean),
                                            out_span_of(posterior_sd), out_span_of(lfsr)),
                  "posteriors");

  py::dict out;
  out["posterior_mean"] = posterior_mean;
  out["posterior_sd"] = posterior_sd;
  out["lfsr"] = lfsr;
  out["pi"] = weights;
  out["grid"] = grid;
  out["iterations"] = iterations;
  return out;
}

// Driver scores for one focal-neighbour pair. Pure arithmetic over per-gene
// vectors, so the same call backs `topDrivers()` in R. `u_raw` is only read on
// a condition fit, where V_S comes from the unshrunken baseline BLUP and V_RxS
// from the shrunken interaction; pass it empty otherwise.
py::dict driver_scores(const ColumnMajor& estimate_shrunk, const ColumnMajor& spec,
                       const ColumnMajor& focal_mean, const ColumnMajor& mu_bar,
                       const ColumnMajor& alpha, const ColumnMajor& u_raw, double var_n,
                       double var_rn, bool has_responder) {
  const std::int64_t n = estimate_shrunk.size();
  for (const auto* a : {&spec, &focal_mean, &mu_bar, &alpha}) {
    if (a->size() != n) throw std::runtime_error("every per-gene input must have the same length");
  }
  if (has_responder && u_raw.size() != n) {
    throw std::runtime_error("a condition fit needs u_raw, one value per gene");
  }
  ColumnMajor score(n), score4(n), v_resid(n), v_s(n), v_total(n), r2_s(n);
  ColumnMajor v_rxs(has_responder ? n : 0), r2_rxs(has_responder ? n : 0);
  raise_if_failed(pace::driver_scores(span_of(estimate_shrunk), span_of(spec), span_of(focal_mean),
                                      span_of(mu_bar), span_of(alpha), span_of(u_raw), var_n,
                                      var_rn, has_responder, n, out_span_of(score),
                                      out_span_of(score4), out_span_of(v_resid), out_span_of(v_s),
                                      out_span_of(v_rxs), out_span_of(v_total), out_span_of(r2_s),
                                      out_span_of(r2_rxs)),
                  "driver scores");
  py::dict out;
  out["driver_score"] = score;
  out["driver_score_spec4"] = score4;
  out["V_resid"] = v_resid;
  out["V_S"] = v_s;
  out["V_RxS"] = v_rxs;
  out["V_total"] = v_total;
  out["R2_S"] = r2_s;
  out["R2_RxS"] = r2_rxs;
  return out;
}

// Pratt attribution of a focal type's spatial variance across its neighbours.
// `sigma` is the n_pairs x n_pairs covariance of the neighbour kernels and `u`
// the n_pairs x n_genes slopes. This is what the pair heatmap reads.
py::dict pair_variance_pratt(const ColumnMajor& sigma, const ColumnMajor& u) {
  if (sigma.ndim() != 2 || sigma.shape(0) != sigma.shape(1)) {
    throw std::runtime_error("sigma must be a square n_pairs x n_pairs matrix");
  }
  if (u.ndim() != 2) throw std::runtime_error("u must be n_pairs x n_genes");
  const std::int64_t n_pairs = sigma.shape(0);
  if (u.shape(0) != n_pairs) throw std::runtime_error("u must have one row per pair");
  const std::int64_t n_genes = u.shape(1);

  ColumnMajor v_pair(n_pairs), v_diag(n_pairs), share_pct(n_pairs), diag_pct(n_pairs), sign(n_pairs);
  double v_total = 0.0, v_diag_total = 0.0, cross_cov_pct = 0.0;
  int n_negative = 0;
  raise_if_failed(pace::pair_variance_pratt(span_of(sigma), span_of(u), n_pairs, n_genes,
                                            out_span_of(v_pair), out_span_of(v_diag),
                                            out_span_of(share_pct), out_span_of(diag_pct),
                                            out_span_of(sign), &v_total, &v_diag_total,
                                            &cross_cov_pct, &n_negative),
                  "pair variance");
  py::dict out;
  out["v_pair"] = v_pair;
  out["v_diag"] = v_diag;
  out["share_pct"] = share_pct;
  out["diag_pct"] = diag_pct;
  out["sign"] = sign;
  out["v_total"] = v_total;
  out["v_diag_total"] = v_diag_total;
  out["cross_cov_pct"] = cross_cov_pct;
  out["n_negative"] = n_negative;
  return out;
}

using IntColumn = py::array_t<int, py::array::c_style | py::array::forcecast>;
// Column-major integers, for the 2-D index matrices the core reads the way R's
// IntegerMatrix is laid out. IntColumn above is c_style, so handing it a
// Fortran array repacks the memory and the core then reads the matrix
// TRANSPOSED -- which for a square one raises nothing at all and simply picks
// the wrong rows.
using IntMatrix = py::array_t<int, py::array::f_style | py::array::forcecast>;

// Gaussian biological and exponential technical neighbour kernels, the first
// stage of a fit. For each cell, the kernel-weighted count of every other cell
// type within the radius.
//
// `cell_type` is 0-based type codes, one per cell. `group` assigns cells to
// images; pass `per_group = true` to build the kernels within each image
// separately, which is what a cohort of separate samples needs. `eps` is the
// radius, defaulting to 3 * h_bio as the method specifies.
py::dict neighbour_kernels(const ColumnMajor& x, const ColumnMajor& y,
                           const IntColumn& cell_type, int n_types, const IntColumn& group,
                           bool per_group, double h_bio, double h_tech, double eps,
                           int n_threads) {
  const std::int64_t n_cells = x.size();
  if (y.size() != n_cells || cell_type.size() != n_cells) {
    throw std::runtime_error("x, y and cell_type must have one entry per cell");
  }
  if (group.size() != n_cells) throw std::runtime_error("group must have one entry per cell");
  if (eps <= 0.0) eps = 3.0 * h_bio;

  ColumnMajor k_bio({n_cells, static_cast<std::int64_t>(n_types)});
  ColumnMajor k_tech({n_cells, static_cast<std::int64_t>(n_types)});
  const pace::Status status = pace::neighbour_kernels(
      span_of(x), span_of(y),
      pace::Span<const int>(cell_type.data(), n_cells), n_types,
      pace::Span<const int>(group.data(), n_cells), per_group,
      h_bio, h_tech, eps, out_span_of(k_bio), out_span_of(k_tech),
      n_threads, pace::InterruptCheck());
  raise_if_failed(status, "neighbour kernels");

  py::dict out;
  out["k_bio"] = k_bio;
  out["k_tech"] = k_tech;
  out["eps"] = eps;
  return out;
}

// Sparse cross-type ambient weight matrix W, the contamination field:
//
//   W[i, j] = exp(-d_ij / h_tech) / edge_fraction[i]
//
// over cells i, j of the same image, different type, within 3 * h_tech.
// Returned in compressed-column form, so
//
//   scipy.sparse.csc_matrix((data, indices, indptr), shape=shape)
//
// rebuilds it without this module depending on scipy.
//
// The edge correction integrates over a quadrature of `n_angles` directions,
// which R fixes at 1000; it is built here rather than asked of the caller.
py::dict ambient_field(const ColumnMajor& x, const ColumnMajor& y, const IntColumn& cell_type,
                       const IntColumn& image, int n_images, double h_tech, bool edge_correct,
                       int n_angles, int n_threads) {
  const std::int64_t n = x.size();
  if (y.size() != n || cell_type.size() != n || image.size() != n) {
    throw std::runtime_error("x, y, cell_type and image must have one entry per cell");
  }
  if (n_angles < 1) throw std::runtime_error("n_angles must be positive");

  // The angles R hands the core, reproduced bit for bit. R writes them as
  //   seq(0, 2 * pi, length.out = n_angles + 1)[-1]
  // and seq() divides ONCE and then multiplies, so theta_k = k * (2 pi / n),
  // with the last element set to 2 pi exactly rather than computed. Writing it
  // the obvious way instead, 2 pi * k / n, moves 352 of the 1000 angles by one
  // ulp, which carries through the edge correction into the ambient field and
  // on into the fit: it was worth 8e-07 in the random effects.
  std::vector<double> cos_theta(static_cast<std::size_t>(n_angles));
  std::vector<double> sin_theta(static_cast<std::size_t>(n_angles));
  const double two_pi = 2.0 * 3.14159265358979323846;
  const double step = two_pi / static_cast<double>(n_angles);
  for (int k = 0; k < n_angles; ++k) {
    const double theta = k + 1 == n_angles ? two_pi : static_cast<double>(k + 1) * step;
    cos_theta[static_cast<std::size_t>(k)] = std::cos(theta);
    sin_theta[static_cast<std::size_t>(k)] = std::sin(theta);
  }

  ColumnMajor edge_fraction(n);
  pace::AmbientFieldBuilder builder;
  raise_if_failed(builder.prepare(span_of(x), span_of(y),
                                  pace::Span<const int>(cell_type.data(), n),
                                  pace::Span<const int>(image.data(), n), n_images, h_tech,
                                  edge_correct,
                                  pace::Span<const double>(cos_theta.data(), n_angles),
                                  pace::Span<const double>(sin_theta.data(), n_angles),
                                  out_span_of(edge_fraction), n_threads, pace::InterruptCheck()),
                  "ambient field");

  const std::int64_t nnz = builder.nnz();
  py::array_t<int> indices(nnz);
  std::fill(indices.mutable_data(), indices.mutable_data() + nnz, 0);
  ColumnMajor data(nnz);
  // int32, because pace::CscView takes Span<const int> and these arrays
  // are meant to feed straight back into the core.
  py::array_t<int> indptr(n + 1);
  const std::vector<std::int64_t>& pointer = builder.column_pointer();
  auto indptr_out = indptr.mutable_unchecked<1>();
  for (std::int64_t j = 0; j <= n; ++j) {
    indptr_out(j) = static_cast<int>(pointer[static_cast<std::size_t>(j)]);
  }

  raise_if_failed(builder.fill(pace::Span<int>(indices.mutable_data(), nnz), out_span_of(data),
                               pace::Span<const double>(edge_fraction.data(), n), n_threads,
                               pace::InterruptCheck()),
                  "ambient field");

  py::dict out;
  out["data"] = data;
  out["indices"] = indices;
  out["indptr"] = indptr;
  out["shape"] = py::make_tuple(n, n);
  out["edge_fraction"] = edge_fraction;
  out["nnz"] = nnz;
  return out;
}

// Zero a freshly allocated integer array, for the same reason out_span_of()
// zeroes a double one.
pace::Span<int> out_int_span(py::array_t<int>& array) {
  int* data = array.mutable_data();
  std::fill(data, data + array.size(), 0);
  return pace::Span<int>(data, static_cast<std::int64_t>(array.size()));
}

// One block of the random-effect design: the cells-by-terms matrix spread
// across groups, so column (group, term) carries that group's cells only.
// Returned in compressed-column form, with the grouping the core computed on
// the same pass -- `cells_by_group` holds cell indices ordered by group and
// `group_start` indexes into it, so group g owns
// cells_by_group[group_start[g]:group_start[g + 1]].
py::dict random_design_block(const ColumnMajor& x_terms, const IntColumn& cell_group,
                             int n_groups) {
  if (x_terms.ndim() != 2) throw std::runtime_error("x_terms must be n_cells x n_terms");
  const std::int64_t n = x_terms.shape(0);
  const std::int64_t n_terms = x_terms.shape(1);
  if (cell_group.size() != n) throw std::runtime_error("cell_group must have one entry per cell");
  if (n_groups < 1) throw std::runtime_error("n_groups must be positive");
  const std::int64_t n_columns = n_terms * static_cast<std::int64_t>(n_groups);
  const std::int64_t nnz = n * n_terms;

  py::array_t<int> column_pointer(n_columns + 1);
  py::array_t<int> row_index(nnz);
  ColumnMajor values(nnz);
  py::array_t<int> cells_by_group(n);
  py::array_t<int> group_start(static_cast<std::int64_t>(n_groups) + 1);

  raise_if_failed(pace::random_design_block(
                      span_of(x_terms), pace::Span<const int>(cell_group.data(), n), n, n_terms,
                      n_groups, out_int_span(column_pointer), out_int_span(row_index),
                      out_span_of(values), out_int_span(cells_by_group),
                      out_int_span(group_start)),
                  "random design block");

  py::dict out;
  out["data"] = values;
  out["indices"] = row_index;
  out["indptr"] = column_pointer;
  out["shape"] = py::make_tuple(n, n_columns);
  out["cells_by_group"] = cells_by_group;
  out["group_start"] = group_start;
  return out;
}

// How the verbose log names the warm-up. R prints "Inf" for an unbounded one
// and the plain number otherwise, so a reported warm-up always matches the one
// in force.
std::string format_warmup(double value) {
  if (!std::isfinite(value)) return value > 0 ? "Inf" : "-Inf";
  if (value == std::floor(value)) return std::to_string(static_cast<long long>(value));
  return std::to_string(value);
}

// ---------------------------------------------------------------------------
// The IRLS loop: the fit itself.
// ---------------------------------------------------------------------------

// A compressed-column matrix handed over from Python. The arrays stay owned by
// numpy, so the holder keeps them alive for as long as the view is used.
struct CscHolder {
  py::array_t<int> indptr, indices;
  ColumnMajor data;
  pace::CscView view;

  CscHolder(const py::dict& d) {
    indptr = d["indptr"].cast<py::array_t<int>>();
    indices = d["indices"].cast<py::array_t<int>>();
    data = d["data"].cast<ColumnMajor>();
    const auto shape = d["shape"].cast<py::tuple>();
    view.column_pointer = pace::Span<const int>(indptr.data(), indptr.size());
    view.row_index = pace::Span<const int>(indices.data(), indices.size());
    view.values = pace::Span<const double>(data.data(), data.size());
    view.n_rows = shape[0].cast<std::int64_t>();
    view.n_cols = shape[1].cast<std::int64_t>();
  }
};

// One block of the random-effect design, as the loop reads it. Mirrors what
// random_design_block() returns, so the two compose directly.
struct BlockHolder {
  ColumnMajor terms;
  py::array_t<int> cell_group, group_start, cells_by_group;
  std::vector<int> group_size;
};

py::dict fit_irls(const ColumnMajor& x_fixed, std::int64_t p, const py::dict& z_design,
                  const py::list& blocks, const py::dict& counts, const py::dict& ambient,
                  const ColumnMajor& offset, const ColumnMajor& alpha_init,
                  const ColumnMajor& mask, const IntColumn& mask_index,
                  std::int64_t n_mask_rows, const ColumnMajor& data_informed_weights,
                  const ColumnMajor& sample_weight,
                  std::int64_t n_cells, std::int64_t q_total, std::int64_t n_genes,
                  const py::dict& options) {
  const std::int64_t n = n_cells, q = q_total, g = n_genes;

  CscHolder z(z_design), y(counts), e_tech(ambient);

  // Blocks: hold every array alive, then build the core's views over them.
  const std::size_t n_blocks = blocks.size();
  std::vector<BlockHolder> held(n_blocks);
  std::vector<pace::SolveBlock> solve_blocks(n_blocks);
  std::vector<pace::TauBlock> tau_blocks(n_blocks);
  for (std::size_t b = 0; b < n_blocks; ++b) {
    const py::dict block = blocks[b].cast<py::dict>();
    BlockHolder& h = held[b];
    h.terms = block["terms"].cast<ColumnMajor>();
    h.cell_group = block["cell_group"].cast<py::array_t<int>>();
    h.group_start = block["group_start"].cast<py::array_t<int>>();
    h.cells_by_group = block["cells_by_group"].cast<py::array_t<int>>();

    const int n_terms = block["n_terms"].cast<int>();
    const int n_groups = block["n_groups"].cast<int>();
    solve_blocks[b].col_offset = block["col_offset"].cast<int>();
    solve_blocks[b].n_terms = n_terms;
    solve_blocks[b].n_groups = n_groups;

    tau_blocks[b].col_offset = solve_blocks[b].col_offset;
    tau_blocks[b].n_terms = n_terms;
    tau_blocks[b].n_groups = n_groups;
    tau_blocks[b].n_cols = block["n_cols"].cast<int>();

    // group sizes, from the offsets the design already computed
    const int* start = h.group_start.data();
    h.group_size.resize(static_cast<std::size_t>(n_groups));
    for (int gi = 0; gi < n_groups; ++gi) h.group_size[static_cast<std::size_t>(gi)] =
        start[gi + 1] - start[gi];
  }
  // Spans only after every vector has stopped growing, so none can dangle.
  for (std::size_t b = 0; b < n_blocks; ++b) {
    BlockHolder& h = held[b];
    solve_blocks[b].terms = pace::Span<const double>(h.terms.data(), h.terms.size());
    solve_blocks[b].cell_group = pace::Span<const int>(h.cell_group.data(), h.cell_group.size());
    solve_blocks[b].group_offsets =
        pace::Span<const int>(h.group_start.data(), h.group_start.size());
    solve_blocks[b].group_cells =
        pace::Span<const int>(h.cells_by_group.data(), h.cells_by_group.size());
    tau_blocks[b].group_size =
        pace::Span<const int>(h.group_size.data(), static_cast<std::int64_t>(h.group_size.size()));
  }

  pace::IrlsLoopInputs inputs;
  // R derives both of these from the fixed design (engine-streaming.R): x1 is
  // its first column, and x1_is_unit is whether that column is all ones. Taking
  // x1_is_unit from the caller instead let a p == 1 fit with a non-constant
  // column build its linear predictor as if the column were ones while the
  // solver used the real values, which is wrong without being an error.
  const bool x_fixed_is_unit =
      p == 1 && x_fixed.size() > 0 &&
      std::all_of(x_fixed.data(), x_fixed.data() + x_fixed.size(),
                  [](double value) { return value == 1.0; });
  inputs.x1 = p == 1 && !x_fixed_is_unit ? span_of(x_fixed) : pace::Span<const double>(nullptr, 0);
  inputs.x1_is_unit = x_fixed_is_unit;
  inputs.x_fixed = p == 1 ? pace::Span<const double>(nullptr, 0) : span_of(x_fixed);
  inputs.solve_x_fixed = span_of(x_fixed);
  inputs.p = p;
  inputs.z = z.view;
  inputs.solve_blocks = &solve_blocks;
  inputs.tau_blocks = &tau_blocks;
  inputs.counts.matrix = y.view;
  inputs.counts.first_gene = 0;
  // Cached, or streamed. Cached takes the whole n x G ambient product; streamed
  // takes W instead and recomputes each chunk's columns from W and the counts,
  // which is what lets a full transcriptome panel fit in memory. The numbers
  // are the same either way: cache mode forms its product with the same
  // sparse_product_csc, so the two differ only in chunking, which is exact.
  // The scratch lives here, one set for the whole fit, because the core takes
  // the source by const reference and writes through these pointers.
  std::vector<int> ambient_scratch_p, ambient_scratch_i;
  std::vector<double> ambient_scratch_x;
  pace::CscView ambient_scratch_view;
  const bool ambient_streamed =
      options.contains("ambient_streamed") && options["ambient_streamed"].cast<bool>();
  if (ambient_streamed) {
    inputs.ambient.streamed = true;
    inputs.ambient.weights = e_tech.view;
    inputs.ambient.counts = y.view;
    inputs.ambient.scratch_column_pointer = &ambient_scratch_p;
    inputs.ambient.scratch_row_index = &ambient_scratch_i;
    inputs.ambient.scratch_values = &ambient_scratch_x;
    inputs.ambient.scratch_view = &ambient_scratch_view;
  } else {
    inputs.ambient.streamed = false;
    inputs.ambient.cached.matrix = e_tech.view;
    inputs.ambient.cached.first_gene = 0;
  }
  inputs.offset = span_of(offset);
  // The anchor mask and the data-informed weights are not optional in practice:
  // without the weights every random effect comes back unshrunk.
  inputs.sample_weight = span_of(sample_weight);
  inputs.mask = span_of(mask);
  inputs.mask_index = pace::Span<const int>(mask_index.data(), mask_index.size());
  inputs.n_mask_rows = n_mask_rows;
  inputs.data_informed_weights = span_of(data_informed_weights);
  inputs.alpha_init = span_of(alpha_init);
  inputs.n = n;
  inputs.q = q;
  inputs.n_genes = g;

  pace::IrlsLoopOptions opt;
  auto get_i = [&](const char* k, std::int64_t d) {
    return options.contains(k) ? options[k].cast<std::int64_t>() : d; };
  auto get_d = [&](const char* k, double d) {
    return options.contains(k) ? options[k].cast<double>() : d; };
  auto get_b = [&](const char* k, bool d) {
    return options.contains(k) ? options[k].cast<bool>() : d; };
  // Every default here mirrors fit_pace_mvpql_streaming()'s in R, so a fit run
  // through either language agrees without the caller restating them. They are
  // not free choices: alpha_warmup decides the iteration at which the
  // dispersion freezes, and moving it moves the fit.
  opt.n_iter = get_i("n_iter", 32);
  opt.min_iter = get_i("min_iter", 12);
  opt.early_stop_tol = get_d("early_stop_tol", 2e-2);
  opt.alpha_warmup = get_d("alpha_warmup", 6);
  opt.alpha_warmup_label = format_warmup(opt.alpha_warmup);
  opt.nb2 = get_b("nb2", false);
  opt.gaussian = get_b("gaussian", false);
  opt.zero_collapse = get_b("zero_collapse", true);
  opt.alpha_max_cells = get_d("alpha_max_cells", std::numeric_limits<double>::infinity());
  opt.fast_density = get_b("fast_density", true);
  opt.interior_precision = static_cast<int>(get_i("interior_precision", 1));
  opt.chunk_size = get_i("chunk_size", 128);
  opt.sub_genes = get_i("sub_genes", 16);
  opt.tau_max = get_d("tau_max", 100);
  opt.d0_min = get_d("d0_min", 1);
  opt.use_reml = get_b("use_reml", false);
  opt.rd_diag = get_b("rd_diag", false);
  // Explicit, because the struct's default order is deliberately not the
  // declaration order and R's own two entry points disagree: pace_fit_streaming
  // overrides fit_pace_mvpql_streaming's "hierarchical" with "adaptive", and
  // "adaptive" is what every PACE fit actually uses.
  opt.tau_shrinkage = static_cast<pace::TauShrinkage>(get_i("tau_shrinkage", 2));
  opt.verbose = get_b("verbose", false);
  opt.n_threads = static_cast<int>(get_i("n_threads", 1));

  ColumnMajor beta({p, g}), u({q, g}), se_beta({p, g}), se_u({q, g});
  ColumnMajor alpha(g), tau_g_array({q, g}), tau_flat(q), rho(n);
  ColumnMajor h_tau({q, opt.n_iter}), h_alpha({g, opt.n_iter}), h_rel(opt.n_iter);
  ColumnMajor h_rel_mean(opt.n_iter), h_nan(opt.n_iter), h_nonfinite(opt.n_iter);
  ColumnMajor h_tau_max(opt.n_iter);
  std::int64_t iterations_run = 0;
  bool converged = false;

  pace::IrlsLoopOutputs out;
  out.beta = out_span_of(beta);
  out.u = out_span_of(u);
  out.se_beta = out_span_of(se_beta);
  out.se_u = out_span_of(se_u);
  out.alpha = out_span_of(alpha);
  out.tau_g_array = out_span_of(tau_g_array);
  out.tau_flat = out_span_of(tau_flat);
  out.rho = out_span_of(rho);
  out.history_tau_flat = out_span_of(h_tau);
  out.history_alpha = out_span_of(h_alpha);
  out.history_rel_delta = out_span_of(h_rel);
  out.history_rel_delta_mean = out_span_of(h_rel_mean);
  out.history_n_nan_genes = out_span_of(h_nan);
  out.history_n_nonfinite = out_span_of(h_nonfinite);
  out.history_tau_max_seen = out_span_of(h_tau_max);
  out.iterations_run = &iterations_run;
  out.converged = &converged;

  // The standard errors are written on the final iteration only, so an
  // unfinished fit must report NaN rather than the zeros out_span_of() leaves.
  std::fill(se_beta.mutable_data(), se_beta.mutable_data() + se_beta.size(),
            std::numeric_limits<double>::quiet_NaN());
  std::fill(se_u.mutable_data(), se_u.mutable_data() + se_u.size(),
            std::numeric_limits<double>::quiet_NaN());

  pace::IrlsLoopReporter reporter;
  std::vector<std::string> warnings;
  // The loop below runs with the GIL released, and the reporter is called from
  // that thread, so touching the Python side needs the GIL back first. Without
  // this a verbose fit segfaults while a quiet one is fine, which is a very
  // misleading symptom.
  reporter.message = [&](const std::string& text) {
    if (!opt.verbose) return;
    py::gil_scoped_acquire hold;
    py::print(text, py::arg("end") = "", py::arg("flush") = true);
  };
  reporter.warn = [&warnings](const std::string& text) { warnings.push_back(text); };
  reporter.warn_nan_genes = [&warnings](std::int64_t iteration, pace::Span<const int>) {
    warnings.push_back("non-finite genes at iteration " + std::to_string(iteration));
  };

  // The density is never reached on the canonical path: NB1 with fast_density
  // uses the core's own likelihood. Anything else needs a real one.
  pace::LogDensity density = [](double, double, double) {
    return std::numeric_limits<double>::quiet_NaN();
  };

  {
    py::gil_scoped_release release;   // the loop is C++ and may use threads
    raise_if_failed(pace::run_irls_loop(inputs, opt, density, pace::trigamma, reporter, out,
                                        pace::InterruptCheck()),
                    "irls loop");
  }

  py::dict result;
  result["B"] = beta;
  result["U"] = u;
  result["se_B"] = se_beta;
  result["se_U"] = se_u;
  result["alpha"] = alpha;
  result["tau_g_array"] = tau_g_array;
  result["tau_flat"] = tau_flat;
  result["rho"] = rho;
  result["n_iter"] = iterations_run;
  result["converged"] = converged;
  result["warnings"] = warnings;
  return result;
}

// ---------------------------------------------------------------------------
// The preprocessing stages between the kernels and the fit.
// ---------------------------------------------------------------------------

// Zero any (focal, neighbour) pair whose effective sample size falls below
// `min_effective`: a slope estimated from almost no contrast is noise. Modifies
// a copy, so the caller's kernel is untouched.
py::dict drop_sparse_kernel(const ColumnMajor& kernel, const IntColumn& cell_type,
                            double min_effective) {
  if (kernel.ndim() != 2) throw std::runtime_error("kernel must be n_cells x n_types");
  const std::int64_t n = kernel.shape(0), n_types = kernel.shape(1);
  if (cell_type.size() != n) throw std::runtime_error("cell_type must have one entry per cell");

  ColumnMajor out({n, n_types});
  std::copy(kernel.data(), kernel.data() + kernel.size(), out.mutable_data());
  std::int64_t n_dropped = 0;
  raise_if_failed(pace::drop_sparse_kernel(
                      pace::Span<double>(out.mutable_data(), out.size()), n, n_types,
                      pace::Span<const int>(cell_type.data(), n), min_effective, &n_dropped),
                  "drop sparse kernel");
  py::dict result;
  result["kernel"] = out;
  result["n_dropped"] = n_dropped;
  return result;
}

// Subtract each (image, cell type) group's mean from every kernel column, so a
// slope is estimated on within-group deviation rather than on composition.
ColumnMajor centre_within_groups(const ColumnMajor& kernel, const IntColumn& image,
                                 const IntColumn& cell_type, int n_images, int n_types) {
  if (kernel.ndim() != 2) throw std::runtime_error("kernel must be n_cells x n_types");
  const std::int64_t n = kernel.shape(0), k = kernel.shape(1);
  if (image.size() != n || cell_type.size() != n) {
    throw std::runtime_error("image and cell_type must have one entry per cell");
  }
  ColumnMajor out({n, k});
  std::copy(kernel.data(), kernel.data() + kernel.size(), out.mutable_data());
  raise_if_failed(pace::centre_within_groups(
                      pace::Span<double>(out.mutable_data(), out.size()), n, k,
                      pace::Span<const int>(image.data(), n),
                      pace::Span<const int>(cell_type.data(), n), n_images, n_types),
                  "centre within groups");
  return out;
}

// Which genes anchor which cell type: a gene owned by some other type, present
// enough there, and near-absent here. These are what pin the contamination
// field, so the fit can tell borrowed signal from real.
py::dict anchor_mask(const ColumnMajor& core_means, double owner_threshold,
                     double core_threshold) {
  if (core_means.ndim() != 2) throw std::runtime_error("core_means must be n_types x n_genes");
  const std::int64_t n_types = core_means.shape(0), n_genes = core_means.shape(1);
  ColumnMajor mask({n_types, n_genes});
  py::array_t<int> n_anchor(n_types);
  raise_if_failed(pace::anchor_mask(span_of(core_means), n_types, n_genes, owner_threshold,
                                    core_threshold, out_span_of(mask), out_int_span(n_anchor)),
                  "anchor mask");
  py::dict out;
  out["mask"] = mask;
  out["n_anchor"] = n_anchor;
  return out;
}

// The per-(random effect, gene) prior weights, from how well each gene is
// detected in the focal type that row belongs to.
ColumnMajor data_informed_weights(const ColumnMajor& detection_rate,
                                  const IntColumn& focal_of_row,
                                  const ColumnMajor& scale_of_row) {
  if (detection_rate.ndim() != 2) {
    throw std::runtime_error("detection_rate must be n_focals x n_genes");
  }
  const std::int64_t n_focals = detection_rate.shape(0), n_genes = detection_rate.shape(1);
  const std::int64_t q = focal_of_row.size();
  if (scale_of_row.size() != q) {
    throw std::runtime_error("scale_of_row must have one entry per random-effect row");
  }
  ColumnMajor weights({q, n_genes});
  raise_if_failed(pace::data_informed_weights(span_of(detection_rate), n_focals, n_genes,
                                              pace::Span<const int>(focal_of_row.data(), q),
                                              span_of(scale_of_row), q, out_span_of(weights)),
                  "data informed weights");
  return weights;
}

// R's scale() of one column, with the divisor it used, so a caller can tell a
// column that carries no slope from one that does.
py::dict standardise(const ColumnMajor& values) {
  const std::int64_t n = values.size();
  ColumnMajor out(n);
  std::copy(values.data(), values.data() + n, out.mutable_data());
  double sd = 0.0;
  raise_if_failed(pace::standardise_column(pace::Span<double>(out.mutable_data(), n), n, &sd),
                  "standardise");
  py::dict result;
  result["values"] = out;
  result["sd"] = sd;
  return result;
}

// Mean of every gene within each cell group, or the detection rate when
// `detection`. Cells coded -1 belong to no group; an empty group gives NaN,
// as colMeans() of no rows does.
ColumnMajor group_column_means(const py::dict& counts, const IntColumn& group, int n_groups,
                               bool detection, int n_threads) {
  CscHolder y(counts);
  const std::int64_t n = y.view.n_rows, n_genes = y.view.n_cols;
  if (group.size() != n) {
    throw std::runtime_error("group must have one entry per cell; note the counts are "
                             "cells by genes here, not genes by cells");
  }
  ColumnMajor means({static_cast<std::int64_t>(n_groups), n_genes});
  raise_if_failed(pace::group_column_means(y.view, pace::Span<const int>(group.data(), n),
                                           n_groups, detection, out_span_of(means), n_threads,
                                           pace::InterruptCheck()),
                  "group column means");
  return means;
}

// Fraction of each cell's neighbours within `radius` that share its type. Used
// to find homotypic-core cells, whose profiles are clean enough to anchor the
// contamination field.
ColumnMajor same_type_fraction(const ColumnMajor& x, const ColumnMajor& y,
                               const IntColumn& cell_type, const IntColumn& image,
                               double radius, int min_image_cells, int n_threads) {
  const std::int64_t n = x.size();
  if (y.size() != n || cell_type.size() != n || image.size() != n) {
    throw std::runtime_error("every input must have one entry per cell");
  }
  ColumnMajor fraction(n);
  raise_if_failed(pace::same_type_fraction(span_of(x), span_of(y),
                                           pace::Span<const int>(cell_type.data(), n),
                                           pace::Span<const int>(image.data(), n), radius,
                                           min_image_cells, out_span_of(fraction), n_threads,
                                           pace::InterruptCheck()),
                  "same type fraction");
  return fraction;
}

// Covariance of the kernel columns within each cell group, as n_groups slices
// of p x p. Its diagonal is the per-focal variance of each neighbour column,
// which is what scales the prior weights.
ColumnMajor group_covariances(const ColumnMajor& values, const IntColumn& group, int n_groups,
                              int n_threads) {
  if (values.ndim() != 2) throw std::runtime_error("values must be n_cells x n_columns");
  const std::int64_t n = values.shape(0), p = values.shape(1);
  if (group.size() != n) throw std::runtime_error("group must have one entry per cell");
  ColumnMajor covariance(static_cast<std::int64_t>(n_groups) * p * p);
  raise_if_failed(pace::group_covariances(span_of(values), n, p,
                                          pace::Span<const int>(group.data(), n), n_groups,
                                          out_span_of(covariance), n_threads,
                                          pace::InterruptCheck()),
                  "group covariances");
  return covariance;
}

// Scale every row to its own maximum and clamp at one, leaving a row with no
// usable maximum at zero. The prior weights are scale-free across focals
// because of this.
ColumnMajor normalise_rows_to_max(const ColumnMajor& values) {
  if (values.ndim() != 2) throw std::runtime_error("values must be 2-D");
  const std::int64_t rows = values.shape(0), cols = values.shape(1);
  ColumnMajor out({rows, cols});
  std::copy(values.data(), values.data() + values.size(), out.mutable_data());
  raise_if_failed(pace::normalise_rows_to_max(
                      pace::Span<double>(out.mutable_data(), out.size()), rows, cols),
                  "normalise rows");
  return out;
}

// Sparse product, in compressed-column form. The ambient field the solver reads
// is W @ Y, and it has to be formed with this rather than with scipy: the
// summation order is part of the answer, and the R package routes the same
// product through here for exactly that reason.
py::dict sparse_product(const py::dict& left, const py::dict& right, int n_threads) {
  CscHolder a(left), b(right);
  if (a.view.n_cols != b.view.n_rows) {
    throw std::runtime_error("left's columns must match right's rows");
  }
  std::vector<int> column_pointer, row_index;
  std::vector<double> values;
  raise_if_failed(pace::sparse_product_csc(a.view, b.view, 0, b.view.n_cols, column_pointer,
                                           row_index, values, n_threads,
                                           pace::InterruptCheck()),
                  "sparse product");

  const std::int64_t nnz = static_cast<std::int64_t>(values.size());
  py::array_t<int> indptr(static_cast<std::int64_t>(column_pointer.size()));
  py::array_t<int> indices(nnz);
  ColumnMajor data(nnz);
  std::copy(column_pointer.begin(), column_pointer.end(), indptr.mutable_data());
  std::copy(row_index.begin(), row_index.end(), indices.mutable_data());
  std::copy(values.begin(), values.end(), data.mutable_data());

  py::dict out;
  out["data"] = data;
  out["indices"] = indices;
  out["indptr"] = indptr;
  out["shape"] = py::make_tuple(a.view.n_rows, b.view.n_cols);
  out["nnz"] = nnz;
  return out;
}

// The per-(focal, gene) variance decomposition. The arguments are in the same
// order as pace_variance_decomposition_cpp()'s in R, so the two assemblies can
// be read side by side. Index vectors are 1-based, as R indexes them, and a
// non-positive entry means the row is absent.
py::dict variance_decomposition(
    const ColumnMajor& ct_means, const IntColumn& group_size, const IntColumn& focal_group,
    const IntColumn& n_focal, const ColumnMajor& u, const ColumnMajor& se_u,
    const IntMatrix& slope_rows, std::int64_t n_kernel, const ColumnMajor& kernel_cov,
    const IntMatrix& responder_rows, std::int64_t n_responder, const IntColumn& responder_keep,
    const ColumnMajor& responder_cov, const IntColumn& intercept_rows,
    const ColumnMajor& toff_var, const ColumnMajor& spill_cov, const ColumnMajor& beta_spill,
    std::int64_t n_spill, const ColumnMajor& mu_mean, const ColumnMajor& alpha, bool nb1,
    int n_threads) {
  pace::DecompositionInput input;
  input.n_genes = ct_means.ndim() == 2 ? ct_means.shape(1) : 0;
  input.n_groups = static_cast<int>(group_size.size());
  input.n_focals = static_cast<int>(focal_group.size());
  input.ct_means = span_of(ct_means);
  input.group_size = pace::Span<const int>(group_size.data(), group_size.size());
  input.focal_group = pace::Span<const int>(focal_group.data(), focal_group.size());
  input.n_focal = pace::Span<const int>(n_focal.data(), n_focal.size());
  input.u = span_of(u);
  input.se_u = span_of(se_u);
  input.q = u.ndim() == 2 ? u.shape(0) : 0;
  input.n_kernel = n_kernel;
  input.slope_rows = pace::Span<const int>(slope_rows.data(), slope_rows.size());
  input.kernel_cov = span_of(kernel_cov);
  input.n_responder = n_responder;
  input.responder_rows = pace::Span<const int>(responder_rows.data(), responder_rows.size());
  input.responder_keep = pace::Span<const int>(responder_keep.data(), responder_keep.size());
  input.responder_cov = span_of(responder_cov);
  input.intercept_rows = pace::Span<const int>(intercept_rows.data(), intercept_rows.size());
  input.toff_var = span_of(toff_var);
  input.n_spill = n_spill;
  input.spill_cov = span_of(spill_cov);
  input.beta_spill = span_of(beta_spill);
  input.mu_mean = span_of(mu_mean);
  input.alpha = span_of(alpha);
  input.nb1 = nb1;

  const std::int64_t rows = static_cast<std::int64_t>(input.n_focals) * input.n_genes;
  ColumnMajor v_state_baseline(rows), v_state_responder(rows), v_spill(rows), v_disp(rows);
  ColumnMajor celltype_offset_sq(rows), focal_mean(rows), max_other_mean(rows), spec(rows);
  ColumnMajor focal_other_ratio(rows), pct_celltype(rows), pct_state(rows);
  ColumnMajor pct_responder(rows), pct_spill(rows), pct_residual(rows);
  py::array_t<int> is_contaminated(rows), keep(rows);
  std::fill(is_contaminated.mutable_data(), is_contaminated.mutable_data() + rows, 0);
  std::fill(keep.mutable_data(), keep.mutable_data() + rows, 0);

  pace::DecompositionOutput output;
  output.v_state_baseline = out_span_of(v_state_baseline);
  output.v_state_responder = out_span_of(v_state_responder);
  output.v_spill = out_span_of(v_spill);
  output.v_disp = out_span_of(v_disp);
  output.celltype_offset_sq = out_span_of(celltype_offset_sq);
  output.focal_mean = out_span_of(focal_mean);
  output.max_other_mean = out_span_of(max_other_mean);
  output.spec = out_span_of(spec);
  output.focal_other_ratio = out_span_of(focal_other_ratio);
  output.is_contaminated = pace::Span<int>(is_contaminated.mutable_data(), rows);
  output.pct_celltype = out_span_of(pct_celltype);
  output.pct_state = out_span_of(pct_state);
  output.pct_responder = out_span_of(pct_responder);
  output.pct_spill = out_span_of(pct_spill);
  output.pct_residual = out_span_of(pct_residual);
  output.keep = pace::Span<int>(keep.mutable_data(), rows);

  raise_if_failed(pace::variance_decomposition(input, output, n_threads,
                                               pace::InterruptCheck()),
                  "variance decomposition");

  py::dict result;
  result["v_state_baseline"] = v_state_baseline;
  result["v_state_responder"] = v_state_responder;
  result["v_spill"] = v_spill;
  result["v_disp"] = v_disp;
  result["celltype_offset_sq"] = celltype_offset_sq;
  result["focal_mean"] = focal_mean;
  result["max_other_mean"] = max_other_mean;
  result["spec"] = spec;
  result["focal_other_ratio"] = focal_other_ratio;
  result["is_contaminated"] = is_contaminated;
  result["pct_celltype"] = pct_celltype;
  result["pct_state"] = pct_state;
  result["pct_responder"] = pct_responder;
  result["pct_spill"] = pct_spill;
  result["pct_residual"] = pct_residual;
  result["keep"] = keep;
  return result;
}

// The aggregate tables over the rows of one focal: the plain mean of each
// block's percentage, its specificity-weighted mean, and the pooled share.
py::dict decomposition_aggregates(const IntColumn& focal_code, int n_focals,
                                  const ColumnMajor& pct, const ColumnMajor& components,
                                  const ColumnMajor& spec, const IntColumn& is_contaminated) {
  const std::int64_t n_rows = focal_code.size();
  ColumnMajor mean_table({static_cast<std::int64_t>(n_focals), std::int64_t{5}});
  ColumnMajor specw_table({static_cast<std::int64_t>(n_focals), std::int64_t{5}});
  ColumnMajor pooled_table({static_cast<std::int64_t>(n_focals), std::int64_t{5}});
  ColumnMajor total(n_rows), total_ss(static_cast<std::int64_t>(n_focals));
  py::array_t<int> n_rows_out(n_focals), n_specific(n_focals);
  std::fill(n_rows_out.mutable_data(), n_rows_out.mutable_data() + n_focals, 0);
  std::fill(n_specific.mutable_data(), n_specific.mutable_data() + n_focals, 0);

  raise_if_failed(pace::decomposition_aggregates(
                      pace::Span<const int>(focal_code.data(), n_rows), n_focals, n_rows,
                      span_of(pct), span_of(components), span_of(spec),
                      pace::Span<const int>(is_contaminated.data(), is_contaminated.size()),
                      out_span_of(mean_table), out_span_of(specw_table),
                      out_span_of(pooled_table), out_span_of(total), out_span_of(total_ss),
                      pace::Span<int>(n_rows_out.mutable_data(), n_focals),
                      pace::Span<int>(n_specific.mutable_data(), n_focals)),
                  "decomposition aggregates");

  py::dict result;
  result["mean"] = mean_table;
  result["specw"] = specw_table;
  result["pooled"] = pooled_table;
  result["total"] = total;
  result["total_ss"] = total_ss;
  result["n_rows"] = n_rows_out;
  result["n_specific"] = n_specific;
  return result;
}

// The observed single-frame statistics on y = log1p(CP10k): each type's mean,
// its within-type sum of squares, and the global mean per gene.
py::dict single_frame_statistics(const py::dict& counts, const ColumnMajor& n_count,
                                 const IntColumn& group, int n_groups, int n_threads) {
  CscHolder y(counts);
  const std::int64_t n_genes = y.view.n_cols;
  ColumnMajor focal_mean({static_cast<std::int64_t>(n_groups), n_genes});
  ColumnMajor within_ss({static_cast<std::int64_t>(n_groups), n_genes});
  ColumnMajor global_mean(n_genes);
  raise_if_failed(pace::single_frame_statistics(
                      y.view, span_of(n_count),
                      pace::Span<const int>(group.data(), group.size()), n_groups,
                      out_span_of(focal_mean), out_span_of(within_ss),
                      out_span_of(global_mean), n_threads, pace::InterruptCheck()),
                  "single frame statistics");
  py::dict out;
  out["focal_mean"] = focal_mean;
  out["within_ss"] = within_ss;
  out["global_mean"] = global_mean;
  return out;
}

// The observed single-frame shares. Mind the two layouts: the inputs are
// n_groups x n_genes column-major, so (type, gene) sits at type + gene *
// n_groups, while the OUTPUTS are grouped by type, so (type, gene) sits at
// type * n_genes + gene. R relies on exactly that, and they are not the same
// order.
py::dict single_frame_shares(const ColumnMajor& focal_mean, const ColumnMajor& within_ss,
                             const ColumnMajor& global_mean, const IntColumn& group_size,
                             const ColumnMajor& v_state, const ColumnMajor& v_responder,
                             const ColumnMajor& v_spill, const ColumnMajor& v_disp,
                             int n_threads) {
  const int n_groups = static_cast<int>(group_size.size());
  const std::int64_t n_genes = focal_mean.ndim() == 2 ? focal_mean.shape(1) : 0;
  const std::int64_t rows = static_cast<std::int64_t>(n_groups) * n_genes;
  const bool has_responder = v_responder.size() != 0;
  ColumnMajor pct_celltype(rows), pct_spatial(rows), pct_spill(rows), pct_residual(rows);
  ColumnMajor pct_responder(has_responder ? rows : 0);
  ColumnMajor ss_lineage(rows), ss_within(rows), denom(rows);
  py::array_t<int> keep(n_groups);
  std::fill(keep.mutable_data(), keep.mutable_data() + n_groups, 0);

  raise_if_failed(pace::single_frame_shares(
                      span_of(focal_mean), span_of(within_ss), span_of(global_mean),
                      pace::Span<const int>(group_size.data(), n_groups), n_groups, n_genes,
                      span_of(v_state), span_of(v_responder), span_of(v_spill),
                      span_of(v_disp), out_span_of(pct_celltype), out_span_of(pct_spatial),
                      out_span_of(pct_responder), out_span_of(pct_spill),
                      out_span_of(pct_residual), out_span_of(ss_lineage),
                      out_span_of(ss_within), out_span_of(denom),
                      pace::Span<int>(keep.mutable_data(), n_groups), n_threads,
                      pace::InterruptCheck()),
                  "single frame shares");
  py::dict out;
  out["pct_celltype"] = pct_celltype;
  out["pct_spatial"] = pct_spatial;
  out["pct_responder"] = pct_responder;
  out["pct_spill"] = pct_spill;
  out["pct_residual"] = pct_residual;
  out["ss_lineage"] = ss_lineage;
  out["ss_within"] = ss_within;
  out["denom"] = denom;
  out["keep"] = keep;
  return out;
}

// Sum-of-squares weighted per-focal blocks of the observed single frame.
ColumnMajor single_frame_focal_blocks(const IntColumn& focal_code, int n_focals,
                                      const ColumnMajor& pct, int n_blocks,
                                      const ColumnMajor& denom) {
  const std::int64_t n_rows = focal_code.size();
  ColumnMajor table({static_cast<std::int64_t>(n_focals), static_cast<std::int64_t>(n_blocks)});
  raise_if_failed(pace::single_frame_focal_blocks(
                      pace::Span<const int>(focal_code.data(), n_rows), n_focals, n_rows,
                      span_of(pct), n_blocks, span_of(denom), out_span_of(table)),
                  "single frame focal blocks");
  return table;
}

}  // namespace



PYBIND11_MODULE(_core, m) {
  m.doc() = "PACE compiled core";
  m.def("mash", &mash, py::arg("bhat"), py::arg("shat"),
        py::arg("v") = ColumnMajor({0, 0}), py::arg("null_weight") = 10.0,
        py::arg("grid_mult") = 1.4142135623730951, py::arg("max_iterations") = 200000,
        py::arg("tolerance") = 1e-13,
        py::arg("condition_names") = std::vector<std::string>(),
        "Canonical multivariate adaptive shrinkage for one slice.");
  m.def("driver_scores", &driver_scores, py::arg("estimate_shrunk"), py::arg("spec"),
        py::arg("focal_mean"), py::arg("mu_bar"), py::arg("alpha"),
        py::arg("u_raw") = ColumnMajor(0), py::arg("var_n") = 0.0, py::arg("var_rn") = 0.0,
        py::arg("has_responder") = false,
        "Driver scores for one focal-neighbour pair.");
  m.def("pair_variance_pratt", &pair_variance_pratt, py::arg("sigma"), py::arg("u"),
        "Pratt attribution of spatial variance across a focal type's neighbours.");
  m.def("neighbour_kernels", &neighbour_kernels, py::arg("x"), py::arg("y"),
        py::arg("cell_type"), py::arg("n_types"), py::arg("group"),
        py::arg("per_group") = false, py::arg("h_bio") = 30.0, py::arg("h_tech") = 5.0,
        py::arg("eps") = 0.0, py::arg("n_threads") = 1,
        "Gaussian biological and exponential technical neighbour kernels.");
  m.def("ambient_field", &ambient_field, py::arg("x"), py::arg("y"), py::arg("cell_type"),
        py::arg("image"), py::arg("n_images"), py::arg("h_tech") = 5.0,
        py::arg("edge_correct") = true, py::arg("n_angles") = 1000, py::arg("n_threads") = 1,
        "Sparse cross-type ambient weight matrix, in compressed-column form.");
  m.def("random_design_block", &random_design_block, py::arg("x_terms"),
        py::arg("cell_group"), py::arg("n_groups"),
        "One block of the random-effect design, in compressed-column form.");
  m.def("fit_irls", &fit_irls, py::arg("x_fixed"), py::arg("p"), py::arg("z_design"),
        py::arg("blocks"), py::arg("counts"), py::arg("ambient"), py::arg("offset"),
        py::arg("alpha_init"), py::arg("mask") = ColumnMajor(0),
        py::arg("mask_index") = IntColumn(0), py::arg("n_mask_rows") = 0,
        py::arg("data_informed_weights") = ColumnMajor(0),
        py::arg("sample_weight") = ColumnMajor(0), py::arg("n_cells"), py::arg("q_total"), py::arg("n_genes"),
        py::arg("options") = py::dict(),
        "The penalised quasi-likelihood IRLS loop: the fit itself.");
  m.def("drop_sparse_kernel", &drop_sparse_kernel, py::arg("kernel"), py::arg("cell_type"),
        py::arg("min_effective"), "Zero pairs with too little effective sample size.");
  m.def("centre_within_groups", &centre_within_groups, py::arg("kernel"), py::arg("image"),
        py::arg("cell_type"), py::arg("n_images"), py::arg("n_types"),
        "Centre each kernel column within every (image, cell type) group.");
  m.def("anchor_mask", &anchor_mask, py::arg("core_means"), py::arg("owner_threshold"),
        py::arg("core_threshold"), "Which genes anchor which cell type.");
  m.def("data_informed_weights", &data_informed_weights, py::arg("detection_rate"),
        py::arg("focal_of_row"), py::arg("scale_of_row"),
        "Per-(random effect, gene) prior weights from detection.");
  m.def("standardise", &standardise, py::arg("values"),
        "R's scale() of one column, with the divisor it used.");
  m.def("group_column_means", &group_column_means, py::arg("counts"), py::arg("group"),
        py::arg("n_groups"), py::arg("detection") = false, py::arg("n_threads") = 1,
        "Per-group gene means, or detection rates.");
  m.def("same_type_fraction", &same_type_fraction, py::arg("x"), py::arg("y"),
        py::arg("cell_type"), py::arg("image"), py::arg("radius") = 30.0,
        py::arg("min_image_cells") = 50, py::arg("n_threads") = 1,
        "Fraction of each cell's neighbours sharing its type.");
  m.def("group_covariances", &group_covariances, py::arg("values"), py::arg("group"),
        py::arg("n_groups"), py::arg("n_threads") = 1,
        "Covariance of the columns within each group, n_groups slices of p x p.");
  m.def("normalise_rows_to_max", &normalise_rows_to_max, py::arg("values"),
        "Scale each row to its own maximum, clamped at one.");
  m.def("sparse_product", &sparse_product, py::arg("left"), py::arg("right"),
        py::arg("n_threads") = 1,
        "Sparse matrix product in compressed-column form.");
  m.def("variance_decomposition", &variance_decomposition,
        py::arg("ct_means"), py::arg("group_size"), py::arg("focal_group"), py::arg("n_focal"),
        py::arg("u"), py::arg("se_u"), py::arg("slope_rows"), py::arg("n_kernel"),
        py::arg("kernel_cov"), py::arg("responder_rows"), py::arg("n_responder"),
        py::arg("responder_keep"), py::arg("responder_cov"), py::arg("intercept_rows"),
        py::arg("toff_var"), py::arg("spill_cov"), py::arg("beta_spill"), py::arg("n_spill"),
        py::arg("mu_mean"), py::arg("alpha"), py::arg("nb1") = true, py::arg("n_threads") = 1,
        "The per-(focal, gene) variance decomposition.");
  m.def("decomposition_aggregates", &decomposition_aggregates,
        py::arg("focal_code"), py::arg("n_focals"), py::arg("pct"), py::arg("components"),
        py::arg("spec"), py::arg("is_contaminated"),
        "Per-focal aggregate tables of the decomposition.");
  m.def("single_frame_statistics", &single_frame_statistics, py::arg("counts"),
        py::arg("n_count"), py::arg("group"), py::arg("n_groups"), py::arg("n_threads") = 1,
        "Observed per-type means and within-type sums of squares on log1p(CP10k).");
  m.def("single_frame_shares", &single_frame_shares, py::arg("focal_mean"),
        py::arg("within_ss"), py::arg("global_mean"), py::arg("group_size"),
        py::arg("v_state"), py::arg("v_responder"), py::arg("v_spill"), py::arg("v_disp"),
        py::arg("n_threads") = 1,
        "The observed single-frame decomposition shares.");
  m.def("single_frame_focal_blocks", &single_frame_focal_blocks, py::arg("focal_code"),
        py::arg("n_focals"), py::arg("pct"), py::arg("n_blocks"), py::arg("denom"),
        "Sum-of-squares weighted per-focal blocks of the observed single frame.");
}
