#! /usr/bin/env python

import sys
import os
import re
import math
import glob
import logging

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
_LOG = logging.getLogger(os.path.basename(__file__))

import pycoevolity
import project_util

import matplotlib as mpl

# Use TrueType (42) fonts rather than Type 3 fonts
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
tex_font_settings = {
        "text.usetex": True,
        "font.family": "sans-serif",
        "text.latex.preamble" : [
                "\\usepackage[T1]{fontenc}",
                "\\usepackage[cm]{sfmath}",
                ]
}

mpl.rcParams.update(tex_font_settings)

import matplotlib.pyplot as plt
from matplotlib import gridspec
import scipy.stats
import numpy



sim_dir_to_priors = {
        "06pops-dpp-root-0005-004-t0002-500k":
            (
                (5.0, 0.04, 0.05),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-dpp-root-0005-004-t0002-500k":
            (
                (5.0, 0.04, 0.05),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-dpp-root-0005-004-t0002-500k-40genomes":
            (
                (5.0, 0.04, 0.05),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-dpp-root-0005-004-t0002-500k-diffuseprior":
            (
                (5.0, 0.04, 0.05),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-dpp-root-0005-004-3_8-t0002-500k-diffuseprior":
            (
                (5.0, 0.04, 3.8),
                (4.0, 0.000475, 0.0001)
            ),
        "03pairs-dpp-root-0005-004-t0002-500k-diffuseprior":
            (
                (5.0, 0.04, 0.05),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-dpp-root-0005-001-t0002-500k-diffuseprior":
            (
                (5.0, 0.01, 0.05),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-dpp-root-0005-001-9_95-t0002-500k-diffuseprior":
            (
                (5.0, 0.01, 9.95),
                (4.0, 0.000475, 0.0001)
            ),
        "03pairs-dpp-root-0005-004-t0002-500k":
            (
                (5.0, 0.04, 0.05),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-dpp-root-0005-004-t0002-500k-0100l":
            (
                (5.0, 0.04, 0.05),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-dpp-root-0005-004-t0002-500k-0100ul":
            (
                (5.0, 0.04, 0.05),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-dpp-root-0005-004-t0002-500k-040s":
            (
                (5.0, 0.04, 0.05),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-dpp-root-0005-004-t0002-500k-060s":
            (
                (5.0, 0.04, 0.05),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-dpp-root-0005-009-t0002-500k":
            (
                (5.0, 0.09, 0.05),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-03pairs-dpp-root-0005-009-t0002-500k":
            (
                (5.0, 0.09, 0.05),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-dpp-root-0005-019-t0002-500k":
            (
                (5.0, 0.19, 0.05),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-dpp-root-0005-079-t0002-500k":
            (
                (5.0, 0.79, 0.05),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-dpp-root-0010-0025-500k":
            (
                (10.0, 0.025, 0.0),
                (1.0, 0.01, 0.0)
            ),
        "03pops-dpp-root-0010-0025-t0001-500k":
            (
                (10.0, 0.025, 0.0),
                (1.0, 0.001, 0.0)
            ),
        "03pops-dpp-root-0010-005-500k":
            (
                (10.0, 0.05, 0.0),
                (1.0, 0.01, 0.0)
            ),
        "03pops-dpp-root-0010-005-t0001-500k":
            (
                (10.0, 0.05, 0.0),
                (1.0, 0.001, 0.0)
            ),
        "03pops-dpp-root-0010-010-500k":
            (
                (10.0, 0.1, 0.0),
                (1.0, 0.01, 0.0)
            ),
        "03pops-dpp-root-0010-020-500k":
            (
                (10.0, 0.2, 0.0),
                (1.0, 0.01, 0.0)
            ),
        "03pops-dpp-root-0050-002-t0002-500k":
            (
                (50.0, 0.02, 0.0),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-dpp-root-0100-001-500k":
            (
                (100.0, 0.01, 0.0),
                (1.0, 0.01, 0.0)
            ),
        "03pops-dpp-root-0500-0002-t0002-500k":
            (
                (500.0, 0.002, 0.0),
                (4.0, 0.000475, 0.0001)
            ),
        "03pops-dpp-root-1000-0001-500k":
            (
                (1000.0, 0.001, 0.0),
                (1.0, 0.01, 0.0)
            ),
}

unique_time_priors = set(p[1] for p in sim_dir_to_priors.values())


def get_prefix_from_sim_dir_name(sim_dir_name):
    conc_prior, time_prior = sim_dir_to_priors[sim_dir_name]
    s = "a-{a_shape}-{a_scale}-{a_offset}-t-{t_shape}-{t_scale}-{t_offset}".format(
            a_shape = str(conc_prior[0]).replace(".", "_"),
            a_scale = str(conc_prior[1]).replace(".", "_"),
            a_offset = str(conc_prior[2]).replace(".", "_"),
            t_shape = str(time_prior[0]).replace(".", "_"),
            t_scale = str(time_prior[1]).replace(".", "_"),
            t_offset = str(time_prior[2]).replace(".", "_"),
            )
    if sim_dir_name.endswith("-0100l"):
        s += "-loci-100-allsites"
    elif sim_dir_name.endswith("-0100ul"):
        s += "-loci-100-unlinkedsnps"
    elif sim_dir_name.endswith("-040s"):
        s += "-singletonprob-0_4"
    elif sim_dir_name.endswith("-060s"):
        s += "-singletonprob-0_6"
    elif sim_dir_name == "03pops-dpp-root-0005-001-t0002-500k-diffuseprior":
        s += "-diffuseprior-10increase"
    elif sim_dir_name == "03pops-dpp-root-0005-001-9_95-t0002-500k-diffuseprior":
        s += "-diffuseprior-10decrease"
    elif sim_dir_name == "03pops-dpp-root-0005-004-t0002-500k-diffuseprior":
        s += "-diffuseprior-4increase"
    elif sim_dir_name == "03pairs-dpp-root-0005-004-t0002-500k-diffuseprior":
        s += "-diffuseprior-pairs-4increase"
    elif sim_dir_name == "03pairs-dpp-root-0005-004-t0002-500k":
        s += "-pairs"
    elif sim_dir_name == "03pops-dpp-root-0005-004-3_8-t0002-500k-diffuseprior":
        s += "-diffuseprior-4decrease"
    elif sim_dir_name == "03pops-dpp-root-0005-004-t0002-500k-40genomes":
        s += "-40genomes"
    elif sim_dir_name == "06pops-dpp-root-0005-004-t0002-500k":
        s += "-06pops"
    if sim_dir_name.startswith("03pops-03pairs"):
        s += "-mixed-comps"
    return s


class BoxData(object):
    def __init__(self, values = [], labels = []):
        assert len(values) == len(labels)
        self._values = values
        self._labels = labels

    @classmethod
    def init(cls, results, parameters, labels = None):
        if labels:
            assert len(parameters) == len(labels)
        bd = cls()
        bd._values = [[] for i in range(len(parameters))]
        if labels:
            bd._labels = list(labels)
        else:
            bd._labels = list(parameters)
        for i, parameter_str in enumerate(parameters):
            x = [float(e) for e in results["{0}".format(parameter_str)]]
            bd._d[i].append(x)
        return bd

    def _get_values(self):
        return self._values

    values = property(_get_values)

    def _get_labels(self):
        return self._labels

    labels = property(_get_labels)

    def _get_number_of_categories(self):
        return len(self._values)

    number_of_categories = property(_get_number_of_categories)


    @classmethod
    def init_time_v_sharing(cls, results,
            estimator_prefix = "mean"):

        bd_abs_error = cls()
        bd_ci_width = cls()

        nreplicates = len(results["true_model"])
        ncomparisons = len(results["true_model"][0])

        labels = list(range(ncomparisons))
        err_vals = [[] for l in labels]
        ci_vals = [[] for l in labels]

        for sim_index in range(nreplicates):
            true_model = [int(x) for x in list(results["true_model"][sim_index])]
            assert len(true_model) == ncomparisons
            # Only use first comparison so as not to multply count the same
            # parameter estimates
            comparison_index = 0
            # for comparison_index in range(ncomparisons):
            num_shared = true_model.count(true_model[comparison_index]) - 1
            height_key_suffix = "root_height_c{0}sp1".format(comparison_index + 1)
            true_height = float(results["true_{0}".format(height_key_suffix)][sim_index])
            est_height = float(results["{0}_{1}".format(estimator_prefix, height_key_suffix)][sim_index])
            ci_lower = float(results["eti_95_lower_{0}".format(height_key_suffix)][sim_index])
            ci_upper = float(results["eti_95_upper_{0}".format(height_key_suffix)][sim_index])
            ci_width = ci_upper - ci_lower
            abs_error = math.fabs(true_height - est_height)

            err_vals[num_shared].append(abs_error)
            ci_vals[num_shared].append(ci_width)

        bd_abs_error._values = err_vals
        bd_ci_width._values = ci_vals
        bd_abs_error._labels = labels
        bd_ci_width._labels = labels
        return bd_abs_error, bd_ci_width


class HistogramData(object):
    def __init__(self, x = []):
        self._x = x

    @classmethod
    def init(cls, results, parameters, parameter_is_discrete):
        d = cls()
        d._x = []
        for parameter_str in parameters:
            if parameter_is_discrete:
                d._x.extend(int(x) for x in results["{0}".format(parameter_str)])
            else:
                d._x.extend(float(x) for x in results["{0}".format(parameter_str)])
        return d

    def _get_x(self):
        return self._x
    x = property(_get_x)


class ScatterData(object):
    def __init__(self,
            x = [],
            y = [],
            y_lower = [],
            y_upper = [],
            highlight_values = [],
            highlight_threshold = None,
            highlight_greater_than = True):
        self._x = x
        self._y = y
        self._y_lower = y_lower
        self._y_upper = y_upper
        self._highlight_values = highlight_values
        self._highlight_threshold = highlight_threshold
        self._highlight_greater_than = highlight_greater_than
        self._vet_data()
        self._highlight_indices = []
        self._populate_highlight_indices()
        self.highlight_color = (184 / 255.0, 90 / 255.0, 13 / 255.0) # pauburn

    @classmethod
    def init(cls, results, parameters,
            highlight_parameter_prefix = None,
            highlight_threshold = None,
            highlight_greater_than = True):
        d = cls()
        d._x = []
        d._y = []
        d._y_lower = []
        d._y_upper = []
        d._highlight_threshold = highlight_threshold
        d._highlight_values = []
        d._highlight_indices = []
        for parameter_str in parameters:
            d._x.extend(float(x) for x in results["true_{0}".format(parameter_str)])
            d._y.extend(float(x) for x in results["mean_{0}".format(parameter_str)])
            d._y_lower.extend(float(x) for x in results["eti_95_lower_{0}".format(parameter_str)])
            d._y_upper.extend(float(x) for x in results["eti_95_upper_{0}".format(parameter_str)])
            if highlight_parameter_prefix:
                d._highlight_values.extend(float(x) for x in results["{0}_{1}".format(
                        highlight_parameter_prefix,
                        parameter_str)])
        d._vet_data()
        d._populate_highlight_indices()
        return d

    @classmethod
    def init_time_v_sharing(cls, results,
            estimator_prefix = "mean",
            highlight_parameter_prefix = "psrf",
            highlight_threshold = 1.2,
            highlight_greater_than = True):
        d_abs_error = cls()
        d_abs_error._x = []
        d_abs_error._y = []
        d_abs_error._y_lower = []
        d_abs_error._y_upper = []
        d_abs_error._highlight_threshold = highlight_threshold
        d_abs_error._highlight_values = []
        d_abs_error._highlight_indices = []
        d_ci_width = cls()
        d_ci_width._x = []
        d_ci_width._y = []
        d_ci_width._y_lower = []
        d_ci_width._y_upper = []
        d_ci_width._highlight_threshold = highlight_threshold
        d_ci_width._highlight_values = []
        d_ci_width._highlight_indices = []
        nreplicates = len(results["true_model"])
        ncomparisons = len(results["true_model"][0])
        for sim_index in range(nreplicates):
            true_model = [int(x) for x in list(results["true_model"][sim_index])]
            assert len(true_model) == ncomparisons
            # Only use first comparison so as not to multply count the same
            # parameter estimates
            comparison_index = 0
            # for comparison_index in range(ncomparisons):
            num_shared = true_model.count(true_model[comparison_index])
            height_key_suffix = "root_height_c{0}sp1".format(comparison_index + 1)
            true_height = float(results["true_{0}".format(height_key_suffix)][sim_index])
            est_height = float(results["{0}_{1}".format(estimator_prefix, height_key_suffix)][sim_index])
            ci_lower = float(results["eti_95_lower_{0}".format(height_key_suffix)][sim_index])
            ci_upper = float(results["eti_95_upper_{0}".format(height_key_suffix)][sim_index])
            ci_width = ci_upper - ci_lower
            abs_error = math.fabs(true_height - est_height)
            d_abs_error._x.append(num_shared)
            d_ci_width._x.append(num_shared)
            d_abs_error._y.append(abs_error)
            d_ci_width._y.append(ci_width)
            if highlight_parameter_prefix:
                d_abs_error._highlight_values.append(float(results["{0}_{1}".format(
                        highlight_parameter_prefix,
                        height_key_suffix)][sim_index]))
                d_ci_width._highlight_values.append(float(results["{0}_{1}".format(
                        highlight_parameter_prefix,
                        height_key_suffix)][sim_index]))
        d_abs_error._vet_data()
        d_ci_width._vet_data()
        d_abs_error._populate_highlight_indices()
        d_ci_width._populate_highlight_indices()
        return d_abs_error, d_ci_width

    def _vet_data(self):
        assert len(self._x) == len(self._y)
        if self._y_lower:
            assert len(self._x) == len(self._y_lower)
        if self._y_upper:
            assert len(self._x) == len(self._y_upper)
        if self._highlight_values:
            assert len(self._x) == len(self._highlight_values)

    def _populate_highlight_indices(self):
        if (self._highlight_values) and (self._highlight_threshold is not None):
            for i in range(len(self._x)):
                if self.highlight(i):
                    self._highlight_indices.append(i)

    def has_y_ci(self):
        return bool(self.y_lower) and bool(self.y_upper)

    def has_highlights(self):
        return bool(self._highlight_indices)

    def _get_x(self):
        return self._x
    x = property(_get_x)

    def _get_y(self):
        return self._y
    y = property(_get_y)

    def _get_y_lower(self):
        return self._y_lower
    y_lower = property(_get_y_lower)

    def _get_y_upper(self):
        return self._y_upper
    y_upper = property(_get_y_upper)

    def _get_highlight_indices(self):
        return self._highlight_indices
    highlight_indices = property(_get_highlight_indices)

    def _get_highlight_x(self):
        return [self._x[i] for i in self._highlight_indices]
    highlight_x = property(_get_highlight_x)

    def _get_highlight_y(self):
        return [self._y[i] for i in self._highlight_indices]
    highlight_y = property(_get_highlight_y)

    def _get_highlight_y_lower(self):
        return [self._y_lower[i] for i in self._highlight_indices]
    highlight_y_lower = property(_get_highlight_y_lower)

    def _get_highlight_y_upper(self):
        return [self._y_upper[i] for i in self._highlight_indices]
    highlight_y_upper = property(_get_highlight_y_upper)

    def highlight(self, index):
        if (not self._highlight_values) or (self._highlight_threshold is None):
            return False

        if self._highlight_greater_than:
            if self._highlight_values[index] > self._highlight_threshold:
                return True
            else:
                return False
        else:
            if self._highlight_values[index] < self._highlight_threshold:
                return True
            else:
                return False
        return False


def get_abs_error(true, estimate):
    return math.fabs(true - estimate)

def get_relative_estimate(true, estimate):
    return estimate / float(true)

def get_relative_error(true, estimate):
    return math.fabs(true - estimate) / true

def get_coal_units_vs_error(results, height_parameters, size_parameters,
        error_func = get_relative_error,
        psrf_as_response = False):
    assert len(height_parameters) == len(size_parameters)
    x = []
    y = []
    psrf = []
    for sp_index in range(len(height_parameters)):
        # Ensure we are getting the height and size for the same population
        assert height_parameters[sp_index].split("_")[-1] == size_parameters[sp_index].split("_")[-1]
        true_height_key = "true_{0}".format(height_parameters[sp_index])
        true_size_key = "true_{0}".format(size_parameters[sp_index])
        mean_height_key = "mean_{0}".format(height_parameters[sp_index])
        psrf_height_key = "psrf_{0}".format(height_parameters[sp_index])
        psrf_size_key = "psrf_{0}".format(size_parameters[sp_index])
        nsims = len(results[true_height_key])
        assert nsims == len(results[true_size_key])
        assert nsims == len(results[mean_height_key])
        for sim_index in range(nsims):
            t = float(results[true_height_key][sim_index])
            t_mean = float(results[mean_height_key][sim_index])
            err = error_func(t, t_mean)
            N = float(results[true_size_key][sim_index])
            coal_units = t / (2.0 * N)
            x.append(coal_units)
            y.append(err)
            psrf_value = max(
                    float(results[psrf_height_key][sim_index]), 
                    float(results[psrf_size_key][sim_index])) 
            psrf.append(psrf_value)
    if psrf_as_response:
        return x, psrf
    return x, y, psrf

def get_true_est_coal_units(results, height_parameters, size_parameters):
    assert len(height_parameters) == len(size_parameters)
    x = []
    y = []
    psrf = []
    for sp_index in range(len(height_parameters)):
        # Ensure we are getting the height and size for the same population
        assert height_parameters[sp_index].split("_")[-1] == size_parameters[sp_index].split("_")[-1]
        true_height_key = "true_{0}".format(height_parameters[sp_index])
        true_size_key = "true_{0}".format(size_parameters[sp_index])
        mean_height_key = "mean_{0}".format(height_parameters[sp_index])
        mean_size_key = "mean_{0}".format(size_parameters[sp_index])
        psrf_height_key = "psrf_{0}".format(height_parameters[sp_index])
        psrf_size_key = "psrf_{0}".format(size_parameters[sp_index])
        nsims = len(results[true_height_key])
        assert nsims == len(results[true_size_key])
        assert nsims == len(results[mean_height_key])
        assert nsims == len(results[mean_size_key])
        for sim_index in range(nsims):
            t = float(results[true_height_key][sim_index])
            t_mean = float(results[mean_height_key][sim_index])
            N = float(results[true_size_key][sim_index])
            N_mean = float(results[mean_size_key][sim_index])
            coal_units = t / (2.0 * N)
            est_coal_units = t_mean / (2.0 * N_mean)
            x.append(coal_units)
            y.append(est_coal_units)
            psrf_value = max(
                    float(results[psrf_height_key][sim_index]), 
                    float(results[psrf_size_key][sim_index])) 
            psrf.append(psrf_value)
    return x, y, psrf


def plot_gamma(shape = 1.0,
        scale = 1.0,
        offset = 0.0,
        x_min = 0.0,
        x_max = None,
        number_of_points = 1000,
        x_label = "Relative root size",
        include_x_label = True,
        include_y_label = True,
        include_title = True,
        curve_line_width = 1.5,
        curve_line_style = '-',
        curve_line_color = (57 / 255.0, 115 / 255.0, 124 / 255.0),
        one_line_width = 1.0,
        one_line_style = '--',
        one_line_color = (184 / 255.0, 90 / 255.0, 13 / 255.0),
        plot_width = 3.5,
        plot_height = 3.0,
        xy_label_size = 16.0,
        title_size = 16.0,
        pad_left = 0.2,
        pad_right = 0.99,
        pad_bottom = 0.18,
        pad_top = 0.9,
        plot_file_prefix = None,
        plot_dir = project_util.PLOT_DIR,
        ):
    if x_max is None:
        x_max = scipy.stats.gamma.ppf(0.999, shape, scale = scale)
        x_max_plot = x_max + offset
    else:
        x_max_plot = x_max
    x = numpy.linspace(x_min, x_max, number_of_points)
    d = scipy.stats.gamma.pdf(x, shape, scale = scale)

    x_plot = [v + offset for v in x]

    if offset > 0.0:
        x_gap = numpy.linspace(x_min, offset, 100)
        d_gap = [0.0 for i in range(100)]
        x_plot = list(x_gap) + list(x_plot)
        d = d_gap + list(d)

    plt.close('all')
    fig = plt.figure(figsize = (plot_width, plot_height))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)
    ax = plt.subplot(gs[0, 0])
    line = ax.plot(x_plot, d)
    ax.set_xlim(x_min, x_max_plot)
    plt.setp(line,
            color = curve_line_color,
            linestyle = curve_line_style,
            linewidth = curve_line_width,
            marker = '',
            zorder = 100)
    ax.axvline(x = 1.0,
            color = one_line_color,
            linestyle = one_line_style,
            linewidth = one_line_width,
            marker = '',
            zorder = 0)
    if include_x_label:
        ax.set_xlabel(
                "{0}".format(x_label),
                fontsize = xy_label_size)
    if include_y_label:
        ax.set_ylabel(
                "Density",
                fontsize = xy_label_size)
    if include_title:
        # col_header = "$\\textrm{{\\sffamily Gamma}}({0:.0f}, \\textrm{{\\sffamily mean}} = {1:.1f})$".format(shape, (shape * scale) + offset)
        if offset == 0.0:
            col_header = "$\\textrm{{\\sffamily Gamma}}({0:.0f}, {1})$\n$\\textrm{{\\sffamily mean}} = {2}$".format(shape, scale, (shape * scale) + offset)
        else:
            col_header = "$\\textrm{{\\sffamily Gamma}}({0:.0f}, {1})$\n$\\textrm{{\\sffamily offset}} = {2}$ $\\textrm{{\\sffamily mean}} = {3}$".format(shape, scale, offset, (shape * scale) + offset)
        ax.set_title(col_header,
                fontsize = title_size)

    gs.update(
            left = pad_left,
            right = pad_right,
            bottom = pad_bottom,
            top = pad_top)

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    plot_path = os.path.join(plot_dir,
            "{0}-gamma.pdf".format(plot_file_prefix))
    plt.savefig(plot_path)
    _LOG.info("Plot written to {0!r}\n".format(plot_path))

def get_nevents_probs(
        results_paths,
        nevents = 1):
    probs = []
    prob_key = "num_events_{0}_p".format(nevents)
    for d in pycoevolity.parsing.spreadsheet_iter(results_paths):
        probs.append((
                float(d[prob_key]),
                int(int(d["true_num_events"]) == nevents)
                ))
    return probs

def bin_prob_correct_tuples(probability_correct_tuples, nbins = 20):
    bin_upper_limits = list(get_sequence_iter(0.0, 1.0, nbins+1))[1:]
    bin_width = (bin_upper_limits[1] - bin_upper_limits[0]) / 2.0
    bins = [[] for i in range(nbins)]
    n = 0
    for (p, t) in probability_correct_tuples:
        n += 1
        binned = False
        for i, l in enumerate(bin_upper_limits):
            if p < l:
                bins[i].append((p, t))
                binned = True
                break
        if not binned:
            bins[i].append((p, t))
    total = 0
    for b in bins:
        total += len(b)
    assert total == n
    assert len(bins) == nbins
    est_true_tups = []
    for i, b in enumerate(bins):
        ests = [p for (p, t) in b]
        est = sum(ests) / float(len(ests))
        correct = [t for (p, t) in b]
        true = sum(correct) / float(len(correct))
        est_true_tups.append((est, true))
    return bins, est_true_tups

def get_nevents_estimated_true_probs(
        results_paths,
        nevents = 1,
        nbins = 20):
    nevent_probs = get_nevents_probs(
            results_paths = results_paths,
            nevents = nevents)
    _LOG.info("\tparsed results for {0} simulations".format(len(nevent_probs)))
    bins, tups = bin_prob_correct_tuples(nevent_probs, nbins = nbins)
    _LOG.info("\tbin sample sizes: {0}".format(
            ", ".join(str(len(b)) for b in bins)
            ))
    return bins, tups

def plot_nevents_estimated_vs_true_probs(
        results_paths,
        nevents = 1,
        nbins = 20,
        plot_file_prefix = "",
        plot_dir = project_util.PLOT_DIR
        ):
    bins, est_true_probs = get_nevents_estimated_true_probs(
            results_paths = results_paths,
            nevents = nevents,
            nbins = nbins)

    plt.close('all')
    fig = plt.figure(figsize = (4.0, 3.5))
    ncols = 1
    nrows = 1
    gs = gridspec.GridSpec(nrows, ncols,
            wspace = 0.0,
            hspace = 0.0)

    ax = plt.subplot(gs[0, 0])
    x = [e for (e, t) in est_true_probs]
    y = [t for (e, t) in est_true_probs]
    sample_sizes = [len(b) for b in bins]
    line, = ax.plot(x, y)
    plt.setp(line,
            marker = 'o',
            markerfacecolor = 'none',
            markeredgecolor = '0.35',
            markeredgewidth = 0.7,
            markersize = 3.5,
            linestyle = '',
            zorder = 100,
            rasterized = False)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    for i, (label, lx, ly) in enumerate(zip(sample_sizes, x, y)):
        if i == 0:
            ax.annotate(
                    str(label),
                    xy = (lx, ly),
                    xytext = (1, 1),
                    textcoords = "offset points",
                    horizontalalignment = "left",
                    verticalalignment = "bottom")
        elif i == len(x) - 1:
            ax.annotate(
                    str(label),
                    xy = (lx, ly),
                    xytext = (-1, -1),
                    textcoords = "offset points",
                    horizontalalignment = "right",
                    verticalalignment = "top")
        else:
            ax.annotate(
                    str(label),
                    xy = (lx, ly),
                    xytext = (-1, 1),
                    textcoords = "offset points",
                    horizontalalignment = "right",
                    verticalalignment = "bottom")
    ylabel_text = ax.set_ylabel("True probability", size = 14.0)
    ax.text(0.5, -0.14,
            "Posterior probability of one divergence",
            horizontalalignment = "center",
            verticalalignment = "top",
            size = 14.0)
    identity_line, = ax.plot(
            [0.0, 1.0],
            [0.0, 1.0])
    plt.setp(identity_line,
            color = '0.8',
            linestyle = '-',
            linewidth = 1.0,
            marker = '',
            zorder = 0)

    gs.update(left = 0.10, right = 0.995, bottom = 0.18, top = 0.91)

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    plot_file_name = "est-vs-true-prob-nevent-1.pdf"
    if plot_file_prefix:
        plot_file_name = plot_file_prefix + "-" + plot_file_name
    plot_path = os.path.join(plot_dir,
            plot_file_name)
    plt.savefig(plot_path)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))

def get_sequence_iter(start = 0.0, stop = 1.0, n = 10):
    assert(stop > start)
    step = (stop - start) / float(n - 1)
    return ((start + (i * step)) for i in range(n))

def truncate_color_map(cmap, min_val = 0.0, max_val = 10, n = 100):
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(
                    n = cmap.name,
                    a = min_val,
                    b = max_val),
            cmap(list(get_sequence_iter(min_val, max_val, n))))
    return new_cmap

def get_root_gamma_parameters(root_shape_string, root_scale_string):
    shape = float(root_shape_string)
    denom = float("1" + ("0" * (len(root_scale_string) - 1)))
    scale = float(root_scale_string) / denom
    return shape, scale

root_gamma_pattern = re.compile(r'root-(?P<alpha_setting>\d+)-(?P<scale_setting>\d+)-')
def get_root_gamma_label(sim_dir):
    m = root_gamma_pattern.search(sim_dir)
    root_shape, root_scale = get_root_gamma_parameters(
            root_shape_string = m.group("alpha_setting"),
            root_scale_string = m.group("scale_setting"))
    mean_decimal_places = 1
    gamma_mean = root_shape * root_scale
    if gamma_mean < 0.49:
        mean_decimal_places = 2
    return "$\\textrm{{\\sffamily Gamma}}({shape}, \\textrm{{\\sffamily mean}} = {mean:.{mean_decimal_places}f})$".format(
            shape = int(root_shape),
            mean = gamma_mean,
            mean_decimal_places = mean_decimal_places)

def get_errors(values, lowers, uppers):
    n = len(values)
    assert(n == len(lowers))
    assert(n == len(uppers))
    return [[values[i] - lowers[i] for i in range(n)],
            [uppers[i] - values[i] for i in range(n)]]

def ci_width_iter(results, parameter_str):
    n = len(results["eti_95_upper_{0}".format(parameter_str)])
    for i in range(n):
        upper = float(results["eti_95_upper_{0}".format(parameter_str)][i])
        lower = float(results["eti_95_lower_{0}".format(parameter_str)][i])
        yield upper - lower

def absolute_error_iter(results, parameter_str):
    n = len(results["true_{0}".format(parameter_str)])
    for i in range(n):
        t = float(results["true_{0}".format(parameter_str)][i])
        e = float(results["mean_{0}".format(parameter_str)][i])
        yield math.fabs(t - e)


def plot_ess_versus_error(
        parameters,
        results_grid,
        column_labels = None,
        row_labels = None,
        parameter_label = "event time",
        plot_file_prefix = None,
        plot_dir = project_util.PLOT_DIR
        ):
    _LOG.info("Generating ESS vs CI scatter plots for {0}...".format(parameter_label))

    assert(len(parameters) == len(set(parameters)))
    if row_labels:
        assert len(row_labels) ==  len(results_grid)
    if column_labels:
        assert len(column_labels) == len(results_grid[0])

    nrows = len(results_grid)
    ncols = len(results_grid[0])

    if not plot_file_prefix:
        plot_file_prefix = parameters[0] 
    plot_file_prefix_ci = plot_file_prefix + "-ess-vs-ci-width"
    plot_file_prefix_error = plot_file_prefix + "-ess-vs-error"

    # Very inefficient, but parsing all results to get min/max for parameter
    ess_min = float('inf')
    ess_max = float('-inf')
    ci_width_min = float('inf')
    ci_width_max = float('-inf')
    error_min = float('inf')
    error_max = float('-inf')
    for row_index, results_grid_row in enumerate(results_grid):
        for column_index, results in enumerate(results_grid_row):
            for parameter_str in parameters:
                ci_widths = tuple(ci_width_iter(results, parameter_str))
                errors = tuple(absolute_error_iter(results, parameter_str))
                ess_min = min(ess_min,
                        min(float(x) for x in results["ess_sum_{0}".format(parameter_str)]))
                ess_max = max(ess_max,
                        max(float(x) for x in results["ess_sum_{0}".format(parameter_str)]))
                ci_width_min = min(ci_width_min, min(ci_widths))
                ci_width_max = max(ci_width_max, max(ci_widths))
                error_min = min(error_min, min(errors))
                error_max = max(error_max, max(errors))
    ess_axis_buffer = math.fabs(ess_max - ess_min) * 0.05
    ess_axis_min = ess_min - ess_axis_buffer
    ess_axis_max = ess_max + ess_axis_buffer
    ci_width_axis_buffer = math.fabs(ci_width_max - ci_width_min) * 0.05
    ci_width_axis_min = ci_width_min - ci_width_axis_buffer
    ci_width_axis_max = ci_width_max + ci_width_axis_buffer
    error_axis_buffer = math.fabs(error_max - error_min) * 0.05
    error_axis_min = error_min - error_axis_buffer
    error_axis_max = error_max + error_axis_buffer

    plt.close('all')
    w = 1.6
    h = 1.5
    fig_width = (ncols * w) + 1.0
    fig_height = (nrows * h) + 0.7
    fig = plt.figure(figsize = (fig_width, fig_height))
    gs = gridspec.GridSpec(nrows, ncols,
            wspace = 0.0,
            hspace = 0.0)

    for row_index, results_grid_row in enumerate(results_grid):
        for column_index, results in enumerate(results_grid_row):

            x = []
            y = []
            for parameter_str in parameters:
                x.extend(float(x) for x in results["ess_sum_{0}".format(parameter_str)])
                y.extend(ci_width_iter(results, parameter_str))

            assert(len(x) == len(y))
            ax = plt.subplot(gs[row_index, column_index])
            line, = ax.plot(x, y)
            plt.setp(line,
                    marker = 'o',
                    markerfacecolor = 'none',
                    markeredgecolor = '0.35',
                    markeredgewidth = 0.7,
                    markersize = 2.5,
                    linestyle = '',
                    zorder = 100,
                    rasterized = True)
            ax.set_xlim(ess_axis_min, ess_axis_max)
            ax.set_ylim(ci_width_axis_min, ci_width_axis_max)
            if column_labels and (row_index == 0):
                col_header = column_labels[column_index]
                ax.text(0.5, 1.015,
                        col_header,
                        horizontalalignment = "center",
                        verticalalignment = "bottom",
                        transform = ax.transAxes)
            if row_labels and (column_index == (ncols - 1)):
                row_label = row_labels[row_index]
                ax.text(1.015, 0.5,
                        row_label,
                        horizontalalignment = "left",
                        verticalalignment = "center",
                        rotation = 270.0,
                        transform = ax.transAxes)

    # show only the outside ticks
    all_axes = fig.get_axes()
    for ax in all_axes:
        if not ax.is_last_row():
            ax.set_xticks([])
        if not ax.is_first_col():
            ax.set_yticks([])

    # show tick labels only for lower-left plot 
    all_axes = fig.get_axes()
    for ax in all_axes:
        if ax.is_last_row() and ax.is_first_col():
            continue
        xtick_labels = ["" for item in ax.get_xticklabels()]
        ytick_labels = ["" for item in ax.get_yticklabels()]
        ax.set_xticklabels(xtick_labels)
        ax.set_yticklabels(ytick_labels)

    # avoid doubled spines
    all_axes = fig.get_axes()
    for ax in all_axes:
        for sp in ax.spines.values():
            sp.set_visible(False)
            sp.set_linewidth(2)
        if ax.is_first_row():
            ax.spines['top'].set_visible(True)
            ax.spines['bottom'].set_visible(True)
        else:
            ax.spines['bottom'].set_visible(True)
        if ax.is_first_col():
            ax.spines['left'].set_visible(True)
            ax.spines['right'].set_visible(True)
        else:
            ax.spines['right'].set_visible(True)

    fig.text(0.5, 0.001,
            "Effective sample size of {0}".format(parameter_label),
            horizontalalignment = "center",
            verticalalignment = "bottom",
            size = 18.0)
    fig.text(0.005, 0.5,
            "CI width {0}".format(parameter_label),
            horizontalalignment = "left",
            verticalalignment = "center",
            rotation = "vertical",
            size = 18.0)

    gs.update(left = 0.08, right = 0.98, bottom = 0.08, top = 0.97)

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    plot_path = os.path.join(plot_dir,
            "{0}-scatter.pdf".format(plot_file_prefix_ci))
    plt.savefig(plot_path, dpi=600)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))


    _LOG.info("Generating ESS vs error scatter plots for {0}...".format(parameter_label))
    plt.close('all')
    w = 1.6
    h = 1.5
    fig_width = (ncols * w) + 1.0
    fig_height = (nrows * h) + 0.7
    fig = plt.figure(figsize = (fig_width, fig_height))
    gs = gridspec.GridSpec(nrows, ncols,
            wspace = 0.0,
            hspace = 0.0)

    for row_index, results_grid_row in enumerate(results_grid):
        for column_index, results in enumerate(results_grid_row):
            x = []
            y = []
            for parameter_str in parameters:
                x.extend(float(x) for x in results["ess_sum_{0}".format(parameter_str)])
                y.extend(absolute_error_iter(results, parameter_str))
                

            assert(len(x) == len(y))
            ax = plt.subplot(gs[row_index, column_index])
            line, = ax.plot(x, y)
            plt.setp(line,
                    marker = 'o',
                    markerfacecolor = 'none',
                    markeredgecolor = '0.35',
                    markeredgewidth = 0.7,
                    markersize = 2.5,
                    linestyle = '',
                    zorder = 100,
                    rasterized = True)
            ax.set_xlim(ess_axis_min, ess_axis_max)
            ax.set_ylim(error_axis_min, error_axis_max)
            if column_labels and (row_index == 0):
                col_header = column_labels[column_index]
                ax.text(0.5, 1.015,
                        col_header,
                        horizontalalignment = "center",
                        verticalalignment = "bottom",
                        transform = ax.transAxes)
            if row_labels and (column_index == (ncols - 1)):
                row_label = row_labels[row_index]
                ax.text(1.015, 0.5,
                        row_label,
                        horizontalalignment = "left",
                        verticalalignment = "center",
                        rotation = 270.0,
                        transform = ax.transAxes)

    # show only the outside ticks
    all_axes = fig.get_axes()
    for ax in all_axes:
        if not ax.is_last_row():
            ax.set_xticks([])
        if not ax.is_first_col():
            ax.set_yticks([])

    # show tick labels only for lower-left plot 
    all_axes = fig.get_axes()
    for ax in all_axes:
        if ax.is_last_row() and ax.is_first_col():
            continue
        xtick_labels = ["" for item in ax.get_xticklabels()]
        ytick_labels = ["" for item in ax.get_yticklabels()]
        ax.set_xticklabels(xtick_labels)
        ax.set_yticklabels(ytick_labels)

    # avoid doubled spines
    all_axes = fig.get_axes()
    for ax in all_axes:
        for sp in ax.spines.values():
            sp.set_visible(False)
            sp.set_linewidth(2)
        if ax.is_first_row():
            ax.spines['top'].set_visible(True)
            ax.spines['bottom'].set_visible(True)
        else:
            ax.spines['bottom'].set_visible(True)
        if ax.is_first_col():
            ax.spines['left'].set_visible(True)
            ax.spines['right'].set_visible(True)
        else:
            ax.spines['right'].set_visible(True)

    fig.text(0.5, 0.001,
            "Effective sample size of {0}".format(parameter_label),
            horizontalalignment = "center",
            verticalalignment = "bottom",
            size = 18.0)
    fig.text(0.005, 0.5,
            "Absolute error of {0}".format(parameter_label),
            horizontalalignment = "left",
            verticalalignment = "center",
            rotation = "vertical",
            size = 18.0)

    gs.update(left = 0.08, right = 0.98, bottom = 0.08, top = 0.97)

    plot_path = os.path.join(plot_dir,
            "{0}-scatter.pdf".format(plot_file_prefix_error))
    plt.savefig(plot_path, dpi=600)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))


def generate_scatter_plots(
        data_grid,
        plot_file_prefix,
        parameter_symbol = "t",
        column_labels = None,
        row_labels = None,
        plot_width = 1.9,
        plot_height = 1.8,
        pad_left = 0.1,
        pad_right = 0.98,
        pad_bottom = 0.12,
        pad_top = 0.92,
        x_label = None,
        x_label_size = 18.0,
        y_label = None,
        y_label_size = 18.0,
        force_shared_x_range = True,
        force_shared_y_range = True,
        force_shared_xy_ranges = True,
        force_shared_spines = True,
        include_coverage = True,
        include_rmse = True,
        include_identity_line = True,
        include_error_bars = True,
        plot_dir = project_util.PLOT_DIR
        ):
    if force_shared_spines or force_shared_xy_ranges:
        force_shared_x_range = True
        force_shared_y_range = True

    if row_labels:
        assert len(row_labels) ==  len(data_grid)
    if column_labels:
        assert len(column_labels) == len(data_grid[0])

    nrows = len(data_grid)
    ncols = len(data_grid[0])

    x_min = float('inf')
    x_max = float('-inf')
    y_min = float('inf')
    y_max = float('-inf')
    for row_index, data_grid_row in enumerate(data_grid):
        for column_index, data in enumerate(data_grid_row):
            x_min = min(x_min, min(data.x))
            x_max = max(x_max, max(data.x))
            y_min = min(y_min, min(data.y))
            y_max = max(y_max, max(data.y))
    if force_shared_xy_ranges:
        mn = min(x_min, y_min)
        mx = max(x_max, y_max)
        x_min = mn
        y_min = mn
        x_max = mx
        y_max = mx
    x_buffer = math.fabs(x_max - x_min) * 0.05
    x_axis_min = x_min - x_buffer
    x_axis_max = x_max + x_buffer
    y_buffer = math.fabs(y_max - y_min) * 0.05
    y_axis_min = y_min - y_buffer
    y_axis_max = y_max + y_buffer


    plt.close('all')
    w = plot_width
    h = plot_height
    fig_width = (ncols * w)
    fig_height = (nrows * h)
    fig = plt.figure(figsize = (fig_width, fig_height))
    if force_shared_spines:
        gs = gridspec.GridSpec(nrows, ncols,
                wspace = 0.0,
                hspace = 0.0)
    else:
        gs = gridspec.GridSpec(nrows, ncols)

    for row_index, data_grid_row in enumerate(data_grid):
        for column_index, data in enumerate(data_grid_row):
            proportion_within_ci = 0.0
            if include_coverage and data.has_y_ci():
                proportion_within_ci = pycoevolity.stats.get_proportion_of_values_within_intervals(
                        data.x,
                        data.y_lower,
                        data.y_upper)
            rmse = 0.0
            if include_rmse:
                rmse = pycoevolity.stats.root_mean_square_error(data.x, data.y)
            ax = plt.subplot(gs[row_index, column_index])
            if include_error_bars and data.has_y_ci():
                line = ax.errorbar(
                        x = data.x,
                        y = data.y,
                        yerr = get_errors(data.y, data.y_lower, data.y_upper),
                        ecolor = '0.65',
                        elinewidth = 0.5,
                        capsize = 0.8,
                        barsabove = False,
                        marker = 'o',
                        linestyle = '',
                        markerfacecolor = 'none',
                        markeredgecolor = '0.35',
                        markeredgewidth = 0.7,
                        markersize = 2.5,
                        zorder = 100,
                        rasterized = True)
                if data.has_highlights():
                    # line = ax.errorbar(
                    #         x = data.highlight_x,
                    #         y = data.highlight_y,
                    #         yerr = get_errors(data.highlight_y,
                    #                 data.highlight_y_lower,
                    #                 data.highlight_y_upper),
                    #         ecolor = data.highlight_color,
                    #         elinewidth = 0.5,
                    #         capsize = 0.8,
                    #         barsabove = False,
                    #         marker = 'o',
                    #         linestyle = '',
                    #         markerfacecolor = 'none',
                    #         markeredgecolor = data.highlight_color,
                    #         markeredgewidth = 0.7,
                    #         markersize = 2.5,
                    #         zorder = 200,
                    #         rasterized = True)
                    line, = ax.plot(data.highlight_x, data.highlight_y)
                    plt.setp(line,
                            marker = 'o',
                            linestyle = '',
                            markerfacecolor = 'none',
                            markeredgecolor = data.highlight_color,
                            markeredgewidth = 0.7,
                            markersize = 2.5,
                            zorder = 200,
                            rasterized = True)
            else:
                line, = ax.plot(data.x, data.y)
                plt.setp(line,
                        marker = 'o',
                        linestyle = '',
                        markerfacecolor = 'none',
                        markeredgecolor = '0.35',
                        markeredgewidth = 0.7,
                        markersize = 2.5,
                        zorder = 100,
                        rasterized = True)
                if data.has_highlights():
                    line, = ax.plot(data.highlight_x, data.highlight_y)
                    plt.setp(line,
                            marker = 'o',
                            linestyle = '',
                            markerfacecolor = 'none',
                            markeredgecolor = data.highlight_color,
                            markeredgewidth = 0.7,
                            markersize = 2.5,
                            zorder = 200,
                            rasterized = True)
            if force_shared_x_range:
                ax.set_xlim(x_axis_min, x_axis_max)
            else:
                ax.set_xlim(min(data.x), max(data.x))
            if force_shared_y_range:
                ax.set_ylim(y_axis_min, y_axis_max)
            else:
                ax.set_ylim(min(data.y), max(data.y))
            if include_identity_line:
                identity_line, = ax.plot(
                        [x_axis_min, x_axis_max],
                        [y_axis_min, y_axis_max])
                plt.setp(identity_line,
                        color = '0.7',
                        linestyle = '-',
                        linewidth = 1.0,
                        marker = '',
                        zorder = 0)
            if include_coverage:
                ax.text(0.02, 0.97,
                        "\\scriptsize\\noindent$p({0:s} \\in \\textrm{{\\sffamily CI}}) = {1:.3f}$".format(
                                parameter_symbol,
                                proportion_within_ci),
                        horizontalalignment = "left",
                        verticalalignment = "top",
                        transform = ax.transAxes,
                        size = 6.0,
                        zorder = 300)
            if include_rmse:
                text_y = 0.97
                if include_coverage:
                    text_y = 0.87
                ax.text(0.02, text_y,
                        "\\scriptsize\\noindent RMSE = {0:.2e}".format(
                                rmse),
                        horizontalalignment = "left",
                        verticalalignment = "top",
                        transform = ax.transAxes,
                        size = 6.0,
                        zorder = 300)
            if column_labels and (row_index == 0):
                col_header = column_labels[column_index]
                ax.text(0.5, 1.015,
                        col_header,
                        horizontalalignment = "center",
                        verticalalignment = "bottom",
                        transform = ax.transAxes)
            if row_labels and (column_index == (ncols - 1)):
                row_label = row_labels[row_index]
                ax.text(1.015, 0.5,
                        row_label,
                        horizontalalignment = "left",
                        verticalalignment = "center",
                        rotation = 270.0,
                        transform = ax.transAxes)

    if force_shared_spines:
        # show only the outside ticks
        all_axes = fig.get_axes()
        for ax in all_axes:
            if not ax.is_last_row():
                ax.set_xticks([])
            if not ax.is_first_col():
                ax.set_yticks([])

        # show tick labels only for lower-left plot 
        all_axes = fig.get_axes()
        for ax in all_axes:
            if ax.is_last_row() and ax.is_first_col():
                continue
            xtick_labels = ["" for item in ax.get_xticklabels()]
            ytick_labels = ["" for item in ax.get_yticklabels()]
            ax.set_xticklabels(xtick_labels)
            ax.set_yticklabels(ytick_labels)

        # avoid doubled spines
        all_axes = fig.get_axes()
        for ax in all_axes:
            for sp in ax.spines.values():
                sp.set_visible(False)
                sp.set_linewidth(2)
            if ax.is_first_row():
                ax.spines['top'].set_visible(True)
                ax.spines['bottom'].set_visible(True)
            else:
                ax.spines['bottom'].set_visible(True)
            if ax.is_first_col():
                ax.spines['left'].set_visible(True)
                ax.spines['right'].set_visible(True)
            else:
                ax.spines['right'].set_visible(True)

    if x_label:
        fig.text(0.5, 0.001,
                x_label,
                horizontalalignment = "center",
                verticalalignment = "bottom",
                size = x_label_size)
    if y_label:
        fig.text(0.005, 0.5,
                y_label,
                horizontalalignment = "left",
                verticalalignment = "center",
                rotation = "vertical",
                size = y_label_size)

    gs.update(left = pad_left,
            right = pad_right,
            bottom = pad_bottom,
            top = pad_top)

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    plot_path = os.path.join(plot_dir,
            "{0}-scatter.pdf".format(plot_file_prefix))
    plt.savefig(plot_path, dpi=600)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))

def generate_specific_scatter_plot(
        data,
        plot_file_prefix,
        parameter_symbol = "t",
        title = None,
        title_size = 16.0,
        x_label = None,
        x_label_size = 16.0,
        y_label = None,
        y_label_size = 16.0,
        plot_width = 3.5,
        plot_height = 3.0,
        pad_left = 0.2,
        pad_right = 0.99,
        pad_bottom = 0.18,
        pad_top = 0.9,
        force_shared_xy_ranges = True,
        xy_limits = None,
        include_coverage = True,
        include_rmse = True,
        include_identity_line = True,
        include_error_bars = True,
        plot_dir = project_util.PLOT_DIR):

    if xy_limits:
        x_axis_min, x_axis_max, y_axis_min, y_axis_max = xy_limits
    else:
        x_min = min(data.x)
        x_max = max(data.x)
        y_min = min(data.y)
        y_max = max(data.y)
        if force_shared_xy_ranges:
            mn = min(x_min, y_min)
            mx = max(x_max, y_max)
            x_min = mn
            y_min = mn
            x_max = mx
            y_max = mx
        x_buffer = math.fabs(x_max - x_min) * 0.05
        x_axis_min = x_min - x_buffer
        x_axis_max = x_max + x_buffer
        y_buffer = math.fabs(y_max - y_min) * 0.05
        y_axis_min = y_min - y_buffer
        y_axis_max = y_max + y_buffer

    plt.close('all')
    fig = plt.figure(figsize = (plot_width, plot_height))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)

    proportion_within_ci = 0.0
    if include_coverage and data.has_y_ci():
        proportion_within_ci = pycoevolity.stats.get_proportion_of_values_within_intervals(
                data.x,
                data.y_lower,
                data.y_upper)
    rmse = 0.0
    if include_rmse:
        rmse = pycoevolity.stats.root_mean_square_error(data.x, data.y)
    ax = plt.subplot(gs[0, 0])
    if include_error_bars and data.has_y_ci():
        line = ax.errorbar(
                x = data.x,
                y = data.y,
                yerr = get_errors(data.y, data.y_lower, data.y_upper),
                ecolor = '0.65',
                elinewidth = 0.5,
                capsize = 0.8,
                barsabove = False,
                marker = 'o',
                linestyle = '',
                markerfacecolor = 'none',
                markeredgecolor = '0.35',
                markeredgewidth = 0.7,
                markersize = 2.5,
                zorder = 100,
                rasterized = True)
        if data.has_highlights():
            # line = ax.errorbar(
            #         x = data.highlight_x,
            #         y = data.highlight_y,
            #         yerr = get_errors(data.highlight_y,
            #                 data.highlight_y_lower,
            #                 data.highlight_y_upper),
            #         ecolor = data.highlight_color,
            #         elinewidth = 0.5,
            #         capsize = 0.8,
            #         barsabove = False,
            #         marker = 'o',
            #         linestyle = '',
            #         markerfacecolor = 'none',
            #         markeredgecolor = data.highlight_color,
            #         markeredgewidth = 0.7,
            #         markersize = 2.5,
            #         zorder = 200,
            #         rasterized = True)
            line, = ax.plot(data.highlight_x, data.highlight_y)
            plt.setp(line,
                    marker = 'o',
                    linestyle = '',
                    markerfacecolor = 'none',
                    markeredgecolor = data.highlight_color,
                    markeredgewidth = 0.7,
                    markersize = 2.5,
                    zorder = 200,
                    rasterized = True)
    else:
        line, = ax.plot(data.x, data.y)
        plt.setp(line,
                marker = 'o',
                linestyle = '',
                markerfacecolor = 'none',
                markeredgecolor = '0.35',
                markeredgewidth = 0.7,
                markersize = 2.5,
                zorder = 100,
                rasterized = True)
        if data.has_highlights():
            line, = ax.plot(x = data.highlight_x, y = data.highlight_y)
            plt.setp(line,
                    marker = 'o',
                    linestyle = '',
                    markerfacecolor = 'none',
                    markeredgecolor = data.highlight_color,
                    markeredgewidth = 0.7,
                    markersize = 2.5,
                    zorder = 200,
                    rasterized = True)
    ax.set_xlim(x_axis_min, x_axis_max)
    ax.set_ylim(y_axis_min, y_axis_max)
    if include_identity_line:
        identity_line, = ax.plot(
                [x_axis_min, x_axis_max],
                [y_axis_min, y_axis_max])
        plt.setp(identity_line,
                color = '0.7',
                linestyle = '-',
                linewidth = 1.0,
                marker = '',
                zorder = 0)
    if include_coverage:
        ax.text(0.02, 0.97,
                "\\normalsize\\noindent$p({0:s} \\in \\textrm{{\\sffamily CI}}) = {1:.3f}$".format(
                        parameter_symbol,
                        proportion_within_ci),
                horizontalalignment = "left",
                verticalalignment = "top",
                transform = ax.transAxes,
                size = 8.0,
                zorder = 300)
    if include_rmse:
        text_y = 0.97
        if include_coverage:
            text_y = 0.87
        ax.text(0.02, text_y,
                "\\normalsize\\noindent RMSE = {0:.2e}".format(
                        rmse),
                horizontalalignment = "left",
                verticalalignment = "top",
                transform = ax.transAxes,
                size = 8.0,
                zorder = 300)
    if x_label is not None:
        ax.set_xlabel(
                x_label,
                fontsize = x_label_size)
    if y_label is not None:
        ax.set_ylabel(
                y_label,
                fontsize = y_label_size)
    if title is not None:
        ax.set_title(plot_title,
                fontsize = title_size)

    gs.update(
            left = pad_left,
            right = pad_right,
            bottom = pad_bottom,
            top = pad_top)

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    plot_path = os.path.join(plot_dir,
            "{0}-scatter.pdf".format(plot_file_prefix))
    plt.savefig(plot_path, dpi=600)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))


def generate_specific_box_plot(
        data,
        plot_file_prefix,
        title = None,
        title_size = 16.0,
        x_label = None,
        x_label_size = 16.0,
        y_label = None,
        y_label_size = 16.0,
        plot_width = 3.5,
        plot_height = 3.0,
        pad_left = 0.2,
        pad_right = 0.99,
        pad_bottom = 0.18,
        pad_top = 0.9,
        jitter = 0.01,
        alpha = 0.4,
        rasterized = False,
        plot_dir = project_util.PLOT_DIR):

    plt.close('all')
    fig = plt.figure(figsize = (plot_width, plot_height))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)

    ax = plt.subplot(gs[0, 0])
    box_dict = ax.boxplot(data.values,
            labels = data.labels,
            notch = False,
            sym = '',
            vert = True,
            showfliers = False,
            whis = 'range',
            zorder = 500)
    # Change median line from default color to black
    plt.setp(box_dict["medians"], color = "black")
    for i in range(data.number_of_categories):
        vals = data.values[i]
        x = numpy.random.uniform(low = i + 1 - jitter, high = i + 1 + jitter, size = len(vals))
        ax.plot(x, vals,
                marker = 'o',
                linestyle = '',
                markerfacecolor = 'none',
                markeredgecolor = '0.35',
                markeredgewidth = 0.7,
                alpha = alpha,
                markersize = 2.5,
                zorder = 100,
                rasterized = rasterized)
    if x_label is not None:
        ax.set_xlabel(
                x_label,
                fontsize = x_label_size)
    if y_label is not None:
        ax.set_ylabel(
                y_label,
                fontsize = y_label_size)
    if title is not None:
        ax.set_title(plot_title,
                fontsize = title_size)

    gs.update(
            left = pad_left,
            right = pad_right,
            bottom = pad_bottom,
            top = pad_top)

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    plot_path = os.path.join(plot_dir,
            "{0}-box.pdf".format(plot_file_prefix))
    if rasterized:
        plt.savefig(plot_path, dpi=600)
    else:
        plt.savefig(plot_path)
    _LOG.info("Box plot written to {0!r}\n".format(plot_path))


def generate_histograms(
        data_grid,
        plot_file_prefix,
        column_labels = None,
        row_labels = None,
        parameter_label = "Number of variable sites",
        range_key = "range",
        number_of_digits = 0,
        plot_width = 1.9,
        plot_height = 1.8,
        pad_left = 0.1,
        pad_right = 0.98,
        pad_bottom = 0.12,
        pad_top = 0.92,
        force_shared_x_range = True,
        force_shared_bins = True,
        force_shared_y_range = True,
        force_shared_spines = True,
        plot_dir = project_util.PLOT_DIR
        ):
    if force_shared_spines:
        force_shared_x_range = True
        force_shared_y_range = True

    if row_labels:
        assert len(row_labels) ==  len(data_grid)
    if column_labels:
        assert len(column_labels) == len(data_grid[0])

    nrows = len(data_grid)
    ncols = len(data_grid[0])

    x_min = float('inf')
    x_max = float('-inf')
    for row_index, data_grid_row in enumerate(data_grid):
        for column_index, data in enumerate(data_grid_row):
            x_min = min(x_min, min(data.x))
            x_max = max(x_max, max(data.x))

    axis_buffer = math.fabs(x_max - x_min) * 0.05
    axis_min = x_min - axis_buffer
    axis_max = x_max + axis_buffer

    plt.close('all')
    w = plot_width
    h = plot_height
    fig_width = (ncols * w)
    fig_height = (nrows * h)
    fig = plt.figure(figsize = (fig_width, fig_height))
    if force_shared_spines:
        gs = gridspec.GridSpec(nrows, ncols,
                wspace = 0.0,
                hspace = 0.0)
    else:
        gs = gridspec.GridSpec(nrows, ncols)

    hist_bins = None
    x_range = None
    if force_shared_x_range:
        x_range = (x_min, x_max)
    for row_index, data_grid_row in enumerate(data_grid):
        for column_index, data in enumerate(data_grid_row):
            summary = pycoevolity.stats.get_summary(data.x)
            _LOG.info("0.025, 0.975 quantiles: {0:.2f}, {1:.2f}".format(
                    summary["qi_95"][0],
                    summary["qi_95"][1]))

            ax = plt.subplot(gs[row_index, column_index])
            n, bins, patches = ax.hist(data.x,
                    weights = [1.0 / float(len(data.x))] * len(data.x),
                    bins = hist_bins,
                    range = x_range,
                    cumulative = False,
                    histtype = 'bar',
                    align = 'mid',
                    orientation = 'vertical',
                    rwidth = None,
                    log = False,
                    color = None,
                    edgecolor = '0.5',
                    facecolor = '0.5',
                    fill = True,
                    hatch = None,
                    label = None,
                    linestyle = None,
                    linewidth = None,
                    zorder = 10,
                    )
            if (hist_bins is None) and force_shared_bins:
                hist_bins = bins
            ax.text(0.98, 0.98,
                    "\\scriptsize {mean:,.{ndigits}f} ({lower:,.{ndigits}f}--{upper:,.{ndigits}f})".format(
                            mean = summary["mean"],
                            lower = summary[range_key][0],
                            upper = summary[range_key][1],
                            ndigits = number_of_digits),
                    horizontalalignment = "right",
                    verticalalignment = "top",
                    transform = ax.transAxes,
                    zorder = 200)

            if column_labels and (row_index == 0):
                col_header = column_labels[column_index]
                ax.text(0.5, 1.015,
                        col_header,
                        horizontalalignment = "center",
                        verticalalignment = "bottom",
                        transform = ax.transAxes)
            if row_labels and (column_index == (ncols - 1)):
                row_label = row_labels[row_index]
                ax.text(1.015, 0.5,
                        row_label,
                        horizontalalignment = "left",
                        verticalalignment = "center",
                        rotation = 270.0,
                        transform = ax.transAxes)

    if force_shared_y_range:
        all_axes = fig.get_axes()
        # y_max = float('-inf')
        # for ax in all_axes:
        #     ymn, ymx = ax.get_ylim()
        #     y_max = max(y_max, ymx)
        for ax in all_axes:
            ax.set_ylim(0.0, 1.0)

    if force_shared_spines:
        # show only the outside ticks
        all_axes = fig.get_axes()
        for ax in all_axes:
            if not ax.is_last_row():
                ax.set_xticks([])
            if not ax.is_first_col():
                ax.set_yticks([])

        # show tick labels only for lower-left plot 
        all_axes = fig.get_axes()
        for ax in all_axes:
            if ax.is_last_row() and ax.is_first_col():
                continue
            xtick_labels = ["" for item in ax.get_xticklabels()]
            ytick_labels = ["" for item in ax.get_yticklabels()]
            ax.set_xticklabels(xtick_labels)
            ax.set_yticklabels(ytick_labels)

        # avoid doubled spines
        all_axes = fig.get_axes()
        for ax in all_axes:
            for sp in ax.spines.values():
                sp.set_visible(False)
                sp.set_linewidth(2)
            if ax.is_first_row():
                ax.spines['top'].set_visible(True)
                ax.spines['bottom'].set_visible(True)
            else:
                ax.spines['bottom'].set_visible(True)
            if ax.is_first_col():
                ax.spines['left'].set_visible(True)
                ax.spines['right'].set_visible(True)
            else:
                ax.spines['right'].set_visible(True)

    fig.text(0.5, 0.001,
            parameter_label,
            horizontalalignment = "center",
            verticalalignment = "bottom",
            size = 18.0)
    fig.text(0.005, 0.5,
            # "Density",
            "Frequency",
            horizontalalignment = "left",
            verticalalignment = "center",
            rotation = "vertical",
            size = 18.0)

    gs.update(left = pad_left,
            right = pad_right,
            bottom = pad_bottom,
            top = pad_top)

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    plot_path = os.path.join(plot_dir,
            "{0}-histograms.pdf".format(plot_file_prefix))
    plt.savefig(plot_path)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))


def generate_specific_histogram(
        data,
        plot_file_prefix,
        title = None,
        title_size = 16.0,
        x_label = None,
        x_label_size = 16.0,
        y_label = None,
        y_label_size = 16.0,
        plot_width = 3.5,
        plot_height = 3.0,
        pad_left = 0.2,
        pad_right = 0.99,
        pad_bottom = 0.18,
        pad_top = 0.9,
        bins = None,
        x_range = None,
        range_key = "range",
        center_key = "mean",
        number_of_digits = 0,
        plot_dir = project_util.PLOT_DIR
        ):

    plt.close('all')
    fig = plt.figure(figsize = (plot_width, plot_height))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)

    summary = pycoevolity.stats.get_summary(data.x)
    _LOG.info("0.025, 0.975 quantiles: {0:.2f}, {1:.2f}".format(
            summary["qi_95"][0],
            summary["qi_95"][1]))

    ax = plt.subplot(gs[0, 0])
    n, b, patches = ax.hist(data.x,
            weights = [1.0 / float(len(data.x))] * len(data.x),
            bins = bins,
            range = x_range,
            cumulative = False,
            histtype = 'bar',
            align = 'mid',
            orientation = 'vertical',
            rwidth = None,
            log = False,
            color = None,
            edgecolor = '0.5',
            facecolor = '0.5',
            fill = True,
            hatch = None,
            label = None,
            linestyle = None,
            linewidth = None,
            zorder = 10,
            )

    ax.text(0.98, 0.98,
            "{mean:,.{ndigits}f} ({lower:,.{ndigits}f}--{upper:,.{ndigits}f})".format(
                    mean = summary[center_key],
                    lower = summary[range_key][0],
                    upper = summary[range_key][1],
                    ndigits = number_of_digits),
            horizontalalignment = "right",
            verticalalignment = "top",
            transform = ax.transAxes,
            zorder = 200)


    if x_label is not None:
        ax.set_xlabel(
                x_label,
                fontsize = x_label_size)
    if y_label is not None:
        ax.set_ylabel(
                y_label,
                fontsize = y_label_size)
    if title is not None:
        ax.set_title(plot_title,
                fontsize = title_size)

    gs.update(
            left = pad_left,
            right = pad_right,
            bottom = pad_bottom,
            top = pad_top)

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    plot_path = os.path.join(plot_dir,
            "{0}-histogram.pdf".format(plot_file_prefix))
    plt.savefig(plot_path)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))


def generate_model_plots(
        results_grid,
        column_labels = None,
        row_labels = None,
        number_of_comparisons = 3,
        plot_width = 1.6,
        plot_height = 1.5,
        pad_left = 0.1,
        pad_right = 0.98,
        pad_bottom = 0.12,
        pad_top = 0.92,
        y_label_size = 18.0,
        y_label = None,
        number_font_size = 12.0,
        plot_file_prefix = None,
        plot_dir = project_util.PLOT_DIR
        ):
    _LOG.info("Generating model plots...")

    cmap = truncate_color_map(plt.cm.binary, 0.0, 0.65, 100)

    if row_labels:
        assert len(row_labels) ==  len(results_grid)
    if column_labels:
        assert len(column_labels) == len(results_grid[0])

    nrows = len(results_grid)
    ncols = len(results_grid[0])

    plt.close('all')
    w = plot_width
    h = plot_height
    fig_width = (ncols * w)
    fig_height = (nrows * h)
    fig = plt.figure(figsize = (fig_width, fig_height))
    gs = gridspec.GridSpec(nrows, ncols,
            wspace = 0.0,
            hspace = 0.0)

    for row_index, results_grid_row in enumerate(results_grid):
        for column_index, results in enumerate(results_grid_row):
            true_map_nevents = []
            true_map_nevents_probs = []
            for i in range(number_of_comparisons):
                true_map_nevents.append([0 for i in range(number_of_comparisons)])
                true_map_nevents_probs.append([[] for i in range(number_of_comparisons)])
            true_nevents = tuple(int(x) for x in results["true_num_events"])
            map_nevents = tuple(int(x) for x in results["map_num_events"])
            true_nevents_cred_levels = tuple(float(x) for x in results["true_num_events_cred_level"])
            # true_model_cred_levels = tuple(float(x) for x in results["true_model_cred_level"])
            assert(len(true_nevents) == len(map_nevents))
            assert(len(true_nevents) == len(true_nevents_cred_levels))
            # assert(len(true_nevents) == len(true_model_cred_levels))

            true_nevents_probs = []
            map_nevents_probs = []
            for i in range(len(true_nevents)):
                true_nevents_probs.append(float(
                    results["num_events_{0}_p".format(true_nevents[i])][i]))
                map_nevents_probs.append(float(
                    results["num_events_{0}_p".format(map_nevents[i])][i]))
            assert(len(true_nevents) == len(true_nevents_probs))
            assert(len(true_nevents) == len(map_nevents_probs))

            mean_true_nevents_prob = sum(true_nevents_probs) / len(true_nevents_probs)
            median_true_nevents_prob = pycoevolity.stats.median(true_nevents_probs)

            nevents_within_95_cred = 0
            # model_within_95_cred = 0
            ncorrect = 0
            for i in range(len(true_nevents)):
                true_map_nevents[map_nevents[i] - 1][true_nevents[i] - 1] += 1
                true_map_nevents_probs[map_nevents[i] - 1][true_nevents[i] - 1].append(map_nevents_probs[i])
                if true_nevents_cred_levels[i] <= 0.95:
                    nevents_within_95_cred += 1
                # if true_model_cred_levels[i] <= 0.95:
                #     model_within_95_cred += 1
                if true_nevents[i] == map_nevents[i]:
                    ncorrect += 1
            p_nevents_within_95_cred = nevents_within_95_cred / float(len(true_nevents))
            # p_model_within_95_cred = model_within_95_cred / float(len(true_nevents))
            p_correct = ncorrect / float(len(true_nevents))

            _LOG.info("p(nevents within CS) = {0:.4f}".format(p_nevents_within_95_cred))
            # _LOG.info("p(model within CS) = {0:.4f}".format(p_model_within_95_cred))
            ax = plt.subplot(gs[row_index, column_index])

            ax.imshow(true_map_nevents,
                    origin = 'lower',
                    cmap = cmap,
                    interpolation = 'none',
                    aspect = 'auto'
                    # extent = [0.5, 3.5, 0.5, 3.5]
                    )
            for i, row_list in enumerate(true_map_nevents):
                for j, num_events in enumerate(row_list):
                    ax.text(j, i,
                            str(num_events),
                            horizontalalignment = "center",
                            verticalalignment = "center",
                            size = number_font_size)
            ax.text(0.98, 0.02,
                    "\\scriptsize$p(k \\in \\textrm{{\\sffamily CS}}) = {0:.3f}$".format(
                            p_nevents_within_95_cred),
                    horizontalalignment = "right",
                    verticalalignment = "bottom",
                    transform = ax.transAxes)
            ax.text(0.02, 0.98,
                    "\\scriptsize$p(\\hat{{k}} = k) = {0:.3f}$".format(
                            p_correct),
                    horizontalalignment = "left",
                    verticalalignment = "top",
                    transform = ax.transAxes)
            ax.text(0.98, 0.98,
                    "\\scriptsize$\\widetilde{{p(k|\\mathbf{{D}})}} = {0:.3f}$".format(
                            median_true_nevents_prob),
                    horizontalalignment = "right",
                    verticalalignment = "top",
                    transform = ax.transAxes)
            if column_labels and (row_index == 0):
                col_header = column_labels[column_index]
                ax.text(0.5, 1.015,
                        col_header,
                        horizontalalignment = "center",
                        verticalalignment = "bottom",
                        transform = ax.transAxes)
            if row_labels and (column_index == (ncols - 1)):
                row_label = row_labels[row_index]
                ax.text(1.015, 0.5,
                        row_label,
                        horizontalalignment = "left",
                        verticalalignment = "center",
                        rotation = 270.0,
                        transform = ax.transAxes)

    # show only the outside ticks
    all_axes = fig.get_axes()
    for ax in all_axes:
        if not ax.is_last_row():
            ax.set_xticks([])
        if not ax.is_first_col():
            ax.set_yticks([])

    # show tick labels only for lower-left plot 
    all_axes = fig.get_axes()
    for ax in all_axes:
        # Make sure ticks correspond only with number of events
        ax.xaxis.set_ticks(range(number_of_comparisons))
        ax.yaxis.set_ticks(range(number_of_comparisons))
        if ax.is_last_row() and ax.is_first_col():
            xtick_labels = [item for item in ax.get_xticklabels()]
            for i in range(len(xtick_labels)):
                xtick_labels[i].set_text(str(i + 1))
            ytick_labels = [item for item in ax.get_yticklabels()]
            for i in range(len(ytick_labels)):
                ytick_labels[i].set_text(str(i + 1))
            ax.set_xticklabels(xtick_labels)
            ax.set_yticklabels(ytick_labels)
        else:
            xtick_labels = ["" for item in ax.get_xticklabels()]
            ytick_labels = ["" for item in ax.get_yticklabels()]
            ax.set_xticklabels(xtick_labels)
            ax.set_yticklabels(ytick_labels)

    # avoid doubled spines
    all_axes = fig.get_axes()
    for ax in all_axes:
        for sp in ax.spines.values():
            sp.set_visible(False)
            sp.set_linewidth(2)
        if ax.is_first_row():
            ax.spines['top'].set_visible(True)
            ax.spines['bottom'].set_visible(True)
        else:
            ax.spines['bottom'].set_visible(True)
        if ax.is_first_col():
            ax.spines['left'].set_visible(True)
            ax.spines['right'].set_visible(True)
        else:
            ax.spines['right'].set_visible(True)

    fig.text(0.5, 0.001,
            "True number of events ($k$)",
            horizontalalignment = "center",
            verticalalignment = "bottom",
            size = 18.0)
    if y_label is None:
        y_label = "Estimated number of events ($\\hat{{k}}$)"
    fig.text(0.005, 0.5,
            y_label,
            horizontalalignment = "left",
            verticalalignment = "center",
            rotation = "vertical",
            size = 18.0)

    gs.update(left = pad_left,
            right = pad_right,
            bottom = pad_bottom,
            top = pad_top)

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    if plot_file_prefix:
        plot_path = os.path.join(plot_dir,
                "{0}-nevents.pdf".format(plot_file_prefix))
    else:
        plot_path = os.path.join(plot_dir,
                "nevents.pdf")
    plt.savefig(plot_path)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))

def generate_specific_model_plots(
        results,
        number_of_comparisons = 3,
        show_all_models = False,
        plot_title = None,
        include_x_label = True,
        include_y_label = True,
        include_median = True,
        include_cs = True,
        include_prop_correct = True,
        plot_width = 3.5,
        plot_height = 3.0,
        xy_label_size = 16.0,
        title_size = 16.0,
        pad_left = 0.2,
        pad_right = 0.99,
        pad_bottom = 0.18,
        pad_top = 0.9,
        lower_annotation_y = 0.02,
        upper_annotation_y = 0.92,
        plot_file_prefix = None,
        plot_dir = project_util.PLOT_DIR,
        model_key = "model",
        num_events_key = "num_events",
        ):
    if show_all_models and (number_of_comparisons != 3):
        raise Exception("show all models only supported for 3 comparisons")
    _LOG.info("Generating model plots...")

    cmap = truncate_color_map(plt.cm.binary, 0.0, 0.65, 100)

    model_to_index = {
            "000": 0,
            "001": 1,
            "010": 2,
            "011": 3,
            "012": 4,
            }
    index_to_model = {}
    for k, v in model_to_index.items():
        index_to_model[v] = k

    plt.close('all')
    fig = plt.figure(figsize = (plot_width, plot_height))
    gs = gridspec.GridSpec(1, 1,
            wspace = 0.0,
            hspace = 0.0)

    true_map_nevents = []
    true_map_model = []
    true_map_nevents_probs = []
    for i in range(number_of_comparisons):
        true_map_nevents.append([0 for i in range(number_of_comparisons)])
        true_map_nevents_probs.append([[] for i in range(number_of_comparisons)])
    for i in range(5):
        true_map_model.append([0 for i in range(5)])
        true_map_nevents_probs.append([[] for i in range(5)])
    true_nevents = tuple(int(x) for x in results["true_{num_events}".format(num_events = num_events_key)])
    map_nevents = tuple(int(x) for x in results["map_{num_events}".format(num_events = num_events_key)])
    true_model = tuple(x for x in results["true_{model}".format(model = model_key)])
    map_model = tuple(x for x in results["map_{model}".format(model = model_key)])
    true_nevents_cred_levels = tuple(float(x) for x in results["true_{num_events}_cred_level".format(num_events = num_events_key)])
    true_model_cred_levels = tuple(float(x) for x in results["true_{model}_cred_level".format(model = model_key)])
    assert(len(true_nevents) == len(map_nevents))
    assert(len(true_nevents) == len(true_nevents_cred_levels))
    assert(len(true_nevents) == len(true_model_cred_levels))
    assert(len(true_nevents) == len(true_model))
    assert(len(true_nevents) == len(map_model))

    true_nevents_probs = []
    map_nevents_probs = []
    for i in range(len(true_nevents)):
        true_nevents_probs.append(float(
            results["{num_events}_{n}_p".format(num_events = num_events_key, n = true_nevents[i])][i]))
        map_nevents_probs.append(float(
            results["{num_events}_{n}_p".format(num_events = num_events_key, n = map_nevents[i])][i]))
    assert(len(true_nevents) == len(true_nevents_probs))
    assert(len(true_nevents) == len(map_nevents_probs))

    mean_true_nevents_prob = sum(true_nevents_probs) / len(true_nevents_probs)
    median_true_nevents_prob = pycoevolity.stats.median(true_nevents_probs)

    true_model_probs = tuple(float(x) for x in results["true_{model}_p".format(model = model_key)])
    assert(len(true_nevents) == len(true_model_probs))

    mean_true_model_prob = sum(true_model_probs) / len(true_model_probs)
    median_true_model_prob = pycoevolity.stats.median(true_model_probs)

    nevents_within_95_cred = 0
    model_within_95_cred = 0
    ncorrect = 0
    model_ncorrect = 0
    for i in range(len(true_nevents)):
        true_map_nevents[map_nevents[i] - 1][true_nevents[i] - 1] += 1
        true_map_nevents_probs[map_nevents[i] - 1][true_nevents[i] - 1].append(map_nevents_probs[i])
        if show_all_models:
            true_map_model[model_to_index[map_model[i]]][model_to_index[true_model[i]]] += 1
        if true_nevents_cred_levels[i] <= 0.95:
            nevents_within_95_cred += 1
        if true_model_cred_levels[i] <= 0.95:
            model_within_95_cred += 1
        if true_nevents[i] == map_nevents[i]:
            ncorrect += 1
        if true_model[i] == map_model[i]:
            model_ncorrect += 1
    p_nevents_within_95_cred = nevents_within_95_cred / float(len(true_nevents))
    p_model_within_95_cred = model_within_95_cred / float(len(true_nevents))
    p_correct = ncorrect / float(len(true_nevents))
    p_model_correct = model_ncorrect /  float(len(true_nevents))

    _LOG.info("p(nevents within CS) = {0:.4f}".format(p_nevents_within_95_cred))
    _LOG.info("p(model within CS) = {0:.4f}".format(p_model_within_95_cred))
    ax = plt.subplot(gs[0, 0])

    if show_all_models:
        ax.imshow(true_map_model,
                origin = 'lower',
                cmap = cmap,
                interpolation = 'none',
                aspect = 'auto'
                )
        for i, row_list in enumerate(true_map_model):
            for j, n in enumerate(row_list):
                ax.text(j, i,
                        str(n),
                        horizontalalignment = "center",
                        verticalalignment = "center",
                        # fontsize = 8,
                        )
    else:
        ax.imshow(true_map_nevents,
                origin = 'lower',
                cmap = cmap,
                interpolation = 'none',
                aspect = 'auto'
                )
        for i, row_list in enumerate(true_map_nevents):
            for j, num_events in enumerate(row_list):
                ax.text(j, i,
                        str(num_events),
                        horizontalalignment = "center",
                        verticalalignment = "center")

    # Currently, we are exclusively showing the model (rather than the number
    # of events) plots. For consistency, let's annotate all plots with summary
    # stats about the models (and not the nevents)
    if include_cs:
        # if show_all_models:
        ax.text(0.98, lower_annotation_y,
                "$p(\\mathcal{{T}} \\in \\textrm{{\\sffamily CS}}) = {0:.3f}$".format(
                        p_model_within_95_cred),
                horizontalalignment = "right",
                verticalalignment = "bottom",
                transform = ax.transAxes)
        # else:
        #     ax.text(0.98, lower_annotation_y,
        #             "$p(k \\in \\textrm{{\\sffamily CS}}) = {0:.3f}$".format(
        #                     p_nevents_within_95_cred),
        #             horizontalalignment = "right",
        #             verticalalignment = "bottom",
        #             transform = ax.transAxes)
    if include_prop_correct:
        # if show_all_models:
        ax.text(0.02, upper_annotation_y,
                "$p(\\hat{{\\mathcal{{T}}}} = \\mathcal{{T}}) = {0:.3f}$".format(
                        p_model_correct),
                horizontalalignment = "left",
                verticalalignment = "bottom",
                transform = ax.transAxes)
        # else:
        #     ax.text(0.02, upper_annotation_y,
        #             "$p(\\hat{{k}} = k) = {0:.3f}$".format(
        #                     p_correct),
        #             horizontalalignment = "left",
        #             verticalalignment = "bottom",
        #             transform = ax.transAxes)
    if include_median:
        # if show_all_models:
        ax.text(0.98, upper_annotation_y,
                "$\\widetilde{{p(\\mathcal{{T}}|\\mathbf{{D}})}} = {0:.3f}$".format(
                        median_true_model_prob),
                horizontalalignment = "right",
                verticalalignment = "bottom",
                transform = ax.transAxes)
        # else:
        #     ax.text(0.98, upper_annotation_y,
        #             "$\\widetilde{{p(k|\\mathbf{{D}})}} = {0:.3f}$".format(
        #                     median_true_nevents_prob),
        #             horizontalalignment = "right",
        #             verticalalignment = "bottom",
        #             transform = ax.transAxes)
    if include_x_label:
        if show_all_models:
            ax.set_xlabel("True model ($\\mathcal{{T}}$)",
                    # labelpad = 8.0,
                    fontsize = xy_label_size)
        else:
            ax.set_xlabel("True \\# of events ($k$)",
                    # labelpad = 8.0,
                    fontsize = xy_label_size)
    if include_y_label:
        if show_all_models:
            ax.set_ylabel("MAP model ($\\hat{{\\mathcal{{T}}}}$)",
                    labelpad = 8.0,
                    fontsize = xy_label_size)
        else:
            ax.set_ylabel("MAP \\# of events ($\\hat{{k}}$)",
                    labelpad = 8.0,
                    fontsize = xy_label_size)
    if plot_title:
        ax.set_title(plot_title,
                fontsize = title_size)

    # Make sure ticks correspond only with number of events
    # ax.xaxis.set_ticks(range(number_of_comparisons))
    # ax.yaxis.set_ticks(range(number_of_comparisons))
    # xtick_labels = [item for item in ax.get_xticklabels()]
    # for i in range(len(xtick_labels)):
    #     xtick_labels[i].set_text(str(i + 1))
    # ytick_labels = [item for item in ax.get_yticklabels()]
    # for i in range(len(ytick_labels)):
    #     ytick_labels[i].set_text(str(i + 1))
    # ax.set_xticklabels(xtick_labels)
    # ax.set_yticklabels(ytick_labels)

    # Make sure ticks correspond only with number of events or model
    if not show_all_models:
        ax.xaxis.set_ticks(range(number_of_comparisons))
        ax.yaxis.set_ticks(range(number_of_comparisons))
    else:
        ax.xaxis.set_ticks(range(5))
        ax.yaxis.set_ticks(range(5))
    xtick_labels = [item for item in ax.get_xticklabels()]
    for i in range(len(xtick_labels)):
        if show_all_models:
            xtick_labels[i].set_text(index_to_model[i])
        else:
            xtick_labels[i].set_text(str(i + 1))
    ytick_labels = [item for item in ax.get_yticklabels()]
    for i in range(len(ytick_labels)):
        if show_all_models:
            ytick_labels[i].set_text(index_to_model[i])
        else:
            ytick_labels[i].set_text(str(i + 1))
    ax.set_xticklabels(xtick_labels)
    ax.set_yticklabels(ytick_labels)



    gs.update(
            left = pad_left,
            right = pad_right,
            bottom = pad_bottom,
            top = pad_top)

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    plot_path = os.path.join(plot_dir,
            "{0}-nevents.pdf".format(plot_file_prefix))
    plt.savefig(plot_path)
    _LOG.info("Plots written to {0!r}\n".format(plot_path))


def parse_results(paths):
    return pycoevolity.parsing.get_dict_from_spreadsheets(
            paths,
            sep = "\t",
            offset = 0)


def main_cli(argv = sys.argv):
    if not os.path.exists(project_util.RESULTS_DIR):
        os.mkdir(project_util.RESULTS_DIR)
    if not os.path.exists(project_util.PLOT_DIR):
        os.mkdir(project_util.PLOT_DIR)

    brooks_gelman_1998_recommended_psrf = 1.2

    pad_left = 0.16
    pad_right = 0.94
    pad_bottom = 0.12
    pad_top = 0.965
    plot_width = 2.8
    plot_height = 2.2

    # Plot relative root priors
    root_gamma_parameters = (
            (5.0, 0.01, 0.05, 8.0, False),
            (5.0, 0.04, 0.05, 8.0, False),
            (5.0, 0.09, 0.05, 8.0, False),
            (5.0, 0.19, 0.05, 8.0, False),
            (5.0, 0.79, 0.05, 8.0, False),
            (50.0, 0.02, 0.0, 8.0, False),
            (500.0, 0.002, 0.0, 8.0, False),
            (10.0, 0.025, 0.0, 4.0, False),
            (10.0, 0.05, 0.0, 4.0, False),
            (10.0, 0.1, 0.0, 4.0, False),
            (10.0, 0.2, 0.0, 4.0, False),
            (100.0, 0.01, 0.0, 4.0, False),
            (1000.0, 0.001, 0.0, 4.0, False),
            (5.0, 0.04, 0.05, 4.5, True),
            (5.0, 0.04, 3.8, 4.5, True),
            (1.0, 2.0, 0.0, 4.5, True),
            (5.0, 0.01, 0.05, 4.5, False),
            (5.0, 0.04, 0.05, 10.5, True),
            (5.0, 0.04, 3.8, 10.5, True),
            (5.0, 0.01, 0.05, 10.5, True),
            (5.0, 0.01, 9.95, 10.5, True),
            )
    # x_max = float("-inf")
    # for shape, scale, offset in root_gamma_parameters:
    #     q = scipy.stats.gamma.ppf(0.975, shape, scale = scale)
    #     q += offset
    #     if q > x_max:
    #         x_max = q
    for i, (shape, scale, offset, x_max, include_x_max_in_name) in enumerate(root_gamma_parameters):
        shape_str = str(shape).replace(".", "_")
        scale_str = str(scale).replace(".", "_")
        offset_str = str(offset).replace(".", "_")
        x_max_str = str(x_max).replace(".", "_")
        if include_x_max_in_name:
            plot_file_prefix = "prior-relative-root-{0}-{1}-{2}-{3}".format(
                    shape_str, scale_str, offset_str, x_max_str)
        else:
            plot_file_prefix = "prior-relative-root-{0}-{1}-{2}".format(
                    shape_str, scale_str, offset_str)
        include_y_label = False
        if i == 0:
            include_y_label = True
        x_min = 0.0
        if x_max == 10.5:
            x_min = -0.5
        plot_gamma(shape = shape,
                scale = scale,
                offset = offset,
                x_min = x_min,
                x_max = x_max,
                number_of_points = 10000,
                x_label = "Relative root size",
                include_x_label = False, # Going to add x-axis in slides
                include_y_label = False,
                include_title = False,
                curve_line_width = 2.5,
                curve_line_style = '-',
                curve_line_color = '0.45',
                one_line_width = 1.5,
                one_line_style = '--',
                one_line_color = '0.7',
                plot_width = plot_width,
                plot_height = plot_height,
                xy_label_size = 16.0,
                title_size = 16.0,
                pad_left = pad_left,
                pad_right = pad_right,
                pad_bottom = pad_bottom,
                pad_top = pad_top,
                plot_file_prefix = plot_file_prefix)

    time_gamma_parameters = (
            (4.0, 0.000475, 0.0001),
            (1.0, 0.01, 0.0),
            (1.0, 0.001, 0.0),
            )
    x_max = float("-inf")
    for shape, scale, offset in time_gamma_parameters:
        q = scipy.stats.gamma.ppf(0.975, shape, scale = scale)
        if q > x_max:
            x_max = q
    for i, (shape, scale, offset) in enumerate(time_gamma_parameters):
        shape_str = str(shape).replace(".", "_")
        scale_str = str(scale).replace(".", "_")
        offset_str = str(offset).replace(".", "_")
        plot_file_prefix = "prior-time-{0}-{1}-{2}".format(shape_str, scale_str, offset_str)
        plot_gamma(shape = shape,
                scale = scale,
                offset = offset,
                x_min = 0.0,
                x_max = x_max,
                number_of_points = 10000,
                x_label = "Time",
                include_x_label = False, # Going to add x-axis in slides
                include_y_label = False,
                include_title = False,
                curve_line_width = 2.5,
                curve_line_style = '-',
                curve_line_color = '0.45',
                one_line_width = 1.5,
                one_line_style = '--',
                one_line_color = '0.7',
                plot_width = plot_width,
                plot_height = plot_height,
                xy_label_size = 16.0,
                title_size = 16.0,
                pad_left = pad_left,
                pad_right = pad_right,
                pad_bottom = pad_bottom,
                pad_top = pad_top,
                plot_file_prefix = plot_file_prefix)


    sim_dirs_opt = [
            "03pops-dpp-root-0005-004-t0002-500k",
            "03pops-dpp-root-0005-009-t0002-500k",
            "03pops-dpp-root-0005-019-t0002-500k",
            "03pops-dpp-root-0005-079-t0002-500k",
    ]
    sim_dirs_opt_centered = [
            "03pops-dpp-root-0005-019-t0002-500k",
            "03pops-dpp-root-0050-002-t0002-500k",
            "03pops-dpp-root-0500-0002-t0002-500k",
    ]
    sim_dirs_opt_all = [
            "03pops-dpp-root-0005-004-t0002-500k",
            "03pops-dpp-root-0005-009-t0002-500k",
            "03pops-dpp-root-0005-079-t0002-500k",
            "03pops-dpp-root-0005-019-t0002-500k",
            "03pops-dpp-root-0050-002-t0002-500k",
            "03pops-dpp-root-0500-0002-t0002-500k",
            "03pops-dpp-root-0005-004-t0002-500k-diffuseprior",
            "03pops-dpp-root-0005-004-3_8-t0002-500k-diffuseprior",
            "03pops-dpp-root-0005-004-t0002-500k-40genomes",
            "03pops-dpp-root-0005-001-t0002-500k-diffuseprior",
            "03pops-dpp-root-0005-001-9_95-t0002-500k-diffuseprior",
    ]

    sim_dirs_opt_6_pops = [
            "06pops-dpp-root-0005-004-t0002-500k",
    ]

    sim_dirs_mixed_comparisons = [
            "03pops-03pairs-dpp-root-0005-009-t0002-500k",
            ]

    sim_dirs_pairs = [
            "03pairs-dpp-root-0005-004-t0002-500k-diffuseprior",
            "03pairs-dpp-root-0005-004-t0002-500k",
            ]

    sim_dirs = [
            "03pops-dpp-root-0005-004-t0002-500k-0100l",
            "03pops-dpp-root-0005-004-t0002-500k-0100ul",
            "03pops-dpp-root-0005-004-t0002-500k-040s",
            "03pops-dpp-root-0005-004-t0002-500k-060s",
            "03pops-dpp-root-0010-0025-500k",
            "03pops-dpp-root-0010-0025-t0001-500k",
            "03pops-dpp-root-0010-005-500k",
            "03pops-dpp-root-0010-005-t0001-500k",
            "03pops-dpp-root-0010-010-500k",
            "03pops-dpp-root-0010-020-500k",
            "03pops-dpp-root-0100-001-500k",
            "03pops-dpp-root-1000-0001-500k",
            ] + sim_dirs_opt_all

    all_sim_dirs = sim_dirs + sim_dirs_mixed_comparisons


    results = {}
    var_only_results = {}
    for sim_dir in all_sim_dirs:
        results[sim_dir] = parse_results(glob.glob(os.path.join(project_util.VAL_DIR,
                        sim_dir,
                        "batch00?",
                        "results.csv.gz")))
        var_only_results_paths = glob.glob(os.path.join(project_util.VAL_DIR,
                        sim_dir,
                        "batch00?",
                        "var-only-results.csv.gz"))
        if var_only_results_paths:
            var_only_results[sim_dir] = parse_results(var_only_results_paths)

    results_mixed_comps = {}
    var_only_results_mixed_comps = {}
    for sim_dir in sim_dirs_mixed_comparisons:
        results_mixed_comps[sim_dir] = results[sim_dir]
        if sim_dir in var_only_results:
            var_only_results_mixed_comps[sim_dir] = var_only_results[sim_dir]

    results_pairs = {}
    var_only_results_pairs = {}
    for sim_dir in sim_dirs_pairs:
        results_pairs[sim_dir] = parse_results(glob.glob(os.path.join(project_util.VAL_DIR,
                        sim_dir,
                        "batch00?",
                        "results.csv.gz")))
        var_only_results_pairs_paths = glob.glob(os.path.join(project_util.VAL_DIR,
                        sim_dir,
                        "batch00?",
                        "var-only-results.csv.gz"))
        if var_only_results_pairs_paths:
            var_only_results_pairs[sim_dir] = parse_results(var_only_results_pairs_paths)

    results_6_pops = {}
    var_only_results_6_pops = {}
    for sim_dir in sim_dirs_opt_6_pops:
        results_6_pops[sim_dir] = parse_results(glob.glob(os.path.join(project_util.VAL_DIR,
                sim_dir,
                "batch00?",
                "results.csv.gz")))
        var_only_results_pairs_paths = glob.glob(os.path.join(project_util.VAL_DIR,
                        sim_dir,
                        "batch00?",
                        "var-only-results.csv.gz"))
        if var_only_results_pairs_paths:
            var_only_results_6_pops[sim_dir] = parse_results(var_only_results_pairs_paths)


    height_parameters = [
            "root_height_c1sp1",
            "root_height_c2sp1",
            "root_height_c3sp1",
    ]
    div_height_parameters = [
            "root_height_c4sp1",
            "root_height_c5sp1",
            "root_height_c6sp1",
    ]
    coal_height_parameters = [
            "coal_root_height_c1sp1",
            "coal_root_height_c2sp1",
            "coal_root_height_c3sp1",
    ]
    div_coal_height_parameters = [
            "coal_root_height_c4sp1",
            "coal_root_height_c5sp1",
            "coal_root_height_c6sp1",
    ]
    root_size_parameters = [
            "pop_size_root_c1sp1",
            "pop_size_root_c2sp1",
            "pop_size_root_c3sp1",
    ]
    div_root_size_parameters = [
            "pop_size_root_c4sp1",
            "pop_size_root_c5sp1",
            "pop_size_root_c6sp1",
    ]
    leaf_size_parameters = [
            "pop_size_c1sp1",
            "pop_size_c2sp1",
            "pop_size_c3sp1",
    ]
    div_leaf_size_parameters = [
            "pop_size_c4sp1",
            "pop_size_c5sp1",
            "pop_size_c6sp1",
            # "pop_size_c4sp2",
            # "pop_size_c5sp2",
            # "pop_size_c6sp2",
    ]


    parameters_to_plot = {
            "event-time": {
                    "headers": height_parameters,
                    "label": "event time",
                    "short_label": "time",
                    "symbol": "t",
                    # "xy_limits": (0.0, 0.008, 0.0, 0.008),
                    "xy_limits": "calc_by_time_prior",
                    "pad_left": pad_left,
                    "mixed_comp_only": False,
                    "exclude_pairs_only": True,
                    "six_pops_extra_headers": div_height_parameters,
            },
            "div-time": {
                    "headers": div_height_parameters,
                    "label": "divergence time",
                    "short_label": "div time",
                    "symbol": "t",
                    # "xy_limits": (0.0, 0.008, 0.0, 0.008),
                    "xy_limits": "calc_by_time_prior",
                    "pad_left": pad_left,
                    "mixed_comp_only": True,
                    "exclude_pairs_only": False,
            },
            "div-demog-time": {
                    "headers": height_parameters + div_height_parameters,
                    "label": "event time",
                    "short_label": "time",
                    "symbol": "t",
                    # "xy_limits": (0.0, 0.008, 0.0, 0.008),
                    "xy_limits": "calc_by_time_prior",
                    "pad_left": pad_left,
                    "mixed_comp_only": True,
                    "exclude_pairs_only": True,
            },
            "event-coal-time": {
                    "headers": coal_height_parameters,
                    "label": "event time in coalescent units",
                    "short_label": "time in coal units",
                    "symbol": "t",
                    "xy_limits": None,
                    "pad_left": pad_left - 0.02,
                    "mixed_comp_only": False,
                    "exclude_pairs_only": True,
                    "six_pops_extra_headers": div_coal_height_parameters,
            },
            "event-div-coal-time": {
                    "headers": div_coal_height_parameters,
                    "label": "divergence time in coalescent units",
                    "short_label": "div time in coal units",
                    "symbol": "t",
                    "xy_limits": None,
                    "pad_left": pad_left - 0.02,
                    "mixed_comp_only": True,
                    "exclude_pairs_only": False,
            },
            "event-div-demog-coal-time": {
                    "headers": coal_height_parameters + div_coal_height_parameters,
                    "label": "event time in coalescent units",
                    "short_label": "time in coal units",
                    "symbol": "t",
                    "xy_limits": None,
                    "pad_left": pad_left - 0.02,
                    "mixed_comp_only": True,
                    "exclude_pairs_only": True,
            },
            "ancestor-size": {
                    "headers": root_size_parameters,
                    "label": "ancestor population size",
                    "short_label": "size",
                    "symbol": "N_e\\mu",
                    "xy_limits": None,
                    "pad_left": pad_left + 0.01,
                    "mixed_comp_only": False,
                    "exclude_pairs_only": True,
                    "six_pops_extra_headers": div_root_size_parameters,
            },
            "div-ancestor-size": {
                    "headers": div_root_size_parameters,
                    "label": "ancestor population size",
                    "short_label": "size",
                    "symbol": "N_e\\mu",
                    "xy_limits": None,
                    "pad_left": pad_left + 0.01,
                    "mixed_comp_only": True,
                    "exclude_pairs_only": False,
            },
            "div-demog-ancestor-size": {
                    "headers": root_size_parameters + div_root_size_parameters,
                    "label": "ancestor population size",
                    "short_label": "size",
                    "symbol": "N_e\\mu",
                    "xy_limits": None,
                    "pad_left": pad_left + 0.01,
                    "mixed_comp_only": True,
                    "exclude_pairs_only": True,
            },
            "descendant-size": {
                    "headers": leaf_size_parameters,
                    "label": "descendant population size",
                    "short_label": "size",
                    "symbol": "N_e\\mu",
                    # "xy_limits": (0.0, 0.008, 0.0, 0.008),
                    "xy_limits": "calc_by_time_prior",
                    "pad_left": pad_left,
                    "mixed_comp_only": False,
                    "exclude_pairs_only": True,
                    "six_pops_extra_headers": div_leaf_size_parameters,
            },
            "div-descendant-size": {
                    "headers": div_leaf_size_parameters,
                    "label": "descendant population size",
                    "short_label": "size",
                    "symbol": "N_e\\mu",
                    # "xy_limits": (0.0, 0.008, 0.0, 0.008),
                    "xy_limits": "calc_by_time_prior",
                    "pad_left": pad_left,
                    "mixed_comp_only": True,
                    "exclude_pairs_only": False,
            },
            "div-demog-descendant-size": {
                    "headers": leaf_size_parameters + div_leaf_size_parameters,
                    "label": "descendant population size",
                    "short_label": "size",
                    "symbol": "N_e\\mu",
                    # "xy_limits": (0.0, 0.008, 0.0, 0.008),
                    "xy_limits": "calc_by_time_prior",
                    "pad_left": pad_left,
                    "mixed_comp_only": True,
                    "exclude_pairs_only": True,
            },
    }


    for parameter, p_info in parameters_to_plot.items():
        data = {}
        var_only_data = {}
        if p_info["mixed_comp_only"]:
            for sim_dir, r in results_mixed_comps.items():
                data[sim_dir] = ScatterData.init(r, p_info["headers"],
                        highlight_parameter_prefix = "psrf",
                        highlight_threshold = brooks_gelman_1998_recommended_psrf,
                        )
            for sim_dir, r in var_only_results_mixed_comps.items():
                var_only_data[sim_dir] = ScatterData.init(r, p_info["headers"],
                        highlight_parameter_prefix = "psrf",
                        highlight_threshold = brooks_gelman_1998_recommended_psrf,
                        )
        else:
            for sim_dir, r in results.items():
                data[sim_dir] = ScatterData.init(r, p_info["headers"],
                        highlight_parameter_prefix = "psrf",
                        highlight_threshold = brooks_gelman_1998_recommended_psrf,
                        )
            for sim_dir, r in var_only_results.items():
                var_only_data[sim_dir] = ScatterData.init(r, p_info["headers"],
                        highlight_parameter_prefix = "psrf",
                        highlight_threshold = brooks_gelman_1998_recommended_psrf,
                        )
            for sim_dir, r in results_6_pops.items():
                data[sim_dir] = ScatterData.init(r,
                        p_info["headers"] + p_info["six_pops_extra_headers"],
                        highlight_parameter_prefix = "psrf",
                        highlight_threshold = brooks_gelman_1998_recommended_psrf,
                        )
            for sim_dir, r in var_only_results_6_pops.items():
                var_only_data[sim_dir] = ScatterData.init(r,
                        p_info["headers"] + p_info["six_pops_extra_headers"],
                        highlight_parameter_prefix = "psrf",
                        highlight_threshold = brooks_gelman_1998_recommended_psrf,
                        )
        if not p_info["exclude_pairs_only"]:
            for sim_dir, r in results_pairs.items():
                data[sim_dir] = ScatterData.init(r, p_info["headers"],
                        highlight_parameter_prefix = "psrf",
                        highlight_threshold = brooks_gelman_1998_recommended_psrf,
                        )
            for sim_dir, r in var_only_results_pairs.items():
                var_only_data[sim_dir] = ScatterData.init(r, p_info["headers"],
                        highlight_parameter_prefix = "psrf",
                        highlight_threshold = brooks_gelman_1998_recommended_psrf,
                        )
        
        x_label = "True {0} (${1}$)".format(
                p_info["label"],
                p_info["symbol"])
        y_label = "Estimated {0} ($\\hat{{{1}}}$)".format(
                p_info["label"],
                p_info["symbol"])

        maximums = {}
        if p_info["xy_limits"] == "calc_by_time_prior":
            maximums = dict(zip(unique_time_priors, (float('-inf') for i in range(len(unique_time_priors)))))
            for sim_dir in data.keys():
                if sim_dir.endswith("-040s") or sim_dir.endswith("-060s"):
                    continue
                time_prior = sim_dir_to_priors[sim_dir][1]
                if sim_dir in var_only_data:
                    maximums[time_prior] = max(maximums[time_prior],
                            max(data[sim_dir].x),
                            max(data[sim_dir].y),
                            max(var_only_data[sim_dir].x),
                            max(var_only_data[sim_dir].y),
                            )
                else:
                    maximums[time_prior] = max(maximums[time_prior],
                            max(data[sim_dir].x),
                            max(data[sim_dir].y),
                            )

        # Generate individual scatters
        for sim_dir in data.keys():
            xy_lim = p_info["xy_limits"]
            if sim_dir.endswith("-040s") or sim_dir.endswith("-060s"):
                xy_lim = None
            if xy_lim == "calc_by_time_prior":
                time_prior = sim_dir_to_priors[sim_dir][1]
                mx = maximums[time_prior]
                xy_lim = (0.0, mx * 1.02, 0.0, mx * 1.02)

            prefix = get_prefix_from_sim_dir_name(sim_dir)
            y_label = "Estimated {0} ($\\hat{{{1}}}$)".format(
                    p_info["short_label"],
                    p_info["symbol"])
            generate_specific_scatter_plot(
                    data = data[sim_dir],
                    plot_file_prefix = parameter + "-" + prefix,
                    parameter_symbol = p_info["symbol"],
                    title = None,
                    title_size = 16.0,
                    x_label = None,
                    x_label_size = 16.0,
                    y_label = None,
                    y_label_size = 16.0,
                    plot_width = plot_width,
                    plot_height = plot_height,
                    pad_left = p_info["pad_left"],
                    pad_right = pad_right,
                    pad_bottom = pad_bottom,
                    pad_top = pad_top,
                    force_shared_xy_ranges = True,
                    xy_limits = xy_lim,
                    include_coverage = True,
                    include_rmse = True,
                    include_identity_line = True,
                    include_error_bars = True,
                    )
            if sim_dir in var_only_data:
                generate_specific_scatter_plot(
                        data = var_only_data[sim_dir],
                        plot_file_prefix = "var-only-" + parameter + "-" + prefix,
                        parameter_symbol = p_info["symbol"],
                        title = None,
                        title_size = 16.0,
                        x_label = None,
                        x_label_size = 16.0,
                        y_label = None,
                        y_label_size = 16.0,
                        plot_width = plot_width,
                        plot_height = plot_height,
                        pad_left = p_info["pad_left"],
                        pad_right = pad_right,
                        pad_bottom = pad_bottom,
                        pad_top = pad_top,
                        force_shared_xy_ranges = True,
                        xy_limits = xy_lim,
                        include_coverage = True,
                        include_rmse = True,
                        include_identity_line = True,
                        include_error_bars = True,
                        )



    # Generate individual model plots
    for sim_dir in sim_dirs + sim_dirs_pairs:
        prefix = get_prefix_from_sim_dir_name(sim_dir)
        r = results.get(sim_dir, results_pairs.get(sim_dir))
        generate_specific_model_plots(
                results = r,
                number_of_comparisons = 3,
                plot_title = None,
                include_x_label = False,
                include_y_label = False,
                include_median = True,
                include_cs = True,
                include_prop_correct = True,
                plot_width = plot_width * 0.94,
                plot_height = plot_height,
                xy_label_size = 16.0,
                title_size = 16.0,
                pad_left = pad_left - 0.04,
                pad_right = 0.985,
                pad_bottom = pad_bottom - 0.015,
                pad_top = pad_top - 0.07,
                lower_annotation_y = 0.01,
                upper_annotation_y = 1.015,
                plot_file_prefix = "nevents-" + prefix)
        generate_specific_model_plots(
                results = r,
                number_of_comparisons = 3,
                show_all_models = True,
                plot_title = None,
                include_x_label = False,
                include_y_label = False,
                include_median = True,
                include_cs = False,
                include_prop_correct = True,
                plot_width = plot_width * 0.96,
                plot_height = plot_height,
                xy_label_size = 16.0,
                title_size = 16.0,
                pad_left = pad_left - 0.03,
                pad_right = 0.985,
                pad_bottom = pad_bottom - 0.015,
                pad_top = pad_top - 0.07,
                lower_annotation_y = 0.01,
                upper_annotation_y = 1.015,
                plot_file_prefix = "model-" + prefix)
        if (sim_dir in var_only_results) or (sim_dir in var_only_results_pairs):
            vor = var_only_results.get(sim_dir, var_only_results_pairs.get(sim_dir))
            generate_specific_model_plots(
                    results = vor,
                    number_of_comparisons = 3,
                    plot_title = None,
                    include_x_label = False,
                    include_y_label = False,
                    include_median = True,
                    include_cs = True,
                    include_prop_correct = True,
                    plot_width = plot_width * 0.94,
                    plot_height = plot_height,
                    xy_label_size = 16.0,
                    title_size = 16.0,
                    pad_left = pad_left - 0.04,
                    pad_right = 0.985,
                    pad_bottom = pad_bottom - 0.015,
                    pad_top = pad_top - 0.07,
                    lower_annotation_y = 0.01,
                    upper_annotation_y = 1.015,
                    plot_file_prefix = "var-only-nevents-" + prefix)
            generate_specific_model_plots(
                    results = vor,
                    number_of_comparisons = 3,
                    show_all_models = True,
                    plot_title = None,
                    include_x_label = False,
                    include_y_label = False,
                    include_median = True,
                    include_cs = False,
                    include_prop_correct = True,
                    plot_width = plot_width * 0.96,
                    plot_height = plot_height,
                    xy_label_size = 16.0,
                    title_size = 16.0,
                    pad_left = pad_left - 0.03,
                    pad_right = 0.985,
                    pad_bottom = pad_bottom - 0.015,
                    pad_top = pad_top - 0.07,
                    lower_annotation_y = 0.01,
                    upper_annotation_y = 1.015,
                    plot_file_prefix = "var-only-model-" + prefix)

    for sim_dir in sim_dirs_mixed_comparisons:
        prefix = get_prefix_from_sim_dir_name(sim_dir)
        generate_specific_model_plots(
                results = results[sim_dir],
                number_of_comparisons = 6,
                plot_title = None,
                include_x_label = False,
                include_y_label = False,
                include_median = True,
                include_cs = True,
                include_prop_correct = True,
                plot_width = plot_width * 0.94,
                plot_height = plot_height,
                xy_label_size = 16.0,
                title_size = 16.0,
                pad_left = pad_left - 0.04,
                pad_right = 0.985,
                pad_bottom = pad_bottom - 0.015,
                pad_top = pad_top - 0.07,
                lower_annotation_y = 0.01,
                upper_annotation_y = 1.015,
                plot_file_prefix = "nevents-" + prefix)
        generate_specific_model_plots(
                results = results[sim_dir],
                number_of_comparisons = 3,
                plot_title = None,
                include_x_label = False,
                include_y_label = False,
                include_median = True,
                include_cs = True,
                include_prop_correct = True,
                plot_width = plot_width * 0.94,
                plot_height = plot_height,
                xy_label_size = 16.0,
                title_size = 16.0,
                pad_left = pad_left - 0.04,
                pad_right = 0.985,
                pad_bottom = pad_bottom - 0.015,
                pad_top = pad_top - 0.07,
                lower_annotation_y = 0.01,
                upper_annotation_y = 1.015,
                plot_file_prefix = "num-div-events-" + prefix,
                model_key = "div_model",
                num_events_key = "num_div_events",
                )
        generate_specific_model_plots(
                results = results[sim_dir],
                number_of_comparisons = 3,
                plot_title = None,
                include_x_label = False,
                include_y_label = False,
                include_median = True,
                include_cs = True,
                include_prop_correct = True,
                plot_width = plot_width * 0.94,
                plot_height = plot_height,
                xy_label_size = 16.0,
                title_size = 16.0,
                pad_left = pad_left - 0.04,
                pad_right = 0.985,
                pad_bottom = pad_bottom - 0.015,
                pad_top = pad_top - 0.07,
                lower_annotation_y = 0.01,
                upper_annotation_y = 1.015,
                plot_file_prefix = "num-demog-events-" + prefix,
                model_key = "demog_model",
                num_events_key = "num_demog_events",
                )
        generate_specific_model_plots(
                results = results[sim_dir],
                number_of_comparisons = 3,
                show_all_models = True,
                plot_title = None,
                include_x_label = False,
                include_y_label = False,
                include_median = True,
                include_cs = False,
                include_prop_correct = True,
                plot_width = plot_width * 0.94,
                plot_height = plot_height,
                xy_label_size = 16.0,
                title_size = 16.0,
                pad_left = pad_left - 0.04,
                pad_right = 0.985,
                pad_bottom = pad_bottom - 0.015,
                pad_top = pad_top - 0.07,
                lower_annotation_y = 0.01,
                upper_annotation_y = 1.015,
                plot_file_prefix = "div-model-" + prefix,
                model_key = "div_model",
                num_events_key = "num_div_events",
                )
        generate_specific_model_plots(
                results = results[sim_dir],
                number_of_comparisons = 3,
                show_all_models = True,
                plot_title = None,
                include_x_label = False,
                include_y_label = False,
                include_median = True,
                include_cs = False,
                include_prop_correct = True,
                plot_width = plot_width * 0.94,
                plot_height = plot_height,
                xy_label_size = 16.0,
                title_size = 16.0,
                pad_left = pad_left - 0.04,
                pad_right = 0.985,
                pad_bottom = pad_bottom - 0.015,
                pad_top = pad_top - 0.07,
                lower_annotation_y = 0.01,
                upper_annotation_y = 1.015,
                plot_file_prefix = "demog-model-" + prefix,
                model_key = "demog_model",
                num_events_key = "num_demog_events",
                )
        if sim_dir in var_only_results:
            generate_specific_model_plots(
                    results = var_only_results[sim_dir],
                    number_of_comparisons = 6,
                    plot_title = None,
                    include_x_label = False,
                    include_y_label = False,
                    include_median = True,
                    include_cs = True,
                    include_prop_correct = True,
                    plot_width = plot_width * 0.94,
                    plot_height = plot_height,
                    xy_label_size = 16.0,
                    title_size = 16.0,
                    pad_left = pad_left - 0.04,
                    pad_right = 0.985,
                    pad_bottom = pad_bottom - 0.015,
                    pad_top = pad_top - 0.07,
                    lower_annotation_y = 0.01,
                    upper_annotation_y = 1.015,
                    plot_file_prefix = "var-only-nevents-" + prefix)
            generate_specific_model_plots(
                    results = var_only_results[sim_dir],
                    number_of_comparisons = 3,
                    plot_title = None,
                    include_x_label = False,
                    include_y_label = False,
                    include_median = True,
                    include_cs = True,
                    include_prop_correct = True,
                    plot_width = plot_width * 0.94,
                    plot_height = plot_height,
                    xy_label_size = 16.0,
                    title_size = 16.0,
                    pad_left = pad_left - 0.04,
                    pad_right = 0.985,
                    pad_bottom = pad_bottom - 0.015,
                    pad_top = pad_top - 0.07,
                    lower_annotation_y = 0.01,
                    upper_annotation_y = 1.015,
                    plot_file_prefix = "var-only-num-div-events-" + prefix,
                    model_key = "div_model",
                    num_events_key = "num_div_events",
                    )
            generate_specific_model_plots(
                    results = var_only_results[sim_dir],
                    number_of_comparisons = 3,
                    plot_title = None,
                    include_x_label = False,
                    include_y_label = False,
                    include_median = True,
                    include_cs = True,
                    include_prop_correct = True,
                    plot_width = plot_width * 0.94,
                    plot_height = plot_height,
                    xy_label_size = 16.0,
                    title_size = 16.0,
                    pad_left = pad_left - 0.04,
                    pad_right = 0.985,
                    pad_bottom = pad_bottom - 0.015,
                    pad_top = pad_top - 0.07,
                    lower_annotation_y = 0.01,
                    upper_annotation_y = 1.015,
                    plot_file_prefix = "var-only-num-demog-events-" + prefix,
                    model_key = "demog_model",
                    num_events_key = "num_demog_events",
                    )
            generate_specific_model_plots(
                    results = var_only_results[sim_dir],
                    number_of_comparisons = 3,
                    show_all_models = True,
                    plot_title = None,
                    include_x_label = False,
                    include_y_label = False,
                    include_median = True,
                    include_cs = False,
                    include_prop_correct = True,
                    plot_width = plot_width * 0.94,
                    plot_height = plot_height,
                    xy_label_size = 16.0,
                    title_size = 16.0,
                    pad_left = pad_left - 0.04,
                    pad_right = 0.985,
                    pad_bottom = pad_bottom - 0.015,
                    pad_top = pad_top - 0.07,
                    lower_annotation_y = 0.01,
                    upper_annotation_y = 1.015,
                    plot_file_prefix = "var-only-div-model-" + prefix,
                    model_key = "div_model",
                    num_events_key = "num_div_events",
                    )
            generate_specific_model_plots(
                    results = var_only_results[sim_dir],
                    number_of_comparisons = 3,
                    show_all_models = True,
                    plot_title = None,
                    include_x_label = False,
                    include_y_label = False,
                    include_median = True,
                    include_cs = False,
                    include_prop_correct = True,
                    plot_width = plot_width * 0.94,
                    plot_height = plot_height,
                    xy_label_size = 16.0,
                    title_size = 16.0,
                    pad_left = pad_left - 0.04,
                    pad_right = 0.985,
                    pad_bottom = pad_bottom - 0.015,
                    pad_top = pad_top - 0.07,
                    lower_annotation_y = 0.01,
                    upper_annotation_y = 1.015,
                    plot_file_prefix = "var-only-demog-model-" + prefix,
                    model_key = "demog_model",
                    num_events_key = "num_demog_events",
                    )

    for sim_dir in sim_dirs_opt_6_pops:
        prefix = get_prefix_from_sim_dir_name(sim_dir)
        r = results_6_pops.get(sim_dir)
        generate_specific_model_plots(
                results = r,
                number_of_comparisons = 6,
                plot_title = None,
                include_x_label = False,
                include_y_label = False,
                include_median = True,
                include_cs = False,
                include_prop_correct = True,
                plot_width = plot_width * 0.96,
                plot_height = plot_height,
                xy_label_size = 16.0,
                title_size = 16.0,
                pad_left = pad_left - 0.03,
                pad_right = 0.985,
                pad_bottom = pad_bottom - 0.015,
                pad_top = pad_top - 0.07,
                lower_annotation_y = 0.01,
                upper_annotation_y = 1.015,
                plot_file_prefix = "nevents-" + prefix)
        if sim_dir in var_only_results_6_pops:
            vor = var_only_results_6_pops.get(sim_dir)
            generate_specific_model_plots(
                    results = vor,
                    number_of_comparisons = 6,
                    plot_title = None,
                    include_x_label = False,
                    include_y_label = False,
                    include_median = True,
                    include_cs = False,
                    include_prop_correct = True,
                    plot_width = plot_width * 0.96,
                    plot_height = plot_height,
                    xy_label_size = 16.0,
                    title_size = 16.0,
                    pad_left = pad_left - 0.03,
                    pad_right = 0.985,
                    pad_bottom = pad_bottom - 0.015,
                    pad_top = pad_top - 0.07,
                    lower_annotation_y = 0.01,
                    upper_annotation_y = 1.015,
                    plot_file_prefix = "var-only-nevents-" + prefix)

    jitter = 0.12
    alpha = 0.5
    for sim_dir in sim_dirs_opt_6_pops:
        prefix = get_prefix_from_sim_dir_name(sim_dir)
        r = results_6_pops.get(sim_dir)
        nshared_v_abs_error, nshared_v_ci_width = BoxData.init_time_v_sharing(
                results = r,
                estimator_prefix = "mean")
        generate_specific_box_plot(
                data = nshared_v_abs_error,
                plot_file_prefix = "shared-v-abs-error-" + prefix,
                title = None,
                title_size = 16.0,
                x_label = None,
                x_label_size = 16.0,
                y_label = None,
                y_label_size = 16.0,
                plot_width = plot_width,
                plot_height = plot_height,
                pad_left = pad_left + 0.02,
                pad_right = pad_right + 0.02,
                pad_bottom = pad_bottom,
                pad_top = pad_top,
                jitter = jitter,
                alpha = alpha,
                rasterized = False)
        generate_specific_box_plot(
                data = nshared_v_ci_width,
                plot_file_prefix = "shared-v-ci-width-" + prefix,
                title = None,
                title_size = 16.0,
                x_label = None,
                x_label_size = 16.0,
                y_label = None,
                y_label_size = 16.0,
                plot_width = plot_width,
                plot_height = plot_height,
                pad_left = pad_left + 0.02,
                pad_right = pad_right + 0.02,
                pad_bottom = pad_bottom,
                pad_top = pad_top,
                jitter = jitter,
                alpha = alpha,
                rasterized = False)
        if sim_dir in var_only_results_6_pops:
            vor = var_only_results_6_pops.get(sim_dir)
            vo_nshared_v_abs_error, vo_nshared_v_ci_width = BoxData.init_time_v_sharing(
                    results = vor,
                    estimator_prefix = "mean")
            generate_specific_box_plot(
                    data = vo_nshared_v_abs_error,
                    plot_file_prefix = "var-only-shared-v-abs-error-" + prefix,
                    title = None,
                    title_size = 16.0,
                    x_label = None,
                    x_label_size = 16.0,
                    y_label = None,
                    y_label_size = 16.0,
                    plot_width = plot_width,
                    plot_height = plot_height,
                    pad_left = pad_left + 0.02,
                    pad_right = pad_right + 0.02,
                    pad_bottom = pad_bottom,
                    pad_top = pad_top,
                    jitter = jitter,
                    alpha = alpha,
                    rasterized = False)
            generate_specific_box_plot(
                    data = vo_nshared_v_ci_width,
                    plot_file_prefix = "var-only-shared-v-ci-width-" + prefix,
                    title = None,
                    title_size = 16.0,
                    x_label = None,
                    x_label_size = 16.0,
                    y_label = None,
                    y_label_size = 16.0,
                    plot_width = plot_width,
                    plot_height = plot_height,
                    pad_left = pad_left + 0.02,
                    pad_right = pad_right + 0.02,
                    pad_bottom = pad_bottom,
                    pad_top = pad_top,
                    jitter = jitter,
                    alpha = alpha,
                    rasterized = False)


    histograms_to_plot = {
            "n-var-sites": {
                    "headers": [
                            "n_var_sites_c1",
                            "n_var_sites_c2",
                            "n_var_sites_c3",
                    ],
                    "label": "Number of variable sites",
                    "short_label": "No. variable sites",
                    "ndigits": 0,
                    "mixed_comp_only": False,
                    "exclude_pairs_only": True,
                    "six_pops_extra_headers": [
                            "n_var_sites_c4",
                            "n_var_sites_c5",
                            "n_var_sites_c6",
                    ],
                    "center_key": "mean",
            },
            "div-n-var-sites": {
                    "headers": [
                            "n_var_sites_c4",
                            "n_var_sites_c5",
                            "n_var_sites_c6",
                    ],
                    "label": "Number of variable sites",
                    "short_label": "No. variable sites",
                    "ndigits": 0,
                    "mixed_comp_only": True,
                    "exclude_pairs_only": False,
                    "six_pops_extra_headers": [],
                    "center_key": "mean",
            },
            "ess-ln-likelihood": {
                    "headers": [
                            "ess_sum_ln_likelihood",
                    ],
                    "label": "Effective sample size of log likelihood",
                    "short_label": "ESS of lnL",
                    "ndigits": 1,
                    "mixed_comp_only": False,
                    "exclude_pairs_only": False,
                    "six_pops_extra_headers": [],
                    "center_key": "mean",
            },
            "ess-event-time": {
                    "headers": [
                            "ess_sum_root_height_c1sp1",
                            "ess_sum_root_height_c2sp1",
                            "ess_sum_root_height_c3sp1",
                    ],
                    "label": "Effective sample size of event time",
                    "short_label": "ESS of time",
                    "ndigits": 1,
                    "mixed_comp_only": False,
                    "exclude_pairs_only": True,
                    "six_pops_extra_headers": [
                            "ess_sum_root_height_c4sp1",
                            "ess_sum_root_height_c5sp1",
                            "ess_sum_root_height_c6sp1",
                    ],
                    "center_key": "mean",
            },
            "ess-div-time": {
                    "headers": [
                            "ess_sum_root_height_c4sp1",
                            "ess_sum_root_height_c5sp1",
                            "ess_sum_root_height_c6sp1",
                    ],
                    "label": "Effective sample size of event time",
                    "short_label": "ESS of time",
                    "ndigits": 1,
                    "mixed_comp_only": True,
                    "exclude_pairs_only": False,
                    "six_pops_extra_headers": [],
                    "center_key": "mean",
            },
            "ess-root-pop-size": {
                    "headers": [
                            "ess_sum_pop_size_root_c1sp1",
                            "ess_sum_pop_size_root_c2sp1",
                            "ess_sum_pop_size_root_c3sp1",
                    ],
                    "label": "Effective sample size of ancestral population size",
                    "short_label": "ESS of ancestral size",
                    "ndigits": 1,
                    "mixed_comp_only": False,
                    "exclude_pairs_only": True,
                    "six_pops_extra_headers": [
                            "ess_sum_pop_size_root_c4sp1",
                            "ess_sum_pop_size_root_c5sp1",
                            "ess_sum_pop_size_root_c6sp1",
                    ],
                    "center_key": "mean",
            },
            "ess-div-root-pop-size": {
                    "headers": [
                            "ess_sum_pop_size_root_c4sp1",
                            "ess_sum_pop_size_root_c5sp1",
                            "ess_sum_pop_size_root_c6sp1",
                    ],
                    "label": "Effective sample size of ancestral population size",
                    "short_label": "ESS of ancestral size",
                    "ndigits": 1,
                    "mixed_comp_only": True,
                    "exclude_pairs_only": False,
                    "six_pops_extra_headers": [],
                    "center_key": "mean",
            },
            "psrf-ln-likelihood": {
                    "headers": [
                            "psrf_ln_likelihood",
                    ],
                    "label": "PSRF of log likelihood",
                    "short_label": "PSRF of lnL",
                    "ndigits": 2,
                    "mixed_comp_only": False,
                    "exclude_pairs_only": False,
                    "six_pops_extra_headers": [],
                    "center_key": "mean",
            },
            "psrf-event-time": {
                    "headers": [
                            "psrf_root_height_c1sp1",
                            "psrf_root_height_c2sp1",
                            "psrf_root_height_c3sp1",
                    ],
                    "label": "PSRF of event time",
                    "short_label": "PSRF of time",
                    "ndigits": 2,
                    "mixed_comp_only": False,
                    "exclude_pairs_only": True,
                    "six_pops_extra_headers": [
                            "psrf_root_height_c4sp1",
                            "psrf_root_height_c5sp1",
                            "psrf_root_height_c6sp1",
                    ],
                    "center_key": "mean",
            },
            "psrf-div-time": {
                    "headers": [
                            "psrf_root_height_c4sp1",
                            "psrf_root_height_c5sp1",
                            "psrf_root_height_c6sp1",
                    ],
                    "label": "PSRF of event time",
                    "short_label": "PSRF of time",
                    "ndigits": 2,
                    "mixed_comp_only": True,
                    "exclude_pairs_only": False,
                    "six_pops_extra_headers": [],
                    "center_key": "mean",
            },
            "run-time": {
                    "headers": [
                            "mean_run_time",
                    ],
                    "label": "Run time (seconds)",
                    "short_label": "Run time (seconds)",
                    "ndigits": 1,
                    "mixed_comp_only": False,
                    "exclude_pairs_only": False,
                    "six_pops_extra_headers": [],
                    "center_key": "median",
            },
    }

    for parameter, p_info in histograms_to_plot.items():
        data = {}
        var_only_data = {}
        if p_info["mixed_comp_only"]:
            for sim_dir, r in results_mixed_comps.items():
                data[sim_dir] = HistogramData.init(r, p_info["headers"], False)
            for sim_dir, r in var_only_results_mixed_comps.items():
                var_only_data[sim_dir] = HistogramData.init(r, p_info["headers"], False)
        else:
            for sim_dir, r in results.items():
                data[sim_dir] = HistogramData.init(r, p_info["headers"], False)
            for sim_dir, r in var_only_results.items():
                var_only_data[sim_dir] = HistogramData.init(r, p_info["headers"], False)

            for sim_dir, r in results_6_pops.items():
                data[sim_dir] = HistogramData.init(r,
                        p_info["headers"] + p_info["six_pops_extra_headers"],
                        False)
            for sim_dir, r in var_only_results_6_pops.items():
                var_only_data[sim_dir] = HistogramData.init(r,
                        p_info["headers"] + p_info["six_pops_extra_headers"],
                        False)

        if not p_info["exclude_pairs_only"]:
            for sim_dir, r in results_pairs.items():
                data[sim_dir] = HistogramData.init(r, p_info["headers"], False)
            for sim_dir, r in var_only_results_pairs.items():
                var_only_data[sim_dir] = HistogramData.init(r, p_info["headers"], False)

        for sim_dir in data.keys():
            prefix = get_prefix_from_sim_dir_name(sim_dir)
            generate_specific_histogram(
                    data = data[sim_dir],
                    plot_file_prefix = parameter + "-" + prefix,
                    title = None,
                    title_size = 16.0,
                    x_label = None,
                    x_label_size = 16.0,
                    y_label = None,
                    y_label_size = 16.0,
                    plot_width = plot_width,
                    plot_height = plot_height,
                    pad_left = pad_left,
                    pad_right = pad_right,
                    pad_bottom = pad_bottom,
                    pad_top = pad_top,
                    bins = None,
                    range_key = "range",
                    center_key = p_info["center_key"],
                    number_of_digits = p_info["ndigits"])
            if sim_dir in var_only_data:
                generate_specific_histogram(
                        data = var_only_data[sim_dir],
                        plot_file_prefix = "var-only-" + parameter + "-" + prefix,
                        title = None,
                        title_size = 16.0,
                        x_label = None,
                        x_label_size = 16.0,
                        y_label = None,
                        y_label_size = 16.0,
                        plot_width = plot_width,
                        plot_height = plot_height,
                        pad_left = pad_left,
                        pad_right = pad_right,
                        pad_bottom = pad_bottom,
                        pad_top = pad_top,
                        bins = None,
                        range_key = "range",
                        center_key = p_info["center_key"],
                        number_of_digits = p_info["ndigits"])


    # root_gamma_labels_opt = []
    # for sim_dir in sim_dirs_opt:
    #     root_gamma_labels_opt.append(get_root_gamma_label(sim_dir))

    # root_gamma_labels_opt_centered = []
    # for sim_dir in sim_dirs_opt_centered:
    #     root_gamma_labels_opt_centered.append(get_root_gamma_label(sim_dir))

    # root_gamma_labels_opt_all = []
    # for sim_dir in sim_dirs_opt_all:
    #     root_gamma_labels_opt_all.append(get_root_gamma_label(sim_dir))

    # root_gamma_labels = []
    # for sim_dir in sim_dirs:
    #     root_gamma_labels.append(get_root_gamma_label(sim_dir))
    
    # results_opt_all = []
    # for sim_dir in sim_dirs_opt_all:
    #     results_opt_all.append(
    #             parse_results(glob.glob(os.path.join(project_util.VAL_DIR,
    #                     sim_dir,
    #                     "batch00?",
    #                     "results.csv.gz")))
    #             )

    # var_only_results_opt_all = []
    # for sim_dir in sim_dirs_opt_all:
    #     var_only_results_opt_all.append(
    #             parse_results(glob.glob(os.path.join(project_util.VAL_DIR,
    #                     sim_dir,
    #                     "batch00?",
    #                     "var-only-results.csv.gz")))
    #             )

    # results_opt = [
    #         results_opt_all[0],
    #         results_opt_all[2],
    #         results_opt_all[1],
    #         ]
    # var_only_results_opt = [
    #         var_only_results_opt_all[0],
    #         var_only_results_opt_all[2],
    #         var_only_results_opt_all[1],
    #         ]
    # results_opt_centered = results_opt_all[2:]
    # var_only_results_opt_centered = var_only_results_opt_all[2:]

    # row_labels = [
    #         "All sites",
    #         "Variable-only sites",
    # ]

    # t_vs_e_opt_all = [ScatterData(
    #         **dict(zip(
    #                 ("x", "y", "highlight_values", "highlight_threshold"),
    #                 list(get_coal_units_vs_error(
    #                         r,
    #                         height_parameters,
    #                         leaf_size_parameters)) + [1.1]
    #                 ))
    #         ) for r in results_opt_all]
    # vo_t_vs_e_opt_all = [ScatterData(
    #         **dict(zip(
    #                 ("x", "y", "highlight_values", "highlight_threshold"),
    #                 list(get_coal_units_vs_error(
    #                         r,
    #                         height_parameters,
    #                         leaf_size_parameters)) + [1.1]
    #                 ))
    #         ) for r in var_only_results_opt_all]

    # t_vs_psrf_opt_all = [ScatterData(
    #         **dict(zip(
    #                 ("x", "y"),
    #                 list(get_coal_units_vs_error(
    #                         r,
    #                         height_parameters,
    #                         leaf_size_parameters,
    #                         psrf_as_response = True))
    #                 ))
    #         ) for r in results_opt_all]
    # vo_t_vs_psrf_opt_all = [ScatterData(
    #         **dict(zip(
    #                 ("x", "y"),
    #                 list(get_coal_units_vs_error(
    #                         r,
    #                         height_parameters,
    #                         leaf_size_parameters,
    #                         psrf_as_response = True))
    #                 ))
    #         ) for r in var_only_results_opt_all]

    # generate_scatter_plots(
    #         data_grid = [t_vs_e_opt_all, vo_t_vs_e_opt_all],
    #         plot_file_prefix = "opt-all-coal-units-vs-error",
    #         column_labels = root_gamma_labels_opt_all,
    #         row_labels = row_labels,
    #         plot_width = 1.9,
    #         plot_height = 1.8,
    #         pad_left = 0.08,
    #         pad_right = 0.98,
    #         pad_bottom = 0.14,
    #         pad_top = 0.94,
    #         x_label = "True event time in coalescent units",
    #         x_label_size = 18.0,
    #         y_label = "Relative error",
    #         y_label_size = 18.0,
    #         force_shared_x_range = False,
    #         force_shared_y_range = False,
    #         force_shared_xy_ranges = False,
    #         force_shared_spines = False,
    #         include_coverage = False,
    #         include_rmse = False,
    #         include_identity_line = False,
    #         include_error_bars = False)

    # generate_scatter_plots(
    #         data_grid = [t_vs_psrf_opt_all, vo_t_vs_psrf_opt_all],
    #         plot_file_prefix = "opt-all-coal-units-vs-psrf",
    #         column_labels = root_gamma_labels_opt_all,
    #         row_labels = row_labels,
    #         plot_width = 1.9,
    #         plot_height = 1.8,
    #         pad_left = 0.08,
    #         pad_right = 0.98,
    #         pad_bottom = 0.14,
    #         pad_top = 0.94,
    #         x_label = "True event time in coalescent units",
    #         x_label_size = 18.0,
    #         y_label = "PSRF",
    #         y_label_size = 18.0,
    #         force_shared_x_range = False,
    #         force_shared_y_range = False,
    #         force_shared_xy_ranges = False,
    #         force_shared_spines = False,
    #         include_coverage = False,
    #         include_rmse = False,
    #         include_identity_line = False,
    #         include_error_bars = False)


if __name__ == "__main__":
    main_cli()
