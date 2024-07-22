#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 17:23:18 2023

@author: David Kaltenbrunner
The following functions are copies from other python packages with
modifications to be usefull for plotting the way we want it.
The original package name is noted in each function.
"""

from bxa.xspec.solver import set_parameters, XSilence
import numpy as np
import scipy.stats
from xspec import Plot
from numpy import log10


def posterior_predictions_plot(analyser, plottype, nsamples=None, group=1):
    """From bxa.xspec.solver
    Internal Routine used by posterior_predictions_unconvolved, posterior_predictions_convolved
    """
    # for plotting, we don't need so many points, and especially the
    # points that barely made it into the analysis are not that interesting.
    # so pick a random subset of at least nsamples points
    posterior = analyser.posterior[:nsamples]

    with XSilence():
        olddevice = Plot.device
        Plot.device = '/null'

        # plot models
        maxncomp = 100 if Plot.add else 0
        for k, row in enumerate(posterior):
            set_parameters(
                values=row, transformations=analyser.transformations)
            Plot.setRebin(0, 1)
            Plot(plottype)
            # get plot data
            if plottype == 'model':
                base_content = np.transpose([
                    Plot.x(group), Plot.xErr(group), Plot.model(group)])
            elif Plot.background:
                base_content = np.transpose([
                    Plot.x(group), Plot.xErr(group), Plot.y(
                        group), Plot.yErr(group),
                    Plot.backgroundVals(group), np.zeros_like(
                        Plot.backgroundVals(group)),
                    Plot.model(group)])
            else:
                base_content = np.transpose([
                    Plot.x(group), Plot.xErr(group), Plot.y(
                        group), Plot.yErr(group),
                    Plot.model(group)])
            # get additive components, if there are any
            comp = []
            for i in range(1, maxncomp):
                try:
                    comp.append(Plot.addComp(i))
                except Exception:
                    print(
                        'The error "***XSPEC Error: Requested array does not exist for this plot." can be ignored.')
                    maxncomp = i
                    break
            content = np.hstack((base_content, np.transpose(
                comp).reshape((len(base_content), -1))))
            yield content
        Plot.device = olddevice


class PredictionBand(object):
    """From ultranest.plot
    Plot bands of model predictions as calculated from a chain.

    call add(y) to add predictions from each chain point

    Example::

        x = numpy.linspace(0, 1, 100)
        band = PredictionBand(x)
        for c in chain:
            band.add(c[0] * x + c[1])
        # add median line
        band.line(color='k')
        # add 1 sigma quantile
        band.shade(color='k', alpha=0.3)
        # add wider quantile
        band.shade(q=0.01, color='gray', alpha=0.1)
        plt.show()

    Parameters
    ----------
    x: array
        The independent variable

    """

    def __init__(self, x, ax, shadeargs={}, lineargs={}):
        """Initialise with independent variable *x*."""
        self.x = x
        self.ys = []
        self.ax = ax
        self.shadeargs = shadeargs
        self.lineargs = lineargs

    def add(self, y):
        """Add a possible prediction *y*."""
        self.ys.append(y)

    def set_shadeargs(self, **kwargs):
        """Set matplotlib style for shading."""
        self.shadeargs = kwargs

    def set_lineargs(self, **kwargs):
        """Set matplotlib style for line."""
        self.lineargs = kwargs

    def get_line(self, q=0.5):
        """Over prediction space x, get quantile *q*. Default is median."""
        if not 0 <= q <= 1:
            raise ValueError("quantile q must be between 0 and 1, not %s" % q)
        assert len(self.ys) > 0, self.ys
        return scipy.stats.mstats.mquantiles(self.ys, q, axis=0)[0]

    def shade(self, q=0.341, **kwargs):
        """Plot a shaded region between 0.5-q and 0.5+q. Default is 1 sigma."""
        if not 0 <= q <= 0.5:
            raise ValueError(
                "quantile distance from the median, q, must be between 0 and 0.5, not %s. For a 99% quantile range, use q=0.48." % q)
        shadeargs = dict(self.shadeargs)
        shadeargs.update(kwargs)
        lo = self.get_line(0.5 - q)
        hi = self.get_line(0.5 + q)
        return self.ax.fill_between(self.x, lo, hi, **shadeargs)

    def line(self, **kwargs):
        """Plot the median curve."""
        lineargs = dict(self.lineargs)
        lineargs.update(kwargs)
        mid = self.get_line(0.5)
        return self.ax.plot(self.x, mid, **lineargs)

def modded_create_uniform_prior_for(model, par):
	"""
    From bxa, change usage of pmin/pmax to pbottom/ptop
	Use for location variables (position)
	The uniform prior gives equal weight in non-logarithmic scale.
	"""
	pval, pdelta, pmin, pbottom, ptop, pmax = par.values
	print('  uniform prior for %s between %f and %f ' % (par.name, pbottom, ptop))
	# TODO: should we use min/max or bottom/top?
	low = float(pbottom)
	spread = float(ptop - pbottom)
	if pbottom > 0 and ptop / pbottom > 100:
		print('   note: this parameter spans several dex. Should it be log-uniform (create_jeffreys_prior_for)?')
	def uniform_transform(x): return x * spread + low
	return dict(model=model, index=par._Parameter__index, name=par.name, 
		transform=uniform_transform, aftertransform=lambda x: x)

def modded_create_jeffreys_prior_for(model, par):
	"""
    From bxa, change usage of pmin/pmax to pbottom/ptop
	Use for scale variables (order of magnitude)
	The Jeffreys prior gives equal weight to each order of magnitude between the
	minimum and maximum value. Flat in logarithmic scale
	"""
	pval, pdelta, pmin, pbottom, ptop, pmax = par.values
	# TODO: should we use min/max or bottom/top?
	#print '  ', par.values
	print('  jeffreys prior for %s between %e and %e ' % (par.name, pbottom, ptop))
	if pbottom == 0:
		raise Exception('You forgot to set reasonable parameter limits on %s' % par.name)
	low = log10(pbottom)
	spread = log10(ptop) - log10(pbottom)
	if spread > 10:
		print('   note: this parameter spans *many* dex. Double-check the limits are reasonable.')
	def log_transform(x): return x * spread + low
	def log_after_transform(x): return 10**x
	return dict(model=model, index=par._Parameter__index, name='log(%s)' % par.name, 
		transform=log_transform, aftertransform=log_after_transform)

def modded_create_gaussian_prior_for(model, par, mean, std):
	"""
	Use for informed variables.
	The Gaussian prior weights by a Gaussian in the parameter.
	"""
	import scipy.stats
	pval, pdelta, pmin, pbottom, ptop, pmax = par.values
	rv = scipy.stats.norm(mean, std)
	def gauss_transform(x): 
		return max(pbottom, min(ptop, rv.ppf(x)))
	print('  gaussian prior for %s of %f +- %f' % (par.name, mean, std))
	return dict(model=model, index=par._Parameter__index, name=par.name, 
		transform=gauss_transform, aftertransform=lambda x: x)