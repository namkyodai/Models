#!/usr/bin/python -tt
# -*- coding: utf-8 -*-
# =============================================================================
# File      : rlstpe.py 
# Creation  : 03 Jul 2015
# Time-stamp: <Don 2015-10-15 17:40 juergen>
#
# Copyright (c) 2015 Jürgen Hackl <hackl@ibi.baug.ethz.ch>
#               http://www.ibi.ethz.ch
# $Id$ 
#
# Description : Restricted Least Square Transition Probability Estimator
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
# =============================================================================

# Import Modules
# ==============

import numpy as np
import cvxpy as opt


# Code
# ====

def load_data(filename):
    """ Load aggregate data from input file

    Aggregate data are loaded from a csv file. There by each column represents
    a condition state, each row a time step and each entry a sample
    proportion. Header and row names are not required.

    +-------+-------+-------+-----+-------+
    |  Time |   S1  |   S2  | ... |   Sr  |
    +=======+=======+=======+=====+=======+
    |  t=0  | y1(0) | y2(0) | ... | yr(0) |
    +-------+-------+-------+-----+-------+
    |  t=1  | y1(1) | y2(1) | ... | yr(1) |
    +-------+-------+-------+-----+-------+
    |  t=t  | y1(t) | y2(t) | ... | yr(t) |
    +-------+-------+-------+-----+-------+
    |  t=T  | y1(T) | y2(T) | ... | yr(T) |
    +-------+-------+-------+-----+-------+

    Parameters
    ----------
    filename : file or string
       File, or filename to read. Here, a file extension IS necessary!

    Returns
    -------
    data : numpy matrix
      Matrix where each row represents a condition states; i.e. data[0] = y1

    Examples
    --------
    >>> rlstpe.load_data('inp_test.csv') # doctest: +SKIP
    """
    # read data from csv file
    data = np.genfromtxt(filename, delimiter=',')

    # return data in transposed format
    # i.e. data[0] is equal to y0
    return data.transpose()

def rlstpe(data):
    """Restricted least square transition probability estimator

    Calculates the restricted least square transition probability estimator
    based on the proposed model in [1]_. Thereby the objective function

    .. math::
       u^T u = (y-Xp)^T (y-Xp)

    is minimized subjected to the constrains

    .. math::
       Gp &= n\\
       p &\leq 0\\
       x_ij &= 0 if i>j

    where :math:`G` is an :math:`(r \times r^2)` known coefficient matrix
    :math:`[I_1,I_2,\dots I_r]` with each :math:`I_i` an :math:`(r \times r)`
    identity matrix and :math:`n` is an :math:`(r \times 1)` column vector
    with all entries equal one [1]_.

    Parameters
    ----------
    data : numpy matrix
      Matrix where each row represents a condition states; i.e. data[0] = y1

    Returns
    -------
    p the estimated transition probability : numpy matrix
      transition probability matrix p where :math:`p_{ij}` is the constant transition
      probability associated wit a change from state :math:`s_i` to :math:`s_j`.

    Note
    ----
    A change of the time period is at the moment not possible. Now the
    calculation goes form :math:`t=0` till :math:`t=T`, thereby :math:`t` is
    fixed by the initial data which are used.

    Examples
    --------
    >>> rlstpe.rlstpe(data) # doctest: +SKIP

    References
    ----------
    .. [1] Lee T. C., Judge G. G. and Zellner A. (1977) Estimating the
    Parameters of the Markov Probability Model from Aggregate Time Series
    Data, 2nd Edition, North-Holland, Amsterdam, Netherlands.

    """

    # get number of Markov stages
    r = data.shape[0]

    # get number of time periods
    T = data.shape[1]

    # control if T is strictly grater than r (Miller, 1952)
    assert T > r, 'T has to be strictly grater than r (Miller, 1952)'

    # setup vector of sample distributions
    y = np.bmat([data[i,1:T] for i in range(r)]).transpose()

    # setup matrix of observed proportions
    Xj = np.column_stack([data[i,0:T-1] for i in range(r)])

    # setup zero matrix
    X0 = np.zeros(Xj.shape)

    # setup identiy matrix
    Ij = np.identity(r)

    # create block matrix X
    block = []
    for i in range(r):
        row = []
        for j in range(r):
            if Ij[i,j] == 1:
                row.append(Xj)
            else:
                row.append(X0)
        block.append(row)
    X = np.bmat(block)

    # create coefficient matrix G
    block = []
    for i in range(r):
        row = []
        for j in range(r):
            row.append(Ij)
        block.append(row)
    G = np.bmat(block)

    # setup unit vector n
    n = np.ones(r*r).reshape(-1,1)

    # set up vector z
    # as a constrain to ensure that the condition state change only to the
    # next lower condition state. 
    z = np.triu(np.ones((r,r))).transpose().reshape(r*r)

    # Optimization
    x = opt.Variable(r*r)
    objective = opt.Minimize(opt.sum_entries(opt.square(y - X*x)))
    constraints = [0 <= x, G*x == n, x<=z]
    prob = opt.Problem(objective, constraints)

    # solve the problem
    # prob.solve()
    prob.solve(verbose=True)

    # format the transition probability matrix
    p = np.round(x.value,4).reshape(r,r).transpose()

    return p


# Example
# =======

# Aggregate data from a simulation study of 1000 individuals
# Table 4.1 p 47
# data = load_data('inp_sim.csv')

# Market shares of Camels, Lucky Strike and Chesterfield
# Table 4.8 p 60
data = load_data('bridge.csv')


# Restricted least square transition probability estimation
p = rlstpe(data)

print('Result of p:')

# The real transition probability matrix of the simulation study can be fond
# on p 44, The estimated solution on p 55.
# The estimated solution of the cigarette example can be found on p 61.

print(p)


# =============================================================================
# eof
#
# Local Variables: 
# mode: python
# End: 

 
