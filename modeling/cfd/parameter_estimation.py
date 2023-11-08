# -*- coding: utf-8 -*-
"""

Summary
-------
Sensitivity analysis and optimization (round-trip) based on an analytical
jet profile.

4 referenece beams are creteate for both analyses.

From an abritrary starting point, the optimizer attempts to recover
the original input parameters:

    1. T0 - jet exit temperature
    2. C0 - jet exit concentration
    3. a  - width growth shape parameter
    4. b  - width growth rate parameter

The absorptions is assumed to composed of the sum of two-gaussians.


Notes
-----
1. the optimization seems to work better using the logrithms of the parameters

TODO
----
1. rename logTau => kappaL for consistency with other files
2. fix unit consistency

"""

import scipy.optimize as optimize
import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt

import jet_profiles
from ray_integration import GaussAbsorption


def linear_to_log_space(ns_in=50, w_in=0.1, do_plot=True):
    """create point distribution over the internval [0,1], the interval
    is split into three internvals [0,-(w-1)/2,(w+1)/2,1].

    Parameters
    ----------
    ns_in : integer
        number of points in the middle region
    w_in : float < 1
        width of the inner region

    In the other outer regio has logrighmic spacing with the number of points
    selected to have continuous spacing at the interface between the intervals

    Returns
    -------
    s : float array
        array of points on the interval [0,1]

    """


    s_in = np.linspace(-w_in, w_in, ns_in)*0.5
    ds_in = np.diff(s_in)

    start = np.log10(s_in[-1])
    end = np.log10(0.5)

    print("start/end",start, end)
    def f_out(ns_out):
        s_out = np.logspace(start, end, ns_out)
        ds_out = np.diff(s_out)

        return s_out, ds_out[0]/ds_in[-1]

    L = (1-w_in)*0.5
    ns_out = np.arange(2,int(L/ds_in[-1]))
    ratio = np.array([f_out(n)[-1] for n in ns_out])

    if do_plot:
        fig, ax = plt.subplots(2,1)
        ax[0].semilogx(ns_out, ratio)
        ax[0].set_xlabel("# outer grid points")
        ax[0].set_ylabel("spacing ratio")
        ax[0].loglog(ns_out, len(ns_out)*[1.0], '--')

    ii = np.argmin( abs(ratio-1.0))
    print(ii, ns_out[ii])
    s_out = f_out( ns_out[ii] )[0]

    s_all = np.concatenate( [-s_out[1:],s_in,s_out[1:]] )
    #s_all = np.unique(s_all)
    s_all = np.sort(s_all)

    if do_plot:
        sc = 0.5*(s_all[1:]+s_all[:-1])
        ax[1].plot(sc, np.diff(s_all), 'o')
        ax[1].set_ylabel("spacing")
        ax[1].set_xlabel("s")
        ax[0].set_title("linear_to_log_space")
        fig.tight_layout()

    return s_all

#s = linear_to_log_space()

#
# 2D plotting of the jet flow-field function
#
jet = jet_profiles.JetProfile()
jet.a = 2.0
jet.b = 6.0
jet.L = 0.5
jet.r0 = 0.006
jet.x0 = 0.0
jet.w0 = 1.0

x_jet = jet.L*np.linspace(0,1,63)
r_jet = jet.r0*np.concatenate([np.linspace(0,1), np.logspace(0,1.0)])

p1, p2 = np.meshgrid(x_jet, r_jet)
pos_jet = np.array([p1, p2])
pos_jet = np.array([p1, p2]).transpose((2,1,0))
f_jet = jet.value(pos_jet[:,:,0], pos_jet[:,:,1])

do_plot = True
if do_plot:
    fig, ax = plt.subplots()
    ax.contour(pos_jet[:,:,0], pos_jet[:,:,1], f_jet[0])

#
# forward calculation
#

x_CL = np.array([0.20, 0.0, 0.0])
x_offset = np.array([-0.1, 0.0, -0.3])
n_dim = 3

def calc_beam_pos(x_CL, x_offset, n_s=100, s_near=0.2):
    x_source = x_CL + x_offset
    x_target = x_CL - x_offset
    s = np.linspace(0, 1, n_s)
    #s = linear_to_log_space(n_s, w_in=0.1, do_plot=False)
    return x_source + (x_target-x_source)*s[:,np.newaxis]

def calc_beam_profile(jet_func, pos, fields={"C":[0,1], "T":[300,2000]}, do_plot=True):
    """calculate flow profile along a beam

    Parameters
    ----------
    jet_func : FlowField function
        calculates flow field at
    pos_beam : 2D float array [i_s,i_dim]
        position vector of beam
    fields : dictionary of 2-tuples
        keys - field names
        values[0] - free stream reference value
        values[1] - jet reference value
    do_plot : boolean
        plot fields

    """

    L = np.sum( (pos[-1,:] - pos[0,:])**2 )**0.5
    s = np.sum( (pos[:,:] - pos[0,:])**2, axis=1 )**0.5
    s /= L
    x_beam = pos[:,0]
    r_beam = (pos[:,1]**2 + pos[:,2]**2)**0.5

    f_beam = jet_func.value(x_beam, r_beam)
    #print("calc_beam_profile",jet_func.a, jet_func.b)
    #
    r1_beam = np.where(pos_beam[:,2] > 0, r_beam, -r_beam)

    if do_plot:
        fig, ax = plt.subplots()
        ax.plot(x_beam, r_beam, 'k--', label="beam")
        ax.plot(x_beam[0], r_beam[0], 'o', label="start")
        ax.plot(x_beam[-1], r_beam[-1], 'o', label="start")
        ax.contour(pos[:,:,0], pos[:,:,1], f_jet[0])
        ax.legend()
        f_beam = jet_func.value(x_beam, r_beam)
        fig, axArray = plt.subplots(2,3,sharey=True)
        ax = axArray.flat
        names = "s x r r* y z".split()
        for i, x in enumerate([s,x_beam,r_beam, r1_beam, pos_beam[:,1], pos_beam[:,2]]):
            ax[i].plot(x, f_beam[0], 'o')
            ax[i].set_xlabel(names[i])
            ax[i].set_ylabel("f")

    d = {"f":f_beam[0], "L":L, "s":s, "x":x_beam, "r":r_beam, "pos":pos}
    #print(fields)
    for name, v_ref in fields.items():
        d[name] = v_ref[0] + (v_ref[1] - v_ref[0])*f_beam[0]

    return d


def calc_absorption(f_abs, flow_field, do_plot=True):
    """

    Parameters
    ----------
    f_abs : Absorption class or function
        DESCRIPTION.
    flow_field : dictionary of float arrays
        s : beam coordinate
        T : temperature flow array
        C : concentration

    Returns
    -------
    dict
        DESCRIPTION.

    """

    T = flow_field["T"]
    C_K = flow_field["C"]
    L = flow_field["L"]
    s = flow_field["s"]

    Q_abs = np.array([f_abs(T, lam) for lam in lam_range])

    #print(kappa_abs)
    kappa_abs = Q_abs*C_K*f_abs.Q0

    if do_plot:
        fig, axArray = plt.subplots(2,3, sharex="col")
        ax = axArray.flat
        ax[0].plot(s, C_K)
        ax[3].plot(s, T)
        ax[3].set_xlabel("beam postion")
        ax[0].set_ylabel("K concentration")
        ax[3].set_ylabel("T")
        n_lam = len(lam_range)
        cmap = plt.colormaps["viridis"]
        for i, lam in enumerate(lam_range):
            color = cmap(i/float(n_lam))
            ax[1].semilogy(s, Q_abs[i,:], color=color)
            ax[4].semilogy(s, kappa_abs[i,:], color=color)

        ax[4].set_xlabel("beam postion")
        ax[1].set_ylabel("kappa")
        ax[4].set_ylabel("Q")

        for i, si in enumerate(s):
            if i % 20 == 0:
                color = cmap(si/s.max())
                ax[2].semilogy(lam_range, Q_abs[:,i], color=color)
                ax[5].semilogy(lam_range, kappa_abs[:,i], color=color)
        ax[5].set_xlabel("frequency")
        ax[5].set_ylabel("kappa")
        ax[2].set_ylabel("Q")
        fig.tight_layout()


    logI_target = L*integrate.trapezoid(kappa_abs, s, axis=1)
    logI =  L*integrate.cumulative_trapezoid(kappa_abs, s, axis=1)
    I = np.exp(-logI)
    I_target = np.exp(-logI_target)

    if do_plot:
        fig, axArray = plt.subplots(2,2,sharex=True)
        ax = axArray.flat

        ax[0].plot(lam_range, logI_target)

        ax[2].plot(lam_range, I_target)

        c = ax[1].contour(lam_range, s[1:], I.T)
        plt.colorbar(c) #, cax=ax[1])
        c = ax[3].contour(lam_range, s[1:], logI.T)
        plt.colorbar(c) #, cax=ax[1])
        fig.tight_layout()
        ax[2].set_xlabel("wavelength [nm]")
        ax[3].set_xlabel("wavelength [nm]")
        labels = ["-logI_target","beam pos.","I_target","beam pos"]
        for i, n in enumerate(labels):
            ax[i].set_ylabel(n)

    return {"tau_t":I_target, "logTau_t":logI_target, "logTau":-logI, "pos":pos_beam}


#
# create a "bundle of" four reference beams
#

pos_beam = calc_beam_pos(x_CL, x_offset, n_s=200)

# x, y spacing
dx = 0.05
dy = 0.01

I_ref = []

f_abs = GaussAbsorption()
f_abs.b[:] = 0.5
f_abs.Q0 = 1e6
lam_range = np.arange(760,775+0.1,0.5)
n_lam = len(lam_range)

beam_list = []

fields0={"C":[0,1e-4], "T":[300,2000]}
for i_x, j_y in [(0,0),(1,0),(0,1),(-1,0)]:
    pos = pos_beam + np.array([i_x*dx, j_y*dy, 0.0])[np.newaxis,:]
    flow = calc_beam_profile(jet, pos, fields=fields0, do_plot=False)
    beam = calc_absorption(f_abs, flow, do_plot=False)
    beam.update({"i_x":i_x,"i_y":j_y,"dx":dx,"dy":dy})
    beam.update(flow)
    beam_list.append(beam)

for i, beam in enumerate(beam_list):
    ax.plot(beam["x"], beam["r"], '--', label="{i_x:} {i_y:}".format_map(beam), alpha=0.5)
ax.legend()
ax.set_aspect(1.0)

n_beam = len(beam_list)
fig, ax = plt.subplots()
for i, beam in enumerate(beam_list):
    pos = beam["pos"]
    ax.semilogy(lam_range, beam["logTau_t"], label="{i_x:} {i_y:}".format_map(beam), alpha=0.5)

ax.legend()
ax.set_ylabel("wavelength [nm]")
ax.set_xlabel("logTau @ target")


#
# sentivitity analysis
#

def calc_new_beams(p_in):
    """calculate new beams based on different jet flow field parameters, *p_in*"""
    C_0, T_0 = p_in["C"], p_in["T"]
    jet.a, jet.b = p_in["a"], p_in["b"]
    new_list = []
    for i, beam in enumerate(beam_list):
        pos = beam["pos"]*1.0
        flow_1 = calc_beam_profile(jet, pos, fields={"C":[0,C_0], "T":[300,T_0]}, do_plot=False)
        beam_1 = calc_absorption(f_abs, flow_1, do_plot=False)
        beam_1.update(flow_1)
        new_list.append(beam_1)
    return new_list


p_names = "T C a b".split()
pref = {k:v[1] for k, v in fields0.items()}
pref.update( {"a":jet.a, "b":jet.b})

figs = []
ax = []
for i in range(3):
    fig, axArray = plt.subplots(2,2)
    a = axArray.flat
    figs.append(fig)
    ax.append(a)

p1 = pref.copy()
eps = 0.1
small = 1e-12
for i, name in enumerate(p_names):
    p1[name] = (1.0 + eps)*pref[name]
    new_beams = calc_new_beams(p1)
    p1[name] = pref[name]
    for j, beam in enumerate(beam_list):

        label = "{i_x:} {i_y:}".format_map(beam)
        #line = ax[i].semilogy(lam_range, beam["logTau_t"], "--") #, label="{i_x:} {i_y:}".format_map(beam), alpha=0.3)
        #color = line[0].get_color()
        val = new_beams[j]["logTau_t"]/beam["logTau_t"]
        ax[0][i].plot(lam_range, val, label=label) #, alpha=0.8) # color=color)
        ax[0][i].set_ylabel("logTau_t ratio")
        ax[0][i].set_xlabel("wavelength [nm]")

        for i_var, var_name in enumerate("T C".split()):
            val = (new_beams[j][var_name] + small)/(beam[var_name] + small)
            ax[i_var+1][i].plot(new_beams[j]["s"], val, label=label)
            ax[i_var+1][i].set_ylabel(var_name)
            ax[i_var+1][i].set_xlabel("beam coordinate")

    label = name + " sentivity"
    for a in ax:
        a[i].set_title(label)

    for i_var in [1,2]:
        ax[i_var][i].set_xlim(0.4,0.6)

    label = name + " sentivity"

for a in ax:
    i = 0
    a[i].legend(title="source")

for f in figs:
    f.tight_layout()


#
# optimization
#
wavelength_mask = np.ones(n_lam)

def opt_func(p, var_name):
    #print(p)
    C_0, T_0 = p[0], p[1]
    a, b = p[2], p[3]
    jet.a = a
    jet.b = b

    f_err = 0.0
    for i, beam in enumerate(beam_list):
        pos = beam["pos"]
        #print(i, pos[0,0])
        flow_1 = calc_beam_profile(jet, pos, fields={"C":[0,C_0], "T":[300,T_0]}, do_plot=False)
        beam_1 = calc_absorption(f_abs, flow_1, do_plot=False)
        f = beam_1[var_name]
        f_ref = beam_list[i][var_name]
        err = (f_ref - f)/f.max()
        f_err += np.sum( err**2*wavelength_mask )

    return f_err*1e8

def log_opt_func(q, var_name):
    return opt_func(np.exp(q), var_name)


p0 = pref.copy()
p0["C"] = 1e-6
p0["T"] = 1500.0
p0["a"] = 1.0
p0["b"] = 1.0

bounds = [(1e-10,1),(500,3500),(1e-3,10),(1e-3,10)]
log_bounds = np.array([np.log(b) for b in bounds])


log_p0_vals = 0.5*(log_bounds[:,0]+log_bounds[:,1])

p0 = dict( zip(p0.keys(), np.exp(log_p0_vals)))

p0_vals = np.array(list(p0.values()))

print("Solving optimization problem")
print("initial paramters", p0)

#
# use log variables for optimization
#
do_log = True
if do_log:
    opt_out = optimize.minimize(log_opt_func, log_p0_vals, args=("tau_t"), bounds=log_bounds, tol=1e-9)
    p1 = dict( zip(p0.keys(), np.exp(opt_out.x)))
else:
    opt_out = optimize.minimize(opt_func, p0_vals, args=("tau_t"), bounds=bounds, tol=1e-9)
    p1 = dict( zip(p0.keys(), opt_out.x))

print("\nDone")
print(opt_out)
#p1["C"] = np.exp(p0["logC"])

print("\nParameters")
print("initial", p0)
print("solution ", p1)
print("reference", pref)

new_beams = calc_new_beams(p1)
for var in ["C","T"]:
    fig, ax = plt.subplots(2,2)
    ax = ax.flat
    for i in range(len(new_beams)):
        label="x{i_x:} y{i_y:}".format_map(beam_list[i])
        ax[i].plot(new_beams[i]["s"], new_beams[i][var], label="opt")
        ax[i].plot(new_beams[i]["s"], beam_list[i][var], '--', label="analytic")
        ax[i].set_title(i)
        ax[i].set_xlim(0.4,0.6)
        ax[i].set_ylabel(var)
        ax[i].set_xlabel("s")
        ax[i].set_title(label)
    ax[0].legend()
    fig.tight_layout()

fig, ax = plt.subplots(2,2)
ax = ax.flat
for i, beam in enumerate(beam_list):
    label="x{i_x:} y{i_y:}".format_map(beam_list[i])
    ax[i].plot(lam_range, new_beams[i]["logTau_t"], label="opt")
    ax[i].plot(lam_range, beam_list[i]["logTau_t"], '--', label="analytic")
    ax[i].set_ylabel("-logTau")
    ax[i].set_xlabel("wavelength [nm]")
ax[0].legend()
fig.tight_layout()
