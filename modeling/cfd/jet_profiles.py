# -*- coding: utf-8 -*-
"""

Summary
-------

Testing of jet profiles and derivative functions


TODO
- add citations/references for profiles

"""
import numpy as np
import matplotlib.pyplot as plt
import sympy


def test_sympy_diff():
    """sequentially calculate sensitivity derivative,

    useful for checking errors in manual function differentiation

    the particular function in question has repeated calculations in
    several of the sensitivities

    """

    def fprint(var, expr):
        print("{:}{:8} {:}".format(2*" ",repr(var),expr))

    A, r, beta = sympy.symbols("A r beta")
    f = A*sympy.exp(-(beta*r)**2)
    print("f", f)
    for p in [A,beta]:
        fprint(p, f.diff(p))
    print()

    c, w = sympy.symbols("c w")
    beta = c/w
    print("beta", beta)
    for p in [w]:
        fprint(p, beta.diff(p))
    print()

    w0 = sympy.symbols("w0")
    beta0 = c/w0
    print("beta0", beta0)
    for p in [w0]:
        fprint(p, beta0.diff(p))
    print()

    a, b, x = sympy.symbols("a b x")
    w = w0*(1 + b*x**a)
    print("w", w)
    for p in [w0, a, b]:
        fprint(p, w.diff(p))



def test_sympy_profile(p_in={"w0":1.0,"a":2,"b":6}):
    """plot function and sensitivity derivatitves using sympy"""

    x, r, x0, L, r0, w0, a, b, c = sympy.symbols("x r x0 L r0 r1 a b c", real=True)

    w = w0*(1 + b*x**a)

    beta = c/w
    beta0 = c/w0

    print("beta0",beta0)

    A = (beta/beta0)**2
    f = A*sympy.exp(-(beta*r)**2)

    df_dp = [f.diff(p) for p in [w0, a, b]]

    c_val = -np.log(0.5)
    rr = np.linspace(0.0,5.0)
    xx = 0.3
    fig, axArray = plt.subplots(2,2, sharey=True)
    ax = axArray.flat
    names = ["df_d" + x for x in "w0 a b".split()]+ ["f"]

    p = sympy.Matrix([w0, a, b])
    df_dp_vec = f.diff(p)

    for i, g in enumerate(df_dp + [f]):
        g_func =  sympy.lambdify([x,r,w0,a,b,c], g)
        vals = g_func(xx,rr,p_in["w0"],p_in["a"],p_in["b"],c_val)
        ax[i].plot(vals, rr)
        ax[i].set_xlabel(names[i])
        ax[i].set_ylabel("r/r0")

        if i < 3:
            h_func =  sympy.lambdify([x,r,w0,a,b,c], df_dp_vec[i])
            print(h_func)
            vals = h_func(xx,rr,p_in["w0"],p_in["a"],p_in["b"],c_val)
            ax[i].plot(vals, rr, '+', label="vec")
            ax[i].legend()

    ax[0].set_title(p_in)
    fig.tight_layout()



def func_generator(x, r, p, s_params=["w0","a","b"]):
    w0, a, b, c = p["w0"], p["a"], p["b"], p["c"]
    w = w0*(1 + b*x**a)
    beta = c/w
    beta0 = c/w0
    A = (beta/beta0)**2
    f = A*sympy.exp(-(beta*r)**2)

    p = sympy.Matrix( [p[v] for v in s_params] )
    df_dp = f.diff(p)

    return f, df_dp


class JetProfile:

    """

    Attributes
    ----------
    a : float
        jet width spreading exponent
    b : float
        jet expansion ratio  (w0/wL)
    w0 : float
        width parameter at x = x0

    r0 : float
        nominal distance [m] in r-direction
    L  : float
        jet length [m]
    x0 : float
        jet starting location

    """

    def __init__(self):

        self.x0 = 0.0
        self.b  = 0.5
        self.a  = 0.5
        self.w0 = 1.0
        self.L  = 1.0
        self.r0 = 0.05


    @property
    def origin(self):
        return self.x0

    @origin.setter
    def origin(self, x0):
        self.x0 = x0

    def value(self, x_in, r_in, deriv=True):

        x0, r0, L = self.x0, self.r0, self.L
        w0, a, b =  self.w0, self.a, self.b
        #print("w0: ", w0, "a:", a, "b:", b)
        #print("x0: ", x0, "L:", L, "r0:", r0)
        x = (x_in - x0)/L
        r = r_in/r0

        w = w0*(1.0 + b*x**a)
        c = (-np.log(0.5))**0.5
        beta = c/w
        beta0 = c/w0
        #print("w0: ", w0, "a:", a, "b:", b)
        A = (beta/beta0)**2.0
        f = A*np.exp(-(beta*r)**2)

        if not deriv:
            return f, beta, w

        df_beta  = f*( 2.0/beta - 2.0*r**2*beta )
        df_beta0 = -f*2.0/beta0
        #print(df_beta.shape)
        #df_beta = -2*beta**3*r**2*np.exp(-beta**2*r**2)/beta0**2 + 2*beta*np.exp(-beta**2*r**2)/beta0**2
        #print(df_beta.shape)
        dbeta_dw = -beta/w
        #print(dbeta_dw.shape)
        dw_dp = np.array([ 1.0 + b*x**a, w0*b*x**a*np.log(x), w0*x**a])
        #print(dw_dp.shape)
        #dw_dp[0] *= r**2
        df_dp = dw_dp*dbeta_dw*df_beta

        dbeta0_dw0 = -beta0/w0
        df_dp[0] += df_beta0*dbeta0_dw0

        return f, df_dp, beta, w, dw_dp
    #@property
    #def params(self):
    #    return self.x0, self.r0, self.a, self.b, self.w0, self.L

    def params_names(self):
        return "w0 a b".split()


def test_jet_profile(p_in={"w0":1,"a":2.0,"b":6}):
    """plot function and sensitivity derivatitves using manual calculations
    and compare with numerial sensitivity derivatives"""

    jet = JetProfile()
    jet.b = p_in["b"]
    jet.a = p_in["a"]
    jet.w0 = p_in["w0"]
    jet.L  = 0.5
    jet.r0 = 0.008

    L, r0 = jet.L, jet.r0
    x = jet.L*np.linspace(0,1,63)
    r = r0*np.concatenate([np.linspace(0,1), np.logspace(0,1.0)])

    pos = np.meshgrid(x, r)

    x1 = pos[0].T
    r1 = pos[1].T

    fig, ax = plt.subplots(2,2)
    val, df, beta, w, dw  = jet.value(x1, r1)
    ax[0,0].contour(x1/L, r1/r0, val )
    ax[0,0].plot(x1[:,0]/L, w[:,0], 'k--', alpha=0.5, label=r"$w_{1/2}/r_0$")
    ax[0,0].set_xlabel("x/L")
    ax[0,0].set_ylabel("r/r0")
    ax[0,0].legend()
    nx, nr = len(x), len(r)

    cmap = plt.colormaps["viridis"]
    jdxs = np.arange(0,nr,20)
    for j in jdxs:
        color = cmap(float(j)/max(jdxs))
        label = None
        if j == 0 or j == jdxs[-1]:
            label="{:3.2f}".format(r1[0,j]/r0)
        ax[1,0].plot(x1[:,j], val[:,j]/val[:,j].max(), color=color, label=label)
    ax[1,0].set_xlabel("x/L")
    ax[1,0].set_ylabel("f(x,r)/f(r)_max")
    ax[1,0].legend(title="r/r0")

    idxs = np.arange(0,nx,8)
    for i in idxs:
        color = cmap(float(i)/idxs.max())
        label = None
        if i == 0 or i == idxs[-1]:
            label="{:3.2f}".format(x1[i,0]/L)
        ax[0,1].plot(val[i,:]/val[i,0],r1[i,:]/r0, color=color, label=label)


    ax[0,1].set_ylabel("r/r0")
    ax[0,1].set_xlabel("f/f_CL")
    ax[0,1].legend(title="x/L")

    ax[1,1].plot(x1[:,0],beta[:,0],'--',label="beta")
    ax[1,1].set_ylabel("beta")
    ax[1,1].legend(loc="center left")
    tw = ax[1,1].twinx()
    tw.plot(x1[:,0],w[:,0], label="w")
    tw.set_ylabel(r"$w_{1/2}/r_0$")
    tw.legend(loc="center right")
    fig.tight_layout()

    fig, axArray = plt.subplots(2,2,sharey=True)
    ax = axArray.flat
    i = 0
    f0 = val*1.0

    # perturbation for numerical derivative calculation
    eps = 1e-5

    for i_p, name in enumerate(jet.params_names()):
        print(i_p, name)
        #for i in range(0,nx,10):
        p0 = getattr(jet, name)*1.0
        setattr(jet, name, p0*(1.0 + eps))
        f1 = jet.value(x1, r1, False)[0]
        setattr(jet, name, p0*(1.0 - eps))
        f2 = jet.value(x1, r1, False)[0]
        df1 = (f1 - f2)/(2.0*p0*eps)
        setattr(jet, name, p0)
        for i in [10]:
            ax[i_p].plot( df[i_p,i,:], (r1/r0)[i,:], label="func" )
            ax[i_p].plot( df1[i,:], (r1/r0)[i,:], '--', label="num" )
            if i_p == 0:
                ax[-1].plot( f0[i,:], (r1/r0)[i,:], 'k--', label="base" )
            ax[-1].plot( f1[i,:], (r1/r0)[i,:], label=name )

        ax[i_p].set_xlabel("df_d"+name)
        ax[i_p].legend()
    ax[-1].legend()
    ax[-1].set_xlabel("f")
    for i in range(2):
        axArray[i,0].set_ylabel("r/r0")
    fig.tight_layout()


if __name__ == "__main__":
    test_sympy_diff()
    test_sympy_profile()
    test_jet_profile()


















