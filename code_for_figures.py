"""
Created on Mon Jul 31 11:20:22 2023

@author: sandy
"""
import numpy as np
import matplotlib.pyplot as plt

i_11 = 1

def m_11(R0):
    m11 = qc_squared(R0)/(R0-1)
    return m11

def m_20(R0):
    m_20 = -R0*qc_squared(R0)*i_11/(2*(R0-1)**2)
    return m_20

def m_22(R0):
    t1 = 4*qc_squared(R0)*(-R0+chi(R0)*(2*qc_squared(R0)+1))-R0
    t2 = m_11(R0)*i_11/(18*(R0-1))
    return t1*t2

def m_1111(R0,phi):
    num = 2*(1+np.cos(phi))*qc_squared(R0)*(-R0+chi(R0)*((1+np.cos(phi))*qc_squared(R0)+1))-R0
    den = (2*(1+np.cos(phi))-1)**2*(R0-1)
    return num/den*m_11(R0)*i_11

def m_111_minus_1(R0,phi):
    num = 2*(1-np.cos(phi))*qc_squared(R0)*(-R0+chi(R0)*((1-np.cos(phi))*qc_squared(R0)+1))-R0
    den = (2*(1-np.cos(phi))-1)**2*(R0-1)
    return num/den*m_11(R0)*i_11

def i_e(R0):
    ie = 1 - 1/R0
    return ie

def i_1111(R0,phi):
    num = R0+(1+np.cos(phi))*qc_squared(R0)*(2*R0+chi(R0)*(R0-1))
    den = (2*(1+np.cos(phi))-1)**2*(R0-1)
    return num/den*m_11(R0)*i_11

def i_111_minus_1(R0,phi):
    num = R0+(1-np.cos(phi))*qc_squared(R0)*(2*R0+chi(R0)*(R0-1))
    den = (2*(1-np.cos(phi))-1)**2*(R0-1)
    return num/den*m_11(R0)*i_11

def i_20(R0):
    return -m_20(R0)

def i_22(R0):
    t1 = R0+2*qc_squared(R0)*(2*R0+chi(R0)*(R0-1))
    t2 = m_11(R0)*i_11/(18*(R0-1))
    return t1*t2

def r1_2101(R0,phi):
    t1 = (chi(R0)*qc_squared(R0)*m_11(R0)/2)*((1+np.cos(phi))*i_1111(R0,phi)+(1-np.cos(phi))*i_111_minus_1(R0,phi))
    t2 = chi(R0)*qc_squared(R0)*i_11*(m_20(R0) + (np.cos(phi)/2)*(m_111_minus_1(R0,phi)-m_1111(R0,phi)))
    t3 = -r2_2101(R0,phi)
    return t1+t2+t3

def r2_2101(R0,phi):
    r2 = R0*(m_11(R0)*(i_20(R0)+(i_1111(R0,phi)+i_111_minus_1(R0,phi))/2)+(m_20(R0)+(m_1111(R0,phi)+m_111_minus_1(R0,phi))/2)*i_11)
    return r2

def r1_31(R0):
    r1 = chi(R0)*qc_squared(R0)*((m_20(R0)-m_22(R0)/2)*i_11+m_11(R0)*i_22(R0))-r2_31(R0)
    return r1
    
def r2_31(R0):
    r2 = R0*(m_11(R0)*(i_20(R0)+i_22(R0)/2)+(m_20(R0)+m_22(R0)/2)*i_11)
    return r2

def chi(R0):
    t1 = R0**2/(R0-1)
    t2 = (2*R0)/(np.sqrt(R0-1))
    return t1 + t2

def qc_squared(R0):
    return np.sqrt(R0-1)

def b1(R0,phi):
    b1 = (-1/(2*qc_squared(R0)+R0))*((R0-1)*r1_2101(R0,phi)+(qc_squared(R0)+R0)*r2_2101(R0,phi))
    return b1
    
    
def a1(R0):
    a1 = (-1/(2*qc_squared(R0)+R0))*((R0-1)*r1_31(R0)+(qc_squared(R0)+R0)*r2_31(R0))
    return a1

def eta(R0,phi):
    eta = b1(R0,phi)/a1(R0)
    return eta

def find_intersection(phis,etas,tol):
    """
    Defines the intersection point of graph based on a specified tolerance.

    Parameters
    ----------
    phis : list of floats
        List of angles for specified function parameters, varying the radius R0.
    etas : list of floats
        List of function outputs for specified function parameters with varying radius R0.
    tol : float
        Tolerance used for finding intersection points.

    Returns
    -------
    phis_at_1 : list of floats
        The angles at a specified radius.
    phis_at_minus_1 : list of floats
        The angles at a specified radius.

    """
    phis_at_1 = []
    phis_at_minus_1 = []
    for i in range(len(phis)):
        if etas[i] > 1-tol and etas[i] < 1+tol:
            phis_at_1.append(phis[i])
        if etas[i] > -1-tol and etas[i] < -1+tol:
            phis_at_minus_1.append(phis[i])
    return phis_at_1, phis_at_minus_1

def find_intersection_chi(R0_list,chis,tol):
    R0_at_10 = []
    for i in range(len(R0_list)):
        if chis[i] > 18-tol and chis[i] < 18+tol:
            R0_at_10.append(R0_list[i])
    return R0_at_10

def new_fun_for_fig(q_squared,R0):
    t1 = R0/(R0-1)
    t2 = R0 + q_squared + (R0-1)/q_squared
    return t1*t2
    
R0_list = np.array((1.762,2,6,10,14,18,22))
phis = np.linspace(0.01,np.pi/3-0.01,100000)
tol=10**(-3)
q_squared_list = np.linspace(0.01,10,100000)
"""
The following commented code was used to recreate Figure 2.2, 
built to test my solutions for all above functions.
"""
#R0_list = np.linspace(1.01,10,100000)
#R0_list = np.ones(100000)*10
# for i in range(len(R0_list)):
#     R0 = R0_list[i]
#     phis_at_1,phis_at_minus_1 = find_intersection(phis,eta(R0,phis),tol)
#     print(f'\n{R0}\n', phis_at_1,phis_at_minus_1)
#plt.plot(q_squared_list,new_fun_for_fig(q_squared_list,R0_list),c='blue')
#labels = ['$R_0$=1.762','$R_0$=2','$R_0$=6','$R_0$=10','$R_0$=14','$R_0$=18','$R_0$=22']
#plt.legend(labels)
# plt.xlabel(r'$q^2$')
# plt.ylabel(r'$\chi$')
# plt.xlim([0.0,10])
# plt.ylim([0,50])
# plt.scatter(3,18,c='black')
#plt.vlines(3,ymin=0,ymax=(17+(7/9)),colors='red',linestyles='--')
#plt.hlines(17+(7/9),xmin=0,xmax=3,colors='red',linestyles='--')
#plt.hlines(0,xmin=0,xmax=1,colors='blue')#,linestyles='--')
#plt.title("Testing Function Accuracy Using Figure 2.2",size="x-large")
# print(new_fun_for_fig(3,10))

"""
Plotting Figure 4.2 for paper.
"""
for i in range(len(R0_list)):
    R0 = R0_list[i]
    phis_at_1,phis_at_minus_1 = find_intersection(phis,eta(R0,phis),tol)
    print(f'\n{R0}\n', phis_at_1,phis_at_minus_1)
    plt.plot(phis,eta(R0,phis))
labels = ['$R_0$=1.762','$R_0$=2','$R_0$=6','$R_0$=10','$R_0$=14','$R_0$=18','$R_0$=22']
plt.legend(labels)
plt.xlabel(r'$\phi$')
plt.ylabel(r'$\eta$')
plt.xlim([0.5,0.7])
plt.ylim([-2,2])
plt.hlines(1,xmin=min(phis),xmax=max(phis),label='Eta=1',linestyles='--')
plt.hlines(-1,xmin=min(phis),xmax=max(phis),label='Eta=-1',linestyles='--')
plt.title("Intersection Points",size="x-large")
