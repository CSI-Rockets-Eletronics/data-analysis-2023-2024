import numpy as np

h = 25
k = 400
L_f = 20*1e-3
L = 5*1e-3
x = 2*1e-3
y = 32*1e-3
P = 2*(x+y)
A_c = x*y
#print(A_c)
A_f = P*L_f
#print(A_f) # good
N = 7
A_base = 26*32*1e-6
A_b = A_base - N*A_c
A_t = N*A_f+A_b
#print(A_t) # good
m = np.sqrt((h*P)/(k*A_c))
# print(m) # good
eta_f = ((np.sqrt(h*P*k*A_c))/(h*A_f))*np.tanh(m*L_f)
# ((np.sinh(m*L_f)+(h/(m*k))*np.cosh(m*L_f))/(np.cosh(m*L_f)+(h/(m*k))*np.sinh(m*L_f)))
#print(eta_f)
eta_0 = 1 - ((N*A_f)/A_t)*(1-eta_f)
print(eta_0)
R_t0 = 1/(eta_0*h*A_t)
print(R_t0)
R_tcpp = 3e-6
theta_b = 75-20
q_c = theta_b*((R_tcpp/A_base)+(L/(k*A_base)+R_t0))**(-1)
print(f'qc = {q_c} W')