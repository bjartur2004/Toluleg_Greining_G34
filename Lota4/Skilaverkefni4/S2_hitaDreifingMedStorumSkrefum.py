import HitadreifingsPlot
import HitaJofnuhneppi_full

P = 5
L = 2 #cm
K = 1.68
H = 0.005
delta = 0.1 #cm
Power_L = 2 #cm
ambient_tempeature = 20 #°C

arguments = [P,L,K,H,delta,Power_L,ambient_tempeature]

n = 10
m = 10
u = HitaJofnuhneppi_full.solve_u(n,m,arguments)
HitadreifingsPlot.plotHitadreifingu(u,n,m, arguments, "Hitadreifing með n=m=10")