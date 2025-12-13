import HitadreifingsPlot
import HitaJofnuhneppi_Sparce

P = 5
L = 2 #cm
K = 1.68
H = 0.005
delta = 0.1 #cm
Power_L = 2 #cm
ambient_tempeature = 20 #°C

arguments = [P,L,K,H,delta,Power_L,ambient_tempeature]

n = 100
m = 100
u = HitaJofnuhneppi_Sparce.solve_u(n,m,arguments)

print(f"u max: {max(u):.4f}")

HitadreifingsPlot.plotHitadreifingu(u,n,m, arguments, "Hitadreifing Leist með Rýrs fylkja aðferð")
