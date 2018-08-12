# Standard Model Neutron Decay Channels

import math

Mn =   939.5654133
sigMn =  0.0000058
Mp =   938.2720813
sigMp =  0.0000058
MK0 =  497.648
MK =   493.677
Mpi =  139.57018
sigMpi = 0.00035
Mpi0 = 134.9766
sigMpi0 =0.0006
Mmu =  105.658
Me =     0.5109989461
sigMe =  0.0000000013
Meta = 547.862
sigMeta =0.018
MetaP =957.78
sigMetaP=0.06 

fwhm=2*math.sqrt(2*math.log(2))
Ry=13.605693009#(84) eV
alpha=0.0072973525664#(17)	2014 CODATA

class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'


def onium(ma, mb, n):
	mu = ma*mb/(ma + mb)
	En = -0.5*mu*(alpha/n)**2
	return ma + mb + En


def antionium(ma, mb, n):
	mu = ma*mb/(ma + mb)
	En = 0.5*mu*(alpha/n)**2
	return ma + mb + En


def oniumspectrum(name, Ma, Mb, max):
	print(name + " (mass "+str(Ma+Mb)+") spectrum:")
	for n in range(1,max):
		print("    "+str(n)+"s: "+str((Ma+Mb-onium(Ma,Mb,n))*10**3) +"keV")


oniumspectrum("Ps", Me, Me, 4)
oniumspectrum("A(pi-mu+)", Mpi, Mmu, 8)
oniumspectrum("A(2pi)", Mpi, Mpi, 8)
oniumspectrum("A([pi0pi+]pi-)", Mpi+Mpi0, Mpi, 4)
oniumspectrum("A([pi0pi-]mu+)", Mpi+Mpi0, Mmu, 4)


particles = {
	"n": {
		'key':		"n",
		'name': 	"neutron",
		'latex': 	"n",
		'mass': 	Mn,
		'charge': 	0,
		'spin':		+0.5,
		'isospin':	-0.5,
		'parity':	+1,
		'baryon':	+1,
		'lepton':	0,
		'weak':		0,
		'width':	sigMn
	}, 
	"p": {
		'key':		"p",
		'name': 	"proton",
		'latex': 	"p",
		'mass': 	Mp,
		'charge': 	+1,
		'spin':		+0.5,
		'isospin':	+0.5,
		'parity':	+1,
		'baryon':	+1,
		'lepton':	0,
		'weak':		0,
		'width':	sigMp
	}, 
	'X(1576)': {
		'key':		"X(1576)",
		'name': 	"bqbq tetraquark",
		'latex': 	"X(",
		'mass': 	1576,
		'charge': 	0,
		'spin':		0,
		'parity':	1,
		'baryon':	0,
		'lepton':	0,
		'weak':		0,
		'width':	2
	},
	'f0(600)': {
		'key':		"f0(600)",
		'name': 	"Light tetraquark canidate",
		'latex': 	"eta'",
		'mass': 	500, #(2*Mpi + 0.7*alpha*(Mpi-Mpi0)),
		'charge': 	0,
		'spin':		0,
		'parity':	-1,
		'baryon':	0,
		'lepton':	0,
		'weak':		0,
		'width':	sigMpi+sigMpi0
	},
	'f(600)[L]': {
		'key':		"f(600)[L]",
		'name': 	"[ud][u*d*] light tetraquark",
		'latex': 	"eta'",
		'mass': 	2./3.*Meta + math.sqrt(2)/3.*(Meta+MetaP) + 1./3.*MetaP - 1./2.*(Mpi0 + 2*Mpi),
		'charge': 	0,
		'spin':		0,
		'parity':	-1,
		'baryon':	0,
		'lepton':	0,
		'weak':		0,
		'width':	sigMpi+sigMpi0
	},
#	"eta'": {
#		'key':		"eta'",
#		'name': 	"eta prime",
#		'latex': 	"\eta'",
#		'mass': 	957.78,
#		'charge': 	0,
#		'spin':		1,
#		'parity':	-1,
#		'baryon':	0,
#		'lepton':	0,
#		'weak':		0,
#		'width':	0.06
#	}, 
	"eta": {
		'key':		"eta",
		'name': 	"eta",
		'latex': 	"\eta",
		'mass': 	547.862,
		'charge': 	0,
		'baryon':	0,
		'lepton':	0,
		'weak':		0,
		'width':	0.018
	}, 
#	"K0": {
#		'key':		"K0",
#		'name': 	"kaon",
#		'latex': 	"K^0",
#		'mass': 	497.648,
#		'spin':		0,
#		'parity':	-1,
#		'charge': 	0,
#		'baryon':	0,
#		'lepton':	0,
#		'strange':  -1,
#		'weak':		0,
#		'width':	0.022
#	}, 
#	"K*": {
#		'key':		"K*",
#		'name': 	"antikaon",
#		'latex': 	"\overline{K}}^{0}",
#		'mass': 	497.648,
#		'spin':		0,
#		'parity':	-1,
#		'charge': 	0,
#		'baryon':	0,
#		'lepton':	0,
#		'strange':  1,
#		'weak':		0,
#		'width':	0.022

#	"K+": {
#		'key':		"K+",
#		'name': 	"charged kaon",
#		'latex': 	"K^{+}",
#		'mass': 	493.667,
#		'spin':		0,
#		'parity':	-1,
#		'charge': 	+1,
#		'baryon':	0,
#		'lepton':	0,
#		'strange':  -1,
#		'weak':		0,
#		'width':	0.013
#	},
#	"K-": {
#		'key':		"K-",
#		'name': 	"charged antikaon",
#		'latex': 	"K^{-}",
#		'mass': 	493.667,
#		'spin':		0,
#		'parity':	-1,
#		'charge': 	-1,
#		'baryon':	0,
#		'lepton':	0,
#		'strange':  +1,
#		'weak':		0,
#		'width':	0.013
#	}, 
	"A(pi-mu+)[1s]": {
		'key':		"A(pi-mu+)[1s]",
		'name': 	"antimuinated pionium",
		'latex': 	"A_{\pi^{-}\mu^+}^{1s}",
		'mass': 	onium(Mpi,Mmu,1),
		'charge': 	0,
		'spin':		0.5,
		'parity':	-1,
		'isospin':	0,
		'baryon':	0,
		'lepton':	-1,
		'weak':		0,
		'width':	2*0.00035
	},
	"A(pi+mu-)[1s]": {
		'key':		"A(pi+mu-)[1s]",
		'name': 	"muinated pionium",
		'latex': 	"A_{\pi^-\mu^-}^{1s}",
		'mass': 	onium(Mpi,Mmu,1),
		'charge': 	0,
		'spin':		0.5,
		'parity':	-1,
		'isospin':	0,
		'baryon':	0,
		'lepton':	+1,
		'weak':		0,
		'width':	2*0.00035
	},
	"A([pi-pi0]mu+)[1s]": {
		'key':		"A([pi-pi0]mu+)[1s]",
		'name': 	"antimuinated pionium",
		'latex': 	"A_{\[\pi^{-}\pi^0\]\mu^+}^{1s}",
		'mass': 	onium(Mpi+Mpi0,Mmu,1),
		'charge': 	0,
		'spin':		0.5,
		'parity':	-1,
		'isospin':	0,
		'baryon':	0,
		'lepton':	-1,
		'weak':		0,
		'width':	2*0.00035
	},
	"A([pi+pi0]mu-)[1s]": {
		'key':		"A([pi+pi0]mu-)[1s]",
		'name': 	"muinated pionium",
		'latex': 	"A_{\[\pi^-\pi^0\]\mu^-}^{1s}",
		'mass': 	onium(Mpi+Mpi0,Mmu,1),
		'charge': 	0,
		'spin':		0.5,
		'parity':	-1,
		'isospin':	0,
		'baryon':	0,
		'lepton':	+1,
		'weak':		0,
		'width':	2*0.00035
	},
	"A(2pi)[1s]": {
		'key':		"A(2pi)[1s]",
		'name': 	"1s pionium",
		'latex': 	"A_{2\pi}^{1s}",
		'mass': 	onium(Mpi,Mpi,1),
		'charge': 	0,
		'spin':		0,
		'parity':	0,
		'isospin':	0,
		'baryon':	0,
		'lepton':	0,
		'weak':		0,
		'width':	2*0.00035
	},
	"A(2pi)[2s]": {
		'key':		"A(2pi)[2s]",
		'name': 	"2s pionium",
		'latex': 	"A_{2\pi}^{2s}",
		'mass': 	onium(Mpi,Mpi,2),
		'charge': 	0,
		'spin':		0,
		'parity':	0,
		'isospin':	0,
		'baryon':	0,
		'lepton':	0,
		'weak':		0,
		'width':	2*0.00035
	},
	"A(2pi)[3s]": {
		'key':		"A(2pi)[3s]",
		'name': 	"3s pionium",
		'latex': 	"A_{2\pi}^{3s}",
		'mass': 	onium(Mpi,Mpi,3),
		'charge': 	0,
		'spin':		0,
		'parity':	0,
		'isospin':	0,
		'baryon':	0,
		'lepton':	0,
		'weak':		0,
		'width':	2*0.00035
	},
	"A(4pi)[1sig]": {
		'key':		"A(4pi)[1sig]",
		'name': 	"dipionium",
		'latex': 	"A_{2\pi}_2^{1\sigma}",
		'mass': 	onium(2*Mpi,2*Mpi,1),
		'charge': 	0,
		'spin':		0,
		'parity':	0,
		'isospin':	0,
		'baryon':	0,
		'lepton':	0,
		'weak':		0,
		'width':	4*0.00035
	},
	"A(4pi)[1sig*]": {
		'key':		"A(4pi)[1sig*]",
		'name': 	"antidipionium",
		'latex': 	"A_{2\pi}_2^{1\sigma^*}",
		'mass': 	antionium(2*Mpi,2*Mpi,1),
		'charge': 	0,
		'spin':		0,
		'parity':	0,
		'isospin':	0,
		'baryon':	0,
		'lepton':	0,
		'weak':		0,
		'width':	4*0.00035
	},
	"pi+": {
		'key':		"pi+",
		'name': 	"pion",
		'latex': 	"\pi^+",
		'mass': 	139.57018,
		'charge': 	+1,
		'spin':		0,
		'parity':	-1,
		'baryon':	0,
		'lepton':	0,
		'weak':		0,
		'width':	0.00035
	},
	"pi-": {
		'key':		"pi-",
		'name': 	"antipion",
		'latex': 	"\pi^{-}",
		'mass': 	139.57018,
		'charge': 	-1,
		'spin':		0,
		'parity':	-1,
		'baryon':	0,
		'lepton':	0,
		'weak':		0,
		'width':	0.00035
	}, 
	"pi0": {
		'key':		"pi0",
		'name': 	"pi zero",
		'latex': 	"\pi^0",
		'mass': 	134.9766,
		'charge': 	0,
		'spin':		0,
		'parity':	-1,
		'baryon':	0,
		'lepton':	0,
		'weak':		0,
		'width':	0.0006
	}, 
	"mu-": {
		'key':		"mu-",
		'name': 	"muon",
		'latex': 	"\mu^{-}",
		'mass': 	105.658,
		'charge': 	-1,
		'isospin':	-0.5,
		'baryon':	0,
		'lepton':	+1,
		'weak':		-1,
		'width':	0
	}, 
	"mu+": {
		'key':		"mu+",
		'name': 	"antimuon",
		'latex': 	"\mu^{+}",
		'mass': 	105.658,
		'charge': 	+1,
		'isospin':	+0.5,
		'baryon':	0,
		'lepton':	-1,
		'weak':		+1,
		'width':	0
	},
	"Ps[1s]": {
		'key':		"Ps[1s]",
		'name': 	"positronium",
		'latex': 	"\text{Ps}^{1s_0}",
		'mass': 	onium(Me,Me,1),
		'charge': 	0,
		'spin':		0,
		'parity':	0,
		'isospin':	0,
		'baryon':	0,
		'lepton':	0,
		'weak':		0,
		'width':	2*Me
	},
	"e-": {
		'key':		"e-",
		'name': 	"electron",
		'latex': 	"e^{-}",
		'mass': 	0.511,
		'charge': 	-1,
		'isospin':	-0.5,
		'baryon':	0,
		'lepton':	+1,
		'weak':		-1,
		'width':	0
	},
	"e+": {
		'key':		"e+",
		'name': 	"positron",
		'latex': 	"e^{+}",
		'mass': 	0.511,
		'charge': 	+1,
		'isospin':	+0.5,
		'baryon':	0,
		'lepton':	-1,
		'weak':		+1,
		'width':	0
	},
	"nu_e+": {
		'key':		"nu_e*",
		'name': 	"antielectron neutrino",
		'latex': 	"\bar{\nu}_e",
		'mass': 	0.000001,
		'charge': 	0,
		'isospin':	-0.5,
		'baryon':	0,
		'lepton':	-1,
		'weak':		+1,
		'width':	0
	},
	"nu_e": {
		'key':		"nu_e",
		'name': 	"electron neutrino",
		'latex': 	"\nu_e",
		'mass': 	0.000001,
		'charge': 	0,
		'isospin':	+0.5,
		'baryon':	0,
		'lepton':	+1,
		'weak':		-1,
		'width':	0
	} 
}


#def eVstr(value, error):
#	if error / 
#	return str("{0:.3f}({0:2f}) MeV ").format(value,error)


def findheaviest(particles, Mmin, Mmax):
	# Find the next heaviest decay channel
	heaviest = ''
	for key, x in particles.items():
		Mx = x["mass"]
		if Mx > Mmin and Mx < Mmax:
			Mmin = Mx
			heaviest = key
	Mx = Mmin
	return heaviest


def massorder(parts):
	keys=[]
	for key, X in parts.items():
		keys.append(key)

	#def bubble_sort(keys):
	changed = True
	while changed:
		changed = False
		for i in range(len(keys) - 1):
			Xa = parts[keys[i]]
			Xb = parts[keys[i+1]]
			if Xa["mass"] < Xb["mass"]:
				keys[i], keys[i+1] = keys[i+1], keys[i]
				changed = True
	return keys



def conserved( parent, coeffs, Emin ):
	M=dM=varM=Q=L=B=S=Yw=I=0
	X0 = particles[parent]
	M0 = X0["mass"] 
	Q0 = X0["charge"] 
	I0 = X0["isospin"] 
	Y0 = X0["weak"] 
	B0 = X0["baryon"]
	L0 = X0["lepton"]
	for key, n in coeffs.items():
		X = particles[key]
		M += n*X["mass"]
		Q += n*X["charge"]
		dM += n*X["width"] 
		varM += n*X["width"]**2
		#I += n*X["isospin"] 
		L += n*X["lepton"]
		B += n*X["baryon"]
		Yw += n*X["weak"]

	#print("M0="+str(M0)+" Q0="+str(Q0)+" B0="+str(B0)+" L0="+str(L0)+" B0-L0="+str(B0-L0))

	if Q == Q0 and ((B-L)-(B0-L0)) == 0 and M > Emin:
		rstr = "M="+str("{0:.3f}").format(M)
		#rstr += "M="+eVstr(M,math.sqrt(varM))
		rstr += str("({0:.0f}) MeV ").format(1000*math.sqrt(varM))
		rstr += " Q="+str(Q)+" B="+str(B)+" L="+str(L)+" B-L="+str(B-L)
		rstr += " Yw="+str(Yw)+" T3="+str(I)
		return rstr
	else:
		return ""


def reactionstr( parent, coeffs, E ):
	rstr = parent+" --> "
	for i in range(len(keys)):
		s = key=keys[i]
		n = coeffs[key]
		if n == 1: 
			rstr += s+" + "
		if n > 1: 
			rstr += str(n)+s+" + "
	if E > 0.100 :
		rstr += str("{0:.3f} MeV").format(E)
	else:
		rstr += bcolors.OKGREEN
		rstr += str("{0:.2f} keV").format(E*1000)
		rstr += bcolors.ENDC
	return rstr

def reactiontex( parent, coeffs, E ):
	rstr = parent+" --> "
	for i in range(len(keys)):
		key=keys[i]
		s = particles[key]["latex"]
		n = coeffs[key]
		if n == 1: 
			rstr += s+" + "
		if n > 1: 
			rstr += str(n)+" "+s+" + "
	rstr += str("{0:.3f} MeV").format(E)
	return rstr

def subbranch( parent, E, i, coeffs, Nmax, Emin ):
	if i==len(keys):
		cstr=conserved(parent, coeffs, Emin)
		if cstr:
			print("{0:90} {1:40}".format(reactionstr(parent, coeffs, E),"("+cstr+")"))
		return
	if Nmax < 1:
		#print("Recursion limit. Abandoning chain")
		return

	key = keys[i]
	X = particles[key]
	Mx = X["mass"]
	Xmax = Nmax
	if X["lepton"] != 0:
		if X["charge"] == 0:
			Xmax = Umax
		else:
			Xmax = Lmax
	Nx = min(int(math.floor(E/Mx)),Xmax)

	for n in range(Nx+1):
		coeffs[key] = n
		subbranch(parent, E-n*Mx, i+1, coeffs, Nmax-n, Emin)	


def newbranch( parent, Nmax, Emin):
	Emax = particles[parent]["mass"];
	subbranch( parent, Emax, 0, {}, Nmax, Emin)


keys=massorder(particles)
cols='{0:30} {1:20} {2:8} {3:3} {4:3} {5:3}'
print(cols.format("name", "key", "mass", "Q", "B", "L"))
for i in range(len(keys)):
	X = particles[keys[i]]
	print(cols.format(X["name"], X["key"], "{:4.3f}".format(X["mass"]), X["charge"], X["baryon"], X["lepton"]))


Nmax=12
Lmax=4
Umax=1
print("Decay branches for neutron...")
newbranch("n", Nmax, Mp)
newbranch("p", Nmax, Mp-1)




"""
print("Old way...")


# Old way
def pstr( n, particle ): 
	if n == 1: 
		return (particle+" + ")
	if n > 1: 
		return (str(n)+particle+" + ")
	else:
		return ""

N=0
def newbranches( Mmin, Mmax, Q0 ):
	for NK0 in range(3):
		N[K0] = NK0
		N = N[K0]
		for NKp in range(3-N):
			N[Kp] = NKp
			N += N[Kp]
			for NKm in range(3-N):
				N += NKm
				for Npip in range(8-N):
					N += Npip
					for Npim in range(8-N):
						N += Npim
						for Npi0 in range(8-N):
							N += Npi0
							for Nmup in range(9-N):
								N += Nmup
								for Nmum in range(9-N):
									N += Nmum
									for Ne in range(6):
										N += Ne
										for Nep in range(6-Ne):
											N += Nep
											M = (NK0*MK0 + (NKp + NKm)*MK 
												+ (Npip+Npim)*Mpi + Npi0*Mpi0 + (Nmup+Nmum)*Mmu + (Nep+Ne)*Me)
											Q = NKp - NKm + Npip - Npim + Nmup - Nmum + Nep - Ne
											L = -Nmup + Nmum - Nep + Ne
											pairs = abs(Nmup + Nmum)/2 + abs(Ne + Nep)/2 + Npi0
											if M <= Mmax and M >= Mmin and Q == Q0 and L <= 1 and L >= -2 and pairs < 4:
												print(
													pstr(NK0,"K0")+
													pstr(NKp,"K+")+
													pstr(NKm,"K-")+
													pstr(Npi0,"pi0")+
													pstr(Npip,"pi+")+
													pstr(Npim,"pi-")+
													pstr(Nmup,"mu+")+
													pstr(Nmum,"mu-")+
													pstr(Nep,"e+")+
													pstr(Ne,"e-")+
													str("{0:.3f} MeV    (B-L=").format(Mmax-M)+
													str(-L)+", "
													"log(a)="+str(pairs)+")")
											N = 0


def branches( Mmin, Mmax, Q0 ):
	for NK0 in range(3):
		N = NK0
		for NKp in range(3-N):
			N += NKp
			for NKm in range(3-N):
				N += NKm
				for Npip in range(8-N):
					N += Npip
					for Npim in range(8-N):
						N += Npim
						for Npi0 in range(8-N):
							N += Npi0
							for Nmup in range(9-N):
								N += Nmup
								for Nmum in range(9-N):
									N += Nmum
									for Ne in range(6):
										N += Ne
										for Nep in range(6-Ne):
											N += Nep
											M = (NK0*MK0 + (NKp + NKm)*MK 
												+ (Npip+Npim)*Mpi + Npi0*Mpi0 + (Nmup+Nmum)*Mmu + (Nep+Ne)*Me)
											Q = NKp - NKm + Npip - Npim + Nmup - Nmum + Nep - Ne
											L = -Nmup + Nmum - Nep + Ne
											pairs = abs(Nmup + Nmum)/2 + abs(Ne + Nep)/2 + Npi0
											if M <= Mmax and M >= Mmin and Q == Q0 and L <= 1 and L >= -2 and pairs < 4:
												print(
													pstr(NK0,"K0")+
													pstr(NKp,"K+")+
													pstr(NKm,"K-")+
													pstr(Npi0,"pi0")+
													pstr(Npip,"pi+")+
													pstr(Npim,"pi-")+
													pstr(Nmup,"mu+")+
													pstr(Nmum,"mu-")+
													pstr(Nep,"e+")+
													pstr(Ne,"e-")+
													str("{0:.3f} MeV    (B-L=").format(Mmax-M)+
													str(-L)+", "
													"log(a)="+str(pairs)+")")
											N = 0

print("n --> ")
branches( Mp, Mn, 0)
print("p --> ")
branches( Mp-10, Mp, 1)
"""
