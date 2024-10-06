import numpy as np
import matplotlib.pyplot as plt
import os
import MSDFuncs as msd


class MSD_ANALYSIS(object):

    def __init__(self):
        
        #######################################BE CAREFUL##############################################################
        # NEMD Data
        #self.NEMDFileName="NEMD_EKV15.txt" # NEMD G* data format: Omega G' +- G" +-
        ###############################################################################################################
        
        #Set parameters for fit and G* calc

        self.filename="JustMSD_InEKV15.txt"       #filename of MSD data   ####################
        self.nmodes = 4                       #number of modes
        
        msd.rho = 9.064479916377107429e+02                     #density of polymers    #########kg/m^3###########
        msd.kBT = ((0.0019872067*300)/(6.02214076*(10**23)))*4184    #kBT set from MD #########J=Kg.m^2.s^-2###########
        
        #N = 1000                            #polymer length
        #Npoly = 100                        #number of chains
        self.M_p = (736*100*(10**-3))/(6.02214076*(10**23))    #particle mass  #########Kg###########
        M_fluid =  (9648750*(10**-3))/(6.02214076*(10**23))   #total mass of fluid  #########Kg###########

        msd.M_frac = M_fluid/(M_fluid+self.M_p)  #Calculate fluid mass fraction

        msd.R = self.R = 10.6*(10**-10)  #radius of particle probe  #########m###########
        self.L = 225*(10**-10)    #box size  #########m###########
        self.RoverL = self.R/self.L     #R/L value
        
        #Time points for MSD (t) and G* (omega) calculation
        tPoints = 1000
        self.tMin = -13  ####################
        self.tMax = -5.5   ####################
        self.tArr=10**(self.tMin+(np.array(range(tPoints), float) + 1.0)/tPoints*(self.tMax-self.tMin))

        omegaPoints = 1000
        self.omegaMin = 5.5  ####################
        self.omegaMax = 13  ####################
        self.omegaArr=10**(self.omegaMin+(np.array(range(omegaPoints), float) + 1.0)/omegaPoints*(self.omegaMax-self.omegaMin))

        #Initialize parameter arrays
        self.tau = []
        self.alpha = []
        self.weights = []

        #Font settings for plots
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['xtick.major.size'] = 12
        plt.rcParams['xtick.major.width'] = 1.2
        plt.rcParams['xtick.minor.size'] = 6
        plt.rcParams['xtick.minor.width'] = 1.2
        plt.rcParams['ytick.major.size'] = 12
        plt.rcParams['ytick.major.width'] = 1.2
        plt.rcParams['ytick.minor.size'] = 6
        plt.rcParams['ytick.minor.width'] = 1.2
        plt.rcParams['xtick.direction'] = 'out'
        plt.rcParams['ytick.direction'] = 'out'

        font = {'family' : 'sans-serif',
                'weight' : 'bold',
                'size'   : 32}
        plt.rc('font', **font)

    def read_fitresults(self):
        #Open file containing fit parameters
        
        with open('MSDfit_params_%dmodes.dat'%(self.nmodes)) as f:
            line = f.readline()
            self.tau.append(float(line.split()[0]))
            self.weights.append(float(line.split()[1]))
            for i in range(0,self.nmodes):
                line = f.readline()
                self.tau.append(float(line.split()[0]))
                self.alpha.append(float(line.split()[1]))
                self.weights.append(float(line.split()[2]))

        #Combine parameters into one array
        self.params = [self.tau,self.alpha]

    def read_MSD(self):
        #Open file containing MSD simulation data
        
        with open(self.filename) as f:
            f.readline()
            lines = f.readlines()
            self.t = np.array([float(line.split()[0]) for line in lines])*(1e-9) #converting ns to s
            self.MSDrawdata = np.array([float(line.split()[1]) for line in lines])*(1e-20) #converting angstrom^2 to m^2 


    def calc_C1(self):
    #Calculate coefficient C1 of MSD in ballistic region (averaging MSD/t^2 at t<0.01)
    
        kmin = 2
        kmax = 5
        C = 0
        for k in range(kmin,kmax+1):
            C += self.MSDrawdata[k]/self.t[k]**2
        msd.C1 = C/(kmax-kmin+1)

    #######################################BE CAREFUL##############################################################
    def write_Gstar(self):
        with open("Gstar_results.txt",'w') as f:
            f.write("%s %14s %18s %18s %18s %18s %18s\n" %("time (omega)","GSER G'","GSER G''", "GSER with HI G'", "GSER with HI G''", "G' by IGSER","G'' by IGSER"))
            f.writelines("%.6E %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n" %(self.omegaArr[i],self.Gs0.real[i],self.Gs0.imag[i],self.Gs.real[i],self.Gs.imag[i],self.Gs_IGSER.real[i],self.Gs_IGSER.imag[i]) for i in range(0,len(self.omegaArr)))
    ##############################################################################################################
        
    ##############################################################################################################
    def write_fittedMSD(self):
        with open("fittedMSD_results.txt",'w') as f:
            f.write("%18s %18s\n" %("time","Fitted MSD"))
            f.writelines("%16.8f %16.8f\n" %(self.tArr[i]*(1e9),self.checkMSD[i]*(1e20)) for i in range(0,len(self.tArr)))
    ##############################################################################################################
    
    #######################################BE CAREFUL##############################################################
    '''
    def Read_NEMD_Data(self):
        with open(self.NEMDFileName) as NEMD_f:
            AllLines=NEMD_f.readlines()
            self.NEMDOmega=np.array([float(line.split()[0]) for line in AllLines])
            self.NEMD_Gp=np.array([float(line.split()[1]) for line in AllLines])
            self.NEMD_Gpp=np.array([float(line.split()[3]) for line in AllLines])
    '''
    ##############################################################################################################
    
    def plot_results(self):
        #Plot MSD fit, G*, and wave propagation results

        fig = plt.figure(figsize=(32, 32))
        fig1 = fig.add_subplot(221)

        fig1.scatter(self.t, self.MSDrawdata, c='b', s=120, label = "Data")
        fig1.plot(self.tArr, self.checkMSD, c='r', lw = 3, label="Fit")
        fig1.axvline(x=self.tau[0], c = 'black', lw = 3, linestyle = "--", alpha = 0.25)
        fig1.axvline(x=self.tau[1], c = 'black', lw = 3, linestyle = "--", alpha = 0.25)
        fig1.axvline(x=self.tau[2], c = 'black', lw = 3, linestyle = "--", alpha = 0.25)
        fig1.axvline(x=self.tau[3], c = 'black', lw = 3, linestyle = "--", alpha = 0.25)
        fig1.axvline(x=self.tau[4], c = 'black', lw = 3, linestyle = "--", alpha = 0.25)
        #fig1.axvline(x=self.tau[5], c = 'black', lw = 3, linestyle = "--", alpha = 0.25)
                
        fig1.set_title("MSD")
        fig1.set_xlabel("t (s)")
        fig1.set_ylabel(r'MSD ($m^2$)')
        fig1.set_xscale('log',)
        fig1.set_yscale('log')
        fig1.legend(frameon=False, prop = { "size": 30 })

        fig1.set_xlim(10**(self.tMin),10**self.tMax)
        fig1.set_ylim((10**-24),(10**-14)) #################### changed -15 to -13 from EKV15 files

        fig2 = fig.add_subplot(222)

        fig2.plot(self.omegaArr, self.Gs0.real,'--', c='r', lw = 3, label="GSER G'")
        fig2.plot(self.omegaArr, self.Gs0.imag, '-',c='r', lw = 3, label="GSER G''")
        fig2.plot(self.omegaArr, self.Gs.real, '--', c='b', lw = 3, label="GSER w/ HI G'")
        fig2.plot(self.omegaArr, self.Gs.imag, '-', c='b', lw = 3, label="GSER w/ HI G''")
        
        #######################################BE CAREFUL##############################################################
        # IGSER added
        fig2.plot(self.omegaArr, self.Gs_IGSER.real,'--', c='k', lw = 3, label="IGSER G'")
        fig2.plot(self.omegaArr, self.Gs_IGSER.imag, '-', c='k', lw = 3, label="IGSER G''")
        # NEMD added
        #fig2.plot(self.NEMDOmega, self.NEMD_Gp,'o', c='k', markersize=16, label="NEMD G'")
        #fig2.plot(self.NEMDOmega, self.NEMD_Gpp,'o', c='k', markersize=16, label="NEMD G''")
        ##############################################################################################################
        
        fig2.set_title("G*")
        fig2.set_xlabel(r'$\omega$ (1/s)')
        fig2.set_ylabel(r'$G^*$ (Pa)')
        fig2.set_xscale('log')
        fig2.set_yscale('log')

        fig2.set_xlim(10**(self.omegaMin),10**self.omegaMax)
        fig2.set_ylim((10**3),(5*10**8)) ####################

        fig2.legend(frameon=False, prop = { "size": 25 })

        fig3 = fig.add_subplot(223)

        fig3.plot(self.omegaArr, self.Delta*(10**10), c='cyan', lw = 3, label=r'$\Delta(\omega)$')
        fig3.plot(self.omegaArr, self.Lambda*(10**10), c='magenta', lw = 3, label=r'$\Lambda(\omega)$')
        fig3.axhline(y=self.L*(10**10), c = 'black', lw = 3, linestyle = "--")
        
        fig3.set_title("Penetration depth & Wave Propagation")
        fig3.set_xlabel(r'$\omega$ (1/s)')
        fig3.set_ylabel(r'$\Delta(\omega), \Lambda(\omega)\,(\AA)$')
        fig3.set_xscale('log')
        fig3.set_yscale('log')

        fig3.set_xlim(10**(self.omegaMin),10**self.omegaMax)
        fig3.set_ylim((10),(5*(10**3))) ####################

        fig3.legend(frameon=False, prop = { "size": 30 })

        fig4 = fig.add_subplot(224)

        fig4.plot(self.omegaArr, self.L/self.Lambda, lw = 3, c='green')
        fig4.axhline(y=1.0, c = 'black', lw = 3, linestyle = "--")

        fig4.set_title("Box length x Wave Propagation")
        fig4.set_ylabel("L"r'$\Lambda^{-1}(\omega)$')
        fig4.set_xlabel(r'$\omega$ (1/s)')
        fig4.set_xscale('log')
        fig4.set_yscale('log')

        fig4.set_xlim(10**(self.omegaMin),10**self.omegaMax)
        fig4.set_ylim(10**(-4),10**4)

        fig.subplots_adjust(hspace=0.4,wspace=0.4)

        plt.show()
        fig.savefig("InEKV15_N250_Loc1.png", dpi = 150)
    def run_analysis(self):
        #Calculate MSD fit from parameters, G*, wave penetration depth, wavelength, and make plots

        #Read in fit results and MSD data
        self.read_fitresults()
        self.read_MSD()

        #Calculate effective particle mass
        self.calc_C1()
        msd.meff0 = self.M_p + ((2/3)*np.pi*((msd.R)**3)*msd.rho) #3.0*msd.kBT/msd.C1
        print("Effective particle mass = %.30f" %(msd.meff0)) ####################
        print("Bare particle mass = %.30f"%(self.M_p)) ####################

        #Calculate MSD fit line
        msd.nmodes = self.nmodes
        self.checkMSD = msd.MSD_t_vec(time=self.tArr,params=self.params,weights=self.weights)

        #Calculate MSD in Fourier space
        r_p = msd.r_p_vec(time=self.omegaArr,params=self.params,weights=self.weights)
        r_dp = msd.r_dp_vec(time=self.omegaArr,params=self.params,weights=self.weights)

        #Calculate Z for IGSER calculation
        Z = msd.calcZ_vec(time=self.omegaArr,rprime=r_p,rdprime=r_dp)
        self.Gs_IGSER = msd.Gstar_IGSER_vec(time=self.omegaArr,z=Z)

        #Calculate GSER
        self.Gs0 = msd.Gstar_GSER_vec(time=self.omegaArr,msd_omega=(r_p+1j*r_dp))

        #Calculate GSER with HI
        self.Gs = msd.Gstar_hydro_vec(time=self.omegaArr,msd_omega=(r_p+1j*r_dp),RoverL=self.RoverL)

        #Write results to file
        self.write_Gstar()

        #Calculate wave propagation properties, returns [penetration depth, wavelength]
        [self.Delta, self.Lambda] = msd.wavefunc(time=self.omegaArr,Gstar=self.Gs)

        self.write_fittedMSD()
        
        #######################################BE CAREFUL##############################################################
        # Read NEMD Data
        #self.Read_NEMD_Data()
        ###############################################################################################################
        
        self.plot_results()


if __name__ == '__main__':

    analysis = MSD_ANALYSIS()
    analysis.run_analysis()


