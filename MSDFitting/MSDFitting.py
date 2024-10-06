import numpy as np
import scipy.optimize as opt
import os
import MSDFuncs as msd
import time
import sys
import itertools
import random

class MSD_FIT(object):

    def __init__(self):

        #Set parameters for fitting procedure
        self.filename="JustMSD_InEKV15"                     #MSD file name with .txt extension ####################
        self.nmodes = 4                             #max number of modes to fit 
        self.tau = [2e-12, 2e-11, 1.75e-10, 3e-09, 7e-07]
        #self.nmodes = 2  
        #self.tau = [0.12, 9, 90]        #array of taus to use for fitting
        msd.diffusive_regime = True						#set to True if the MSD data reaches the diffusive regime (set to false for plateau)
        self.stdev_include = False                      #set to True if the MSD data file contains a third column with standard deviation 

        #Check if tau values are equal to nmodes+1
        self.check_taus()           

        #Make labels for alpha values
        self.make_labels()          

        #Load MSD data to fit
        self.read_msd()             
        
        #Calculate C1 coefficient
        msd.C1 = self.calc_C()

        #Initialize arrays and values for optimizer
        self.fit_alphas = np.ones(len(self.tau)-1,float)        #initialize array of alphas as 1.0's
        self.weights = np.zeros(len(self.tau),float)            #initialize array for weights (g_j)
        self.Iter = 0                                           #set iteration number
        self.cost_values = []                                   #list of chi2 values
        self.bnds = [(-9, 9) for j in range(0,self.nmodes)]     #boundaries for alpha values (if known)


    def check_taus(self):
        #Check if number of tau values = nmodes+1

        if (len(self.tau)-1) != self.nmodes:
            print("Not all tau values set.")
            sys.exit()


    def make_labels(self): 
        #Make labels for command line print
        
        self.labels = []
        for i in range(1,self.nmodes+1):
            self.labels.append("alpha%d" %i)
        self.labels.append("ChiSq")
        self.frmt = "{:9s}   "*len(self.labels)


    def find_nearest(self, array, value):
        #Find nearest value

        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx


    def read_msd(self): 
        #Read in MSD data to fit and shorten to 100 points and write to _short.txt file
        
        with open(self.filename + ".txt") as f:
            #f.readline() #Comment out if there is no header in file
            lines = f.readlines()
            self.tdata = np.array([float(line.split()[0]) for line in lines])*(1e-9) #converting ns to s
            self.MSDrawdata = np.array([float(line.split()[1]) for line in lines])*(1e-20) #converting angstrom^2 to m^2
            if self.stdev_include:
                self.stdev = np.array([float(line.split()[2]) for line in lines])*(1e-20) #converting angstrom^2 to m^2
            else:
                self.stdev = np.ones(len(self.tdata),float)

        cut = np.logspace(-13,-5.5,2500) #################### used -7 for try 1 and try 2 and it worked
        cut_tdata = []
        cut_MSDdata = []
        cut_stdev=[]
        skip = 1
        for j in range(0,len(cut)):
            idx = self.find_nearest(array=self.tdata,value=cut[j])
            
            if j == 0:
                cut_tdata.append(self.tdata[idx])
                cut_MSDdata.append(self.MSDrawdata[idx])
                cut_stdev.append(self.stdev[idx])

            
            elif (j > 0) and (self.tdata[idx] != cut_tdata[j-skip]):
                cut_tdata.append(self.tdata[idx])
                cut_MSDdata.append(self.MSDrawdata[idx])
                cut_stdev.append(self.stdev[idx])

            else:
                skip+=1
                continue

        self.tdata = cut_tdata
        self.MSDrawdata = cut_MSDdata
        self.stdev = cut_stdev

        nfile = open(self.filename + "_short.txt",'w')
        for i in range(0,len(self.tdata)):
            nfile.write("%16.8E %16.8E %16.8E\n"%(self.tdata[i],self.MSDrawdata[i],self.stdev[i]))
                

    def calc_C(self): 
        #Calculate coefficient C1 of MSD in ballistic region (averaging MSD/t^2 at t<0.01) 
        
        kmin = 2
        kmax = 5
        C = 0
        for k in range(kmin,kmax+1):
            C += self.MSDrawdata[k]/self.tdata[k]**2
        C1 = C/(kmax-kmin+1)
        #print "C = %0.8f" %(C1) 

        return C1


    def set_newparams(self, params):
        #Propogates minimizer's parameter shifts

        shifts = 0
        for i in range(len(params)):
            if params[i] != self.fit_alphas[i]: shifts += 1
            self.fit_alphas[i] = params[i]

        return shifts


    def get_status(self, shifts, chi2):
        #Monitor minimization procedure

        nparams = len(self.fit_alphas)

        if len(self.cost_values) == 0: self.cost_values.append(chi2)

        if shifts > 2 and chi2 < self.cost_values[-1]:
            t1 = time.time()
            print('Walltime=%0.2f(mins)\n'%((t1-self.t0)/60.0))
            self.cost_values.append(chi2)
            self.Iter += 1
            x = self.fit_alphas
            frmt = "{: 3.6f}   "*len(x)
            print('{: 4d}'.format(self.Iter), frmt.format(*x), '{: 3.6f}'.format(chi2),'\n')


    def chi2_func(self, params, xdata,ydata,stdev):
        #Calculate and store chi2 value for tracking minimization

        shifts = self.set_newparams(params)
        fit_params = [self.tau, params]
        
        if self.stdev_include:
            stdev = np.divide(stdev,ydata)

        logres = np.divide(msd.logresidual(xdata,ydata,fit_params),stdev[1:])
        current_chi2 = np.real(np.sum(np.square(np.divide(logres,stdev[1:]))))
        
        self.get_status(shifts, current_chi2)

        return logres


    def run_fit(self):
        #Begin fitting procedure

        print("\nMinimizing ChiSq...\n")
        print('{:4s}'.format('Iter'),' ', self.frmt.format(*self.labels), '\n')
        self.t0 = time.time()
        res = opt.leastsq(self.chi2_func,self.fit_alphas,args=(self.tdata,self.MSDrawdata,self.stdev))
        t1 = time.time()
        print("\n")
        print("Done. Walltime = %0.2f(mins)\n"%((t1-self.t0)/60.0))

        fit_params = [self.tau,res[0]]

        #Calculate chi squared from the fit
        if self.stdev_include:
            stdev = np.divide(self.stdev,self.MSDrawdata)
        else:
            stdev = self.stdev

        ChiSq = np.real(np.sum(np.square(np.divide(msd.logresidual(self.tdata,self.MSDrawdata,fit_params),stdev[1:]))))
        print("Chi squared = %.6f" %(ChiSq))

        #Write fitting parameters to file
        with open('MSDfit_params_%dmodes.dat'%(self.nmodes),'w') as f:
            f.write('%16.20f  \t \t    %16.20f\n' %(fit_params[0][0], msd.weights[0])) ####################
            for i in range(1,self.nmodes+1):
                f.write('%16.20f  %16.20f  %16.20f\n' %(fit_params[0][i], fit_params[1][i-1], msd.weights[i])) ####################
            f.write(' \n')
            f.write('ChiSq = %0.4f' %(ChiSq))


if __name__ == '__main__':

    fit_msd = MSD_FIT()
    fit_msd.run_fit()
