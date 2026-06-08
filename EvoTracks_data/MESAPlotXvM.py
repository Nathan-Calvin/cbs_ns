"""
MESAPlotXvM.py

MESA summary files created by MESAZAMS.py, which reads digests of MESA profiles made by MESAoutput3.py.
Summaries of current interest combine models of varying masses but with the same core value of X.

Read in one or more MESA summary files, noting number of entries and range of masses.
Make the following plots for each summary table:
1) log(time) vs. log(M) for nuclear, age, thermal, and dynamic timescales.
2) log(time) vs. log(M) for nuclear, thermal, and dynamic timescales.
3) zeta_thermal vs. log(M), inferred from adjacent rows of data. (Also make csv of this.)
4) Delta(log10(M)) vs. log(M), inferred from adjacent rows of data.
5) Convection zones vs. log(M), scaled to fractional radius
6) log(R) vs. log(M)

For a pair of MESA summary files, plot the following comparisons:
1) log(time) vs. log(M) for nuclear, age, thermal, and dynamic timescales. (Dashed lines for second file.)
2) zeta_thermal vs. log(M), inferred from adjacent rows of data. (Different color for second file.)
3) Ratio R(second file)/R(first file) vs. log(M).

Written 2/27/2023 by L. Molnar
Revised 5/7/2023 by L. Molnar to add log(R) vs. log(M) plot
Revised 5/7/2023 by L. Molnar to optionally add binaryMESA values to zeta plot
"""
from readHeader import readHeader
from wc import readData
import matplotlib.pyplot as plt
from numpy import array, exp, zeros, log, log10, interp
from math import pi, sqrt, ceil, cos
from scipy.optimize import fsolve, minimize, minimize_scalar
import numpy as np

def MakeMRKTable(MRKinFile,case):
    """Make a table of stellar parameters from a MESA summary file.
    INPUTS:
        MRKinFile: directory of MRK file
        case: string name for summary file
    OUTPUTS:
        MRKtable: lists of 0) M, 1) R, 2) k2, 3) kAM, 4) logT, 5) logL, 6) t_age, 7) t_dyn, 
        8) t_KH, 9) t_nuc, 10) logM, 11) Nconv, 12) Rcz,min, 13) Rcz,max, 14) Rcz,min, 15) Rcz,max
    """
    if MRKinFile == None:
        raise Exception('No MRK file found.') 
    else:
        MRKFile = open(MRKinFile)
        (header, HasEOH) = readHeader(MRKFile)
        data = readData(MRKFile,',')
        MRKtable = [list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),\
                    list(),list(),list(),list(),list()]
        for i in range(len(data)):
            MRKtable[0].append(float(data[i][4])) # [Msun] mass
            MRKtable[1].append(float(data[i][5])) # [Rsun] radius
            MRKtable[2].append(float(data[i][9])) # [k2_ave] average moment of inertia
            MRKtable[3].append(float(data[i][12])) # [kAM] apsidal motion coefficient
            MRKtable[4].append(log10(float(data[i][6]))) # log10(temperature)           
            MRKtable[5].append(log10(float(data[i][7]))) # log10(luminosity)
            MRKtable[6].append(float(data[i][0])) # t_age (year)
            MRKtable[7].append(float(data[i][1])) # t_dyn (year)
            MRKtable[8].append(float(data[i][2])) # t_KH (year)
            MRKtable[9].append(float(data[i][3])) # t_nuc (year)
            MRKtable[10].append(log10(float(data[i][4]))) # log10(M/Msun)
            MRKtable[11].append(int(data[i][17])) # Number of convection zones
            if MRKtable[11][-1] == 0:
                MRKtable[12].append(0.)
                MRKtable[13].append(0.)
                MRKtable[14].append(0.)
                MRKtable[15].append(0.)
            elif MRKtable[11][-1] == 1:
                MRKtable[12].append(float(data[i][18]))
                MRKtable[13].append(float(data[i][19]))
                MRKtable[14].append(0.)
                MRKtable[15].append(0.)
            elif MRKtable[11][-1] == 2:
                MRKtable[12].append(float(data[i][18]))
                MRKtable[13].append(float(data[i][19]))
                MRKtable[14].append(float(data[i][20]))
                MRKtable[15].append(float(data[i][21]))
            else:
                print('Error on number of convection zones:N,mass:',MRKtable[11][-1],MRKtable[0][-1])
    print('MESA case {:s}, {:d} Rows, mass range: ({:.4f},{:.4f})'.\
          format(case,len(MRKtable[0]),min(MRKtable[0]),max(MRKtable[0])))
    return MRKtable

def FilterMRKTable(Tin,DMmin):
    """Remove excess data rows in MRKtable
    INPUTS:
    Tin: raw MRKtable (column 10 assumed to be log10(M)
    DMmin: minimum allowed value of Delta(log10(M))
    OUTPUTS:
    Tout: filtered MRKtable
    """
    # Generate empty list of lists
    Tout = list()
    for i in range(len(Tin)):
        Tout.append(list())
    # Populate with selected mass values
    for i in range(len(Tin[0])):
        Add = False
        if i == 0: # Always add the first line
            Add = True
        else:
            Delta = Tin[10][i] - Tout[10][-1]
            if Delta >= DMmin: # Only include interior lines with proper spacing
                Add = True
            else:
                if i == (len(Tin[10])-1): # Always add the last line
                    Add = True
        if Add:
            for j in range(len(Tin)):
                Tout[j].append(Tin[j][i])
    print('Rows in filtered table:',len(Tout[0]))
    return Tout

def PlotMRKtable(MRKtable,title,prefix,outFolder,output,BM):
    """Plot key properties of one MRK table"""
    # Plot three timescales plus age vs. mass
    plotAgeM(MRKtable[0],MRKtable[6],MRKtable[7],MRKtable[8],MRKtable[9],title,prefix,outFolder,output,True)
    # Plot three timescales only vs. mass
    plotAgeM(MRKtable[0],MRKtable[6],MRKtable[7],MRKtable[8],MRKtable[9],title,prefix,outFolder,output,False)
    # Plot zeta_thermal vs. mass, also save as a csv file
    plotAlphaM(MRKtable[0],MRKtable[1],title,prefix,outFolder,output,BM)
    csvAlphaM(MRKtable[0],MRKtable[1],title,prefix,outFolder,output)
    # Plot Delta(log10(M)) of table vs. mass
    plotDlmM(MRKtable[0],MRKtable[10],title,prefix,outFolder,output)
    # Plot convection zones of table vs. mass
    plotConvection(MRKtable,title,prefix,outFolder,output)
    # Plot log(R) vs log(M)
    plotRM(MRKtable[0],MRKtable[1],title,prefix,outFolder,output)
    return

def PlotMRKtables(MRKtable,Cases,title,prefix,outFolder,output):
    """Plot comparisons of MRK tables"""
    # Plot three timescales plus age vs. mass
    plotAgeMCompare(MRKtable,title,prefix,outFolder,output)
    # Plot zeta_thermal vs. mass, also save as a csv file
    plotAlphaMCompare(MRKtable,Cases,prefix,outFolder,output)
    # Plot R_case1/R_case0 vs. mass
    plotRCompare(MRKtable,title,prefix,outFolder,output)
    return

def plotAgeM(mass,tage,tdyn,tkh,tnuc,title,prefix,outFolder,output,Age):
    """ Plot three timescales (optionally adding age) as a function of log10(M).
    INPUTS:
    Age: boolean to also plot age
    """
    plt.clf()
    ax1 = plt.gca()
    #
    if Age:
        Label = [ 'nuclear','age','thermal','dynamic' ]
    else:
        Label = [ 'nuclear','thermal','dynamic' ]
    plt.axis([0.1,max(mass),1.E3,1.E12])
    plt.loglog(mass,tnuc,'r-',label='nuclear')
    if Age: plt.loglog(mass,tage,'-',color='orange',label='age')
    plt.loglog(mass,tkh,'g-',label='thermal')
    plt.loglog(mass,tdyn,'b-',label='dynamic')
    #
    plt.xlabel(r'Mass ($M_{\odot}$)', fontweight='bold', fontsize='large')
    plt.ylabel('Time (y)', fontweight='bold', fontsize='large')
    leg = plt.legend(Label,loc=(.05,0.2), shadow=True, numpoints=1,prop={'size':15})
    plt.title(title)
    if output == 'screen':                      # Final graph
        plt.show()
    elif output == 'png' or 'eps':
        if Age:
            fname = outFolder+prefix+'AgeM.'+output
        else:
            fname = outFolder+prefix+'TimeM.'+output
        plt.savefig(fname,dpi=180,format=output)
        print('Save AgeM plot.')
    else:
        print('Unknown plot option', output)
    return

def plotAgeMCompare(MRKtable,title,prefix,outFolder,output):
    """ Plot three timescales and age as a function of log10(M) for two values of X.
    INPUTS:
    """
    plt.clf()
    ax1 = plt.gca()
    #
    Label = [ 'nuclear','age','thermal','dynamic' ]
    plt.axis([0.1,max(MRKtable[0][0]),1.E3,1.E12])
    # Plot first Xcore case, include labels
    plt.loglog(MRKtable[0][0],MRKtable[0][9],'r-',label='nuclear')
    plt.loglog(MRKtable[0][0],MRKtable[0][6],'-',color='orange',label='age')
    plt.loglog(MRKtable[0][0],MRKtable[0][8],'g-',label='thermal')
    plt.loglog(MRKtable[0][0],MRKtable[0][7],'b-',label='dynamic')
    leg = plt.legend(Label,loc=(.05,0.2), shadow=True, numpoints=1,prop={'size':15})

    # Plot second Xcore case, use dash lines
    plt.loglog(MRKtable[0][0],MRKtable[1][9],'r--')
    plt.loglog(MRKtable[0][0],MRKtable[1][6],'--',color='orange')
    plt.loglog(MRKtable[0][0],MRKtable[1][8],'g--')
    plt.loglog(MRKtable[0][0],MRKtable[1][7],'b--')
    
    plt.xlabel(r'Mass ($M_{\odot}$)', fontweight='bold', fontsize='large')
    plt.ylabel('Time (y)', fontweight='bold', fontsize='large')

    plt.title(title)
    if output == 'screen':                      # Final graph
        plt.show()
    elif output == 'png' or 'eps':
        fname = outFolder+prefix+'TimeM.'+output
        plt.savefig(fname,dpi=180,format=output)
        print('Save Time comparison plot.')
    else:
        print('Unknown plot option', output)
    return

def plotAlphaM(mass,radius,title,prefix,outFolder,output,BM):
    """Plot zeta adiabated vs. log10(M)
    INPUTS:
    BM = True implies also plot binaryMESA result"""
    # Binary MESA values
    MESAm = [0.6,0.7,0.8,0.9]
    MESAq = [0.900,0.843,0.738,0.689]
    MESAz = [0.697,0.859,1.222,1.429]
    #
    plt.clf()
    ax1 = plt.gca()
    # Compute alpha
    massave = list()
    alpha = list()
    for i in range(len(mass)-1):
        massave.append(sqrt(mass[i]*mass[i+1]))
        alpha.append( log(radius[i+1]/radius[i])/log(mass[i+1]/mass[i]) )
    #
    plt.axis([0.1,max(mass),0.,1.5])
    plt.semilogx(massave,alpha,'r-',label='Thermal')
    if BM: plt.semilogx(MESAm,MESAz,'b-',label='Adiabatic')
    #
    plt.xlabel(r'Mass ($M_{\odot}$)', fontweight='bold', fontsize='large')
    plt.ylabel(r'$\zeta\equiv \frac{d{\rm ln}(R)}{d{\rm ln}(M)}$', fontweight='bold', fontsize='large')
    plt.title(title)
    leg = plt.legend(loc=(.05,0.1), shadow=True, numpoints=1,prop={'size':15})
    if output == 'screen':                      # Final graph
        plt.show()
    elif output == 'png' or 'eps':
        fname = outFolder+prefix+'AlphaM.'+output
        plt.savefig(fname,dpi=180,format=output)
        print('Save AlphaM plot.')
    else:
        print('Unknown plot option', output)
    return

def plotAlphaMCompare(MRKtable,Cases,prefix,outFolder,output):
    """Plot zeta adiabated vs. log10(M)"""
    plt.clf()
    ax1 = plt.gca()
    # Compute average masses
    massave = list()
    for i in range(len(MRKtable[0][0])-1):
        massave.append(sqrt(MRKtable[0][0][i]*MRKtable[0][0][i+1]))
    # Compute zetas
    alpha = [list(),list()]
    for i in range(len(MRKtable[0][0])-1):
        for j in range(2):
            alpha[j].append( log(MRKtable[j][1][i+1]/MRKtable[j][1][i])/\
                             log(MRKtable[0][0][i+1]/MRKtable[0][0][i]) )
    #
    plt.axis([0.1,max(MRKtable[0][0]),0.,1.5])
    plt.semilogx(massave,alpha[0],'r-',label=Cases[0])
    plt.semilogx(massave,alpha[1],'b-',label=Cases[1])
    leg = plt.legend(loc=(.05,0.1), shadow=True, numpoints=1,prop={'size':15})
    #
    plt.xlabel(r'Mass ($M_{\odot}$)', fontweight='bold', fontsize='large')
    plt.ylabel(r'$\zeta_{\rm thermal}\equiv \frac{d{\rm ln}(R)}{d{\rm ln}(M)}$', fontweight='bold', fontsize='large')
    if output == 'screen':                      # Final graph
        plt.show()
    elif output == 'png' or 'eps':
        fname = outFolder+prefix+'zeta.'+output
        plt.savefig(fname,dpi=180,format=output)
        print('Save zetaCompare plot.')
    else:
        print('Unknown plot option', output)
    return

def csvAlphaM(mass,radius,title,prefix,outFolder,output):
    """Make csv table of  power law index dln(R)/dln(M) vs. log10(M)"""
    if output != 'screen':
        # Compute alpha
        massave = list()
        alpha = list()
        for i in range(len(mass)-1):
            massave.append(sqrt(mass[i]*mass[i+1]))
            alpha.append( log(radius[i+1]/radius[i])/log(mass[i+1]/mass[i]) )
        # Save to csv
        fname = outFolder+prefix+'AlphaM.csv'
        report = open(fname,'w')
        report.write('MESAPlotXvM.py\n'+title+'\n')
        report.write('M/Msun,alpha\nEOH\n')
        for i in range(len(massave)):
            report.write('{:.5f},{:.5f}\n'.format(massave[i],alpha[i]))
        report.close()
        print('Save AlphaM csv.')
    return

def plotDlmM(mass,lm,title,prefix,outFolder,output):
    """Plot Delta(log10(M)) vs. log10(M)"""
    plt.clf()
    ax1 = plt.gca()
    # Compute Delta(log10(M))
    massave = list()
    dlm = list()
    for i in range(len(mass)-1):
        massave.append(sqrt(mass[i]*mass[i+1]))
        dlm.append(lm[i+1]-lm[i])
    #
    plt.axis([0.1,max(mass),0.1,max(dlm)])
    print(mass)
    print(massave)
    print(dlm)
    plt.semilogx(massave,dlm,'r-')
    #
    plt.xlabel(r'Mass ($M_{\odot}$)', fontweight='bold', fontsize='large')
    plt.ylabel(r'$\Delta\log_{10}(M)$', fontweight='bold', fontsize='large')
    plt.title(title)
    if output == 'screen':                      # Final graph
        plt.show()
    elif output == 'png' or 'eps':
        fname = outFolder+prefix+'DlmM.'+output
        plt.savefig(fname,dpi=180,format=output)
        print('Save DeltaLogM plot.')
    else:
        print('Unknown plot option', output)
    return

def plotConvection(MRKtable,title,prefix,outFolder,output):
    """Plot convection zones for one MRKtable vs. mass
    INPUTS:
    MRKtable: 0) Mass, 11) Nconv, 12) Rcz1,min, 13) Rcz1,max, 14) Rcz2,min, 15) Rcz2,max
    """
    # First, generate lists with extent of zones
    mass = MRKtable[0]
    nconv = MRKtable[11]
    mc = list() # core convective zone
    rcmin = list()
    rcmax = list()
    cindex = list() 
    me = list() # envelope convective zone
    remin = list()
    remax = list()
    for i in range(len(mass)):
        for j in range(nconv[i]):
            rmin = MRKtable[12+j*2][i]
            rmax = MRKtable[13+j*2][i]
            if rmin+rmax<0.9: # A core convective zone
                if len(mc) != 0:
                    if i != cindex[-1]+1: # If there is a gap, put zeroes at boundaries
                        mc.append(mass[cindex[-1]+1])
                        rcmin.append(0.)
                        rcmax.append(0.)
                        cindex.append(cindex[-1]+1)
                        mc.append(mass[i-1])
                        rcmin.append(0.)
                        rcmax.append(0.)
                        cindex.append(i-1)
                else: # For first value, add in a zero on the left
                    mc.append(mass[i-1])
                    rcmin.append(0.)
                    rcmax.append(0.)
                    cindex.append(i-1)
                mc.append(mass[i])
                rcmin.append(rmin)
                rcmax.append(rmax)
                cindex.append(i)
            else: # An envelope convective zone
                me.append(mass[i])
                remin.append(rmin)
                remax.append(rmax)
    # Second, plot the zones
    plt.clf()
    ax1 = plt.gca()
    #
    plt.axis([0.1,max(mass),0.,1.])
    plt.fill_between(mc,rcmin,rcmax,color='lightgray')
    plt.fill_between(me,remin,remax,color='lightgray')
    plt.xscale('log')
    plt.xlabel(r'Mass ($M_{\odot}$)', fontweight='bold', fontsize='large')
    plt.ylabel('Fractional Radius', fontweight='bold', fontsize='large')
    plt.title('Convective zones, '+title)
    if output == 'screen':                      # Final graph
        plt.show()
    elif output == 'png' or 'eps':
        fname = outFolder+prefix+'Convection.'+output
        plt.savefig(fname,dpi=180,format=output)
        print('Save convection zone plot.')
    else:
        print('Unknown plot option', output)
    return

def plotRM(mass,radius,title,prefix,outFolder,output):
    """Plot log10(R) vs. log10(M)"""
    plt.clf()
    ax1 = plt.gca()
    #
    plt.axis([0.1,max(mass),0.1,max(radius)])
    plt.loglog(mass,radius,'r-')
    #
    plt.xlabel(r'Mass ($M_{\odot}$)', fontweight='bold', fontsize='large')
    plt.ylabel(r'Radius ($R_\odot$)', fontweight='bold', fontsize='large')
    plt.title(title)
    if output == 'screen':                      # Final graph
        plt.show()
    elif output == 'png' or 'eps':
        fname = outFolder+prefix+'RM.'+output
        plt.savefig(fname,dpi=180,format=output)
        print('Save RM plot.')
    else:
        print('Unknown plot option', output)
    return

def plotRCompare(MRKtable,title,prefix,outFolder,output):
    """Plot R_Case[1]/R_Case[0] vs. M"""
    plt.clf()
    ax1 = plt.gca()
    # Compute ratios
    r = list()
    for i in range(len(MRKtable[0][1])):
        r.append(MRKtable[1][1][i]/MRKtable[0][1][i])

    
    plt.axis([0.1,max(MRKtable[0][0]),min(r),max(r)])
    plt.semilogx(MRKtable[0][0],r,'r-')
    #
    plt.xlabel(r'Mass ($M_{\odot}$)', fontweight='bold', fontsize='large')
    plt.ylabel(r'$R_1/R_0$', fontweight='bold', fontsize='large')
    plt.title(title)
    if output == 'screen':                      # Final graph
        plt.show()
    elif output == 'png' or 'eps':
        fname = outFolder+prefix+'Radius.'+output
        plt.savefig(fname,dpi=180,format=output)
        print('Save Rcompare plot.')
    else:
        print('Unknown plot option', output)
    return

''' Plots below this point '''

def printoutput(outFolder,outFile,output):
    ''' Generates output for plots
    INPUT:
        outFolder, outFile, output: output information
    '''
    if output == 'screen':
        plt.show()
    elif output == 'png' or 'eps':
        plt.savefig(outFolder+outFile+'.'+output,dpi=300,format=output)
        print('Write out',output)
    else:
        print('Unknown plot option', output)

if __name__ == '__main__':

    #Run "MESAPlotXvM.py"
    # User input parameters

    Cases = [ '0.6', '0.697' ] # Choose '0.6', '0.69696'
    outFolder = '/storage/wumas/nls/EvoTracks_data/plots'
    output = 'png' # 'png', 'screen', 'eps'

    # Read in MESA summary files, save as MRK tables
    MRKtable = list()
    for i in range(len(Cases)):
        MESA = '/storage/wumas/nls/EvoTracks_data/Summary/Sum-X='+Cases[i]+'.csv'
        MRKtable.append(MakeMRKTable(MESA,Cases[i]))

    # Explore range of log10(M) evaluated
    old = -0.82
    print('  logA,   logB,  A-B, A-A(i-1)')
    for i in range(len(MRKtable[0][0])):
        a = log10(MRKtable[0][0][i])
        b = log10(MRKtable[1][0][i])
        c = a - b
        d = a - old
        print('{:.4f},{:.4f},{:.4f},{:.4f}'.format(a,b,c,d))
        old = a
        
    # Make a variety of plots based on them.

    # First, make plots based on a single value of X
    for i in range(len(Cases)):
        title = r'$X_{\rm core} = $'+Cases[i]
        prefix = 'X'+Cases[i]+'-'
        # Plot age,zeta_thermal, and Delta log10(M) vs. M
        BM = False
        PlotMRKtable(MRKtable[i],title,prefix,outFolder,output,BM)

    # Next, make plots with two values of X.
    title = r'$X_{\rm core} = $'+Cases[0] +', '+Cases[1] 
    prefix = 'Xcomparison-'
    # Plot R1/R0, overlay zeta_thermal, overlay timescales vs. M
    PlotMRKtables(MRKtable,Cases,title,prefix,outFolder,output)
  
    print('\nMESAPlotXvM finished')

"""
Description of results:

Valued of mass evaluated:
Plan: log10(M/Msun) = [-0.80,-0.78, ..., 1.20 ] (N=101)
Actual: three missing values (in both runs): log10(M) = -0.66, +1.18, and +1.20,
(corresponding to 0.2188, 15.1356, 15.8489 Msun.)
Actual min,max mass: 0.1585, 14.4544 Msun, 98 rows.

Values of radius computed:
Between ZAMS (Xcore=0.69696) and Xcore=0.6, radii increased by 3 to 13%, hovering around 4% in mass
range 0.3-1.5 Msun. As the changes are gradual, zeta_thermal changed little. The larger changes in
radius occur where there are convective cores. This is not a surprise as they must burn more to achieve the
same Xcore.

Convection zones mass ranges:
                     Xcore
Convection type:   0.69696        0.6
Full               0.1585-0.2884 0.1585-0.2512
Split              0.3020-0.3631 0.2630-0.3311
End of envelope    1.6596        1.7378
Beginning of core  1.5849        1.4454

"""
