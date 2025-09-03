# -*- coding: utf-8 -*-

"""
Created on Tues Jan 23 2024
@author: skwon
"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

def setup(ANALYSISPATH):
    if not os.path.isdir(ANALYSISPATH + "/EDSMINMAXDATA"):
        os.mkdir(ANALYSISPATH + "/EDSMINMAXDATA")
        
    if not os.path.isdir(ANALYSISPATH + "/ELEMETALEXCESSDATA"):
        os.mkdir(ANALYSISPATH + "/ELEMETALEXCESSDATA")
        
    if not os.path.isdir(ANALYSISPATH + "/ELEMETALEXCESSMINMAXDATA"):
        os.mkdir(ANALYSISPATH + "/ELEMETALEXCESSMINMAXDATA")
        
    if not os.path.isdir(ANALYSISPATH + "/ELEMETALEXCESSPLOTS"):
        os.mkdir(ANALYSISPATH + "/ELEMETALEXCESSPLOTS")
    
    if not os.path.isdir(ANALYSISPATH + "/EDSPLOTS"):
        os.mkdir(ANALYSISPATH + "/EDSPLOTS")

  
    return()

# "IMPORT THE DATA" FUNCTION
def importthedata(datafolder, filename, linenumber, element, atomicvolume):
    """
    # Imports a linescan text file from the datafolder
    # Data folder = Path to the folder which holds the text file
    # filename = ID/name of the text file for a specific element line scan. 
        ASSUMED LINE SCAN STRUCTURE: column 1 (position in units of m) 
                                     column 2 (concentration in units of atomic fraction)
    # linenumber = the linescan number for the filename ID. 
    # element = name of the element asscociated with the line scan
    
    #### PATH to the data would then be: datafolder + '/' + filename + "_line" + str(linenumber) + "-" + element
    #### Example: EDSScanData/0001_line1_Cr.txt
            
    # atomicvolume = the atomic volume of the material in units of atoms/nm3
    """
    PATH = datafolder + '/' + filename + "_line" + str(linenumber) + "-" + element
    

    # IMPORT
    Ps = np.loadtxt(PATH + ".txt",delimiter=";",usecols=[0],  dtype = np.float64)*1e9
    X_atper = np.loadtxt(PATH  + ".txt",delimiter=";",usecols=[1], dtype = np.float64)*100
    X_frac = np.loadtxt(PATH  + ".txt",delimiter=";",usecols=[1], dtype = np.float64)
    
    initialmatrix = [Ps, X_atper]
    
    # CONVERT THE DATA
    X_atnm3 = X_frac*atomicvolume
    # OUTPUTS = THE  DATA CENTERED & CONVERTED 
    matrix = [Ps, X_atnm3]
    

    return initialmatrix, matrix


def findingtheindexofminmax(matrix, L_guess_nm_alongline, weight, expectedminormax):
    """
    Centers the line scan using the minimum or maximum of the element, X
    intialmatrix = from importthedata(), the original matrix with X in units of atomic percent
    matrix = from importthedata(), the converted matrix with X in untis of atoms/nm3
    L_guess_nm_alongline = the estimated length in nm where the line scan intercepts the GB
    weight = a fraction that helps set the range before and after the  L_guess_nm_alongline self estimate
    expectedminormax = either "min" or "max", dictates how the line scan will be centered
    """
    Ps = matrix[0]
    X_atnm3 = matrix[1]
    
    #### Based on the visual estimation of where the GB is - you get a L_nm
    # need to find the index of where the position matrix is closest to this length
    index_GBguess = (np.abs(Ps - L_guess_nm_alongline)).argmin()
    lb_idx = round(index_GBguess - index_GBguess*weight) # limits the region of interest
    rb_idx = round(index_GBguess + index_GBguess*weight) # limits the region of interest

 # CENTERS THE POSITION DATA BASED ON MIN/MAX IF DEPLETED/ENRICHED INPUT
        # finding the index of Position where an element is most segreated
    
    if expectedminormax == "min":
        # minimum
        indexingminmax = lb_idx + np.argmin(X_atnm3[lb_idx:rb_idx])
    if expectedminormax == "max":
        # maximum
        indexingminmax = lb_idx + np.argmax(X_atnm3[lb_idx:rb_idx])
    
    return indexingminmax


def findingtheindexofminmax_specifyingelement(datafolder, filename, linenumber, elementforminmaxcenter, atomicvolume, L_guess_nm_alongline, weight, elementexpectedminormax):
    """
    Centers the line scan using the minimum or maximum of the element, X
    intialmatrix = from importthedata(), the original matrix with X in units of atomic percent
    matrix = from importthedata(), the converted matrix with X in untis of atoms/nm3
    L_guess_nm_alongline = the estimated length in nm where the line scan intercepts the GB
    weight = a fraction that helps set the range before and after the  L_guess_nm_alongline self estimate
    expectedminormax = either "min" or "max", dictates how the line scan will be centered
    """
    initialmatrix, matrix = importthedata(datafolder, filename, linenumber, elementforminmaxcenter, atomicvolume)
    
    Ps = matrix[0]
    X_atnm3 = matrix[1]
    
    #### Based on the visual estimation of where the GB is - you get a L_nm
    # need to find the index of where the position matrix is closest to this length
    index_GBguess = (np.abs(Ps - L_guess_nm_alongline)).argmin()
    lb_idx = round(index_GBguess - index_GBguess*weight) # limits the region of interest
    rb_idx = round(index_GBguess + index_GBguess*weight) # limits the region of interest

 # CENTERS THE POSITION DATA BASED ON MIN/MAX IF DEPLETED/ENRICHED INPUT
        # finding the index of Position where an element is most segreated
    
    if elementexpectedminormax == "min":
        # minimum
        indexingminmax = lb_idx + np.argmin(X_atnm3[lb_idx:rb_idx])
    if elementexpectedminormax == "max":
        # maximum
        indexingminmax = lb_idx + np.argmax(X_atnm3[lb_idx:rb_idx])
    
    return indexingminmax

def centerthedatausingelement(initialmatrix, matrix, indexingminmax):

    Ps = matrix[0]
    X_atnm3 = matrix[1]
    X_atper = initialmatrix[1]
    
    Ps = Ps - Ps[indexingminmax]
    matrixshifted = [Ps, X_atnm3]
    initalmatrixshifted = [Ps, X_atper]

    return initalmatrixshifted, matrixshifted



def recordingminmaxEDSdata(ANALYSISPATH, filename, linenumber, element, initialmatrix, indexingminmax):
    '''
    This saves the position of the element min/max and the element's minimum or maximum in units of atomic percent
    This is saved in a text file corresponding to the filename. 
    useful if you want to extract and plot this information
    '''
    Ps = initialmatrix[0]
    Xatper = initialmatrix[1]
    Psminmax = Ps[indexingminmax]
    Xatperminmaxat = Xatper[indexingminmax]
    
    data = pd.DataFrame({'Position (nm)': Psminmax,'Conecentration Extrema (atomic percent)': Xatperminmaxat}, index=[0])
    data.to_csv(ANALYSISPATH + "/" + "EDSMINMAXDATA" + "/"  + filename + "_line" + str(linenumber) + "-" + element + "_EDSMinMaxData" + ".txt", index=False)

    
def determineaveragingdistance(initialmatrix, fraction):
    Ps = initialmatrix[0]
    averaginglength = round(fraction*Ps[-1])
    return averaginglength



# "DECIDING THE AVERAGING BOUNDS INDEX" FUNCTION
def decideavgbound(matrix, averaginglength):
    Ps = matrix[0]
    stepsize = Ps[1]-Ps[0]
    Nstep =  averaginglength/stepsize
    left = round(Nstep)
    right = len(Ps) - 1 -left
    return left, right




# "SINGLE LINE AVERAGE" FUNCTION
def singlelineaverage(initalmatrix, matrix, averaginglength):
    ''' 
    Takes an average of the line scan a distance away from the GB. 
    initalmatrixshifted = from centerthedatausingelement(), the matrix with X in untis of atomic percent, shifted so the GB is at x = 0 
    matrixshifted = from centerthedatausingelement(), the converted matrix with X in untis of atoms/nm3, shifted so the GB is at x = 0 
    averaginglength = the distance from each side of the line scan which is used to take the average concentration 
    (example: N = 50 nm means the first 50 nm and the last 50 nm will be used)
    '''
    

    left, right = decideavgbound(matrix, averaginglength)

    # FINDS AVERAGE BASED ON DISTANCE BOUNDS SPECIFED BY USER
    # USING THE x CONCENTRATION IN UNITS OF ATOMIC PERCENT
    X_atper = initalmatrix[1]
    ending = len(X_atper)- 1 # the index of the end of the data 
    
    # Aquiring the average of the "flat"/"bulk" area 
    summation_atper = sum(X_atper[0:left]) + sum(X_atper[right:ending])
    numberofpoints_atper = len(X_atper[0:left]) + len(X_atper[right:ending])
    average_atper_onelinescan = summation_atper/numberofpoints_atper
    

    # FINDS AVERAGE BASED ON DISTANCE BOUNDS SPECIFED BY USER
    # USING THE x CONCENTRATION IN UNITS OF ATOMS/NM3
    X_atnm3 = matrix[1]
    ending = len(X_atnm3)- 1 # the index of the end of the data 
    
    # Aquiring the average of the "flat"/"bulk" area 
    summation_atnm3 = sum(X_atnm3[0:left]) + sum(X_atnm3[right:ending])
    numberofpoints_atnm3 = len(X_atnm3[0:left]) + len(X_atnm3[right:ending])
    average_atnm3_onelinescan = summation_atnm3/numberofpoints_atnm3
    

    return average_atper_onelinescan, average_atnm3_onelinescan




# "GETTING BULK AVERAGE" FUNCTION
def gettingbulkaverage(datafolder, filename, totallinenumber, element, atomicvolume, L_guess_nm_alongline, weight, averaginglength ):
    
    
    
    # INITALIZE MATRIXS OF IMPORTANCE  
    average_atper_multiplelinescans = np.zeros(totallinenumber)
    average_atnm3_multiplelinescans = np.zeros(totallinenumber)
    
    # GET THE AVERAGE FOR ALL OF THE LINES USING THE "SINGLE LINE AVERAGE" FUNCTION
        # THIS SHOULD REQUIRE THE USE OF A FOR LOOP FOR ALL THE LINES
        # WILL GIVE AN ARRAY OF AVERAGES
    for M in range(0,totallinenumber):
        M = M + 1

        desiredline = M
        #print(basetextnamefront + str(M) + basetextback)
       # initialmatrix, matrix, textname = importthedata(filepathtodata, scannumber, basetextnamefront, basetextback,expectedminormax, desiredline, indexingminmax)
        initialmatrix, matrix = importthedata(datafolder, filename, desiredline, element, atomicvolume)
       
        
        # AVERAGE THE AVERAGE ARRAY - GIVING ONE VALUE
        average_atper_multiplelinescans[M-1], average_atnm3_multiplelinescans[M-1] = singlelineaverage(initialmatrix, matrix, averaginglength)

    bulkaverage_atnm3 = np.mean(average_atnm3_multiplelinescans)
    bulkinitialaverage_atper = np.mean(average_atper_multiplelinescans)
    # OUTPUTS = BULK AVERAGE FOR CONVERTED CASE AND INITIAL CASE
    return bulkinitialaverage_atper, bulkaverage_atnm3
    









# "SUBTRACTING THE BULK AVERAGE" FUNCTION
def subtractingbulkavg(matrixshifted, bulkaverage_atnm3):
    '''
    subtracts the bulk average from the marix shifted 
    matrixshifted = from centerthedatausingelement(), the converted matrix with X in untis of atoms/nm3, shifted so the GB is at x = 0 
    bulkaverage = the average concentration of the bulk away from the GB, defined either with singlelineaverage() or gettingbulkaverage()
    this is done for a specific element
    '''
    Ps = matrixshifted[0] 
    X_atnm3 = matrixshifted[1] 
    # SUBTRACT THE BULK AVERAGE FROM THE COMPOSITION DATA
    X_minusbulkavg = X_atnm3 - bulkaverage_atnm3
    
    matrixshifted_minusbulkavg = [Ps, X_minusbulkavg]
    # OUTPUTS = THE DATA NEWLY SUBTRACTED
    


    
    return matrixshifted_minusbulkavg



# "INTREGRATE THE DATA" FUNCTION
def integratingthedata(ANALYSISPATH, filename, linenumber, element, matrixshifted_minusbulkavg, expectedminormax, linescanhalfdistance, integratecutoff, plotcolor, viewplots):
    # WRITTEN FOR ONE LINE EACH 
    # INTEGRATES THE DATA IN A STEP-PER-STEP FASHION

    import scipy.integrate as integrate
    from scipy.signal import argrelextrema
    
    Ps = matrixshifted_minusbulkavg[0]
    X = matrixshifted_minusbulkavg[1]
    integratestepsize = 1
    

    maxorminlocation = 0 # 0 nm, because the GB is shifted to zero
    #A = 120 # an estimation of the distance half the size of the line scan linescanhalfdistance
    C = int(linescanhalfdistance/integratestepsize) # the amount of steps in the for loop
    D = C + 1 # needed to help predefine an array
    
    # pre-defining the arrays
    left_x = np.zeros((D), dtype=int)
    right_x = np.zeros((D), dtype=int)
    var_lx = np.zeros((C), dtype=int)
    var_rx = np.zeros((C), dtype=int)
    RIS = np.zeros(C)
    Dist = np.zeros(C)
    

    # preparing for the for loop 
    
    x = 0
    left_x[x] = maxorminlocation
    right_x[x] = maxorminlocation
    
    for x in range(0,C):
        x = x + 1
        left_x[x] = left_x[x-1] - integratestepsize
        right_x[x] = right_x[x-1] + integratestepsize
        
        [var_lx[x-1], var_rx[x-1]]= findnearest(Ps, left_x[x], right_x[x])
        
        Dist[x-1] = Ps[var_rx[x-1]] - Ps[ var_lx[x-1]]
        RIS[x-1] = integrate.simps(X[var_lx[x-1]:var_rx[x-1]],Ps[var_lx[x-1]:var_rx[x-1]], even='avg')

         
    
    ## determinging a cut off of when the local min/max elemeental excess 
    #cut = 150 # nm integratecutoff
    cutoffindex = (np.abs(Dist - integratecutoff)).argmin()
 
    if expectedminormax == "min":
        # minimum
        q = argrelextrema(RIS[0:cutoffindex], np.less_equal, order = 5) 
        
        if np.size(q) == 0:
            print('Unable to find min in ' + filename + "_line" + str(linenumber) + "-" + element + '-Smoothed')
            minmaxRIS = 0
            minmaxRIS_positions = 0
            leftbound = 0
            rightbound = 0
        else:
            minmaxRIS = np.array(RIS[q])[0]
            minmaxRIS_positions = np.array(Dist[q])[0]
            leftbound = var_lx[q][0]
            rightbound = var_rx[q][0]

    if expectedminormax == "max":
        # maximum
        q = argrelextrema(RIS[0:cutoffindex], np.greater_equal, order = 5)
        
        if np.size(q) == 0:
            print('Unable to find max in ' +  filename + "_line" + str(linenumber) + "-" + element + '-Smoothed')
            minmaxRIS = 0
            minmaxRIS_positions = 0
            leftbound = 0
            rightbound = 0
        
        else:
            minmaxRIS = np.array(RIS[q])[0]
            minmaxRIS_positions = np.array(Dist[q])[0]
            leftbound = var_lx[q][0]
            rightbound = var_rx[q][0]
    
    # everything is indexed with 0 to denote choosing the FIRST extrema in the elemental excess profile

    fig, ax = plt.subplots()
    ax.plot(Ps, X, label=filename + "_line" + str(linenumber) + "-" + element, color=plotcolor, marker='o', markersize=5)
    plt.axvline(x = Ps[leftbound], color='k', linestyle='--')
    plt.axvline(x = Ps[rightbound], color='k', linestyle='--')
    plt.axhline(y = 0, color='k', linestyle='--')
    plt.xlabel('Position [nm]', fontsize=20)
    plt.ylabel("Concentration [$nm^{3}$]", fontsize=20)
    plt.title( filename + "_line" + str(linenumber) + "-" + element, fontsize=25, y=1.08)
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18) 
    ax.tick_params(which='both', width=2)
    ax.tick_params(which='major', length=8)
    ax.tick_params(which='minor', length=4)
    figurestring = ANALYSISPATH +  "/EDSPLOTS/" +  filename + "_line" + str(linenumber) + "-" + element + "-4_profile_converted_GB_centered_minusbulkavg_integrated_bounds_at_maxormin" + '.png'
    plt.savefig(figurestring, bbox_inches='tight')
    if viewplots == "no":
        plt.close()
    
    
    fig, ax = plt.subplots()
    ax.scatter(Dist,RIS,label= filename + "_line" + str(linenumber) + "-" + element, color=plotcolor, alpha=1)
    plt.scatter(minmaxRIS_positions, minmaxRIS,  alpha = 1)
    plt.axhline(y=0, color='k', linestyle='--')
    plt.xlabel('Separation of Bounds [nm]', fontsize=20)
    plt.ylabel('Elemental Excess [$nm^{-2}$]', fontsize=20)
    plt.title('RIS(bounds) for ' + filename + "_line" + str(linenumber) + "-" + element, fontsize=25, y=1.08)
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18) 
    ax.tick_params(which='both', width=2)
    ax.tick_params(which='major', length=8)
    ax.tick_params(which='minor', length=4)
    figurestring = ANALYSISPATH + "/ELEMETALEXCESSPLOTS/" +  filename + "_line" + str(linenumber) + "-" + element + '.png'
    plt.savefig(figurestring, bbox_inches='tight')
    if viewplots == "no":
        plt.close()
    
    
    ######### SAVE THE DATA

    data = pd.DataFrame({'Distance of Bounds (nm)': Dist,'Elemental Excess (atoms/nm2)':RIS})
    data.to_csv(ANALYSISPATH + "/ELEMETALEXCESSDATA/" + filename + "_line" + str(linenumber) + "-" + element + "_ELEMENTALEXCESSData" + ".txt", index=False,)
    data2 = pd.DataFrame({'Boundary Distance for Extrema Elemental Excess': minmaxRIS_positions,'Elemental Excess Extrema': minmaxRIS}, index=[0])
    data2.to_csv(ANALYSISPATH + "/ELEMETALEXCESSMINMAXDATA/" + filename + "_line" + str(linenumber) + "-" + element + "_ELEMENTALEXCESSMinMaxData" + ".txt", index=False,)
    return 





def findnearest(array, valuel, valuem):
    import numpy as np
    array = np.asarray(array)
    idlx = (np.abs(array - valuel)).argmin()
    idmx = (np.abs(array - valuem)).argmin()
    
    return idlx, idmx
    
    
    
    

def plotinitallinescan(ANALYSISPATH, filename, linenumber, element, initialmatrix, matrix, plotcolor, viewplots):
    

    fig, ax = plt.subplots()
    ax.plot(initialmatrix[0],initialmatrix[1], color=plotcolor, marker='o', markersize=5)
    plt.xlabel('Position [nm]', fontsize=20)
    plt.ylabel( "Concentration [at%]", fontsize=20)
    plt.title(filename + "_line" + str(linenumber) + "-" + element, fontsize=25, y=1.08)
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18) 
    ax.tick_params(which='both', width=2)
    ax.tick_params(which='major', length=8)
    ax.tick_params(which='minor', length=4)

    figurestring =  ANALYSISPATH + "/EDSPLOTS/"  +  filename + "_line" + str(linenumber) + "-" + element + "-1_initial_profile" + '.png'
    plt.savefig(figurestring, bbox_inches='tight')
    if viewplots == "no":
        plt.close()
    
    fig, ax = plt.subplots()
    ax.plot(matrix[0],matrix[1], color=plotcolor, marker='o', markersize=5)
    plt.xlabel('Position [nm]', fontsize=20)
    plt.ylabel( "Concentration [$nm^{-3}$]", fontsize=20)
    plt.title(filename + "_line" + str(linenumber) + "-" + element, fontsize=25, y=1.08)
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18) 
    ax.tick_params(which='both', width=2)
    ax.tick_params(which='major', length=8)
    ax.tick_params(which='minor', length=4)

    figurestring =  ANALYSISPATH + "/EDSPLOTS/"  +  filename + "_line" + str(linenumber) + "-" + element + "-2_profile_converted" + '.png'
    plt.savefig(figurestring, bbox_inches='tight')
    if viewplots == "no":
        plt.close()
    
    
    
    
    
def plotcenteredlinescans(ANALYSISPATH, filename, linenumber, element, matrixshifted, averaginglength, bulkaverage_atnm3, plotcolor, viewplots):
   
    left, right = decideavgbound(matrixshifted, averaginglength)

    fig, ax = plt.subplots()
    ax.plot(matrixshifted[0],matrixshifted[1], color=plotcolor, marker='o', markersize=5)
    plt.axvline(x = matrixshifted[0][left], color='k', linestyle='--')
    plt.axvline(x =  matrixshifted[0][right], color='k', linestyle='--')
    plt.axhline(y =  bulkaverage_atnm3, color='k', linestyle='--')  
    plt.xlabel('Position [nm]', fontsize=20)
    plt.ylabel( "Concentration [$nm^{-3}$]", fontsize=20)
    plt.title(filename + "_line" + str(linenumber) + "-" + element, fontsize=25, y=1.08)
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18) 
    ax.tick_params(which='both', width=2)
    ax.tick_params(which='major', length=8)
    ax.tick_params(which='minor', length=4)

    figurestring =  ANALYSISPATH + "/EDSPLOTS/"  +  filename + "_line" + str(linenumber) + "-" + element + "-3_profile_converted_GB_centered" + '.png'
    plt.savefig(figurestring, bbox_inches='tight')
    if viewplots == "no":
        plt.close()
    

    
    
    
    
    


def analysis_singlelinescan(ANALYSISPATH, datafolder, filename, linenumber, mainelement, expectedminormax, centerelement, centerelementexpectedminormax,  plotcolor,  atomicvolume, L_guess_nm_alongline, weight, fraction, linescanhalfdistance, integratecutoff, viewplots):
    setup(ANALYSISPATH)
    initialmatrix, matrix = importthedata(datafolder, filename, linenumber, mainelement, atomicvolume)
    
    plotinitallinescan(ANALYSISPATH, filename, linenumber, mainelement, initialmatrix, matrix, plotcolor, viewplots)
    #indexingminmax = findingtheindexofminmax(matrix, L_guess_nm_alongline, weight, expectedminormax)
    indexingminmax = findingtheindexofminmax_specifyingelement(datafolder, filename, linenumber, centerelement, atomicvolume, L_guess_nm_alongline, weight, centerelementexpectedminormax)
    initialmatrixshifted, matrixshifted  = centerthedatausingelement(initialmatrix, matrix, indexingminmax)
    recordingminmaxEDSdata(ANALYSISPATH, filename, linenumber, mainelement, initialmatrix, indexingminmax)
    averaginglength = determineaveragingdistance(initialmatrix, fraction)
    average_atper_onelinescan, average_atnm3_onelinescan = singlelineaverage(initialmatrix, matrix, averaginglength)
    plotcenteredlinescans(ANALYSISPATH, filename, linenumber, mainelement, matrixshifted, averaginglength, average_atnm3_onelinescan, plotcolor, viewplots)
    matrixshifted_minusbulkavg = subtractingbulkavg(matrixshifted, average_atnm3_onelinescan)
    integratingthedata(ANALYSISPATH, filename, linenumber, mainelement, matrixshifted_minusbulkavg, expectedminormax, linescanhalfdistance, integratecutoff, plotcolor, viewplots)
      


def analysis_usingmultiplelinescans(ANALYSISPATH, datafolder, filename, mainelement, mainelementexpectedminormax, centerelement ,centerelementexpectedminormax, totallinenumber, plotcolor,  atomicvolume, L_guess_nm_alongline, weight, fraction, linescanhalfdistance, integratecutoff, viewplots):
    setup(ANALYSISPATH)
    
    for M in range(0,totallinenumber):
        M = M + 1
        desiredline = M 
        initialmatrix, matrix = importthedata(datafolder, filename, desiredline, mainelement, atomicvolume)
        plotinitallinescan(ANALYSISPATH, filename, desiredline, mainelement, initialmatrix, matrix, plotcolor, viewplots)
        indexingminmax = findingtheindexofminmax_specifyingelement(datafolder, filename, desiredline, centerelement, atomicvolume, L_guess_nm_alongline, weight, centerelementexpectedminormax)
        initalmatrixshifted, matrixshifted = centerthedatausingelement(initialmatrix, matrix, indexingminmax)
        recordingminmaxEDSdata(ANALYSISPATH, filename, desiredline, mainelement, initialmatrix, indexingminmax)
        averaginglength = determineaveragingdistance(initialmatrix, fraction)

        if M == 1:
            bulkinitialaverage_atper, bulkaverage_atnm3 = gettingbulkaverage(datafolder, filename, totallinenumber, mainelement, atomicvolume, L_guess_nm_alongline, weight, averaginglength )

        plotcenteredlinescans(ANALYSISPATH, filename,desiredline, mainelement, matrixshifted, averaginglength, bulkaverage_atnm3, plotcolor,viewplots)
        matrixshifted_minusbulkavg = subtractingbulkavg(matrixshifted, bulkaverage_atnm3)
        integratingthedata(ANALYSISPATH, filename, desiredline, mainelement, matrixshifted_minusbulkavg, mainelementexpectedminormax, linescanhalfdistance, integratecutoff, plotcolor, viewplots)
          











    

