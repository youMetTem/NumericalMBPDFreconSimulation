import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.spatial.distance import pdist, squareform
from scipy.stats import gaussian_kde
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
from sklearn.metrics import mean_squared_error

# setting up parameters
mlcName = "Helium"
mlcMW = 6.646 * 10**-27 # kg
mlcRadius = 3.1 * 10**(-11) # m
N = 2000 # number of particles
temperature = 373 # K
boundaryLength = 10**(-9) # m
boltzmannConstant = 1.38 * 10**-23 # J/K

# initialize random position, velocity matrix, v_rms, and dt
def initSys(N, temperature, mlcMW, boundaryLength):

    # setting up initial positions in a grid layout
    # preventing overlapping placed particle
    particlePerSide = int(np.ceil(N**(1/3)))
    particleDistance = boundaryLength / particlePerSide

    if particleDistance < 2 * mlcRadius: # checking whether the box is large enough
        print("Overcrowded, box is too small.")
        return

    posMat = np.zeros((N, 3))
    placingCount = 0

    # placing particle
    for x in range(particlePerSide):
        for y in range(particlePerSide):
            for z in range(particlePerSide):
                if placingCount >= N: break

                posMat[placingCount] = np.array([
                    x*particleDistance + particleDistance/2,
                    y*particleDistance + particleDistance/2,
                    z*particleDistance + particleDistance/2
                    ])
                
                # adding tiny randomness
                posMat[placingCount] += np.random.uniform(-particleDistance/10, -particleDistance/10, 3)
                placingCount += 1
            if placingCount >= N: break
        if placingCount >= N: break

    # setting up initial velocities
    velocityMat = np.random.uniform(-1, 1, (N, 3))

    # scaling velocities to temperature parameter
    Vsqaure = np.sum(velocityMat**2, axis = 1)
    totKE = mlcMW * np.sum(Vsqaure) / 2
    avgKE = totKE/N
    currentTemp = (2/3) * (avgKE/boltzmannConstant)
    scalingFactor = np.sqrt(temperature/currentTemp)

    velocityMat = velocityMat * scalingFactor

    # setting up v_rms
    vRMS = np.sqrt(3 * boltzmannConstant * temperature / mlcMW)

    # setting up safe dt for calculation, avoiding tunneling
    maxV = np.max(np.linalg.norm(velocityMat, axis=1))
    dt = (mlcRadius / maxV) * 0.2

    return posMat, velocityMat, vRMS, dt

# theoretical maxwell-boltzmann PDF
def maxwellBoltzmannPDF(v, m, T,):
    A = (m / (2 * np.pi * boltzmannConstant * T))**(3/2)
    exponent = - (m * v**2) / (2 * boltzmannConstant * T)
    return 4 * np.pi * A * (v**2) * np.exp(exponent)

# visualizing PDF compared to theoretical PDF
def initPDF(velocity):
    plt.figure(figsize = (9, 6))

    # setting up axis limit and calculating speed
    speedMat = np.linalg.norm(velocity, axis = 1)
    maxPlotSpeed = vRMS * 3

    # plotting theoretical PDF
    vPoints = np.linspace(0, maxPlotSpeed, 500)
    PDFTheory = maxwellBoltzmannPDF(vPoints, mlcMW, temperature)
    plt.plot(vPoints, PDFTheory, "r-", lw = 3, alpha = 0.8, label = "Theoretical Distribution")

    # plotting initial state
    bins = np.linspace(0, maxPlotSpeed, 100)
    plt.hist(speedMat, bins = bins, density = True, color = "skyblue", alpha = 0.75, edgecolor = "black", linewidth = 0.5, label = "Initial Simulation Data")

    plt.xlim(0, maxPlotSpeed)
    plt.ylim(0, np.max(PDFTheory) * 1.2)
    plt.xlabel("Speed [m/s]")
    plt.ylabel("Probability Density")
    plt.title(f"Initial Speed Distribution of {mlcName} at T = {temperature} K")
    plt.legend(loc = "upper right")
    plt.grid(True, alpha = 0.2)

    plt.show()

    return

# update function, output the new positions and velocities matrix
def update(pos, velocity, dt, boundaryLength, mlcRadius):

    #update positions
    pos += velocity * dt

    # handle wall collisions
    leftCollide = pos[:, 0] < mlcRadius
    rightCollide = pos[:, 0] > (boundaryLength - mlcRadius)
    velocity[leftCollide | rightCollide, 0] *= -1

    bottomCollide = pos[:, 1] < mlcRadius
    topCollide = pos[:, 1] > (boundaryLength - mlcRadius)
    velocity[bottomCollide | topCollide, 1] *= -1

    frontCollide = pos[:, 2] < mlcRadius
    backCollide = pos[:, 2] > (boundaryLength - mlcRadius)
    velocity[frontCollide | backCollide, 2] *= -1

    # handle particle collisions
    distanceMat = squareform(pdist(pos))
    collidePair1, collidePair2 = np.where((distanceMat < 2 * mlcRadius))
    uniqueFilter = collidePair1 < collidePair2

    for i, j in zip(collidePair1[uniqueFilter], collidePair2[uniqueFilter]):
        relativePos = pos[i] - pos[j]
        relativeVelocity = velocity[i] - velocity[j]
        distanceSquare = np.sum(np.square(relativePos))

        # checking if particle is moving toward each other or not
        if np.dot(relativePos, relativeVelocity) < 0:
            deltaVelocity = relativePos * (np.dot(relativePos, relativeVelocity)/distanceSquare)

            velocity[i] -= deltaVelocity
            velocity[j] += deltaVelocity

    # prevent wall clipping (double check)
    pos[pos < mlcRadius] = mlcRadius
    pos[pos > boundaryLength - mlcRadius] = boundaryLength - mlcRadius
    
    return pos, velocity

# visualizing simulated PDF
def runSimulationHist():
    fig, ax = plt.subplots(figsize = (9, 6))

    # setting up axis limit to be 3 times the vRMS
    maxPlotSpeed = vRMS * 3

    # plotting theoretical maxwell-boltzmann PDF
    vPoints = np.linspace(0, maxPlotSpeed, 500)
    PDFTheory = maxwellBoltzmannPDF(vPoints, mlcMW, temperature)
    ax.plot(vPoints, PDFTheory, "r-", lw = 3, alpha = 0.8, label = "Theoretical Distribution")

    # setting up PDF histogram
    bins = np.linspace(0, maxPlotSpeed, 45)
    binWidth = bins[1] - bins[0]

    m, n, barContainer = ax.hist([], bins = bins, density = True, color = "skyblue", alpha = 0.6, ec = "black", label = "Simulation Data")

    ax.set_xlim(9, maxPlotSpeed)
    ax.set_ylim(0, np.max(PDFTheory) * 1.2) 
    ax.set_xlabel("Speed [m/s]")
    ax.set_ylabel("Probability Density")
    ax.set_title(f"He Atom Speed PDF at {temperature}K")
    ax.legend()

    # animation loop
    def animate(frame):
        global posMat, velocityMat

        # showing initial state
        if frame > 0:
            for n in range(20):
                posMat, velocityMat = update(posMat, velocityMat, dt, boundaryLength, mlcRadius)
        else: pass

        speeds = np.linalg.norm(velocityMat, axis = 1)
        counts, n = np.histogram(speeds, bins = bins)

        # Normalizing counts to density
        densityValue = counts / (N * binWidth)

        # updating rectagle height
        for rect, h in zip(barContainer.patches, densityValue):
            rect.set_height(h)

        return barContainer.patches

    anim = FuncAnimation(fig, animate, frames = 300, interval = 50, blit = False)
    plt.show()

    return

# visualizing curve of PDF instead of histogram
def runSimulationCurve():
    fig, ax = plt.subplots(figsize = (9, 6))
    maxPlotSpeed = vRMS * 3

    vPoints = np.linspace(0, maxPlotSpeed, 500)
    PDFTheory = maxwellBoltzmannPDF(vPoints, mlcMW, temperature)
    ax.plot(vPoints, PDFTheory, "r-", lw = 3, alpha = 0.6, label = "Theoretical Distribution")

    simulationLine, = ax.plot([], [], color='#1f77b4', lw = 2, alpha = 0.8, markersize = 5, label = "Simulation Data")
    fillPoly = ax.fill_between(vPoints, 0, 0, color='#1f77b4', alpha =0.2)

    ax.set_xlim(0, maxPlotSpeed)
    ax.set_ylim(0, np.max(PDFTheory) * 1.2)
    ax.set_xlabel("Speed [m/s]")
    ax.set_ylabel("Probability Density")
    ax.set_title(f"He Atom Speed Distribution (Gaussian KDE) at T={temperature}K")
    ax.legend(loc = "upper right")
    ax.grid(True, alpha = 0.3)

    def animate(frame):
        global posMat, velocityMat
        nonlocal fillPoly

        for n in range(20):
            posMat, velocityMat = update(posMat, velocityMat, dt, boundaryLength, mlcRadius)
        
        speeds = np.linalg.norm(velocityMat, axis = 1)
        kde = gaussian_kde(speeds, bw_method = "scott")
        densityValue = kde(vPoints)

        simulationLine.set_data(vPoints, densityValue)
        fillPoly.remove()
        fillPoly = ax.fill_between(vPoints, densityValue, 0, color='#1f77b4', alpha = 0.2)

        return simulationLine, fillPoly

    anim = FuncAnimation(fig, animate, frames = 300, interval = 50, blit = False)
    plt.show()

    return

def hypoFuncFit(velocityMat):

    # definining hypothesized/generic form function
    def hypoFunc(v, a, b):
        return a * v**2 * np.exp(-b * v**2)

    speeds = np.linalg.norm(velocityMat, axis = 1)
    
    maxV = np.max(speeds)
    vAxis = np.linspace(0, maxV * 1.2, 200)

    # sorting data into histogram to obtain y value for the data fitting
    histDensity, binEdges = np.histogram(speeds, bins = "auto", density = True)
    binCenters = (binEdges[:-1] + binEdges[1:])/2

    # guessing initial value for A, B
    peakIDX = np.argmax(histDensity)
    vPeak = binCenters[peakIDX]
    if vPeak == 0: vPeak = 0.1

    guessedB = 1 / (vPeak**2)
    guessedA = np.max(histDensity/(vPeak**2 * np.exp(-1)))

    # fitting data with hypothesized function
    params, cov = curve_fit(hypoFunc, binCenters, histDensity, p0 =[guessedA, guessedB])
    predictedA, predictedB = params

    # calculating theoretical A, B values
    theoryA = 4*np.pi*((mlcMW)/(2*np.pi*boltzmannConstant*temperature))**(3/2)
    theoryB = mlcMW/(2*boltzmannConstant*temperature)

    # comparing theory and predicted A, B values
    errorA = ((theoryA - predictedA)/theoryA)*100
    errorB = ((theoryB - predictedB)/theoryB)*100

    # calculating residuals
    yhat = hypoFunc(vAxis, predictedA, predictedB)
    yTheory = maxwellBoltzmannPDF(vAxis, mlcMW, temperature)
    residuals = yhat - yTheory

    # calculating MSE
    mse = mean_squared_error(yTheory, yhat)

    print(f"---------------------------\nHypothesized Generic Function depicted as : A * v**2 * exp(-B * v**2)\n\npredictedA = {predictedA}\npredictedB = {predictedB}\n\ntheoryA = {theoryA}\ntheoryB = {theoryB}\n\nPercentage Error A = {errorA} %\nPercentage Error B = {errorB} %\n\nMSE = {mse}\n---------------------------\n")

    # plotting data
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (10, 8), sharex = True, gridspec_kw = {"height_ratios":[3, 1]})

    # distribution plot
    ax1.hist(speeds, bins = 50, density = True, alpha = 0.5, color = "gray", label = "Raw Simulation Data")
    ax1.scatter(binCenters, histDensity, color = "black", s = 10, marker = "x", alpha =0.6, label = "Histogram Bin Centers")
    ax1.plot(vAxis, yTheory, "b--", linewidth = 4, alpha = 0.8, label = "Theoretical Distribution")
    ax1.plot(vAxis, yhat, "r-", linewidth = 2, alpha = 0.7, label = "Estimated Function\n(From Hypothesized Generic Function)")

    ax1.set_ylabel("Probability Density")
    ax1.set_title("Maxwell-Boltzmann Distribution: Generic Function Curve Fitted")
    ax1.legend()
    ax1.grid(True, alpha = 0.3)

    # residuals plot
    ax2.scatter(vAxis, residuals, color = "crimson", s = 10, alpha = 0.6)
    ax2.axhline(0, color = "black", linestyle = "--", linewidth = 1)

    ax2.set_ylim(np.min(residuals)*1.5, np.max(residuals)*1.5)
    ax2.set_ylabel("Residuals")
    ax2.set_title("Residual Plot")
    ax2.grid(True, alpha =0.3)

    plt.tight_layout
    plt.show()

    return predictedA, predictedB

# fitting data with cubic spline after sorting data into histogram
def histCubicSplineFit(velocityMat):

    speeds = np.linalg.norm(velocityMat, axis = 1)
    maxV = np.max(speeds)
    vAxis = np.linspace(0, maxV*1.2, 200)

    # sorting data into histogram to obtain y value for the interpolation
    counts, binEdges = np.histogram(speeds, bins = 50, density = True)
    binCenters = (binEdges[:-1] + binEdges[1:])/2
    yObserved = counts

    # fitting data points into cubic spline
    cubicSplineFunc = CubicSpline(binCenters, yObserved)
    yhat = cubicSplineFunc(vAxis)

    # cutting off out of range data (below v_min, above v_max)
    yhat[vAxis < np.min(speeds)] = 0
    yhat[vAxis > np.max(speeds)] = 0
    yhat[yhat < 0] = 0

    # setting up theoretical curve
    yTheory = maxwellBoltzmannPDF(vAxis, mlcMW, temperature)

    # calculating residuals
    residuals = yhat - yTheory

    # calculating MSE
    mse = mean_squared_error(yTheory, yhat)

    print(f"---------------------------\nMSE = {mse}\n\nEvaluate probability of particle being in speed v using:\nreturned function [cubicSplineFunc]\n---------------------------\n")

    # plotting
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (10, 8), sharex = True, gridspec_kw = {"height_ratios":[3, 1]})

    ax1.hist(speeds, bins = 50, density = True, alpha = 0.3, color = "gray", label = "Raw Simulation Data")
    ax1.plot(vAxis, yTheory, "b--", linewidth = 2, alpha =0.8, label = "Theoretical Distribition")
    ax1.plot(vAxis, yhat, "r-", linewidth = 2.5, alpha = 0.7, label = "Cubic Spline Interpolation")
    ax1.scatter(binCenters, yObserved, color = "black", s = 10, marker = "x", alpha = 0.6, label = " Histogarm Bin Centers")

    ax1.set_ylabel("Probability Density")
    ax1.set_title(f"Maxwell-Boltzmann: Histogram's Bin Midpoint Spline Interpolated Curve")
    ax1.legend()
    ax1.grid(True, alpha = 0.3)

    ax2.scatter(vAxis, residuals, color = "crimson", s = 10, alpha = 0.6)
    ax2.axhline(0, color = "black", linestyle = "--", linewidth = 1)

    ax2.set_ylim(np.min(residuals)*1.5, np.max(residuals)*1.5)
    ax2.set_ylabel("Residuals")
    ax2.set_xlabel("Speed [m/s]")
    ax2.set_title("Residual Plot")
    ax2.grid(True, alpha = 0.3)

    plt.tight_layout()
    plt.show()

    return cubicSplineFunc

# main
def main():

    while True:
        print(f"====={mlcName} Molecules Simulation=====\n------------------------------------------\nPlease Select Mode of Computation\n------------------------------------------\n[1] Visualizing Initial PDF Compared to Theoretical PDF\n[2] Visualize Dynamic Simulated PDF Histogram\n[3] Visualize Dynamic Simulated PDF (KDE curve)\n[4] Visualize Estimated Function From Hypothesized Generic Function\n[5] Visualize Histogram's Bin Midpoint Spline Interpolated Curve\n[6] Exit\n<*> Default mode is set to Helium, change this in main.py, parameters block\n<*> Mode [4], [5] might take a few second\n------------------------------------------")
        
        mode = int(input("Selected Mode: ").strip())

        if mode not in [1, 2, 3, 4, 5, 6]:
            print("Enter Provided Mode [1] - [5]\n\n\n")
            continue

        if mode == 1: initPDF(velocityMat); print("Successfully Executed")
        if mode == 2: runSimulationHist(); print("Successfully Executed")
        if mode == 3: runSimulationCurve(); print("Successfully Executed")
        if mode in [4, 5]:

            print("Loading ...")

            for i in range(1000):
                update(posMat, velocityMat, dt, boundaryLength, mlcRadius)

            stackedVelocities = []
            Nsampling = 50
            samplingStep = 50

            for i in range(Nsampling):
                for n in range(samplingStep):
                    update(posMat, velocityMat, dt, boundaryLength, mlcRadius)
                stackedVelocities.append(velocityMat.copy())

            fullVelocitiesData = np.vstack(stackedVelocities)
            
            if mode == 4:
                predictedA, predictedB = hypoFuncFit(fullVelocitiesData)

            if mode == 5:
                cubicSplineFunc = histCubicSplineFit(fullVelocitiesData)
            
            print("Successfully Executed")

        if mode == 6: print("Closing ..."); break

        print("\n\n\n")

# execute code
if __name__ == "__main__":
    posMat, velocityMat, vRMS, dt = initSys(N, temperature, mlcMW, boundaryLength)
    main()