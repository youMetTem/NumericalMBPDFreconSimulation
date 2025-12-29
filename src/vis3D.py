from vpython import canvas, box, sphere, vector, color, rate, label
import numpy as np
import colorsys
import main

# setting up parameters
mlcName = "Helium"
mlcMW = 6.646 * 10**-27 # kg
mlcRadius = 3.1 * 10**(-11) # m
N = 500 # number of particles (500 for optimized visualization)
temperature = 373 # K
boundaryLength = 10**(-9) # m
boltzmannConstant = 1.38 * 10**-23 # J/K

posMat, velocityMat, vRMS, dt = main.initSys(N, temperature, mlcMW, boundaryLength)
dt = 7 * 10**(-16) # defaulting tiny dt for visualization

# setting up color spectrum assigning function
def getSpectrumColor(speed):
    normalized = speed/(2.5*vRMS)
    if normalized > 1: normalized = 1

    hue = normalized * 0.8
    r, g, b = colorsys.hsv_to_rgb(hue, 1, 1)

    return vector(r, g, b)

# setting up scene
scene = canvas(
    title = f"3D {mlcName} Collision Simulation",
    width = 1000,
    height =800,
    center = vector(boundaryLength/2, boundaryLength/2, boundaryLength/2),
    background = color.black
)

# setting up main box
container = box(
    pos = vector(boundaryLength/2, boundaryLength/2, boundaryLength/2),
    length = boundaryLength, height = boundaryLength, width = boundaryLength,
    opacity = 0.3,
    color = color.white
)

# setting up spectrum bar
barX = boundaryLength * 1.3
steps = 100; barH = boundaryLength

for i in range(steps):
    fract = i/steps
    yPos = (i*(barH/steps)) + (barH/steps)/2

    box(
        pos = vector(barX, yPos, boundaryLength/2),
        length = boundaryLength/15,
        height = barH/steps,
        width = boundaryLength/15,
        color = getSpectrumColor(fract*2.5*vRMS),
        opacity = 1
    )

label(pos = vector(barX, -boundaryLength/10, boundaryLength/2), text = "min speed", box = False, height = 20)
label(pos = vector(barX, boundaryLength*1.1, boundaryLength/2), text = "max speed", box = False, height = 20)

# setting up particles
visualParticles = []
for i in range(N):
    p = sphere(
        pos = vector(posMat[i, 0], posMat[i, 1], posMat[i, 2]),
        radius = mlcRadius,
        color = color.white,
        made_trail = False
    )
    visualParticles.append(p)

stepsPerFrame = 2

# main activation loop
while True:
    rate(60)

    for n in range(stepsPerFrame):
        result = main.update(posMat, velocityMat, dt, boundaryLength, mlcRadius)

        if result is not None:
            posMat, velocityMat = result

    for i in range(N):
        visualParticles[i].pos = vector(posMat[i, 0], posMat[i, 1], posMat[i, 2])

        speed = np.linalg.norm(velocityMat[i])
        visualParticles[i].color = getSpectrumColor(speed)