import matplotlib.pyplot as plt

V = []
theta = []
theta_b = []

with open('DEBUG_FIND_OUTPUT.txt', 'r') as f:
    for l in f.readlines():
        things = l.split(' ')
        V.append(things[0])
        theta.append(things[1])
        theta_b.append(things[2])


plt.plot(V, 'r-')
plt.plot(theta, 'b--')
plt.plot(theta_b, 'r--')
plt.title('Interpolated voltage + boundaries')
plt.show()

##

Qdiff = []
vMovingAvg = []
vGlobalMovingAvg = []
baseline = []
variability = []
vGlobal = []
Qmin = []

with open('DEBUG_OUTPUT.txt', 'r') as f:
    for l in f.readlines():
        things = l.split(' ')
        vMovingAvg.append(things[0])
        vGlobalMovingAvg.append(things[1])
        baseline.append(things[2])
        variability.append(things[3])
        Qdiff.append(things[4])
        vGlobal.append(things[5])
        Qmin.append(things[6])


plt.plot(vGlobalMovingAvg, 'r-')
plt.plot(vGlobal, 'b-')
plt.title('vGlobalMovingAvg and vGlobal')
plt.show()

plt.plot(baseline, 'r-')
plt.plot(vMovingAvg, 'b-')
plt.title('Baseline and vMovingAvg')
plt.show()

plt.plot(variability, 'r-')
plt.title('variability')
plt.show()


plt.plot(Qdiff, 'r-')
plt.title('Qdiff')
plt.show()

plt.plot(Qmin, 'r-')
plt.title('Qmin')
plt.show()