baseline = []
variability = []
interp5 = []
vm = []

with open('out.txt', 'r') as f:
    L = f.readlines()
    for l in L[3:]:
        things = l.split(' ')
        baseline.append(things[0])
        variability.append(things[1])
        interp5.append(things[2])
        vm.append(things[3])
import matplotlib.pyplot as plt

plt.plot(baseline, 'r-')
plt.plot(variability, 'g-')
plt.plot(interp5, 'b-')
plt.plot(vm, 'k-')
plt.show()