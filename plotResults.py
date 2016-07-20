import matplotlib.pyplot as plt

V = []
theta = []
theta_b = []


with open('spikes.txt', 'r') as f:
    L = f.readlines()
    i = 0
    while i < len(L):
        l = L[i]
        things = l.split(' ')

        shape = []
        t = things[0]
        ch = things[1]
        shape = things[2:]

        plt.plot(shape, 'k-')
        plt.title('Channel: ' + ch)
        plt.show()

        i += 1
        
        f, axarr = plt.subplots(3, 3)
        for j in range(9):
            l = L[i]
            things = l.split(' ')
            shape = things[1:]
            
            axarr[j/3, j%3].plot(shape)
            axarr[j/3, j%3].set_xlim([0, len(shape)])
            axarr[j/3, j%3].set_ylim([1900, 2200])
            i += 1
            
        plt.show()
        



