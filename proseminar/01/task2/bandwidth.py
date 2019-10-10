import matplotlib.pyplot as plt
import csv

def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

def getDataFromFile(file):
    size = []
    latency = []
    with open(file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        for row in csv_reader:
            size.append(num(row[0]))
            latency.append(num(row[1]))
    return size, latency
    
    
sockets = getDataFromFile("bw/2sockets")
cores = getDataFromFile("bw/2cores")
hosts = getDataFromFile("bw/2hosts")

plt.ylim(0, 6500)
plt.xlim(0, 5000000)

plt.xlabel('Bandwitdh in MB/s')
plt.ylabel('time in us')
plt.plot(hosts[0], hosts[1], label = '2 hosts')
plt.plot(cores[0], cores[1], label = '2 cores')
plt.plot(sockets[0], sockets[1], label = '2 sockets')
plt.legend()
#plt.show()
plt.savefig('bandwidthTest.png')
