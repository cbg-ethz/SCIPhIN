import sys

monovarIn = sys.argv[1]
sciphiIn = sys.argv[2]
sciphinIn = sys.argv[3]

def getBinData(tool, filename, binData):

    inFile = open(filename, 'r')
    for line in inFile:
        if line[0] == '#':
            continue

        lineSplit = line.split('\t')
        bin = []
        for i in range(9, 25):
            gt = lineSplit[i].split(':')[0]
            if gt == './.':
                bin.append(None)
            elif gt == '0/0' or gt == '0|0':
                bin.append(0)
            elif gt == '0/1' or gt == '1/1':
                bin.append(1)
            else:
                print(lineSplit[i])

            pos = lineSplit[0] + "_" + lineSplit[1]
            if not (pos in binData):
                binData[pos] = {}
            binData[pos][tool] = bin

    return binData
        

binData = {}
binData = getBinData('monovar', monovarIn, binData)
binData = getBinData('sciphi', sciphiIn, binData)
binData = getBinData('sciphin', sciphinIn, binData)

counter = {'all':[0]*16, 'monovar_sciphi': [0]*16, 'monovar_sciphin': [0]*16, 'monovar': [0]*16, 'sciphi_sciphin':[0]*16, 'sciphi': [0]*16, 'sciphin': [0]*16}

for pos in binData:
    if len(binData[pos]) == 3:
        values = binData[pos]
        for i in range(16):
            if values['monovar'][i] == None:
                continue
            if values['monovar'][i] == 1 and values['sciphi'][i] == 1 and values['sciphin'][i] == 1:
                counter['all'][i] += 1
            elif values['monovar'][i] == 1 and values['sciphi'][i] == 1:
                counter['monovar_sciphi'][i] += 1
            elif values['monovar'][i] == 1 and values['sciphin'][i] == 1:
                counter['monovar_sciphin'][i] += 1
            elif values['sciphi'][i] == 1 and values['sciphin'][i] == 1:
                counter['sciphi_sciphin'][i] += 1
            elif values['monovar'][i] == 1:
                counter['monovar'][i] += 1
            elif values['sciphi'][i] == 1:
                counter['sciphi'][i] += 1
            elif values['sciphin'][i] == 1:
                counter['sciphin'][i] += 1

print('\ta1\ta2\ta3\ta4\ta5\ta6\ta7\ta8\th1\th2\th3\th4\th5\th6\th7\th8')
for mix in counter:
    string_ints = [str(int) for int in counter[mix]]
    print(mix + '\t' + '\t'.join(string_ints))
