import subprocess
from pathlib import *
import matplotlib.pyplot as plt
import os
from xlwt import Workbook

def saveExperimentsToExelTable(sheet1, listSizes):
    sheet1.write(0, 0, 'Size of matrix')

    sheet1.write(0, 1, 'Simple version')
    sheet1.write(1, 1, 'Time of decomposition')
    sheet1.write(1, 2, 'Time of calculation Q')
    sheet1.write(1, 3, '||A - Q * R||')

    sheet1.write(0, 4, 'RowHouse version')
    sheet1.write(1, 4, 'Time of decomposition')
    sheet1.write(1, 5, 'Time of calculation Q')
    sheet1.write(1, 6, '||A - Q * R||')

    sheet1.write(0, 7, 'Reverse accumulation version')
    sheet1.write(1, 7, 'Time of decomposition')
    sheet1.write(1, 8, 'Time of calculation Q')
    sheet1.write(1, 9, '||A - Q * R||')

    sheet1.write(0, 10, 'Final version')
    sheet1.write(1, 10, 'Time of decomposition')
    sheet1.write(1, 11, 'Time of calculation Q')
    sheet1.write(1, 12, '||A - Q * R||')

    for i in range(2, len(listSizes) + 2, 1):
        sheet1.write(i, 0, listSize[i - 2])


def saveExperimentForOneMode(listTimeDecompositionForOneMode, listTimeSelectionQForOneMode, listAccuracyForOneMode):
    a = 1

def savePng(name='', fmt='png'):
    pwd = os.getcwd()
    iPath = './pictures/{}'.format(fmt)
    if not os.path.exists(iPath):
        os.mkdir(iPath)
    os.chdir(iPath)
    plt.savefig('{}.{}'.format(name, fmt), fmt='png')
    os.chdir(pwd)


if __name__ == "__main__":
    path = PurePath(r'C:\Users\Julia\source\repos\Householder\x64\Release\Householder.exe')
    listSize = [3, 5, 10, 50, 100, 300, 500, 800, 1000, 1500, 2000, 2500, 3000]
    listTimeDecompositionForOneMode = []
    listTimeSelectionQForOneMode = []
    listAccuracyForOneMode = []

    # Workbook is created
    wb = Workbook()

    # add_sheet is used to create sheet.
    sheet1 = wb.add_sheet('May 2020 - MS')
    saveExperimentsToExelTable(sheet1, listSize)

    for mode in range(1, 5, 1):
        for size in listSize:
            if mode < 3 and size > 1000:
                continue
            print("---------------------------------------------------------------- \n")
            print(f"mode = {mode}; size = {size} \n")
            command = f"{mode} {size}"
            output = subprocess.check_output(f"{path} {command}",shell=True, universal_newlines=True)
            _, outputDecomposition = output.split("Time for simple version of decomposition: ")
            outputDecomposition, outputSelectionQ = outputDecomposition.split("Time for simple version of Q: ")

            print(f"Time for decomposition = {outputDecomposition}")
            listTimeDecompositionForOneMode.append(float(outputDecomposition))

            outputSelectionQ, outputMaxAbs = outputSelectionQ.split("Check result matrices...")
            print(f"Time for selection = {outputSelectionQ}")
            listTimeSelectionQForOneMode.append(float(outputSelectionQ))

            _, outputMaxAbs = outputMaxAbs.split("max abs = ")
            outputMaxAbs, _ = outputMaxAbs.split(" OK", maxsplit=1)
            print(f"Accuracy = {outputMaxAbs}")
            listAccuracyForOneMode.append(float(outputMaxAbs))

        saveExperimentForOneMode(listTimeDecompositionForOneMode, listTimeSelectionQForOneMode, listAccuracyForOneMode)



    wb.save('xlwt experiments.xls')
    '''   fig = plt.figure()
    grid1 = plt.grid(True)
    graph1 = plt.plot(listSize, listTimeDecomposition, label="Insert in avl tree", color='red')


    #save("Insert", 'png')
    plt.show()
    
    '''