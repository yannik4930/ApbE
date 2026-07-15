# standard libraries
import sys
import os
import re
import shutil
# requires package installment
import numpy as np

# if set to True, maintains exactly the same input structure as the original FRESCO scripts
# if set to False, strand names and number of strands do not need to be specified
LEGACY = False
# provides additional information about each step if set to True
DEBUG = False

# ============================================Define functions =========================================================
# ======================================================================================================================


def ExampleUse(legacymode):
    print('This is a python implementation of the DistributeFoldx script written by Hein Wijma.\n'
          'The script requires less inputs than the original script. by setting LEGACY to True inside the script,\n'
          'it wil instead use the original input structure. '
          'Legacy mode is currently {}.'. format({True:'ON', False: 'OFF'}[LEGACY]))
    if legacymode:
        print('This program is intended to distribute the calculation of a large number of mutations by FoldX.\n'
              'For the first phase of preparing the files for calculations, use a command like:\n\n'
              '\tpython DistributeFoldX.py Phase1 MyPrecious.pdb 2 A B MyPreciousSelection.tab '
              '100 /home/wijma/Fold_X/FoldX.linux64 \n\n'
              'This should PREPARE mutations of MyPrecious.pdb in the 2 '
              'subunits A and B of the residues in the File MyPreciousSelection.tab and\n'
              'distribute them over directories with each a 100 different mutations. '
              
              'For the second phase of collecting the results into a list, use a command like:\n\n'
              '\tpython DistributeFoldX.py Phase2 MyPrecious 2 A B MyPreciousSelection.tab 100 -5 \n\n'
              'This should COLLECT the mutations of MyPrecious.pdb in the 2 subunits A and B '
              'of the residues in the File MyPreciousSelection.tab and make FOUR lists:\n'
              '- a complete list of all mutations\n'
              '- a list of the mutations that are less than -5 kJ mol-1\n'
              '- a list with the best mutation per position\n'
              '- a list with the best mutation per position that are less than -5 kJ mol-1\n\n'
              )
    else:
        print('This program is intended to distribute the calculation of a large number of mutations by FoldX.\n'
              'For the first phase of preparing the files for calculations, use a command like:\n\n'
              '\tpython DistributeFoldX.py Phase1 MyPrecious.pdb MyPreciousSelection.tab '
              '100 /home/wijma/Fold_X/FoldX.linux64 \n\n'
              'This should PREPARE mutations of MyPrecious.pdb for residues specified '
              'in the File MyPreciousSelection.tab and\n'
              'distribute them over directories with each a 100 different mutations. '
              
              'For the second phase of collecting the results into a list, use a command like:\n\n'
              '\tpython DistributeFoldX.py Phase2 MyPrecious.pdb -5 \n\n'
              'This should COLLECT the mutations of MyPrecious.pdb assuming a cutoff of -5 KJ mol-1,'
              'resulting in FOUR lists:\n'
              '- a complete list of all mutations\n'
              '- a list of the mutations that are less than -5 kJ mol-1\n'
              '- a list with the best mutation per position\n'
              '- a list with the best mutation per position that are less than -5 kJ mol-1\n\n'
              )


def CheckError(conditional, errormessage):
    '''
    Takes in conditional statement (Bool) and error message (str).
    Prints error message and quits if conditional is True
    '''
    # print message and end script if condition is met
    if conditional:
        print('-------- ERROR -------- '+errormessage)
        sys.exit()


def CheckFileExtension(name, extension):
    '''
    takes in a filename (str) and list of expected extensions (list),
    adds correct extension if filename misses extension,
    and checks for correct extension if the filename has an extension
    '''
    # only split the last occurrence of '.', i.e the file extension
    name = name.rsplit('.', 1)
    # if no extension, use first in list of allowed extensions
    if len(name) == 1:
        return '.'.join([name[0],extension[0]])
    # if allowed extension, use that extension
    elif len(name) == 2 and name[1] in extension:
        return '.'.join(name)
    # if no allowed extension, give error
    else:
        CheckError(True, '.{} file {} does not have any of the allowed '
                         'extensions {}'.format(extension[0], '.'.join(name), extension))


def ReadTabFile(tabfile):
    '''
    reads specified .tab file (str),
    removes blank first line and 'END' on last line,
    and returns all data in file as a numpy array
    '''
    with open(tabfile) as tab:
        # skip first line to get all data
        TabContents = tab.read().splitlines()
        # check if it starts with empty line and ends with END
        CheckError((TabContents[0] != ''), '.tab file did not start with empty line')
        CheckError((not 'END' in TabContents), '.tab file did not terminate with END')
    # shorten to just the sequence, skipping first line and stopping at first END
    TabContents = TabContents[1:TabContents.index('END')]
    # split entries to nested list, then turn into 2d array
    return np.array([res.strip().split(' ') for res in TabContents])

def PurgeSequence(sequence, selection1, selection2):
    '''
    Takes in a sequence read from .tab using ReadTabFile(), (np.array)
    and two selections from this sequence given as an array of indices. (np.array)
    Removes all differences between the two selections from the sequence and returns purged sequence
    '''
    # Combine first column (res nr) and last column (amino acid) to compare both at the same time
    # i.e. for strand B, [22, B, V,] becomes '22V',
    strand1 = np.char.add(sequence[selection1, 0], sequence[selection1, 2])
    strand2 = np.char.add(sequence[selection2, 0], sequence[selection2, 2])
    # get indices of residues in strand1 but not in strand2 and vice versa
    delsel1 = selection1[~np.isin(strand1, strand2)]
    delsel2 = selection2[~np.isin(strand2, strand1)]
    # remove all rows with selected indices
    purgedsequence = np.delete(sequence, np.append(delsel1, delsel2), axis=0)
    if DEBUG:
        print('Purged Residues:\n', sequence[np.append(delsel1, delsel2), :])
    return purgedsequence

def SaturationScanning(sequence):
    '''
    takes in a sequence from a .tab file (np.array),
    and does saturation scanning mutagenesis,
    outputs an array with all possible mutations for that sequence
    '''
    residues = np.array([*'ADEFGHIKLMNPQRSTVWY'])
    satmuts = np.array([])
    for res in sequence:
        # copy residue 19 times and reorder to make sure output will be [aminoacid, resnr, mutation]
        muts = np.repeat(res[np.newaxis, ...], 19, axis=0)[:,[0, 2, 1]]
        # replace strandname with mutation
        muts[:, 2] = residues
        # add mutations to final array
        if len(satmuts) == 0:
            satmuts = muts
        else:
            satmuts = np.vstack((satmuts, muts))
    return satmuts

def FindBestPerPos(mutationlist, energylist):
    '''
    takes in a list of mutations and a corresponding list of energies
    and returns a list of indices with the best mutation per residue
    '''
    # find indices of best mutation per position
    currentres = ''
    currentbest = ''
    bestperpos_out = []
    for ind, (mut, kj) in enumerate(zip(mutationlist, energylist)):
        # only take relevant part
        mut = mut[0:-1]
        # if this is new residue or last entry, save indice of current best
        if mut != currentres or ind == len(energylist) - 1:
            # skip if no best residue has yet been saved
            if not (currentres == '' or currentbest == ''):
                bestperpos_out.append(currentbest[0])
            currentres = mut
            currentbest = [ind, kj]
        # if this is the same residue, check if new mutation is better than best
        elif kj < currentbest[-1]:
            currentbest = [ind, kj]
    # return list of best per residue
    return bestperpos_out

def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    """ 
    Turn a string into a list of string and number chunks,
    but skips first value.
        "z23a" -> [23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ][1:]

# ====================check and assign all arguments in the command line to variables===================================
# ======================================================================================================================
# print example use
ExampleUse(LEGACY)

# some sanity checks before assigning to variables:
CheckError(len(sys.argv) == 1, 'Arguments are missing, see example above')
CheckError(sys.argv[1] not in ['Phase1', 'Phase2'], 'Phase has to be Phase1 or Phase2')
WhichPhaseAreWeIn = sys.argv[1]

# more sanity checks
if LEGACY:
    CheckError((len(sys.argv) < 3), 'Arguments are missing, see example above')
    try:
        int(sys.argv[3])
    except:
        CheckError(True, 'Given number of subunits is not an integer')
    CheckError((len(sys.argv) < 7+int(sys.argv[3])), 'Arguments are missing, see example above')
    CheckError((len(sys.argv) > 7+int(sys.argv[3])), 'Too many arguments, see example above')
    try:
        int(sys.argv[5+int(sys.argv[3])])
    except:
        CheckError(True, 'Given number of mutations per directory is not an integer')
elif WhichPhaseAreWeIn == 'Phase1':
    CheckError((len(sys.argv) < 6), 'Arguments are missing, see example above')
    CheckError((len(sys.argv) > 6), 'Too many arguments, see example above')
    try:
        int(sys.argv[4])
    except:
        CheckError(True, 'Given number of mutations per directory is not an integer')
elif WhichPhaseAreWeIn == 'Phase2':
    CheckError((len(sys.argv) > 4), 'Too many arguments, see example above')
    CheckError((len(sys.argv) < 4), 'Arguments are missing, see example above')
    try:
        float(sys.argv[3])
    except:
        CheckError(True, 'Given energy cutoff is not a number')


# assign arguments to variables:

if LEGACY:
    NamePDBFile = sys.argv[2]
    NumberOfStrands = int(sys.argv[3])
    Strands = sys.argv[4:4+NumberOfStrands]
    NameTableFile = sys.argv[4+NumberOfStrands]
    MutationsPerDirectory = int(sys.argv[5+NumberOfStrands])
    if WhichPhaseAreWeIn == 'Phase1':
        Foldxlocation = sys.argv[6+NumberOfStrands]
    if WhichPhaseAreWeIn == 'Phase2':
        EnergyCutoff = sys.argv[6+NumberOfStrands]
else:
    if WhichPhaseAreWeIn == 'Phase1':
        NamePDBFile = sys.argv[2]
        NameTableFile = sys.argv[3]
        MutationsPerDirectory = int(sys.argv[4])
        Foldxlocation = sys.argv[5]
    if WhichPhaseAreWeIn == 'Phase2':
        NamePDBFile = sys.argv[2]
        EnergyCutoff = float(sys.argv[3])

# make sure the file extension on all files is correct and all files exist
NamePDBFile = CheckFileExtension(NamePDBFile, ['pdb', 'ent', 'brk'])
CheckError(not os.path.exists(NamePDBFile), '.pdb file does not exist: ' + NamePDBFile)

if WhichPhaseAreWeIn == 'Phase1':
    NameTableFile = CheckFileExtension(NameTableFile, ['tab'])
    CheckError(not os.path.exists(NameTableFile), '.tab file does not exist: ' + NameTableFile)
    CheckError(not os.path.exists(Foldxlocation), 'Foldx executable does not exist: ' + Foldxlocation)
    #CheckError(not os.path.exists('rotabase.txt'), 'rotabase.txt file does not exist.')
elif LEGACY:
    NameTableFile = CheckFileExtension(NameTableFile, ['tab'])
    CheckError(not os.path.exists(NameTableFile), '.tab file does not exist: ' + NameTableFile)
    CheckError((int(NumberOfStrands) > 4), 'More than 4 subunits')
    CheckError((int(NumberOfStrands) != len(Strands)), '# of subunits does not correspond to # of subunit names')

# print info about variables
print('Now performing: ', WhichPhaseAreWeIn)
if DEBUG:
    print('PDB: ', NamePDBFile)
    if LEGACY:
        print('Tab file: ', NameTableFile)
        print('Number of strands: ', NumberOfStrands)
        print('Strands: ', Strands)
        print('Files per directory: ', MutationsPerDirectory)
if WhichPhaseAreWeIn == 'Phase1':
    print('Tab file: ', NameTableFile)
    print('Files per directory: ', MutationsPerDirectory)
    print('The name and location of FoldX is: ', Foldxlocation)
if WhichPhaseAreWeIn == 'Phase2':
    print('Mutation energy Cutoff: ', EnergyCutoff)

# ===========read .tab file and if it is a multimer, ensure the residues are found in every strand =====================
# ======================================================================================================================

if WhichPhaseAreWeIn == 'Phase1':
    # read table file into np array
    # this will return an str array with [resnr, strandname, aminoacid] for each residue
    RawSequence = ReadTabFile(NameTableFile)
    if not LEGACY:
        # get strands by finding all unique entries in 2nd column
        Strands = list(set(RawSequence[:, 1]))
        Strands.sort()
        NumberOfStrands = len(Strands)
        if DEBUG:
            print('not in LEGACY mode, so NumberofStrands and Strands were found automatically:')
            print('Found ', NumberOfStrands, ' Strand(s): ', Strands)

    # purge sequence from residues not shared by all strands
    PurgedSequence = RawSequence.copy()
    # go over each possible combination of the strands pairwise
    for strand_i in Strands:
        for strand_j in Strands:
            # get residues belonging to each strand as array of indices
            Strandind_i = np.where(PurgedSequence[:, 1] == strand_i)
            Strandind_j = np.where(PurgedSequence[:, 1] == strand_j)
            # purge residues that are not shared between the two strands
            PurgedSequence = PurgeSequence(PurgedSequence, Strandind_i[0], Strandind_j[0])
    # after purging differences, take just one of the strands
    PurgedSequence_One = PurgedSequence[np.where(PurgedSequence[:, 1] == Strands[0])]

    # print information about purging
    if DEBUG:
        print('before purging non-matching residues:', len(RawSequence))
        print('after purging non-matching residues:', len(PurgedSequence))
        print('purged ', len(RawSequence) - len(PurgedSequence), ' residues.')
        print('ending up with {} residues per strand'.format(len(PurgedSequence_One)))
    print('Finished eliminating residues that are not present in all of the {} strands'.format(NumberOfStrands))

# ===============  make list of mutations and write them to files that can be used by FoldX ============================
# ======================================================================================================================

if WhichPhaseAreWeIn == 'Phase1':
    # get list of all mutations and its length
    MutatedProteinList = SaturationScanning(PurgedSequence_One)
    NumberOfMutations = len(MutatedProteinList)
    print('The number of mutations is: ', NumberOfMutations)

    # decide on how many subdirectories to make
    NumberOfSubdirectories = int(NumberOfMutations/int(MutationsPerDirectory))
    if NumberOfMutations % MutationsPerDirectory != 0:
        NumberOfSubdirectories +=1
    print('With maximum {} mutations per subdirectory this requires '
          '{} directories'.format(MutationsPerDirectory, NumberOfSubdirectories))

    # for loop over each subdirectory
    todolines = []
    print('Preparing directories...')
    for directory_nr in range(1, NumberOfSubdirectories+1):
        # define name and range of mutations
        Subdirectory_name = 'Subdirectory{}'.format(directory_nr)
        StartMutationRange = (directory_nr-1)*MutationsPerDirectory
        EndMutationRange = directory_nr*MutationsPerDirectory
        # create new directory if it does not exist
        try:
            os.mkdir(Subdirectory_name)
        except:
            print('{} already exists, skipping this directory'.format(Subdirectory_name))
            continue

        # create 'list.txt' file with name of original pdb
        with open(os.path.join(Subdirectory_name, 'list.txt'), "w") as pdblist:
            pdblist.write(NamePDBFile)
        # create 'individual_list.txt file with list of mutants
        with open(os.path.join(Subdirectory_name, 'individual_list.txt'), "w") as indivlist:
            for m in MutatedProteinList[StartMutationRange:EndMutationRange]:
                # line should be structured like "KA4E,KB4E,KC4E;" for mutation K4E on strands A,B,C
                allstrands = [m[0]+s+m[1]+m[2] for s in Strands]
                allstrands = ','.join(allstrands)+';\n'
                indivlist.write(allstrands)
        # create 'List_Mutations_readable.txt' with list of mutants
        with open(os.path.join(Subdirectory_name, 'List_Mutations_readable.txt'), "w") as readlist:
            for i, m in enumerate(MutatedProteinList[StartMutationRange:EndMutationRange]):
                # line should be structured like "3 K 4 E" for mutation K4E as the 3rd entry in file.
                readableline = str(i+1)+' '+' '.join(m)+'\n'
                readlist.write(readableline)

        # copy pdb file and rotabase.txt into new directory
        shutil.copy2(NamePDBFile, Subdirectory_name)
        # shutil.copy2('rotabase.txt', Subdirectory_name)
        # append new lines to the todolist file
        todolines.append('cd {}\n'.format(Subdirectory_name))
        todolines.append('{} --command=BuildModel --pdb={}  --mutant-file=individual_list.txt'
                         ' --numberOfRuns=5 > LOG&\n'.format(Foldxlocation, NamePDBFile))
        todolines.append('cd ..\n')

    # write todolist and make it executable
    with open('todolist', "w") as todo:
        todo.writelines(todolines)
    os.chmod("todolist", 0o777)

# =====Collect the data, analyse it and write it to the right files, write it in kJ mol instead of the kcal mol-1=======
# ======================================================================================================================
DirectoriesPresent = False
if WhichPhaseAreWeIn == 'Phase2':
    # get list of directories
    Directory_list = os.listdir('.')
    Directory_list = [i for i in Directory_list if re.match('Subdirectory[1-9]+', i)]
    Directory_list.sort(key=alphanum_key)
    print('Found {} Subdirectories'.format(len(Directory_list)))
    if len(Directory_list) > 0:
        DirectoriesPresent = True

if WhichPhaseAreWeIn == 'Phase2' and DirectoriesPresent:
    # initialize list with results
    Mutlist_full = np.array([])
    Energylist_full = np.array([])
    SDlist_full = np.array([])
    Total_expected = 0

    # for loop going over each directory
    for directory_name in Directory_list:
        # get filename and location
        NamePDBnoExtension = NamePDBFile.rsplit('.', 1)[0]
        FoldxEnergyFile = 'Average_{}.fxout'.format(NamePDBnoExtension)
        FoldxEnergyFileLoc = os.path.join(directory_name, FoldxEnergyFile)
        MutlistLoc = os.path.join(directory_name, 'List_Mutations_readable.txt')

        # check if neccesary files exists
        print("\nCurrently looking for {}....".format(FoldxEnergyFileLoc))
        if not os.path.exists(FoldxEnergyFileLoc):
            print('{} does not exist, cant read energies'.format(FoldxEnergyFileLoc))
            continue
        if not os.path.exists(MutlistLoc):
            print('{} does not exist, cant read mutants'.format(MutlistLoc))
            continue

        # read in energy file and mutation list
        with open(MutlistLoc, "r") as Mutfile:
            Mutlist = [line.rstrip('\n') for line in Mutfile]
        if len(Mutlist) == 0:
            print('{} is empty'.format(Mutfile))
        with open(FoldxEnergyFileLoc, "r") as Energyfile:
            Energylist = [line.rstrip('\n') for line in Energyfile]
        if len(Energylist) == 0:
            print('{} is empty'.format(Energyfile))
            continue

        # obtain just the SD and Energy entries from .fxout file, and convert energy from kcal/mol to kj/mol
        Energylist = np.array([line.split('\t', 3)[1:3] for line in Energylist if NamePDBnoExtension in line])
        Energylist = Energylist.astype(float)*4.1840
        # format entries in mutation list from '1 K 4 E to 'K4E'
        Mutlist = np.array([''.join(line.split(' ')[1:]) for line in Mutlist])
        # keep track of total amount of expected outputs
        Total_expected += len(Mutlist)

        # give warning if the two are not the same length
        print('Found {} mutations'.format(len(Energylist)))
        if len(Mutlist) != len(Energylist):
            print('WARNING: the number of energies in {} '
                  'do not match the number of mutations in {}:'.format(FoldxEnergyFile, MutlistLoc))
            print('Only the mutations for which energy values have been calculated will be reported')
            print('Input number of mutations: {}'.format(len(Mutlist)))
            print('Output number of mutations: {}'.format(len(Energylist)))

        # keep only mutation list entries with matching energy output, append both to full list
        Mutlist = Mutlist[0:len(Energylist)]
        Mutlist_full = np.append(Mutlist_full, Mutlist)
        Energylist_full = np.append(Energylist_full, Energylist[:, 1])
        SDlist_full = np.append(SDlist_full, Energylist[:, 0])

    print('\n')
    # find indices where energy is below cutoff, and find infices of best per position
    Completelist = np.where(Energylist_full == Energylist_full)[0]
    BelowCutoff = np.where(np.less(Energylist_full, EnergyCutoff))[0]
    BestPerPos = np.array(FindBestPerPos(Mutlist_full, Energylist_full))
    BestPerPosBelowCutoff = np.array([ind for ind in BestPerPos if np.isin(ind, BelowCutoff)])

    # go over each of the selections and create a sorted and unsorted .tab file
    for name, indices in zip(['CompleteList', 'BelowCutOff', 'BestPerPosition', 'BestPerPositionBelowCutOff'],
                             [Completelist, BelowCutoff, BestPerPos, BestPerPosBelowCutoff]):
        # create empty file and skip if indices are empty
        if len(indices) == 0:
            open('MutationEnergies_' + name + '.tab', "a").close()
            open('MutationEnergies_' + name + '_SortedByEnergy.tab', "a").close()
            continue
        # create array with mut name, energy, and SD
        fullarray = np.vstack((Mutlist_full[indices],
                               np.round(Energylist_full[indices], 3),
                               np.round(SDlist_full[indices], 3))).T
        # obtain sorted array using some numpy trickery
        fullarraysorted = fullarray[Energylist_full[indices].argsort()]
        # write both to file
        header = 'Below are the mutations and the change in stability and SD from 5 calculations, all in kJ mol -1'
        np.savetxt('MutationEnergies_'+name+'.tab', fullarray,
                   fmt='%s%9s%9s', header=header, comments='')
        np.savetxt('MutationEnergies_'+name+'_SortedByEnergy.tab', fullarraysorted,
                   fmt='%s%9s%9s', header=header, comments='')

    NumberOfMutationsCollectedBelowCutOff = len(BelowCutoff)
    NumberOfMutationsCollectedPerPosition = len(BestPerPos)
    NumberOfMutationsCollectedPerPositionBelowCutOff = len(BestPerPosBelowCutoff)
    print("\n")
    print("The total number of expected mutations from Foldx is"
          "..................................................... {}".format(Total_expected))
    print("The total number of mutations collected from FoldX is "
          "................................................... {}".format(len(Energylist_full)))
    print("Of those the following number of mutations passes the energy cut off"
          "..................................... {}".format(len(BelowCutoff)))
    print("Of those the following number of mutations was the best at that position"
          "................................. {}".format(len(BestPerPos)))
    print("Of those the following number of mutations was the best at that position "
          "and passed the energy cut off... {}".format(len(BestPerPosBelowCutoff)))
    print("The used energy cut off was ............................................"
          "................................. {} kJ mol-1".format(EnergyCutoff))
    print("\nThe results were written to files called:\n"
          " - MutationsEnergies_CompleteList.tab\n"
          " - MutationsEnergies_BelowCutOff.tab\n"
          " - MutationsEnergies_BestPerPosition.tab\n"
          " - MutationsEnergies_BestPerPositionBelowCutOff.tab\n")
    print("The same results sorted by best energy were written to files called:\n"
          " - MutationsEnergies_CompleteList_SortedByEnergy.tab\n"
          " - MutationsEnergies_BelowCutOff_SortedByEnergy.tab\n"
          " - MutationsEnergies_BestPerPosition_SortedByEnergy.tab\n"
          " - MutationsEnergies_BestPerPositionBelowCutOff_SortedByEnergy.tab\n")

print("\nFinished {} succesfully, exiting...\n".format(WhichPhaseAreWeIn))
