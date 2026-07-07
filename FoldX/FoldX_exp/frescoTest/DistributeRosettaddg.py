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
    print('This is a python implementation of the DistributeRosettaddg script written by Hein Wijma.\n'
          'By setting LEGACY to True inside the script, it wil use the original input structure.\n'
          'Legacy mode is currently {}.\n'. format({True:'ON', False: 'OFF'}[LEGACY]))
    if legacymode:
        print('This program is intended to distribute the calculation of a large number of mutations by Rosetta.\n'
              'For the first phase of preparing the files for calculations, use a command like:\n\n'
              '\tpython DistributeRosettaddg.py Phase1 DimerSelection.tab 2 A 4 B 149 DimerParsedForRosetta.pdb 100 '
              'FLAGrow3 /home/wijma/mini20101104/mini/bin/fix_bb_monomer_ddg.linuxgccrelease \n\n'
              'This should PREPARE mutations of DimerParsedForRosetta.pdb in the 2 subunits A and B '
              'of the residues in the File MyPreciousSelection.tab \nand distribute them over directories with each '
              'a 100 different mutations, in the pdb file the subunits start at residue 4 and 149. \n FLAGrow3 is the '
              'FLAG file to be used. \n/home/wijma/mini/bin/fix_bb_monomer_ddg.linuxgccrelease '
              'is the name and location of the ddg software.\n'
              
              'For the second phase of collecting the results into a list, use a command like:\n\n'
              '\tpython DistributeRosettaddg.py Phase2 MyPrecious.pdb -5 \n\n'
              'This should COLLECT the mutations of MyPrecious.pdb assuming a cutoff of -5 KJ mol-1,'
              'resulting in FOUR lists:\n'
              '- a complete list of all mutations\n'
              '- a list of the mutations that are less than -5 kJ mol-1\n'
              '- a list with the best mutation per position\n'
              '- a list with the best mutation per position that are less than -5 kJ mol-1\n\n'
              )
    else:
        print('This program is intended to distribute the calculation of a large number of mutations by FoldX.\n'
              'For the first phase of preparing the files for calculations, use a command like:\n\n'
              '\tpython DistributeRosettaddg.py Phase1 DimerParsedForRosetta.pdb DimerSelection.tab 100 '
              'FLAGrow3 /home/wijma/mini20101104/mini/bin/fix_bb_monomer_ddg.linuxgccrelease \n\n'
              'This should PREPARE mutations of DimerParsedForRosetta.pdb of the residues in the File '
              'MyPreciousSelection.tab and \ndistribute them over directories with each a 100 different mutations. '
              'FLAGrow3 is the FLAG file to be used. \n/home/wijma/mini/bin/fix_bb_monomer_ddg.linuxgccrelease '
              'is the name and location of the ddg software.\n'
              
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


def ReadStrandStart(pdbname, strands):
    '''
    given a path to an existing .pdb and a list of strands,
    will read the .pdb and find the first residue number for each strand
    '''
    with open(pdbname) as pdb:
        pdbcontents = pdb.read().splitlines()
    # read only lines starting with ATOM. since its a fixed with format, directly index strand and resnr
    pdbcontents = np.array([[line[21], line[22:26]] for line in pdbcontents if line[0:4]=='ATOM'])
    # seperate strand and resnr colums
    pdb_strand = pdbcontents[:,0]
    pdb_resnr = pdbcontents[:,1].astype(int)
    # get lowest resnr for each strand
    strandstarts = []
    for strand in strands:
        strand_resnr = pdb_resnr[np.where(pdb_strand == strand)]
        if len(strand_resnr) >= 1:
            strandstarts.append(min(strand_resnr))
        else:
            strandstarts.append(np.NaN)
    return strandstarts


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
        strandname1, strandname2 = sequence[selection1, 1][0], sequence[selection2, 1][0]
        print('Checking Strands: {} and {}'.format(strandname1, strandname2))
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
    for i in sequence:
        # copy residue 19 times and reorder to make sure output will be [aminoacid, resnr, mutation]
        muts = np.repeat(i[np.newaxis, ...], 19, axis=0)[:,[0, 2, 1]]
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
    # find indices of the best mutation per position
    currentres = ''
    currentbest = ''
    bestperpos_out = []
    for ind, (mut, kj) in enumerate(zip(mutationlist, energylist)):
        # only take relevant part
        mut = mut[0:-1]
        # if this is new residue or last entry, save index of current best
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
    ''' 
    Turn a string into a list of string and number chunks,
    but skips first value.
        "z23a" -> [23, "a"]
    '''
    return [ tryint(c) for c in re.split('([0-9]+)', s) ][1:]

# ====================check and assign all arguments in the command line to variables===================================
# ======================================================================================================================
# print example use
print('\n')
ExampleUse(LEGACY)

# some sanity checks before assigning to variables:
CheckError(len(sys.argv) == 1, 'Arguments are missing, see example above')
CheckError(sys.argv[1] not in ['Phase1', 'Phase2'], 'Phase has to be Phase1 or Phase2: '+sys.argv[1])
WhichPhaseAreWeIn = sys.argv[1]

# more sanity checks
if LEGACY:
    CheckError((len(sys.argv) < 4), 'Arguments are missing, see example above')
    try:
        int(sys.argv[3])
    except:
        CheckError(True, 'Given number of subunits is not an integer')
    CheckError((len(sys.argv) < 8+2*int(sys.argv[3])), 'Arguments are missing, see example above')
    CheckError((len(sys.argv) > 8+2*int(sys.argv[3])), 'Too many arguments, see example above')
    try:
        int(sys.argv[5+2*int(sys.argv[3])])
    except:
        CheckError(True, 'Given number of mutations per directory is not an integer')
    try:
        [int(start) for start in sys.argv[5:4 + 2 * int(sys.argv[3]):2]]
    except:
        CheckError(True, 'One of the given start residues is not an integer')
elif WhichPhaseAreWeIn == 'Phase1':
    CheckError((len(sys.argv) > 7), 'Too many arguments, see example above')
    CheckError((len(sys.argv) < 7), 'Arguments are missing, see example above')
    try:
        int(sys.argv[4])
    except:
        CheckError(True, 'Given number of mutations per directory is not an integer: ' + sys.argv[4])
elif WhichPhaseAreWeIn == 'Phase2':
    CheckError((len(sys.argv) > 4), 'Too many arguments, see example above')
    CheckError((len(sys.argv) < 4), 'Arguments are missing, see example above')
    try:
        float(sys.argv[3])
    except:
        CheckError(True, 'Given energy cutoff is not a number')


# assign arguments to variables:
if LEGACY:
    NameTableFile = sys.argv[2]
    NumberOfStrands = int(sys.argv[3])
    Strands = sys.argv[4:4+2*NumberOfStrands:2]
    Strandstarts = [int(start) for start in sys.argv[5:4+2*NumberOfStrands:2]]
    NamePDBFile = sys.argv[4+2*NumberOfStrands]
    MutationsPerDirectory = int(sys.argv[5+2*NumberOfStrands])
    if WhichPhaseAreWeIn == 'Phase1':
        Flagname = sys.argv[6+2*NumberOfStrands]
        Rosettalocation = sys.argv[7+2*NumberOfStrands]
    if WhichPhaseAreWeIn == 'Phase2':
        EnergyCutoff = sys.argv[6+2*NumberOfStrands]
else:
    if WhichPhaseAreWeIn == 'Phase1':
        NamePDBFile = sys.argv[2]
        NameTableFile = sys.argv[3]
        MutationsPerDirectory = int(sys.argv[4])
        Flagname = sys.argv[5]
        Rosettalocation = sys.argv[6]
    if WhichPhaseAreWeIn == 'Phase2':
        NamePDBFile = sys.argv[2]
        EnergyCutoff = float(sys.argv[3])

# make sure the file extension on all files is correct and all files exist
NamePDBFile = CheckFileExtension(NamePDBFile, ['pdb', 'ent', 'brk'])
CheckError(not os.path.exists(NamePDBFile), '.pdb file does not exist: ' + NamePDBFile)

if WhichPhaseAreWeIn == 'Phase1':
    NameTableFile = CheckFileExtension(NameTableFile, ['tab'])
    CheckError(not os.path.exists(NameTableFile), '.tab file does not exist: ' + NameTableFile)
    CheckError(not os.path.exists(Rosettalocation), 'Rosetta executable does not exist: ' + Rosettalocation)
    CheckError(not os.path.exists(Flagname), 'Rosetta Flag file does not exist.')
elif LEGACY:
    NameTableFile = CheckFileExtension(NameTableFile, ['tab'])
    CheckError(not os.path.exists(NameTableFile), '.tab file does not exist: ' + NameTableFile)
    CheckError((int(NumberOfStrands) > 4), 'More than 4 subunits')
    CheckError((int(NumberOfStrands) != len(Strands)), 'Number of subunits does not correspond to subunit names')

# print info about variables
print('Now performing: ', WhichPhaseAreWeIn)
if DEBUG:
    print('PDB: ', NamePDBFile)
    if LEGACY:
        print('Tab file: ', NameTableFile)
        print('Number of strands: ', NumberOfStrands)
        print('Strands: ', Strands)
        print('Start of strands: ', Strandstarts)
        print('Files per directory: ', MutationsPerDirectory)
    if WhichPhaseAreWeIn == 'Phase1':
        print('Tab file: ', NameTableFile)
        print('Files per directory: ', MutationsPerDirectory)
        print('The name and location of the Rosetta script is: ', Rosettalocation)
    if WhichPhaseAreWeIn == 'Phase2':
        print('Mutation energy Cutoff: ', EnergyCutoff)
    print('')

# ============read .tab file and if it is a multimer, ensure the residues are found in every strand ====================
# ======================================================================================================================

if WhichPhaseAreWeIn == 'Phase1':
    # read table file into np array
    # this will return an str array with [resnr, strandname, aminoacid] for each residue
    RawSequence = ReadTabFile(NameTableFile)
    if not LEGACY:
        # get strands by finding all unique entries in 2nd column
        Strands = list(set(RawSequence[:, 1]))
        Strands.sort()
        Strandstarts = ReadStrandStart(NamePDBFile, Strands)
        NumberOfStrands = len(Strands)
        if DEBUG:
            print('not in LEGACY mode, so NumberofStrands Strands, and Strandstarts were found automatically:')
            print('Found ', NumberOfStrands, ' Strand(s): ', Strands, ' Starting at: ', Strandstarts)
        CheckError(np.isnan(np.array(Strandstarts)).any(),
                   'One or more of the strands in the .tab file do not exist in the .pdb')

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

    print('Finished eliminating residues that are not present in all of the {} strands\n'.format(NumberOfStrands))

# ===============  make list of mutations and write them to files that can be used by Rosetta ==========================
# ======================================================================================================================

if WhichPhaseAreWeIn == 'Phase1':
    #get list of all mutations and its length
    MutatedProteinList = SaturationScanning(PurgedSequence_One)
    NumberOfMutations = len(MutatedProteinList)
    print('The number of mutations is: ', NumberOfMutations)

    #decide on how many subdirectories to make
    NumberOfSubdirectories = int(NumberOfMutations/int(MutationsPerDirectory))
    NumberInLastDirectory = NumberOfMutations % MutationsPerDirectory
    if NumberInLastDirectory != 0:
        NumberOfSubdirectories +=1
    print('With maximum {} mutations per subdirectory this requires {} directories'.format(MutationsPerDirectory,
           NumberOfSubdirectories))

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
            #continue

        # create 'list.txt' file with name of original pdb
        with open(os.path.join(Subdirectory_name, 'list.txt'), "w") as pdblist:
            pdblist.write(NamePDBFile)

        # create 'RosettaFormatMutations.mut file with list of mutants
        with open(os.path.join(Subdirectory_name, 'RosettaFormatMutations.mut'), "w") as mutlist:
            # write total number of mutations to file
            if directory_nr == NumberOfSubdirectories:
                mutlist.write('total {}\n'.format(len(Strands)*NumberInLastDirectory))
            else:
                mutlist.write('total {}\n'.format(len(Strands)*MutationsPerDirectory))
            #for loop going over each mutation
            for mut in MutatedProteinList[StartMutationRange:EndMutationRange]:
                # subtract the start from the residue number to get new residue numbers
                newres = [str(int(mut[1])-int(start)+1) for start in Strandstarts]
                # for mutation 'K4E' with 3 subunits starting at 2, make [3, 'K 3 E', 'K 3 E', 'K 3 E']
                mutlines = [str(len(Strands))+'\n']
                mutlines += [' '.join([mut[0], res, mut[2], '\n']) for res in newres]
                mutlist.writelines(mutlines)

        # create 'List_Mutations_readable.txt' with list of mutants
        with open(os.path.join(Subdirectory_name, 'List_Mutations_readable.txt'), "w") as readlist:
            #for loop going over each mutation
            for mut in MutatedProteinList[StartMutationRange:EndMutationRange]:
                # subtract the start from the residue number to get new residue numbers
                newres = [str(int(mut[1])-int(start)+1) for start in Strandstarts]
                # for mutation 'K4E' with 3 subunits starting at 2, make 'K3EK3EK3E is K 4 E'
                newmuts = ''.join([mut[0]+start+mut[2] for start in newres])
                oldmut = ' '.join(mut)
                readlist.write('{} is {}\n'.format(newmuts, oldmut))

        # copy pdb file and rotabase.txt into new directory
        shutil.copy2(NamePDBFile, Subdirectory_name)
        shutil.copy2(Flagname, Subdirectory_name)
        # append new lines to the todolist file
        todolines.append('cd {}\n'.format(Subdirectory_name))
        todolines.append('{} @{} -in:file:s {} -ddg::mut_file RosettaFormatMutations.mut '
                         '>LOG&\n'.format(Rosettalocation, Flagname, NamePDBFile))
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
        RosettaEnergyFileLoc = os.path.join(directory_name, 'ddg_predictions.out')
        MutlistLoc = os.path.join(directory_name, 'List_Mutations_readable.txt')

        # check if neccesary files exists
        print("\nCurrently looking for {}....".format(RosettaEnergyFileLoc))
        if not os.path.exists(RosettaEnergyFileLoc):
            print('{} does not exist, cant read energies'.format(RosettaEnergyFileLoc))
            continue
        if not os.path.exists(MutlistLoc):
            print('{} does not exist, cant read mutants'.format(MutlistLoc))
            continue

        # read in energy file and mutation list
        with open(MutlistLoc, "r") as Mutfile:
            Mutlist = [line.rstrip('\n') for line in Mutfile]
        if len(Mutlist) == 0:
            print('{} is empty'.format(Mutfile))
        with open(RosettaEnergyFileLoc, "r") as Energyfile:
            Energylist = [line.rstrip('\n') for line in Energyfile]
        if len(Energylist) == 0:
            print('{} is empty'.format(Energyfile))
            continue

        # skip first line, get Energy entries from the energy file, and convert energy from kcal/mol to kj/mol
        Energylist = np.array([line.split()[2] for line in Energylist[1:] if 'ddG' in line])
        Energylist = Energylist.astype(float)*4.1840
        # format entries in mutation list from '1 K 4 E to 'K4E'
        Mutlist = np.array([''.join(line.split(' ')[-3:]) for line in Mutlist])
        # keep track of total amount of expected outputs
        Total_expected += len(Mutlist)

        print('energy', Energylist, 'muts', Mutlist)

        # give warning if the two are not the same length
        print('Found {} mutations'.format(len(Energylist)))
        if len(Mutlist) != len(Energylist):
            print('WARNING: the number of energies in ddg_predictions.out '
                  'do not match the number of mutations in {}:'.format(MutlistLoc))
            print('Only the mutations for which energy values have been calculated will be reported')
            print('Input number of mutations: {}'.format(len(Mutlist)))
            print('Output number of mutations: {}'.format(len(Energylist)))
        # keep only mutation list entries with matching energy output, append both to full list
        Mutlist = Mutlist[0:len(Energylist)]
        Mutlist_full = np.append(Mutlist_full, Mutlist)
        Energylist_full = np.append(Energylist_full, Energylist)
        SDlist_full = np.append(SDlist_full, np.zeros(len(Energylist)))

    print('\n')
    # find indices where energy is below cutoff, and find indices of best per position
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
          "................................................... {}".format(len(Completelist)))
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
