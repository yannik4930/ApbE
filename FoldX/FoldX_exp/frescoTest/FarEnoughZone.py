import MDAnalysis

def get_chainIDs(selection):
    '''
    get a single chainID for each residue in the selection
    '''
    chainids = []
    for res in selection.residues.chainIDs:
        if len(set(res)) != 1:
            raise Error
        else:
            chainid = list(set(res))[0]
        chainids.append(chainid)
    return chainids

def switch_AA_code(input_string, mode='to_one'):
    '''
    change amino acid letter code from one to three letter or vice versa
    use mode=to_one to go from three letter code to one letter code
    use mode=to_three to go from one letter code to three letter code
    '''
    switch_AA = {'R':'ARG', 'H':'HIS', 'K':'LYS', 'D':'ASP', 'E':'GLU',
                 'S':'SER', 'T':'THR', 'N':'ASN', 'Q':'GLN', 'C':'CYS',
                 'G':'GLY', 'P':'PRO', 'A':'ALA', 'V':'VAL', 'I':'ILE',
                 'L':'LEU', 'M':'MET', 'F':'PHE', 'Y':'TYR', 'W':'TRP'}
    rswitch_AA = {v: k for k, v in switch_AA.items()}
    if mode == 'to_three':
        return switch_AA[input_string]
    if mode == 'to_one':
        return rswitch_AA[input_string]
    else:
        return input_string


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog='FarEnoughZone',
                    description='create the input list for the DistributeFoldX/DistributeRosetta scripts using MDAnalysis')
    parser.add_argument('-p', '--pdb', required=True, type=str, help='the input structure')
    parser.add_argument('-r', '--residue', default='', type=str, help='residue name around which residues are excluded from the selection')
    parser.add_argument('-d', '--distance',default=5, help='distance around which residues are excluded from the selection')
    parser.add_argument('-c', '--customsel', default='', type=str, help='instead of distance and resname, use a custom selection in MDAnalysis selection language' )
    parser.add_argument('-o', '--output', default='output.tab')
    args = parser.parse_args()
    inputpdb = args.pdb
    distance = args.distance
    selection= args.customsel
    output = args.output
    residue = args.residue

    # get the residues that fit criteria
    u = MDAnalysis.Universe(inputpdb)
    if len(selection) == 0 and len(residue) == 0:
        u_sel = u.select_atoms('protein')
        print('no selection specified, used all protein residues')
    elif len(selection) > 0:
        u_sel = u.select_atoms(selection)
        u_sel = u_sel.select_atoms('protein')
        print('selected {} residues using custom selection: {}'.format(len(u_sel.residues), selection))
    else:
        u_sel = u.select_atoms('(protein and not (byres around {} resname {}))'.format(distance, residue))
        print('selected {} residues with distance {} from {}'.format(len(u_sel.residues), distance, residue))

    # generate list
    list_out = ['']
    chainids = get_chainIDs(u_sel)
    for chain, res in zip(chainids, u_sel.residues):
        aa = switch_AA_code(res.resname, mode='to_one')
        resnr = res.resid
        list_out.append('{} {} {}'.format(aa, chain, resnr))
    list_out.append('END')

    with open(output, 'w') as f:
        f.writelines([i+'\n' for i in list_out])
