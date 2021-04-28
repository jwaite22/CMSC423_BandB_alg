file = open("input", "r")
string = file.read()
data = string.split(" ")
spectrum = []
for i in data:
    spectrum.append(int(i))
spectrum = spectrum[1:]


AA_Table =[57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186]


def can_table(AA_Table,spectrum):
    cans = []
    for curr in AA_Table:
        if curr in spectrum:
            cans.append([curr])
    return cans


def table(AA_Table,spectrum):
    cans = []
    for curr in AA_Table:
        if curr in spectrum:
            cans.append(curr)
    return cans


can_peptides = can_table(AA_Table, spectrum)
table = table(AA_Table, spectrum)
fin_peptides = []


def Expand(can_peptides):
    pep = can_peptides[:]
    for i in can_peptides:
        for amino in table:
            new = i[:]  
            new.append(amino)
            pep.append(new)
        pep.remove(i)
    return pep

test = Expand(can_peptides)

def Mass(peptide):
    return sum(peptide)



def ParentMass(spectrum):
    return spectrum[-1]


def lin_spectrum(peptide):
    final = []
    plst = peptide
    if peptide[0] == 0:
        plst = plst[1:]
        final.append(0)
    for i, A in enumerate(peptide):
        
        mass = 0
        
        for x in range(i, len(plst)):
            mass += plst[x]
            print(mass)
            final.append(mass)
    final.sort()
    return final




def consistent(cyclo, spectrum):
    spec = spectrum[:]
    for curr in cyclo:
        if curr not in spec:
            return False
        spec.remove(curr)
    return True


def Cyclospectrum(peptide):
    final = []
    plst = peptide + peptide
    for i, A in enumerate(peptide):
        for j in range(i,i+len(peptide)):
            mass = 0
            for x in range(i,j+1):
                mass += plst[x]
            final.append(mass)
    final.sort()
    final = final[:(-len(peptide)+1)]
    return final


while can_peptides:
    can_peptides = Expand(can_peptides)
    test_peptides = can_peptides[:]
    for pep in can_peptides:
        cyc = Cyclospectrum(pep)
        if Mass(pep) == ParentMass(spectrum):
            if ( cyc == spectrum and (pep not in fin_peptides)):
                fin_peptides.append(pep)
            test_peptides.remove(pep)
        elif not consistent(lin_spectrum(pep),spectrum): 
            test_peptides.remove(pep)
    can_peptides = test_peptides



outfile= open("output","w+")

for j, pep in enumerate(fin_peptides):
    for i,nuc in enumerate(pep):
        if  i < len(pep)-1:
            outfile.write(str(nuc) + "-")
        elif i == len(pep)-1 and j < len(fin_peptides):
            outfile.write(str(nuc) + " ")
        else:
            outfile.write(str(nuc))
outfile.close()
