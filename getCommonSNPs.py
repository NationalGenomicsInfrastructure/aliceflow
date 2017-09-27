# - reading dbSNP rsXXX IDs from a file (one line each) and storing in a set
# - reading a VCF line, extract its rs ID
#   - if rsXXX is in the previous set, write out line

import click


@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--rsids',    '-r', type=str, help='List of rsXXX IDs from dbSNP, one a line ', required=True)
@click.option('--vcf',      '-v', type=str, help='VCF to extract from', required=True)
def selectSNPs(rsids,vcf):
    rsIDs = getIDs(rsids)
    printIDsOnly(vcf,rsIDs)

def getIDs(rsids):
    idSet = set()
    with open(rsids,'r') as rsidsFile:
        for line in rsidsFile:
            rsID = line.rstrip()
            idSet.add(rsID)
    return idSet

def printIDsOnly(vcf,rsIDs):
    with open(vcf,'r') as vcfFile:
        for line in vcfFile:
            line = line.rstrip()
            if line.startswith("#"):
                print line
            else:
                cols = line.split()
                if cols[2] in rsIDs:
                    print line

if __name__ == '__main__':
    selectSNPs()
