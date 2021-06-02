import sys

UMI_format = sys.argv[1]

for line in sys.stdin:
    line = line.replace('\n', '')
    if line[0] == '@':
        print(line)
        #GATCCCCTGC:K00180:234:HCTHMBBXX:1:1106:29244:17315 99  chr1    13199   255 32M =   13308   138 TGGATACCACTTTCCCACGAAGGCAGGGCCAT    JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ    NH:i:1  HI:i:1  AS:i:41 nM:i:9  NM:i:9  MD:Z:0C0T0T1G0C1T0C0T0G21   jM:B:c,-1   jI:B:i,-1   RG:Z:foo
    else:
        line = line.split('\t')
        UMI = line[0].split(':')[0]
        OG_read_name = ':'.join(line[0].split(':')[1:])
        if UMI_format == 'pureclip':
            new_name = OG_read_name + '_' + UMI
        elif UMI_format == 'iCount':
            new_name = OG_read_name + ':' + UMI
        line[0] = new_name
        print('\t'.join(line))

