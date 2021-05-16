import sys
import csv
import click

import ms_deisotope

@click.command("collision-energy-extractor")
@click.argument("infile", type=click.Path(readable=True))
@click.argument("outfile", type=click.Path(writable=True))
def main(infile, outfile):
    '''Read the real HCD collision energy from a Thermo RAW file for each
    scan.
    '''
    reader = ms_deisotope.MSFileLoader(infile)
    columns = ['scan_id', 'ms_level', 'charge', 'precursor_mz', 'activation_name',
               'energy', 'Thermo_Trailer_Extra_HCD_Energy_eV',
               'Thermo_Trailer_Extra_HCD_Energy']

    if sys.version_info.major == 3:
        stream = open(outfile, 'wt', newline='')
    else:
        stream = open(outfile, 'wb')
    with stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(columns)
        for scan in reader.make_iterator(grouped=False):
            if scan.ms_level == 1:
                continue
            activation = scan.activation
            annotations = scan.annotations
            fields = [
                scan.id, scan.ms_level,
                scan.precursor_information.charge,
                scan.precursor_information.mz,
                activation.method, activation.energy,
                annotations.get('[Thermo Trailer Extra]HCD Energy eV'),
                annotations.get('[Thermo Trailer Extra]HCD Energy'),
            ]
            writer.writerow(fields)


if __name__ == "__main__":
    main.main()
