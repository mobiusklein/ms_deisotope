#!/usr/bin/env python
import click
import ms_deisotope

from ms_deisotope.output import MzMLSerializer
from ms_deisotope.data_source import query


@click.command('faims-split')
@click.argument("source_file", type=click.Path(exists=True, readable=True))
@click.argument("output_prefix")
def main(source_file, output_prefix):
    """Read in `source_file` and split it based upon FAIMS compensation voltage into separate
    mzML files whose path prefix matches `output_prefix` and ends with the compensation voltage
    dedicated to that stream.
    """
    reader = ms_deisotope.MSFileLoader(source_file)

    sinks = {}
    n = len(reader)
    i = 0
    last = 0
    interval = 1000
    demultiplexer = query.FAIMSDemultiplexingIterator(reader, 30)

    for heads in demultiplexer:
        empty_channels = 0
        for channel, value in heads.items():
            if value is not None:
                if channel not in sinks:
                    click.echo("Opening new channel for CV %r" % (channel, ))
                    writer = MzMLSerializer(
                        open(output_prefix + "_%r.mzML" % (channel, ), 'wb'),
                        len(reader), deconvoluted=False)
                    writer.copy_metadata_from(reader)
                    method = writer.build_processing_method(1, False, False, False, ["ion mobility seperation"])
                    writer.add_data_processing(method)
                    sinks[channel] = writer
                else:
                    writer = sinks[channel]
                writer.save_scan(value)
                i += 1
            else:
                empty_channels += 1
        if empty_channels == len(heads):
            click.echo("All channels empty, finishing")
            break
        if i - last >= interval:
            click.echo("Processed %d spectra (%0.3f%%)" % (i, i * 100.0 / n))
            last = i

    click.echo("Closing buffers.")
    for sink in sinks.values():
        click.echo("Closing %r" % (sink.handle.name, ))
        sink.close()

if __name__ == "__main__":
    main.main()



