import click
import ms_deisotope
from ms_deisotope import output

from ms_deisotope.output import MzMLSerializer
from ms_deisotope.data_source import query

from ms_deisotope.tools.utils import register_debug_hook


@click.command('faims-split')
@click.argument("source_file", type=click.Path(exists=True, readable=True))
@click.argument("output_prefix")
def main(source_file, output_prefix):
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
        click.echo("Closing %r" % (sink, ))
        sink.close()

if __name__ == "__main__":
    main.main()



