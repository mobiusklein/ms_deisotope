import click
import ms_deisotope

from ms_deisotope.output import ProcessedMzMLDeserializer, MzMLSerializer
from ms_deisotope.data_source.scan import ScanBunch, ScanBase
from ms_deisotope.data_source.scan.mobility_frame import IonMobilityFrame


class TimeOrderMergingIterator(object):
    def __init__(self, sources):
        self.sources = list(sources)
        self._init_heads()

    def _init_heads(self):
        n = len(self.sources)
        self.heads = [None] * n
        for i in range(n):
            self._advance(i)

    def __next__(self):
        result = self.get_next_head()
        if result is None:
            raise StopIteration()
        return result

    next = __next__

    def __iter__(self):
        return self

    def count_exhausted(self):
        n = 0
        for h in self.heads:
            n += h is None
        return n

    def _advance(self, i):
        try:
            self.heads[i] = next(self.sources[i])
        except StopIteration:
            self.heads[i] = None

    def get_next_head(self):
        best = None
        best_i = None
        best_time = float('inf')
        for i, head in enumerate(self.heads):
            if head is None:
                continue
            time = self.get_time(head)
            assert time is not None
            if time < best_time:
                best = head
                best_i = i
                best_time = time
        if best is None:
            return best
        i = best_i
        self._advance(i)
        return best

    def get_time(self, obj):
        if isinstance(obj, ScanBunch):
            if obj.precursor is not None:
                return self.get_time(obj.precursor)
            elif obj.products:
                return self.get_time(obj.products[0])
        elif isinstance(obj, ScanBase):
            return obj.scan_time
        elif isinstance(obj, IonMobilityFrame):
            return obj.time
        else:
            # Just guessing at this point
            return obj.scan_time


@click.command("merge-mzml")
@click.argument("source_paths", nargs=-1, required=True, type=click.Path(readable=True))
@click.argument("output_path", required=True, type=click.Path(writable=True))
def main(source_paths, output_path):
    '''Combine multiple processed mzML files together into a single file sorted by time.
    '''
    sources = []
    for source_path in source_paths:
        click.echo("Reading %r" % source_path)
        sources.append(ProcessedMzMLDeserializer(source_path))
    total_n = sum(map(len, sources))
    writer = MzMLSerializer(open(output_path, 'wb'), total_n)
    iterator = TimeOrderMergingIterator(sources)
    writer.copy_metadata_from(sources[0])
    with writer:
        i = 0
        for bunch in iterator:
            i += 1
            if i % 100 == 0:
                click.echo("Processed %d batches. %d sources depeted" % (i, iterator.count_exhausted()))
            writer.save(bunch)


if __name__ == "__main__":
    main.main()
