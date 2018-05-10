File Description
----------------

A mass spectrometry data file contains heterogenous types of information
derived from one or more provenance sources. Some formats, like `mzML
<https://doi.org/10.1074/mcp.R110.000133>`_, track this
information. This information can be queried to make the object select
behavior appropriate to its contents and to provide a clear chai of
evidence back to the original raw data.


.. automodule:: ms_deisotope.data_source.metadata.file_information

    .. autoclass:: FileInformation
        :members:

    .. autoclass:: SourceFile
        :members:

    .. data:: id_formats
        
        These are the recognized formats for encoding scan identifiers from raw and
        processed mass spectrometry data files, derived from the :title-reference:`HUPO PSI-MS`
        controlled vocabulary.

        .. exec::

            from ms_deisotope.data_source.metadata.file_information import id_formats
            from rst_table import as_rest_table

            def xsd_detection(line, directive='``%s``'):
                import re
                xsd_pattern = re.compile(r"(\S+=xsd:\S+)")
                matches = list(xsd_pattern.finditer(line))
                parts = []
                start = 0
                for match in matches:
                    parts.append(line[start:match.start()])
                    start = match.start()
                    parts.append(directive % line[start:match.end()])
                    start = match.end()
                parts.append(line[start:])
                return ''.join(parts)

            data = sorted(id_formats)
            rows = [("Format Name", "ID", ' ', "Description")]
            for desc in data:
                rows.append((desc.name, desc.id, ' ', xsd_detection(desc.description)))
            print(as_rest_table(rows))

    .. data:: file_formats

        These are the recognized formats for storing raw and processed mass spectrometry
        data in, derived from the :title-reference:`HUPO PSI-MS` controlled vocabulary.

        .. exec::

            from ms_deisotope.data_source.metadata.file_information import file_formats
            from rst_table import as_rest_table

            rows = [('Format Name', ' ', 'ID')]

            for file_format in sorted(file_formats):
                rows.append((file_format.name, ' ', file_format.id))

            print(as_rest_table(rows))

    .. data:: content_keys

        These are commonly used to describe the contents of a mass spectrometry data file,
        derived from the :title-reference:`HUPO PSI-MS` controlled vocabulary.

        .. exec::

            from ms_deisotope.data_source.metadata.file_information import content_keys
            from rst_table import as_rest_table

            rows = [('Content Type', ' ', 'ID')]

            for content in sorted(content_keys):
                rows.append((content.name, ' ', content.id))

            print(as_rest_table(rows))
