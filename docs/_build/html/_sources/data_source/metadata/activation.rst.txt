Activation
----------


.. automodule:: ms_deisotope.data_source.metadata.activation

    .. autoclass:: ActivationInformation
        :members:

    .. autoclass:: MultipleActivationInformation
        :members:

    .. autoclass:: DissociationMethod
        :members:

    .. data:: dissociation_methods
        
        These are the recognized dissociation method names derived from
        the :title-reference:`HUPO PSI-MS` controlled vocabulary.

        .. exec::

            from ms_deisotope.data_source.metadata.activation import dissociation_methods_map
            from rst_table import as_rest_table

            rows = [("Dissociation Method", "ID", "Category")]
            for term in sorted(set(dissociation_methods_map.values())):
                rows.append((term.name, term.id, term.specialization[0]))

            print(as_rest_table(rows))
