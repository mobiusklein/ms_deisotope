mzML
----

.. automodule:: ms_deisotope.data_source.mzml

    .. autoclass:: MzMLLoader
        :show-inheritance:

        Scan Access and Iteration Control
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        .. automethod:: get_scan_by_id
        .. automethod:: get_scan_by_index
        .. automethod:: get_scan_by_time

        .. automethod:: __getitem__
        .. automethod:: __len__

        .. automethod:: start_from_scan
        .. automethod:: make_iterator
        .. automethod:: __iter__

