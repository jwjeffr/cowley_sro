.. _OVITO: https://www.youtube.com/watch?v=zbBx78u3gJI

A couple of things need to be specified in the script:

- A first cutoff
- A second cutoff
- A type map

To calculate the Cowley SRO parameter, we need to search for neighbors within some cutoff distance. The first and second cutoff are respectively the minimum and maximum distance to search for neighbors with. You can find this from the radial distribution function in `OVITO`_.

The type map (line 108) is optional, but specifying it labels the atom-atom pairs in the plot with their atom names. If type 1 means iron and type 2 means nickel in your script, you can write:

  type_map = {1: 'Fe', 2: 'Ni'}

on line 108. Then, run:

  python sro_param.py input.dump output.png first_cutoff second_cutoff

which will perform the SRO calculation using the dump file ``input.dump`` with cutoffs ``first_cutoff`` and ``second_cutoff``, saving the plot to ``output.png``.

If you have multiple dump files labelled like ``input.*.dump``, it should be fine to write ``input.*.dump`` in the command above. I haven't tested this, though.
