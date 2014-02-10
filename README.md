nodegroups - Node group structures
==================================

Implementation of the *node group extraction framework* enables the exploration of node group structures of different networks, such as communities, modules, core/periphery, hubs&spokes, or similar structures. Description of the algorithm can be found in:

- L. Šubelj, N. Blagus, and M. Bajec, "Group extraction for real-world networks: The case of communities, modules, and hubs and spokes," in Proc. of NetSci '13, 2013, p. 152.
- L. Šubelj, S. Žitnik, N. Blagus, and M. Bajec, "Node mixing and group structure of complex software networks," Adv. Complex Syst., 2014. (in review)


Compile
-------

The code should work on any platform using any standard C++ compiler (eg. GCC, Visual Studio), but it has only been extensively tested in the following environment.

Testing environment:

- Debian 7.3
- *build-essential* (~11.5)
- *g++* (~4:4.7.2-1)
- *SNAP* (~2.1) (newest from <https://github.com/snap-stanford/snap>)

Compile *SNAP* from console using GCC:

    wget http://snap.stanford.edu/releases/Snap-2.1.zip
    unzip Snap-2.1.zip
    mv Snap-2.1 snap
    cd ./snap
    make all
    cd ..

Compile *nodegroups* from console using GCC:

    cd ./src
    make all


Usage
-----

Parameters:

    -i: Input graph (list of undirected edges) (default:'graph.edgelist')
    -l: Optional input node labels (node ID, node label)  (default:'graph.labels')
    -o: Output extracted group assignments (default:'graph.groups')

Example:

    ./nodegroups -i:graph.edgelist


Feedback
--------

If you encounter any bugs or have feature requests, please file them in the [issue tracker](https://github.com/gw0/nodegroups/issues), or even develop it yourself and submit a pull request over [GitHub](https://github.com/gw0/nodegroups).


License
-------

Copyright (c) 2014 - *gw0* [<http://gw.tnode.com/>] &lt;<gw.2014@tnode.com>&gt;

This library is licensed under the [GNU Affero General Public License 3.0+](LICENSE_AGPL-3.0.txt) (AGPL-3.0+). Note that it is mandatory to make all modifications and complete source code of this library publicly available to any user.
