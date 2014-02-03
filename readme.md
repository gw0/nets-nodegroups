Node group structure
====================

Implementation of the *node group extraction framework* enables the exploration of node group structure of different networks, such as communities, modules, hubs and similar. Description of the algorithm can be found in:

- L. Šubelj, N. Blagus, and M. Bajec, "Group extraction for real-world networks: The case of communities, modules, and hubs and spokes," in Proc. of NetSci '13, 2013, p. 152.
- L. Šubelj, S. Žitnik, N. Blagus, and M. Bajec, "Node mixing and group structure of complex software networks," Adv. Complex Syst., 2014.


Compile
-------

The code should work on any platform using any standard C++ compiler (eg. GCC, Visual Studio), but it has only been extensively tested in the following environment.

Testing environment:

- Debian 7.3
- *build-essential* (~11.5)
- *g++* (~4:4.7.2-1)
- SNAP (~2.1) (download from <http://snap.stanford.edu/snap/download.html>)

Compile from console using GCC:

    cd src
    make all


Usage
-----

Parameters:

    -i: Input graph (list of undirected edges) (default:'graph.edgelist')
    -l: Optional input node labels (node ID, node label)
    -o: Output extracted group assignments (default:'graph.groups')

Example:

    ./ngs -i:bullet.edgelist
