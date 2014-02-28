nodegroups - Node group structures
==================================

Implementation of the *node group extraction framework* (into S, T) enables the exploration of node group structures of different networks, such as communities, modules, core/periphery, hubs&spokes, or similar structures. Description of the algorithm can be found in:

- L. Šubelj, N. Blagus, and M. Bajec, "Group extraction for real-world networks: The case of communities, modules, and hubs and spokes," in Proc. of NetSci '13, 2013, p. 152.
- L. Šubelj, S. Žitnik, N. Blagus, and M. Bajec, "Node mixing and group structure of complex software networks," Adv. Complex Syst., 2014. (in review)

The adopted group extraction framework extracts groups from a simple undirected graph sequentially. An optimization method (currently random-restart hill climbing) is used to maximize the group criterion *W(S,T)* and extract group *S* with the corresponding linking pattern *T*. After extraction edges between *S* and *T* are removed and the whole process repeated on the largest weakly-connected component until the group criterion *W* is larger than expected on a Erdös-Rényi random graph.


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

    -o: Input and output file name prefix (can be overriden) (default:'graph')
    -i: Input file with graph edges (undirected edge per line) (default:'graph.edgelist')
    -l: Optional input file with node labels (node ID, node label) (default:'graph.labels')
    -og: Output file with group assignments (for S and T) (default:'graph.groups')
    -or: Output file with only group reports (for S and T) (default:'graph.groupsreport')
    -n: Number of restarts of the optimization algorithm (default:2000)
    -sm: Maximal number of steps in each optimization run (default:100000)
    -sw: Stop optimization if no W improvement in steps (default:1000)
    -ss: Initial random-sample size of S ant T (0=random) (default:1)
    -fn: Finish after extracting so many groups (default:0)
    -fw: Finish if W smaller than top percentile on random graphs (default:1.0)
    -rg: Number of different Erdos-Renyi random graphs (default:500)
    -rn: Number of restarts on each Erdos-Renyi random graph (default:10)
    -rf: Force W recomputation on random graphs when relative difference smaller (default:inf)

Example command:

    ./nodegroups -o:graph

Example input `graph.edgelist` (`-i:`):

    0 1
    0 2
    1 2
    2 3
    3 4
    3 5
    ...

Example output `graph.groups` (`-og:`):

    # NId GroupS GroupT NLabel
    0     0      0      foo
    1     0      0      bar
    2     0      0      foobar
    3     1      -1     -
    2     -1     1      -
    ...

Example output `graph.groupsreport` (`-or:`):

    # Graphs: 12  Nodes: 115  Edges: 613
    N   M   N_S M_S N_T M_T N_ST    M_ST    L_ST    L_STc   W   Tau Mod_S   Mod_T   Type
    115 613 9   36  9   36  9   36  72  25  823.0000    1.0000  0.1352  0.1352  COM
    ...

Description of group type names:

- "COM": community (S = T)
- "MOD": module (S inters T = 0)
- "HSD": hub&spokes module (module and |T| = 1)
- "MIX": mixture (else)
- "CPX": core/periphery mixture (S subset T or T subset S)


Feedback
--------

If you encounter any bugs or have feature requests, please file them in the [issue tracker](https://github.com/gw0/nodegroups/issues), or even develop it yourself and submit a pull request over [GitHub](https://github.com/gw0/nodegroups).


License
-------

Copyright (c) 2014 - *gw0* [<http://gw.tnode.com/>] &lt;<gw.2014@tnode.com>&gt;

This library is licensed under the [GNU Affero General Public License 3.0+](LICENSE_AGPL-3.0.txt) (AGPL-3.0+). Note that it is mandatory to make all modifications and complete source code of this library publicly available to any user.
