nets-nodegroups - Research node group structure of networks
===========================================================

To research and explore node group structure of various networks, we introduced the *group type parameter Tau* and implemented a *node group extraction framework*. The framework is capable of identifying groups, such as communities, modules, core/periphery, hubs&spokes, or similar structures. Detailed description of the algorithm can be found in:

- L. Šubelj, N. Blagus, and M. Bajec, "Group extraction for real-world networks: The case of communities, modules, and hubs and spokes," in Proc. of NetSci '13, 2013, p. 152.
- L. Šubelj, S. Žitnik, N. Blagus, and M. Bajec, “Node mixing and group structure of complex software networks,” Advs. Complex Syst., vol. 17, 2014.

The adopted network node group extraction framework extracts groups from a simple undirected graph sequentially. An optimization method (currently random-restart hill climbing) is used to maximize the group criterion *W(S,T)* and extract group *S* with the corresponding linking pattern *T*. After extraction edges between *S* and *T* are removed and the whole process repeated on the largest weakly-connected component until the group criterion *W* is larger than expected on a Erdös-Rényi random graph.

Cite this code repository by using DOI:

[![DOI:10.5281/zenodo.11589](https://zenodo.org/badge/doi/10.5281/zenodo.11589.png)](http://dx.doi.org/10.5281/zenodo.11589)


Compile
-------

This code should work on any platform using any standard C++ compiler (eg. GCC, Visual Studio), but it has only been extensively tested in one environment.

- Debian 7.3
- *build-essential* (~11.5)
- *g++* (~4:4.7.2-1)
- *SNAP* (~2.1) (newest from <https://github.com/snap-stanford/snap>)

First you need to get the *source code* (newest from <https://github.com/gw0/nets-nodegroups>):

    git clone http://github.com/gw0/nets-nodegroups.git
    cd ./nets-nodegroups


Compile *SNAP* from console using GCC:

    wget http://snap.stanford.edu/releases/Snap-2.1.zip
    unzip Snap-2.1.zip
    mv Snap-2.1 snap
    cd ./snap
    make all
    cd ..

Compile *net-nodegroups* from console using GCC:

    cd ./src
    make all

Executable file is put in `./src/nodegroups`.


Usage
-----

Parameters:

    -o: Input and output file name prefix (can be overriden) (default:'graph')
    -i: Input file with graph edges (undirected edge per line) (default:'graph.edgelist')
    -l: Optional input file with node labels (node ID, node label) (default:'graph.labels')
    -og: Output file with ST-group assignments (default:'graph.groups')
    -os: Output file with only ST-group extraction summary (default:'graph.groupssum')
    -n: Number of restarts of the optimization algorithm (default:2000)
    -sm: Maximal number of steps in each optimization run (default:100000)
    -sw: Stop optimization if no W improvement in steps (default:1000)
    -ss: Initial random-sample size of S ant T (0=random) (default:1)
    -fn: Finish after extracting so many groups (default:0)
    -fw: Finish if W smaller than top percentile on random graphs (default:1.0)
    -rg: Number of different Erdos-Renyi random graphs (default:500)
    -rn: Number of restarts on each Erdos-Renyi random graph (default:10)
    -rf: Force W recomputation on random graphs when relative difference smaller (default:inf)

Example command (read from `graph.edgelist` and output to `graph.groups` and `graph.groupsreport`):

    ./nodegroups -o:graph

Input file `graph.edgelist` with graph edges (`-i:`):

    0 1
    0 2
    1 2
    2 3
    3 4
    3 5
    ...

Output file `graph.groups` with extracted groups *S* and linking patterns *T* (`-og:`):

    # NId GroupS GroupT NLabel
    0     0      0      foo
    1     0      0      bar
    2     0      0      foobar
    3     1      -1     -
    2     -1     1      -
    ...

Output file `graph.groupsreport` with report on extracted groups (`-or:`):

    # Graphs: 12  Nodes: 115  Edges: 613
    N   M   N_S M_S N_T M_T N_ST    M_ST    L_ST    L_STc   W   Tau Mod_S   Mod_T   Type
    115 613 9   36  9   36  9   36  72  25  823.0000    1.0000  0.1352  0.1352  COM
    ...

Description of group type parameter Tau names (column `Type`):

- `COM`: community (*S = T*)
- `MOD`: module (*S intersection with T = 0*)
- `HSD`: hub&spokes module (module and *|T| = 1*)
- `MIX`: mixture (otherwise)
- `CPX`: core/periphery mixture (*S subset of T* or *T subset of S*)


Feedback
--------

If you encounter any bugs or have feature requests, please file them in the [issue tracker](https://github.com/gw0/nets-nodegroups/issues), or even develop it yourself and submit a pull request over [GitHub](https://github.com/gw0/nets-nodegroups).


License
-------

Copyright (c) 2014 *gw0* [<http://gw.tnode.com/>] &lt;<gw.2014@tnode.com>&gt;

This library is licensed under the [GNU Affero General Public License 3.0+](LICENSE_AGPL-3.0.txt) (AGPL-3.0+). Note that it is mandatory to make all modifications and complete source code of this library publicly available to any user.

