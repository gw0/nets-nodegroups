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

- Debian 7.3, 8.0
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

    -o:Prefix for all file names (simplified usage) (default:'graph')
    -i:Input file with graph edges (undirected edge per line) (default:'graph.edgelist')
    -l:Input file (optional) with node labels (node ID, node label) (default:'graph.labels')
    -og:Output file with ST-group assignments (default:'graph.groups')
    -os:Output file with only ST-group extraction summary (default:'graph.groupssum')
    -n:Number of restarts of the optimization algorithm (default:2000)
    -sm:Maximal number of steps in each optimization run (default:100000)
    -sw:Stop optimization if no W improvement in steps (default:1000)
    -ss:Initial random-sample size of S ant T (0=random) (default:1)
    -fn:Finish after extracting so many groups (turn off random graphs) (default:0)
    -fw:Finish if W smaller than top percentile on random graphs (default:1)
    -rg:Random graphs (Erdos-Renyi) to construct for estimating W (default:500)
    -rn:Random graph restarts of the optimization algorithm (default:10)
    -rf:Random graph re-estimation of W if relative difference smaller (default:inf)

Example to read from `graph.edgelist` and output to `graph.groups` and `graph.groupsreport`:

    ./nodegroups -o:graph

Same example but extract only first 12 groups (ignoring estimated W on random graphs):

    ./nodegroups -o:graph -fn:12

Input file `graph.edgelist` contains undirected graph edges (`-i:`):

    0 1
    0 2
    1 2
    2 3
    3 4
    3 5
    ...

Output file `graph.groups` contains extracted node groups *S* and linking patterns *T* (`-og:`):

    # NId GroupS GroupT NLabel
    0     0      0      foo
    1     0      0      bar
    2     0      0      foobar
    3     1      -1     -
    2     -1     1      -
    ...

Output file `graph.groupsreport` contains a summary of extracted node groups (`-or:`):

    # Graphs: 12  Nodes: 115  Edges: 613
    N   M   N_S M_S N_T M_T N_ST M_ST L_ST L_STc W        Tau    Mod_S  Mod_T  Type
    115 613 9   36  9   36  9    36   72   25    823.0000 1.0000 0.1352 0.1352 COM
    115 582 9   36  9   36  9    36   72   30    818.0000 1.0000 0.1164 0.1164 COM
    115 550 10  40  10  40  10   40   80   30    810.0000 1.0000 0.1191 0.1191 COM
    ...

Description of columns:

- `N`: number of nodes left in graph
- `M`: number of edges left in graph
- `N_S`: number of nodes in subgraph on group *S*
- `M_S`: number of edges in subgraph on group *S*
- `N_T`: number of nodes in subgraph on linking pattern *T*
- `M_T`: number of edges in subgraph on linking pattern *T*
- `N_ST`: number of nodes in subgraph on intersection of *S* and *T*
- `M_ST`: number of edges in subgraph on intersection of *S* and *T*
- `L_ST`: number of edges *L(S,T)* between groups *S* and *T*
- `L_STc`: number of edges *L(S,Tc)* between groups *S* and complement of *T*
- `W`: group critetion *W(S,T)*
- `Tau`: group type parameter *Tau(S,T)*
- `Mod_S`: modularity measure on group *S*
- `Mod_T`: modularity measure on linking pattern *T*
- `Type`: human name for group type parameter *Tau*:
  - `COM`: community (*S = T*)
  - `MOD`: module (*S intersection with T = 0*)
  - `HSD`: hub&spokes module (module and *|T| = 1*)
  - `MIX`: mixture (otherwise)
  - `CPX`: core/periphery mixture (*S subset of T* or *T subset of S*)

Exit status codes:

- `-1`: error, file not found or crash
- `0` : success, groups extracted
- `1` : no groups extracted


Feedback
--------

If you encounter any bugs or have feature requests, please file them in the [issue tracker](https://github.com/gw0/nets-nodegroups/issues), or even develop it yourself and submit a pull request over [GitHub](https://github.com/gw0/nets-nodegroups).


License
-------

Copyright (c) 2014 *gw0* [<http://gw.tnode.com/>] &lt;<gw.2014@tnode.com>&gt;

This library is licensed under the [GNU Affero General Public License 3.0+](LICENSE_AGPL-3.0.txt) (AGPL-3.0+). Note that it is mandatory to make all modifications and complete source code of this library publicly available to any user.

