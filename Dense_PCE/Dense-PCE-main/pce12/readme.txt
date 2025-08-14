######################################################################
 PCE: Pseudo Clique Enumerator, ver. 1.0
    10/July/2006 Takeaki Uno    e-mail:uno@nii.jp, 
    homepage:   http://research.nii.ac.jp/~uno/index.html
######################################################################

** This program is available for only academic use, basically.   **
** Anyone can modify this program, but he/she has to write down  **
** the change of the modification on the top of the source code. **
** Neither contact nor appointment to Takeaki Uno is needed.     **
** If one wants to re-distribute this code, do not forget to     **
** refer the newest code, and show the link to homepage          **
** of Takeaki Uno, to notify the news about codes for the users. **
** For the commercial use, please make a contact to Takeaki Uno. **

##################################
####    Problem Definition    ####
##################################

For a graph G=(V,E) of vertex set V and edge set E, a clique is a subset
of V such that any two its vertices are connected by an edge of E.
A clique is called maximal if it is included in no other clique.
Here our aim is to enumerate vertex sets which induce clique like subgraphs,
which we call pseudo cliques. Let e(n) denote the number of edges in a
clique of n vertices. For a vertex set S, we define the density of S by 
 |E(S)|/e(n),  where E(S) denotes the edge set of the subgraph induced by S.
Intuitively, the density is the ratio of the number of edges compared to
a clique. For a threshold value theta, 0<=theta<=1, we call S a pseudo
clique if its density is no less than theta.
A pseudo clique is called locally maximal if there is no pseudo clique
obtained by one vertex to the pseudo clique.
This program enumerates all pseudo cliques or locally maximal pseudo cliques 
of the given graph and threshold value theta.
If the threshold value is 1.0, then the problem is to enumerates all cliques.
For this task, mace has a better performance, about 10 times faster than PCE.


#####################
####    Usage    ####
#####################

The program is written in C code (gcc). It uses only the basic library,
so you can compile it in any environment. To compile the program, first
put all the source files in a directly, and just execute 
 
 % make
 
Then, you can execute "pce". The format of the input parameter is,

 pce MCcSs [options] input-filename threshold [output-filename]

PCE is executed with at least two parameters. The first parameter is
composed of commands, given by a combination of letters. The meaning
of the letters are:

 _: do not output the reports of the execution, such as size of input graph,
      to the standard output
 %: show progress, by outputting "-----" for each 100,000 cliques
 +: if the output file exists, append the solutions to the output file
 C: enumerate pseudo cliques of the give graph
 M: enumerate locally maximal cliques of the given graph
 s: terminate after finding 1,000,000 pseudo cliques to be output

The second parameter is the name of the input file, the third is the 
threshold value for the density, and the fourth parameter is the name of
output file. The threshold value is given by real number ranging from 0 to 1.
You can omit the output file name. In this case, the program
only counts the numbers of cliques to be output, classified by their sizes.
Note that the name of the input file can not start with '-'. If the
output file name is "-", the solutions will be output to the standard
output.

Between the first parameter and the second parameter, you can give some
options as follows.

 -# [num]:stop after outputting [num] solutions
 -, [char]: give the separator of the numbers in the output
    the numbers in the output file are separated by the given
    character [char].
 -Q [filename]: replace the output numbers 
    according to the permutation table written in the file of
    [filename], replace the numbers in the output. The numbers in the
    file can be separated by any non-numeric character such as newline
    character. 
 -l [num]  :output cliques with size at least [num] 
 -u [num]  :output cliques with size at most [num]

For example, by executing

   pce C -l 5 -u 7 g100.grh 0.9 clq.out

the program finds all pseudo cliques of sizes from 5 to 7 (including both
size 5 and size 7) of density at least 90% in the "g100.grh", 
and output to the file "clq.out"

Example of the execution)

   pce MS g15.grh 75 clq.out 
(output locally maximal pseudo cliques with at least 75% of densities in
"g15.grh" to "clq.out". Show the progress during the execution)



########################################
####    Input/Output File Format    ####
########################################

At first PCE outputs the information and statistics of the input graph
to standard error. The input graph, with n vertices, is considered to
be composed of vertices from 0 to n-1. The pseudo cliques found are
output to the output file specified by user, represented by a sequence
of numbers corresponding to the vertices of the pseudo clique. One
line of the output file is for one pseudo clique. The numbers in each
line is separated by " ", and by giving "-," option we can change the
separator. At the termination of the program, it outputs the number of
pseudo cliques found, and the number of pseudo cliques classified by
their sizes. For example, if there are pseudo cliques {0,1}, {2},
{0,1,3}, {1,2}, the output to the standard output will be

-----------
#vertices=??? #edges ????    <= numbers of vertices and edges
4     <= total number
0     <= number of cliques of size 0
1     <= number of cliques of size 1
2     <= number of cliques of size 2
1     <= number of cliques of size 3
-----------

and the output file will be

-----------
0 1
2
0 1 3
1 2
[EOF]
-----------

If q is given in the parameter, then "#vertices, #edges" is not printed.
If output file name is not given, then no output file is generated.

The output pseudo cliques are not sorted. If you want to sort it, use the
script "sortout.pl". The usage is just,

% sortout.pl < input-file > output-file

"input-file" is the name of file to which pce output, and the sorted output
will be written in the file of the name "output-file".
The vertices of each pseudo clique will be sorted in the increasing order of
indices, and all the pseudo cliques (lines) will be also sorted, by the
lexicographical order (considered as a string).
(Actually, you can specify separator like sortout.pl ",").

====   The format of the input file   =====

Each following i-th line is the list of vertices adjacent to vertex
i-1, so each line is: (vertex), (vertex), ... Any non-numeric letter
(except for newline and end-of-file) is allowed to be used for the
separator. Each vertex has to be ranged from 0 to (#vertices-1).
For an edge connecting vertices u and v, we do not need to write both
"v in the (u-1)th line" and "u in the (v-1)th line".
We need just one, "v in the (u-1)th line" or "u in the (v-1)th line".


Example) a graph with edges (1,0), (2,0), (1,3), (2,3), (3,4), (0,4):

1,2   
3

2 4
0
[EOF]


###############################################
###  Transforming Other Graph File Format   ###
###############################################

For the use of other formats for input graph files, we have several scripts.
We explain the functions of the Perl scripts in the following.

-- compgrh.pl [b] [separator] < input-file > output-file
Write to output-file the complement graph of the graph read from the
input-file. If you specify "b" option, then the input-file is regarded as
a bipartite graph.

Ex.)
% compgrh.pl b "," < test.grh > test2.grh

-- transnum.pl table-file [separator] < input-file > output-file
Read file from standard input, and give a unique number to each name written
by general strings (such as ABC, ttt), and transform every string name to 
a number, and output it to standard output. The mapping from string names to
numbers is output to table-file. The default character for the separator of 
the string names is " "(space). It can be changed by giving a character for
the option [separator]. For example, A,B is a string name, if the separator 
is space, but if we set the separator to ",", it is regarded as two names 
A and B. This is executed by "transnum.pl table-file "," < input-file...".

-- untransnum.pl table-file < input-file > output-file
According to the table-file output by transnum.pl, un-transform numbers to 
string names. The output of the program is composed of numbers, thus 
it is used if we want to transform to the original string names.
It reads file from standard output, and output to the standard output.

-- transgrh.pl [Bb] [separator] < input-file > output-file
Transform a file in the format of that every line writes two end vertices 
of an edge, to the format of this program. For example, the file

0,2
0,1
1,4
3,4
1,3
[EOF]

, representing the graph with edge (0,2), edge (0,1) and..., is transformed to

1,2
3,4

4
[EOF]

If the name of vertices are written in general strings, transform them to 
numbers by transnum.pl before the execution.
  When we give parameter D or d, the input graph is regarded as a directed
graph. If we give d, the first number is the origin of a directed edge, 
and if we give D, the second number is the origin.
  Bipartite graph can be transformed by giving options b or B. Suppose that
we have a bipartite graph with vertex sets A and B both indices start from 0.
If we transform it by the above way, vertex i in A and i in B are considered
to be the same vertex. By giving b, the second number in each row is added
by a constant so that no two vertices have the same index. In the case of B, 
the first number is added.
When we transform the above graph by ( transgrh.pl b "," < input > output ),
we have a graph in the format

5,6
7,8

8
[EOF]

-- sortgrhid
Transform a graph file format such that each row corresponds to a vertex.
The first number of each row is the ID of the vertex, and the following 
numbers are the vertex ID's to which the vertex adjacent.
The script sorts the row according to the ID's, and remove the ID from the file.
For example, the file,

-------------
4,2,3
1,2
0,1,3,4
3
2,3
[EOF]
-------------

is transformed to

-------------
1,3,4
2
3

2,3
[EOF]
-------------

If the ID's are given in general strings, use transnum.pl before the execution.

-- example of the usage: when transform the file test.grh which is edge list
 format with general string vertex names with separator ",", and enumerate
 chordless cycles:

transnum.pl table "," < test.grh > tmp.grh
transgrh.pl < tmp.grh > tmp2.grh
pce C tmp2.grh out
untransnum.pl table < out > out2


###########################
####    Performance    ####
###########################

This algorithm is polynomial delay, taking at most O(|V|) time for each pseudo
clique in the case of clique enumeration, where |V| is the number of
vertices in the input graph. In the case of locally maximal
pseudo cliques, it also enumerates all pseudo cliques and output only
locally maximal ones, so the computation time may not be output sensitive,
i.e., it may take much time for one locally maximal pseudo clique.
Practical computation time is constant for each pseudo clique clique, and
finds about 100,000 pseudo cliques per second, if the graph is sparse.
As increase the density of the graph, the computation time increases
almost linearly.


##################################################
####    Algorithm and Implementation Issue    ####
##################################################

Let G=(V,E) be a graph with vertex set V=(1,...,n) and edge set E defined
on V. To enumerate all pseudo cliques, we look at the adjacency of pseudo 
cliques, and use a useful traversing route defined on the pseudo cliques.
For a vertex subset U of V, we denote by G[U] the subgraph induced by U.
Let K be a pseudo clique which is not an empty vertex set, and u be a minimum
degree vertex in G[U]. Then, we can see that K-{u} is also a pseudo clique.
Let u*(K) be the minimum vertex among the minimum degree vertices in G[U].
We define the parent of K by K-u*(K). Since the parent has vertices one less 
than K, thus the binary relation parent-child induces a rooted tree whose 
root is an empty set (we regard an empty set as a pseudo clique).
By traversing this rooted tree, starting from the root, by depth first
search, we can enumerate all pseudo cliques. The depth first search needs
an algorithm for finding all children of a given pseudo clique, and it is
enough, since what we have to do is recursively find children of pseudo 
cliques.

  ====  updating adjacency to K ====
Now the important matter for efficient computation is how to find children
speedy. For a vertex not in K, let d(v,K) be the number of vertices in K
adjacent to K. K\cup {v} is a pseudo clique if and only if
(|E(K)|+d(v,K))/e(n+1) is no less than theta. This condition is written as
d(v,K) >= theta*e(n+1)-|E(K)|. Thus, we have to check only the vertices
satisfy this condition. To find the vertices satisfying this condition,
we update d(v,K) for all vertices when we add/remove a vertex to/from K,
and classify vertices into groups according to values of d(v,K). Then,
for each l>=theta*e(n+1)-|E(K)|, we can get the vertices satisfying
d(v,K) = l by looking at the group. When we add/remove a vertex u to from
K, d(v,K) changes if and only if v is adjacent to K. Thus, we
increment/decrement d(v,K) for vertices v adjacent to u. This can be done
in linear time of the degree of v.


################################################
####    Introduction: Graphs and Cliques    ####
################################################

A graph is an object, composed of vertices, and edges such that each
edge connects two vertices. The following is an example of a graph,
composed of vertices {0,1,2,3,4,5} and edges connecting 0 and 1, 1 
and 2, 0 and 3, 1 and 3, 1 and 4, 3 and 4.

0-1-2
|/|
3-4

A clique is a set of vertex in the graph. For example, the vertices 
0, 1, and 3 compose a clique. A vertex itself is a clique, and two
vertices connected by an edge is also a clique. A maximal clique is
a clique included in no other clique. For example, vertices 0 and 1 
compose a clique, but it is not maximal clique. vertices 0, 1, and 3
compose a clique, and it is also a maximal clique.

Clique is used widely for modeling a kind of related objects.
Enumerating all cliques is a popular method in database analysis
and solving the models. However, to be a clique is in some sense
a very strict condition. If there are some errors or incompleteness
in the database, some objects which should be considered as a group will
not be a clique. Moreover, in some problems, we just want to find dense
parts of the input graph. In such cases we should enumerate pseudo cliques.
A pseudo clique, often called dense subgraph, is a subgraph having many 
edges. With a threshold value theta, a vertex set S is a pseudo clique
if at least theta-ratio pair of vertices are connected.
For example, if the threshold value theta is 0.8, vertex subset {0,1,3,4}
of the above graph is a pseudo clique, but {1,2,4} is not.
In the following, we show some applications of pseudo cliques 
to real-world problems.

==== Finding Web Communities ===
Consider a web network graph, in which vertices are web sites,
and two vertices are connected by an edge if the corresponding 
sites are linked. Then, a clique of this web network graph can 
be considered to compose a community, since they have link in 
any pair of its members. Finding web communities by hand is not 
an easy task. By enumerating the cliques in the web network graph,
we can automatically find candidates and seeds of communities.

==== Finding Clusters ====
Suppose that we have data of collection of objects, and we want to 
find groups which we can consider that the members of the group 
are similar, or have common features. Such a group is called a cluster.
Consider a relation graph such that two objects similar to each other
is connected. Then, pseudo cliques are considered as seeds of clusters.
By joining cliques.


###############################
####    Acknowledgments    ####
###############################

We would like to thank Koji Tsuda of AIST, Japan and Elizabeth ??? for
their comments and advices. We also thank Stephanie Henrichs of
University of South Carolina for a bug report.



