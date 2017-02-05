List of scripts
===============

getmyancestrydna.py
-------------------

getmyancestrydna.py is a python3 script that downloads DNA matches sharing information from AncestryDNA

This script requires the python3 requests module to work. To install this module, run this in your terminal: "python3 -m pip install requests" (or "python3 -m pip install --user requests" if you don't have admin rights on your machine)

getmy23andme.py
---------------

getmy23andme.py is a python3 script that downloads DNA match sharing information from 23andMe

This script requires the python3 requests module to work. To install this module, run this in your terminal: "python3 -m pip install requests" (or "python3 -m pip install --user requests" if you don't have admin rights on your machine)

ancestry2graph.py
-----------------

ancestry2graph.py is a python3 script that converts the output of getmyancestrydna.py into a graph file

ibdview2graph.py
----------------

ibdview2graph.py is a python3 script that converts the output of getmy23andme.py into a graph file

To obtain distances in centiMorgans it requires a genetic map for the GRCh37 genome

graph2matrix.py
---------------

graph2matrix.py is a python3 script that converts a graph file into matrix format which can be subsequently loaded into Gephi

graph2plot.py
-------------

graph2plot.py is a python script that generates sharing plots from graph files

This script requires the networkx module to work. To install this module, run this in your terminal: "python -m pip install networkx" (or "python -m pip install --user networkx" if you don't have admin rights on your machine)

This script requires the pydot module to work. To install this module, run this in your terminal: "python -m pip install pydot" (or "python -m pip install --user pydot" if you don't have admin rights on your machine)

Unfortunately this script requires python 2.7 to run due to some portability issues within the pydot code

matches2plot.py
-------------

matches2plot is a python script that shows relative sharing of DNA matches with two separate individuals in your account

Examples
========

download and unpack the genetic map
-----------------------------------

wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip

unzip -d . -o plink.GRCh37.map.zip

download DNA match sharing information from 23andMe
---------------------------------------------------

./getmy23andme.py -u %USERNAME% -p %PASSWORD% -x -v

convert your 23andMe information into a graph file
--------------------------------------------------

./ibd2graph.py -h %ACCOUNT_ID%.inheritance.tsv -i %ACCOUNT_ID%.ibdview.tsv -o %ACCOUNT_ID%.graph.tsv

convert your 23andMe information into a graph file and add centimorgan information
----------------------------------------------------------------------------------

./ibd2graph.py -h %ACCOUNT_ID%.inheritance.tsv -i %ACCOUNT_ID%.ibdview.tsv -c {1..22} X -g plink.chr{{1..22},X}.GRCh37.map -o %ACCOUNT_ID%.graph.tsv

Notice that on Windows the above line will not work as the default Windows shell will not expand the curly braces. Therefore instead of {1..22} X you will have to write:

1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X

And instead of plink.chr{{1..22},X}.GRCh37.map you will have to write:

plink.chr1.GRCh37.map plink.chr2.GRCh37.map plink.chr3.GRCh37.map plink.chr4.GRCh37.map plink.chr5.GRCh37.map plink.chr6.GRCh37.map plink.chr7.GRCh37.map plink.chr8.GRCh37.map plink.chr9.GRCh37.map plink.chr10.GRCh37.map plink.chr11.GRCh37.map plink.chr12.GRCh37.map plink.chr13.GRCh37.map plink.chr14.GRCh37.map plink.chr15.GRCh37.map plink.chr16.GRCh37.map plink.chr17.GRCh37.map plink.chr18.GRCh37.map plink.chr19.GRCh37.map plink.chr20.GRCh37.map plink.chr21.GRCh37.map plink.chr22.GRCh37.map plink.chrX.GRCh37.map

plot your 23andMe graph file
----------------------------

./graph2plot.py -i %ACCOUNT_ID%.graph.tsv -r %EHID% -rel %ACCOUNT_ID%.%EHID%.relfinder.tsv -o %EHID%.pdf

download DNA match sharing information from AncestryDNA
-------------------------------------------------------

./getmyancestrydna.py -u %USERNAME% -p %PASSWORD% -x -v

convert your AncestryDNA information into a graph file
------------------------------------------------------

./ancestry2graph.py -i %UCDMID%.%GUID%.tsv -o %GUID%.graph.tsv

plot your AncestryDNA graph file
--------------------------------

./graph2plot.py -r %GUID% -i %GUID%.graph.tsv -anc %UCDMID%.%GUID%.tsv -o %GUID%.pdf

convert your 23andMe graph file into a matrix that you can open with Gephi
--------------------------------------------------------------------------

./graph2matrix.py -t \; -l -v -c -g -i %ACCOUNT_ID%.graph.tsv -h %ACCOUNT_ID%.inheritance.tsv -o %ACCOUNT_ID%.matrix.csv

plot the amount of sharing from DNAmatches with a parent and a child
--------------------------------------------------------------------

./matches2plot -a %UCDMID%.%GUID1%.tsv -b %UCDMID%.%GUID2%.tsv

Support
=======

To learn more about the data being retrieved by this tool, see <a href="http://apol1.blogspot.com/2015/09/visualize-dna-matches-graph.html">here</a>

This set of programs is still in beta phase, and bugs are still present. Features will be added on request. It is provided as is

These scripts require both python 3.4 to run due to some novel features in the argparse module (https://docs.python.org/3/whatsnew/3.4.html#argparse) and python 2.7 due to current incompabilities within the python3 pydot package

Current version was updated on June 26th 2016

A special thank goes to Vincenzo Palleschi for the numerous discussions and acting as the first beta tester. Part of the inspiration of this project came from his <a href="https://peerj.com/articles/cs-27/">work</a> on pairwise genome sharing and the Shared Matches <a href="http://blogs.ancestry.com/ancestry/2015/08/26/see-your-dna-matches-in-a-whole-new-way/">tool</a> from AncestryDNA.

Send questions, suggestions, or feature requests to giulio.genovese@gmail.com
