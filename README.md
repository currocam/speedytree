## speedytree


Speedytree is a Rust implementation of Neighbor-Joining for building phylogenetic trees from large Phylip distance matrices. There are two strategies: the Canonical algorithm (as [QuickTree](https://github.com/khowe/quicktree)) and something in the spirit of [RapidNJ](https://birc.au.dk/software/rapidnj/) but with B-trees.


You can read more about Neighbor-Joining [here](https://en.wikipedia.org/wiki/Neighbor_joining). The RapidNJ algorithm should be faster for very big problems at the cost of a larger memory overhead. 

You may also find interesting the [slides](speedytree_slides.pdf) I made about this project.

### Benchmark


You can run the benchmark with `make benchmark`. As an example, a small sth file is included, but you probably want to add your own. The alignments I used are included in the release with the binary executable file.


### Testing


Speedytree has a test suite that can be run with `cargo test`. There are three types of tests:


- Unit tests, which test individual functions.
- Integration tests, which test the whole program with small matrices.
- Property-based tests, which test the program with randomly additive binary trees. These tests rely on the fact that the program should output the same tree that it was given as input.




### About the command line application

#### Installation
You can either download the binary from the releases page or build it yourself. To build it yourself you need to have Rust installed. Then you can run `cargo build --release` to build the binary. The binary will be located in `target/release/speedytree`.


#### Usage
Speedytree takes a Phylip distance matrix as input and outputs a Newick tree. It can be used like this:


```
speedytree < input.phy > output.nwk
```


Speedytree has a few options that can be used to tweak the output. You can see them by running `speedytree --help`. The most important options are:


- `-c` to set the number of threads to use. By default it will use 1.
- `--naive` to use the canonical implementation. This algorithm is equivalent to QuickTree, and it's fast in practice for small matrices.
- `--rapidnj` to use the RapidNJ heuristics, but implemented with BTrees.
- `--hybrid` to use a mix of the two algorithms.



