# Vertex Unfolding

This project implements vertex unfolding for an arbitrary polyhedron.
It is using the algorithm from Eric Demaine, David Eppstein, Jeff Erickson, George Hart and Joseph O'Rourke introduced in *Vertex-Unfoldings of Simplicial Manifolds*.

## Input
Input is a STL file containing a single polyhedron.

The stanford bunny as polyhedra:
![stanford bunny polyhedra](images/stanford-bunny.png)

## Output
Output is a PNG showing the unfolded representation of the given polyhedra.

Which looks like this for a simple cube:
![vertex unfolding for a cube](outputs/cube-vertex-unfolding.png)

Or for the stanford bunny showed before:
![vertex unfolding for a stanford bunny](outputs/stanford_bunny-vertex-unfolding.png)


## Computation of the Vertex Unfolding
### Spanning Tree
The first step to create the vertex unfolding is creating the spanning tree.
With the help of this spanning tree, the levels for each face of the polyhedra is calculated.
The level of a face is a number which describes how many faces are until a leave node of the spanning tree comes.

The created spanning tree for a cube with colored nodes depending on its level can look like this:
![spanning tree for a cube](outputs/cube-spanning-tree.png)

The more complex example of the stanford bunny looks like this:
![spanning tree for a stanford bunny](outputs/stanford_bunny-spanning-tree.png)

### Facet Path / Eulerian Path
With the help of the spanning tree and the levels of the facet, a facet path can be created.
This path is used to create an eulerian path.

To create a valid facet graph which has an eulerian path is a lot of magic due to many edge cases.

The created facet path for a cube looks like this (red nods are faces and blue one are vertices):
![facet path for a cube](outputs/cube-facet-path.png "Mesh")

And for the stanford bunny
![facet path for a stanford bunny](outputs/stanford_bunny-facet-path.png "Mesh")

### Vertex Unfolding
The eulerian path of the facet path is used to traverse the polyhedra for the vertex unfolding.
The triangles of the polyhedra are laid out in the sequence of the path and are connected at the specified vertices.
In order to make sure that the unfolding is valid (contains no overlaps) the triangles are rotated so that each has there one vertical space.

![vertex-unfolding-space.png](images/vertex-unfolding-layed-out.png)

## How to Use the Script:

### Install dependencies
Dependencies can be installed by running the following commands:
```commandline
python3 -m venv venv
pip install -r requirements.txt
```

### Run script
Create a script in which you want to use the vertex unfolding and import `unfold` from the `vertex_unfolding` file:

```python
import vertex_unfolding

if __name__ == "__main__":
    name = "my_stl"
    vertex_unfolding.unfold(
        stl_file=f"{name}.stl",
        output=f"{name}_vertex_unfolding.png",
        spanning_tree_output=f"{name}-spanning-tree.png",
        facet_path_output=f"{name}-facet-path.png"
    )
```

Parameter:
* `stl_file`: Path to the stl file containing the polyhedra
* `output`: Name of the output file

### CONTRIBUTING:

Install extra dependencies:

```commandline
pip install -r requirements-dev.txt
```

Run linting:

```commandline
pre-commit run --all-files
```

Run test:

```commandline
pytest
```
