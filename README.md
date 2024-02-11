Vertex Unfolding

From this:

![standford bunny](shapes/stanford_bunny.stl "Mesh")

To this:

![standford bunny](outputs/stanford_bunny-spanning-tree.png "Mesh")

To this:

![standford bunny](outputs/stanford_bunny-facet-path.png "Mesh")

To this:

![standford bunny](outputs/stanford_bunny-vertex-unfolding.png "Mesh")

RUNNING:

```commandline
python3 -m venv venv
. venv/bin/activate
pip install -r requirements.txt
```

Create your script and import `unfold` from the `vertex_unfolding` file:

```python
import vertex_unfolding

vertex_unfolding.unfold(stl_file="my_stl.stl", output="my_vertex_unfolding.png")
```

CONTRIBUTING:

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
