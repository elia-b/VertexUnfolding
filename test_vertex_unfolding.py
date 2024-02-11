import pathlib

import vertex_unfolding
import os

SHAPES_REPO = pathlib.Path(__file__).parent / "shapes"
OUTPUT_REPO = pathlib.Path(__file__).parent / "outputs"


def test_vertex_unfolding():
    for shape in os.listdir(SHAPES_REPO):
        shape_name = shape[:-4]
        output_name = f"{shape_name}-vertex-unfolding.png"
        spanning_tree_name = f"{shape_name}-spanning-tree.png"
        facet_path_name = f"{shape_name}-facet-path.png"
        vertex_unfolding.unfold(
            stl_file=SHAPES_REPO / shape,
            output=OUTPUT_REPO / output_name,
            spanning_tree_output=OUTPUT_REPO / spanning_tree_name,
            facet_path_output=OUTPUT_REPO / facet_path_name,
        )
        assert output_name in os.listdir(OUTPUT_REPO)
        assert spanning_tree_name in os.listdir(OUTPUT_REPO)
        assert facet_path_name in os.listdir(OUTPUT_REPO)
