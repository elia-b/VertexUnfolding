import pathlib
import igl
import networkx as nx
import numpy as np
import drawsvg as draw
import matplotlib.pyplot as plt
import scipy.spatial.transform
import math


def unfold(
    stl_file: str | pathlib.Path,
    output: str | pathlib.Path,
    spanning_tree_output: None | str | pathlib.Path = None,
    facet_path_output: None | str | pathlib.Path = None,
):
    if isinstance(stl_file, pathlib.Path):
        stl_file = stl_file.as_posix()
    raw_vertices, raw_faces = igl.read_triangle_mesh(
        filename=stl_file,
        dtypef=np.float64,
    )
    vertices, _, _, faces = igl.remove_duplicate_vertices(
        raw_vertices,
        raw_faces,
        0.00001,
    )
    graph = _graph_from_mesh(faces)
    spanning_tree = nx.minimum_spanning_tree(graph)
    level_array = _level_array(spanning_tree)

    if spanning_tree_output:
        plt.figure(1, figsize=(24, 24))  # compute size based on number of nodes
        max_level = max(level_array)
        node_color = [level_array[node] / max_level for node in spanning_tree.nodes]
        nx.draw(spanning_tree, with_labels=True, node_color=node_color)
        plt.savefig(spanning_tree_output)
        plt.clf()

    facet_graph = generate_facet_graph(spanning_tree, faces, level_array=level_array)

    _remove_single_vertices_connected_components(facet_graph)

    if facet_path_output:
        plt.figure(2, figsize=(14, 14))  # compute size based on number of nodes
        node_color = [
            "#115f9a" if str(node).startswith("v") else "#991f17"
            for node in facet_graph.nodes
        ]
        nx.draw(facet_graph, with_labels=True, node_color=node_color)
        plt.savefig(facet_path_output)
        plt.clf()

    eulerian_path = list(nx.eulerian_path(facet_graph))

    polygons, offsets = _get_polygons_and_offsets(
        eulerian_path=eulerian_path,
        faces=faces,
        vertices=vertices,
    )

    zoom = 5

    drawing = draw.Drawing(
        (offsets[-1][0] + 100) * zoom,
        1000,
        origin=(-50 * zoom, -500),
    )
    for polygon in polygons:
        drawing.append(
            draw.Lines(
                *np.array(polygon).flatten() * zoom,
                close=True,
                fill="#eeee00",
                stroke="#000",
                stroke_width=0.1,
            ),
        )

    # for offset in offsets:
    # drawing.append(draw.Circle(*np.array(offset) * zoom, 1 * zoom, fill='lime'))
    # drawing.append(draw.Text(x=offset[0]* zoom, y=offset[1]* zoom, font_size=7 * zoom, text=f"{round(offset[0], 1)}, {round(offset[1], 1)}", fill='red'))

    drawing.rasterize()
    if isinstance(output, pathlib.Path):
        output = output.as_posix()
    drawing.save_png(output)


"""
VERTEX UNFOLDING UTILS
This set of functions streamlines the vertex unfolding process
"""


def _level_array(spanning_tree: nx.Graph) -> np.ndarray:
    """
    Assigns levels to nodes in a spanning tree based on their distance from leaf nodes.

    Args:
    spanning_tree: A NetworkX graph representing the spanning tree.

    Returns:
    A NumPy array where each element represents the level of the corresponding node in the tree.
    Leaf nodes have level 1, their parents have level 2, and so on.
    """
    current_tree = spanning_tree.copy()
    next_tree = spanning_tree.copy()

    level_array = np.zeros(len(spanning_tree.nodes))

    iteration_no = 1
    while next_tree.number_of_nodes() >= 1:
        for node_id in list(current_tree.nodes):
            if current_tree.degree[node_id] <= 1:
                next_tree.remove_node(node_id)
                level_array[node_id] = iteration_no
        current_tree = next_tree.copy()
        iteration_no = iteration_no + 1

    return level_array


def generate_facet_graph(
    spanning_tree: nx.Graph,
    faces: np.ndarray,
    level_array: np.ndarray,
) -> nx.Graph:
    facet_graph = nx.Graph()

    cycles = []

    for even_iterator in range(2, int(max(level_array)) + 1, 2):
        for face_id in range(len(level_array)):
            if level_array[face_id] == even_iterator:
                neighbour_faces = list(spanning_tree.adj[face_id])
                n_adj_faces_lower_level = []
                for adj_face in neighbour_faces:
                    if (
                        level_array[adj_face] < even_iterator
                        and level_array[adj_face] % 2 == 1
                    ):
                        n_adj_faces_lower_level.append(adj_face)
                if len(n_adj_faces_lower_level) == 1:
                    adj_face = n_adj_faces_lower_level[0]
                    common_vertices = np.intersect1d(faces[face_id], faces[adj_face])
                    facet_graph.add_edge(face_id, "v" + str(common_vertices[0]))
                    facet_graph.add_edge(face_id, "v" + str(common_vertices[1]))
                    facet_graph.add_edge(adj_face, "v" + str(common_vertices[0]))
                    facet_graph.add_edge(adj_face, "v" + str(common_vertices[1]))
                    cycles.append([face_id, adj_face])
                if len(n_adj_faces_lower_level) == 2:
                    adj_face0 = n_adj_faces_lower_level[0]
                    adj_face1 = n_adj_faces_lower_level[1]
                    common_vertices0 = np.intersect1d(faces[face_id], faces[adj_face0])
                    common_vertices1 = np.intersect1d(faces[face_id], faces[adj_face1])
                    facet_graph.add_edge(adj_face0, "v" + str(common_vertices0[0]))
                    facet_graph.add_edge(adj_face0, "v" + str(common_vertices0[1]))
                    facet_graph.add_edge(adj_face1, "v" + str(common_vertices1[0]))
                    facet_graph.add_edge(adj_face1, "v" + str(common_vertices1[1]))
                    for v in faces[face_id]:
                        if v not in common_vertices1:
                            facet_graph.add_edge(face_id, "v" + str(v))
                    for v in faces[face_id]:
                        if v not in common_vertices0:
                            facet_graph.add_edge(face_id, "v" + str(v))
                    cycles.append([face_id, adj_face0, adj_face1])
                if len(n_adj_faces_lower_level) == 3:
                    adj_face0 = n_adj_faces_lower_level[0]
                    adj_face1 = n_adj_faces_lower_level[1]
                    adj_face2 = n_adj_faces_lower_level[2]
                    common_vertices0 = np.intersect1d(faces[face_id], faces[adj_face0])
                    common_vertices1 = np.intersect1d(faces[face_id], faces[adj_face1])
                    common_vertices01 = np.intersect1d(
                        faces[adj_face0],
                        faces[adj_face1],
                    )
                    common_vertices12 = np.intersect1d(
                        faces[adj_face1],
                        faces[adj_face2],
                    )
                    common_vertices02 = np.intersect1d(
                        faces[adj_face0],
                        faces[adj_face2],
                    )
                    facet_graph.add_edge(face_id, "v" + str(common_vertices0[0]))
                    facet_graph.add_edge(face_id, "v" + str(common_vertices0[1]))
                    if common_vertices0[0] in common_vertices01:
                        facet_graph.add_edge(adj_face0, "v" + str(common_vertices0[0]))
                    else:
                        facet_graph.add_edge(adj_face0, "v" + str(common_vertices0[1]))
                    facet_graph.add_edge(adj_face1, "v" + str(common_vertices1[0]))
                    facet_graph.add_edge(adj_face1, "v" + str(common_vertices1[1]))
                    assert len(common_vertices12) == 1
                    facet_graph.add_edge(adj_face2, "v" + str(common_vertices12[0]))
                    # Add extra edges, not in slide
                    for v in faces[adj_face0]:
                        if v not in common_vertices0:
                            facet_graph.add_edge(adj_face0, "v" + str(v))

                    facet_graph.add_edge(adj_face2, "v" + str(common_vertices02[0]))
                    cycles.append([face_id, adj_face0, adj_face1, adj_face2])

    # second phase
    for node in spanning_tree:
        if node not in facet_graph.nodes:
            neighbour_faces = list(spanning_tree.adj[node])
            assert level_array[node] % 2 == 1
            assert len(neighbour_faces) > 0
            if not level_array[node] == max(level_array):
                n_adj_faces_lower_level = []
                for adj_face in neighbour_faces:
                    if level_array[adj_face] % 2 == 1 and level_array[adj_face] != max(
                        level_array,
                    ):
                        n_adj_faces_lower_level.append(adj_face)
                if len(n_adj_faces_lower_level) < 1:
                    continue

                if len(n_adj_faces_lower_level) > 1:
                    pruned_n_adj_faces_lower_level = []
                    for f in n_adj_faces_lower_level:
                        if f not in facet_graph.nodes:
                            pruned_n_adj_faces_lower_level.append(f)

                    n_adj_faces_lower_level = pruned_n_adj_faces_lower_level
                assert len(n_adj_faces_lower_level) == 1

                neighbour = n_adj_faces_lower_level[0]
                cycle_length = None
                cycle_found = None
                cycle_index = None
                for index, cycle in enumerate(cycles):
                    if neighbour in cycle:
                        cycle_found = cycle
                        cycle_index = index
                        cycle_length = len(cycle)

                if cycle_length is None:
                    common_vertices = np.intersect1d(faces[node], faces[neighbour])
                    facet_graph.add_edge(neighbour, "v" + str(common_vertices[1]))
                    facet_graph.add_edge(neighbour, "v" + str(common_vertices[0]))
                    facet_graph.add_edge(node, "v" + str(common_vertices[1]))
                    facet_graph.add_edge(node, "v" + str(common_vertices[0]))
                else:
                    assert cycle_length is not None
                    if cycle_length == 2:
                        common_vertices = np.intersect1d(faces[node], faces[neighbour])
                        if facet_graph.has_edge(
                            neighbour,
                            "v" + str(common_vertices[1]),
                        ):
                            facet_graph.remove_edge(
                                neighbour,
                                "v" + str(common_vertices[1]),
                            )
                            facet_graph.add_edge(
                                neighbour,
                                "v" + str(common_vertices[0]),
                            )
                        elif facet_graph.has_edge(
                            neighbour,
                            "v" + str(common_vertices[0]),
                        ):
                            facet_graph.remove_edge(
                                neighbour,
                                "v" + str(common_vertices[0]),
                            )
                            facet_graph.add_edge(
                                neighbour,
                                "v" + str(common_vertices[1]),
                            )
                        else:
                            facet_graph.add_edge(
                                neighbour,
                                "v" + str(common_vertices[1]),
                            )
                            facet_graph.add_edge(
                                neighbour,
                                "v" + str(common_vertices[0]),
                            )
                        facet_graph.add_edge(node, "v" + str(common_vertices[1]))
                        facet_graph.add_edge(node, "v" + str(common_vertices[0]))
                        assert cycle_index is not None
                        cycles[cycle_index].append(node)
                    elif cycle_length == 3:
                        assert cycle_found is not None
                        cycle_leader = cycle_found[0]
                        common_vertices = np.intersect1d(faces[node], faces[neighbour])
                        cycle_neighbour = (
                            cycle_found[2]
                            if cycle_found[1] == neighbour
                            else cycle_found[1]
                        )
                        common_vertices_leader_neighbour = np.intersect1d(
                            faces[cycle_leader],
                            faces[neighbour],
                        )
                        common_vertices_cycle_neighbour_neighbour = np.intersect1d(
                            faces[cycle_neighbour],
                            faces[neighbour],
                        )
                        assert len(common_vertices_cycle_neighbour_neighbour) == 1
                        facet_graph.remove_edge(
                            neighbour,
                            "v" + str(common_vertices_leader_neighbour[0]),
                        )
                        facet_graph.remove_edge(
                            neighbour,
                            "v" + str(common_vertices_leader_neighbour[1]),
                        )
                        for v in common_vertices_leader_neighbour:
                            if v not in common_vertices_cycle_neighbour_neighbour:
                                facet_graph.remove_edge(cycle_leader, "v" + str(v))
                        facet_graph.add_edge(
                            cycle_leader,
                            "v" + str(common_vertices_cycle_neighbour_neighbour[0]),
                        )
                        facet_graph.add_edge(neighbour, "v" + str(common_vertices[1]))
                        facet_graph.add_edge(neighbour, "v" + str(common_vertices[0]))
                        facet_graph.add_edge(node, "v" + str(common_vertices[1]))
                        facet_graph.add_edge(node, "v" + str(common_vertices[0]))
                        cycles[cycle_index].remove(neighbour)
                        cycles.append([node, neighbour])
                    else:
                        raise ValueError("Cycle edgecase")

    # Deal with max values
    for node in spanning_tree:
        if node not in facet_graph.nodes:
            neighbour_faces = list(spanning_tree.adj[node])
            assert level_array[node] % 2 == 1
            assert len(neighbour_faces) > 0
            if level_array[node] == max(level_array):
                assert all(
                    [level_array[node] >= level_array[n] for n in neighbour_faces],
                )
                n_adj_faces_lower_level = []
                for adj_face in neighbour_faces:
                    if (
                        level_array[adj_face] % 2 == 1
                        and adj_face not in facet_graph.nodes
                    ):
                        n_adj_faces_lower_level.append(adj_face)
                extra_nodes = []
                for adj_face in n_adj_faces_lower_level:
                    if level_array[adj_face] == max(level_array):
                        adj_neighbour_faces = list(spanning_tree.adj[adj_face])
                        for adj_neighbour_face in adj_neighbour_faces:
                            if (
                                level_array[adj_neighbour_face] % 2 == 1
                                and adj_neighbour_face not in facet_graph.nodes
                            ):
                                extra_nodes.append(adj_neighbour_face)
                n_adj_faces_lower_level += extra_nodes
                face_id = node
                if len(n_adj_faces_lower_level) == 0:
                    facet_graph.add_edge(face_id, "v" + str(faces[face_id][0]))
                    facet_graph.add_edge(face_id, "v" + str(faces[face_id][1]))
                if len(n_adj_faces_lower_level) == 1:
                    adj_face = n_adj_faces_lower_level[0]
                    common_vertices = np.intersect1d(faces[face_id], faces[adj_face])
                    facet_graph.add_edge(face_id, "v" + str(common_vertices[0]))
                    facet_graph.add_edge(face_id, "v" + str(common_vertices[1]))
                    facet_graph.add_edge(adj_face, "v" + str(common_vertices[0]))
                    facet_graph.add_edge(adj_face, "v" + str(common_vertices[1]))
                if len(n_adj_faces_lower_level) == 2:
                    adj_face0 = n_adj_faces_lower_level[0]
                    adj_face1 = n_adj_faces_lower_level[1]
                    common_vertices0 = np.intersect1d(faces[face_id], faces[adj_face0])
                    common_vertices1 = np.intersect1d(faces[face_id], faces[adj_face1])
                    facet_graph.add_edge(adj_face0, "v" + str(common_vertices0[0]))
                    facet_graph.add_edge(adj_face0, "v" + str(common_vertices0[1]))
                    facet_graph.add_edge(adj_face1, "v" + str(common_vertices1[0]))
                    facet_graph.add_edge(adj_face1, "v" + str(common_vertices1[1]))
                    for v in faces[face_id]:
                        if v not in common_vertices1:
                            facet_graph.add_edge(face_id, "v" + str(v))
                    for v in faces[face_id]:
                        if v not in common_vertices0:
                            facet_graph.add_edge(face_id, "v" + str(v))
                if len(n_adj_faces_lower_level) == 3:
                    adj_face0 = n_adj_faces_lower_level[0]
                    adj_face1 = n_adj_faces_lower_level[1]
                    adj_face2 = n_adj_faces_lower_level[2]
                    common_vertices0 = np.intersect1d(faces[face_id], faces[adj_face0])
                    common_vertices1 = np.intersect1d(faces[face_id], faces[adj_face1])
                    common_vertices01 = np.intersect1d(
                        faces[adj_face0],
                        faces[adj_face1],
                    )
                    common_vertices12 = np.intersect1d(
                        faces[adj_face1],
                        faces[adj_face2],
                    )
                    common_vertices02 = np.intersect1d(
                        faces[adj_face0],
                        faces[adj_face2],
                    )
                    facet_graph.add_edge(face_id, "v" + str(common_vertices0[0]))
                    facet_graph.add_edge(face_id, "v" + str(common_vertices0[1]))
                    if common_vertices0[0] in common_vertices01:
                        facet_graph.add_edge(adj_face0, "v" + str(common_vertices0[0]))
                    else:
                        facet_graph.add_edge(adj_face0, "v" + str(common_vertices0[1]))
                    facet_graph.add_edge(adj_face1, "v" + str(common_vertices1[0]))
                    facet_graph.add_edge(adj_face1, "v" + str(common_vertices1[1]))
                    assert len(common_vertices12) == 1
                    facet_graph.add_edge(adj_face2, "v" + str(common_vertices12[0]))
                    # Add extra edges, not in slide
                    for v in faces[adj_face0]:
                        if v not in common_vertices0:
                            facet_graph.add_edge(adj_face0, "v" + str(v))

                    facet_graph.add_edge(adj_face2, "v" + str(common_vertices02[0]))

    for face_id in range(len(faces) - 1):
        nei_nodes = list(spanning_tree.adj[face_id])
        for adj_face_id in nei_nodes:
            if not nx.has_path(facet_graph, face_id, adj_face_id):
                common_vertices = np.intersect1d(faces[face_id], faces[adj_face_id])

                if facet_graph.has_edge(face_id, "v" + str(common_vertices[0])):
                    facet_graph.remove_edge(face_id, "v" + str(common_vertices[0]))
                else:
                    facet_graph.add_edge(face_id, "v" + str(common_vertices[0]))
                if facet_graph.has_edge(face_id, "v" + str(common_vertices[1])):
                    facet_graph.remove_edge(face_id, "v" + str(common_vertices[1]))
                else:
                    facet_graph.add_edge(face_id, "v" + str(common_vertices[1]))

                if facet_graph.has_edge(adj_face_id, "v" + str(common_vertices[0])):
                    facet_graph.remove_edge(adj_face_id, "v" + str(common_vertices[0]))
                else:
                    facet_graph.add_edge(adj_face_id, "v" + str(common_vertices[0]))
                if facet_graph.has_edge(adj_face_id, "v" + str(common_vertices[1])):
                    facet_graph.remove_edge(adj_face_id, "v" + str(common_vertices[1]))
                else:
                    facet_graph.add_edge(adj_face_id, "v" + str(common_vertices[1]))

    return facet_graph


def _remove_single_vertices_connected_components(facet_graph: nx.Graph) -> None:
    for component in list(nx.connected_components(facet_graph)):
        component = list(component)
        if len(component) == 1 and component[0].startswith("v"):
            facet_graph.remove_node(component[0])
        else:
            ValueError("The graph is not fully connected")


def _get_polygons_and_offsets(
    eulerian_path: list[tuple],
    faces: np.ndarray,
    vertices: np.ndarray,
) -> tuple:
    face_normals = igl.per_face_normals(vertices, faces, np.ones((1, 3)))

    polygons = []
    offsets = [[0, 0]]
    x_offset = 0
    y_offset = 0

    face_order = []
    vertex_order = []

    vertex_start = False

    if "v" in str(eulerian_path[0][0]):
        vertex_start = True

    for first, second in eulerian_path:
        if "v" not in str(first):
            face_order.append(first)
        else:
            vertex_order.append(int(first[1:]))

    last = eulerian_path[len(eulerian_path) - 1][1]
    if eulerian_path[0][0] != last:
        if "v" not in str(last):
            face_order.append(last)
        else:
            vertex_order.append(int(last[1:]))

    for n, face_index in enumerate(face_order):
        face_coordinates = [vertices[vertex_id] for vertex_id in faces[face_index]]
        for i in range(3):
            face_coordinates[i] = get_2d_projection(face_normals[face_index]).dot(
                face_coordinates[i],
            )
        v1 = None
        v2 = None
        if vertex_start:
            for i, vertex_index in enumerate(faces[face_index]):
                if vertex_index == vertex_order[n]:
                    v1 = i
                if (
                    vertex_index
                    == vertex_order[0 if n + 1 >= len(vertex_order) else n + 1]
                ):
                    v2 = i
        else:
            for i, vertex_index in enumerate(faces[face_index]):
                if (
                    vertex_index
                    == vertex_order[len(vertex_order) - 1 if n - 1 < 0 else n - 1]
                ):
                    v1 = i
                if vertex_index == vertex_order[n]:
                    v2 = i
        assert v1 is not None
        assert v2 is not None
        face_coordinates = rotate_triangle_to_line(face_coordinates, v1, v2)
        face_coordinates = translate_to_origin(face_coordinates, v1)
        for i in range(3):
            face_coordinates[i] = face_coordinates[i] + np.array([x_offset, y_offset])
        polygon = face_coordinates
        x_offset = polygon[v2][0]
        y_offset = polygon[v2][1]
        polygons.append(polygon)
        offsets.append([x_offset, y_offset])

    return polygons, offsets


"""
GRAPH UTILS
This collection of functions provides convenient tools for working with graphs.
"""


def _is_adjacent(current_face: np.ndarray, face: np.ndarray) -> bool:
    """
    Determines whether two faces in a mesh are adjacent.

    Args:
    current_face: A NumPy array representing the first face, containing vertex indices.
    face: A NumPy array representing the second face, containing vertex indices.

    Returns:
    True if the faces are adjacent (share at least two vertices), False otherwise.
    """
    common_vertices = np.intersect1d(current_face, face)
    if len(common_vertices) >= 2:
        return True
    return False


def _add_edges_from_mesh(graph: nx.Graph, faces: np.ndarray) -> None:
    """
    Adds edges to a graph based on face adjacency in a mesh.

    Args:
    graph: A NetworkX graph to which edges will be added.
    faces: A NumPy array of faces, each represented as a list of vertex indices.

    Returns:
    None. Modifies the graph in-place by adding edges representing face adjacency.
    """
    for current_face_id, current_face in enumerate(faces):
        for face_id, face in enumerate(faces):
            if _is_adjacent(current_face, face) and current_face_id != face_id:
                graph.add_edge(current_face_id, face_id)


def _graph_from_mesh(faces: np.ndarray) -> nx.Graph:
    """
    Creates a graph representation of a mesh based on face adjacency.

    Args:
    faces: A NumPy array of faces, each represented as a list of vertex indices.

    Returns:
    A NetworkX graph where nodes represent faces and edges represent adjacency between faces.
    """
    graph = nx.Graph()
    _add_edges_from_mesh(graph, faces)
    return graph


"""
TRANSFORMATION UTILS
This set of functions offers user-friendly tools for handling transformations like translation and rotation.
"""


def get_rotation_matrix(axis: np.ndarray, angle: np.float64) -> np.ndarray:
    return scipy.spatial.transform.Rotation.from_rotvec(
        axis / np.linalg.norm(axis) * angle,
    ).as_matrix()


def get_2d_projection(face_normal: np.ndarray) -> np.ndarray:
    xy_plane_normal = np.array([0, 0, 1])
    rotation_axis = np.cross(face_normal, xy_plane_normal)
    angle = np.arccos(np.clip(np.dot(xy_plane_normal, face_normal), -1.0, 1.0))

    discard_z_matrix = np.array(
        [
            [1, 0, 0],
            [0, 1, 0],
        ],
    )
    if np.linalg.norm(rotation_axis) == 0:
        return discard_z_matrix
    rotation_matrix = get_rotation_matrix(rotation_axis, angle)
    return discard_z_matrix.dot(rotation_matrix)


def translate_to_origin(
    face_coordinates: np.ndarray,
    index_to_origin: int,
) -> np.ndarray:
    translation_vector = face_coordinates[index_to_origin]
    for i in range(3):
        face_coordinates[i] = face_coordinates[i] - translation_vector
    return face_coordinates


def rotate_triangle_to_line(
    face_coordinates,
    common_vertex1,
    common_vertex2,
):
    c1 = face_coordinates[common_vertex1][1] - face_coordinates[common_vertex2][1]
    c2 = face_coordinates[common_vertex1][0] - face_coordinates[common_vertex2][0]
    if c2 == 0:
        for i in range(3):
            face_coordinates[i] = np.array(
                [
                    [0, 1],
                    [-1, 0],
                ],
            ).dot(face_coordinates[i])
    else:
        rotation_matrix = get_rotation_matrix(
            np.array([0.0, 0.0, 1.0]),
            -math.atan(c1 / c2),
        )
        discard_z_matrix = np.array(
            [
                [1, 0, 0],
                [0, 1, 0],
            ],
        )
        for i in range(3):
            face_coordinates[i] = np.delete(
                discard_z_matrix.dot(rotation_matrix),
                2,
                1,
            ).dot(face_coordinates[i])

    if face_coordinates[common_vertex1][0] > face_coordinates[common_vertex2][0]:
        for i in range(3):
            face_coordinates[i] = np.array(
                [
                    [-1, 0],
                    [0, -1],
                ],
            ).dot(face_coordinates[i])

    for index in [0, 1, 2]:
        if index not in [common_vertex1, common_vertex2]:
            not_common_vertex = index

    assert not_common_vertex is not None

    if face_coordinates[not_common_vertex][0] < face_coordinates[common_vertex1][0]:
        angle = (
            face_coordinates[not_common_vertex][1] - face_coordinates[common_vertex1][1]
        ) / (
            face_coordinates[not_common_vertex][0] - face_coordinates[common_vertex1][0]
        )
        side = (
            1.5
            if face_coordinates[not_common_vertex][1]
            > face_coordinates[common_vertex1][1]
            else 0.5
        )
        rotation_matrix = get_rotation_matrix(
            np.array([0.0, 0.0, 1.0]),
            math.pi * side - math.atan(angle),
        )
        discard_z_matrix = np.array(
            [
                [1, 0, 0],
                [0, 1, 0],
            ],
        )
        for i in range(3):
            face_coordinates[i] = np.delete(
                discard_z_matrix.dot(rotation_matrix),
                2,
                1,
            ).dot(face_coordinates[i])

    if face_coordinates[not_common_vertex][0] > face_coordinates[common_vertex2][0]:
        angle = (
            face_coordinates[not_common_vertex][1] - face_coordinates[common_vertex2][1]
        ) / (
            face_coordinates[not_common_vertex][0] - face_coordinates[common_vertex2][0]
        )
        side = (
            0.5
            if face_coordinates[not_common_vertex][1]
            > face_coordinates[common_vertex2][1]
            else 1.5
        )
        rotation_matrix = get_rotation_matrix(
            np.array([0.0, 0.0, 1.0]),
            math.pi * side - math.atan(angle),
        )
        discard_z_matrix = np.array(
            [
                [1, 0, 0],
                [0, 1, 0],
            ],
        )
        for i in range(3):
            face_coordinates[i] = np.delete(
                discard_z_matrix.dot(rotation_matrix),
                2,
                1,
            ).dot(face_coordinates[i])

    return face_coordinates
