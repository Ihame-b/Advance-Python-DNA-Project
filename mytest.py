import json

import networkx as nx
import pandas as pd
from pytest import mark

from project import (
    clean_data,
    construct_dna_sequence,
    construct_graph,
    generate_sequences,
    is_valid_graph,
)
COLS = ['SegmentNr', 'Position', 'A', 'C', 'G', 'T']
def _df(rows):
    return pd.DataFrame(rows, columns=COLS)

''' Testing Clean data'''

@mark.parametrize(
    'raw, expected_rows',
    [
        (
            _df([[4, 1, 0, 0, 1, 1], [4, 2, 1, 0, 0, 0]]),
            [],
        ),
        (
            _df([
                [9, 1, 1, 0, 0, 0], [9, 2, 0, 1, 0, 0], [9, 4, 0, 0, 0, 1],
                [10, 1, 0, 0, 1, 0], [10, 2, 0, 0, 0, 1],
            ]),
            [[10, 1, 0, 0, 1, 0], [10, 2, 0, 0, 0, 1]],
        ),
        (
            _df([
                [6, 1, 1, 0, 0, 0], [6, 1, 1, 0, 0, 0], [6, 2, 0, 0, 1, 0],
            ]),
            [[6, 1, 1, 0, 0, 0], [6, 2, 0, 0, 1, 0]],
        ),
        (
            _df([
                [11, 1, 0, 0, 1, 0], [11, 2, 0, 0, 0, 1],
                [22, 1, 0, 0, 1, 0], [22, 2, 0, 0, 0, 1],
            ]),
            [[11, 1, 0, 0, 1, 0], [11, 2, 0, 0, 0, 1]],
        ),
    ],
)
def test_clean_data_custom(raw, expected_rows):
    got = clean_data(raw).reset_index(drop=True)
    want = _df(expected_rows).reset_index(drop=True)
    if want.empty:
        assert got.empty
    else:
        assert got.equals(want)

''' Generate the sequence test'''
@mark.parametrize(
    'raw, expected_json',
    [
        (_df([[5, 1, 0, 0, 0, 1], [5, 2, 0, 0, 1, 0]]), {"5": "TG"}),
        (
            _df([
                [1, 1, 0, 1, 0, 0], [1, 2, 0, 1, 0, 0],
                [2, 1, 0, 0, 1, 0], [2, 2, 0, 0, 1, 0],
            ]),
            {"1": "CC", "2": "GG"},
        ),
        (
            _df([
                [8, 1, 0, 0, 1, 0], [8, 2, 1, 0, 0, 0],
                [8, 3, 0, 1, 0, 0], [8, 4, 0, 0, 0, 1],
            ]),
            {"8": "GACT"},
        ),
    ],
)
def test_generate_sequences_custom(raw, expected_json):
    assert json.loads(generate_sequences(raw)) == expected_json


''' Test construction of graph'''
@mark.parametrize(
    'payload, k, edges',
    [
        ({"1": "CAGCT"}, 4, [("AGC", "GCT"), ("CAG", "AGC")]),
        ({"1": "TAG", "2": "GAT"}, 2, [("A", "G"), ("A", "T"), ("G", "A"), ("T", "A")]),
        ({"1": "CGCAT"}, 3, [("CA", "AT"), ("CG", "GC"), ("GC", "CA")]),
    ],
)
def test_construct_graph_custom(payload, k, edges):
    g = construct_graph(json.dumps(payload), k)
    assert sorted(g.edges()) == sorted(edges)

''' Testing if graph is valid '''
@mark.parametrize(
    'edges, ok',
    [
        ([("GT", "TG"), ("TG", "GT")], True),                      # 2-node cycle
        ([("CA", "AC"), ("AC", "CG"), ("CG", "GC")], True),        # euler path
        ([("AA", "AB"), ("CD", "DE")], False),                     # disconnected
        ([("M", "N"), ("M", "O"), ("P", "M")], False),             # 3 imbalance nodes
    ],
)
def test_is_valid_graph_custom(edges, ok):
    g = nx.MultiDiGraph()
    g.add_edges_from(edges)
    assert is_valid_graph(g) is ok

''' Testing construction of sequence'''
@mark.parametrize(
    'edges, expected',
    [
        ([("CA", "AC"), ("AC", "CG"), ("CG", "GC")], "CACGC"),
        ([("GT", "TG"), ("TG", "GT")], "GTG"),
        ([("XY", "YZ")], "XYZ"),
        ([("NO", "PE"), ("QR", "ST")], "DNA sequence can not be constructed."),
    ],
)
def test_construct_dna_sequence_custom(edges, expected):
    g = nx.MultiDiGraph()
    g.add_edges_from(edges)
    assert construct_dna_sequence(g) == expected