"""
Tests for project.py — DNA Sequence Assembly
"""

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


@mark.parametrize(
    'dna_df, expected',
    [
        (
            pd.DataFrame(
                data=[
                    [1, 1, 1, 0, 0, 1],
                    [1, 2, 0, 0, 0, 1],
                    [2, 1, 1, 0, 0, 0],
                    [2, 2, 0, 1, 0, 0],
                ],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'],
            ),
            pd.DataFrame(
                data=[
                    [2, 1, 1, 0, 0, 0],
                    [2, 2, 0, 1, 0, 0],
                ],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'],
            ),
        ),
        (
            pd.DataFrame(
                data=[
                    [1, 1, 1, 0, 0, 0],
                    [1, 3, 0, 0, 0, 1],
                    [2, 1, 0, 1, 0, 0],
                    [2, 2, 0, 0, 1, 0],
                ],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'],
            ),
            pd.DataFrame(
                data=[
                    [2, 1, 0, 1, 0, 0],
                    [2, 2, 0, 0, 1, 0],
                ],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'],
            ),
        ),
        (
            pd.DataFrame(
                data=[
                    [1, 1, 1, 0, 0, 0],
                    [1, 2, 0, 0, 0, 1],
                    [1, 2, 0, 1, 0, 0],
                ],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'],
            ),
            pd.DataFrame(
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'],
            ),
        ),
    ],
)
def test_clean_data(dna_df: pd.DataFrame, expected: pd.DataFrame) -> None:
    result = clean_data(dna_df).reset_index(drop=True)
    expected = expected.reset_index(drop=True)
    if expected.empty:
        assert result.empty
    else:
        assert result.equals(expected)


@mark.parametrize(
    'dna_df, expected_json_str',
    [
        (
            pd.DataFrame(
                data=[
                    [1, 1, 1, 0, 0, 0],
                    [1, 2, 0, 0, 0, 1],
                ],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'],
            ),
            json.dumps({"1": "AT"}),
        ),
        (
            pd.DataFrame(
                data=[
                    [1, 1, 1, 0, 0, 0],
                    [1, 2, 0, 0, 0, 1],
                    [2, 1, 0, 1, 0, 0],
                    [2, 2, 0, 0, 1, 0],
                ],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'],
            ),
            json.dumps({"1": "AT", "2": "CG"}),
        ),
        (
            pd.DataFrame(
                data=[
                    [1, 1, 1, 0, 0, 0],
                    [1, 2, 0, 0, 0, 1],
                    [1, 3, 0, 1, 0, 0],
                ],
                columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'],
            ),
            json.dumps({"1": "ATC"}),
        ),
    ],
)
def test_generate_sequences(dna_df: pd.DataFrame, expected_json_str: str) -> None:
    assert generate_sequences(dna_df) == expected_json_str


@mark.parametrize(
    'json_data, k, expected_edge_list',
    [
        (
            json.dumps({"1": "ATTACTC"}),
            5,
            [('ATTA', 'TTAC'), ('TTAC', 'TACT'), ('TACT', 'ACTC')],
        ),
        (
            json.dumps({"1": "ATGAA", "2": "CTGAATGA"}),
            5,
            [
                ('ATGA', 'TGAA'),
                ('CTGA', 'TGAA'),
                ('TGAA', 'GAAT'),
                ('GAAT', 'AATG'),
                ('AATG', 'ATGA'),
            ],
        ),
        (
            json.dumps({"1": "ATCG"}),
            3,
            [('AT', 'TC'), ('TC', 'CG')],
        ),
    ],
)
def test_construct_graph(json_data: str, k: int, expected_edge_list: list) -> None:
    graph = construct_graph(json_data, k)
    assert sorted(list(graph.edges())) == sorted(expected_edge_list)


@mark.parametrize(
    'DNA_edge_list, expected_validity',
    [
        (
            [('ATTA', 'TTAC'), ('TTAC', 'TACT'), ('TACT', 'ACTC'), ('ACTC', 'ATTA')],
            True,
        ),
        (
            [('AAA', 'AAC'), ('AAC', 'ACA'), ('ACA', 'CAC')],
            True,
        ),
        (
            [('ATTA', 'TTAC'), ('TTAC', 'TACT'), ('ACTC', 'CTCA')],
            False,
        ),
        (
            [('A', 'B'), ('A', 'C'), ('A', 'D')],
            False,
        ),
    ],
)
def test_is_valid_graph(DNA_edge_list: list, expected_validity: bool) -> None:
    graph = nx.MultiDiGraph()
    for edge in DNA_edge_list:
        graph.add_edge(edge[0], edge[1])
    assert is_valid_graph(graph) is expected_validity


@mark.parametrize(
    'DNA_edge_list, possible_dna_sequence',
    [
        (
            [('AAA', 'AAC'), ('AAC', 'ACA'), ('ACA', 'CAC')],
            ["AAACAC"],
        ),
        (
            [('ATTA', 'TTAC'), ('TTAC', 'TACT'), ('TACT', 'ACTC'), ('ACTC', 'ATTA')],
            ["ATTACTC"],
        ),
        (
            [('AT', 'TG'), ('TG', 'GA')],
            ["ATGA"],
        ),
    ],
)
def test_construct_dna_sequence(DNA_edge_list: list, possible_dna_sequence: list) -> None:
    graph = nx.MultiDiGraph()
    for edge in DNA_edge_list:
        graph.add_edge(edge[0], edge[1])
    result = construct_dna_sequence(graph)
    assert result in possible_dna_sequence
