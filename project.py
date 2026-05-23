"""
author: Ihame Gilbert
"""

import json
import os
import re
import sys

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

''' DNA Sequence Assembly Project of Advanced Python Programming project. '''

''' Load CSV into DataFrame. '''
def read_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(
        name,
        header=None,
        names=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']
    )

'''Converting Segments to JSON.'''
def _segment_to_sequence(DataFramSeg: pd.DataFrame) -> str:
    """Build DNA string from one segment's rows (sorted by Position)."""
    DataFramSeg = DataFramSeg.sort_values('Position')
    chars = []
    for _, row in DataFramSeg.iterrows():
        for base in ['A', 'C', 'G', 'T']:
            if row[base] == 1:
                chars.append(base)
                break
    return ''.join(chars)

'''Clean Data by removing invalid segments.'''
def clean_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    remove_seg = set()
    for segment in df['SegmentNr'].unique():
        seg = df[df['SegmentNr'] == segment]

        '''Wrong position (not exactly one base = 1) -> drop segment'''
        if not ((seg[['A', 'C', 'G', 'T']].sum(axis=1) == 1).all()):
            remove_seg.add(segment)
            continue
        
        '''Missing position (not 1 .. max) -> drop segment'''
        positions = set(seg['Position'].astype(int))
        max_pos = int(seg['Position'].max())
        if positions != set(range(1, max_pos + 1)):
            remove_seg.add(segment)
            continue

        '''Duplicate position: same values -> keep one row , conflict -> drop segment.'''
        for pos in seg['Position'].unique():
            rows_at_pos = seg[seg['Position'] == pos]
            if len(rows_at_pos) > 1:
                if len(rows_at_pos.drop_duplicates(subset=['A', 'C', 'G', 'T'])) > 1:
                    remove_seg.add(segment)
                    break
    df = df[~df['SegmentNr'].isin(remove_seg)]

    '''Duplicate rows at same position with same values — keep one'''
    df = df.drop_duplicates(subset=['SegmentNr', 'Position'], keep='first')

    '''Duplicate segments (same sequence) -> keep one segment'''
    seen_sequences = {}
    duplicated_seg = set()
    for segment in df['SegmentNr'].unique():
        seg = df[df['SegmentNr'] == segment]
        seq = _segment_to_sequence(seg)
        if seq in seen_sequences.values():
            duplicated_seg.add(segment)
        else:
            seen_sequences[segment] = seq

    df = df[~df['SegmentNr'].isin(duplicated_seg)]
    return df.reset_index(drop=True)

'''Convert one segment's rows into a DNA string'''
def _segment_to_sequence(DataFramSeg: pd.DataFrame) -> str:
    
    DataFramSeg = DataFramSeg.sort_values('Position')
    chars = []
    for _, row in DataFramSeg.iterrows():
        for base in ['A', 'C', 'G', 'T']:
            if row[base] == 1:
                chars.append(base)
                break
    return ''.join(chars)

'''Convert cleaned dataframe to JSON( eg: segment id ->DNA string {"1": "ATGAA", "2": "CTGAATGA"})'''
def generate_sequences(df: pd.DataFrame) -> str:

    sequences = {}
    for segment in sorted(df['SegmentNr'].unique()):
        seg = df[df['SegmentNr'] == segment]
        sequences[str(int(segment))] = _segment_to_sequence(seg)
    return json.dumps(sequences)

'''Creating a Bruijn graph as a NetworkX MultiDiGraph'''
def construct_graph(json_data: str, k: int) -> nx.MultiDiGraph:
    
    G = nx.MultiDiGraph()
    sequences = json.loads(json_data)
    for seq in sequences.values():
        if len(seq) < k:
            continue 
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            left = kmer[:-1]  
            right = kmer[1:]  
            G.add_edge(left, right)

    return G

'''Plot de Bruijn graph and save as PNG'''
def plot_graph(graph: nx.MultiDiGraph, filename: str) -> None:

    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(graph)
    nx.draw(
        graph, pos,
        with_labels=True,
        node_color='lightblue',
        node_size=500,
        font_size=8,
    )
    plt.savefig(filename)
    plt.close()

'''Check if de Bruijn graph can be sequenced'''
def is_valid_graph(graph: nx.MultiDiGraph) -> bool:
    
    if graph.number_of_nodes() == 0:
        return False
    if not nx.is_weakly_connected(graph):
        return False

    in_deg = dict(graph.in_degree())
    out_deg = dict(graph.out_degree())
    diff_nodes = [n for n in graph.nodes if in_deg[n] != out_deg[n]]

    if len(diff_nodes) not in (0, 2):
        return False
    if len(diff_nodes) == 2:
        start_ok = any(out_deg[n] - in_deg[n] == 1 for n in diff_nodes)
        end_ok = any(in_deg[n] - out_deg[n] == 1 for n in diff_nodes)
        return start_ok and end_ok
    return True    

'''Build DNA string from Euler path'''
def construct_dna_sequence(graph: nx.MultiDiGraph) -> str:
    
    if not is_valid_graph(graph):
        return "DNA sequence can not be constructed."
    G = graph.copy()
    in_d, out_d = dict(G.in_degree()), dict(G.out_degree())
    diff = [n for n in G.nodes if in_d[n] != out_d[n]]

    if len(diff) == 2:
        start = next(n for n in diff if out_d[n] - in_d[n] == 1)
        G.add_edge(next(n for n in diff if in_d[n] - out_d[n] == 1), start)
    else:
        start = next(iter(G.nodes))

    stack, path = [start], []
    while stack:
        v = stack[-1]
        if G.out_degree(v):
            w = next(iter(G.successors(v)))
            stack.append(w)
            G.remove_edge(v, w)
        else:
            path.append(stack.pop())

    path = path[::-1]
    if len(path) > 1 and path[0] == path[-1]:
        path = path[:-1]

    print(" -> ".join(path))
    if len(path) <= 1:
        return "DNA sequence can not be constructed."
    return path[0] + "".join(n[-1] for n in path[1:])

def save_output(s: str, filename: str) -> None:
    with open(filename, 'w') as f:
        f.write(s)

'''Mai Run the full DNA assembly pipeline: read -> clean -> JSON -> graph -> plot -> validate -> DNA -> save'''
def main() -> None:
    
    if len(sys.argv) != 2:
        print("Usage: python project.py DNA_x_k.csv")
        return
    csv_file = sys.argv[1]
    basename = os.path.basename(csv_file)
    m = re.match(r'DNA_(\d+)_(\d+)\.csv', basename)
    if not m:
        print("Invalid filename format")
        return
    x, k = int(m.group(1)), int(m.group(2))
    df = clean_data(read_csv(csv_file))  # full path still works
    json_data = generate_sequences(df)
    graph = construct_graph(json_data, k)
    plot_graph(graph, f"DNA_{x}.png")
    if is_valid_graph(graph):
        sequence = construct_dna_sequence(graph)
    else:
        sequence = "DNA sequence can not be constructed."
    save_output(sequence, f"DNA_{x}.txt")

if __name__ == "__main__":
    main()