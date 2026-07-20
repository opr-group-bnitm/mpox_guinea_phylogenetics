import argparse
import pandas as pd
from collections import defaultdict, Counter
from Bio import SeqIO
import baltic as bt


# Assign SNPs to branches
# This is based on code extracted from squirrel


def get_state(nt):
    """Get A, C, G or T if the position calls one of these, else 'N'
    """
    nt_upper = nt.upper()
    if nt_upper in 'ACGT':
        return nt_upper
    return 'N'


def read_node_states(fname_alignment, fname_states):
    """Get a nested dictionairy with node + site -> state (=Nucleobase)
    """
    node_states = defaultdict(dict)

    # read the state file, skip the comments and the header
    with open(fname_states) as f:
        entries = [
            l.rstrip("\n").split("\t")[:3]
            for l in f
            if not l.startswith('#') and l.strip()
        ][1:]

    for node, site, state in entries:
        node_states[int(site)-1][node] = get_state(state)

    # now the tips
    for record in SeqIO.parse(fname_alignment, 'fasta'):
        for i, nt in enumerate(record.seq):
            node_states[i][record.id] = get_state(nt)

    return node_states


def get_snp_states(node_states):
    """Get all columns in the ancestral MSA that has unambiguous variance.
    """
    return {
        site: states
        for site, states in node_states.items()
        if len(set(state for state in states.values() if state != 'N')) > 1
    }


def get_context(node_states, node, site):
    pre =  'N' if site == 0 else node_states[site - 1][node]
    nt = node_states[site][node]
    post = node_states[site+1][node] if (site + 1) in node_states else 'N'
    return f'{pre}{nt}{post}'


def get_snps_two_branches(parent, child, node_states, snp_states):
    return [
        (parent, child, i, states[parent], states[child], get_context(node_states, parent, i), get_context(node_states, child, i))
        for i, states in snp_states.items()
        if states[parent] != states[child] and states[parent] != 'N' and states[child] != 'N'
    ]


def get_snps_along_branches(fname_tree, node_states, snp_states):
    snps = []
    tree = bt.loadNewick(fname_tree, absoluteTime=False)
    for parent in tree.getBranches():
        if not parent.is_leaf():
            parent_label = parent.traits['label']
            for child in parent.children:
                child_label = child.name if child.is_leaf() else child.traits['label']
                snps.extend(get_snps_two_branches(parent_label, child_label, node_states, snp_states))
    return snps


def get_apobec_label(parent_state, child_state, context_parent, context_child):
    snp = (parent_state, child_state)
    if snp == ('C', 'T'):
        if context_parent[0] == 'T':
            if context_child[0] == 'T':
                return 'apobec'
            else:
                return 'unclear'
        elif context_parent[0] == 'N':
            return 'unclear'
    elif snp == ('G', 'A'):
        if context_parent[2] == 'A':
            if context_child[2] == 'A':
                return 'apobec'
            else:
                return 'unclear'
        elif context_parent[2] == 'N':
            return 'unclear'
    return 'non-apobec'


def count_apobec_mutations(branch_snps):
    """Iterate through all branches of a tree and collect the SNPs with contexts.
    """
    snps_with_labels = [
        (parent_id, child_id, site, parent_state, child_state, context_parent, context_child,
         get_apobec_label(parent_state, child_state, context_parent, context_child))
        for parent_id, child_id, site, parent_state, child_state, context_parent, context_child in branch_snps
    ]
    apobec_counter = Counter(
        (parent_id, child_id, label)
        for parent_id, child_id, _, _, _, _, _, label in snps_with_labels
    )
    apobec_labels = ['apobec', 'non-apobec', 'unclear']
    node_pairs = sorted(set((parent, child) for parent, child, _ in apobec_counter.keys()))
    apobed_counts = [
        {
            'parent': parent_id,
            'child': child_id,
            **{label: apobec_counter[(parent_id, child_id, label)] for label in apobec_labels},
        }
        for parent_id, child_id in node_pairs
    ]
    return apobed_counts, snps_with_labels


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--msa', required=True, help='Multile sequence alignment fasta file')
    parser.add_argument('--ancestral', required=True, help='IQTREE ancestral reconstruction state file.')
    parser.add_argument('--tree', required=True, help='IQTREE phylogenetic tree in newick format.')
    parser.add_argument('--output', required=True, help='Output .csv file.')
    parser.add_argument('--table', required=True, help='Filename for .csv file with per mutation table.')
    args = parser.parse_args()

    node_states = read_node_states(args.msa, args.ancestral)
    snp_states = get_snp_states(node_states)
    snps_on_branches = get_snps_along_branches(args.tree, node_states, snp_states)
    apobec_counts, apobec_snps = count_apobec_mutations(snps_on_branches)

    if args.table:
        cols = [
            'parent',
            'child',
            'site',
            'ref',
            'alt',
            'context_parent',
            'context_child',
            'apobec_state',
        ]
        pd.DataFrame(apobec_snps, columns=cols).to_csv(args.table, index=False)

    pd.DataFrame(apobec_counts).to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
