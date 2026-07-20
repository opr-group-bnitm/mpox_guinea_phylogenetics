import argparse
import baltic as bt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input',
        required=True,
        help='Newick tree with labels in format NODE_BOOTSTRAPSUPPORT'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Name for output tree that can be read by baltic and has traits label and bootstrap'
    )
    args = parser.parse_args()

    tree = bt.loadNewick(args.input)
    for node in tree.getBranches():
        if not node.is_leaf() and 'label' in node.traits and '_' in node.traits['label']:
            node.traits['label'], node.traits['bootstrap'] = node.traits['label'].split('_')

    with open(args.output, 'w') as f_out:
        f_out.write(tree.toString())


if __name__ == '__main__':
    main()
